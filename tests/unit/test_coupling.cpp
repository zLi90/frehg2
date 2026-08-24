/// \file test_coupling.cpp
/// \brief Unit tests for the surface-subsurface coupler (plan §8.1
///        test_massbalance_coupled and §10 P3): miniature coupled
///        configurations pinning the per-step volume closure, the §5.7
///        exchange resolution (infiltration limited to available surface
///        water, no porosity factor), the dry-cell seepage hold, the
///        saturation bounce-back, the reallocation vent, the coupled
///        evaporation correction, and the sync-mode adaptive lockstep.

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Types.hpp"
#include "coupling/Coupler.hpp"
#include "gw/RichardsSolver.hpp"
#include "gw/TerrainMetric.hpp"
#include "swe/SurfaceSolver.hpp"

#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <memory>
#include <vector>

namespace {

using frehg::real_t;

/// A single-rank miniature coupled run assembled from an in-code
/// configuration, mirroring the driver's construction order.
struct MiniCoupled {
  frehg::FrehgConfig cfg;
  std::unique_ptr<frehg::Grid> grid;
  std::unique_ptr<frehg::HaloExchanger> halo;
  std::unique_ptr<frehg::gw::TerrainMetric> mesh;
  std::unique_ptr<frehg::BoundarySet> boundaries;
  std::unique_ptr<frehg::swe::SurfaceSolver> surface;
  std::unique_ptr<frehg::gw::RichardsSolver> gw;
  std::unique_ptr<frehg::coupling::Coupler> coupler;
  real_t t = 0.0;

  void build() {
    grid = std::make_unique<frehg::Grid>(MPI_COMM_WORLD, cfg.domain);
    halo = std::make_unique<frehg::HaloExchanger>(*grid, false);
    mesh = std::make_unique<frehg::gw::TerrainMetric>(*grid, cfg, *halo);
    grid->buildGlobalIds(mesh->ktop());
    boundaries =
        std::make_unique<frehg::BoundarySet>(*grid, cfg.boundaryConditions, cfg.configDir);
    surface = std::make_unique<frehg::swe::SurfaceSolver>(*grid, cfg, *boundaries, *halo);
    gw = std::make_unique<frehg::gw::RichardsSolver>(*grid, cfg, *mesh, *boundaries, *halo);
    coupler = std::make_unique<frehg::coupling::Coupler>(*grid, cfg, *surface, *gw);
  }

  /// One coupled step at the coupler's own step size; returns the dt taken.
  real_t step() {
    const real_t dt = coupler->nextDt();
    t += dt;
    coupler->step(t, dt);
    return dt;
  }

  real_t surfaceVolume() const { return surface->ownedVolume(); }
  real_t gwVolume() const { return gw->ownedVolume(); }

  /// The coupled per-step closure identity (plan §5.7 / §8.1): every term is
  /// measured, so the residual must sit at rounding level.
  real_t closureResidual(real_t vSurfBefore, real_t vGwBefore, real_t accumBefore) const {
    const frehg::coupling::CouplingAudit& c = coupler->audit();
    const frehg::swe::SurfaceStepAudit& s = surface->audit();
    const real_t lhs = (surfaceVolume() - vSurfBefore) + (gwVolume() - vGwBefore) +
                       (c.accumVolume - accumBefore);
    const real_t rhs = (s.rainVolume - s.evapVolume - s.boundaryOutflow + s.bcInflow +
                        s.clampVolume) +
                       (c.gw.boundaryIn + c.gw.cplExchanged + c.gw.cplVent) + c.gw.cplBounce -
                       c.gw.cplEvap + c.discarded - c.gw.ssStorage + c.gw.reallocAdjust -
                       c.gw.vloss;
    return lhs - rhs;
  }
};

/// The soft test soil of test_gw.cpp.
frehg::SoilType testSoil(real_t ks = 1.0e-5) {
  frehg::SoilType soil;
  soil.name = "test";
  soil.ksx = ks;
  soil.ksy = ks;
  soil.ksz = ks;
  soil.thetaS = 0.46;
  soil.thetaR = 0.04;
  soil.vgAlpha = 5.9;
  soil.vgN = 2.68;
  soil.aev = -0.02;
  return soil;
}

/// A coupled base configuration: nx columns of nz cells under a flat bed at
/// z = 0, fixed lockstep dt (dt_min = dt_max = dt), quiescent surface.
frehg::FrehgConfig coupledBaseConfig(int nx, int nz, real_t dz, real_t dt) {
  frehg::FrehgConfig cfg;
  cfg.simulation.id = "coupled-unit";
  cfg.domain.nx = nx;
  cfg.domain.ny = 1;
  cfg.domain.nz = nz;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 1.0;
  cfg.domain.dz = dz;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = dt;
  cfg.time.tEnd = 1000.0 * dt;
  cfg.time.outputInterval = 1000.0 * dt;
  cfg.modules.surfaceWater = true;
  cfg.modules.groundwater = true;
  cfg.coupling.mode = frehg::CouplingConfig::Mode::Sync;
  cfg.surfaceWater.gravity = 9.81;
  cfg.surfaceWater.friction.coefficient.constant = 0.02;
  cfg.surfaceWater.minDepth = 1.0e-6;
  cfg.surfaceWater.wettingFaceDepth = 1.0e-6;
  cfg.surfaceWater.friction.thinLayerDepth = 0.01;
  cfg.groundwater.timestep.dtInit = dt;
  cfg.groundwater.timestep.dtMin = dt;
  cfg.groundwater.timestep.dtMax = dt;
  cfg.groundwater.specificStorage = 0.0;
  cfg.soil.types = {testSoil()};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.surface.eta.constant = 0.0;
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  cfg.initialConditions.groundwater.value.constant = 0.2;
  return cfg;
}

std::vector<std::array<real_t, 2>> wholeFootprint(int nx, int ny) {
  const real_t x1 = static_cast<real_t>(nx) + 0.1;
  const real_t y1 = static_cast<real_t>(ny) + 0.1;
  return {{-0.1, -0.1}, {x1, -0.1}, {x1, y1}, {-0.1, y1}};
}

frehg::BoundaryConditionConfig fluxBc(frehg::BcTarget target, real_t value, int nx, int ny) {
  frehg::BoundaryConditionConfig bc;
  bc.name = "flux";
  bc.polygon = wholeFootprint(nx, ny);
  bc.target = target;
  bc.kind = frehg::BcKind::Flux;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = value;
  return bc;
}

real_t interior2(const frehg::Field2<real_t>& field, int j, int i) {
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  return host(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
}

real_t interior3(const frehg::Field3<real_t>& field, int j, int i, int k) {
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  return host(static_cast<std::size_t>(j), static_cast<std::size_t>(i),
              static_cast<std::size_t>(k));
}

// ---------------------------------------------------------------------------
// The coupled mass audit (plan §5.7 enforcement, §8.1
// test_massbalance_coupled).
// ---------------------------------------------------------------------------

TEST(CoupledModule, PondedInfiltrationClosesMassBudget) {
  // The §8.1 toy: ponded columns infiltrating, no rain, no boundary
  // conditions. Per-step global closure of the coupled budget at 1e-10 of
  // the total volume, and the exchange itself cancels exactly: the volume
  // the surface lost is the volume the subsurface's top faces carried.
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(3, 10, 0.05, 0.5);
  sim.cfg.initialConditions.surface.eta.constant = 0.02;  // 2 cm pond
  sim.build();
  const real_t vTotal = sim.surfaceVolume() + sim.gwVolume();
  EXPECT_GT(sim.surfaceVolume(), 0.0);
  real_t vs = sim.surfaceVolume();
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  for (int n = 0; n < 60; ++n) {
    sim.step();
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    // §5.7: no porosity factor anywhere in the exchange — what left the
    // surface arrived below (no bounce/vent/hold in this toy).
    EXPECT_EQ(c.gw.cplBounce, 0.0);
    EXPECT_EQ(c.gw.cplVent, 0.0);
    EXPECT_NEAR(c.applied, c.gw.cplExchanged, 1.0e-12 * vTotal);
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  // Water genuinely infiltrated.
  EXPECT_LT(sim.surfaceVolume(), 0.06);
  EXPECT_GT(sim.surfaceVolume(), 0.0);
  EXPECT_GT(interior3(sim.gw->waterContent(), 1, 2, 0), 0.2);
}

TEST(CoupledModule, SubcycledWindowClosesMassBudget) {
  // The same toy in subcycled mode with dtg < dt: several subsurface steps
  // per surface step, the exchange accumulated as volume (amendment A9).
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(3, 10, 0.05, 1.0);
  sim.cfg.coupling.mode = frehg::CouplingConfig::Mode::Subcycled;
  sim.cfg.groundwater.timestep.dtInit = 0.2;
  sim.cfg.groundwater.timestep.dtMin = 0.2;
  sim.cfg.groundwater.timestep.dtMax = 0.2;
  sim.cfg.initialConditions.surface.eta.constant = 0.02;
  sim.build();
  const real_t vTotal = sim.surfaceVolume() + sim.gwVolume();
  real_t vs = sim.surfaceVolume();
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  for (int n = 0; n < 30; ++n) {
    EXPECT_EQ(sim.step(), 1.0);  // the surface dt stays fixed in this mode
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    EXPECT_EQ(c.substeps, 5);
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  EXPECT_GT(interior3(sim.gw->waterContent(), 1, 2, 0), 0.2);
}

TEST(CoupledModule, InfiltrationLimitedToAvailableWater) {
  // A thin pond over conductive soil: the unrestricted Darcy flux would
  // drain several times the pond in one step. The §5.7-resolved limit
  // (|q| dtg <= available depth, no wcs factor) hands the subsurface
  // exactly the pond and lands the surface on the bed — nothing is created
  // or discarded.
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(1, 10, 0.05, 0.5);
  sim.cfg.soil.types = {testSoil(1.0e-2)};
  sim.cfg.initialConditions.surface.eta.constant = 1.0e-3;
  sim.build();
  const real_t pond = sim.surfaceVolume();
  const real_t vg0 = sim.gwVolume();
  const real_t vTotal = pond + vg0;
  ASSERT_NEAR(pond, 1.0e-3, 1.0e-12);
  sim.step();
  const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
  EXPECT_NEAR(sim.surfaceVolume(), 0.0, 1.0e-12);
  // The exchange handed the subsurface exactly the pond: the limited
  // top-face flux equals the surface loss (the consistency restore may add
  // separately audited θ, so the closure identity is the conservation
  // statement).
  EXPECT_NEAR(c.applied, -pond, 1.0e-12);
  EXPECT_NEAR(c.gw.cplExchanged, -pond, 1.0e-12);
  EXPECT_LE(std::fabs(c.discarded), 1.0e-12);
  EXPECT_LE(std::fabs(sim.closureResidual(pond, vg0, 0.0)), 1.0e-10 * vTotal);
}

TEST(CoupledModule, SeepageHeldOnDrySurfaceUntilMinDepth) {
  // Exfiltration onto a dry surface: a high prescribed bottom head drives
  // upward flow through a nearly saturated column. The per-step seepage is
  // below min_depth, so the accumulator holds (legacy reset_seepage,
  // shallowwater.c:668-675) until the accumulated depth exceeds it.
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(1, 5, 0.05, 0.5);
  sim.cfg.surfaceWater.minDepth = 2.0e-4;
  sim.cfg.initialConditions.groundwater.value.constant = 0.459;
  frehg::BoundaryConditionConfig head;
  head.name = "bottom-head";
  head.polygon = wholeFootprint(1, 1);
  head.target = frehg::BcTarget::GroundwaterBottom;
  head.kind = frehg::BcKind::Head;
  head.value.form = frehg::BcValueConfig::Form::Constant;
  head.value.constant = 2.0;
  sim.cfg.boundaryConditions = {head};
  sim.build();
  const real_t vTotal = sim.surfaceVolume() + sim.gwVolume();
  bool held = false;
  bool released = false;
  real_t vs = sim.surfaceVolume();
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  for (int n = 0; n < 200 && !released; ++n) {
    sim.step();
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    if (c.accumVolume > 0.0 && sim.surfaceVolume() == 0.0) {
      held = true;  // seepage in transit, surface still dry
    }
    if (sim.surfaceVolume() > 0.0) {
      released = true;
      EXPECT_GT(vs + accum, 0.0);
      EXPECT_GT(sim.surfaceVolume(), sim.cfg.surfaceWater.minDepth);
    }
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  EXPECT_TRUE(held);
  EXPECT_TRUE(released);
}

TEST(CoupledModule, BounceBackReturnsInjectionOnFullColumn) {
  // A configured top flux injects into an almost saturated column under a
  // dry surface (legacy qtop < 0): the top cell has no pore room, so the
  // exchange returns the volume to the surface at once
  // (groundwater.c:844-851) and the budget still closes.
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(1, 5, 0.05, 0.5);
  sim.cfg.initialConditions.groundwater.value.constant = 0.4599;
  sim.cfg.boundaryConditions = {fluxBc(frehg::BcTarget::GroundwaterTop, -1.0e-4, 1, 1)};
  sim.build();
  const real_t vTotal = sim.surfaceVolume() + sim.gwVolume();
  real_t vs = sim.surfaceVolume();
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  real_t bounced = 0.0;
  for (int n = 0; n < 5; ++n) {
    sim.step();
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    bounced += c.gw.cplBounce;
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  EXPECT_GT(bounced, 0.0);
  EXPECT_GT(interior2(sim.surface->eta(), 1, 1), 0.0);
}

TEST(CoupledModule, ReallocationVentDepositsOnWetSurface) {
  // Bottom injection faster than the tight soil can pass upward through a
  // saturated column: the corrector over-saturates the bottom cell, the
  // post-allocation send walk finds no room above and vents onto the wet
  // surface (allocate_send, groundwater.c:1185-1196).
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(1, 5, 0.05, 0.5);
  sim.cfg.soil.types = {testSoil(1.0e-7)};
  sim.cfg.groundwater.specificStorage = 1.0e-4;
  sim.cfg.initialConditions.surface.eta.constant = 0.01;  // wet surface
  sim.cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  sim.cfg.initialConditions.groundwater.value.constant = 0.46;
  sim.cfg.boundaryConditions = {fluxBc(frehg::BcTarget::GroundwaterBottom, 1.0e-4, 1, 1)};
  sim.build();
  const real_t vTotal = sim.surfaceVolume() + sim.gwVolume();
  const real_t eta0 = interior2(sim.surface->eta(), 1, 1);
  real_t vs = sim.surfaceVolume();
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  real_t vented = 0.0;
  for (int n = 0; n < 5; ++n) {
    sim.step();
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    vented += c.gw.cplVent;
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  EXPECT_GT(vented, 0.0);
  EXPECT_GT(interior2(sim.surface->eta(), 1, 1), eta0);
}

TEST(CoupledModule, EvaporationCorrectionKeepsSubsurfaceLossOffTheSurface) {
  // A positive configured top flux (evaporation) on a moist column under a
  // dry surface: the volume leaves through the top face but must not pond
  // (groundwater.c:858-863) — the correction removes it from the seepage
  // accumulator and the audit books it as evaporated.
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(1, 5, 0.05, 0.5);
  sim.cfg.initialConditions.groundwater.value.constant = 0.3;
  sim.cfg.boundaryConditions = {fluxBc(frehg::BcTarget::GroundwaterTop, 5.0e-6, 1, 1)};
  sim.build();
  const real_t vTotal = sim.surfaceVolume() + sim.gwVolume();
  real_t vs = sim.surfaceVolume();
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  for (int n = 0; n < 10; ++n) {
    sim.step();
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    EXPECT_GT(c.gw.cplEvap, 0.0);
    EXPECT_EQ(sim.surfaceVolume(), 0.0);
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  EXPECT_LT(sim.gwVolume(), vg + 1.0e-15);
}

TEST(CoupledModule, SyncModeMarchesOnTheCommonAdaptiveStep) {
  // Sync coupling is the legacy adaptive lockstep (solve.c:37,193 —
  // amendment A10): the marching step starts at time.dt and follows the
  // subsurface controller within [dt_min, dt_max].
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(1, 10, 0.05, 0.5);
  sim.cfg.groundwater.timestep.dtInit = 0.5;
  sim.cfg.groundwater.timestep.dtMin = 0.1;
  sim.cfg.groundwater.timestep.dtMax = 2.0;
  sim.build();
  EXPECT_EQ(sim.coupler->nextDt(), 0.5);  // starts from time.dt
  const real_t first = sim.step();
  EXPECT_EQ(first, 0.5);
  // A quiescent column keeps dq below dq_grow: the common step grows by
  // the legacy factor 1.25 until the dt_max clamp.
  real_t expected = 0.5;
  for (int n = 0; n < 10; ++n) {
    expected = std::min(expected * 1.25, 2.0);
    EXPECT_NEAR(sim.step(), expected, 1.0e-12) << "step " << n;
  }
}

}  // namespace
