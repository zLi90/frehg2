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
#include "transport/ScalarSolver.hpp"

#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
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
  /// Optional temperature-transport instance (item: the coupled thermal
  /// side-Dirichlet admission); built by buildTemperature() after build().
  std::unique_ptr<frehg::transport::ScalarSolver> temperature;
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

  /// Wire the temperature ScalarSolver exactly as the driver does
  /// (Simulation::buildTransport with modules.temperature).
  void buildTemperature() {
    frehg::transport::SurfaceWiring sw;
    sw.active = true;
    sw.minDepth = surface->minDepth();
    sw.dept = surface->depth();
    sw.etan = surface->etaStart();
    sw.bottom = surface->bottom();
    sw.uu = surface->uu();
    sw.vv = surface->vv();
    sw.fu = surface->flowRateX();
    sw.fv = surface->flowRateY();
    sw.asx = surface->faceAreaX();
    sw.asy = surface->faceAreaY();
    sw.rainMask = surface->rainApplyMask();
    frehg::transport::SubsurfaceWiring gwW;
    gwW.active = true;
    gwW.wc = gw->waterContent();
    gwW.wcn = gw->waterContentStart();
    gwW.qx = gw->fluxXVolumetric();
    gwW.qy = gw->fluxYVolumetric();
    gwW.qzF = gw->fluxZFaceVolumetric();
    gwW.kx = gw->faceConductivityX();
    gwW.ky = gw->faceConductivityY();
    gwW.kzF = gw->faceConductivityZFace();
    gwW.wcs = gw->soilThetaS();
    gwW.ksz = gw->soilKsz();
    gwW.dz3d = mesh->dz3d();
    gwW.ax = mesh->areaX();
    gwW.ay = mesh->areaY();
    gwW.cosx = mesh->cosX();
    gwW.cosy = mesh->cosY();
    gwW.az = mesh->areaZ();
    gwW.topCode = gw->topBcCode();
    gwW.topValue = gw->topBcValue();
    gwW.sideCodeYp = gw->sideBcCodeYp();
    gwW.sideCodeYm = gw->sideBcCodeYm();
    gwW.sideCodeXm = gw->sideBcCodeXm();
    gwW.sideCodeXp = gw->sideBcCodeXp();
    frehg::transport::CouplingWiring cw;
    cw.active = true;
    cw.qss = coupler->seepageRate();
    temperature = std::make_unique<frehg::transport::ScalarSolver>(
        *grid, cfg, *boundaries, *halo, sw, gwW, cw,
        frehg::transport::ScalarSpec::temperature(cfg));
  }

  /// One coupled step at the coupler's own step size; returns the dt taken.
  real_t step() {
    const real_t dt = coupler->nextDt();
    t += dt;
    coupler->step(t, dt);
    if (temperature) {
      temperature->step(t, dt, gw->lastDtg(), surface->currentRain(),
                        surface->currentEvaporation());
    }
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

/// Writes a flat-list data file into a fresh temp directory; returns the
/// dir (the test_gw.cpp pattern).
std::string writeFlatFile(const std::string& name, const std::vector<real_t>& values) {
  const std::filesystem::path dir =
      std::filesystem::temp_directory_path() /
      ("frehg_coupling_test_" + std::to_string(::getpid()) + "_" + name);
  std::filesystem::create_directories(dir);
  std::ofstream out(dir / name);
  out.precision(17);
  for (const real_t v : values) {
    out << v << "\n";
  }
  return dir.string();
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

// ---------------------------------------------------------------------------
// v2 Q7 §8.2 backfill: coupled surface boundary conditions under an active
// exchange, and the coupled branch of the hydrostatic side ghost.
// ---------------------------------------------------------------------------

TEST(CoupledModule, DischargeInflowClosesBudgetUnderExchange) {
  // A ponded coupled column set fed by a surface discharge inflow while the
  // pond infiltrates: the injected volume is measured into bcInflow
  // (Q dt / global member count per cell, FreeSurface.cpp assembleRhs) and
  // the coupled per-step closure identity holds with the exchange active.
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(3, 10, 0.05, 0.5);
  sim.cfg.initialConditions.surface.eta.constant = 0.02;  // pond: exchange on
  frehg::BoundaryConditionConfig inlet;
  inlet.name = "inlet";
  inlet.polygon = {{-0.1, -0.1}, {1.1, -0.1}, {1.1, 1.1}, {-0.1, 1.1}};  // cell i = 0
  inlet.target = frehg::BcTarget::Surface;
  inlet.kind = frehg::BcKind::Discharge;
  inlet.value.form = frehg::BcValueConfig::Form::Constant;
  inlet.value.constant = 2.0e-4;  // total inflow [m^3/s] over the region
  sim.cfg.boundaryConditions = {inlet};
  sim.build();
  const real_t vTotal = sim.surfaceVolume() + sim.gwVolume();
  real_t vs = sim.surfaceVolume();
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  real_t inflowSum = 0.0;
  real_t exchangedSum = 0.0;
  for (int n = 0; n < 40; ++n) {
    sim.step();
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    const frehg::swe::SurfaceStepAudit& s = sim.surface->audit();
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    // The discharge volume is exact: Q dt into the single member cell.
    EXPECT_NEAR(s.bcInflow, 2.0e-4 * 0.5, 1.0e-15) << "step " << n;
    inflowSum += s.bcInflow;
    exchangedSum += c.gw.cplExchanged;
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  EXPECT_NEAR(inflowSum, 2.0e-4 * 0.5 * 40, 1.0e-12);
  // Infiltration really ran alongside the inflow (net top-face flux down).
  EXPECT_LT(exchangedSum, -1.0e-4);
  EXPECT_GT(interior3(sim.gw->waterContent(), 1, 2, 0), 0.2);
  // The injected volume stayed in the coupled system (closed edges).
  EXPECT_GT(sim.surfaceVolume() + sim.gwVolume(), vTotal + 1.0e-3);
}

TEST(CoupledModule, OutflowEdgeDrainsCoupledPondConservatively) {
  // A coupled pond over a bed sloping down to a WEST free-outflow edge
  // (kind outflow continues the bed slope across the boundary face; the
  // west/south faces read the halo-column face geometry — the V2-A11 face
  // fix). The budget closes per step, the pond genuinely drains through
  // the edge, and the below-bed clamp mints nothing.
  const std::string dir =
      writeFlatFile("bath.dat", {0.0, 0.05, 0.10, 0.0, 0.05, 0.10, 0.0, 0.05, 0.10});
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(3, 6, 0.05, 0.5);
  sim.cfg.domain.ny = 3;
  sim.cfg.configDir = dir;
  sim.cfg.domain.bottomElevation.fromFile = true;
  sim.cfg.domain.bottomElevation.file = "bath.dat";
  sim.cfg.initialConditions.surface.eta.constant = 0.15;  // 5 cm over the crest
  // The closure identity is asserted at 1e-10 of the coupled volume; the
  // outflow correction is an assembly-side booking, so the eta solve must
  // converge past the assertion floor (the default rtol 1e-8 leaves
  // solver-residual-sized closure defects in this actively draining toy,
  // unlike the quiescent ponds above).
  sim.cfg.solver.surface.rtol = 1.0e-13;
  frehg::BoundaryConditionConfig out;
  out.name = "outlet";
  out.polygon = {{-0.1, -0.1}, {0.9, -0.1}, {0.9, 3.1}, {-0.1, 3.1}};  // west column
  out.target = frehg::BcTarget::Surface;
  out.kind = frehg::BcKind::Outflow;
  sim.cfg.boundaryConditions = {out};
  sim.build();
  const real_t vSurf0 = sim.surfaceVolume();
  const real_t vTotal = vSurf0 + sim.gwVolume();
  real_t vs = vSurf0;
  real_t vg = sim.gwVolume();
  real_t accum = 0.0;
  real_t outflowSum = 0.0;
  real_t clampSum = 0.0;
  // 40 steps drain most of the pond while every cell stays wet; the final
  // dry-out (around step 50 at this slope) carries a ~4e-10 absolute
  // closure spike at the wet/dry transition steps, which is drying-front
  // bookkeeping, not the outflow condition under test.
  for (int n = 0; n < 40; ++n) {
    sim.step();
    const frehg::coupling::CouplingAudit& c = sim.coupler->audit();
    const frehg::swe::SurfaceStepAudit& s = sim.surface->audit();
    EXPECT_LE(std::fabs(sim.closureResidual(vs, vg, accum)), 1.0e-10 * vTotal) << "step " << n;
    outflowSum += s.boundaryOutflow;
    clampSum += std::fabs(s.clampVolume);
    vs = sim.surfaceVolume();
    vg = sim.gwVolume();
    accum = c.accumVolume;
  }
  // The edge really released the pond (measured 0.851 of the 0.9 m^3)...
  EXPECT_GT(outflowSum, 0.5);
  // ...the pond drained (measured 0.025 m^3 remaining)...
  EXPECT_LT(sim.surfaceVolume(), 0.1 * vSurf0);
  // ...and the below-bed clamp minted nothing.
  EXPECT_EQ(clampSum, 0.0);
}

TEST(CoupledModule, HydrostaticSideHeadFollowsLiveSurfaceWhenCoupled) {
  // The coupled branch of ghostHead (gw/Predictor.cpp): a hydrostatic side
  // stage ABOVE the edge bed is a ponded/tidal boundary, and in a coupled
  // run its ghost follows the LIVE local surface — (bed - zc) + dept —
  // rather than the configured stage (the uncoupled fixed form
  // (value - zc), pinned by GwModule.HydrostaticSideHoldsEquilibrium).
  // Bed 0, nz = 5 x dz = 0.1; west-middle cell; mid-depth cell k = 2 has
  // zc = -0.25.
  const auto makeSim = [](real_t eta0, real_t stage) {
    MiniCoupled sim;
    sim.cfg = coupledBaseConfig(3, 5, 0.1, 0.5);
    sim.cfg.domain.ny = 3;
    sim.cfg.initialConditions.surface.eta.constant = eta0;
    frehg::BoundaryConditionConfig side;
    side.name = "sea";
    side.polygon = {{-0.1, 1.1}, {0.9, 1.1}, {0.9, 1.9}, {-0.1, 1.9}};  // west mid cell
    side.target = frehg::BcTarget::GroundwaterSide;
    side.kind = frehg::BcKind::Head;
    side.value.form = frehg::BcValueConfig::Form::Hydrostatic;
    side.value.hydrostaticEta = stage;  // above the bed: the coupled branch
    sim.cfg.boundaryConditions = {side};
    sim.build();
    // Re-enforce the ghost heads now that the coupler is attached (the gw
    // constructor ran before attachCoupling); no physics has run, so dept
    // is exactly the initial pond depth.
    sim.gw->refreshDerivedState();
    return sim;
  };
  const real_t zc = -0.25;

  MiniCoupled shallow = makeSim(0.05, 0.5);
  const real_t ghostShallow = interior3(shallow.gw->head(), 2, 0, 2);
  // The ghost is the live-surface hydrostatic column, not the fixed stage.
  EXPECT_NEAR(ghostShallow, (0.0 - zc) + 0.05, 1.0e-12);
  EXPECT_GT(std::fabs(ghostShallow - (0.5 - zc)), 0.1);

  // A deeper pond moves the ghost by exactly the depth change...
  MiniCoupled deep = makeSim(0.15, 0.5);
  const real_t ghostDeep = interior3(deep.gw->head(), 2, 0, 2);
  EXPECT_NEAR(ghostDeep - ghostShallow, 0.10, 1.0e-12);

  // ...while the configured stage value is inert once above the bed.
  MiniCoupled otherStage = makeSim(0.05, 0.9);
  EXPECT_DOUBLE_EQ(interior3(otherStage.gw->head(), 2, 0, 2), ghostShallow);

  // And the side flux follows the live surface the way the formula says:
  // the deeper pond's larger ghost head drives more inflow through the
  // west boundary face (positive q points toward -x, so inflow is
  // negative at the west face slot i = 0).
  shallow.step();
  deep.step();
  const real_t qxShallow = interior3(shallow.gw->fluxXVolumetric(), 2, 0, 2);
  const real_t qxDeep = interior3(deep.gw->fluxXVolumetric(), 2, 0, 2);
  EXPECT_LT(qxShallow, 0.0);
  EXPECT_LT(qxDeep, qxShallow);
}

TEST(CoupledModule, SideThermalDirichletAdmittedUnderCoupling) {
  // The coupled twin of the Q5 thermal side-admission fix (V2-A17): a
  // temperature side scalar_value on a NON-y+ side (west here) is
  // limiter-admitted in a COUPLED run — the thermal extrema branches admit
  // prescribed side ghosts on ALL FOUR sides, where the salinity spec
  // admits y+ only (TransportModule.SalinitySideValueIsLimiterClippedOffYPlus
  // pins that asymmetry uncoupled). Until now the thermal admission was
  // exercised only uncoupled (the heat 8-orientation battery).
  MiniCoupled sim;
  sim.cfg = coupledBaseConfig(4, 6, 0.1, 2.0);
  sim.cfg.domain.ny = 3;
  sim.cfg.groundwater.specificStorage = 1.0e-5;
  sim.cfg.soil.types = {testSoil(1.0e-4)};
  sim.cfg.initialConditions.surface.eta.constant = 0.0;  // dry surface
  sim.cfg.initialConditions.groundwater.value.constant = 0.42;  // conductive
  sim.cfg.modules.temperature = true;
  sim.cfg.temperature.thermalConductivity = 1.0;
  sim.cfg.temperature.heatCapacitySolid = 2.2e6;
  sim.cfg.initialConditions.temperature.groundwater.constant = 10.0;
  sim.cfg.initialConditions.temperature.surface.constant = 0.0;
  frehg::BoundaryConditionConfig head;
  head.name = "side-head";
  head.polygon = {{-0.1, 1.1}, {0.9, 1.1}, {0.9, 1.9}, {-0.1, 1.9}};  // west mid cell
  head.target = frehg::BcTarget::GroundwaterSide;
  head.kind = frehg::BcKind::Head;
  head.value.form = frehg::BcValueConfig::Form::Constant;
  head.value.constant = 0.5;  // pressurized side: inflow
  frehg::BoundaryConditionConfig heat = head;
  heat.name = "side-heat";
  heat.kind = frehg::BcKind::ScalarValue;
  heat.scalarField = frehg::BcScalar::Temperature;
  heat.value.constant = 30.0;
  sim.cfg.boundaryConditions = {head, heat};
  sim.build();
  sim.buildTemperature();

  const real_t heatBefore = sim.temperature->ownedSubsurfaceMass();
  for (int n = 0; n < 40; ++n) {
    sim.step();
  }
  // The prescribed ghost is written on the west side...
  EXPECT_DOUBLE_EQ(interior3(sim.temperature->subsurfaceScalar(), 2, 0, 3), 30.0);
  // ...and — unlike the salinity y+-only rule — the limiter ADMITS it: the
  // edge cell warms measurably above the wet-stencil maximum (10, where a
  // clipped update would be pinned exactly; measured 10.127), and the
  // boundary heat enters the ledger (measured 0.081 K m^3 over 40 steps).
  EXPECT_GT(interior3(sim.temperature->subsurfaceScalar(), 2, 1, 3), 10.05);
  EXPECT_GT(sim.temperature->ownedSubsurfaceMass(), heatBefore + 0.02);
  EXPECT_GT(sim.temperature->audit().subsBoundary, 0.0);
}

}  // namespace
