/// \file test_transport.cpp
/// \brief Unit tests for the scalar-transport module (plan §8.1
///        test_limiters and §10 P4): the superbee limiter values, square-wave
///        advection without new extrema, closed-domain scalar conservation on
///        both grids, the coupled scalar exchange, the configurable bounds,
///        the scalar boundary conditions, the dispersion-tensor diagonal, and
///        the baroclinic activation.

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Types.hpp"
#include "coupling/Coupler.hpp"
#include "gw/RichardsSolver.hpp"
#include "gw/TerrainMetric.hpp"
#include "swe/SurfaceSolver.hpp"
#include "transport/Limiters.hpp"
#include "transport/ScalarSolver.hpp"

#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <memory>
#include <vector>

namespace {

using frehg::real_t;

// ---------------------------------------------------------------------------
// Limiter values (plan §8.1 test_limiters).
// ---------------------------------------------------------------------------

TEST(Limiters, SuperbeePhiMatchesAnalyticValues) {
  // φ(r) = max(0, min(2r, 1), min(r, 2)) at the plan §8.1 sample ratios.
  const real_t r[6] = {-1.0, 0.0, 0.5, 1.0, 2.0, 4.0};
  const real_t expected[6] = {0.0, 0.0, 1.0, 1.0, 2.0, 2.0};
  for (int n = 0; n < 6; ++n) {
    EXPECT_DOUBLE_EQ(frehg::transport::superbeePhi(r[n]), expected[n]) << "r = " << r[n];
  }
}

TEST(Limiters, TvdSuperbeeFaceValue) {
  // tvd_superbee (subroutines.c:541-561): sc + φ(r)/2 (1 - |u| dt/δ)(sp - sc)
  // with r = (sc - sm)/(sp - sc).
  const real_t u = 0.5;
  const real_t delta = 1.0;
  const real_t dt = 0.4;  // coef = 0.2
  // Smooth data: sp = 2, sc = 1, sm = 0 → r = 1, φ = 1.
  EXPECT_NEAR(frehg::transport::tvdSuperbee(2.0, 1.0, 0.0, u, delta, dt),
              1.0 + 0.5 * 1.0 * 0.8 * 1.0, 1.0e-15);
  // Extremum: sp = 2, sc = 1, sm = 2 → r = -1, φ = 0 → donor value.
  EXPECT_DOUBLE_EQ(frehg::transport::tvdSuperbee(2.0, 1.0, 2.0, u, delta, dt), 1.0);
  // Flat pair: sp = sc → donor value regardless of sm.
  EXPECT_DOUBLE_EQ(frehg::transport::tvdSuperbee(1.0, 1.0, 0.3, u, delta, dt), 1.0);
  // Steep upwind gradient: sp = 1.1, sc = 1, sm = 0 → r = 10, φ = 2.
  EXPECT_NEAR(frehg::transport::tvdSuperbee(1.1, 1.0, 0.0, u, delta, dt),
              1.0 + 0.5 * 2.0 * 0.8 * 0.1, 1.0e-15);
}

// ---------------------------------------------------------------------------
// Fixtures.
// ---------------------------------------------------------------------------

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

void setInterior3(const frehg::Field3<real_t>& field, int j, int i, int k, real_t value) {
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  host(static_cast<std::size_t>(j), static_cast<std::size_t>(i), static_cast<std::size_t>(k)) =
      value;
  Kokkos::deep_copy(field, host);
}

/// A single-rank miniature run with the transport module, assembled in the
/// driver's construction and stepping order.
struct MiniTransport {
  frehg::FrehgConfig cfg;
  std::unique_ptr<frehg::Grid> grid;
  std::unique_ptr<frehg::HaloExchanger> halo;
  std::unique_ptr<frehg::gw::TerrainMetric> mesh;
  std::unique_ptr<frehg::BoundarySet> boundaries;
  std::unique_ptr<frehg::swe::SurfaceSolver> surface;
  std::unique_ptr<frehg::gw::RichardsSolver> gw;
  std::unique_ptr<frehg::coupling::Coupler> coupler;
  std::unique_ptr<frehg::transport::ScalarSolver> transport;
  real_t t = 0.0;

  void build() {
    grid = std::make_unique<frehg::Grid>(MPI_COMM_WORLD, cfg.domain);
    halo = std::make_unique<frehg::HaloExchanger>(*grid, false);
    if (cfg.modules.groundwater) {
      mesh = std::make_unique<frehg::gw::TerrainMetric>(*grid, cfg, *halo);
      grid->buildGlobalIds(mesh->ktop());
    } else {
      frehg::HostField2<int> ktop("test_ktop", static_cast<std::size_t>(grid->nyLocal()),
                                  static_cast<std::size_t>(grid->nxLocal()));
      Kokkos::deep_copy(ktop, 0);
      grid->buildGlobalIds(ktop);
    }
    boundaries =
        std::make_unique<frehg::BoundarySet>(*grid, cfg.boundaryConditions, cfg.configDir);
    if (cfg.modules.surfaceWater) {
      surface = std::make_unique<frehg::swe::SurfaceSolver>(*grid, cfg, *boundaries, *halo);
    }
    if (cfg.modules.groundwater) {
      gw = std::make_unique<frehg::gw::RichardsSolver>(*grid, cfg, *mesh, *boundaries, *halo);
    }
    if (surface && gw) {
      coupler = std::make_unique<frehg::coupling::Coupler>(*grid, cfg, *surface, *gw);
    }
    frehg::transport::SurfaceWiring sw;
    if (surface) {
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
    }
    frehg::transport::SubsurfaceWiring gwW;
    if (gw) {
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
    }
    frehg::transport::CouplingWiring cw;
    if (coupler) {
      cw.active = true;
      cw.qss = coupler->seepageRate();
    }
    transport = std::make_unique<frehg::transport::ScalarSolver>(*grid, cfg, *boundaries, *halo,
                                                                 sw, gwW, cw);
    if (gw) {
      gw->attachScalar(transport->subsurfaceScalar(),
                       surface ? transport->surfaceScalar() : frehg::Field2<real_t>(),
                       cfg.groundwater.densityCoupling.enabled);
    }
  }

  /// One flow + transport step in the driver's order; returns the dt taken.
  real_t step() {
    real_t dt = cfg.time.dt;
    if (coupler) {
      dt = coupler->nextDt();
      t += dt;
      coupler->step(t, dt);
      transport->step(t, dt, gw->lastDtg(), surface->currentRain(),
                      surface->currentEvaporation());
    } else if (surface) {
      t += dt;
      surface->beginStep(t);
      surface->solveFreeSurface();
      surface->updateVelocity();
      transport->step(t, dt, 0.0, surface->currentRain(), surface->currentEvaporation());
    } else {
      t += dt;
      gw->step(t, dt);
      transport->step(t, dt, gw->lastDtg(), 0.0, 0.0);
    }
    return dt;
  }
};

/// Surface-only channel: nx cells of deep frictionless water with a uniform
/// initial velocity — the 1D advection stage for the square-wave tests.
frehg::FrehgConfig channelConfig(int nx, real_t depth, real_t u0, real_t dt) {
  frehg::FrehgConfig cfg;
  cfg.simulation.id = "transport-channel";
  cfg.domain.nx = nx;
  cfg.domain.ny = 1;
  cfg.domain.nz = 1;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 1.0;
  cfg.domain.dz = 1.0;
  cfg.domain.bottomElevation.constant = -depth;
  cfg.time.dt = dt;
  cfg.time.tEnd = 1000.0 * dt;
  cfg.time.outputInterval = 1000.0 * dt;
  cfg.modules.surfaceWater = true;
  cfg.modules.transport = true;
  cfg.surfaceWater.gravity = 9.81;
  cfg.surfaceWater.friction.coefficient.constant = 0.0;
  cfg.surfaceWater.minDepth = 1.0e-6;
  cfg.surfaceWater.wettingFaceDepth = 1.0e-6;
  cfg.surfaceWater.friction.thinLayerDepth = 0.01;
  cfg.initialConditions.surface.eta.constant = 0.0;
  cfg.initialConditions.surface.hasUu = true;
  cfg.initialConditions.surface.uu.constant = u0;
  cfg.transport.scheme.advection = frehg::TransportSchemeConfig::Advection::Superbee;
  cfg.transport.surfaceDiffusivityX = 0.0;
  cfg.transport.surfaceDiffusivityY = 0.0;
  cfg.initialConditions.transport.surface.constant = 0.0;
  return cfg;
}

/// Set a square wave s = 1 on i in [lo, hi) of the (single-row) surface
/// scalar, then refresh the ghost/halo state through the module.
void setSquareWave(MiniTransport& mini, int lo, int hi) {
  const frehg::Field2<real_t>& s = mini.transport->surfaceScalar();
  auto host = Kokkos::create_mirror_view(s);
  Kokkos::deep_copy(host, s);
  for (int i = lo; i < hi; ++i) {
    host(1, static_cast<std::size_t>(i) + 1) = 1.0;
  }
  Kokkos::deep_copy(s, host);
  mini.transport->refreshDerivedState(0.0);
}

real_t surfaceScalarMass(const MiniTransport& mini) {
  return mini.transport->ownedSurfaceMass();
}

// ---------------------------------------------------------------------------
// Square-wave advection (plan §10 P4 exit): no new extrema beyond 1e-12, and
// the monotonicity clipping preserves the local extrema of the stencil.
// ---------------------------------------------------------------------------

TEST(TransportModule, SquareWaveAdvectionAddsNoNewExtrema) {
  MiniTransport mini;
  mini.cfg = channelConfig(60, 10.0, 1.0, 0.3);  // CFL = 0.3
  mini.build();
  setSquareWave(mini, 10, 20);

  real_t com0 = 0.0;
  real_t mass0 = 0.0;
  for (int i = 0; i < 60; ++i) {
    const real_t v = interior2(mini.transport->surfaceScalar(), 1, i + 1);
    com0 += v * static_cast<real_t>(i);
    mass0 += v;
  }
  com0 /= mass0;

  for (int n = 0; n < 40; ++n) {
    mini.step();
    for (int i = 0; i < 60; ++i) {
      const real_t v = interior2(mini.transport->surfaceScalar(), 1, i + 1);
      ASSERT_GE(v, -1.0e-12) << "step " << n << " cell " << i;
      ASSERT_LE(v, 1.0 + 1.0e-12) << "step " << n << " cell " << i;
    }
  }

  real_t com1 = 0.0;
  real_t mass1 = 0.0;
  for (int i = 0; i < 60; ++i) {
    const real_t v = interior2(mini.transport->surfaceScalar(), 1, i + 1);
    com1 += v * static_cast<real_t>(i);
    mass1 += v;
  }
  com1 /= mass1;
  // The wave actually advected (the closed basin decelerates the initial
  // uniform flow, so the displacement is a few cells, not u t / dx).
  EXPECT_GT(com1 - com0, 1.5);
}

TEST(TransportModule, UpwindSquareWaveStaysMonotone) {
  MiniTransport mini;
  mini.cfg = channelConfig(60, 10.0, 1.0, 0.3);
  mini.cfg.transport.scheme.advection = frehg::TransportSchemeConfig::Advection::Upwind;
  mini.build();
  setSquareWave(mini, 10, 20);
  for (int n = 0; n < 40; ++n) {
    mini.step();
    for (int i = 0; i < 60; ++i) {
      const real_t v = interior2(mini.transport->surfaceScalar(), 1, i + 1);
      ASSERT_GE(v, -1.0e-12);
      ASSERT_LE(v, 1.0 + 1.0e-12);
    }
  }
}

// ---------------------------------------------------------------------------
// Closed-domain scalar conservation (plan §10 P4 exit: <= 1e-8 per step).
// ---------------------------------------------------------------------------

TEST(TransportModule, StillWaterDiffusionConservesSurfaceScalarMass) {
  // Quiescent closed basin, pure lateral diffusion: the plan §10 P4
  // closed-domain conservation criterion at full strictness (no flow, so
  // the legacy flux-volume lag is inert and Σ s dept A is exact).
  MiniTransport mini;
  mini.cfg = channelConfig(40, 10.0, 0.0, 0.2);
  mini.cfg.transport.surfaceDiffusivityX = 0.01;
  mini.build();
  setSquareWave(mini, 8, 16);

  real_t before = surfaceScalarMass(mini);
  for (int n = 0; n < 30; ++n) {
    mini.step();
    const real_t after = surfaceScalarMass(mini);
    EXPECT_NEAR(after, before, 1.0e-8 * before) << "step " << n;
    before = after;
  }
  // The blob actually spread.
  EXPECT_GT(interior2(mini.transport->surfaceScalar(), 1, 7 + 1), 1.0e-4);
}

TEST(TransportModule, SloshingBasinClosesTheScalarBudget) {
  // Accelerating flow: the legacy flux volume lags the velocity update by
  // one step (volume_by_flux before update_velocity, solve.c/shallowwater.c
  // order), so Σ s dept A drifts by the measured anchor term — the closure
  // identity with every audited term must still hold to rounding.
  MiniTransport mini;
  mini.cfg = channelConfig(40, 10.0, 0.5, 0.2);
  mini.build();
  setSquareWave(mini, 8, 16);

  real_t before = surfaceScalarMass(mini);
  for (int n = 0; n < 30; ++n) {
    mini.step();
    const real_t after = surfaceScalarMass(mini);
    const frehg::transport::TransportAudit& a = mini.transport->audit();
    const real_t residual = (after - before) - (a.exchange + a.surfSource + a.surfBoundary +
                                                a.surfAdjust + a.surfAnchor);
    EXPECT_NEAR(residual, 0.0, 1.0e-8 * before) << "step " << n;
    before = after;
  }
}

TEST(TransportModule, DrainingColumnClosesTheScalarBudget) {
  // A draining unsaturated column, closed on every boundary: the scalar
  // rides the vertical Darcy fluxes. The limiter clips at the moving front
  // and the ledger re-anchors against the post-reallocation θ — both
  // measured, so the closure identity holds to rounding per step.
  MiniTransport mini;
  frehg::FrehgConfig& cfg = mini.cfg;
  cfg.simulation.id = "transport-column";
  cfg.domain.nx = 1;
  cfg.domain.ny = 1;
  cfg.domain.nz = 24;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 1.0;
  cfg.domain.dz = 0.05;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = 2.0;
  cfg.time.tEnd = 1.0e4;
  cfg.time.outputInterval = 1.0e4;
  cfg.modules.groundwater = true;
  cfg.modules.transport = true;
  cfg.groundwater.timestep.dtInit = 2.0;
  cfg.groundwater.timestep.dtMin = 2.0;
  cfg.groundwater.timestep.dtMax = 2.0;
  cfg.groundwater.specificStorage = 0.0;
  cfg.soil.types = {testSoil(1.0e-5)};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  cfg.initialConditions.groundwater.value.constant = 0.3;
  cfg.modules.transport = true;
  cfg.transport.scheme.advection = frehg::TransportSchemeConfig::Advection::Superbee;
  cfg.transport.dispersionLongitudinal = 0.01;
  cfg.transport.dispersionTransverse = 0.002;
  cfg.transport.dispersionMolecular = 1.0e-9;
  cfg.initialConditions.transport.groundwater.constant = 0.0;
  mini.build();

  // A smooth scalar bump mid-column.
  for (int k = 8; k < 16; ++k) {
    const real_t z = (static_cast<real_t>(k) - 11.5) / 4.0;
    setInterior3(mini.transport->subsurfaceScalar(), 1, 1, k, std::exp(-4.0 * z * z));
  }
  mini.transport->refreshDerivedState(0.0);

  real_t before = mini.transport->ownedSubsurfaceMass();
  ASSERT_GT(before, 0.0);
  for (int n = 0; n < 40; ++n) {
    mini.step();
    const real_t after = mini.transport->ownedSubsurfaceMass();
    const frehg::transport::TransportAudit& a = mini.transport->audit();
    const real_t residual =
        (after - before) - (a.subsBoundary + a.subsAdjust + a.subsAnchor - a.exchange);
    EXPECT_NEAR(residual, 0.0, 1.0e-8 * before) << "step " << n;
    before = after;
  }
}

TEST(TransportModule, WaterTableColumnConservesScalarMass) {
  // Hydrostatic-equilibrium column (water-table initial state): the flow is
  // at rest, the dispersive fluxes move the blob, and Σ s θ V holds the
  // plan §10 P4 closed-domain criterion at full strictness.
  MiniTransport mini;
  frehg::FrehgConfig& cfg = mini.cfg;
  cfg.simulation.id = "transport-wt-column";
  cfg.domain.nx = 1;
  cfg.domain.ny = 1;
  cfg.domain.nz = 24;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 1.0;
  cfg.domain.dz = 0.05;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = 2.0;
  cfg.time.tEnd = 1.0e4;
  cfg.time.outputInterval = 1.0e4;
  cfg.modules.groundwater = true;
  cfg.modules.transport = true;
  cfg.groundwater.timestep.dtInit = 2.0;
  cfg.groundwater.timestep.dtMin = 2.0;
  cfg.groundwater.timestep.dtMax = 2.0;
  cfg.groundwater.specificStorage = 0.0;
  cfg.soil.types = {testSoil(1.0e-5)};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::WaterTable;
  cfg.initialConditions.groundwater.value.constant = -0.6;  // mid-box table
  cfg.transport.scheme.advection = frehg::TransportSchemeConfig::Advection::Superbee;
  cfg.transport.dispersionLongitudinal = 0.01;
  cfg.transport.dispersionTransverse = 0.002;
  cfg.transport.dispersionMolecular = 1.0e-7;
  cfg.initialConditions.transport.groundwater.constant = 0.0;
  mini.build();

  // A smooth scalar bump in the saturated zone, clear of the box bottom:
  // the legacy limiter pins a single-conductive-neighbor cell (the bottom
  // cell) to its neighbor's value, which is a measured adjustment, not
  // conservation — keep it inert here.
  for (int k = 13; k < 19; ++k) {
    const real_t z = (static_cast<real_t>(k) - 15.5) / 3.0;
    setInterior3(mini.transport->subsurfaceScalar(), 1, 1, k, std::exp(-4.0 * z * z));
  }
  mini.transport->refreshDerivedState(0.0);

  real_t before = mini.transport->ownedSubsurfaceMass();
  ASSERT_GT(before, 0.0);
  for (int n = 0; n < 40; ++n) {
    mini.step();
    const real_t after = mini.transport->ownedSubsurfaceMass();
    EXPECT_NEAR(after, before, 1.0e-8 * before) << "step " << n;
    before = after;
  }
}

TEST(TransportModule, CoupledExchangeMovesScalarConservatively) {
  // Ponded saline water infiltrating a fresh column (the §8.1 coupled toy
  // with a scalar): the infiltration carries the surface concentration into
  // the subsurface; the combined scalar mass closes per step.
  MiniTransport mini;
  frehg::FrehgConfig& cfg = mini.cfg;
  cfg.simulation.id = "transport-coupled";
  cfg.domain.nx = 3;
  cfg.domain.ny = 1;
  cfg.domain.nz = 12;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 1.0;
  cfg.domain.dz = 0.1;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = 1.0;
  cfg.time.tEnd = 1000.0;
  cfg.time.outputInterval = 1000.0;
  cfg.modules.surfaceWater = true;
  cfg.modules.groundwater = true;
  cfg.modules.transport = true;
  cfg.coupling.mode = frehg::CouplingConfig::Mode::Sync;
  cfg.surfaceWater.gravity = 9.81;
  cfg.surfaceWater.friction.coefficient.constant = 0.02;
  cfg.surfaceWater.minDepth = 1.0e-6;
  cfg.surfaceWater.wettingFaceDepth = 1.0e-6;
  cfg.surfaceWater.friction.thinLayerDepth = 0.01;
  cfg.groundwater.timestep.dtInit = 1.0;
  cfg.groundwater.timestep.dtMin = 1.0;
  cfg.groundwater.timestep.dtMax = 1.0;
  cfg.groundwater.specificStorage = 0.0;
  cfg.soil.types = {testSoil(1.0e-5)};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.surface.eta.constant = 0.2;  // 0.2 m pond
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  cfg.initialConditions.groundwater.value.constant = 0.2;
  cfg.transport.scheme.advection = frehg::TransportSchemeConfig::Advection::Superbee;
  cfg.transport.dispersionLongitudinal = 0.002;
  cfg.transport.dispersionTransverse = 0.0004;
  cfg.transport.dispersionMolecular = 1.0e-10;
  cfg.initialConditions.transport.surface.constant = 35.0;
  cfg.initialConditions.transport.groundwater.constant = 0.0;
  mini.build();

  real_t before = mini.transport->ownedSurfaceMass() + mini.transport->ownedSubsurfaceMass();
  ASSERT_GT(before, 0.0);
  real_t subsGain = 0.0;
  for (int n = 0; n < 30; ++n) {
    mini.step();
    const real_t after =
        mini.transport->ownedSurfaceMass() + mini.transport->ownedSubsurfaceMass();
    // Combined closure: the exchange term cancels between the two grids;
    // the measured source/adjust/anchor terms account for the rest.
    const frehg::transport::TransportAudit& a = mini.transport->audit();
    const real_t residual = (after - before) -
                            (a.surfSource + a.surfBoundary + a.surfAdjust + a.surfAnchor +
                             a.subsBoundary + a.subsAdjust + a.subsAnchor);
    EXPECT_NEAR(residual, 0.0, 1.0e-8 * before) << "step " << n;
    before = after;
    subsGain = mini.transport->ownedSubsurfaceMass();
  }
  // The exchange really moved scalar into the subsurface.
  EXPECT_GT(subsGain, 1.0e-4);
  EXPECT_LT(mini.transport->audit().exchange, 0.0);
}

// ---------------------------------------------------------------------------
// Configurable bounds (plan §3.2: the legacy hard [0, 200] clamp replaced).
// ---------------------------------------------------------------------------

TEST(TransportModule, ConfiguredUpperBoundClampsTheUpdate) {
  MiniTransport mini;
  mini.cfg = channelConfig(20, 10.0, 0.5, 0.2);
  mini.cfg.transport.hasBoundMax = true;
  mini.cfg.transport.boundMax = 0.5;
  mini.build();
  setSquareWave(mini, 5, 10);  // initial s = 1 exceeds the configured max
  mini.step();
  for (int i = 0; i < 20; ++i) {
    EXPECT_LE(interior2(mini.transport->surfaceScalar(), 1, i + 1), 0.5 + 1.0e-12);
  }
}

// ---------------------------------------------------------------------------
// Scalar boundary conditions.
// ---------------------------------------------------------------------------

TEST(TransportModule, EtaPairedScalarValuePrescribesStageSalinity) {
  MiniTransport mini;
  mini.cfg = channelConfig(20, 10.0, 0.0, 0.2);
  frehg::BoundaryConditionConfig etaBc;
  etaBc.name = "stage";
  etaBc.polygon = {{-0.1, -0.1}, {0.9, -0.1}, {0.9, 1.1}, {-0.1, 1.1}};  // cell i = 0
  etaBc.target = frehg::BcTarget::Surface;
  etaBc.kind = frehg::BcKind::Eta;
  etaBc.value.form = frehg::BcValueConfig::Form::Constant;
  etaBc.value.constant = 0.0;
  frehg::BoundaryConditionConfig sBc = etaBc;
  sBc.name = "stage-salinity";
  sBc.kind = frehg::BcKind::ScalarValue;
  sBc.value.constant = 35.0;
  mini.cfg.boundaryConditions = {etaBc, sBc};
  mini.build();
  mini.step();
  EXPECT_DOUBLE_EQ(interior2(mini.transport->surfaceScalar(), 1, 1), 35.0);
  // The Dirichlet reset is measured as a source, not silently created.
  EXPECT_GT(mini.transport->audit().surfSource, 0.0);
}

TEST(TransportModule, SideScalarValueSetsSubsurfaceGhosts) {
  // Groundwater-only box with a head condition and a scalar Dirichlet on
  // the y+ side: the ghost carries the configured concentration, and the
  // inflow advects it into the domain.
  MiniTransport mini;
  frehg::FrehgConfig& cfg = mini.cfg;
  cfg.simulation.id = "transport-side";
  cfg.domain.nx = 1;
  cfg.domain.ny = 6;
  cfg.domain.nz = 8;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 0.5;
  cfg.domain.dz = 0.1;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = 2.0;
  cfg.time.tEnd = 1.0e4;
  cfg.time.outputInterval = 1.0e4;
  cfg.modules.groundwater = true;
  cfg.modules.transport = true;
  cfg.groundwater.timestep.dtInit = 2.0;
  cfg.groundwater.timestep.dtMin = 2.0;
  cfg.groundwater.timestep.dtMax = 2.0;
  cfg.groundwater.specificStorage = 1.0e-5;
  cfg.soil.types = {testSoil(1.0e-4)};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  cfg.initialConditions.groundwater.value.constant = 0.3;
  cfg.transport.scheme.advection = frehg::TransportSchemeConfig::Advection::Upwind;
  cfg.initialConditions.transport.groundwater.constant = 0.0;

  const real_t yEdge = 6.0 * 0.5;
  frehg::BoundaryConditionConfig head;
  head.name = "side-head";
  head.polygon = {{-0.1, yEdge - 0.5}, {1.1, yEdge - 0.5}, {1.1, yEdge + 0.1},
                  {-0.1, yEdge + 0.1}};
  head.target = frehg::BcTarget::GroundwaterSide;
  head.kind = frehg::BcKind::Head;
  head.value.form = frehg::BcValueConfig::Form::Hydrostatic;
  head.value.hydrostaticEta = 0.8;  // pressurized side: inflow
  frehg::BoundaryConditionConfig salt = head;
  salt.name = "side-salt";
  salt.kind = frehg::BcKind::ScalarValue;
  salt.value.form = frehg::BcValueConfig::Form::Constant;
  salt.value.constant = 10.0;
  cfg.boundaryConditions = {head, salt};
  mini.build();

  real_t massBefore = mini.transport->ownedSubsurfaceMass();
  for (int n = 0; n < 20; ++n) {
    mini.step();
  }
  // Ghost slot carries the Dirichlet value; scalar mass entered the domain.
  EXPECT_DOUBLE_EQ(interior3(mini.transport->subsurfaceScalar(), 6 + 1, 1, 4), 10.0);
  EXPECT_GT(mini.transport->ownedSubsurfaceMass(), massBefore + 1.0e-6);
  EXPECT_GT(mini.transport->audit().subsBoundary, 0.0);
}

// ---------------------------------------------------------------------------
// Dispersion tensor diagonal (through the carried top-cell snapshot).
// ---------------------------------------------------------------------------

TEST(TransportModule, DispersionTensorDiagonalOnVerticalFlux) {
  // A column draining under gravity: the top-cell Dzz snapshot equals
  // molecular θs + α_L |q_z| with the legacy *volumetric* face flux
  // (dispersion_tensor, scalar.c:965-994 — the face-area factor stays in).
  MiniTransport mini;
  frehg::FrehgConfig& cfg = mini.cfg;
  cfg.simulation.id = "transport-tensor";
  cfg.domain.nx = 1;
  cfg.domain.ny = 1;
  cfg.domain.nz = 10;
  cfg.domain.dx = 2.0;
  cfg.domain.dy = 3.0;  // cell area 6 m^2: the quirk is visible
  cfg.domain.dz = 0.1;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = 2.0;
  cfg.time.tEnd = 1.0e4;
  cfg.time.outputInterval = 1.0e4;
  cfg.modules.groundwater = true;
  cfg.modules.transport = true;
  cfg.groundwater.timestep.dtInit = 2.0;
  cfg.groundwater.timestep.dtMin = 2.0;
  cfg.groundwater.timestep.dtMax = 2.0;
  cfg.groundwater.specificStorage = 0.0;
  cfg.soil.types = {testSoil(1.0e-5)};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  cfg.initialConditions.groundwater.value.constant = 0.3;
  cfg.transport.dispersionLongitudinal = 0.01;
  cfg.transport.dispersionTransverse = 0.002;
  cfg.transport.dispersionMolecular = 1.0e-9;
  cfg.initialConditions.transport.groundwater.constant = 0.0;
  mini.build();
  mini.step();

  // Expected from the module's own flux observables: |q| at the top cell's
  // lower face (per-area output times the face area).
  const real_t qzPerArea = interior3(mini.gw->fluxZPerArea(), 1, 1, 0);
  const real_t qzVolumetric = qzPerArea * 6.0;
  const real_t expected = 1.0e-9 * 0.46 + 0.01 * std::fabs(qzVolumetric);
  const real_t dzzTop = interior2(mini.transport->dispersionTopSnapshot(), 1, 1);
  EXPECT_NEAR(dzzTop, expected, 1.0e-12 + 1.0e-9 * expected);
  EXPECT_GT(std::fabs(qzVolumetric), 0.0);
}

// ---------------------------------------------------------------------------
// Baroclinic activation (plan §10 P4: r_rho / r_visc from the scalar).
// ---------------------------------------------------------------------------

TEST(Baroclinic, RatiosFollowTheScalarAndAverageAtFaces) {
  MiniTransport mini;
  frehg::FrehgConfig& cfg = mini.cfg;
  cfg.simulation.id = "transport-baroclinic";
  cfg.domain.nx = 1;
  cfg.domain.ny = 4;
  cfg.domain.nz = 6;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 1.0;
  cfg.domain.dz = 0.1;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = 1.0;
  cfg.time.tEnd = 100.0;
  cfg.time.outputInterval = 100.0;
  cfg.modules.groundwater = true;
  cfg.modules.transport = true;
  cfg.groundwater.timestep.dtInit = 1.0;
  cfg.groundwater.timestep.dtMin = 1.0;
  cfg.groundwater.timestep.dtMax = 1.0;
  cfg.groundwater.specificStorage = 1.0e-5;
  cfg.groundwater.densityCoupling.enabled = true;
  cfg.soil.types = {testSoil(1.0e-5)};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  cfg.initialConditions.groundwater.value.constant = 0.3;
  cfg.transport.scheme.advection = frehg::TransportSchemeConfig::Advection::Upwind;
  cfg.initialConditions.transport.groundwater.constant = 0.0;
  mini.build();

  // A split field: rows j = 0, 1 fresh, rows j = 2, 3 at 35 psu.
  for (int j = 2; j < 4; ++j) {
    for (int k = 0; k < 6; ++k) {
      setInterior3(mini.transport->subsurfaceScalar(), j + 1, 1, k, 35.0);
    }
  }
  mini.transport->refreshDerivedState(0.0);
  mini.step();

  const real_t rFresh = 1.0;
  const real_t rSalt = 1.0 + 35.0 * frehg::gw::kBaroclinicBetaRho;
  EXPECT_NEAR(interior3(mini.gw->densityRatio(), 1, 1, 2), rFresh, 1.0e-12);
  EXPECT_NEAR(interior3(mini.gw->densityRatio(), 3, 1, 2), rSalt, 1.0e-12);
  EXPECT_NEAR(interior3(mini.gw->viscosityRatio(), 3, 1, 2), 1.0 / (1.0 + 35.0 * 0.0022),
              1.0e-12);
  // Interior y face between rows 1 and 2: the arithmetic mean.
  EXPECT_NEAR(interior3(mini.gw->densityRatioFaceY(), 2, 1, 2), 0.5 * (rFresh + rSalt),
              1.0e-12);
  // Interior z face inside the salty rows: the mean of equal values.
  EXPECT_NEAR(interior3(mini.gw->densityRatioFaceZ(), 3, 1, 3), rSalt, 1.0e-12);
}

}  // namespace
