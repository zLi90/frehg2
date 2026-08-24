/// \file test_swe.cpp
/// \brief Unit tests for the surface-water module (plan §8.1): closure
///        formulas against hand-computed values, and miniature
///        configurations that pin the source terms, boundary-condition
///        kinds, wet/dry behavior, and volume bookkeeping.

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Types.hpp"
#include "swe/SurfaceSolver.hpp"
#include "swe/SweFormulas.hpp"

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

real_t interiorValue(const frehg::Field2<real_t>& field, int j, int i) {
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  return host(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
}

/// A single-rank miniature surface-water run assembled from an in-code
/// configuration.
struct MiniSim {
  frehg::FrehgConfig cfg;
  std::unique_ptr<frehg::Grid> grid;
  std::unique_ptr<frehg::HaloExchanger> halo;
  std::unique_ptr<frehg::BoundarySet> boundaries;
  std::unique_ptr<frehg::swe::SurfaceSolver> swe;

  void build() {
    grid = std::make_unique<frehg::Grid>(MPI_COMM_WORLD, cfg.domain);
    frehg::HostField2<int> ktop("test_ktop", static_cast<std::size_t>(grid->nyLocal()),
                                static_cast<std::size_t>(grid->nxLocal()));
    Kokkos::deep_copy(ktop, 0);
    grid->buildGlobalIds(ktop);
    halo = std::make_unique<frehg::HaloExchanger>(*grid, false);
    boundaries =
        std::make_unique<frehg::BoundarySet>(*grid, cfg.boundaryConditions, cfg.configDir);
    swe = std::make_unique<frehg::swe::SurfaceSolver>(*grid, cfg, *boundaries, *halo);
  }

  void step(real_t t) {
    swe->beginStep(t);
    swe->solveFreeSurface();
    swe->updateVelocity();
  }
};

frehg::FrehgConfig baseConfig(int nx, int ny, real_t dx, real_t dt) {
  frehg::FrehgConfig cfg;
  cfg.simulation.id = "swe-unit";
  cfg.domain.nx = nx;
  cfg.domain.ny = ny;
  cfg.domain.nz = 1;
  cfg.domain.dx = dx;
  cfg.domain.dy = dx;
  cfg.domain.dz = 1.0;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = dt;
  cfg.time.tEnd = 1000.0 * dt;
  cfg.time.outputInterval = 1000.0 * dt;
  cfg.modules.surfaceWater = true;
  cfg.surfaceWater.friction.coefficient.constant = 0.01;
  cfg.surfaceWater.minDepth = 1.0e-8;
  cfg.surfaceWater.wettingFaceDepth = 1.0e-8;
  cfg.initialConditions.surface.eta.constant = 0.5;
  return cfg;
}

/// Writes a flat-list raster into a fresh temp directory; returns the dir.
std::string writeRaster(const std::string& name, const std::vector<real_t>& values) {
  const std::filesystem::path dir =
      std::filesystem::temp_directory_path() /
      ("frehg_swe_test_" + std::to_string(::getpid()) + "_" + name);
  std::filesystem::create_directories(dir);
  std::ofstream out(dir / name);
  out.precision(17);
  for (const real_t v : values) {
    out << v << "\n";
  }
  return dir.string();
}

// ---------------------------------------------------------------------------
// Closure formulas vs hand-computed values (legacy provenance in
// SweFormulas.hpp).
// ---------------------------------------------------------------------------

TEST(SweFormulas, AdvectionRampMatchesLegacyBand) {
  // Below the band: untouched. Inside: linear. Above or at rest: zero.
  EXPECT_DOUBLE_EQ(frehg::swe::cflDampedAdvection(2.0, 1.0, 0.4), 2.0);
  EXPECT_DOUBLE_EQ(frehg::swe::cflDampedAdvection(2.0, 1.0, 0.6), 2.0 * (0.7 - 0.6) / 0.2);
  EXPECT_DOUBLE_EQ(frehg::swe::cflDampedAdvection(2.0, 1.0, 0.65), 2.0 * (0.7 - 0.65) / 0.2);
  EXPECT_DOUBLE_EQ(frehg::swe::cflDampedAdvection(2.0, 1.0, 0.71), 0.0);
  EXPECT_DOUBLE_EQ(frehg::swe::cflDampedAdvection(2.0, 0.0, 0.1), 0.0);
}

TEST(SweFormulas, ManningDragSwitchesExponentAtThinLayerDepth) {
  const real_t g = 9.81;
  const real_t n = 0.019;
  // Below hD the exponent is 2/3, above it 1/3 (shallowwater.c:988-989).
  EXPECT_DOUBLE_EQ(frehg::swe::manningDrag(g, n, 0.05, 0.1),
                   g * n * n / std::pow(0.05, 2.0 / 3.0));
  EXPECT_DOUBLE_EQ(frehg::swe::manningDrag(g, n, 0.2, 0.1),
                   g * n * n / std::pow(0.2, 1.0 / 3.0));
}

TEST(SweFormulas, ChezyDragHasNoDepthDependence) {
  EXPECT_DOUBLE_EQ(frehg::swe::chezyDrag(9.81, 1.767), 9.81 / (1.767 * 1.767));
}

TEST(SweFormulas, PointImplicitFactor) {
  EXPECT_DOUBLE_EQ(frehg::swe::pointImplicitFactor(5.0, 0.03, 0.4, 2.0),
                   1.0 / (0.5 * 5.0 * 0.03 * 0.4 * 2.0 + 1.0));
  EXPECT_DOUBLE_EQ(frehg::swe::pointImplicitFactor(5.0, 0.03, 0.0, 2.0), 1.0);
}

TEST(SweFormulas, WindStressAndThinLayerAttenuation) {
  const real_t omega = 0.5;
  const real_t tau = frehg::swe::windStress(0.0013, 10.0, 0.2, 0.1, omega);
  const real_t relative = 10.0 - 0.2 * std::cos(omega) - 0.1 * std::sin(omega);
  EXPECT_DOUBLE_EQ(tau, frehg::swe::kAirDensity * 0.0013 * relative * relative);
  // Full stress at or above hD; exponential attenuation below; zero below
  // hD/2 (shallowwater.c:271-286).
  EXPECT_DOUBLE_EQ(frehg::swe::windStressAttenuated(tau, 0.2, 0.1, 5.0), tau);
  EXPECT_DOUBLE_EQ(frehg::swe::windStressAttenuated(tau, 0.08, 0.1, 5.0),
                   tau * std::exp(5.0 * (0.08 - 0.1) / 0.1));
  EXPECT_DOUBLE_EQ(frehg::swe::windStressAttenuated(tau, 0.04, 0.1, 5.0), 0.0);
}

// ---------------------------------------------------------------------------
// Sources on miniature domains.
// ---------------------------------------------------------------------------

TEST(SweModule, RainAccumulatesExactlyOnAClosedCell) {
  // A 1x1 domain has zero face areas in both directions (the legacy nx==1 /
  // ny==1 handling, shallowwater.c:1056-1057), so it is genuinely closed.
  // Larger flat basins are NOT still under rain: legacy refreshes the
  // physical-edge ghost stage only at the start of solve_shallowwater, so
  // update_velocity sees the pre-rain ghost and rain induces a small
  // outward boundary-face velocity — a preserved call-order artifact
  // (docs/theory/surface-water.md).
  MiniSim sim;
  sim.cfg = baseConfig(1, 1, 1.0, 1.0);
  sim.cfg.surfaceWater.rainfall.constant = 1.0e-4;
  sim.build();
  const real_t v0 = sim.swe->ownedVolume();
  real_t rainTotal = 0.0;
  for (int n = 1; n <= 5; ++n) {
    sim.step(static_cast<real_t>(n));
    rainTotal += sim.swe->audit().rainVolume;
  }
  EXPECT_NEAR(interiorValue(sim.swe->eta(), 1, 1), 0.5 + 5.0e-4, 1.0e-12);
  EXPECT_NEAR(sim.swe->ownedVolume() - v0, rainTotal, 1.0e-9);
  EXPECT_NEAR(rainTotal, 5.0 * 1.0e-4, 1.0e-12);
}

TEST(SweModule, RainExclusionRegionReceivesNothing) {
  MiniSim sim;
  sim.cfg = baseConfig(3, 3, 1.0, 1.0);
  sim.cfg.surfaceWater.rainfall.constant = 1.0e-4;
  // Exclude the center cell (center at (1.5, 1.5)).
  sim.cfg.surfaceWater.rainfallExcludePolygon = {{1.1, 1.1}, {1.9, 1.1}, {1.9, 1.9}, {1.1, 1.9}};
  sim.build();
  sim.step(1.0);
  EXPECT_NEAR(interiorValue(sim.swe->eta(), 2, 2), 0.5, 1.0e-12);        // excluded center
  EXPECT_NEAR(interiorValue(sim.swe->eta(), 1, 1), 0.5 + 1.0e-4, 1.0e-12);
  EXPECT_NEAR(sim.swe->audit().rainVolume, 1.0e-4 * 8.0, 1.0e-12);
}

TEST(SweModule, EvaporationDriesToTheBedAndNoFurther) {
  MiniSim sim;
  sim.cfg = baseConfig(2, 2, 1.0, 1.0);
  sim.cfg.initialConditions.surface.eta.constant = 1.5e-4;  // 0.15 mm of water
  sim.cfg.surfaceWater.evaporation.constant = 1.0e-4;
  sim.build();
  sim.step(1.0);
  EXPECT_NEAR(interiorValue(sim.swe->depth(), 1, 1), 0.5e-4, 1.0e-15);
  sim.step(2.0);
  // The clamp absorbs over-drying (shallowwater.c:628-635).
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->eta(), 1, 1), 0.0);
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->depth(), 1, 1), 0.0);
}

// ---------------------------------------------------------------------------
// Boundary-condition kinds.
// ---------------------------------------------------------------------------

TEST(SweModule, EtaConditionHoldsTheStageAndFillsTheChannel) {
  MiniSim sim;
  sim.cfg = baseConfig(5, 1, 1.0, 0.1);
  sim.cfg.initialConditions.surface.eta.constant = 0.2;
  frehg::BoundaryConditionConfig bc;
  bc.name = "stage";
  bc.polygon = {{-0.1, -0.1}, {1.0, -0.1}, {1.0, 1.1}, {-0.1, 1.1}};  // cell i = 0
  bc.target = frehg::BcTarget::Surface;
  bc.kind = frehg::BcKind::Eta;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = 0.4;
  sim.cfg.boundaryConditions.push_back(bc);
  sim.build();
  for (int n = 1; n <= 400; ++n) {
    sim.step(0.1 * static_cast<real_t>(n));
  }
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->eta(), 1, 1), 0.4);  // held exactly
  // The stage propagates down the channel; the weakly damped seiche and the
  // open east end keep the far cell near — not at — the stage.
  const real_t far = interiorValue(sim.swe->eta(), 1, 5);
  EXPECT_GT(far, 0.35);
  EXPECT_LT(far, 0.47);
}

TEST(SweModule, DischargeConditionInjectsExactVolume) {
  MiniSim sim;
  sim.cfg = baseConfig(3, 3, 1.0, 1.0);
  frehg::BoundaryConditionConfig bc;
  bc.name = "inflow";
  bc.polygon = {{1.1, 1.1}, {1.9, 1.1}, {1.9, 1.9}, {1.1, 1.9}};  // center cell
  bc.target = frehg::BcTarget::Surface;
  bc.kind = frehg::BcKind::Discharge;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = 0.01;  // m^3/s
  sim.cfg.boundaryConditions.push_back(bc);
  sim.build();
  const real_t v0 = sim.swe->ownedVolume();
  real_t leaked = 0.0;
  for (int n = 1; n <= 3; ++n) {
    sim.step(static_cast<real_t>(n));
    EXPECT_NEAR(sim.swe->audit().bcInflow, 0.01, 1.0e-12);
    leaked += sim.swe->audit().boundaryOutflow;
  }
  EXPECT_NEAR(sim.swe->ownedVolume() - v0, 3.0 * 0.01 - leaked, 1.0e-6);
}

TEST(SweModule, VelocityConditionDrivesAndForcesTheEdgeFace) {
  MiniSim sim;
  sim.cfg = baseConfig(5, 1, 1.0, 0.1);
  frehg::BoundaryConditionConfig bc;
  bc.name = "west_inflow";
  bc.polygon = {{-0.1, -0.1}, {1.0, -0.1}, {1.0, 1.1}, {-0.1, 1.1}};  // west edge cell
  bc.target = frehg::BcTarget::Surface;
  bc.kind = frehg::BcKind::Velocity;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = 0.05;  // m/s into the domain through the west face
  sim.cfg.boundaryConditions.push_back(bc);
  sim.build();
  const real_t v0 = sim.swe->ownedVolume();
  real_t injected = 0.0;
  real_t leaked = 0.0;
  for (int n = 1; n <= 10; ++n) {
    sim.step(0.1 * static_cast<real_t>(n));
    injected += sim.swe->audit().bcInflow;
    leaked += sim.swe->audit().boundaryOutflow;
  }
  // The prescribed face velocity is applied to the halo face.
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->uu(), 1, 0), 0.05);
  EXPECT_GT(injected, 0.0);
  EXPECT_NEAR(sim.swe->ownedVolume() - v0, injected - leaked, 1.0e-6);
}

TEST(SweModule, OutflowConditionDrainsASlopedChannelConservatively) {
  // Bed sloping down toward the east edge; free outflow at the east column.
  std::vector<real_t> bed = {0.4, 0.3, 0.2, 0.1, 0.0};
  const std::string dir = writeRaster("bed.dat", bed);
  MiniSim sim;
  sim.cfg = baseConfig(5, 1, 1.0, 0.1);
  sim.cfg.configDir = dir;
  sim.cfg.domain.bottomElevation.fromFile = true;
  sim.cfg.domain.bottomElevation.file = "bed.dat";
  sim.cfg.initialConditions.surface.eta.constant = 0.6;
  frehg::BoundaryConditionConfig bc;
  bc.name = "outlet";
  bc.polygon = {{4.0, -0.1}, {5.1, -0.1}, {5.1, 1.1}, {4.0, 1.1}};  // east edge cell
  bc.target = frehg::BcTarget::Surface;
  bc.kind = frehg::BcKind::Outflow;
  sim.cfg.boundaryConditions.push_back(bc);
  sim.build();
  const real_t v0 = sim.swe->ownedVolume();
  real_t released = 0.0;
  for (int n = 1; n <= 600; ++n) {
    sim.step(0.1 * static_cast<real_t>(n));
    released += sim.swe->audit().boundaryOutflow;
  }
  const real_t vEnd = sim.swe->ownedVolume();
  EXPECT_LT(vEnd, 0.25 * v0);                      // the channel drains
  EXPECT_NEAR(v0 - vEnd, released, 1.0e-4 * v0);   // and the budget closes
}

// ---------------------------------------------------------------------------
// Wet/dry machinery.
// ---------------------------------------------------------------------------

TEST(SweModule, HigherOfTwoBottomsBlocksFlowAcrossADam) {
  std::vector<real_t> bed = {0.0, 1.0, 0.0};
  const std::string dir = writeRaster("dam.dat", bed);
  MiniSim sim;
  sim.cfg = baseConfig(3, 1, 1.0, 0.1);
  sim.cfg.configDir = dir;
  sim.cfg.domain.bottomElevation.fromFile = true;
  sim.cfg.domain.bottomElevation.file = "dam.dat";
  sim.cfg.initialConditions.surface.eta.constant = 0.5;  // below the dam crest
  sim.build();
  const real_t v0 = sim.swe->ownedVolume();
  for (int n = 1; n <= 50; ++n) {
    sim.step(0.1 * static_cast<real_t>(n));
  }
  // Face depths use the higher of the two bottoms (initialize.c:937-953):
  // the dam face area is zero, so nothing moves.
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->uu(), 1, 1), 0.0);
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->uu(), 1, 2), 0.0);
  EXPECT_NEAR(sim.swe->ownedVolume(), v0, 1.0e-12);
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->depth(), 1, 2), 0.0);  // crest cell dry
}

TEST(SweModule, ElevationOffsetLiftsNegativeBeds) {
  std::vector<real_t> bed = {-0.4, -0.2, 0.2};
  const std::string dir = writeRaster("lift.dat", bed);
  MiniSim sim;
  sim.cfg = baseConfig(3, 1, 1.0, 0.1);
  sim.cfg.configDir = dir;
  sim.cfg.domain.bottomElevation.fromFile = true;
  sim.cfg.domain.bottomElevation.file = "lift.dat";
  sim.cfg.initialConditions.surface.eta.constant = -0.1;
  sim.build();
  // offset = -min(bottom) (initialize.c:142-147); internal frame is shifted,
  // the output view is not.
  EXPECT_DOUBLE_EQ(sim.swe->elevationOffset(), 0.4);
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->bottom(), 1, 1), 0.0);
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->eta(), 1, 1), 0.3);
  EXPECT_DOUBLE_EQ(interiorValue(sim.swe->etaAbsolute(), 1, 1), -0.1);
}

}  // namespace
