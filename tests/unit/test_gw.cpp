/// \file test_gw.cpp
/// \brief Unit tests for the groundwater module (plan §8.1, §10 P2):
///        miniature configurations pinning conservation, the audit identity,
///        boundary-condition kinds, use_full3d, the post-allocation modes,
///        the adaptive-step controller, and the subsurface mesh geometry.

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Types.hpp"
#include "gw/RichardsSolver.hpp"
#include "gw/TerrainMetric.hpp"
#include "gw/VanGenuchten.hpp"

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

/// A single-rank miniature groundwater run assembled from an in-code
/// configuration.
struct MiniGw {
  frehg::FrehgConfig cfg;
  std::unique_ptr<frehg::Grid> grid;
  std::unique_ptr<frehg::HaloExchanger> halo;
  std::unique_ptr<frehg::gw::TerrainMetric> mesh;
  std::unique_ptr<frehg::BoundarySet> boundaries;
  std::unique_ptr<frehg::gw::RichardsSolver> gw;
  real_t t = 0.0;

  void build() {
    grid = std::make_unique<frehg::Grid>(MPI_COMM_WORLD, cfg.domain);
    halo = std::make_unique<frehg::HaloExchanger>(*grid, false);
    mesh = std::make_unique<frehg::gw::TerrainMetric>(*grid, cfg, *halo);
    grid->buildGlobalIds(mesh->ktop());
    boundaries =
        std::make_unique<frehg::BoundarySet>(*grid, cfg.boundaryConditions, cfg.configDir);
    gw = std::make_unique<frehg::gw::RichardsSolver>(*grid, cfg, *mesh, *boundaries, *halo);
  }

  /// One adaptive step; returns the dtg that was taken.
  real_t step() {
    const real_t dtg = gw->nextDt();
    t += dtg;
    gw->step(t, dtg);
    return dtg;
  }

  real_t volume() const { return gw->ownedVolume(); }
};

real_t interior3(const frehg::Field3<real_t>& field, int j, int i, int k) {
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  return host(static_cast<std::size_t>(j), static_cast<std::size_t>(i),
              static_cast<std::size_t>(k));
}

/// A soft test soil: sand-shaped retention with a slow saturated
/// conductivity so miniature columns evolve gently at fixed dtg.
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

frehg::FrehgConfig gwBaseConfig(int nx, int ny, int nz, real_t dz, real_t dtg) {
  frehg::FrehgConfig cfg;
  cfg.simulation.id = "gw-unit";
  cfg.domain.nx = nx;
  cfg.domain.ny = ny;
  cfg.domain.nz = nz;
  cfg.domain.dx = 1.0;
  cfg.domain.dy = 1.0;
  cfg.domain.dz = dz;
  cfg.domain.bottomElevation.constant = 0.0;
  cfg.time.dt = dtg;
  cfg.time.tEnd = 1000.0 * dtg;
  cfg.time.outputInterval = 1000.0 * dtg;
  cfg.modules.groundwater = true;
  cfg.groundwater.timestep.dtInit = dtg;
  cfg.groundwater.timestep.dtMin = dtg;
  cfg.groundwater.timestep.dtMax = dtg;
  cfg.groundwater.specificStorage = 0.0;
  cfg.soil.types = {testSoil()};
  cfg.soil.map.constantName = "test";
  cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
  cfg.initialConditions.groundwater.value.constant = 0.2;
  return cfg;
}

/// Polygon covering the full (single-cell-wide) footprint.
std::vector<std::array<real_t, 2>> wholeFootprint(int nx, int ny) {
  const real_t x1 = static_cast<real_t>(nx) + 0.1;
  const real_t y1 = static_cast<real_t>(ny) + 0.1;
  return {{-0.1, -0.1}, {x1, -0.1}, {x1, y1}, {-0.1, y1}};
}

frehg::BoundaryConditionConfig topBc(frehg::BcKind kind, real_t value, int nx, int ny) {
  frehg::BoundaryConditionConfig bc;
  bc.name = "top";
  bc.polygon = wholeFootprint(nx, ny);
  bc.target = frehg::BcTarget::GroundwaterTop;
  bc.kind = kind;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = value;
  return bc;
}

/// Writes a flat-list data file into a fresh temp directory; returns the dir.
std::string writeFlatFile(const std::string& name, const std::vector<real_t>& values) {
  const std::filesystem::path dir =
      std::filesystem::temp_directory_path() /
      ("frehg_gw_test_" + std::to_string(::getpid()) + "_" + name);
  std::filesystem::create_directories(dir);
  std::ofstream out(dir / name);
  out.precision(17);
  for (const real_t v : values) {
    out << v << "\n";
  }
  return dir.string();
}

// ---------------------------------------------------------------------------
// Conservation and the audit identity (plan §10 P2 exit criteria).
// ---------------------------------------------------------------------------

TEST(GwModule, ClosedColumnConservesMassToOneInEight) {
  // No boundary conditions at all: every face is no-flux. Uniform moisture
  // drains internally under gravity; the total volume must not drift by
  // more than 1e-8 relative per step (plan §10 P2).
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 20, 0.05, 0.5);
  sim.build();
  const real_t v0 = sim.volume();
  real_t prev = v0;
  for (int n = 0; n < 50; ++n) {
    sim.step();
    const real_t v = sim.volume();
    EXPECT_LE(std::fabs(v - prev), 1.0e-8 * v0) << "step " << n;
    EXPECT_EQ(sim.gw->audit().boundaryIn, 0.0);
    prev = v;
  }
  // The state genuinely evolves (water moves down; heads change).
  EXPECT_GT(interior3(sim.gw->waterContent(), 1, 1, 19), 0.2);
  EXPECT_LT(interior3(sim.gw->waterContent(), 1, 1, 0), 0.2);
}

TEST(GwModule, AuditIdentityClosesWithCompressibleStorage) {
  // With Ss > 0 the θ update moves water into compressible storage; the
  // audit identity dV = boundary_in - ss_storage + realloc - vloss closes
  // to rounding every step.
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 20, 0.05, 0.5);
  sim.cfg.groundwater.specificStorage = 1.0e-4;
  sim.build();
  real_t prev = sim.volume();
  for (int n = 0; n < 30; ++n) {
    sim.step();
    const real_t v = sim.volume();
    const frehg::gw::GwStepAudit& a = sim.gw->audit();
    const real_t residual =
        (v - prev) - (a.boundaryIn - a.ssStorage + a.reallocAdjust - a.vloss);
    EXPECT_LE(std::fabs(residual), 1.0e-12 * std::max(v, 1.0)) << "step " << n;
    prev = v;
  }
}

// ---------------------------------------------------------------------------
// Boundary-condition kinds (plan §10 P2 deliverable).
// ---------------------------------------------------------------------------

TEST(GwModule, PrescribedHeadTopInfiltrates) {
  // b2 in miniature: ponded head 0 on the column top. Water enters through
  // the saturated top face, the top cell saturates, and the audit sees the
  // inflow.
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 10, 0.05, 0.5);
  sim.cfg.initialConditions.groundwater.value.constant = 0.06;
  sim.cfg.boundaryConditions = {topBc(frehg::BcKind::Head, 0.0, 1, 1)};
  sim.build();
  const real_t v0 = sim.volume();
  real_t inflow = 0.0;
  for (int n = 0; n < 40; ++n) {
    sim.step();
    inflow += sim.gw->audit().boundaryIn;
  }
  EXPECT_GT(inflow, 0.0);
  EXPECT_GT(sim.volume(), v0);
  // The top cell wets steadily toward saturation through the saturated
  // (K = Ks) top face; full saturation takes longer than this miniature run.
  EXPECT_GT(interior3(sim.gw->waterContent(), 1, 1, 0), 0.12);
}

TEST(GwModule, PrescribedFluxTopDeliversExactVolume) {
  // The flux override replaces the Darcy value exactly (darcy_flux
  // subroutines.c:183): total inflow is q A t to rounding. Sign follows
  // legacy qtop: negative = downward (into the ground).
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 10, 0.1, 0.5);
  sim.cfg.boundaryConditions = {topBc(frehg::BcKind::Flux, -1.0e-6, 1, 1)};
  sim.build();
  const real_t v0 = sim.volume();
  real_t elapsed = 0.0;
  for (int n = 0; n < 40; ++n) {
    elapsed += sim.step();
  }
  const real_t expected = 1.0e-6 * elapsed;  // area = dx dy = 1
  EXPECT_NEAR(sim.volume() - v0, expected, 1.0e-12 * std::max(expected, 1.0));
}

TEST(GwModule, GravityBottomDrainsAtConductivityRate) {
  // Free drainage (kind flux, value gravity): the bottom face releases
  // K(h) r_visc r_rho A (darcy_flux subroutines.c:190, with the bottom-code
  // fix). A saturated column drains at the saturated conductivity.
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 10, 0.1, 0.5);
  // Ss regularizes the fully saturated column: with zero storage everywhere
  // and only Neumann conditions the head system is singular (legacy
  // included: C(h > 0) = 0), so saturated cases carry the benchmark Ss.
  sim.cfg.groundwater.specificStorage = 1.0e-5;
  sim.cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::WaterTable;
  sim.cfg.initialConditions.groundwater.value.constant = 0.0;  // table at the surface
  frehg::BoundaryConditionConfig bc;
  bc.name = "drain";
  bc.polygon = wholeFootprint(1, 1);
  bc.target = frehg::BcTarget::GroundwaterBottom;
  bc.kind = frehg::BcKind::Flux;
  bc.value.form = frehg::BcValueConfig::Form::Gravity;
  sim.cfg.boundaryConditions = {bc};
  sim.build();
  const real_t v0 = sim.volume();
  const real_t dtg = sim.step();
  // The flux override is exact: q = -K r_visc r_rho A with K = Ks in the
  // saturated column, so the audit records Ks A dtg outflow to rounding.
  // The volume change itself is dominated by the vadose zone the falling
  // pressure creates through the consistency restore (an audited PCA
  // behavior, not boundary flux), so the assertion is the audit identity.
  EXPECT_NEAR(sim.gw->audit().boundaryIn, -1.0e-5 * dtg, 1.0e-18);
  const frehg::gw::GwStepAudit& a = sim.gw->audit();
  const real_t residual = (sim.volume() - v0) -
                          (a.boundaryIn - a.ssStorage + a.reallocAdjust - a.vloss);
  EXPECT_LE(std::fabs(residual), 1.0e-12);
  EXPECT_LT(sim.volume(), v0);
}

TEST(GwModule, PrescribedFluxBottomInjects) {
  // Bottom flux, positive = upward = into the domain (legacy qbot).
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 10, 0.1, 0.5);
  frehg::BoundaryConditionConfig bc;
  bc.name = "recharge";
  bc.polygon = wholeFootprint(1, 1);
  bc.target = frehg::BcTarget::GroundwaterBottom;
  bc.kind = frehg::BcKind::Flux;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = 2.0e-6;
  sim.cfg.boundaryConditions = {bc};
  sim.build();
  const real_t v0 = sim.volume();
  real_t elapsed = 0.0;
  for (int n = 0; n < 20; ++n) {
    elapsed += sim.step();
  }
  EXPECT_NEAR(sim.volume() - v0, 2.0e-6 * elapsed, 1.0e-12);
}

TEST(GwModule, HydrostaticSideHoldsEquilibrium) {
  // A hydrostatic side ghost (enforce_head_bc:769-788) at the water-table
  // stage of the initial condition leaves the column in equilibrium: heads
  // do not drift.
  MiniGw sim;
  sim.cfg = gwBaseConfig(3, 3, 10, 0.1, 0.5);
  sim.cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::WaterTable;
  sim.cfg.initialConditions.groundwater.value.constant = -0.35;
  frehg::BoundaryConditionConfig bc;
  bc.name = "sea";
  // Middle cell of the west edge only (it touches exactly one edge face).
  bc.polygon = {{-0.1, 1.0 + 0.1}, {0.9, 1.1}, {0.9, 1.9}, {-0.1, 1.9}};
  bc.target = frehg::BcTarget::GroundwaterSide;
  bc.kind = frehg::BcKind::Head;
  bc.value.form = frehg::BcValueConfig::Form::Hydrostatic;
  bc.value.hydrostaticEta = -0.35;
  sim.cfg.boundaryConditions = {bc};
  sim.build();
  auto headBefore = Kokkos::create_mirror_view(sim.gw->head());
  Kokkos::deep_copy(headBefore, sim.gw->head());
  for (int n = 0; n < 10; ++n) {
    sim.step();
  }
  auto headAfter = Kokkos::create_mirror_view(sim.gw->head());
  Kokkos::deep_copy(headAfter, sim.gw->head());
  real_t drift = 0.0;
  for (int j = 1; j <= 3; ++j) {
    for (int i = 1; i <= 3; ++i) {
      for (int k = 0; k < 10; ++k) {
        drift = std::max(drift, std::fabs(headAfter(j, i, k) - headBefore(j, i, k)));
      }
    }
  }
  EXPECT_LE(drift, 1.0e-6);
}

TEST(GwModule, SideFluxInjectsWithInflowPositive) {
  // frehg2 side-flux convention: positive value = into the domain
  // (documented in parameters.md; the legacy qyp/qym signs disagreed
  // between predictor and corrector and were never exercised).
  MiniGw sim;
  sim.cfg = gwBaseConfig(3, 3, 5, 0.1, 0.5);
  frehg::BoundaryConditionConfig bc;
  bc.name = "inflow";
  bc.polygon = {{-0.1, 1.1}, {0.9, 1.1}, {0.9, 1.9}, {-0.1, 1.9}};
  bc.target = frehg::BcTarget::GroundwaterSide;
  bc.kind = frehg::BcKind::Flux;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = 4.0e-6;
  sim.cfg.boundaryConditions = {bc};
  sim.build();
  const real_t v0 = sim.volume();
  real_t elapsed = 0.0;
  for (int n = 0; n < 20; ++n) {
    elapsed += sim.step();
  }
  // One column's west faces: area = dy * total depth = 1 * 0.5.
  EXPECT_NEAR(sim.volume() - v0, 4.0e-6 * 0.5 * elapsed, 1.0e-12);
}

// ---------------------------------------------------------------------------
// Scheme behaviors: use_full3d, post-allocation modes, adaptive stepping.
// ---------------------------------------------------------------------------

TEST(GwModule, UseFull3dSwitchGatesLateralFlow) {
  // Two columns with different moisture; lateral exchange happens only with
  // use_full3d (compute_K_face:295-307 seals unsaturated lateral faces).
  const std::vector<real_t> moisture = {0.10, 0.10, 0.10, 0.10, 0.30, 0.30, 0.30, 0.30};
  const std::string dir = writeFlatFile("moisture.dat", moisture);
  for (const bool full3d : {false, true}) {
    MiniGw sim;
    sim.cfg = gwBaseConfig(2, 1, 4, 0.1, 0.5);
    sim.cfg.groundwater.useFull3d = full3d;
    sim.cfg.configDir = dir;
    sim.cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::Moisture;
    sim.cfg.initialConditions.groundwater.value.fromFile = true;
    sim.cfg.initialConditions.groundwater.value.file = "moisture.dat";
    sim.build();
    sim.step();
    const real_t qInterface = interior3(sim.gw->fluxXPerArea(), 1, 1, 1);
    if (full3d) {
      EXPECT_NE(qInterface, 0.0);
    } else {
      EXPECT_EQ(qInterface, 0.0);
    }
  }
}

TEST(GwModule, ReallocationSurplusModesDropOrConserve) {
  // The b3 mechanism in miniature: ponded infiltration through a sand layer
  // onto dry, flat-retention clay. The clay cells behind the saturation
  // front carry corrector moisture above theta(h_predicted) — the predictor
  // barely moves their head, since the interface conductivity is the clay's
  // (darcy z faces use the lower cell's Ks). The legacy sweep discards that
  // surplus (drop; groundwater.c:1026 left the transfer disabled);
  // redistribute walks it into the column's pore room (amendment A7). Both
  // modes must close the audit identity every step.
  const std::vector<real_t> ids = {0, 0, 0, 1, 1, 1, 1, 1, 1, 1};
  const std::string dir = writeFlatFile("soil_id.dat", ids);
  real_t droppedTotal[2] = {0.0, 0.0};
  for (const bool redistribute : {false, true}) {
    MiniGw sim;
    sim.cfg = gwBaseConfig(1, 1, 10, 0.05, 0.5);
    frehg::SoilType sand = testSoil(3.0e-4);
    sand.name = "sand";
    frehg::SoilType clay = testSoil(1.52e-6);
    clay.name = "clay";
    clay.thetaS = 0.4686;
    clay.thetaR = 0.106;
    clay.vgAlpha = 1.04;
    clay.vgN = 1.395;
    clay.aev = 0.0;
    sim.cfg.soil.types = {sand, clay};
    sim.cfg.soil.map.fromFile = true;
    sim.cfg.soil.map.file = "soil_id.dat";
    sim.cfg.soil.map.constantName.clear();
    sim.cfg.configDir = dir;
    sim.cfg.initialConditions.groundwater.value.constant = 0.14;
    sim.cfg.groundwater.reallocationSurplus =
        redistribute ? frehg::GroundwaterConfig::ReallocationSurplus::Redistribute
                     : frehg::GroundwaterConfig::ReallocationSurplus::Drop;
    sim.cfg.boundaryConditions = {topBc(frehg::BcKind::Head, 0.0, 1, 1)};
    sim.build();
    const real_t v0 = sim.volume();
    real_t prev = v0;
    for (int n = 0; n < 200; ++n) {
      sim.step();
      const real_t v = sim.volume();
      const frehg::gw::GwStepAudit& a = sim.gw->audit();
      const real_t residual =
          (v - prev) - (a.boundaryIn - a.ssStorage + a.reallocAdjust - a.vloss);
      ASSERT_LE(std::fabs(residual), 1.0e-12) << "mode " << redistribute << " step " << n;
      droppedTotal[redistribute ? 1 : 0] += a.reallocDropped;
      prev = v;
    }
    EXPECT_GT(sim.volume(), v0);
  }
  // Drop mode discards the clay front's surplus outright; redistribute
  // places it in the column (a 1D column has no lateral fraction to lose).
  EXPECT_GT(droppedTotal[0], 1.0e-5);
  EXPECT_LE(droppedTotal[1], 1.0e-9);
}

TEST(GwModule, AdaptiveStepGrowsInQuietFlowAndStaysClamped) {
  // Equilibrium column: dq = 0, so dtg grows by the legacy factor 1.25 per
  // step until dt_max clamps it (adaptive_time_step:1709-1719).
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 8, 0.1, 0.1);
  sim.cfg.groundwater.timestep.dtInit = 0.1;
  sim.cfg.groundwater.timestep.dtMin = 0.1;
  sim.cfg.groundwater.timestep.dtMax = 1.0;
  sim.cfg.groundwater.specificStorage = 1.0e-5;  // regularizes the saturated column
  sim.cfg.initialConditions.groundwater.form = frehg::GroundwaterInitialConfig::Form::WaterTable;
  sim.cfg.initialConditions.groundwater.value.constant = 0.0;  // saturated equilibrium
  sim.build();
  EXPECT_EQ(sim.gw->nextDt(), 0.1);
  sim.step();
  EXPECT_NEAR(sim.gw->nextDt(), 0.125, 1.0e-12);
  real_t dtg = 0.0;
  for (int n = 0; n < 30; ++n) {
    dtg = sim.step();
    EXPECT_LE(sim.gw->nextDt(), 1.0);
    EXPECT_GE(sim.gw->nextDt(), 0.1);
  }
  EXPECT_EQ(dtg, 1.0);  // reached and held dt_max
}

TEST(GwModule, ReallocateRestoresIsolatedHeadConsistency) {
  // After a step, isolated unsaturated cells carry h = h(θ)
  // (reallocate_water_content:995-999).
  MiniGw sim;
  sim.cfg = gwBaseConfig(1, 1, 10, 0.05, 0.5);
  sim.build();
  sim.step();
  for (int k = 2; k < 8; ++k) {
    const real_t wc = interior3(sim.gw->waterContent(), 1, 1, k);
    const real_t h = interior3(sim.gw->head(), 1, 1, k);
    const frehg::gw::VgSoil soil{5.9, 2.68, 0.46, 0.04, -0.02};
    EXPECT_NEAR(h, frehg::gw::headFromWaterContent(soil, wc), 1.0e-12) << "k=" << k;
  }
}

// ---------------------------------------------------------------------------
// Subsurface mesh geometry (TerrainMetric; legacy map.c:196-616).
// ---------------------------------------------------------------------------

TEST(GwMesh, RegularMeshAppliesLegacyPartialCellRules) {
  // Three columns against a box anchored at the highest bed: full, a
  // shortened crossing layer (>= dz/4 remains), and a merged crossing layer
  // (< dz/4) — legacy map.c:357-393.
  const std::string dir = writeFlatFile("bath.dat", {0.0, -0.05, -0.09});
  frehg::FrehgConfig cfg = gwBaseConfig(3, 1, 3, 0.1, 0.5);
  cfg.configDir = dir;
  cfg.domain.bottomElevation.fromFile = true;
  cfg.domain.bottomElevation.file = "bath.dat";
  frehg::Grid grid(MPI_COMM_WORLD, cfg.domain);
  frehg::HaloExchanger halo(grid, false);
  frehg::gw::TerrainMetric mesh(grid, cfg, halo);

  EXPECT_EQ(mesh.ktop()(0, 0), 0);
  EXPECT_EQ(mesh.ktop()(0, 1), 0);  // 0.05 of the first layer remains
  EXPECT_EQ(mesh.ktop()(0, 2), 1);  // 0.01 < dz/4: merged into the layer below
  EXPECT_EQ(mesh.elevationOffset(), 0.09);
  EXPECT_EQ(mesh.boxTop(), 0.09);
  EXPECT_NEAR(interior3(mesh.dz3d(), 1, 1, 0), 0.1, 1.0e-14);
  EXPECT_NEAR(interior3(mesh.dz3d(), 1, 2, 0), 0.05, 1.0e-14);
  EXPECT_NEAR(interior3(mesh.dz3d(), 1, 3, 1), 0.11, 1.0e-14);  // absorbed 0.01
  // zcell reports centers in the un-shifted datum.
  EXPECT_NEAR(mesh.zCellHost()(0, 0, 0), -0.05, 1.0e-14);
  EXPECT_NEAR(mesh.zCellHost()(0, 1, 0), -0.075, 1.0e-14);
}

TEST(GwMesh, TerrainMeshScalesColumnsAndTiltsFaces) {
  // follow_terrain: each column spans [bed, box bottom] with nz layers;
  // lateral faces carry slope metrics (legacy map.c:316-353, :530-581).
  const std::string dir = writeFlatFile("bath.dat", {0.0, -0.1});
  frehg::FrehgConfig cfg = gwBaseConfig(2, 1, 4, 0.1, 0.5);
  cfg.configDir = dir;
  cfg.domain.followTerrain = true;
  cfg.domain.bottomElevation.fromFile = true;
  cfg.domain.bottomElevation.file = "bath.dat";
  frehg::Grid grid(MPI_COMM_WORLD, cfg.domain);
  frehg::HaloExchanger halo(grid, false);
  frehg::gw::TerrainMetric mesh(grid, cfg, halo);

  // Box: top at max bath (offset frame 0.1), bottom 0.4 below.
  EXPECT_EQ(mesh.elevationOffset(), 0.1);
  EXPECT_EQ(mesh.ktop()(0, 0), 0);
  EXPECT_EQ(mesh.ktop()(0, 1), 0);
  EXPECT_NEAR(interior3(mesh.dz3d(), 1, 1, 0), 0.1, 1.0e-14);    // (0.1+0.3)/4
  EXPECT_NEAR(interior3(mesh.dz3d(), 1, 2, 0), 0.075, 1.0e-14);  // (0.0+0.3)/4
  // Face metric between the columns: cell centers at z = 0.05 and 0.0375
  // (offset frame 0.05... first layer centers 0.1-0.05 and 0.0-0.0375+...).
  const real_t zc0 = 0.1 - 0.5 * 0.1;
  const real_t zc1 = 0.0 - 0.5 * 0.075;
  const real_t hdiff = std::fabs(zc1 - zc0);
  const real_t dist = std::sqrt(hdiff * hdiff + 1.0);
  EXPECT_NEAR(interior3(mesh.sinX(), 1, 1, 0), hdiff / dist, 1.0e-12);
  EXPECT_NEAR(interior3(mesh.cosX(), 1, 1, 0), 1.0 / dist, 1.0e-12);
  EXPECT_NEAR(interior3(mesh.areaX(), 1, 1, 0), 0.5 * (0.1 + 0.075) * (1.0 / dist) * 1.0,
              1.0e-12);
}

TEST(GwMesh, UniformTerrainLayersStackBelowEachBed) {
  // terrain_layers: uniform (amendment A12): every column carries the
  // configured dz * dz_stretch^k profile below its own bed — the
  // terrain-parallel slab of the b5 reference (SERGHEI GwInit.h:109-117) —
  // while the slope metrics stay the legacy terrain-following forms.
  const std::string dir = writeFlatFile("bath.dat", {0.0, -0.1});
  frehg::FrehgConfig cfg = gwBaseConfig(2, 1, 4, 0.1, 0.5);
  cfg.configDir = dir;
  cfg.domain.followTerrain = true;
  cfg.domain.terrainLayers = frehg::DomainConfig::TerrainLayers::Uniform;
  cfg.domain.bottomElevation.fromFile = true;
  cfg.domain.bottomElevation.file = "bath.dat";
  frehg::Grid grid(MPI_COMM_WORLD, cfg.domain);
  frehg::HaloExchanger halo(grid, false);
  frehg::gw::TerrainMetric mesh(grid, cfg, halo);

  EXPECT_EQ(mesh.ktop()(0, 0), 0);
  EXPECT_EQ(mesh.ktop()(0, 1), 0);
  for (int i = 1; i <= 2; ++i) {
    for (int k = 0; k < 4; ++k) {
      EXPECT_NEAR(interior3(mesh.dz3d(), 1, i, k), 0.1, 1.0e-14) << i << "," << k;
    }
  }
  // Column bottoms follow the local bed: bath - (k+1) dz in the offset
  // frame (offset 0.1 lifts the beds to 0.1 and 0.0).
  EXPECT_NEAR(interior3(mesh.bot3d(), 1, 1, 3), 0.1 - 0.4, 1.0e-14);
  EXPECT_NEAR(interior3(mesh.bot3d(), 1, 2, 3), 0.0 - 0.4, 1.0e-14);
  // The slab is terrain-parallel, so the face slope equals the bed slope.
  const real_t hdiff = 0.1;
  const real_t dist = std::sqrt(hdiff * hdiff + 1.0);
  EXPECT_NEAR(interior3(mesh.sinX(), 1, 1, 0), hdiff / dist, 1.0e-12);
  // zcell reports centers in the un-shifted datum, below each local bed.
  EXPECT_NEAR(mesh.zCellHost()(0, 0, 0), -0.05, 1.0e-14);
  EXPECT_NEAR(mesh.zCellHost()(0, 1, 0), -0.15, 1.0e-14);
}

}  // namespace
