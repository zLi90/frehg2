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
#include "swe/WindForcing.hpp"

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

// ---------------------------------------------------------------------------
// v2 Q6 wind: Cd(U10) laws (hand tables per the published formulas, plan
// §5.3 unit battery) and the WindForcing direction conventions including
// the 360->0 wrap.
// ---------------------------------------------------------------------------

TEST(WindForcing, DragLawTablesMatchHandValues) {
  using frehg::WindConfig;
  using frehg::swe::windDragCoefficient;
  const frehg::real_t cap = 3.5e-3;
  struct Row {
    frehg::real_t u10;
    frehg::real_t garratt, smith, wu, largePond;
  };
  // Hand-computed from the published formulas (1e-12 relative, plan §5.3):
  //   garratt   = (0.75 + 0.067 U) 1e-3   (cap 3.5e-3 above ~41.0 m/s)
  //   smith     = (0.63 + 0.066 U) 1e-3
  //   wu        = (0.80 + 0.065 U) 1e-3
  //   large-pond= 1.2e-3 (U < 11) | (0.49 + 0.065 U) 1e-3, both capped
  const Row rows[] = {
      {0.0, 0.75e-3, 0.63e-3, 0.80e-3, 1.2e-3},
      {5.0, 1.085e-3, 0.96e-3, 1.125e-3, 1.2e-3},
      {10.0, 1.42e-3, 1.29e-3, 1.45e-3, 1.2e-3},
      {20.0, 2.09e-3, 1.95e-3, 2.10e-3, 1.79e-3},
      {30.0, 2.76e-3, 2.61e-3, 2.75e-3, 2.44e-3},
      {40.0, 3.43e-3, 3.27e-3, 3.40e-3, 3.09e-3},
      {60.0, 3.5e-3, 3.5e-3, 3.5e-3, 3.5e-3},  // every law capped
  };
  for (const Row& r : rows) {
    EXPECT_NEAR(windDragCoefficient(WindConfig::DragLaw::Garratt, r.u10, cap, 0.0),
                r.garratt, 1.0e-12 * r.garratt) << "garratt U=" << r.u10;
    EXPECT_NEAR(windDragCoefficient(WindConfig::DragLaw::SmithBanke, r.u10, cap, 0.0),
                r.smith, 1.0e-12 * r.smith) << "smith-banke U=" << r.u10;
    EXPECT_NEAR(windDragCoefficient(WindConfig::DragLaw::Wu, r.u10, cap, 0.0),
                r.wu, 1.0e-12 * r.wu) << "wu U=" << r.u10;
    EXPECT_NEAR(windDragCoefficient(WindConfig::DragLaw::LargePond, r.u10, cap, 0.0),
                r.largePond, 1.0e-12 * r.largePond) << "large-pond U=" << r.u10;
  }
  // The Large & Pond breakpoint: the published law's 5e-6 jump at 11 m/s.
  EXPECT_DOUBLE_EQ(windDragCoefficient(WindConfig::DragLaw::LargePond, 10.999, cap, 0.0),
                   1.2e-3);
  EXPECT_NEAR(windDragCoefficient(WindConfig::DragLaw::LargePond, 11.0, cap, 0.0),
              1.205e-3, 1.0e-15);
  // The constant law is the uncapped legacy Cw.
  EXPECT_DOUBLE_EQ(windDragCoefficient(WindConfig::DragLaw::Constant, 60.0, cap, 1.0277551),
                   1.0277551);
}

TEST(WindForcing, DirectionSeriesInterpolatesAcrossTheWrap) {
  // A direction series stepping 350 -> 10 degrees must pass through
  // 0/360, never through 180 (the documented circle convention).
  const std::string dir = ::testing::TempDir() + "/wind_dir_wrap.dat";
  {
    std::ofstream out(dir);
    out << "0.0 350.0\n100.0 10.0\n";
  }
  frehg::WindConfig cfg;
  cfg.enabled = true;
  cfg.speed.constant = 10.0;
  cfg.direction.fromSeries = true;
  cfg.direction.file = "wind_dir_wrap.dat";
  const frehg::swe::WindForcing forcing(
      cfg, [](const std::string& f) { return ::testing::TempDir() + "/" + f; });
  // Midpoint lands on 0/360 (the short arc's centre); the truncated
  // legacy pi literal bounds the round-trip at ~5e-8 rad.
  const frehg::real_t omegaMid = forcing.sample(50.0).omega;
  EXPECT_NEAR(std::sin(omegaMid), 0.0, 1.0e-6);
  EXPECT_NEAR(std::cos(omegaMid), 1.0, 1.0e-6);
  // The wrap guarantee: every interpolated direction stays inside the
  // short arc [350, 10] (chord interpolation of the unit vectors —
  // within 0.04 deg of constant rate for this 20-degree step; it NEVER
  // takes the 180-degree long way the legacy linear-in-degrees
  // interpolation took).
  for (frehg::real_t tt = 0.0; tt <= 100.0; tt += 5.0) {
    const frehg::real_t omega = forcing.sample(tt).omega;
    EXPECT_GT(std::cos(omega), std::cos(10.5 * frehg::swe::kLegacyPi / 180.0))
        << "t = " << tt;
  }
  const frehg::real_t omegaQ = forcing.sample(25.0).omega;  // ~355 deg (chord)
  EXPECT_NEAR(omegaQ * 180.0 / frehg::swe::kLegacyPi, -5.0, 0.05);
  std::remove(dir.c_str());
}

TEST(WindForcing, ComponentFormMatchesCompassForm) {
  // (u10, v10) = W (cos a, sin a) must give the same sample as the
  // compass form with direction a and north_angle 0 (both feed the
  // momentum as omega-from-+x), and the law sees the same |U10|.
  frehg::WindConfig compass;
  compass.enabled = true;
  compass.law = frehg::WindConfig::DragLaw::Garratt;
  compass.speed.constant = 13.0;
  compass.direction.constant = 37.0;
  const frehg::swe::WindForcing a(compass, [](const std::string& f) { return f; });

  frehg::WindConfig comp;
  comp.enabled = true;
  comp.law = frehg::WindConfig::DragLaw::Garratt;
  comp.componentForm = true;
  const frehg::real_t rad = 37.0 * frehg::swe::kLegacyPi / 180.0;
  comp.u10.constant = 13.0 * std::cos(rad);
  comp.v10.constant = 13.0 * std::sin(rad);
  const frehg::swe::WindForcing b(comp, [](const std::string& f) { return f; });

  const auto sa = a.sample(0.0);
  const auto sb = b.sample(0.0);
  EXPECT_NEAR(sb.speed, sa.speed, 1.0e-12);
  EXPECT_NEAR(sb.omega, sa.omega, 1.0e-9);
  EXPECT_NEAR(sb.dragCd, sa.dragCd, 1.0e-15);
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
