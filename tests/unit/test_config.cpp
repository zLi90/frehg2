/// \file test_config.cpp
/// \brief Configuration schema tests: acceptance of the six benchmark
///        configurations and rejection of a battery of crafted-invalid
///        configurations with the correct messages (plan §8.1).

#include "core/Config.hpp"
#include "core/Types.hpp"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

using frehg::ValidationResult;

/// A minimal valid surface-water-only configuration used as the mutation
/// base for the invalid battery. Files it references are created by the
/// fixture.
const char* kBaseConfig = R"(simulation: {id: unit-test}
domain:
  nx: 4
  ny: 3
  nz: 2
  dx: 1.0
  dy: 1.0
  dz: 0.5
  bottom_elevation: {constant: -1.0}
time: {dt: 0.5, t_end: 100.0, output_interval: 10.0}
modules: {surface_water: true}
surface_water:
  friction: {coefficient: {constant: 0.02}}
  min_depth: 1.0e-6
  wetting_face_depth: 1.0e-6
initial_conditions: {surface: {eta: {constant: 0.0}}}
output: {filename: out/output.h5}
)";

class ConfigTest : public ::testing::Test {
 protected:
  std::filesystem::path dir_;

  void SetUp() override {
    dir_ = std::filesystem::path(::testing::TempDir()) / "frehg_config_test";
    std::filesystem::create_directories(dir_);
  }
  void TearDown() override { std::filesystem::remove_all(dir_); }

  std::string writeConfig(const std::string& text) {
    const std::filesystem::path path = dir_ / "config.yaml";
    std::ofstream out(path);
    out << text;
    out.close();
    return path.string();
  }

  ValidationResult validate(const std::string& text) {
    return frehg::validateConfigFile(writeConfig(text));
  }

  static bool hasError(const ValidationResult& result, const std::string& needle) {
    for (const std::string& error : result.errors) {
      if (error.find(needle) != std::string::npos) {
        return true;
      }
    }
    return false;
  }

  static std::string joined(const ValidationResult& result) {
    std::string all;
    for (const std::string& error : result.errors) {
      all += error + "\n";
    }
    return all;
  }
};

TEST_F(ConfigTest, MinimalConfigIsValid) {
  const ValidationResult result = validate(kBaseConfig);
  EXPECT_TRUE(result.ok()) << joined(result);
}

TEST_F(ConfigTest, LoadMaterializesDefaults) {
  const frehg::FrehgConfig cfg = frehg::loadConfig(writeConfig(kBaseConfig));
  EXPECT_EQ(cfg.simulation.id, "unit-test");
  EXPECT_DOUBLE_EQ(cfg.surfaceWater.gravity, 9.81);
  EXPECT_DOUBLE_EQ(cfg.surfaceWater.friction.thinLayerDepth, 0.1);
  EXPECT_EQ(cfg.surfaceWater.friction.law, frehg::FrictionConfig::Law::Manning);
  EXPECT_DOUBLE_EQ(cfg.surfaceWater.viscosityX, 1.0e-6);
  EXPECT_EQ(cfg.coupling.mode, frehg::CouplingConfig::Mode::Sync);
  EXPECT_EQ(cfg.runtime.gpuAwareMpi, frehg::RuntimeConfig::GpuAwareMpi::Auto);
  EXPECT_EQ(cfg.domain.decomposition.mpiNx, 0);  // auto
  EXPECT_DOUBLE_EQ(cfg.time.tStart, 0.0);
  const std::string description = frehg::describeConfig(cfg);
  EXPECT_NE(description.find("unit-test"), std::string::npos);
  EXPECT_NE(description.find("manning"), std::string::npos);
}

TEST_F(ConfigTest, AllSixBenchmarkConfigsValidate) {
  const std::vector<std::string> cases = {
      "benchmarks/b1-sw/b1-sw.yaml",
      "benchmarks/b2-gw/b2-gw.yaml",
      "benchmarks/b3-kirkland/b3-kirkland.yaml",
      "benchmarks/b4-govindaraju/b4-govindaraju.yaml",
      "benchmarks/b5-vcatchment/b5-vcatchment.yaml",
      "benchmarks/b6-kuan/b6-kuan-ss.yaml",
      "benchmarks/b6-kuan/b6-kuan-td.yaml",
  };
  for (const std::string& rel : cases) {
    const std::string path = std::string(FREHG_REPO_DIR) + "/" + rel;
    const ValidationResult result = frehg::validateConfigFile(path);
    EXPECT_TRUE(result.ok()) << rel << ":\n" << joined(result);
    if (result.ok()) {
      const frehg::FrehgConfig cfg = frehg::loadConfig(path);
      EXPECT_FALSE(cfg.simulation.id.empty()) << rel;
      EXPECT_FALSE(frehg::describeConfig(cfg).empty()) << rel;
    }
  }
}

// --------------------------------------------------------------------------
// Invalid battery (>= 10 distinct rejections with message checks)
// --------------------------------------------------------------------------

TEST_F(ConfigTest, RejectsUnknownKeyWithSuggestion) {
  std::string text = kBaseConfig;
  text += "extra:\n  value: 1\n";
  ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "unknown key 'extra'")) << joined(result);

  // Misspelled nested key gets a nearest-key suggestion.
  std::string misspelled = kBaseConfig;
  const std::string needle = "min_depth:";
  misspelled.replace(misspelled.find(needle), needle.size(), "min_dept:");
  result = validate(misspelled);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "unknown key 'min_dept'")) << joined(result);
  EXPECT_TRUE(hasError(result, "did you mean 'min_depth'?")) << joined(result);
}

TEST_F(ConfigTest, RejectsMissingRequiredKey) {
  std::string text = kBaseConfig;
  const std::string needle = "time: {dt: 0.5, t_end: 100.0, output_interval: 10.0}";
  text.replace(text.find(needle), needle.size(),
               "time: {t_end: 100.0, output_interval: 10.0}");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "time: missing required key 'dt'")) << joined(result);
}

TEST_F(ConfigTest, RejectsWrongType) {
  std::string text = kBaseConfig;
  const std::string needle = "nx: 4";
  text.replace(text.find(needle), needle.size(), "nx: four");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "domain.nx: expected an integer")) << joined(result);
}

TEST_F(ConfigTest, RejectsOutOfRangeValue) {
  std::string text = kBaseConfig;
  const std::string needle = "dx: 1.0";
  text.replace(text.find(needle), needle.size(), "dx: -1.0");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "domain.dx: must be > 0")) << joined(result);
}

TEST_F(ConfigTest, RejectsInvalidEnumValue) {
  std::string text = kBaseConfig;
  const std::string needle = "friction: {coefficient: {constant: 0.02}}";
  text.replace(text.find(needle), needle.size(),
               "friction: {law: cheezy, coefficient: {constant: 0.02}}");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "invalid value 'cheezy'")) << joined(result);
  EXPECT_TRUE(hasError(result, "manning, chezy")) << joined(result);
}

TEST_F(ConfigTest, RejectsGroundwaterTimestepOrdering) {
  std::string text = kBaseConfig;
  const std::string needle = "modules: {surface_water: true}";
  text.replace(text.find(needle), needle.size(),
               "modules: {surface_water: true, groundwater: true}");
  text += R"(groundwater:
  timestep: {dt_init: 5.0, dt_min: 0.1, dt_max: 1.0}
  specific_storage: 1.0e-5
soil:
  types:
    - {name: loam, ksx: 1.0e-5, ksy: 1.0e-5, ksz: 1.0e-5,
       theta_s: 0.4, theta_r: 0.08, vg_alpha: 6.0, vg_n: 2.0}
  map: {constant: loam}
)";
  ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "dt_min <= dt_init <= dt_max")) << joined(result);
}

TEST_F(ConfigTest, RejectsNonIntegerOutputInterval) {
  std::string text = kBaseConfig;
  const std::string needle = "output_interval: 10.0";
  text.replace(text.find(needle), needle.size(), "output_interval: 10.5");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "whole seconds")) << joined(result);
}

TEST_F(ConfigTest, RejectsMissingReferencedFile) {
  std::string text = kBaseConfig;
  const std::string needle = "bottom_elevation: {constant: -1.0}";
  text.replace(text.find(needle), needle.size(),
               "bottom_elevation: {file: does_not_exist.dat}");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "does_not_exist.dat")) << joined(result);
  EXPECT_TRUE(hasError(result, "does not exist")) << joined(result);
}

TEST_F(ConfigTest, RejectsBothConstantAndFile) {
  std::string text = kBaseConfig;
  const std::string needle = "bottom_elevation: {constant: -1.0}";
  text.replace(text.find(needle), needle.size(),
               "bottom_elevation: {constant: -1.0, file: also.dat}");
  {
    std::ofstream touch(dir_ / "also.dat");
    touch << "0\n";
  }
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "exactly one of")) << joined(result);
}

TEST_F(ConfigTest, RejectsNoModulesEnabled) {
  std::string text = kBaseConfig;
  const std::string needle = "modules: {surface_water: true}";
  text.replace(text.find(needle), needle.size(), "modules: {surface_water: false}");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "at least one of surface_water/groundwater")) << joined(result);
}

TEST_F(ConfigTest, RejectsDensityCouplingWithoutTransport) {
  std::string text = kBaseConfig;
  const std::string needle = "modules: {surface_water: true}";
  text.replace(text.find(needle), needle.size(),
               "modules: {surface_water: true, groundwater: true}");
  text += R"(groundwater:
  timestep: {dt_init: 0.1, dt_min: 0.1, dt_max: 1.0}
  specific_storage: 1.0e-5
  density_coupling: {enabled: true}
soil:
  types:
    - {name: loam, ksx: 1.0e-5, ksy: 1.0e-5, ksz: 1.0e-5,
       theta_s: 0.4, theta_r: 0.08, vg_alpha: 6.0, vg_n: 2.0}
  map: {constant: loam}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "density coupling requires modules.transport")) << joined(result);
}

TEST_F(ConfigTest, RejectsKindTargetMismatch) {
  std::string text = kBaseConfig;
  text += R"(boundary_conditions:
  - name: bad
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: groundwater_top
    kind: eta
    value: {constant: 0.0}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "kind eta applies to target surface only")) << joined(result);
  EXPECT_TRUE(hasError(result, "requires modules.groundwater")) << joined(result);
}

TEST_F(ConfigTest, RejectsUndefinedSoilMapName) {
  std::string text = kBaseConfig;
  const std::string needle = "modules: {surface_water: true}";
  text.replace(text.find(needle), needle.size(),
               "modules: {surface_water: true, groundwater: true}");
  text += R"(groundwater:
  timestep: {dt_init: 0.1, dt_min: 0.1, dt_max: 1.0}
  specific_storage: 1.0e-5
soil:
  types:
    - {name: loam, ksx: 1.0e-5, ksy: 1.0e-5, ksz: 1.0e-5,
       theta_s: 0.4, theta_r: 0.08, vg_alpha: 6.0, vg_n: 2.0}
  map: {constant: sand}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "soil type 'sand' is not defined")) << joined(result);
}

TEST_F(ConfigTest, RejectsMonitorOutsideDomain) {
  std::string text = kBaseConfig;
  const std::string needle = "output: {filename: out/output.h5}";
  text.replace(text.find(needle), needle.size(),
               "output:\n  filename: out/output.h5\n  monitors:\n"
               "    - {name: probe, i: 99, j: 0, variables: [depth]}");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "outside the domain (nx = 4)")) << joined(result);
}

TEST_F(ConfigTest, RejectsMalformedYaml) {
  const ValidationResult result = validate("simulation: {id: broken\n  domain: [\n");
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "YAML parse error")) << joined(result);
}

TEST_F(ConfigTest, RejectsMissingConfigFile) {
  const ValidationResult result =
      frehg::validateConfigFile((dir_ / "absent.yaml").string());
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "cannot open configuration file")) << joined(result);
}

TEST_F(ConfigTest, RejectsRestartWithoutFile) {
  std::string text = kBaseConfig;
  text += "restart: {enabled: true}\n";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "restart.file: required when restart.enabled")) << joined(result);
}

TEST_F(ConfigTest, RejectsGravityValueOffTheBottom) {
  // P2 cross-check: legacy free drainage (bctype_GW code 3) is a
  // bottom-face behavior.
  std::string text = kBaseConfig;
  text += R"(boundary_conditions:
  - name: drain
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: groundwater_side
    kind: flux
    value: {gravity: true}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "free drainage applies to target groundwater_bottom only"))
      << joined(result);
}

TEST_F(ConfigTest, RejectsHydrostaticValueOffTheSides) {
  // P2 cross-check: the hydrostatic ghost-head form reproduces the legacy
  // side condition; top/bottom heads are prescribed pressure heads.
  std::string text = kBaseConfig;
  text += R"(boundary_conditions:
  - name: pond
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: groundwater_top
    kind: head
    value: {hydrostatic: {eta: 0.5}}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "hydrostatic values apply to target groundwater_side only"))
      << joined(result);
}

TEST_F(ConfigTest, RejectsUnknownReallocationSurplus) {
  std::string text = kBaseConfig;
  const std::string needle = "modules: {surface_water: true}";
  text.replace(text.find(needle), needle.size(),
               "modules: {surface_water: true, groundwater: true}");
  text += R"(groundwater:
  timestep: {dt_init: 0.1, dt_min: 0.1, dt_max: 1.0}
  specific_storage: 1.0e-5
  reallocation_surplus: conserve
soil:
  types:
    - {name: loam, ksx: 1.0e-5, ksy: 1.0e-5, ksz: 1.0e-5,
       theta_s: 0.4, theta_r: 0.08, vg_alpha: 6.0, vg_n: 2.0}
  map: {constant: loam}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "invalid value 'conserve'")) << joined(result);
}

TEST_F(ConfigTest, RejectsGravityValueOnNonFluxKind) {
  std::string text = kBaseConfig;
  text += R"(boundary_conditions:
  - name: tide
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: surface
    kind: eta
    value: {gravity: true}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "free drainage applies to kind flux only")) << joined(result);
}

TEST_F(ConfigTest, RejectsRainfallWithBothConstantAndSeries) {
  {
    std::ofstream series(dir_ / "rain.dat");
    series << "0 1e-6\n100 1e-6\n";
  }
  std::string text = kBaseConfig;
  const std::string needle = "  wetting_face_depth: 1.0e-6\n";
  text.insert(text.find(needle) + needle.size(),
              "  rainfall: {constant: 1.0e-6, series: {file: rain.dat}}\n");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "exactly one of {constant, series}")) << joined(result);
}

TEST_F(ConfigTest, RejectsRainfallExclusionWithTooFewVertices) {
  std::string text = kBaseConfig;
  const std::string needle = "  wetting_face_depth: 1.0e-6\n";
  text.insert(text.find(needle) + needle.size(),
              "  rainfall:\n    constant: 1.0e-6\n"
              "    exclude: {polygon: [[0.0, 0.0], [1.0, 1.0]]}\n");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "rainfall.exclude.polygon")) << joined(result);
}

TEST_F(ConfigTest, RejectsOutflowWithValue) {
  std::string text = kBaseConfig;
  text += R"(boundary_conditions:
  - name: outlet
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: surface
    kind: outflow
    value: {constant: 0.0}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "kind outflow takes no value")) << joined(result);
}

TEST_F(ConfigTest, RejectsMissingValueOnNonOutflowKind) {
  std::string text = kBaseConfig;
  text += R"(boundary_conditions:
  - name: tide
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: surface
    kind: eta
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "missing required key 'value'")) << joined(result);
}

TEST_F(ConfigTest, RejectsOutflowOnGroundwaterTarget) {
  std::string text = kBaseConfig;
  text += R"(boundary_conditions:
  - name: outlet
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: groundwater_side
    kind: outflow
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "kind outflow applies to target surface only")) << joined(result);
}

TEST_F(ConfigTest, RejectsSeepageOutputWithoutGroundwater) {
  std::string text = kBaseConfig;
  const std::string needle = "output: {filename: out/output.h5}";
  text.replace(text.find(needle), needle.size(),
               "output:\n  filename: out/output.h5\n"
               "  variables: {surface: [eta, seepage]}");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "'seepage' requires modules.groundwater")) << joined(result);
}

TEST_F(ConfigTest, LoadConfigOnInvalidFileIsFatal) {
  const std::string path = writeConfig("modules: {}\n");
  EXPECT_THROW(frehg::loadConfig(path), frehg::FatalError);
}

TEST_F(ConfigTest, RejectsUniformTerrainLayersWithoutFollowTerrain) {
  // P3 cross-check (amendment A12): the vertical-extent rule only selects
  // between the two terrain-following meshes.
  std::string text = kBaseConfig;
  const std::string needle = "bottom_elevation: {constant: -1.0}";
  text.replace(text.find(needle), needle.size(),
               "bottom_elevation: {constant: -1.0}\n  terrain_layers: uniform");
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "requires domain.follow_terrain: true")) << joined(result);
}

TEST_F(ConfigTest, AcceptsUniformTerrainLayersWithFollowTerrain) {
  std::string text = kBaseConfig;
  const std::string needle = "bottom_elevation: {constant: -1.0}";
  text.replace(text.find(needle), needle.size(),
               "bottom_elevation: {constant: -1.0}\n"
               "  follow_terrain: true\n  terrain_layers: uniform");
  const ValidationResult result = validate(text);
  EXPECT_TRUE(result.ok()) << joined(result);
  const frehg::FrehgConfig cfg = frehg::loadConfig(writeConfig(text));
  EXPECT_EQ(cfg.domain.terrainLayers, frehg::DomainConfig::TerrainLayers::Uniform);
}

TEST_F(ConfigTest, RejectsCoupledGroundwaterTopHead) {
  // P3 cross-check (amendment A11): the coupler owns the coupled top head;
  // only flux-kind conditions (the legacy qtop source) may be configured.
  std::string text = kBaseConfig;
  const std::string needle = "modules: {surface_water: true}";
  text.replace(text.find(needle), needle.size(),
               "modules: {surface_water: true, groundwater: true}");
  const std::string ic = "initial_conditions: {surface: {eta: {constant: 0.0}}}";
  text.replace(text.find(ic), ic.size(),
               "initial_conditions:\n"
               "  surface: {eta: {constant: 0.0}}\n"
               "  groundwater: {moisture: {constant: 0.2}}");
  text += R"(groundwater:
  timestep: {dt_init: 0.1, dt_min: 0.1, dt_max: 1.0}
  specific_storage: 1.0e-5
soil:
  types:
    - {name: loam, ksx: 1.0e-5, ksy: 1.0e-5, ksz: 1.0e-5,
       theta_s: 0.4, theta_r: 0.08, vg_alpha: 6.0, vg_n: 2.0}
  map: {constant: loam}
boundary_conditions:
  - name: pond
    region: {polygon: [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]}
    target: groundwater_top
    kind: head
    value: {constant: 0.1}
)";
  const ValidationResult result = validate(text);
  ASSERT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "must be kind flux in coupled runs")) << joined(result);
  // The flux kind (the legacy qtop source) stays accepted.
  const std::string headKind = "kind: head\n    value: {constant: 0.1}";
  text.replace(text.find(headKind), headKind.size(),
               "kind: flux\n    value: {constant: 0.1}");
  EXPECT_TRUE(validate(text).ok()) << joined(validate(text));
}

// v2 Q1 (plan §2.2): the solver block parses per-system preconditioner
// selection and reuse policy, defaults reproduce v1 (bjacobi-icc, gw cap
// 1000), and the schema rejects unknown preconditioners.
TEST_F(ConfigTest, SolverDefaultsReproduceV1) {
  const frehg::FrehgConfig cfg = frehg::loadConfig(writeConfig(kBaseConfig));
  EXPECT_EQ(cfg.solver.surface.preconditioner, "bjacobi-icc");
  EXPECT_EQ(cfg.solver.groundwater.preconditioner, "bjacobi-icc");
  EXPECT_DOUBLE_EQ(cfg.solver.surface.rtol, 1.0e-8);
  EXPECT_DOUBLE_EQ(cfg.solver.groundwater.atol, 1.0e-14);
  EXPECT_EQ(cfg.solver.surface.maxIterations, 500);
  EXPECT_EQ(cfg.solver.groundwater.maxIterations, 1000);
  EXPECT_EQ(cfg.solver.surface.reuseMaxSolves, 50);
  EXPECT_DOUBLE_EQ(cfg.solver.groundwater.reuseIterationFactor, 1.5);
}

TEST_F(ConfigTest, SolverBlockParsesSelections) {
  std::string text = kBaseConfig;
  text += R"(solver:
  surface: {preconditioner: amg, rtol: 1.0e-10, reuse_max_solves: 20}
  groundwater: {preconditioner: gamg, max_iterations: 2000,
                reuse_iteration_factor: 2.0}
)";
  const ValidationResult result = validate(text);
  EXPECT_TRUE(result.ok()) << joined(result);
  const frehg::FrehgConfig cfg = frehg::loadConfig(writeConfig(text));
  EXPECT_EQ(cfg.solver.surface.preconditioner, "amg");
  EXPECT_DOUBLE_EQ(cfg.solver.surface.rtol, 1.0e-10);
  EXPECT_EQ(cfg.solver.surface.reuseMaxSolves, 20);
  EXPECT_EQ(cfg.solver.groundwater.preconditioner, "gamg");
  EXPECT_EQ(cfg.solver.groundwater.maxIterations, 2000);
  EXPECT_DOUBLE_EQ(cfg.solver.groundwater.reuseIterationFactor, 2.0);
}

TEST_F(ConfigTest, SolverRejectsUnknownPreconditioner) {
  std::string text = kBaseConfig;
  text += "solver: {surface: {preconditioner: ilu}}\n";
  const ValidationResult result = validate(text);
  EXPECT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "preconditioner")) << joined(result);
}

TEST_F(ConfigTest, SolverRejectsUnknownKeys) {
  std::string text = kBaseConfig;
  text += "solver: {surface: {krylov: gmres}}\n";
  const ValidationResult result = validate(text);
  EXPECT_FALSE(result.ok());
  EXPECT_TRUE(hasError(result, "krylov")) << joined(result);
}

}  // namespace
