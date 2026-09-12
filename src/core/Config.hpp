/// \file Config.hpp
/// \brief Typed Frehg2 configuration (YAML v2 schema, plan §5.4 and §6).
///
/// The configuration pipeline has two stages:
///  1. structural validation against the declarative schema in
///     ConfigSchema.cpp — type/range/allowed-value checks, hard unknown-key
///     rejection with nearest-key suggestions, and cross-field checks;
///  2. extraction into the plain structs below, with every defaulted key
///     materialized so the effective configuration can be echoed to the log
///     and embedded into the HDF5 output.
///
/// `frehg --validate config.yaml` runs stage 1 only and reports all errors.

#ifndef FREHG_CORE_CONFIG_HPP
#define FREHG_CORE_CONFIG_HPP

#include "core/Types.hpp"

#include <array>
#include <string>
#include <vector>

namespace YAML {
class Node;  // forward declaration; yaml-cpp is an implementation detail
}

namespace frehg {

/// A scalar quantity given either as a constant or as a per-cell input file.
struct FileOrConstant {
  bool fromFile = false;   ///< true if \c file names the source
  std::string file;        ///< path relative to the configuration file
  real_t constant = 0.0;   ///< used when \c fromFile is false
};

/// A time-dependent quantity given either as a constant or as a series file.
struct SeriesOrConstant {
  bool fromSeries = false;  ///< true if \c file names a time-series file
  std::string file;         ///< path relative to the configuration file
  real_t constant = 0.0;    ///< used when \c fromSeries is false
};

/// simulation: run identification.
struct SimulationConfig {
  std::string id;     ///< short run identifier
  std::string title;  ///< free-text description (optional)
};

/// domain.decomposition: MPI block layout; 0 means "auto" (MPI_Dims_create).
struct DecompositionConfig {
  int mpiNx = 0;  ///< ranks along i, 0 = auto
  int mpiNy = 0;  ///< ranks along j, 0 = auto
};

/// domain: global grid extents, spacings, and bathymetry source.
struct DomainConfig {
  /// Vertical extent rule of the terrain-following mesh (plan amendment
  /// A12). Scaled: every column's nz layers are scaled to span from the
  /// local bed to the common box bottom (the legacy follow_terrain wedge,
  /// map.c:316-353 — b6's tank). Uniform: every column carries the
  /// configured dz * dz_stretch^k layer profile below its own bed (the
  /// terrain-parallel slab of the b5 reference, SERGHEI GwInit.h:109-117).
  enum class TerrainLayers { Scaled, Uniform };
  int nx = 0;                      ///< global cells in i
  int ny = 0;                      ///< global cells in j
  int nz = 0;                      ///< global cells in k (subsurface layers)
  real_t dx = 0.0;                 ///< cell size in i [m]
  real_t dy = 0.0;                 ///< cell size in j [m]
  real_t dz = 0.0;                 ///< top-layer thickness [m]
  real_t dzStretch = 1.0;          ///< geometric layer growth factor (legacy dz_incre)
  FileOrConstant bottomElevation;  ///< bed elevation [m]
  bool followTerrain = false;      ///< terrain-following subsurface mesh
  TerrainLayers terrainLayers = TerrainLayers::Scaled;  ///< vertical extent rule
  DecompositionConfig decomposition;  ///< MPI block layout request
};

/// time: stepping and output cadence.
struct TimeConfig {
  /// Surface time step [s]. Fixed in single-module and subcycled-coupling
  /// runs; the initial value of the common adaptive step in sync-coupled
  /// runs (legacy solve.c:37,193 — plan amendment A10).
  real_t dt = 0.0;
  real_t tStart = 0.0;          ///< start time [s], integer-valued (§7)
  real_t tEnd = 0.0;            ///< end time [s]
  real_t outputInterval = 0.0;  ///< spatial output cadence [s], integer-valued (§7)
};

/// modules: which physics modules run.
struct ModulesConfig {
  bool surfaceWater = false;  ///< run the SWE module
  bool groundwater = false;   ///< run the Richards module
  bool transport = false;     ///< run scalar transport
};

/// surface_water.friction: drag law and roughness (plan §5.8).
struct FrictionConfig {
  /// Available drag laws; Chezy is the plan §5.8 additive extension.
  enum class Law { Manning, Chezy };
  Law law = Law::Manning;  ///< selected drag law
  FileOrConstant coefficient;    ///< Manning n or Chezy C
  real_t thinLayerDepth = 0.1;   ///< legacy hD: exponent-switch depth [m]
};

/// surface_water.wind: quadratic wind stress with thin-layer attenuation.
struct WindConfig {
  bool enabled = false;            ///< apply wind stress
  real_t cd = 0.0013;              ///< legacy Cw drag coefficient
  real_t attenuationDepth = 5.0;   ///< legacy CwT thin-layer attenuation depth [m]
  real_t northAngle = 0.0;         ///< grid-to-north rotation [deg]
  SeriesOrConstant speed;          ///< wind speed [m/s]
  SeriesOrConstant direction;      ///< wind direction [deg]
};

/// surface_water: SWE module parameters.
struct SurfaceWaterConfig {
  real_t gravity = 9.81;           ///< [m/s^2]
  FrictionConfig friction;        ///< drag law and roughness
  real_t viscosityX = 1.0e-6;      ///< eddy viscosity in i [m^2/s]
  real_t viscosityY = 1.0e-6;      ///< eddy viscosity in j [m^2/s]
  real_t minDepth = 0.0;           ///< wet/dry threshold [m] (legacy min_dept)
  real_t wettingFaceDepth = 0.0;   ///< face-wetting threshold [m] (legacy wtfh)
  WindConfig wind;                 ///< wind stress forcing
  SeriesOrConstant rainfall;       ///< rain rate [m/s]
  /// Region receiving no rainfall (empty = rain everywhere). Expresses the
  /// legacy hardcoded skip of the last global row (shallowwater.c:596),
  /// which the b1 goldens embed for their outlet row.
  std::vector<std::array<real_t, 2>> rainfallExcludePolygon;
  SeriesOrConstant evaporation;    ///< evaporation rate [m/s]
};

/// groundwater.timestep: adaptive dtg controller (legacy semantics, §3.1).
struct GroundwaterTimestepConfig {
  real_t dtInit = 0.0;      ///< initial subsurface step [s]
  real_t dtMin = 0.0;       ///< lower clamp [s]
  real_t dtMax = 0.0;       ///< upper clamp [s]
  real_t dqGrow = 0.01;     ///< max flux change below which dtg grows
  real_t dqShrink = 0.02;   ///< max flux change above which dtg shrinks
  real_t courantMax = 2.0;  ///< legacy Co_max: Courant limit on dK/dtheta
};

/// groundwater.density_coupling: baroclinic feedback (activated in P4).
struct DensityCouplingConfig {
  bool enabled = false;  ///< activate r_rho / r_visc face ratios (P4)
};

/// groundwater: Richards module parameters.
struct GroundwaterConfig {
  /// How the post-allocation step handles the surplus water of unsaturated
  /// cells adjacent to saturation, θ_corrector - θ(h). The legacy sweep
  /// computes the surplus and its gradient split but leaves the transfer
  /// disabled (groundwater.c:1020-1027), discarding the volume — the b2
  /// golden embeds that behavior. Redistribute delivers the vertical
  /// fractions into the column's pore room instead (mass-conserving; the
  /// b3 physics gate requires it). Plan amendment A7.
  enum class ReallocationSurplus { Drop, Redistribute };
  std::string scheme = "pca";     ///< sole allowed value; documents the design
  bool useFull3d = true;          ///< false: zero lateral K in unsaturated cells
  GroundwaterTimestepConfig timestep;  ///< adaptive dtg controller
  real_t specificStorage = 0.0;   ///< Ss [1/m]
  /// Surplus handling in the post-allocation step (legacy default: drop).
  ReallocationSurplus reallocationSurplus = ReallocationSurplus::Drop;
  DensityCouplingConfig densityCoupling;  ///< baroclinic feedback switch
};

/// soil.types[]: van Genuchten–Mualem soil parameters.
struct SoilType {
  std::string name;      ///< unique type name referenced by soil.map
  real_t ksx = 0.0;      ///< saturated conductivity in i [m/s]
  real_t ksy = 0.0;      ///< saturated conductivity in j [m/s]
  real_t ksz = 0.0;      ///< saturated conductivity in k [m/s]
  real_t thetaS = 0.0;   ///< saturated water content
  real_t thetaR = 0.0;   ///< residual water content
  real_t vgAlpha = 0.0;  ///< van Genuchten alpha [1/m]
  real_t vgN = 0.0;      ///< van Genuchten n
  real_t aev = 0.0;      ///< saturation-cutoff (air-entry) head [m], <= 0
};

/// soil.map: soil-type assignment, by uniform type name or id raster.
struct SoilMapConfig {
  bool fromFile = false;     ///< true if \c file names the id raster
  std::string file;          ///< integer ids indexing soil.types (0-based)
  std::string constantName;  ///< uniform soil type by name
};

/// soil: soil-type table and spatial map.
struct SoilConfig {
  std::vector<SoilType> types;  ///< the soil-type table
  SoilMapConfig map;            ///< spatial assignment of types
};

/// coupling: surface–subsurface synchronization mode (plan §5.7, §10 P3).
struct CouplingConfig {
  /// Lockstep (dtg = dt) or adaptive-dtg subcycling inside the surface step.
  enum class Mode { Sync, Subcycled };
  Mode mode = Mode::Sync;  ///< synchronization mode
};

/// initial_conditions.surface
struct SurfaceInitialConfig {
  FileOrConstant eta;          ///< free-surface elevation [m]
  bool hasUu = false;          ///< true when an initial u-velocity is given
  FileOrConstant uu;           ///< initial u-velocity [m/s]
  bool hasVv = false;          ///< true when an initial v-velocity is given
  FileOrConstant vv;           ///< initial v-velocity [m/s]
};

/// initial_conditions.groundwater — exactly one of the three forms.
struct GroundwaterInitialConfig {
  /// Which of the three exclusive forms the configuration used.
  enum class Form { WaterTable, Head, Moisture };
  Form form = Form::Head;  ///< selected form
  FileOrConstant value;    ///< the initial field for that form
};

/// initial_conditions.transport
struct TransportInitialConfig {
  FileOrConstant surface;      ///< initial surface concentration
  FileOrConstant groundwater;  ///< initial subsurface concentration
};

/// initial_conditions
struct InitialConditionsConfig {
  SurfaceInitialConfig surface;          ///< SWE initial state
  GroundwaterInitialConfig groundwater;  ///< Richards initial state
  TransportInitialConfig transport;      ///< scalar initial state
};

/// boundary_conditions[].target: which sub-boundary the condition applies to.
enum class BcTarget { Surface, GroundwaterTop, GroundwaterBottom, GroundwaterSide };

/// boundary_conditions[].kind (plan §5.6): covers every legacy behavior.
/// Outflow is the free (transmissive) outflow decided by the P1 b4 gate
/// (plan §5.6 amendment): the prescribed-stage sink mapping kept the outlet
/// column dry, so SERGHEI's bctype-9 free outflow gets its own kind.
enum class BcKind { Eta, Discharge, Velocity, Outflow, Head, Flux, ScalarValue };

/// boundary_conditions[].value — exactly one form.
struct BcValueConfig {
  /// The exclusive value forms a boundary condition may take.
  enum class Form {
    Constant,     ///< fixed value
    Series,       ///< time series file
    Gravity,      ///< free drainage (kind flux only; legacy bctype_GW code 3)
    Hydrostatic   ///< head = eta - z (kind head only)
  };
  Form form = Form::Constant;   ///< selected value form
  real_t constant = 0.0;        ///< value for Form::Constant
  std::string seriesFile;       ///< series path for Form::Series
  real_t hydrostaticEta = 0.0;  ///< reference stage for Form::Hydrostatic
};

/// boundary_conditions[]: one polygon-region boundary condition.
struct BoundaryConditionConfig {
  std::string name;   ///< unique condition name
  std::vector<std::array<real_t, 2>> polygon;  ///< region vertices (x, y) [m]
  BcTarget target = BcTarget::Surface;  ///< sub-boundary the condition acts on
  BcKind kind = BcKind::Eta;            ///< what quantity is prescribed
  BcValueConfig value;                  ///< the prescribed value
};

/// transport.scheme
struct TransportSchemeConfig {
  /// Advection scheme: first-order upwind or TVD superbee (plan §3.1).
  enum class Advection { Upwind, Superbee };
  Advection advection = Advection::Upwind;  ///< selected scheme
};

/// transport: scalar transport parameters.
struct TransportConfig {
  TransportSchemeConfig scheme;          ///< advection scheme selection
  real_t surfaceDiffusivityX = 1.0e-10;  ///< [m^2/s]
  real_t surfaceDiffusivityY = 1.0e-10;  ///< [m^2/s]
  real_t dispersionLongitudinal = 0.0;   ///< [m]
  real_t dispersionTransverse = 0.0;     ///< [m]
  real_t dispersionMolecular = 1.0e-10;  ///< [m^2/s]
  real_t boundMin = 0.0;                 ///< lower scalar bound (plan §3.2)
  bool hasBoundMax = false;              ///< true if an upper bound is set
  real_t boundMax = 0.0;                 ///< upper scalar bound when present
};

/// output.monitors[]: point-probe time series.
struct MonitorConfig {
  std::string name;  ///< unique monitor name (HDF5 table name)
  int i = 0;  ///< global cell index in i
  int j = 0;  ///< global cell index in j
  std::vector<std::string> variables;  ///< variables sampled at the point
};

/// output: HDF5 file, variable selection, monitors, checkpoint cadence.
struct OutputConfig {
  std::string filename;  ///< single HDF5 output file (plan §7)
  std::vector<std::string> surfaceVariables;      ///< /surface/* selection
  std::vector<std::string> groundwaterVariables;  ///< /groundwater/* selection
  std::vector<std::string> transportVariables;    ///< /transport/* selection
  std::vector<MonitorConfig> monitors;            ///< point probes
  real_t checkpointInterval = 0.0;  ///< 0 = off; final checkpoint always written
};

/// restart
struct RestartConfig {
  bool enabled = false;  ///< resume from a checkpoint
  std::string file;      ///< checkpointed output file
  real_t time = 0.0;     ///< checkpoint time to resume from [s]
};

/// solver.\<system\>: per-system linear-solver selection (v2 plan §2.2).
/// Defaults reproduce the v1 hardcoded behavior (CG + block-Jacobi/ICC(0));
/// `amg` selects hypre BoomerAMG, `gamg` PETSc's built-in smoothed
/// aggregation. Every PETSc detail stays overridable through the options
/// database (the constructor calls KSPSetFromOptions last).
struct SolverSystemConfig {
  std::string preconditioner = "bjacobi-icc";  ///< bjacobi-icc | amg | gamg
  /// PETSc matrix/vector backend (v2 plan §2B.2 B1): "aij" (host, the
  /// default) or "aijkokkos" (Kokkos Kernels on the build's execution
  /// space — the threaded-CPU and GPU solve path; requires a Kokkos-enabled
  /// PETSc). Device builds force "aijkokkos" regardless of this value.
  std::string matType = "aij";
  real_t rtol = 1.0e-8;    ///< relative tolerance (legacy SetRTCAccuracy)
  real_t atol = 1.0e-14;   ///< absolute tolerance
  int maxIterations = 500; ///< iteration cap
  /// AMG hierarchy reuse cadence: rebuild the preconditioner at least every
  /// this many solves (0 = rebuild every solve, the bjacobi-icc behavior).
  /// Ignored for bjacobi-icc, whose per-solve refresh is cheap.
  int reuseMaxSolves = 50;
  /// Early-rebuild trigger: rebuild when an iteration count exceeds this
  /// factor times the count measured right after the last rebuild.
  real_t reuseIterationFactor = 1.5;
};

/// solver
struct SolverConfig {
  SolverConfig() { groundwater.maxIterations = 1000; }
  SolverSystemConfig surface;      ///< the fs_ free-surface system
  SolverSystemConfig groundwater;  ///< the gw_ Richards predictor system
  std::string petscOptionsFile;    ///< optional PETSc options file
};

/// runtime
struct RuntimeConfig {
  /// Whether MPI receives device pointers directly (plan §5.3).
  enum class GpuAwareMpi { Auto, On, Off };
  GpuAwareMpi gpuAwareMpi = GpuAwareMpi::Auto;  ///< staging mode
};

/// The complete, materialized Frehg2 configuration.
struct FrehgConfig {
  SimulationConfig simulation;    ///< run identification
  DomainConfig domain;            ///< grid extents and bathymetry
  TimeConfig time;                ///< stepping and cadence
  ModulesConfig modules;          ///< enabled physics
  SurfaceWaterConfig surfaceWater;///< SWE parameters
  GroundwaterConfig groundwater;  ///< Richards parameters
  SoilConfig soil;                ///< soil table and map
  CouplingConfig coupling;        ///< coupling mode
  InitialConditionsConfig initialConditions;  ///< initial state
  std::vector<BoundaryConditionConfig> boundaryConditions;  ///< BC list
  TransportConfig transport;      ///< scalar-transport parameters
  OutputConfig output;            ///< output selection
  RestartConfig restart;          ///< restart request
  SolverConfig solver;            ///< PETSc options hookup
  RuntimeConfig runtime;          ///< runtime toggles

  std::string configDir;  ///< directory of the YAML file (for relative paths)
  std::string rawText;    ///< full configuration text (embedded in HDF5 output)

  /// Resolve a config-relative path against configDir.
  std::string resolvePath(const std::string& relative) const;
};

/// Outcome of configuration validation.
struct ValidationResult {
  std::vector<std::string> errors;  ///< empty when the configuration is valid
  /// \return true when no errors were recorded.
  bool ok() const { return errors.empty(); }
};

/// Validate \p path against the v2 schema, including cross-field checks and
/// referenced-file existence. Never throws for configuration mistakes; all
/// problems are collected into the result. I/O failure to read the file
/// itself is also reported as an error entry.
ValidationResult validateConfigFile(const std::string& path);

/// Load, validate, and materialize a configuration. Fatal (throws
/// frehg::FatalError via the logger) if validation reports any error.
FrehgConfig loadConfig(const std::string& path);

/// Render the effective (defaults-materialized) configuration as a
/// human-readable listing for the log header (plan §5.4).
std::string describeConfig(const FrehgConfig& config);

/// Serialize the resolved configuration back to schema-valid YAML (v2 plan
/// §2A). This is the exact inverse of extraction for every materialized
/// field: `loadConfig(resolvedConfigYaml(loadConfig(f)))` is a fixed point,
/// and `frehg --resolve f` prints it so the run record's embedded
/// configuration can be verified byte-for-byte against a re-resolve of the
/// original input (gate r1). Blocks are emitted per enabled module
/// (mirroring the extraction guards), so the output revalidates.
std::string resolvedConfigYaml(const FrehgConfig& config);

namespace detail {
/// Structural + cross-field validation used by both entry points above.
/// \p configDir anchors relative file references. Declared here so
/// ConfigSchema.cpp and Config.cpp share it without a public yaml-cpp
/// dependency.
void validateRoot(const YAML::Node& root, const std::string& configDir,
                  std::vector<std::string>& errors);
}  // namespace detail

}  // namespace frehg

#endif  // FREHG_CORE_CONFIG_HPP
