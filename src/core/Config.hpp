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
  /// Stage-3 T3 (plan §13.B): periodic lateral topology per axis — a
  /// domain property for EVERY module (the surface momentum/continuity,
  /// the groundwater matrix and fluxes, transport, baroclinic faces).
  bool periodicX = false;  ///< x axis wraps (seam is an interior face)
  bool periodicY = false;  ///< y axis wraps (seam is an interior face)
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
  bool transport = false;     ///< run scalar (salinity) transport
  bool temperature = false;   ///< run temperature transport (v2 Q5)
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
  /// Selectable Cd(U10) laws (v2 Q6, plan §5.2). Constant is the legacy
  /// Cw (no cap); the published laws are capped at `cap`.
  enum class DragLaw { Constant, Garratt, SmithBanke, Wu, LargePond };
  bool enabled = false;            ///< apply wind stress
  DragLaw law = DragLaw::Constant; ///< wind.law
  real_t cd = 0.0013;              ///< legacy Cw drag coefficient (constant law)
  real_t cap = 3.5e-3;             ///< Cd cap for the U10 laws (wind.cap)
  real_t attenuationDepth = 5.0;   ///< legacy CwT thin-layer attenuation depth [m]
  real_t northAngle = 0.0;         ///< grid-to-north rotation [deg]
  /// Input form (v2 Q6): grid-frame components (u10, v10) — the wrap-free
  /// form — or the legacy compass speed/direction pair.
  bool componentForm = false;
  SeriesOrConstant speed;          ///< wind speed [m/s]
  SeriesOrConstant direction;      ///< wind direction [deg]
  SeriesOrConstant u10;            ///< wind +x component [m/s]
  SeriesOrConstant v10;            ///< wind +y component [m/s]
};

/// atmosphere: met forcing for the bulk-aerodynamic module (v2 Q4, plan
/// §3.2). One vapor-physics module, several consumers (open-water and soil
/// evaporation now; the Q5 surface heat exchange later). All previously
/// hardcoded legacy constants are required inputs here — that is what
/// makes the legacy aerodynamic formulation shippable again (v1 removed
/// it for its hardcoded 20 C air/pressure/humidity, not its physics).
struct AtmosphereConfig {
  bool present = false;                 ///< true when the block is configured
  SeriesOrConstant airTemperature;      ///< [C]
  /// T_s [C] for the Q4 evaporation consumers (surface bulk / soil bulk).
  /// Optional since v2 Q5: the surface HEAT exchange never reads it (it
  /// evaluates the bulk chain at the local water temperature, V2-A17);
  /// the schema requires it exactly when a Q4 consumer is configured.
  bool hasSurfaceTemperature = false;   ///< true when the key was given
  SeriesOrConstant surfaceTemperature;  ///< T_s [C]
  SeriesOrConstant pressure;            ///< P0 [kPa]
  /// Exactly one humidity form (schema cross-check).
  bool humidityIsRelative = false;      ///< which form was configured
  SeriesOrConstant specificHumidity;    ///< q_a [-]
  SeriesOrConstant relativeHumidity;    ///< [-] of q_sat(air temperature)
  SeriesOrConstant windSpeed;           ///< U [m/s] (shared with the Q6 wind)
  /// v2 Q5 (plan V2-A17): prescribed absorbed shortwave and incident
  /// longwave radiation [W/m^2] for the bulk surface heat exchange
  /// (radiation SCHEMES stay out of scope — plan §10), and the GLM
  /// still-air lower bound on the wind speed entering the bulk transfer
  /// functions (applied at MetForcing::sample for every bulk consumer).
  SeriesOrConstant shortwave;           ///< absorbed SW [W/m^2] (default 0)
  SeriesOrConstant longwaveIn;          ///< incident LW [W/m^2] (default 0)
  /// Stage-3 T7 (plan §11 W7.1): INCIDENT shortwave R_s [W/m^2] for the
  /// Penman-Monteith potential, whose channels apply their own albedo.
  /// Distinct from the Q5 `shortwave`, which stays the ABSORBED SW of
  /// the bulk surface heat exchange (its meaning is unchanged).
  SeriesOrConstant shortwaveIncoming;   ///< incident SW [W/m^2] (default 0)
  real_t windSpeedFloor = 0.5;          ///< [m/s] still-air lower bound
  /// Stage-3 T5 (plan §9 W5.1): the aerodynamic-resistance law of the ET
  /// consumers. The Q4/Q5 bulk consumers (surface bulk evaporation, the
  /// heat module's latent/sensible terms, groundwater.evaporation) stay
  /// pinned to liu2006 regardless (their goldens).
  enum class AeroResistance { Liu2006, Fao56, LogProfile };
  AeroResistance aeroResistance = AeroResistance::Liu2006;  ///< selected r_a law
  real_t aeroRefHeight = 2.0;     ///< measurement height z_ref [m]
  real_t aeroZ0m = 0.01;          ///< momentum roughness z0m [m] (log_profile)
  real_t aeroZ0h = 0.001;         ///< scalar roughness z0h [m] (log_profile)
  real_t aeroDisplacement = 0.0;  ///< displacement height d [m] (log_profile)
};

/// evapotranspiration.energy (Stage-3 T7, plan §11): the Penman-Monteith
/// energy-budget settings. The available energy is
///   A = R_n - G,  R_n = (1 - albedo) R_s + eps (LW_in - sigma T_K^4),
///   G = fraction * R_n,
/// with eps = kWaterEmissivity and T_K = T_air + 273.15. The albedo and
/// ground-heat fraction are per channel: water uses albedo_water and the
/// water fraction; the soil and canopy channels take their albedo from
/// the soil / vegetation type and their fraction from here.
struct EnergyConfig {
  real_t albedoWater = 0.06;     ///< open-water albedo [-]
  real_t gFractionWater = 0.0;   ///< ground-heat fraction, water channel
  real_t gFractionSoil = 0.1;    ///< ground-heat fraction, soil channel
  real_t gFractionCanopy = 0.1;  ///< ground-heat fraction, canopy channel
};

/// evapotranspiration: the T5 ET block (plan §9) — replaces BOTH legacy
/// evaporation keys (surface_water.evaporation / groundwater.evaporation;
/// the schema rejects mixing the forms). Open-water and soil channels
/// share one potential mode; complementarity holds on a time basis by
/// construction (theory/evapotranspiration.md).
struct EvapotranspirationConfig {
  bool enabled = false;           ///< true when the block is configured
  /// The shared potential mode of both channels (PM parses in T7).
  enum class Potential { Bulk, PenmanMonteith, Prescribed };
  Potential potential = Potential::Bulk;  ///< selected mode
  SeriesOrConstant prescribed;    ///< E_pot [m/s] (potential prescribed)
  bool openWater = true;          ///< the surface channel
  bool soil = false;              ///< the soil channel (block present)
  /// Soil-channel region (empty = every column).
  std::vector<std::array<real_t, 2>> soilPolygon;
  /// Soil-rate limiter (alpha_barton acts inside the bulk humidity
  /// difference; beta_lee_pielke multiplies; none passes the potential).
  enum class Limiter { AlphaBarton, BetaLeePielke, None };
  Limiter limiter = Limiter::AlphaBarton;  ///< selected limiter
  EnergyConfig energy;            ///< PM energy budget (T7)
};

/// surface_water: SWE module parameters.
struct SurfaceWaterConfig {
  /// Open-water evaporation mode (v2 Q4 §3.2). Prescribed applies the
  /// configured rate unconditionally with the legacy dry clamp (b1 golden
  /// fidelity; optional exclude region masks cells, the rain symmetry
  /// fix). Bulk computes the rate per step from the atmosphere block
  /// (water-surface q_g = q_sat(T_s)) and removes at most the available
  /// depth per wet cell — no volume creation.
  enum class EvapMode { Prescribed, Bulk };
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
  SeriesOrConstant evaporation;    ///< evaporation rate [m/s] (prescribed)
  EvapMode evapMode = EvapMode::Prescribed;  ///< evaporation.mode
  /// Region receiving no evaporation (prescribed mode; empty = everywhere).
  std::vector<std::array<real_t, 2>> evaporationExcludePolygon;
};

/// groundwater.evaporation: bulk-aerodynamic soil evaporation (v2 Q4
/// §3.2). Potential rate from the atmosphere block; actual rate limited
/// through the surface-layer soil relative humidity alpha_1 (Geng &
/// Boufadel Eq. 6), applied as the top-face flux over the region.
/// Uncoupled groundwater runs only in v2.0 (schema-enforced): the coupled
/// top exchange is owned by the coupler and no gate exercises the coupled
/// combination yet.
struct GwEvaporationConfig {
  bool enabled = false;                        ///< bulk mode configured
  std::vector<std::array<real_t, 2>> polygon;  ///< evaporation zone (x, y)
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

/// groundwater.density_coupling: baroclinic feedback (activated in P4;
/// coefficients made configuration in v2 Q5 with the legacy compile-time
/// constants as defaults — golden-neutral, V2-A17):
///   r_rho  = 1 + beta_saline s - thermal_expansion (T - reference_temperature)
///   r_visc = 1 / (1 + beta_saline_viscosity s)   (mu(T) out of scope)
struct DensityCouplingConfig {
  bool enabled = false;  ///< activate r_rho / r_visc face ratios (P4)
  real_t betaSaline = 0.000744;          ///< legacy r_rho slope [L/g]
  real_t betaSalineViscosity = 0.0022;   ///< legacy r_visc slope [L/g]
  real_t thermalExpansion = 0.0;         ///< beta_T [1/K] (Q5; 0 = off)
  real_t referenceTemperature = 20.0;    ///< T0 [C] for the thermal term
  /// groundwater.density_coupling.reallocation_gradient: density (post-v2.0
  /// opt-in). The post-allocation head-gradient split (legacy
  /// check_head_gradient, groundwater.c:1074-1150) measures the vertical
  /// gradient against freshwater gravity, dh/dz - 1. Under density coupling
  /// a hydrostatic saline column has dh/dz = r_rho, so the legacy test sees
  /// a spurious upward gradient and routes restore surpluses upward. The
  /// density form subtracts the face density ratio instead — the gravity
  /// term of the Darcy flux itself. Off by default (`freshwater`, legacy).
  bool densityReallocationGradient = false;
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
  /// groundwater.reallocation (Stage-3 T1, owner decision O9; S3-A4).
  /// Legacy: the full v2.0 post-allocation sweep (sends, theta := theta(h)
  /// restores, drops). Off (YAML token "disabled"): the same sweep minus
  /// every non-conservative action — no sends, no drops, no theta :=
  /// theta(h) creation; the conservative isolated-cell h := h(theta)
  /// restore is kept, and the realloc audit columns are exactly zero.
  enum class Reallocation { Legacy, Off };
  std::string scheme = "pca";     ///< sole allowed value; documents the design
  bool useFull3d = true;          ///< false: zero lateral K in unsaturated cells
  GroundwaterTimestepConfig timestep;  ///< adaptive dtg controller
  real_t specificStorage = 0.0;   ///< Ss [1/m]
  /// Surplus handling in the post-allocation step (legacy default: drop).
  ReallocationSurplus reallocationSurplus = ReallocationSurplus::Drop;
  /// THE STAGE-3 T8 DEFAULT-FLIP CONSTANT (plan §2.7): the T8 evaluation
  /// changes exactly this initializer to Reallocation::Off and re-runs
  /// every test with unchanged tolerances; it stays flipped only if all
  /// pass. Do not flip it anywhere else.
  Reallocation reallocation = Reallocation::Legacy;
  DensityCouplingConfig densityCoupling;  ///< baroclinic feedback switch
  GwEvaporationConfig evaporation;        ///< bulk soil evaporation (v2 Q4)
};

/// transport.scheme / temperature.scheme
struct TransportSchemeConfig {
  /// Advection scheme: first-order upwind or TVD superbee (plan §3.1).
  enum class Advection { Upwind, Superbee };
  Advection advection = Advection::Upwind;  ///< selected scheme
};

/// temperature.surface_exchange: the Q5 surface heat-flux term.
struct SurfaceExchangeConfig {
  /// None: no surface heat exchange. Equilibrium: dT/dt =
  /// -K_e (T - T_e)/((rho c)_w h) (Edinger). Bulk: Q_net = Q_sw +
  /// eps (LW_in - sigma T_K^4) - Q_lat - Q_sens on the shared Q4 atm
  /// module chain, evaluated at the local water temperature.
  enum class Mode { None, Equilibrium, Bulk };
  Mode mode = Mode::None;         ///< selected mode
  real_t equilibriumTemperature = 0.0;  ///< T_e [C] (equilibrium mode)
  real_t equilibriumCoefficient = 0.0;  ///< K_e [W/m^2/K] (equilibrium mode)
};

/// temperature: the second registered scalar (v2 Q5, plan §4.1/V2-A17).
/// The transport machinery is shared with salinity; these are the
/// per-scalar parameters and the thermal physics switches.
struct TemperatureConfig {
  TransportSchemeConfig scheme;          ///< advection scheme selection
  real_t surfaceDiffusivityX = 1.0e-10;  ///< [m^2/s]
  real_t surfaceDiffusivityY = 1.0e-10;  ///< [m^2/s]
  /// Subsurface effective thermal conduction lambda_eff [W/m/K]; enters
  /// the dispersion tensor's molecular slot as lambda/(rho c)_w (a bulk
  /// property — NOT multiplied by theta_s, unlike solute molecular
  /// diffusion). Required with the groundwater module.
  real_t thermalConductivity = 0.0;
  real_t heatCapacityWater = 4.184e6;    ///< (rho c)_w [J/m^3/K]
  /// Solid volumetric heat capacity (rho c)_s [J/m^3/K]: the thermal
  /// retardation kappa = (1 - theta_s) (rho c)_s / (rho c)_w joins the
  /// subsurface mass basis (theta + kappa) — the SEAWAT
  /// retardation-as-scalar form. Required with the groundwater module.
  real_t heatCapacitySolid = 0.0;
  real_t dispersivityLongitudinal = 0.0; ///< thermal dispersivity [m]
  real_t dispersivityTransverse = 0.0;   ///< thermal dispersivity [m]
  /// Optional temperature bounds; both default OPEN (unlike salinity,
  /// whose lower bound defaults to 0).
  bool hasBoundMin = false;              ///< true if a lower bound is set
  real_t boundMin = 0.0;                 ///< lower bound when present
  bool hasBoundMax = false;              ///< true if an upper bound is set
  real_t boundMax = 0.0;                 ///< upper bound when present
  SurfaceExchangeConfig surfaceExchange; ///< surface heat-flux term
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
  /// Stage-3 T5: field capacity for the beta_lee_pielke ET limiter.
  /// Absent -> theta(h = -3.3 m) from the van Genuchten curve.
  bool hasThetaFc = false;  ///< theta_fc was configured
  real_t thetaFc = 0.0;     ///< configured field capacity [-]
  /// Stage-3 T7: surface albedo for the Penman-Monteith soil channel
  /// (the schema requires it on every type exactly when PM runs the
  /// soil channel).
  bool hasAlbedo = false;  ///< albedo was configured
  real_t albedo = 0.0;     ///< soil surface albedo [-]
};

/// vegetation.types[] (Stage-3 T6, plan §10 W6.1): one static canopy /
/// root description. Beer–Lambert transmission tau = exp(-k LAI) scales
/// the ground (open-water and soil) ET rates; the canopy transpires
/// T_pot = (1 - tau) E_x(r_s = r_c) with r_c = r_l / (0.5 LAI), taken
/// from the root zone as a Feddes-weighted volumetric sink (O7: it runs
/// regardless of inundation).
struct VegetationType {
  std::string name;            ///< unique type name
  real_t lai = 0.0;            ///< leaf area index [-]
  real_t extinction = 0.5;     ///< Beer–Lambert k [-]
  real_t leafResistance = 0.0; ///< r_l [s/m]
  real_t rootDepth = 0.0;      ///< L [m] below the column's top surface
  /// Root-density shape: uniform 1/L, or linear 2(L-d)/L^2 (triangular).
  enum class RootProfile { Uniform, Linear };
  RootProfile profile = RootProfile::Uniform;  ///< root-density shape
  /// Feddes stress heads [m], h1 >= h2 >= h3 >= h4. h1/h2 (the anaerobic
  /// cutoff) are optional: wetland species omit them (alpha = 1 for
  /// h >= h3 — O7's transpiration-while-inundated).
  bool hasH12 = false;         ///< h1/h2 configured
  real_t h1 = 0.0;             ///< alpha = 0 at/above h1 [m]
  real_t h2 = 0.0;             ///< alpha = 1 at h2 [m]
  real_t h3 = 0.0;             ///< alpha = 1 down to h3 [m]
  real_t h4 = 0.0;             ///< alpha = 0 at/below h4 [m]
  /// Canopy albedo (penman_monteith, T7; the schema requires it on every
  /// type exactly when PM runs with vegetation).
  bool hasAlbedo = false;      ///< albedo was configured
  real_t albedo = 0.0;         ///< canopy albedo [-]
};

/// vegetation: the static per-cell map (ids: 0 = bare, k >= 1 =
/// types[k-1]) with an optional per-cell LAI override raster.
struct VegetationConfig {
  bool enabled = false;               ///< block configured
  std::vector<VegetationType> types;  ///< type table
  bool mapFromFile = false;           ///< map is an id raster
  std::string mapFile;                ///< id raster (j*nx + i)
  std::string mapConstantName;        ///< uniform type name
  bool hasLaiFile = false;            ///< per-cell LAI override given
  std::string laiFile;                ///< LAI raster (j*nx + i)
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
  TransportInitialConfig temperature;    ///< temperature initial state (Q5)
};

/// boundary_conditions[].target: which sub-boundary the condition applies to.
enum class BcTarget { Surface, GroundwaterTop, GroundwaterBottom, GroundwaterSide };

/// boundary_conditions[].kind (plan §5.6): covers every legacy behavior.
/// Outflow is the free (transmissive) outflow decided by the P1 b4 gate
/// (plan §5.6 amendment): the prescribed-stage sink mapping kept the outlet
/// column dry, so SERGHEI's bctype-9 free outflow gets its own kind.
/// ScalarCauchy is the v2 Q4 zero-total-scalar-flux top condition (v2 plan
/// §3.2, Geng & Boufadel 2015 Eq. (7)): water crosses the subsurface top
/// face, scalar mass does not, and the top-cell limiter admits the exact
/// evaporative concentration/dilution the water loss implies. Homogeneous
/// (takes no value); target groundwater_top in groundwater-only runs.
enum class BcKind { Eta, Discharge, Velocity, Outflow, Head, Flux, ScalarValue, ScalarCauchy };

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

/// boundary_conditions[].scalar: which registered scalar a scalar_value
/// condition prescribes (v2 Q5, V2-A17). Salinity is the default; the
/// selector is only valid on kind scalar_value (schema cross-check).
enum class BcScalar { Salinity, Temperature };

/// boundary_conditions[].inflow.closure (Stage-3 T2, owner decision O1/O2):
/// the momentum state of the water BEYOND a prescribed-stage (kind eta)
/// boundary, read by the advection stencils at the reservoir faces.
/// ZeroGradient (default) is the v2.0 behavior, bitwise: the inflow
/// imports the interior momentum through the mass-consistent boundary
/// reconstruction — valid where the boundary cuts a uniform flow, a free
/// energy source where the outer water is still (the 2026-10 campaign
/// gyre, report_V1). External treats the stage cells as an external
/// reservoir moving at inflow.velocity (default still): self-advection
/// across a reservoir face takes the conservative d(u^2/2) form, so a
/// steady entrance costs exactly (u^2 - u_ext^2)/2g of head, and coupled
/// deposits on the stage cells leave through the reservoir (booked as
/// stage_deposit) instead of kicking the velocity update.
enum class BcInflowClosure { ZeroGradient, External };

/// boundary_conditions[]: one polygon-region boundary condition.
struct BoundaryConditionConfig {
  std::string name;   ///< unique condition name
  std::vector<std::array<real_t, 2>> polygon;  ///< region vertices (x, y) [m]
  BcTarget target = BcTarget::Surface;  ///< sub-boundary the condition acts on
  BcKind kind = BcKind::Eta;            ///< what quantity is prescribed
  BcScalar scalarField = BcScalar::Salinity;  ///< scalar_value selector (Q5)
  BcValueConfig value;                  ///< the prescribed value
  /// Stage-3 T2 (kind eta only): the inflow closure and, under External,
  /// the external reservoir's velocity components (constant or series;
  /// both default to a still reservoir).
  BcInflowClosure inflowClosure = BcInflowClosure::ZeroGradient;
  bool hasInflowBlock = false;          ///< the inflow key was given
  BcValueConfig inflowU;                ///< u_ext [m/s] along +x (External)
  BcValueConfig inflowV;                ///< v_ext [m/s] along +y (External)
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
  /// v2 Q4 (plan §3.2): keep the legacy hardcoded +0.01 evaporative-
  /// concentration limiter allowance on coupled dry evaporating columns
  /// instead of the exact in-step concentration factor. Off by default;
  /// exists so a b6-class golden can be pinned to the legacy behavior.
  bool legacyEvapAllowance = false;
  /// Stage-3 T5 (O4): what happens to a drying surface cell's remaining
  /// scalar mass. Delete (v2.0, golden-pinned) books it into surf_adjust
  /// and discards it; Store moves it into the per-cell salt store,
  /// released at the first wet step (the T8 DF evaluation may flip the
  /// default). Salinity only — temperature keeps the legacy rule.
  enum class SurfaceDryout { Delete, Store };
  SurfaceDryout surfaceDryout = SurfaceDryout::Store;  ///< selected mode (token: discard/store)
  /// transport.surface_update: consistent (v2.0 post-release opt-in): the
  /// surface concentration divides the transported mass by the volume of
  /// the *current* step's flow rates (the rates that moved the mass), and
  /// the local min/max limiter applies in every wet cell, still ones
  /// included. Off by default (`legacy`): the one-step-lagged flux volume
  /// and the moving-cell-only limiter are golden-pinned (docs/theory/
  /// transport.md quirk 2).
  bool consistentSurfaceUpdate = false;
  /// transport.subsurface_update: consistent (post-v2.0 opt-in): the
  /// subsurface scalar moves with every volume the groundwater module moves
  /// over the surface step. The face volumes are accumulated over all
  /// subcycled substeps (legacy: the last substep's fluxes times the
  /// surface dt), the post-allocation transfers are recorded on every face
  /// they cross (legacy: on the receiver's face only), the mass basis is the
  /// water content at the previous transport step, the surface-interface
  /// diffusion acts on every wet step (legacy: only on steps with a nonzero
  /// seepage rate) and the legacy one-sided interface dispersive gain is
  /// dropped. Off by default (`legacy`): the b6 goldens pin the legacy update.
  bool consistentSubsurfaceUpdate = false;
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
  std::vector<std::string> temperatureVariables;  ///< /temperature/* selection
  std::vector<std::string> atmosphereVariables;   ///< /atmosphere/* (T5 ET rates)
  std::vector<MonitorConfig> monitors;            ///< point probes
  /// Stage-3 T4 (plan §8 W4.3): monitor/audit ROW-WRITE decimation in
  /// whole seconds. 0 (the default) writes a row every step (the v2.0
  /// behavior, bitwise). Accumulation is always per step; rows are also
  /// forced at checkpoint instants and at t_end.
  real_t monitorInterval = 0.0;
  real_t checkpointInterval = 0.0;  ///< 0 = off (none written); > 0 also writes one at t_end
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
  AtmosphereConfig atmosphere;    ///< met forcing (v2 Q4)
  SurfaceWaterConfig surfaceWater;///< SWE parameters
  GroundwaterConfig groundwater;  ///< Richards parameters
  SoilConfig soil;                ///< soil table and map
  CouplingConfig coupling;        ///< coupling mode
  InitialConditionsConfig initialConditions;  ///< initial state
  std::vector<BoundaryConditionConfig> boundaryConditions;  ///< BC list
  EvapotranspirationConfig evapotranspiration;  ///< ET block (Stage-3 T5)
  VegetationConfig vegetation;    ///< vegetation block (Stage-3 T6)
  TransportConfig transport;      ///< scalar-transport parameters
  TemperatureConfig temperature;  ///< temperature-transport parameters (Q5)
  OutputConfig output;            ///< output selection
  RestartConfig restart;          ///< restart request
  SolverConfig solver;            ///< PETSc options hookup
  RuntimeConfig runtime;          ///< runtime toggles

  std::string configDir;  ///< directory of the YAML file (for relative paths)
  std::string rawText;    ///< full configuration text (embedded in HDF5 output)
  /// Non-fatal validation findings (Stage-3 W0.9 channel): keys that are
  /// legal but inert under the current mode, and similar advisories. Logged
  /// at run start and written to the run record's `warnings:` list; a
  /// configuration with warnings and no errors is VALID.
  std::vector<std::string> validationWarnings;

  /// Resolve a config-relative path against configDir.
  std::string resolvePath(const std::string& relative) const;
};

/// Outcome of configuration validation.
struct ValidationResult {
  std::vector<std::string> errors;  ///< empty when the configuration is valid
  /// Non-fatal advisories (Stage-3 W0.9): printed as WARNING by --validate
  /// and carried into FrehgConfig::validationWarnings; never invalidate.
  std::vector<std::string> warnings;
  /// \return true when no errors were recorded (warnings do not invalidate).
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
                  std::vector<std::string>& errors,
                  std::vector<std::string>& warnings);
}  // namespace detail

}  // namespace frehg

#endif  // FREHG_CORE_CONFIG_HPP
