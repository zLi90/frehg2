/// \file ConfigSchema.cpp
/// \brief Declarative YAML v2 schema and validation (plan §5.4, §6).
///
/// The schema is a value-semantic tree of Spec nodes built once per
/// validation. Validation walks the YAML document against the tree and
/// collects errors: type/range/allowed-value violations, missing required
/// keys, hard unknown-key rejection with a nearest-key suggestion, exactly-
/// one-of groups, referenced-file existence, and the cross-field rules listed
/// in plan §5.4.

#include "core/Config.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace frehg::detail {

namespace {

// ---------------------------------------------------------------------------
// Schema model
// ---------------------------------------------------------------------------

struct Spec;
struct KeyedSpec;
/// The named children of a map node. A named aggregate element rather than
/// std::pair<std::string, Spec>: instantiating pair's constructor traits
/// re-enters the enclosing, still-incomplete Spec through the children
/// member — ill-formed, and diagnosed as an incomplete-type error by
/// clang 18 against libstdc++ 14's completeness asserts (older libstdc++
/// and libc++ accepted the cycle by instantiation-order luck). std::vector
/// itself supports an incomplete element type, so vector<KeyedSpec> is
/// conformant with KeyedSpec defined once Spec is complete, below.
using KeyedSpecs = std::vector<KeyedSpec>;

/// One node of the declarative schema tree.
struct Spec {
  enum class Type { Map, Bool, Int, Real, String, Enum, IntOrAuto, Sequence, Point2 };

  Type type = Type::Map;
  bool required = false;

  // Scalar constraints.
  double minVal = -std::numeric_limits<double>::infinity();
  double maxVal = std::numeric_limits<double>::infinity();
  bool minExclusive = false;
  bool integerValued = false;           ///< Real that must land on whole numbers (§7)
  bool mustExist = false;               ///< String naming a file that must exist
  std::vector<std::string> allowed;     ///< Enum values

  // Map behavior.
  KeyedSpecs children;
  bool exactlyOne = false;              ///< exactly one child key must be present

  // Sequence behavior.
  std::vector<Spec> element;            ///< single-element vector = element spec
  std::size_t minItems = 0;
};

/// One map child: the YAML key and the spec of the node beneath it.
struct KeyedSpec {
  std::string key;   ///< the YAML map key
  Spec spec;         ///< the schema node the key maps to
};

Spec boolean(bool required = false) {
  Spec s;
  s.type = Spec::Type::Bool;
  s.required = required;
  return s;
}

Spec integer(bool required, double minVal) {
  Spec s;
  s.type = Spec::Type::Int;
  s.required = required;
  s.minVal = minVal;
  return s;
}

Spec real(bool required, double minVal = -std::numeric_limits<double>::infinity(),
          double maxVal = std::numeric_limits<double>::infinity(), bool minExclusive = false) {
  Spec s;
  s.type = Spec::Type::Real;
  s.required = required;
  s.minVal = minVal;
  s.maxVal = maxVal;
  s.minExclusive = minExclusive;
  return s;
}

Spec realPositive(bool required) { return real(required, 0.0, std::numeric_limits<double>::infinity(), true); }

Spec realNonNegative(bool required) { return real(required, 0.0); }

Spec integerSeconds(bool required, double minVal) {
  Spec s = real(required, minVal);
  s.integerValued = true;
  return s;
}

Spec text(bool required = false) {
  Spec s;
  s.type = Spec::Type::String;
  s.required = required;
  return s;
}

Spec path(bool required = false) {
  Spec s = text(required);
  s.mustExist = true;
  return s;
}

Spec enumeration(bool required, std::vector<std::string> allowed) {
  Spec s;
  s.type = Spec::Type::Enum;
  s.required = required;
  s.allowed = std::move(allowed);
  return s;
}

Spec intOrAuto() {
  Spec s;
  s.type = Spec::Type::IntOrAuto;
  s.minVal = 1.0;
  return s;
}

Spec map(KeyedSpecs children, bool required = false) {
  Spec s;
  s.type = Spec::Type::Map;
  s.required = required;
  s.children = std::move(children);
  return s;
}

Spec oneOf(KeyedSpecs children, bool required = false) {
  Spec s = map(std::move(children), required);
  s.exactlyOne = true;
  return s;
}

Spec sequence(Spec element, bool required = false, std::size_t minItems = 0) {
  Spec s;
  s.type = Spec::Type::Sequence;
  s.required = required;
  s.element.push_back(std::move(element));
  s.minItems = minItems;
  return s;
}

Spec point2() {
  Spec s;
  s.type = Spec::Type::Point2;
  return s;
}

/// {constant: <real>} | {file: <path>}
Spec fileOrConstant(bool required = false) {
  return oneOf({{"constant", real(false)}, {"file", path(false)}}, required);
}

/// {constant: <real>} | {series: {file: <path>}}
Spec seriesOrConstant(bool required = false) {
  return oneOf({{"constant", real(false)}, {"series", map({{"file", path(true)}})}}, required);
}

/// Allowed output variables per group (plan §7).
const std::vector<std::string>& surfaceVariables() {
  static const std::vector<std::string> v = {"eta",     "depth",       "uu", "vv",
                                             "seepage", "seepage_rate",
                                             // Stage-3 T9 (W9.1/W9.2).
                                             "film_removal_rate", "water_store"};
  return v;
}
const std::vector<std::string>& groundwaterVariables() {
  static const std::vector<std::string> v = {"hydraulic_head", "water_content",
                                             "qx", "qy", "qz",
                                             // Stage-3 T1 per-cell cumulative
                                             // post-allocation diagnostics
                                             // (plan W1.3) [m^3]
                                             "realloc_created", "realloc_dropped",
                                             "vloss",
                                             // Stage-3 T6: interval-mean
                                             // root-uptake sink [1/s].
                                             "root_uptake"};
  return v;
}
const std::vector<std::string>& transportVariables() {
  static const std::vector<std::string> v = {"concentration", "concentration_surface",
                                             // Stage-3 T5 (O4): the per-cell
                                             // surface salt store [psu m^3].
                                             "salt_store"};
  return v;
}
// Stage-3 T5/T6 (plan §13.C): the /atmosphere 2-D group — interval-mean
// ET rates (the T4 IntervalRates mechanism) and the cumulative clocks.
const std::vector<std::string>& atmosphereVariables() {
  static const std::vector<std::string> v = {
      "evaporation_rate", "soil_evaporation_rate",
      "soil_evaporation_unrealised_rate", "wet_time", "soil_evaporation_time",
      // Stage-3 T6: the transpiration channel.
      "transpiration_rate", "transpiration_unrealised_rate"};
  return v;
}
const std::vector<std::string>& temperatureVariables() {
  static const std::vector<std::string> v = {"temperature", "temperature_surface"};
  return v;
}

std::vector<std::string> allVariables() {
  std::vector<std::string> v = surfaceVariables();
  v.insert(v.end(), groundwaterVariables().begin(), groundwaterVariables().end());
  v.insert(v.end(), transportVariables().begin(), transportVariables().end());
  v.insert(v.end(), temperatureVariables().begin(), temperatureVariables().end());
  v.insert(v.end(), atmosphereVariables().begin(), atmosphereVariables().end());
  return v;
}

/// The complete v2 schema (plan §6).
Spec buildRootSchema() {
  Spec soilType = map({{"name", text(true)},
                       {"ksx", realNonNegative(true)},
                       {"ksy", realNonNegative(true)},
                       {"ksz", realNonNegative(true)},
                       {"theta_s", real(true, 0.0, 1.0, true)},
                       {"theta_r", real(true, 0.0, 1.0)},
                       {"vg_alpha", realPositive(true)},
                       {"vg_n", real(true, 1.0, std::numeric_limits<double>::infinity(), true)},
                       {"aev", real(false, -std::numeric_limits<double>::infinity(), 0.0)},
                       {"theta_fc", real(false, 0.0, 1.0)},
                       {"albedo", real(false, 0.0, 1.0)}});

  Spec bcEntry = map(
      {{"name", text(true)},
       {"region", map({{"polygon", sequence(point2(), true, 3)}}, true)},
       {"target", enumeration(true, {"surface", "groundwater_top", "groundwater_bottom",
                                     "groundwater_side"})},
       {"kind",
        enumeration(true, {"eta", "discharge", "velocity", "outflow", "head", "flux",
                           "scalar_value", "scalar_cauchy"})},
       // scalar_value selector (v2 Q5): which registered scalar the
       // condition prescribes; valid on kind scalar_value only
       // (cross-field check), default salinity.
       {"scalar", enumeration(false, {"salinity", "temperature"})},
       // value is required for every kind except outflow and scalar_cauchy
       // (cross-field check).
       {"value", oneOf({{"constant", real(false)},
                        {"series", map({{"file", path(true)}})},
                        {"gravity", boolean(false)},
                        {"hydrostatic", map({{"eta", real(true)}})}},
                       false)},
       // Stage-3 T2 (plan §13.A): the inflow closure of a prescribed-stage
       // condition — the momentum state of the water beyond the boundary.
       // kind eta only; velocity only with closure external (cross-field
       // checks).
       {"inflow",
        map({{"closure", enumeration(true, {"zero_gradient", "external"})},
             {"velocity", map({{"u", seriesOrConstant(false)},
                               {"v", seriesOrConstant(false)}})}})}});

  Spec monitorEntry = map({{"name", text(true)},
                           {"type", enumeration(false, {"point"})},
                           {"i", integer(true, 0.0)},
                           {"j", integer(true, 0.0)},
                           {"variables", sequence(enumeration(false, allVariables()), true, 1)}});

  return map({
      {"simulation", map({{"id", text(true)}, {"title", text(false)}}, true)},
      {"domain", map({{"nx", integer(true, 1.0)},
                      {"ny", integer(true, 1.0)},
                      {"nz", integer(true, 1.0)},
                      {"dx", realPositive(true)},
                      {"dy", realPositive(true)},
                      {"dz", realPositive(true)},
                      {"dz_stretch", realPositive(false)},
                      {"bottom_elevation", fileOrConstant(true)},
                      {"follow_terrain", boolean(false)},
                      {"terrain_layers", enumeration(false, {"scaled", "uniform"})},
                      // Stage-3 T3 (plan §13.B): periodic lateral topology
                      // per axis — a DOMAIN property (every module wraps);
                      // nx/ny >= 3 and the edge-BC rules are cross-field
                      // checks, the BC face-attachment rule is a setup
                      // fatal (BoundarySet).
                      {"periodic", map({{"x", boolean(false)}, {"y", boolean(false)}})},
                      {"decomposition", map({{"mpi_nx", intOrAuto()}, {"mpi_ny", intOrAuto()}})}},
                     true)},
      {"time", map({{"dt", realPositive(true)},
                    {"t_start", integerSeconds(false, 0.0)},
                    {"t_end", realPositive(true)},
                    {"output_interval", integerSeconds(true, 1.0)}},
                   true)},
      {"modules", map({{"surface_water", boolean(false)},
                       {"groundwater", boolean(false)},
                       {"transport", boolean(false)},
                       {"temperature", boolean(false)}},
                      true)},
      {"surface_water",
       map({{"gravity", realPositive(false)},
            {"friction", map({{"law", enumeration(false, {"manning", "chezy"})},
                              {"coefficient", fileOrConstant(true)},
                              {"thin_layer_depth", realPositive(false)}},
                             true)},
            {"viscosity", map({{"x", realNonNegative(false)}, {"y", realNonNegative(false)}})},
            {"min_depth", realPositive(true)},
            {"wetting_face_depth", realPositive(true)},
            // wind (v2 Q6 §5.2): selectable Cd(U10) law with cap, and
            // either the legacy compass speed/direction pair or the
            // grid-frame (u10, v10) component pair — exclusivity and the
            // cap-vs-constant rule are cross-field checks.
            {"wind", map({{"enabled", boolean(false)},
                          {"law", enumeration(false, {"constant", "garratt", "smith-banke",
                                                      "wu", "large-pond"})},
                          {"cd", realNonNegative(false)},
                          {"cap", realPositive(false)},
                          {"attenuation_depth", realPositive(false)},
                          {"north_angle", real(false)},
                          {"speed", seriesOrConstant(false)},
                          {"direction", seriesOrConstant(false)},
                          {"u10", seriesOrConstant(false)},
                          {"v10", seriesOrConstant(false)}})},
            // rainfall additionally accepts an exclusion region (no rain
            // inside; the legacy hardcoded last-row skip, made explicit);
            // the constant/series exclusivity is a cross-field check.
            {"rainfall", map({{"constant", real(false)},
                              {"series", map({{"file", path(true)}})},
                              {"exclude", map({{"polygon", sequence(point2(), true, 3)}})}})},
            // evaporation (v2 Q4 §3.2): prescribed rate (constant/series,
            // optional exclusion region — the rain symmetry) or the
            // bulk-aerodynamic mode; mode/value exclusivity and the
            // atmosphere-block dependency are cross-field checks.
            {"evaporation", map({{"mode", enumeration(false, {"bulk"})},
                                 {"constant", real(false)},
                                 {"series", map({{"file", path(true)}})},
                                 {"exclude", map({{"polygon", sequence(point2(), true, 3)}})}})}})},
      {"atmosphere",
       // Met forcing for the bulk-aerodynamic module (v2 Q4 §3.2); exactly
       // one humidity form (cross-field check). surface_temperature became
       // optional in v2 Q5 (required exactly when a Q4 evaporation consumer
       // is configured — cross-field check; the Q5 heat exchange evaluates
       // at the local water temperature instead, V2-A17). shortwave /
       // longwave_in / wind_speed_floor feed the Q5 bulk heat exchange.
       map({{"air_temperature", seriesOrConstant(true)},
            {"surface_temperature", seriesOrConstant(false)},
            {"pressure", seriesOrConstant(true)},
            {"specific_humidity", seriesOrConstant(false)},
            {"relative_humidity", seriesOrConstant(false)},
            {"wind_speed", seriesOrConstant(true)},
            {"shortwave", seriesOrConstant(false)},
            {"longwave_in", seriesOrConstant(false)},
            // Stage-3 T7 (plan §11 W7.1): INCIDENT R_s for the PM
            // potential (distinct from the Q5 absorbed `shortwave`).
            {"shortwave_incoming", seriesOrConstant(false)},
            {"wind_speed_floor", realNonNegative(false)},
            // Stage-3 T5 (plan §13.C): the ET consumers' r_a law.
            {"aerodynamic_resistance",
             map({{"law", enumeration(false, {"liu2006", "fao56", "log_profile"})},
                  {"reference_height", realPositive(false)},
                  {"z0m", realPositive(false)},
                  {"z0h", realPositive(false)},
                  {"displacement", realNonNegative(false)}})}})},
      // Stage-3 T5 (plan §9/§13.C): the ET block — replaces BOTH legacy
      // evaporation keys (cross-field rule). The soil sub-block's
      // presence enables the soil channel; its region defaults to all.
      // Stage-3 T7 (plan §11): penman_monteith validates (energy keys).
      {"evapotranspiration",
       map({{"potential", enumeration(true, {"bulk", "penman_monteith", "prescribed"})},
            {"prescribed", seriesOrConstant(false)},
            {"open_water", boolean(false)},
            {"soil",
             map({{"region", map({{"polygon", sequence(point2(), true, 3)}})},
                  {"limiter",
                   enumeration(false, {"alpha_barton", "beta_lee_pielke", "none"})}})},
            // T7 (plan §11 W7.1): the PM energy budget. Per-channel
            // albedo comes from here (water) / the soil type / the
            // vegetation type; G = fraction * R_n per channel.
            {"energy",
             map({{"albedo_water", real(false, 0.0, 1.0)},
                  {"ground_heat_fraction",
                   map({{"water", real(false, 0.0, 1.0)},
                        {"soil", real(false, 0.0, 1.0)},
                        {"canopy", real(false, 0.0, 1.0)}})}})}})},
      // Stage-3 T6 (plan §10/§13.C): static vegetation. Map ids: 0 = bare,
      // k >= 1 = types[k-1]; feddes h1/h2 optional (wetland species).
      {"vegetation",
       map({{"types",
             sequence(map({{"name", text(true)},
                           {"lai", realNonNegative(true)},
                           {"extinction", realPositive(true)},
                           {"leaf_resistance", realNonNegative(true)},
                           {"root_depth", realPositive(true)},
                           {"root_profile", enumeration(false, {"uniform", "linear"})},
                           {"feddes",
                            map({{"h1", real(false)},
                                 {"h2", real(false)},
                                 {"h3", real(true)},
                                 {"h4", real(true)}},
                                true)},
                           {"albedo", real(false, 0.0, 1.0)}}),
                      true, 1)},
            {"map", oneOf({{"constant", text(false)}, {"file", path(false)}}, true)},
            {"lai", map({{"file", text(true)}})}})},
      {"groundwater",
       map({{"scheme", enumeration(false, {"pca"})},
            {"use_full3d", boolean(false)},
            {"timestep", map({{"dt_init", realPositive(true)},
                              {"dt_min", realPositive(true)},
                              {"dt_max", realPositive(true)},
                              {"dq_grow", realPositive(false)},
                              {"dq_shrink", realPositive(false)},
                              {"courant_max", realPositive(false)}},
                             true)},
            {"specific_storage", realNonNegative(true)},
            // Stage-3 T1 (owner decision O9): run the post-allocation
            // step (legacy, the v2.0 default) or skip it entirely
            // (disabled: flux-conservative theta, retention-curve heads,
            // exact-zero realloc audit). The T8 default-flip evaluation may
            // change the default (plan §2.7). The token is "disabled", not
            // "off": YAML 1.1 parsers (PyYAML — the gate harness, the
            // run-record checker, user postprocessing) read off as the
            // boolean false, which corrupts every round trip (S3-A3).
            {"reallocation", enumeration(false, {"legacy", "disabled"})},
            {"predictor_corrector", enumeration(false, {"legacy", "consistent"})},
            {"face_conductivity", enumeration(false, {"arithmetic", "harmonic"})},
            {"reallocation_surplus", enumeration(false, {"drop", "redistribute"})},
            {"density_coupling",
             // Coefficients configurable since v2 Q5 (defaults = the
             // legacy compile-time constants; thermal term default 0).
             map({{"enabled", boolean(false)},
                  {"beta_saline", realNonNegative(false)},
                  {"beta_saline_viscosity", realNonNegative(false)},
                  {"thermal_expansion", realNonNegative(false)},
                  {"reference_temperature", real(false)},
                  {"reallocation_gradient", enumeration(false, {"freshwater", "density"})}})},
            // Bulk-aerodynamic soil evaporation over a region (v2 Q4
            // §3.2); dependencies (atmosphere block, uncoupled runs) are
            // cross-field checks.
            {"evaporation",
             map({{"mode", enumeration(true, {"bulk"})},
                  {"region", map({{"polygon", sequence(point2(), true, 3)}}, true)}})}})},
      {"soil", map({{"types", sequence(soilType, true, 1)},
                    {"map", oneOf({{"constant", text(false)}, {"file", path(false)}}, true)}})},
      {"coupling", map({{"mode", enumeration(false, {"sync", "subcycled"})}})},
      {"initial_conditions",
       map({{"surface", map({{"eta", fileOrConstant(true)},
                             {"uu", fileOrConstant(false)},
                             {"vv", fileOrConstant(false)}})},
            {"groundwater", oneOf({{"water_table", fileOrConstant(false)},
                                   {"head", fileOrConstant(false)},
                                   {"moisture", fileOrConstant(false)}})},
            {"transport", map({{"surface", fileOrConstant(false)},
                               {"groundwater", fileOrConstant(false)}})},
            {"temperature", map({{"surface", fileOrConstant(false)},
                                 {"groundwater", fileOrConstant(false)}})}})},
      {"boundary_conditions", sequence(bcEntry)},
      {"transport",
       map({{"scheme", map({{"advection", enumeration(false, {"upwind", "superbee"})}})},
            {"surface_diffusivity",
             map({{"x", realNonNegative(false)}, {"y", realNonNegative(false)}})},
            {"dispersion", map({{"longitudinal", realNonNegative(false)},
                                {"transverse", realNonNegative(false)},
                                {"molecular", realNonNegative(false)}})},
            {"bounds", map({{"min", real(false)}, {"max", real(false)}})},
            {"legacy_evap_allowance", boolean(false)},
            // Stage-3 T10 (§12B): conservative = the solve-consistent
            // surface update (salt moves with the solve's face volumes).
            {"surface_update", enumeration(false, {"legacy", "consistent", "conservative"})},
            {"subsurface_update", enumeration(false, {"legacy", "consistent"})},
            // Stage-3 T5 (O4): drying-cell scalar fate; the T8 DF constant.
            {"surface_dryout", enumeration(false, {"discard", "store"})}})},
      {"temperature",
       // The second registered scalar (v2 Q5, plan §4.1/V2-A17). The
       // thermal parameters' module dependencies are cross-field checks.
       map({{"scheme", map({{"advection", enumeration(false, {"upwind", "superbee"})}})},
            {"surface_diffusivity",
             map({{"x", realNonNegative(false)}, {"y", realNonNegative(false)}})},
            {"thermal_conductivity", realPositive(false)},
            {"heat_capacity_water", realPositive(false)},
            {"heat_capacity_solid", realPositive(false)},
            {"dispersivity", map({{"longitudinal", realNonNegative(false)},
                                  {"transverse", realNonNegative(false)}})},
            {"bounds", map({{"min", real(false)}, {"max", real(false)}})},
            {"surface_exchange",
             map({{"mode", enumeration(true, {"equilibrium", "bulk"})},
                  {"equilibrium", map({{"temperature", real(false)},
                                       {"coefficient", realPositive(false)}})}})}})},
      {"output",
       map({{"filename", text(true)},
            {"variables",
             map({{"surface", sequence(enumeration(false, surfaceVariables()), false, 1)},
                  {"groundwater",
                   sequence(enumeration(false, groundwaterVariables()), false, 1)},
                  {"transport", sequence(enumeration(false, transportVariables()), false, 1)},
                  {"temperature",
                   sequence(enumeration(false, temperatureVariables()), false, 1)},
                  {"atmosphere",
                   sequence(enumeration(false, atmosphereVariables()), false, 1)}})},
            {"monitors", sequence(monitorEntry)},
            {"monitor_interval", integerSeconds(false, 0.0)},
            {"checkpoint", map({{"interval", integerSeconds(false, 0.0)}})}},
           true)},
      {"restart", map({{"enabled", boolean(false)}, {"file", text(false)}, {"time", real(false)}})},
      {"solver", map({{"surface", map({{"preconditioner",
                                        enumeration(false, {"bjacobi-icc", "amg", "gamg"})},
                                       {"mat_type", enumeration(false, {"aij", "aijkokkos"})},
                                       {"rtol", realPositive(false)},
                                       {"atol", realNonNegative(false)},
                                       {"max_iterations", integer(false, 1.0)},
                                       {"reuse_max_solves", integer(false, 0.0)},
                                       {"reuse_iteration_factor", real(false, 1.0)}})},
                      {"groundwater", map({{"preconditioner",
                                            enumeration(false, {"bjacobi-icc", "amg", "gamg"})},
                                           {"mat_type", enumeration(false, {"aij", "aijkokkos"})},
                                           {"rtol", realPositive(false)},
                                           {"atol", realNonNegative(false)},
                                           {"max_iterations", integer(false, 1.0)},
                                           {"reuse_max_solves", integer(false, 0.0)},
                                           {"reuse_iteration_factor", real(false, 1.0)}})},
                      {"petsc_options_file", text(false)}})},
      {"runtime", map({{"gpu_aware_mpi", enumeration(false, {"auto", "on", "off"})}})},
  });
}

// ---------------------------------------------------------------------------
// Validation walker
// ---------------------------------------------------------------------------

struct Context {
  std::string configDir;
  std::vector<std::string>* errors = nullptr;
  std::vector<std::string>* warnings = nullptr;

  void addError(const std::string& where, const std::string& what) {
    errors->push_back((where.empty() ? "(document)" : where) + ": " + what);
  }

  /// Non-fatal advisory (Stage-3 W0.9): a legal-but-inert key or similar.
  /// Printed by --validate as WARNING, logged at run start, and written to
  /// the run record's `warnings:` list. Never invalidates.
  void addWarning(const std::string& where, const std::string& what) {
    warnings->push_back((where.empty() ? "(document)" : where) + ": " + what);
  }
};

/// Join a parent path and a key ("" + "domain" -> "domain").
std::string joinPath(const std::string& where, const std::string& key) {
  return where.empty() ? key : where + "." + key;
}

/// Levenshtein distance for the unknown-key suggestion.
std::size_t editDistance(const std::string& a, const std::string& b) {
  const std::size_t na = a.size();
  const std::size_t nb = b.size();
  std::vector<std::size_t> prev(nb + 1);
  std::vector<std::size_t> curr(nb + 1);
  for (std::size_t jb = 0; jb <= nb; ++jb) {
    prev[jb] = jb;
  }
  for (std::size_t ia = 1; ia <= na; ++ia) {
    curr[0] = ia;
    for (std::size_t jb = 1; jb <= nb; ++jb) {
      const std::size_t sub = prev[jb - 1] + (a[ia - 1] == b[jb - 1] ? 0 : 1);
      curr[jb] = std::min({prev[jb] + 1, curr[jb - 1] + 1, sub});
    }
    std::swap(prev, curr);
  }
  return prev[nb];
}

std::string nearestKey(const std::string& key, const KeyedSpecs& children) {
  std::string best;
  std::size_t bestDist = std::numeric_limits<std::size_t>::max();
  for (const auto& [name, spec] : children) {
    const std::size_t d = editDistance(key, name);
    if (d < bestDist) {
      bestDist = d;
      best = name;
    }
  }
  if (!best.empty() && bestDist <= std::max<std::size_t>(2, key.size() / 3)) {
    return best;
  }
  return std::string();
}

std::string joinAllowed(const std::vector<std::string>& allowed) {
  std::ostringstream out;
  for (std::size_t n = 0; n < allowed.size(); ++n) {
    out << (n ? ", " : "") << allowed[n];
  }
  return out.str();
}

bool nodeAsDouble(const YAML::Node& node, double& out) {
  try {
    out = node.as<double>();
    return true;
  } catch (const YAML::Exception&) {
    return false;
  }
}

bool nodeAsBool(const YAML::Node& node, bool& out) {
  try {
    out = node.as<bool>();
    return true;
  } catch (const YAML::Exception&) {
    return false;
  }
}

void checkNode(const YAML::Node& node, const Spec& spec, const std::string& where, Context& ctx);

void checkMap(const YAML::Node& node, const Spec& spec, const std::string& where, Context& ctx) {
  if (!node.IsMap()) {
    ctx.addError(where, "expected a mapping");
    return;
  }
  // Unknown-key rejection with suggestion.
  for (const auto& item : node) {
    const std::string key = item.first.as<std::string>();
    const bool known = std::any_of(spec.children.begin(), spec.children.end(),
                                   [&key](const auto& kv) { return kv.key == key; });
    if (!known) {
      const std::string suggestion = nearestKey(key, spec.children);
      std::string message = "unknown key '" + key + "'";
      if (!suggestion.empty()) {
        message += " (did you mean '" + suggestion + "'?)";
      }
      ctx.addError(where, message);
    }
  }
  // Required keys, exactly-one groups, recursion.
  std::size_t presentCount = 0;
  for (const auto& [key, child] : spec.children) {
    const YAML::Node childNode = node[key];
    if (childNode.IsDefined()) {
      ++presentCount;
      checkNode(childNode, child, joinPath(where, key), ctx);
    } else if (child.required && !spec.exactlyOne) {
      ctx.addError(where, "missing required key '" + key + "'");
    }
  }
  if (spec.exactlyOne && presentCount != 1) {
    std::vector<std::string> names;
    names.reserve(spec.children.size());
    for (const auto& [key, child] : spec.children) {
      names.push_back(key);
    }
    ctx.addError(where, "exactly one of {" + joinAllowed(names) + "} must be given (found " +
                            std::to_string(presentCount) + ")");
  }
}

void checkNode(const YAML::Node& node, const Spec& spec, const std::string& where, Context& ctx) {
  switch (spec.type) {
    case Spec::Type::Map:
      checkMap(node, spec, where, ctx);
      return;

    case Spec::Type::Bool: {
      bool value = false;
      if (!node.IsScalar() || !nodeAsBool(node, value)) {
        ctx.addError(where, "expected a boolean (true/false)");
      }
      return;
    }

    case Spec::Type::Int: {
      double value = 0;
      if (!node.IsScalar() || !nodeAsDouble(node, value) ||
          value != std::floor(value)) {
        ctx.addError(where, "expected an integer");
        return;
      }
      if (value < spec.minVal) {
        ctx.addError(where, "must be >= " + std::to_string(static_cast<long>(spec.minVal)) +
                                " (got " + node.as<std::string>() + ")");
      }
      return;
    }

    case Spec::Type::IntOrAuto: {
      if (node.IsScalar() && node.as<std::string>() == "auto") {
        return;
      }
      double value = 0;
      if (!node.IsScalar() || !nodeAsDouble(node, value) || value != std::floor(value) ||
          value < spec.minVal) {
        ctx.addError(where, "expected 'auto' or an integer >= 1");
      }
      return;
    }

    case Spec::Type::Real: {
      double value = 0;
      if (!node.IsScalar() || !nodeAsDouble(node, value)) {
        ctx.addError(where, "expected a number");
        return;
      }
      if (spec.minExclusive ? !(value > spec.minVal) : !(value >= spec.minVal)) {
        std::ostringstream out;
        out << "must be " << (spec.minExclusive ? "> " : ">= ") << spec.minVal << " (got " << value
            << ")";
        ctx.addError(where, out.str());
      }
      if (!(value <= spec.maxVal)) {
        std::ostringstream out;
        out << "must be <= " << spec.maxVal << " (got " << value << ")";
        ctx.addError(where, out.str());
      }
      if (spec.integerValued && std::abs(value - std::round(value)) > 1.0e-9) {
        std::ostringstream out;
        out << "must land on whole seconds (got " << value
            << "); the HDF5 layout keys times by integer seconds";
        ctx.addError(where, out.str());
      }
      return;
    }

    case Spec::Type::String: {
      if (!node.IsScalar()) {
        ctx.addError(where, "expected a string");
        return;
      }
      if (spec.mustExist) {
        const std::filesystem::path p =
            std::filesystem::path(ctx.configDir) / node.as<std::string>();
        if (!std::filesystem::exists(p)) {
          ctx.addError(where, "referenced file '" + p.string() + "' does not exist");
        }
      }
      return;
    }

    case Spec::Type::Enum: {
      if (!node.IsScalar()) {
        ctx.addError(where, "expected one of {" + joinAllowed(spec.allowed) + "}");
        return;
      }
      const std::string value = node.as<std::string>();
      if (std::find(spec.allowed.begin(), spec.allowed.end(), value) == spec.allowed.end()) {
        ctx.addError(where,
                     "invalid value '" + value + "' (allowed: " + joinAllowed(spec.allowed) + ")");
      }
      return;
    }

    case Spec::Type::Sequence: {
      if (!node.IsSequence()) {
        ctx.addError(where, "expected a list");
        return;
      }
      if (node.size() < spec.minItems) {
        ctx.addError(where, "needs at least " + std::to_string(spec.minItems) + " entr" +
                                (spec.minItems == 1 ? "y" : "ies"));
      }
      for (std::size_t n = 0; n < node.size(); ++n) {
        checkNode(node[n], spec.element.front(), where + "[" + std::to_string(n) + "]", ctx);
      }
      return;
    }

    case Spec::Type::Point2: {
      double value = 0;
      if (!node.IsSequence() || node.size() != 2 || !nodeAsDouble(node[0], value) ||
          !nodeAsDouble(node[1], value)) {
        ctx.addError(where, "expected an [x, y] coordinate pair");
      }
      return;
    }
  }
}

// ---------------------------------------------------------------------------
// Cross-field checks (plan §5.4)
// ---------------------------------------------------------------------------

bool moduleOn(const YAML::Node& root, const char* name) {
  const YAML::Node modules = root["modules"];
  if (!modules.IsDefined() || !modules.IsMap()) {
    return false;
  }
  const YAML::Node flag = modules[name];
  bool value = false;
  return flag.IsDefined() && flag.IsScalar() && nodeAsBool(flag, value) && value;
}

/// Safe nested lookup: returns an undefined node instead of throwing when
/// any parent along the way is absent or not a mapping (yaml-cpp rejects
/// operator[] on undefined nodes).
YAML::Node sub(const YAML::Node& node, const std::string& key) {
  if (!node.IsDefined() || !node.IsMap()) {
    return YAML::Node(YAML::NodeType::Undefined);
  }
  return node[key];
}

/// \copydoc sub
YAML::Node sub(const YAML::Node& node, const std::string& key1, const std::string& key2) {
  return sub(sub(node, key1), key2);
}

double realOr(const YAML::Node& node, double fallback) {
  double value = fallback;
  if (node.IsDefined() && node.IsScalar() && nodeAsDouble(node, value)) {
    return value;
  }
  return fallback;
}

void crossChecks(const YAML::Node& root, Context& ctx) {
  const bool sw = moduleOn(root, "surface_water");
  const bool gw = moduleOn(root, "groundwater");
  const bool tr = moduleOn(root, "transport");
  const bool tp = moduleOn(root, "temperature");

  if (!sw && !gw) {
    ctx.addError("modules", "at least one of surface_water/groundwater must be enabled");
  }
  if (tr && !sw && !gw) {
    ctx.addError("modules", "transport requires a flow module (surface_water or groundwater)");
  }
  if (tp && !sw && !gw) {
    ctx.addError("modules",
                 "temperature requires a flow module (surface_water or groundwater)");
  }
  if (sw && !root["surface_water"].IsDefined()) {
    ctx.addError("surface_water", "section is required when modules.surface_water is true");
  }
  if (gw && !root["groundwater"].IsDefined()) {
    ctx.addError("groundwater", "section is required when modules.groundwater is true");
  }
  if (gw && !root["soil"].IsDefined()) {
    ctx.addError("soil", "section is required when modules.groundwater is true");
  }
  if (tr && !root["transport"].IsDefined()) {
    ctx.addError("transport", "section is required when modules.transport is true");
  }
  if (tp && !root["temperature"].IsDefined()) {
    ctx.addError("temperature", "section is required when modules.temperature is true");
  }

  // v2 Q5 temperature cross-field rules (plan §4.1, V2-A17).
  {
    const YAML::Node temp = root["temperature"];
    if (temp.IsDefined() && temp.IsMap()) {
      if (gw && tp) {
        // The subsurface heat physics has no defaults: conduction and the
        // solid heat capacity are the case's thermal identity.
        if (!temp["thermal_conductivity"].IsDefined()) {
          ctx.addError("temperature.thermal_conductivity",
                       "required when modules.groundwater is true");
        }
        if (!temp["heat_capacity_solid"].IsDefined()) {
          ctx.addError("temperature.heat_capacity_solid",
                       "required when modules.groundwater is true");
        }
      }
      if (!gw) {
        for (const char* key : {"thermal_conductivity", "heat_capacity_solid",
                                "dispersivity"}) {
          if (temp[key].IsDefined()) {
            ctx.addError(std::string("temperature.") + key,
                         "subsurface thermal parameter requires modules.groundwater: true");
          }
        }
      }
      const YAML::Node exchange = temp["surface_exchange"];
      if (exchange.IsDefined() && exchange.IsMap()) {
        if (!sw) {
          ctx.addError("temperature.surface_exchange",
                       "requires modules.surface_water: true");
        }
        const std::string mode =
            (exchange["mode"].IsDefined() && exchange["mode"].IsScalar())
                ? exchange["mode"].as<std::string>()
                : std::string();
        const YAML::Node eq = exchange["equilibrium"];
        if (mode == "equilibrium") {
          if (!eq.IsDefined() || !sub(eq, "temperature").IsDefined() ||
              !sub(eq, "coefficient").IsDefined()) {
            ctx.addError("temperature.surface_exchange.equilibrium",
                         "mode equilibrium requires equilibrium.temperature and "
                         "equilibrium.coefficient");
          }
        }
        if (mode == "bulk") {
          if (eq.IsDefined()) {
            ctx.addError("temperature.surface_exchange.equilibrium",
                         "mode bulk takes no equilibrium block (the flux comes from "
                         "the atmosphere block)");
          }
          if (!root["atmosphere"].IsDefined()) {
            ctx.addError("temperature.surface_exchange",
                         "mode bulk requires the atmosphere block");
          }
        }
      }
      if (temp["surface_diffusivity"].IsDefined() && !sw) {
        ctx.addError("temperature.surface_diffusivity",
                     "requires modules.surface_water: true");
      }
    }
  }

  // time ordering.
  {
    const YAML::Node time = root["time"];
    if (time.IsDefined() && time.IsMap()) {
      const double tStart = realOr(time["t_start"], 0.0);
      const double tEnd = realOr(time["t_end"], 0.0);
      if (time["t_end"].IsDefined() && !(tEnd > tStart)) {
        ctx.addError("time.t_end", "must be greater than time.t_start");
      }
    }
  }

  // rainfall: exactly one of constant/series (the map also carries the
  // optional exclude region, so the generic exactly-one group cannot apply).
  {
    const YAML::Node rainfall = sub(root, "surface_water", "rainfall");
    if (rainfall.IsDefined() && rainfall.IsMap()) {
      const int present = (rainfall["constant"].IsDefined() ? 1 : 0) +
                          (rainfall["series"].IsDefined() ? 1 : 0);
      if (present != 1) {
        ctx.addError("surface_water.rainfall",
                     "exactly one of {constant, series} must be given (found " +
                         std::to_string(present) + ")");
      }
    }
  }

  // v2 Q6 wind cross-field rules (plan §5.2).
  {
    const YAML::Node wind = sub(root, "surface_water", "wind");
    if (wind.IsDefined() && wind.IsMap()) {
      const bool components = wind["u10"].IsDefined() || wind["v10"].IsDefined();
      const bool compass = wind["speed"].IsDefined() || wind["direction"].IsDefined();
      if (components && compass) {
        ctx.addError("surface_water.wind",
                     "give either the (u10, v10) component pair or the "
                     "speed/direction pair, not both");
      }
      if (components && (!wind["u10"].IsDefined() || !wind["v10"].IsDefined())) {
        ctx.addError("surface_water.wind",
                     "the component form requires both u10 and v10");
      }
      const std::string law =
          (wind["law"].IsDefined() && wind["law"].IsScalar())
              ? wind["law"].as<std::string>()
              : std::string("constant");
      if (law == "constant" && wind["cap"].IsDefined()) {
        ctx.addError("surface_water.wind.cap",
                     "the cap applies to the U10 laws only (the constant law is "
                     "the uncapped legacy Cw)");
      }
      if (law != "constant" && wind["cd"].IsDefined()) {
        ctx.addError("surface_water.wind.cd",
                     "cd is the constant law's coefficient; the U10 laws compute "
                     "Cd from the wind speed");
      }
    }
  }

  // v2 Q4 evaporation / atmosphere cross-field rules (plan §3.2).
  {
    const YAML::Node atm = root["atmosphere"];
    const bool atmPresent = atm.IsDefined();
    if (atmPresent && atm.IsMap()) {
      const int humid = (atm["specific_humidity"].IsDefined() ? 1 : 0) +
                        (atm["relative_humidity"].IsDefined() ? 1 : 0);
      if (humid != 1) {
        ctx.addError("atmosphere",
                     "exactly one of {specific_humidity, relative_humidity} must be "
                     "given (found " + std::to_string(humid) + ")");
      }
    }
    const YAML::Node evap = sub(root, "surface_water", "evaporation");
    const bool surfBulkEvap = evap.IsDefined() && evap.IsMap() && evap["mode"].IsDefined();
    const bool gwBulkEvap = sub(root, "groundwater", "evaporation").IsDefined();
    if (atmPresent && atm.IsMap() && (surfBulkEvap || gwBulkEvap) &&
        !atm["surface_temperature"].IsDefined()) {
      ctx.addError("atmosphere.surface_temperature",
                   "required when a bulk evaporation consumer is configured "
                   "(surface_water.evaporation.mode bulk or "
                   "groundwater.evaporation)");
    }
    if (evap.IsDefined() && evap.IsMap()) {
      const bool bulk = evap["mode"].IsDefined();
      const int present = (evap["constant"].IsDefined() ? 1 : 0) +
                          (evap["series"].IsDefined() ? 1 : 0);
      if (bulk && (present > 0 || evap["exclude"].IsDefined())) {
        ctx.addError("surface_water.evaporation",
                     "mode: bulk takes no constant/series/exclude (the rate comes "
                     "from the atmosphere block, wet cells only)");
      }
      if (!bulk && present != 1) {
        ctx.addError("surface_water.evaporation",
                     "exactly one of {constant, series} must be given (found " +
                         std::to_string(present) + ")");
      }
      if (bulk && !atmPresent) {
        ctx.addError("surface_water.evaporation",
                     "mode: bulk requires the atmosphere block");
      }
    }
    const YAML::Node gwEvap = sub(root, "groundwater", "evaporation");
    if (gwEvap.IsDefined()) {
      if (!atmPresent) {
        ctx.addError("groundwater.evaporation",
                     "bulk soil evaporation requires the atmosphere block");
      }
      if (sw && gw) {
        // Same restriction (and reason) as scalar_cauchy: the coupled top
        // exchange is owned by the coupler, and no gate exercises the
        // coupled combination — §8.3 forbids shipping it unexercised.
        ctx.addError("groundwater.evaporation",
                     "bulk soil evaporation applies to uncoupled groundwater runs "
                     "only (the coupler owns the top exchange) — coupled soil "
                     "evaporation is configured through the evapotranspiration "
                     "block");
      }
    }
  }

  // Stage-3 T5 (plan §13.C): the evapotranspiration block's rules.
  {
    const YAML::Node et = root["evapotranspiration"];
    if (et.IsDefined() && et.IsMap()) {
      const std::string mode =
          et["potential"].IsDefined() ? et["potential"].as<std::string>("") : "bulk";
      const bool atmPresent = root["atmosphere"].IsDefined();
      const YAML::Node atm = root["atmosphere"];
      if (sub(root, "surface_water", "evaporation").IsDefined() ||
          sub(root, "groundwater", "evaporation").IsDefined()) {
        ctx.addError("evapotranspiration",
                     "replaces surface_water.evaporation and "
                     "groundwater.evaporation; give one form");
      }
      if ((mode == "bulk" || mode == "penman_monteith") && !atmPresent) {
        ctx.addError("evapotranspiration.potential",
                     mode + ": requires the atmosphere block");
      }
      // Stage-3 T7 (plan §13.C): PM needs the incident radiation pair;
      // the albedo of each configured channel must be declared (water in
      // the energy block's default; soil / canopy on their types).
      if (mode == "penman_monteith" && atmPresent && atm.IsMap()) {
        for (const char* key : {"shortwave_incoming", "longwave_in"}) {
          if (!atm[key].IsDefined()) {
            ctx.addError("evapotranspiration.potential",
                         std::string("penman_monteith: requires atmosphere.") + key);
          }
        }
      }
      if (mode == "penman_monteith") {
        if (et["soil"].IsDefined()) {
          const YAML::Node soilTypes = sub(root, "soil", "types");
          if (soilTypes.IsDefined() && soilTypes.IsSequence()) {
            for (std::size_t i = 0; i < soilTypes.size(); ++i) {
              if (!soilTypes[i]["albedo"].IsDefined()) {
                ctx.addError("soil.types[" + std::to_string(i) + "]",
                             "penman_monteith with the soil channel needs every "
                             "soil type's albedo");
              }
            }
          }
        }
        const YAML::Node vegTypes = sub(root, "vegetation", "types");
        if (vegTypes.IsDefined() && vegTypes.IsSequence()) {
          for (std::size_t i = 0; i < vegTypes.size(); ++i) {
            if (!vegTypes[i]["albedo"].IsDefined()) {
              ctx.addError("vegetation.types[" + std::to_string(i) + "]",
                           "penman_monteith with vegetation needs every "
                           "vegetation type's albedo");
            }
          }
        }
      }
      if (mode == "bulk" && atmPresent && atm.IsMap() &&
          !atm["surface_temperature"].IsDefined()) {
        ctx.addError("evapotranspiration.potential",
                     "bulk: requires atmosphere.surface_temperature");
      }
      if (mode == "prescribed" && !et["prescribed"].IsDefined()) {
        ctx.addError("evapotranspiration.potential",
                     "prescribed: requires evapotranspiration.prescribed");
      }
      const YAML::Node soilBlock = et["soil"];
      if (soilBlock.IsDefined()) {
        if (!gw) {
          ctx.addError("evapotranspiration.soil", "needs modules.groundwater");
        }
        const std::string limiter = soilBlock.IsMap() && soilBlock["limiter"].IsDefined()
                                        ? soilBlock["limiter"].as<std::string>("")
                                        : "alpha_barton";
        if (limiter == "alpha_barton" && mode != "bulk") {
          ctx.addError("evapotranspiration.soil.limiter",
                       "alpha_barton acts inside the bulk humidity difference; "
                       "use beta_lee_pielke or none with potential " + mode);
        }
      }
      const bool openWater = !et["open_water"].IsDefined() ||
                             et["open_water"].as<bool>(true);
      if (openWater && !sw) {
        ctx.addError("evapotranspiration.open_water",
                     "needs modules.surface_water");
      }
      if (mode == "penman_monteith" && root["modules"].IsDefined() &&
          root["modules"]["temperature"].IsDefined() &&
          root["modules"]["temperature"].as<bool>(false)) {
        ctx.addError("evapotranspiration",
                     "penman_monteith with modules.temperature is not available "
                     "(the heat module's latent flux uses the bulk rate at the "
                     "simulated temperature; the combination is ungated)");
      }
      const YAML::Node aero = sub(root, "atmosphere", "aerodynamic_resistance");
      const std::string raLaw = aero.IsDefined() && aero.IsMap() && aero["law"].IsDefined()
                                    ? aero["law"].as<std::string>("")
                                    : "liu2006";
      if (mode == "bulk" && raLaw != "liu2006" && root["modules"].IsDefined() &&
          root["modules"]["temperature"].IsDefined() &&
          root["modules"]["temperature"].as<bool>(false)) {
        ctx.addError("atmosphere.aerodynamic_resistance",
                     "with modules.temperature the evaporated volume and the heat "
                     "module's latent flux must share the liu2006 resistance");
      }
      // ET + salinity: the exact salt bookkeeping needs the consistent
      // updates of the modules present, and (with a surface) the store.
      // Since the S3-A18 default flip an ABSENT surface_update resolves to
      // conservative (and an absent subsurface_update beside it to
      // consistent), so only EXPLICIT legacy values violate the rule.
      if (tr) {
        const YAML::Node trans = root["transport"];
        const auto explicitLegacy = [&trans](const char* key) {
          return trans.IsDefined() && trans.IsMap() && trans[key].IsDefined() &&
                 trans[key].as<std::string>("") == "legacy";
        };
        const auto surfaceExplicit = [&trans]() -> std::string {
          if (!(trans.IsDefined() && trans.IsMap() &&
                trans["surface_update"].IsDefined())) {
            return "";  // absent: resolves to conservative
          }
          return trans["surface_update"].as<std::string>("");
        };
        const bool surfaceOk = sw ? surfaceExplicit() != "legacy" : true;
        // The subsurface default is consistent only under the resolved
        // conservative surface scheme; an explicit consistent/absent-under-
        // conservative passes, an explicit legacy (or an absent key under an
        // explicit consistent surface) fails.
        const bool subAbsent = !(trans.IsDefined() && trans.IsMap() &&
                                 trans["subsurface_update"].IsDefined());
        const bool subOk =
            !gw || (subAbsent
                        ? (surfaceExplicit().empty() ||
                           surfaceExplicit() == "conservative")
                        : !explicitLegacy("subsurface_update"));
        if (!surfaceOk || !subOk) {
          ctx.addError("evapotranspiration",
                       "with salinity transport requires transport.surface_update: "
                       "consistent (with modules.surface_water) and "
                       "transport.subsurface_update: consistent (with "
                       "modules.groundwater)");
        }
        // Stage-3 T8 (W8.2/W8.3): `store` is the shipped default now, so
        // only an EXPLICIT `discard` violates the rule (an absent key
        // resolves to store).
        if (sw && trans.IsDefined() && trans.IsMap() &&
            trans["surface_dryout"].IsDefined() &&
            trans["surface_dryout"].as<std::string>("") == "discard") {
          ctx.addError("evapotranspiration",
                       "with salinity transport requires transport.surface_dryout: "
                       "store (salt left by evaporating water must be kept)");
        }
      }
    }
  }

  // Stage-3 T3 (plan §13.B): periodic-axis cross-field rules. A periodic
  // axis needs >= 3 cells (the 5/7-point stencils would otherwise read the
  // same cell on both sides), and an open surface condition (eta /
  // discharge / outflow) listed on a periodic axis's edge has no boundary
  // to act on — the polygon test here is conservative (bounding-box
  // against the edge strips); the exact member-face rule is the BoundarySet
  // setup fatal.
  {
    const YAML::Node per = sub(root, "domain", "periodic");
    if (per.IsDefined() && per.IsMap()) {
      const bool perX = per["x"].IsDefined() && per["x"].as<bool>(false);
      const bool perY = per["y"].IsDefined() && per["y"].as<bool>(false);
      const YAML::Node dom = root["domain"];
      if (perX && dom["nx"].IsDefined() && dom["nx"].as<int>(0) < 3) {
        ctx.addError("domain.periodic.x",
                     "needs nx >= 3 (got " + std::to_string(dom["nx"].as<int>(0)) + ")");
      }
      if (perY && dom["ny"].IsDefined() && dom["ny"].as<int>(0) < 3) {
        ctx.addError("domain.periodic.y",
                     "needs ny >= 3 (got " + std::to_string(dom["ny"].as<int>(0)) + ")");
      }
    }
  }

  // Stage-3 T1 (W1.1): with reallocation: disabled the surplus and gradient
  // knobs steer a step that never runs. WARNINGS, not errors, for two
  // load-bearing reasons: the resolved configuration materializes both keys
  // unconditionally (an error would fail every off run's r1 round-trip),
  // and existing case files carrying them must stay valid if the T8
  // default-flip evaluation makes disabled the default (plan §2.7).
  {
    const YAML::Node re = sub(root, "groundwater", "reallocation");
    if (re.IsDefined() && re.IsScalar() && re.as<std::string>() == "disabled") {
      if (sub(root, "groundwater", "reallocation_surplus").IsDefined()) {
        ctx.addWarning("groundwater.reallocation_surplus",
                       "has no effect with reallocation: disabled (nothing "
                       "is transferred)");
      }
      if (sub(root, "groundwater", "density_coupling").IsDefined() &&
          root["groundwater"]["density_coupling"]["reallocation_gradient"].IsDefined()) {
        ctx.addWarning("groundwater.density_coupling.reallocation_gradient",
                       "has no effect with reallocation: disabled");
      }
      const YAML::Node vars = sub(root, "output", "variables");
      if (vars.IsDefined() && vars.IsMap() && vars["groundwater"].IsDefined() &&
          vars["groundwater"].IsSequence()) {
        for (std::size_t n = 0; n < vars["groundwater"].size(); ++n) {
          const std::string var = vars["groundwater"][n].as<std::string>("");
          if (var == "realloc_created" || var == "realloc_dropped") {
            ctx.addWarning("output.variables.groundwater",
                           var + " is identically zero with reallocation: disabled");
          }
        }
      }
    }
  }

  // T11 (O17b): the harmonic face rule exists only in the consistent
  // predictor–corrector pair — the legacy corrector keeps its one-sided
  // quirk and cannot take a different mean. Since the S3-A18 default flip
  // an ABSENT predictor_corrector resolves to consistent, so only an
  // explicit `legacy` conflicts.
  {
    const YAML::Node fc = sub(root, "groundwater", "face_conductivity");
    const YAML::Node pcm = sub(root, "groundwater", "predictor_corrector");
    const bool harmonic = fc.IsDefined() && fc.IsScalar() &&
                          fc.as<std::string>() == "harmonic";
    const bool legacyPc = pcm.IsDefined() && pcm.IsScalar() &&
                          pcm.as<std::string>() == "legacy";
    if (harmonic && legacyPc) {
      ctx.addError("groundwater.face_conductivity",
                   "harmonic needs groundwater.predictor_corrector: consistent "
                   "(the legacy corrector keeps its one-sided face rule)");
    }
  }

  // transport.subsurface_update: consistent rewires the subsurface scalar
  // update onto the groundwater module's window volumes; without that
  // module there is no subsurface scalar to update.
  {
    const YAML::Node update = sub(root, "transport", "subsurface_update");
    if (update.IsDefined() && update.IsScalar() && update.as<std::string>() == "consistent" &&
        !gw) {
      ctx.addError("transport.subsurface_update",
                   "subsurface_update: consistent requires modules.groundwater: true");
    }
  }

  // transport.surface_update: conservative advects with the solve face
  // volumes, whose coupled-deposit book (T10 C3) reads the consistent-
  // exchange machinery (applied gain, vent windows); with groundwater the
  // subsurface update must therefore be consistent too (a legacy
  // subsurface exchange would deliver the water's volume without its
  // scalar — the sx5 heat-drain defect). Since the S3-A18 default flip an
  // ABSENT surface_update resolves to conservative (and an absent
  // subsurface_update beside it to consistent), so the rule fires on the
  // RESOLVED conservative scheme with an EXPLICIT legacy subsurface.
  {
    const YAML::Node update = sub(root, "transport", "surface_update");
    const YAML::Node subUpdate = sub(root, "transport", "subsurface_update");
    const bool conservative =
        !(update.IsDefined() && update.IsScalar()) ||
        update.as<std::string>("") == "conservative";
    const bool subLegacy = subUpdate.IsDefined() && subUpdate.IsScalar() &&
                           subUpdate.as<std::string>() == "legacy";
    if (conservative && gw && subLegacy) {
      ctx.addError("transport.surface_update",
                   "surface_update: conservative (the default) with "
                   "modules.groundwater requires transport.subsurface_update: "
                   "consistent (the coupled-deposit salt/heat books read the "
                   "consistent exchange machinery); to keep the legacy "
                   "subsurface update, set surface_update: legacy or "
                   "consistent explicitly");
    }
  }

  // density coupling requires a scalar to couple to (salinity and/or
  // temperature since v2 Q5); the per-term coefficients require their
  // scalar's module.
  {
    const YAML::Node dc = sub(root, "groundwater", "density_coupling");
    if (dc.IsDefined() && dc.IsMap()) {
      bool enabled = false;
      if (dc["enabled"].IsDefined() && dc["enabled"].IsScalar()) {
        nodeAsBool(dc["enabled"], enabled);
      }
      if (enabled && !tr && !tp) {
        ctx.addError("groundwater.density_coupling.enabled",
                     "density coupling requires modules.transport or "
                     "modules.temperature");
      }
      if (!tr && (dc["beta_saline"].IsDefined() ||
                  dc["beta_saline_viscosity"].IsDefined())) {
        ctx.addError("groundwater.density_coupling",
                     "the saline coefficients require modules.transport: true");
      }
      if (!tp && (dc["thermal_expansion"].IsDefined() ||
                  dc["reference_temperature"].IsDefined())) {
        ctx.addError("groundwater.density_coupling",
                     "the thermal coefficients require modules.temperature: true");
      }
    }
  }

  // terrain_layers selects between the two terrain-following vertical rules
  // (amendment A12); it has no meaning on the regular mesh.
  {
    const YAML::Node layers = sub(root, "domain", "terrain_layers");
    bool follow = false;
    const YAML::Node followNode = sub(root, "domain", "follow_terrain");
    if (followNode.IsDefined() && followNode.IsScalar()) {
      nodeAsBool(followNode, follow);
    }
    if (layers.IsDefined() && layers.IsScalar() && layers.as<std::string>() == "uniform" &&
        !follow) {
      ctx.addError("domain.terrain_layers",
                   "'uniform' requires domain.follow_terrain: true");
    }
  }

  // dt_min <= dt_init <= dt_max.
  {
    const YAML::Node ts = sub(root, "groundwater", "timestep");
    if (ts.IsDefined() && ts.IsMap() && ts["dt_init"].IsDefined() && ts["dt_min"].IsDefined() &&
        ts["dt_max"].IsDefined()) {
      const double dtInit = realOr(ts["dt_init"], 0.0);
      const double dtMin = realOr(ts["dt_min"], 0.0);
      const double dtMax = realOr(ts["dt_max"], 0.0);
      if (!(dtMin <= dtInit && dtInit <= dtMax)) {
        ctx.addError("groundwater.timestep",
                     "requires dt_min <= dt_init <= dt_max (got " + std::to_string(dtMin) +
                         " / " + std::to_string(dtInit) + " / " + std::to_string(dtMax) + ")");
      }
    }
  }

  // initial conditions per enabled module.
  {
    const YAML::Node ic = root["initial_conditions"];
    if (sw && !sub(ic, "surface", "eta").IsDefined()) {
      ctx.addError("initial_conditions.surface.eta",
                   "required when modules.surface_water is true");
    }
    if (gw && !sub(ic, "groundwater").IsDefined()) {
      ctx.addError("initial_conditions.groundwater",
                   "required when modules.groundwater is true");
    }
    if (tr) {
      if (sw && !sub(ic, "transport", "surface").IsDefined()) {
        ctx.addError("initial_conditions.transport.surface",
                     "required when transport and surface_water are enabled");
      }
      if (gw && !sub(ic, "transport", "groundwater").IsDefined()) {
        ctx.addError("initial_conditions.transport.groundwater",
                     "required when transport and groundwater are enabled");
      }
    }
    if (tp) {
      if (sw && !sub(ic, "temperature", "surface").IsDefined()) {
        ctx.addError("initial_conditions.temperature.surface",
                     "required when temperature and surface_water are enabled");
      }
      if (gw && !sub(ic, "temperature", "groundwater").IsDefined()) {
        ctx.addError("initial_conditions.temperature.groundwater",
                     "required when temperature and groundwater are enabled");
      }
    }
  }

  // soil.types name uniqueness and soil.map reference.
  {
    const YAML::Node soil = root["soil"];
    if (soil.IsDefined() && soil.IsMap()) {
      std::set<std::string> names;
      const YAML::Node types = soil["types"];
      if (types.IsDefined() && types.IsSequence()) {
        for (std::size_t n = 0; n < types.size(); ++n) {
          const YAML::Node name = sub(types[n], "name");
          if (name.IsDefined() && name.IsScalar()) {
            if (!names.insert(name.as<std::string>()).second) {
              ctx.addError("soil.types[" + std::to_string(n) + "].name",
                           "duplicate soil type name '" + name.as<std::string>() + "'");
            }
          }
        }
      }
      const YAML::Node mapConst = sub(soil, "map", "constant");
      if (mapConst.IsDefined() && mapConst.IsScalar() &&
          names.find(mapConst.as<std::string>()) == names.end()) {
        ctx.addError("soil.map.constant", "soil type '" + mapConst.as<std::string>() +
                                              "' is not defined in soil.types");
      }
    }
  }

  // Stage-3 T6 (plan §13.C): vegetation cross-field rules — the
  // transpiration sink lives in the root zone (groundwater) and its
  // potential comes from the ET block; with salinity the exact sink
  // factor needs the consistent subsurface update; feddes head ordering
  // h1 >= h2 >= h3 >= h4 with h1/h2 together or not at all.
  {
    const YAML::Node veg = root["vegetation"];
    if (veg.IsDefined() && veg.IsMap()) {
      if (!gw) {
        ctx.addError("vegetation", "needs modules.groundwater (the root zone)");
      }
      if (!root["evapotranspiration"].IsDefined()) {
        ctx.addError("vegetation",
                     "needs the evapotranspiration block (it supplies the "
                     "potential T_pot = (1 - tau) E_x(r_c))");
      }
      if (tr && gw) {
        const YAML::Node trans = root["transport"];
        const bool consistentSub =
            trans.IsDefined() && trans.IsMap() &&
            trans["subsurface_update"].IsDefined() &&
            trans["subsurface_update"].as<std::string>("") == "consistent";
        if (!consistentSub) {
          ctx.addError("vegetation",
                       "with salinity transport requires "
                       "transport.subsurface_update: consistent (the root "
                       "sink's exact salt factor rides the window volumes)");
        }
      }
      std::set<std::string> vnames;
      const YAML::Node types = veg["types"];
      if (types.IsDefined() && types.IsSequence()) {
        for (std::size_t n = 0; n < types.size(); ++n) {
          const std::string where = "vegetation.types[" + std::to_string(n) + "]";
          const YAML::Node name = sub(types[n], "name");
          if (name.IsDefined() && name.IsScalar() &&
              !vnames.insert(name.as<std::string>()).second) {
            ctx.addError(where + ".name", "duplicate vegetation type name '" +
                                              name.as<std::string>() + "'");
          }
          const YAML::Node fed = sub(types[n], "feddes");
          if (fed.IsDefined() && fed.IsMap()) {
            const bool has1 = fed["h1"].IsDefined();
            const bool has2 = fed["h2"].IsDefined();
            if (has1 != has2) {
              ctx.addError(where + ".feddes",
                           "h1 and h2 come together or not at all");
            }
            if (fed["h3"].IsDefined() && fed["h4"].IsDefined()) {
              const double h3 = fed["h3"].as<double>(0.0);
              const double h4 = fed["h4"].as<double>(0.0);
              double h2v = h3, h1v = h3;
              if (has1 && has2) {
                h1v = fed["h1"].as<double>(0.0);
                h2v = fed["h2"].as<double>(0.0);
              }
              if (!(h1v >= h2v && h2v >= h3 && h3 > h4)) {
                ctx.addError(where + ".feddes",
                             "needs h1 >= h2 >= h3 > h4 (heads in m)");
              }
            }
          }
        }
      }
      const YAML::Node mapConst = sub(veg, "map", "constant");
      if (mapConst.IsDefined() && mapConst.IsScalar() &&
          vnames.find(mapConst.as<std::string>()) == vnames.end()) {
        ctx.addError("vegetation.map.constant",
                     "vegetation type '" + mapConst.as<std::string>() +
                         "' is not defined in vegetation.types");
      }
    }
  }

  // Boundary conditions: uniqueness, kind/target/value compatibility, deps.
  {
    const YAML::Node bcs = root["boundary_conditions"];
    if (bcs.IsDefined() && bcs.IsSequence()) {
      std::set<std::string> names;
      for (std::size_t n = 0; n < bcs.size(); ++n) {
        const YAML::Node bc = bcs[n];
        const std::string where = "boundary_conditions[" + std::to_string(n) + "]";
        if (!bc.IsMap()) {
          continue;
        }
        const YAML::Node nameNode = bc["name"];
        if (nameNode.IsDefined() && nameNode.IsScalar() &&
            !names.insert(nameNode.as<std::string>()).second) {
          ctx.addError(where + ".name",
                       "duplicate boundary condition name '" + nameNode.as<std::string>() + "'");
        }
        const std::string target =
            (bc["target"].IsDefined() && bc["target"].IsScalar()) ? bc["target"].as<std::string>()
                                                                  : std::string();
        const std::string kind =
            (bc["kind"].IsDefined() && bc["kind"].IsScalar()) ? bc["kind"].as<std::string>()
                                                              : std::string();
        const bool surfaceTarget = (target == "surface");
        const bool gwTarget = !target.empty() && !surfaceTarget;

        if (surfaceTarget && !sw) {
          ctx.addError(where, "target surface requires modules.surface_water: true");
        }
        if (gwTarget && !gw) {
          ctx.addError(where, "target " + target + " requires modules.groundwater: true");
        }
        if ((kind == "eta" || kind == "discharge" || kind == "velocity" || kind == "outflow") &&
            gwTarget) {
          ctx.addError(where, "kind " + kind + " applies to target surface only");
        }
        if ((kind == "head" || kind == "flux") && surfaceTarget) {
          ctx.addError(where, "kind " + kind + " applies to groundwater targets only");
        }
        if (kind == "head" && target == "groundwater_top" && sw && gw) {
          // In coupled runs the coupler owns the top boundary: wet columns
          // take the Dirichlet surface depth, dry columns are seepage faces
          // (plan §10 P3, amendment A11). A configured flux adds the legacy
          // qtop source; a configured head has no coupled meaning.
          ctx.addError(where, "target groundwater_top must be kind flux in coupled "
                              "runs (the coupler owns the top head)");
        }
        const std::string scalarSel =
            (bc["scalar"].IsDefined() && bc["scalar"].IsScalar())
                ? bc["scalar"].as<std::string>()
                : std::string("salinity");
        if (bc["scalar"].IsDefined() && kind != "scalar_value") {
          ctx.addError(where + ".scalar",
                       "the scalar selector applies to kind scalar_value only");
        }
        if (kind == "scalar_value") {
          if (scalarSel == "temperature" && !tp) {
            ctx.addError(where,
                         "scalar: temperature requires modules.temperature: true");
          }
          if (scalarSel == "salinity" && !tr) {
            ctx.addError(where, "kind scalar_value requires modules.transport: true");
          }
          // Cell-pinning Dirichlet tops/bottoms are the Q5 temperature
          // form (g6/g8); the salinity rows of the §8.2 matrix stay
          // schema-rejected until a gate exercises them.
          if ((target == "groundwater_top" || target == "groundwater_bottom") &&
              scalarSel != "temperature") {
            ctx.addError(where,
                         "scalar_value on " + target + " is temperature-only in "
                         "v2.0 (no salinity gate exercises the pinned-cell form)");
          }
          // The coupled top face is coupler-owned (same reasoning as the
          // head-on-top rejection above): a pinned top cell under a live
          // exchange has no gated meaning, so the §8.2 matrix rejects the
          // cell loudly rather than accepting it unverified (Q7).
          if (target == "groundwater_top" && sw && gw) {
            ctx.addError(where,
                         "scalar_value on groundwater_top is rejected in coupled "
                         "runs (the coupler owns the top exchange; no gate "
                         "exercises a pinned top cell under coupling)");
          }
        }
        if (kind == "scalar_cauchy") {
          // The v2 Q4 zero-total-scalar-flux top condition (Geng & Boufadel
          // 2015 Eq. (7)). §8.2 matrix: groundwater_top × uncoupled is the
          // one verified cell; every other target/mode is rejected loudly
          // rather than accepted unverified.
          if (!tr) {
            ctx.addError(where, "kind scalar_cauchy requires modules.transport: true");
          }
          if (target != "groundwater_top") {
            ctx.addError(where, "kind scalar_cauchy applies to target groundwater_top only "
                                "(the zero-total-flux salt condition rides the top-face "
                                "water flux)");
          }
          if (sw && gw) {
            // In coupled runs the coupler owns the top scalar exchange (the
            // sseepage path); a Cauchy top has no coupled meaning yet.
            ctx.addError(where, "kind scalar_cauchy applies to uncoupled groundwater runs "
                                "only (the coupler owns the top scalar exchange)");
          }
        }
        // Stage-3 T2 (plan §13.A): inflow-closure cross-field rules.
        {
          const YAML::Node inflow = bc["inflow"];
          if (inflow.IsDefined()) {
            if (kind != "eta") {
              ctx.addError(where + ".inflow",
                           "applies to kind eta only (it sets the momentum of "
                           "the external reservoir behind a prescribed stage)");
            }
            if (inflow.IsMap() && inflow["velocity"].IsDefined() &&
                inflow["closure"].IsDefined() && inflow["closure"].IsScalar() &&
                inflow["closure"].as<std::string>() != "external") {
              ctx.addError(where + ".inflow.velocity",
                           "only valid with closure external (zero_gradient "
                           "takes the inflow momentum from the interior)");
            }
          }
        }
        const YAML::Node value = bc["value"];
        if (kind == "outflow" && value.IsDefined()) {
          ctx.addError(where + ".value",
                       "kind outflow takes no value (free outflow is transmissive)");
        }
        if (kind == "scalar_cauchy" && value.IsDefined()) {
          ctx.addError(where + ".value",
                       "kind scalar_cauchy takes no value (the condition is the homogeneous "
                       "zero-total-flux relation)");
        }
        if (kind != "outflow" && kind != "scalar_cauchy" && !kind.empty() &&
            !value.IsDefined()) {
          ctx.addError(where, "missing required key 'value'");
        }
        if (value.IsDefined() && value.IsMap()) {
          if (value["gravity"].IsDefined() && kind != "flux") {
            ctx.addError(where + ".value.gravity", "free drainage applies to kind flux only");
          }
          if (value["gravity"].IsDefined() && target != "groundwater_bottom") {
            // Legacy free drainage (bctype_GW code 3) is a bottom-face
            // behavior; a gravity outflow through the top or a side has no
            // defined physics (P2).
            ctx.addError(where + ".value.gravity",
                         "free drainage applies to target groundwater_bottom only");
          }
          bool gravityOn = false;
          if (value["gravity"].IsDefined() && value["gravity"].IsScalar() &&
              nodeAsBool(value["gravity"], gravityOn) && !gravityOn) {
            ctx.addError(where + ".value.gravity", "must be true when given");
          }
          if (value["hydrostatic"].IsDefined() && kind != "head") {
            ctx.addError(where + ".value.hydrostatic",
                         "hydrostatic values apply to kind head only");
          }
          if (value["hydrostatic"].IsDefined() && target != "groundwater_side") {
            // The hydrostatic ghost-head form reproduces the legacy side
            // condition (enforce_head_bc, groundwater.c:769-788); top and
            // bottom heads are prescribed pressure heads (P2).
            ctx.addError(where + ".value.hydrostatic",
                         "hydrostatic values apply to target groundwater_side only");
          }
        }
      }
    }
  }

  // Output variable groups require their module.
  {
    const YAML::Node vars = sub(root, "output", "variables");
    if (vars.IsDefined() && vars.IsMap()) {
      if (vars["surface"].IsDefined() && !sw) {
        ctx.addError("output.variables.surface", "requires modules.surface_water: true");
      }
      if (vars["groundwater"].IsDefined() && !gw) {
        ctx.addError("output.variables.groundwater", "requires modules.groundwater: true");
      }
      if (vars["transport"].IsDefined() && !tr) {
        ctx.addError("output.variables.transport", "requires modules.transport: true");
      }
      // seepage and seepage_rate are surface variables produced by the
      // coupling (seepage_rate: Stage-3 T4 W4.2).
      const YAML::Node surfaceVars = vars["surface"];
      if (surfaceVars.IsDefined() && surfaceVars.IsSequence() && !gw) {
        for (std::size_t n = 0; n < surfaceVars.size(); ++n) {
          if (!surfaceVars[n].IsScalar()) {
            continue;
          }
          const std::string name = surfaceVars[n].as<std::string>();
          if (name == "seepage" || name == "seepage_rate") {
            ctx.addError("output.variables.surface",
                         "variable '" + name + "' requires modules.groundwater: true");
          }
        }
      }
      // Stage-3 T9 (plan §13.C): the paired water store exists only with
      // the store (surface water + transport + surface_dryout: store; an
      // absent key resolves to store since T8).
      if (surfaceVars.IsDefined() && surfaceVars.IsSequence()) {
        const YAML::Node trans = root["transport"];
        const bool storeMode =
            tr && !(trans.IsDefined() && trans.IsMap() && trans["surface_dryout"].IsDefined() &&
                    trans["surface_dryout"].as<std::string>("") == "discard");
        for (std::size_t n = 0; n < surfaceVars.size(); ++n) {
          if (surfaceVars[n].IsScalar() &&
              surfaceVars[n].as<std::string>() == "water_store" && !storeMode) {
            ctx.addError("output.variables.surface",
                         "water_store needs transport.surface_dryout: store "
                         "(with modules.surface_water and modules.transport)");
          }
        }
      }
    }
  }

  // Monitors: unique names, in-range indices, variables of enabled modules.
  {
    const YAML::Node monitors = sub(root, "output", "monitors");
    if (monitors.IsDefined() && monitors.IsSequence()) {
      const long nx = static_cast<long>(realOr(sub(root, "domain", "nx"), 0.0));
      const long ny = static_cast<long>(realOr(sub(root, "domain", "ny"), 0.0));
      std::set<std::string> names;
      for (std::size_t n = 0; n < monitors.size(); ++n) {
        const YAML::Node mon = monitors[n];
        const std::string where = "output.monitors[" + std::to_string(n) + "]";
        if (!mon.IsMap()) {
          continue;
        }
        const YAML::Node nameNode = mon["name"];
        if (nameNode.IsDefined() && nameNode.IsScalar() &&
            !names.insert(nameNode.as<std::string>()).second) {
          ctx.addError(where + ".name",
                       "duplicate monitor name '" + nameNode.as<std::string>() + "'");
        }
        const long mi = static_cast<long>(realOr(mon["i"], -1.0));
        const long mj = static_cast<long>(realOr(mon["j"], -1.0));
        if (nx > 0 && (mi < 0 || mi >= nx)) {
          ctx.addError(where + ".i", "outside the domain (nx = " + std::to_string(nx) + ")");
        }
        if (ny > 0 && (mj < 0 || mj >= ny)) {
          ctx.addError(where + ".j", "outside the domain (ny = " + std::to_string(ny) + ")");
        }
        const YAML::Node vars = mon["variables"];
        if (vars.IsDefined() && vars.IsSequence()) {
          for (std::size_t m = 0; m < vars.size(); ++m) {
            if (!vars[m].IsScalar()) {
              continue;
            }
            const std::string var = vars[m].as<std::string>();
            const auto& sv = surfaceVariables();
            const auto& gv = groundwaterVariables();
            const auto& tv = transportVariables();
            const bool isSurface = std::find(sv.begin(), sv.end(), var) != sv.end();
            const bool isGw = std::find(gv.begin(), gv.end(), var) != gv.end();
            const bool isTr = std::find(tv.begin(), tv.end(), var) != tv.end();
            if ((isSurface && !sw) || (isGw && !gw) || (isTr && !tr)) {
              ctx.addError(where + ".variables", "variable '" + var +
                                                     "' belongs to a module that is not enabled");
            }
          }
        }
      }
    }
  }

  // restart.enabled requires an existing file.
  {
    const YAML::Node restart = root["restart"];
    bool enabled = false;
    if (restart.IsDefined() && restart.IsMap() && restart["enabled"].IsDefined() &&
        restart["enabled"].IsScalar() && nodeAsBool(restart["enabled"], enabled) && enabled) {
      const YAML::Node file = restart["file"];
      const std::string fileValue =
          (file.IsDefined() && file.IsScalar()) ? file.as<std::string>() : std::string();
      if (fileValue.empty()) {
        ctx.addError("restart.file", "required when restart.enabled is true");
      } else {
        const std::filesystem::path p = std::filesystem::path(ctx.configDir) / fileValue;
        if (!std::filesystem::exists(p)) {
          ctx.addError("restart.file", "restart file '" + p.string() + "' does not exist");
        }
      }
    }
  }

  // solver.petsc_options_file must exist when set.
  {
    const YAML::Node file = sub(root, "solver", "petsc_options_file");
    if (file.IsDefined() && file.IsScalar()) {
      const std::string fileValue = file.as<std::string>();
      if (!fileValue.empty()) {
        const std::filesystem::path p = std::filesystem::path(ctx.configDir) / fileValue;
        if (!std::filesystem::exists(p)) {
          ctx.addError("solver.petsc_options_file",
                       "options file '" + p.string() + "' does not exist");
        }
      }
    }
  }
}

}  // namespace

void validateRoot(const YAML::Node& root, const std::string& configDir,
                  std::vector<std::string>& errors,
                  std::vector<std::string>& warnings) {
  Context ctx;
  ctx.configDir = configDir;
  ctx.errors = &errors;
  ctx.warnings = &warnings;

  if (!root.IsMap()) {
    ctx.addError("(document)", "top level must be a mapping");
    return;
  }
  const Spec schema = buildRootSchema();
  checkNode(root, schema, std::string(), ctx);
  crossChecks(root, ctx);
}

}  // namespace frehg::detail
