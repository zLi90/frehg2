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
using KeyedSpecs = std::vector<std::pair<std::string, Spec>>;

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
  static const std::vector<std::string> v = {"eta", "depth", "uu", "vv", "seepage"};
  return v;
}
const std::vector<std::string>& groundwaterVariables() {
  static const std::vector<std::string> v = {"hydraulic_head", "water_content", "qx", "qy", "qz"};
  return v;
}
const std::vector<std::string>& transportVariables() {
  static const std::vector<std::string> v = {"concentration", "concentration_surface"};
  return v;
}

std::vector<std::string> allVariables() {
  std::vector<std::string> v = surfaceVariables();
  v.insert(v.end(), groundwaterVariables().begin(), groundwaterVariables().end());
  v.insert(v.end(), transportVariables().begin(), transportVariables().end());
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
                       {"aev", real(false, -std::numeric_limits<double>::infinity(), 0.0)}});

  Spec bcEntry = map(
      {{"name", text(true)},
       {"region", map({{"polygon", sequence(point2(), true, 3)}}, true)},
       {"target", enumeration(true, {"surface", "groundwater_top", "groundwater_bottom",
                                     "groundwater_side"})},
       {"kind",
        enumeration(true,
                    {"eta", "discharge", "velocity", "outflow", "head", "flux", "scalar_value"})},
       // value is required for every kind except outflow (cross-field check).
       {"value", oneOf({{"constant", real(false)},
                        {"series", map({{"file", path(true)}})},
                        {"gravity", boolean(false)},
                        {"hydrostatic", map({{"eta", real(true)}})}},
                       false)}});

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
                      {"decomposition", map({{"mpi_nx", intOrAuto()}, {"mpi_ny", intOrAuto()}})}},
                     true)},
      {"time", map({{"dt", realPositive(true)},
                    {"t_start", integerSeconds(false, 0.0)},
                    {"t_end", realPositive(true)},
                    {"output_interval", integerSeconds(true, 1.0)}},
                   true)},
      {"modules", map({{"surface_water", boolean(false)},
                       {"groundwater", boolean(false)},
                       {"transport", boolean(false)}},
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
            {"wind", map({{"enabled", boolean(false)},
                          {"cd", realNonNegative(false)},
                          {"attenuation_depth", realPositive(false)},
                          {"north_angle", real(false)},
                          {"speed", seriesOrConstant(false)},
                          {"direction", seriesOrConstant(false)}})},
            // rainfall additionally accepts an exclusion region (no rain
            // inside; the legacy hardcoded last-row skip, made explicit);
            // the constant/series exclusivity is a cross-field check.
            {"rainfall", map({{"constant", real(false)},
                              {"series", map({{"file", path(true)}})},
                              {"exclude", map({{"polygon", sequence(point2(), true, 3)}})}})},
            {"evaporation", seriesOrConstant(false)}})},
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
            {"reallocation_surplus", enumeration(false, {"drop", "redistribute"})},
            {"density_coupling", map({{"enabled", boolean(false)}})}})},
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
                               {"groundwater", fileOrConstant(false)}})}})},
      {"boundary_conditions", sequence(bcEntry)},
      {"transport",
       map({{"scheme", map({{"advection", enumeration(false, {"upwind", "superbee"})}})},
            {"surface_diffusivity",
             map({{"x", realNonNegative(false)}, {"y", realNonNegative(false)}})},
            {"dispersion", map({{"longitudinal", realNonNegative(false)},
                                {"transverse", realNonNegative(false)},
                                {"molecular", realNonNegative(false)}})},
            {"bounds", map({{"min", real(false)}, {"max", real(false)}})}})},
      {"output",
       map({{"filename", text(true)},
            {"variables",
             map({{"surface", sequence(enumeration(false, surfaceVariables()), false, 1)},
                  {"groundwater",
                   sequence(enumeration(false, groundwaterVariables()), false, 1)},
                  {"transport", sequence(enumeration(false, transportVariables()), false, 1)}})},
            {"monitors", sequence(monitorEntry)},
            {"checkpoint", map({{"interval", integerSeconds(false, 0.0)}})}},
           true)},
      {"restart", map({{"enabled", boolean(false)}, {"file", text(false)}, {"time", real(false)}})},
      {"solver", map({{"petsc_options_file", text(false)}})},
      {"runtime", map({{"gpu_aware_mpi", enumeration(false, {"auto", "on", "off"})}})},
  });
}

// ---------------------------------------------------------------------------
// Validation walker
// ---------------------------------------------------------------------------

struct Context {
  std::string configDir;
  std::vector<std::string>* errors = nullptr;

  void addError(const std::string& where, const std::string& what) {
    errors->push_back((where.empty() ? "(document)" : where) + ": " + what);
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
                                   [&key](const auto& kv) { return kv.first == key; });
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

  if (!sw && !gw) {
    ctx.addError("modules", "at least one of surface_water/groundwater must be enabled");
  }
  if (tr && !sw && !gw) {
    ctx.addError("modules", "transport requires a flow module (surface_water or groundwater)");
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

  // density coupling requires transport.
  {
    const YAML::Node dc = sub(sub(root, "groundwater", "density_coupling"), "enabled");
    bool enabled = false;
    if (dc.IsDefined() && dc.IsScalar() && nodeAsBool(dc, enabled) && enabled && !tr) {
      ctx.addError("groundwater.density_coupling.enabled",
                   "density coupling requires modules.transport: true");
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
        if (kind == "scalar_value" && !tr) {
          ctx.addError(where, "kind scalar_value requires modules.transport: true");
        }
        const YAML::Node value = bc["value"];
        if (kind == "outflow" && value.IsDefined()) {
          ctx.addError(where + ".value",
                       "kind outflow takes no value (free outflow is transmissive)");
        }
        if (kind != "outflow" && !kind.empty() && !value.IsDefined()) {
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
      // seepage is a surface variable produced by the coupling.
      const YAML::Node surfaceVars = vars["surface"];
      if (surfaceVars.IsDefined() && surfaceVars.IsSequence() && !gw) {
        for (std::size_t n = 0; n < surfaceVars.size(); ++n) {
          if (surfaceVars[n].IsScalar() && surfaceVars[n].as<std::string>() == "seepage") {
            ctx.addError("output.variables.surface",
                         "variable 'seepage' requires modules.groundwater: true");
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
                  std::vector<std::string>& errors) {
  Context ctx;
  ctx.configDir = configDir;
  ctx.errors = &errors;

  if (!root.IsMap()) {
    ctx.addError("(document)", "top level must be a mapping");
    return;
  }
  const Spec schema = buildRootSchema();
  checkNode(root, schema, std::string(), ctx);
  crossChecks(root, ctx);
}

}  // namespace frehg::detail
