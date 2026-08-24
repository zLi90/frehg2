/// \file Config.cpp
/// \brief Configuration loading, extraction, and effective-config echo.

#include "core/Config.hpp"

#include "core/Logger.hpp"

#include <yaml-cpp/yaml.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <utility>

namespace frehg {

namespace {

FileOrConstant extractFileOrConstant(const YAML::Node& node) {
  FileOrConstant out;
  if (node["file"].IsDefined()) {
    out.fromFile = true;
    out.file = node["file"].as<std::string>();
  } else {
    out.constant = node["constant"].as<real_t>();
  }
  return out;
}

SeriesOrConstant extractSeriesOrConstant(const YAML::Node& node) {
  SeriesOrConstant out;
  if (node["series"].IsDefined()) {
    out.fromSeries = true;
    out.file = node["series"]["file"].as<std::string>();
  } else {
    out.constant = node["constant"].as<real_t>();
  }
  return out;
}

template <class T>
T valueOr(const YAML::Node& node, const T& fallback) {
  return node.IsDefined() ? node.as<T>() : fallback;
}

int intOrAutoValue(const YAML::Node& node) {
  if (!node.IsDefined() || (node.IsScalar() && node.as<std::string>() == "auto")) {
    return 0;
  }
  return node.as<int>();
}

SimulationConfig extractSimulation(const YAML::Node& node) {
  SimulationConfig out;
  out.id = node["id"].as<std::string>();
  out.title = valueOr<std::string>(node["title"], std::string());
  return out;
}

DomainConfig extractDomain(const YAML::Node& node) {
  DomainConfig out;
  out.nx = node["nx"].as<int>();
  out.ny = node["ny"].as<int>();
  out.nz = node["nz"].as<int>();
  out.dx = node["dx"].as<real_t>();
  out.dy = node["dy"].as<real_t>();
  out.dz = node["dz"].as<real_t>();
  out.dzStretch = valueOr<real_t>(node["dz_stretch"], 1.0);
  out.bottomElevation = extractFileOrConstant(node["bottom_elevation"]);
  out.followTerrain = valueOr<bool>(node["follow_terrain"], false);
  out.terrainLayers =
      valueOr<std::string>(node["terrain_layers"], "scaled") == "uniform"
          ? DomainConfig::TerrainLayers::Uniform
          : DomainConfig::TerrainLayers::Scaled;
  if (node["decomposition"].IsDefined()) {
    out.decomposition.mpiNx = intOrAutoValue(node["decomposition"]["mpi_nx"]);
    out.decomposition.mpiNy = intOrAutoValue(node["decomposition"]["mpi_ny"]);
  }
  return out;
}

TimeConfig extractTime(const YAML::Node& node) {
  TimeConfig out;
  out.dt = node["dt"].as<real_t>();
  out.tStart = valueOr<real_t>(node["t_start"], 0.0);
  out.tEnd = node["t_end"].as<real_t>();
  out.outputInterval = node["output_interval"].as<real_t>();
  return out;
}

ModulesConfig extractModules(const YAML::Node& node) {
  ModulesConfig out;
  out.surfaceWater = valueOr<bool>(node["surface_water"], false);
  out.groundwater = valueOr<bool>(node["groundwater"], false);
  out.transport = valueOr<bool>(node["transport"], false);
  return out;
}

SurfaceWaterConfig extractSurfaceWater(const YAML::Node& node) {
  SurfaceWaterConfig out;
  out.gravity = valueOr<real_t>(node["gravity"], 9.81);
  const YAML::Node friction = node["friction"];
  out.friction.law = valueOr<std::string>(friction["law"], "manning") == "chezy"
                         ? FrictionConfig::Law::Chezy
                         : FrictionConfig::Law::Manning;
  out.friction.coefficient = extractFileOrConstant(friction["coefficient"]);
  out.friction.thinLayerDepth = valueOr<real_t>(friction["thin_layer_depth"], 0.1);
  if (node["viscosity"].IsDefined()) {
    out.viscosityX = valueOr<real_t>(node["viscosity"]["x"], 1.0e-6);
    out.viscosityY = valueOr<real_t>(node["viscosity"]["y"], 1.0e-6);
  }
  out.minDepth = node["min_depth"].as<real_t>();
  out.wettingFaceDepth = node["wetting_face_depth"].as<real_t>();
  const YAML::Node wind = node["wind"];
  if (wind.IsDefined()) {
    out.wind.enabled = valueOr<bool>(wind["enabled"], false);
    out.wind.cd = valueOr<real_t>(wind["cd"], 0.0013);
    out.wind.attenuationDepth = valueOr<real_t>(wind["attenuation_depth"], 5.0);
    out.wind.northAngle = valueOr<real_t>(wind["north_angle"], 0.0);
    if (wind["speed"].IsDefined()) {
      out.wind.speed = extractSeriesOrConstant(wind["speed"]);
    }
    if (wind["direction"].IsDefined()) {
      out.wind.direction = extractSeriesOrConstant(wind["direction"]);
    }
  }
  if (node["rainfall"].IsDefined()) {
    out.rainfall = extractSeriesOrConstant(node["rainfall"]);
    const YAML::Node exclude = node["rainfall"]["exclude"];
    if (exclude.IsDefined()) {
      const YAML::Node polygon = exclude["polygon"];
      for (std::size_t n = 0; n < polygon.size(); ++n) {
        out.rainfallExcludePolygon.push_back(
            {polygon[n][0].as<real_t>(), polygon[n][1].as<real_t>()});
      }
    }
  }
  if (node["evaporation"].IsDefined()) {
    out.evaporation = extractSeriesOrConstant(node["evaporation"]);
  }
  return out;
}

GroundwaterConfig extractGroundwater(const YAML::Node& node) {
  GroundwaterConfig out;
  out.scheme = valueOr<std::string>(node["scheme"], "pca");
  out.useFull3d = valueOr<bool>(node["use_full3d"], true);
  const YAML::Node ts = node["timestep"];
  out.timestep.dtInit = ts["dt_init"].as<real_t>();
  out.timestep.dtMin = ts["dt_min"].as<real_t>();
  out.timestep.dtMax = ts["dt_max"].as<real_t>();
  out.timestep.dqGrow = valueOr<real_t>(ts["dq_grow"], 0.01);
  out.timestep.dqShrink = valueOr<real_t>(ts["dq_shrink"], 0.02);
  out.timestep.courantMax = valueOr<real_t>(ts["courant_max"], 2.0);
  out.specificStorage = node["specific_storage"].as<real_t>();
  const std::string surplus =
      valueOr<std::string>(node["reallocation_surplus"], "drop");
  out.reallocationSurplus = (surplus == "redistribute")
                                ? GroundwaterConfig::ReallocationSurplus::Redistribute
                                : GroundwaterConfig::ReallocationSurplus::Drop;
  if (node["density_coupling"].IsDefined()) {
    out.densityCoupling.enabled = valueOr<bool>(node["density_coupling"]["enabled"], false);
  }
  return out;
}

SoilConfig extractSoil(const YAML::Node& node) {
  SoilConfig out;
  for (const YAML::Node& typeNode : node["types"]) {
    SoilType type;
    type.name = typeNode["name"].as<std::string>();
    type.ksx = typeNode["ksx"].as<real_t>();
    type.ksy = typeNode["ksy"].as<real_t>();
    type.ksz = typeNode["ksz"].as<real_t>();
    type.thetaS = typeNode["theta_s"].as<real_t>();
    type.thetaR = typeNode["theta_r"].as<real_t>();
    type.vgAlpha = typeNode["vg_alpha"].as<real_t>();
    type.vgN = typeNode["vg_n"].as<real_t>();
    type.aev = valueOr<real_t>(typeNode["aev"], 0.0);
    out.types.push_back(std::move(type));
  }
  const YAML::Node mapNode = node["map"];
  if (mapNode["file"].IsDefined()) {
    out.map.fromFile = true;
    out.map.file = mapNode["file"].as<std::string>();
  } else {
    out.map.constantName = mapNode["constant"].as<std::string>();
  }
  return out;
}

InitialConditionsConfig extractInitialConditions(const YAML::Node& node, const ModulesConfig& mods) {
  InitialConditionsConfig out;
  const YAML::Node surface = node["surface"];
  if (surface.IsDefined()) {
    out.surface.eta = extractFileOrConstant(surface["eta"]);
    if (surface["uu"].IsDefined()) {
      out.surface.hasUu = true;
      out.surface.uu = extractFileOrConstant(surface["uu"]);
    }
    if (surface["vv"].IsDefined()) {
      out.surface.hasVv = true;
      out.surface.vv = extractFileOrConstant(surface["vv"]);
    }
  }
  const YAML::Node gw = node["groundwater"];
  if (mods.groundwater && gw.IsDefined()) {
    if (gw["water_table"].IsDefined()) {
      out.groundwater.form = GroundwaterInitialConfig::Form::WaterTable;
      out.groundwater.value = extractFileOrConstant(gw["water_table"]);
    } else if (gw["head"].IsDefined()) {
      out.groundwater.form = GroundwaterInitialConfig::Form::Head;
      out.groundwater.value = extractFileOrConstant(gw["head"]);
    } else {
      out.groundwater.form = GroundwaterInitialConfig::Form::Moisture;
      out.groundwater.value = extractFileOrConstant(gw["moisture"]);
    }
  }
  const YAML::Node transport = node["transport"];
  if (transport.IsDefined()) {
    if (transport["surface"].IsDefined()) {
      out.transport.surface = extractFileOrConstant(transport["surface"]);
    }
    if (transport["groundwater"].IsDefined()) {
      out.transport.groundwater = extractFileOrConstant(transport["groundwater"]);
    }
  }
  return out;
}

BcTarget targetFromString(const std::string& s) {
  if (s == "surface") {
    return BcTarget::Surface;
  }
  if (s == "groundwater_top") {
    return BcTarget::GroundwaterTop;
  }
  if (s == "groundwater_bottom") {
    return BcTarget::GroundwaterBottom;
  }
  return BcTarget::GroundwaterSide;
}

BcKind kindFromString(const std::string& s) {
  if (s == "eta") {
    return BcKind::Eta;
  }
  if (s == "discharge") {
    return BcKind::Discharge;
  }
  if (s == "velocity") {
    return BcKind::Velocity;
  }
  if (s == "outflow") {
    return BcKind::Outflow;
  }
  if (s == "head") {
    return BcKind::Head;
  }
  if (s == "flux") {
    return BcKind::Flux;
  }
  return BcKind::ScalarValue;
}

std::vector<BoundaryConditionConfig> extractBoundaryConditions(const YAML::Node& node) {
  std::vector<BoundaryConditionConfig> out;
  if (!node.IsDefined()) {
    return out;
  }
  for (const YAML::Node& bcNode : node) {
    BoundaryConditionConfig bc;
    bc.name = bcNode["name"].as<std::string>();
    for (const YAML::Node& pt : bcNode["region"]["polygon"]) {
      bc.polygon.push_back({pt[0].as<real_t>(), pt[1].as<real_t>()});
    }
    bc.target = targetFromString(bcNode["target"].as<std::string>());
    bc.kind = kindFromString(bcNode["kind"].as<std::string>());
    const YAML::Node value = bcNode["value"];
    if (!value.IsDefined()) {
      // Only the outflow kind takes no value (schema cross-check).
      bc.value.form = BcValueConfig::Form::Constant;
    } else if (value["constant"].IsDefined()) {
      bc.value.form = BcValueConfig::Form::Constant;
      bc.value.constant = value["constant"].as<real_t>();
    } else if (value["series"].IsDefined()) {
      bc.value.form = BcValueConfig::Form::Series;
      bc.value.seriesFile = value["series"]["file"].as<std::string>();
    } else if (value["gravity"].IsDefined()) {
      bc.value.form = BcValueConfig::Form::Gravity;
    } else {
      bc.value.form = BcValueConfig::Form::Hydrostatic;
      bc.value.hydrostaticEta = value["hydrostatic"]["eta"].as<real_t>();
    }
    out.push_back(std::move(bc));
  }
  return out;
}

TransportConfig extractTransport(const YAML::Node& node) {
  TransportConfig out;
  if (node["scheme"].IsDefined()) {
    out.scheme.advection = valueOr<std::string>(node["scheme"]["advection"], "upwind") == "superbee"
                               ? TransportSchemeConfig::Advection::Superbee
                               : TransportSchemeConfig::Advection::Upwind;
  }
  if (node["surface_diffusivity"].IsDefined()) {
    out.surfaceDiffusivityX = valueOr<real_t>(node["surface_diffusivity"]["x"], 1.0e-10);
    out.surfaceDiffusivityY = valueOr<real_t>(node["surface_diffusivity"]["y"], 1.0e-10);
  }
  if (node["dispersion"].IsDefined()) {
    out.dispersionLongitudinal = valueOr<real_t>(node["dispersion"]["longitudinal"], 0.0);
    out.dispersionTransverse = valueOr<real_t>(node["dispersion"]["transverse"], 0.0);
    out.dispersionMolecular = valueOr<real_t>(node["dispersion"]["molecular"], 1.0e-10);
  }
  if (node["bounds"].IsDefined()) {
    out.boundMin = valueOr<real_t>(node["bounds"]["min"], 0.0);
    if (node["bounds"]["max"].IsDefined()) {
      out.hasBoundMax = true;
      out.boundMax = node["bounds"]["max"].as<real_t>();
    }
  }
  return out;
}

OutputConfig extractOutput(const YAML::Node& node) {
  OutputConfig out;
  out.filename = node["filename"].as<std::string>();
  const YAML::Node vars = node["variables"];
  if (vars.IsDefined()) {
    if (vars["surface"].IsDefined()) {
      out.surfaceVariables = vars["surface"].as<std::vector<std::string>>();
    }
    if (vars["groundwater"].IsDefined()) {
      out.groundwaterVariables = vars["groundwater"].as<std::vector<std::string>>();
    }
    if (vars["transport"].IsDefined()) {
      out.transportVariables = vars["transport"].as<std::vector<std::string>>();
    }
  }
  const YAML::Node monitors = node["monitors"];
  if (monitors.IsDefined()) {
    for (const YAML::Node& monNode : monitors) {
      MonitorConfig mon;
      mon.name = monNode["name"].as<std::string>();
      mon.i = monNode["i"].as<int>();
      mon.j = monNode["j"].as<int>();
      mon.variables = monNode["variables"].as<std::vector<std::string>>();
      out.monitors.push_back(std::move(mon));
    }
  }
  if (node["checkpoint"].IsDefined()) {
    out.checkpointInterval = valueOr<real_t>(node["checkpoint"]["interval"], 0.0);
  }
  return out;
}

RestartConfig extractRestart(const YAML::Node& node) {
  RestartConfig out;
  if (!node.IsDefined()) {
    return out;
  }
  out.enabled = valueOr<bool>(node["enabled"], false);
  out.file = valueOr<std::string>(node["file"], std::string());
  out.time = valueOr<real_t>(node["time"], 0.0);
  return out;
}

}  // namespace

std::string FrehgConfig::resolvePath(const std::string& relative) const {
  const std::filesystem::path p(relative);
  if (p.is_absolute()) {
    return relative;
  }
  return (std::filesystem::path(configDir) / p).string();
}

ValidationResult validateConfigFile(const std::string& path) {
  ValidationResult result;
  std::ifstream in(path);
  if (!in) {
    result.errors.push_back("cannot open configuration file '" + path + "'");
    return result;
  }
  std::stringstream buffer;
  buffer << in.rdbuf();

  YAML::Node root;
  try {
    root = YAML::Load(buffer.str());
  } catch (const YAML::Exception& e) {
    result.errors.push_back("YAML parse error: " + std::string(e.what()));
    return result;
  }
  const std::string dir = std::filesystem::path(path).parent_path().string();
  detail::validateRoot(root, dir.empty() ? "." : dir, result.errors);
  return result;
}

FrehgConfig loadConfig(const std::string& path) {
  const ValidationResult validation = validateConfigFile(path);
  if (!validation.ok()) {
    std::ostringstream all;
    all << "invalid configuration '" << path << "':";
    for (const std::string& err : validation.errors) {
      all << "\n  - " << err;
    }
    log::fatal(all.str());
  }

  std::ifstream in(path);
  std::stringstream buffer;
  buffer << in.rdbuf();
  const YAML::Node root = YAML::Load(buffer.str());

  FrehgConfig cfg;
  const std::string dir = std::filesystem::path(path).parent_path().string();
  cfg.configDir = dir.empty() ? "." : dir;
  cfg.rawText = buffer.str();

  cfg.simulation = extractSimulation(root["simulation"]);
  cfg.domain = extractDomain(root["domain"]);
  cfg.time = extractTime(root["time"]);
  cfg.modules = extractModules(root["modules"]);
  if (root["surface_water"].IsDefined()) {
    cfg.surfaceWater = extractSurfaceWater(root["surface_water"]);
  }
  if (root["groundwater"].IsDefined()) {
    cfg.groundwater = extractGroundwater(root["groundwater"]);
  }
  if (root["soil"].IsDefined()) {
    cfg.soil = extractSoil(root["soil"]);
  }
  if (root["coupling"].IsDefined()) {
    cfg.coupling.mode = valueOr<std::string>(root["coupling"]["mode"], "sync") == "subcycled"
                            ? CouplingConfig::Mode::Subcycled
                            : CouplingConfig::Mode::Sync;
  }
  if (root["initial_conditions"].IsDefined()) {
    cfg.initialConditions = extractInitialConditions(root["initial_conditions"], cfg.modules);
  }
  cfg.boundaryConditions = extractBoundaryConditions(root["boundary_conditions"]);
  if (root["transport"].IsDefined()) {
    cfg.transport = extractTransport(root["transport"]);
  }
  cfg.output = extractOutput(root["output"]);
  cfg.restart = extractRestart(root["restart"]);
  if (root["solver"].IsDefined()) {
    cfg.solver.petscOptionsFile =
        valueOr<std::string>(root["solver"]["petsc_options_file"], std::string());
  }
  if (root["runtime"].IsDefined()) {
    const std::string mode = valueOr<std::string>(root["runtime"]["gpu_aware_mpi"], "auto");
    cfg.runtime.gpuAwareMpi = mode == "on"    ? RuntimeConfig::GpuAwareMpi::On
                              : mode == "off" ? RuntimeConfig::GpuAwareMpi::Off
                                              : RuntimeConfig::GpuAwareMpi::Auto;
  }
  return cfg;
}

std::string describeConfig(const FrehgConfig& cfg) {
  std::ostringstream out;
  out << "effective configuration for '" << cfg.simulation.id << "'";
  if (!cfg.simulation.title.empty()) {
    out << " (" << cfg.simulation.title << ")";
  }
  out << "\n";
  out << "  domain: " << cfg.domain.nx << " x " << cfg.domain.ny << " x " << cfg.domain.nz
      << ", dx/dy/dz = " << cfg.domain.dx << "/" << cfg.domain.dy << "/" << cfg.domain.dz
      << ", dz_stretch = " << cfg.domain.dzStretch
      << (cfg.domain.followTerrain
              ? (cfg.domain.terrainLayers == DomainConfig::TerrainLayers::Uniform
                     ? ", terrain-following (uniform layers)"
                     : ", terrain-following")
              : "")
      << "\n";
  out << "  time: dt = " << cfg.time.dt << " s, t = [" << cfg.time.tStart << ", " << cfg.time.tEnd
      << "] s, output every " << cfg.time.outputInterval << " s\n";
  out << "  modules: surface_water=" << (cfg.modules.surfaceWater ? "on" : "off")
      << " groundwater=" << (cfg.modules.groundwater ? "on" : "off")
      << " transport=" << (cfg.modules.transport ? "on" : "off") << "\n";
  if (cfg.modules.surfaceWater) {
    out << "  surface_water: gravity = " << cfg.surfaceWater.gravity << ", friction = "
        << (cfg.surfaceWater.friction.law == FrictionConfig::Law::Chezy ? "chezy" : "manning")
        << ", thin_layer_depth = " << cfg.surfaceWater.friction.thinLayerDepth
        << ", min_depth = " << cfg.surfaceWater.minDepth
        << ", wetting_face_depth = " << cfg.surfaceWater.wettingFaceDepth
        << ", viscosity = " << cfg.surfaceWater.viscosityX << "/" << cfg.surfaceWater.viscosityY
        << ", wind = " << (cfg.surfaceWater.wind.enabled ? "on" : "off") << "\n";
  }
  if (cfg.modules.groundwater) {
    out << "  groundwater: scheme = " << cfg.groundwater.scheme
        << ", use_full3d = " << (cfg.groundwater.useFull3d ? "true" : "false")
        << ", Ss = " << cfg.groundwater.specificStorage << ", dtg = ["
        << cfg.groundwater.timestep.dtMin << ", " << cfg.groundwater.timestep.dtMax
        << "] (init " << cfg.groundwater.timestep.dtInit << ", dq " << cfg.groundwater.timestep.dqGrow
        << "/" << cfg.groundwater.timestep.dqShrink << ", Co_max "
        << cfg.groundwater.timestep.courantMax << ")"
        << ", density_coupling = " << (cfg.groundwater.densityCoupling.enabled ? "on" : "off")
        << ", reallocation_surplus = "
        << (cfg.groundwater.reallocationSurplus ==
                    GroundwaterConfig::ReallocationSurplus::Redistribute
                ? "redistribute"
                : "drop")
        << "\n";
    out << "  soil: " << cfg.soil.types.size() << " type(s), map = "
        << (cfg.soil.map.fromFile ? cfg.soil.map.file : cfg.soil.map.constantName) << "\n";
  }
  if (cfg.modules.surfaceWater && cfg.modules.groundwater) {
    out << "  coupling: mode = "
        << (cfg.coupling.mode == CouplingConfig::Mode::Sync ? "sync" : "subcycled") << "\n";
  }
  if (cfg.modules.transport) {
    out << "  transport: advection = "
        << (cfg.transport.scheme.advection == TransportSchemeConfig::Advection::Superbee
                ? "superbee"
                : "upwind")
        << ", dispersion l/t/m = " << cfg.transport.dispersionLongitudinal << "/"
        << cfg.transport.dispersionTransverse << "/" << cfg.transport.dispersionMolecular << "\n";
  }
  out << "  boundary_conditions: " << cfg.boundaryConditions.size() << "\n";
  out << "  output: " << cfg.output.filename << ", " << cfg.output.monitors.size()
      << " monitor(s), checkpoint interval = " << cfg.output.checkpointInterval << " s\n";
  out << "  runtime: gpu_aware_mpi = "
      << (cfg.runtime.gpuAwareMpi == RuntimeConfig::GpuAwareMpi::On    ? "on"
          : cfg.runtime.gpuAwareMpi == RuntimeConfig::GpuAwareMpi::Off ? "off"
                                                                       : "auto");
  return out.str();
}

}  // namespace frehg
