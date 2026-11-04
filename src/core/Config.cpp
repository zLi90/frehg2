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
  if (node["periodic"].IsDefined()) {
    out.periodicX = valueOr<bool>(node["periodic"]["x"], false);
    out.periodicY = valueOr<bool>(node["periodic"]["y"], false);
  }
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
  out.temperature = valueOr<bool>(node["temperature"], false);
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
    const std::string law = valueOr<std::string>(wind["law"], "constant");
    out.wind.law = law == "garratt"       ? WindConfig::DragLaw::Garratt
                   : law == "smith-banke" ? WindConfig::DragLaw::SmithBanke
                   : law == "wu"          ? WindConfig::DragLaw::Wu
                   : law == "large-pond"  ? WindConfig::DragLaw::LargePond
                                          : WindConfig::DragLaw::Constant;
    out.wind.cd = valueOr<real_t>(wind["cd"], 0.0013);
    out.wind.cap = valueOr<real_t>(wind["cap"], 3.5e-3);
    out.wind.attenuationDepth = valueOr<real_t>(wind["attenuation_depth"], 5.0);
    out.wind.northAngle = valueOr<real_t>(wind["north_angle"], 0.0);
    out.wind.componentForm =
        wind["u10"].IsDefined() || wind["v10"].IsDefined();
    if (wind["u10"].IsDefined()) {
      out.wind.u10 = extractSeriesOrConstant(wind["u10"]);
    }
    if (wind["v10"].IsDefined()) {
      out.wind.v10 = extractSeriesOrConstant(wind["v10"]);
    }
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
    const YAML::Node evap = node["evaporation"];
    if (valueOr<std::string>(evap["mode"], "prescribed") == "bulk") {
      out.evapMode = SurfaceWaterConfig::EvapMode::Bulk;
    } else {
      out.evaporation = extractSeriesOrConstant(evap);
      const YAML::Node exclude = evap["exclude"];
      if (exclude.IsDefined()) {
        const YAML::Node polygon = exclude["polygon"];
        for (std::size_t n = 0; n < polygon.size(); ++n) {
          out.evaporationExcludePolygon.push_back(
              {polygon[n][0].as<real_t>(), polygon[n][1].as<real_t>()});
        }
      }
    }
  }
  return out;
}

AtmosphereConfig extractAtmosphere(const YAML::Node& node) {
  AtmosphereConfig out;
  out.present = true;
  out.airTemperature = extractSeriesOrConstant(node["air_temperature"]);
  if (node["surface_temperature"].IsDefined()) {
    out.hasSurfaceTemperature = true;
    out.surfaceTemperature = extractSeriesOrConstant(node["surface_temperature"]);
  }
  out.pressure = extractSeriesOrConstant(node["pressure"]);
  if (node["relative_humidity"].IsDefined()) {
    out.humidityIsRelative = true;
    out.relativeHumidity = extractSeriesOrConstant(node["relative_humidity"]);
  } else {
    out.specificHumidity = extractSeriesOrConstant(node["specific_humidity"]);
  }
  out.windSpeed = extractSeriesOrConstant(node["wind_speed"]);
  if (node["shortwave"].IsDefined()) {
    out.shortwave = extractSeriesOrConstant(node["shortwave"]);
  }
  if (node["longwave_in"].IsDefined()) {
    out.longwaveIn = extractSeriesOrConstant(node["longwave_in"]);
  }
  if (node["shortwave_incoming"].IsDefined()) {
    out.shortwaveIncoming = extractSeriesOrConstant(node["shortwave_incoming"]);
  }
  out.windSpeedFloor = valueOr<real_t>(node["wind_speed_floor"], 0.5);
  const YAML::Node aero = node["aerodynamic_resistance"];
  if (aero.IsDefined()) {
    const std::string law = valueOr<std::string>(aero["law"], "liu2006");
    out.aeroResistance = law == "fao56" ? AtmosphereConfig::AeroResistance::Fao56
                         : law == "log_profile"
                             ? AtmosphereConfig::AeroResistance::LogProfile
                             : AtmosphereConfig::AeroResistance::Liu2006;
    out.aeroRefHeight = valueOr<real_t>(aero["reference_height"], 2.0);
    out.aeroZ0m = valueOr<real_t>(aero["z0m"], 0.01);
    out.aeroZ0h = valueOr<real_t>(aero["z0h"], 0.001);
    out.aeroDisplacement = valueOr<real_t>(aero["displacement"], 0.0);
  }
  return out;
}

VegetationConfig extractVegetation(const YAML::Node& node) {
  VegetationConfig out;
  out.enabled = true;
  for (const YAML::Node& typeNode : node["types"]) {
    VegetationType type;
    type.name = typeNode["name"].as<std::string>();
    type.lai = typeNode["lai"].as<real_t>();
    type.extinction = typeNode["extinction"].as<real_t>();
    type.leafResistance = typeNode["leaf_resistance"].as<real_t>();
    type.rootDepth = typeNode["root_depth"].as<real_t>();
    type.profile = valueOr<std::string>(typeNode["root_profile"], "uniform") == "linear"
                       ? VegetationType::RootProfile::Linear
                       : VegetationType::RootProfile::Uniform;
    const YAML::Node fed = typeNode["feddes"];
    type.h3 = fed["h3"].as<real_t>();
    type.h4 = fed["h4"].as<real_t>();
    if (fed["h1"].IsDefined()) {
      type.hasH12 = true;
      type.h1 = fed["h1"].as<real_t>();
      type.h2 = fed["h2"].as<real_t>();
    }
    if (typeNode["albedo"].IsDefined()) {
      type.hasAlbedo = true;
      type.albedo = typeNode["albedo"].as<real_t>();
    }
    out.types.push_back(std::move(type));
  }
  const YAML::Node mapNode = node["map"];
  if (mapNode["file"].IsDefined()) {
    out.mapFromFile = true;
    out.mapFile = mapNode["file"].as<std::string>();
  } else {
    out.mapConstantName = mapNode["constant"].as<std::string>();
  }
  if (node["lai"].IsDefined()) {
    out.hasLaiFile = true;
    out.laiFile = node["lai"]["file"].as<std::string>();
  }
  return out;
}

EvapotranspirationConfig extractEvapotranspiration(const YAML::Node& node) {
  EvapotranspirationConfig out;
  out.enabled = true;
  const std::string mode = node["potential"].as<std::string>();
  out.potential = mode == "penman_monteith"
                      ? EvapotranspirationConfig::Potential::PenmanMonteith
                  : mode == "prescribed"
                      ? EvapotranspirationConfig::Potential::Prescribed
                      : EvapotranspirationConfig::Potential::Bulk;
  if (node["prescribed"].IsDefined()) {
    out.prescribed = extractSeriesOrConstant(node["prescribed"]);
  }
  out.openWater = valueOr<bool>(node["open_water"], true);
  const YAML::Node soil = node["soil"];
  if (soil.IsDefined()) {
    out.soil = true;
    const YAML::Node region = soil["region"];
    if (region.IsDefined() && region["polygon"].IsDefined()) {
      const YAML::Node polygon = region["polygon"];
      for (std::size_t n = 0; n < polygon.size(); ++n) {
        out.soilPolygon.push_back(
            {polygon[n][0].as<real_t>(), polygon[n][1].as<real_t>()});
      }
    }
    const std::string limiter = valueOr<std::string>(soil["limiter"], "alpha_barton");
    out.limiter = limiter == "beta_lee_pielke"
                      ? EvapotranspirationConfig::Limiter::BetaLeePielke
                  : limiter == "none" ? EvapotranspirationConfig::Limiter::None
                                      : EvapotranspirationConfig::Limiter::AlphaBarton;
  }
  const YAML::Node energy = node["energy"];
  if (energy.IsDefined()) {
    out.energy.albedoWater = valueOr<real_t>(energy["albedo_water"], 0.06);
    const YAML::Node g = energy["ground_heat_fraction"];
    if (g.IsDefined()) {
      out.energy.gFractionWater = valueOr<real_t>(g["water"], 0.0);
      out.energy.gFractionSoil = valueOr<real_t>(g["soil"], 0.1);
      out.energy.gFractionCanopy = valueOr<real_t>(g["canopy"], 0.1);
    }
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
  // Stage-3 T8 (plan §2.7): the struct initializer in Config.hpp is THE
  // default-flip constant — an absent key keeps it, so the DF evaluation
  // changes exactly one line.
  if (node["reallocation"].IsDefined()) {
    out.reallocation = node["reallocation"].as<std::string>() == "disabled"
                           ? GroundwaterConfig::Reallocation::Off
                           : GroundwaterConfig::Reallocation::Legacy;
  }
  // Stage-3 T11 (O17/S3-A16): same default-flip convention as reallocation.
  if (node["predictor_corrector"].IsDefined()) {
    out.predictorCorrector =
        node["predictor_corrector"].as<std::string>() == "consistent"
            ? GroundwaterConfig::PredictorCorrector::Consistent
            : GroundwaterConfig::PredictorCorrector::Legacy;
  }
  if (node["face_conductivity"].IsDefined()) {
    out.faceConductivity = node["face_conductivity"].as<std::string>() == "harmonic"
                               ? GroundwaterConfig::FaceConductivity::Harmonic
                               : GroundwaterConfig::FaceConductivity::Arithmetic;
  }
  const std::string surplus =
      valueOr<std::string>(node["reallocation_surplus"], "drop");
  out.reallocationSurplus = (surplus == "redistribute")
                                ? GroundwaterConfig::ReallocationSurplus::Redistribute
                                : GroundwaterConfig::ReallocationSurplus::Drop;
  if (node["density_coupling"].IsDefined()) {
    const YAML::Node dc = node["density_coupling"];
    out.densityCoupling.enabled = valueOr<bool>(dc["enabled"], false);
    out.densityCoupling.betaSaline = valueOr<real_t>(dc["beta_saline"], 0.000744);
    out.densityCoupling.betaSalineViscosity =
        valueOr<real_t>(dc["beta_saline_viscosity"], 0.0022);
    out.densityCoupling.thermalExpansion = valueOr<real_t>(dc["thermal_expansion"], 0.0);
    out.densityCoupling.referenceTemperature =
        valueOr<real_t>(dc["reference_temperature"], 20.0);
    out.densityCoupling.densityReallocationGradient =
        valueOr<std::string>(dc["reallocation_gradient"], "freshwater") == "density";
  }
  if (node["evaporation"].IsDefined()) {
    out.evaporation.enabled = true;  // mode "bulk" is the sole schema value
    const YAML::Node polygon = node["evaporation"]["region"]["polygon"];
    for (std::size_t n = 0; n < polygon.size(); ++n) {
      out.evaporation.polygon.push_back(
          {polygon[n][0].as<real_t>(), polygon[n][1].as<real_t>()});
    }
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
    if (typeNode["theta_fc"].IsDefined()) {
      type.hasThetaFc = true;
      type.thetaFc = typeNode["theta_fc"].as<real_t>();
    }
    if (typeNode["albedo"].IsDefined()) {
      type.hasAlbedo = true;
      type.albedo = typeNode["albedo"].as<real_t>();
    }
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
  const YAML::Node temperature = node["temperature"];
  if (temperature.IsDefined()) {
    if (temperature["surface"].IsDefined()) {
      out.temperature.surface = extractFileOrConstant(temperature["surface"]);
    }
    if (temperature["groundwater"].IsDefined()) {
      out.temperature.groundwater = extractFileOrConstant(temperature["groundwater"]);
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
  if (s == "scalar_cauchy") {
    return BcKind::ScalarCauchy;
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
    bc.scalarField = valueOr<std::string>(bcNode["scalar"], "salinity") == "temperature"
                         ? BcScalar::Temperature
                         : BcScalar::Salinity;
    const YAML::Node value = bcNode["value"];
    if (!value.IsDefined()) {
      // Only the outflow and scalar_cauchy kinds take no value (schema
      // cross-check).
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
    // Stage-3 T2 (plan §13.A): the inflow closure of a prescribed-stage
    // condition, with the external reservoir's velocity components.
    if (bcNode["inflow"].IsDefined()) {
      bc.hasInflowBlock = true;
      bc.inflowClosure =
          bcNode["inflow"]["closure"].as<std::string>() == "external"
              ? BcInflowClosure::External
              : BcInflowClosure::ZeroGradient;
      const auto extractComponent = [&](const char* key, BcValueConfig& outValue) {
        const YAML::Node comp = bcNode["inflow"]["velocity"][key];
        if (!comp.IsDefined()) {
          return;  // still-reservoir default (constant 0)
        }
        if (comp["constant"].IsDefined()) {
          outValue.form = BcValueConfig::Form::Constant;
          outValue.constant = comp["constant"].as<real_t>();
        } else {
          outValue.form = BcValueConfig::Form::Series;
          outValue.seriesFile = comp["series"]["file"].as<std::string>();
        }
      };
      if (bcNode["inflow"]["velocity"].IsDefined()) {
        extractComponent("u", bc.inflowU);
        extractComponent("v", bc.inflowV);
      }
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
  out.legacyEvapAllowance = valueOr<bool>(node["legacy_evap_allowance"], false);
  // Stage-3 T8 (plan §2.7-5): the struct initializer in Config.hpp is
  // the salt-store default-flip constant — an absent key keeps it.
  if (node["surface_dryout"].IsDefined()) {
    out.surfaceDryout = node["surface_dryout"].as<std::string>() == "store"
                            ? TransportConfig::SurfaceDryout::Store
                            : TransportConfig::SurfaceDryout::Delete;
  }
  {
    // Stage-3 T10: `conservative` implies the consistent semantics (one
    // more change on top: the solve's face volumes drive the advection).
    // Default flipped legacy -> conservative at the T10/T11 close (owner
    // decision, 2026-11-04; plan S3-A18); `legacy`/`consistent` stay the
    // golden-pinned opt-outs.
    const std::string su = valueOr<std::string>(node["surface_update"], "conservative");
    out.consistentSurfaceUpdate = su == "consistent" || su == "conservative";
    out.conservativeSurfaceUpdate = su == "conservative";
  }
  // The conservative surface update's coupled books need the consistent
  // subsurface exchange (the T10 schema rule), so when `conservative` is
  // the resolved surface scheme an ABSENT subsurface_update resolves to
  // consistent too; an explicit `legacy` beside it stays a schema error.
  out.consistentSubsurfaceUpdate =
      valueOr<std::string>(node["subsurface_update"],
                           out.conservativeSurfaceUpdate ? "consistent" : "legacy") ==
      "consistent";
  return out;
}

TemperatureConfig extractTemperature(const YAML::Node& node) {
  TemperatureConfig out;
  if (node["scheme"].IsDefined()) {
    out.scheme.advection =
        valueOr<std::string>(node["scheme"]["advection"], "upwind") == "superbee"
            ? TransportSchemeConfig::Advection::Superbee
            : TransportSchemeConfig::Advection::Upwind;
  }
  if (node["surface_diffusivity"].IsDefined()) {
    out.surfaceDiffusivityX = valueOr<real_t>(node["surface_diffusivity"]["x"], 1.0e-10);
    out.surfaceDiffusivityY = valueOr<real_t>(node["surface_diffusivity"]["y"], 1.0e-10);
  }
  out.thermalConductivity = valueOr<real_t>(node["thermal_conductivity"], 0.0);
  out.heatCapacityWater = valueOr<real_t>(node["heat_capacity_water"], 4.184e6);
  out.heatCapacitySolid = valueOr<real_t>(node["heat_capacity_solid"], 0.0);
  if (node["dispersivity"].IsDefined()) {
    out.dispersivityLongitudinal =
        valueOr<real_t>(node["dispersivity"]["longitudinal"], 0.0);
    out.dispersivityTransverse = valueOr<real_t>(node["dispersivity"]["transverse"], 0.0);
  }
  if (node["bounds"].IsDefined()) {
    if (node["bounds"]["min"].IsDefined()) {
      out.hasBoundMin = true;
      out.boundMin = node["bounds"]["min"].as<real_t>();
    }
    if (node["bounds"]["max"].IsDefined()) {
      out.hasBoundMax = true;
      out.boundMax = node["bounds"]["max"].as<real_t>();
    }
  }
  const YAML::Node exchange = node["surface_exchange"];
  if (exchange.IsDefined()) {
    out.surfaceExchange.mode =
        valueOr<std::string>(exchange["mode"], "equilibrium") == "bulk"
            ? SurfaceExchangeConfig::Mode::Bulk
            : SurfaceExchangeConfig::Mode::Equilibrium;
    if (exchange["equilibrium"].IsDefined()) {
      out.surfaceExchange.equilibriumTemperature =
          valueOr<real_t>(exchange["equilibrium"]["temperature"], 0.0);
      out.surfaceExchange.equilibriumCoefficient =
          valueOr<real_t>(exchange["equilibrium"]["coefficient"], 0.0);
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
    if (vars["temperature"].IsDefined()) {
      out.temperatureVariables = vars["temperature"].as<std::vector<std::string>>();
    }
    if (vars["atmosphere"].IsDefined()) {
      out.atmosphereVariables = vars["atmosphere"].as<std::vector<std::string>>();
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
  out.monitorInterval = valueOr<real_t>(node["monitor_interval"], 0.0);
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
  detail::validateRoot(root, dir.empty() ? "." : dir, result.errors, result.warnings);
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
  // Carried, not logged here: --resolve must keep its sentinel-fenced
  // stdout clean, so the run entry point (Simulation) logs them and the
  // run record stores them (Stage-3 W0.9).
  cfg.validationWarnings = validation.warnings;

  cfg.simulation = extractSimulation(root["simulation"]);
  cfg.domain = extractDomain(root["domain"]);
  cfg.time = extractTime(root["time"]);
  cfg.modules = extractModules(root["modules"]);
  if (root["atmosphere"].IsDefined()) {
    cfg.atmosphere = extractAtmosphere(root["atmosphere"]);
  }
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
  if (root["evapotranspiration"].IsDefined()) {
    cfg.evapotranspiration = extractEvapotranspiration(root["evapotranspiration"]);
  }
  if (root["vegetation"].IsDefined()) {
    cfg.vegetation = extractVegetation(root["vegetation"]);
  }
  if (root["transport"].IsDefined()) {
    cfg.transport = extractTransport(root["transport"]);
  } else if (cfg.modules.transport) {
    // An absent transport block takes the same resolved defaults the
    // extractor applies (S3-A18: conservative surface update, consistent
    // subsurface update beside it). A heat-only configuration (no
    // transport block) keeps the v2 Q5 legacy surface update its goldens
    // pin — the heat scalar opts in through an explicit transport block.
    cfg.transport = extractTransport(YAML::Node(YAML::NodeType::Map));
  }
  // The consistent subsurface update exists only with the groundwater
  // module (the schema rejects an EXPLICIT consistent without it); the
  // S3-A18 default must not synthesize it on surface-only runs, or the
  // resolved configuration embedded in the run record would not
  // revalidate.
  if (!cfg.modules.groundwater) {
    cfg.transport.consistentSubsurfaceUpdate = false;
  }
  if (root["temperature"].IsDefined()) {
    cfg.temperature = extractTemperature(root["temperature"]);
  }
  cfg.output = extractOutput(root["output"]);
  cfg.restart = extractRestart(root["restart"]);
  if (root["solver"].IsDefined()) {
    const auto parseSystem = [](const YAML::Node& node, SolverSystemConfig& out) {
      if (!node.IsDefined()) {
        return;
      }
      out.preconditioner = valueOr<std::string>(node["preconditioner"], out.preconditioner);
      out.matType = valueOr<std::string>(node["mat_type"], out.matType);
      out.rtol = valueOr<real_t>(node["rtol"], out.rtol);
      out.atol = valueOr<real_t>(node["atol"], out.atol);
      out.maxIterations = valueOr<int>(node["max_iterations"], out.maxIterations);
      out.reuseMaxSolves = valueOr<int>(node["reuse_max_solves"], out.reuseMaxSolves);
      out.reuseIterationFactor =
          valueOr<real_t>(node["reuse_iteration_factor"], out.reuseIterationFactor);
    };
    parseSystem(root["solver"]["surface"], cfg.solver.surface);
    parseSystem(root["solver"]["groundwater"], cfg.solver.groundwater);
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
      << " transport=" << (cfg.modules.transport ? "on" : "off")
      << " temperature=" << (cfg.modules.temperature ? "on" : "off") << "\n";
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
        << ", reallocation = "
        << (cfg.groundwater.reallocation == GroundwaterConfig::Reallocation::Off
                ? "disabled"
                : "legacy")
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
  if (cfg.evapotranspiration.enabled) {
    const EvapotranspirationConfig& et = cfg.evapotranspiration;
    out << "  evapotranspiration: potential = "
        << (et.potential == EvapotranspirationConfig::Potential::PenmanMonteith
                ? "penman_monteith"
            : et.potential == EvapotranspirationConfig::Potential::Prescribed
                ? "prescribed"
                : "bulk")
        << ", open_water = " << (et.openWater ? "on" : "off") << ", soil = ";
    if (et.soil) {
      out << "on (limiter = "
          << (et.limiter == EvapotranspirationConfig::Limiter::BetaLeePielke
                  ? "beta_lee_pielke"
              : et.limiter == EvapotranspirationConfig::Limiter::None ? "none"
                                                                      : "alpha_barton")
          << (et.soilPolygon.empty() ? ", region = all" : ", region = polygon") << ")";
    } else {
      out << "off";
    }
    if (cfg.vegetation.enabled) {
      out << ", vegetation = " << cfg.vegetation.types.size() << " type(s)";
    }
    if (cfg.atmosphere.present) {
      out << ", r_a = "
          << (cfg.atmosphere.aeroResistance == AtmosphereConfig::AeroResistance::Fao56
                  ? "fao56"
              : cfg.atmosphere.aeroResistance ==
                      AtmosphereConfig::AeroResistance::LogProfile
                  ? "log_profile"
                  : "liu2006");
    }
    out << "\n";
  }
  if (cfg.modules.transport) {
    out << "  transport: advection = "
        << (cfg.transport.scheme.advection == TransportSchemeConfig::Advection::Superbee
                ? "superbee"
                : "upwind")
        << ", dispersion l/t/m = " << cfg.transport.dispersionLongitudinal << "/"
        << cfg.transport.dispersionTransverse << "/" << cfg.transport.dispersionMolecular
        << ", surface_dryout = "
        << (cfg.transport.surfaceDryout == TransportConfig::SurfaceDryout::Store
                ? "store"
                : "discard")
        << ", surface_update = "
        << (cfg.transport.conservativeSurfaceUpdate ? "conservative"
            : cfg.transport.consistentSurfaceUpdate ? "consistent"
                                                    : "legacy")
        << "\n";
  }
  if (cfg.modules.temperature) {
    out << "  temperature: advection = "
        << (cfg.temperature.scheme.advection == TransportSchemeConfig::Advection::Superbee
                ? "superbee"
                : "upwind");
    if (cfg.modules.groundwater) {
      out << ", lambda_eff = " << cfg.temperature.thermalConductivity
          << " W/m/K, (rho c)_w/s = " << cfg.temperature.heatCapacityWater << "/"
          << cfg.temperature.heatCapacitySolid;
    }
    out << ", surface_exchange = "
        << (cfg.temperature.surfaceExchange.mode == SurfaceExchangeConfig::Mode::Bulk
                ? "bulk"
                : cfg.temperature.surfaceExchange.mode ==
                          SurfaceExchangeConfig::Mode::Equilibrium
                      ? "equilibrium"
                      : "none")
        << "\n";
  }
  out << "  boundary_conditions: " << cfg.boundaryConditions.size() << "\n";
  out << "  output: " << cfg.output.filename << ", " << cfg.output.monitors.size()
      << " monitor(s), monitor rows every "
      << (cfg.output.monitorInterval > 0.0
              ? std::to_string(static_cast<long>(cfg.output.monitorInterval)) + " s"
              : std::string("step"))
      << ", checkpoint interval = " << cfg.output.checkpointInterval << " s\n";
  out << "  runtime: gpu_aware_mpi = "
      << (cfg.runtime.gpuAwareMpi == RuntimeConfig::GpuAwareMpi::On    ? "on"
          : cfg.runtime.gpuAwareMpi == RuntimeConfig::GpuAwareMpi::Off ? "off"
                                                                       : "auto");
  return out.str();
}

namespace {

// ---------------------------------------------------------------------------
// Resolved-config serialization (v2 plan §2A): each emit* mirrors the
// corresponding extract* above exactly, so the pair is a fixed point.
// ---------------------------------------------------------------------------

YAML::Node emitFileOrConstant(const FileOrConstant& v) {
  YAML::Node n;
  if (v.fromFile) {
    n["file"] = v.file;
  } else {
    n["constant"] = v.constant;
  }
  return n;
}

YAML::Node emitSeriesOrConstant(const SeriesOrConstant& v) {
  YAML::Node n;
  if (v.fromSeries) {
    n["series"]["file"] = v.file;
  } else {
    n["constant"] = v.constant;
  }
  return n;
}

YAML::Node emitPolygon(const std::vector<std::array<real_t, 2>>& polygon) {
  YAML::Node n(YAML::NodeType::Sequence);
  for (const std::array<real_t, 2>& pt : polygon) {
    YAML::Node vertex(YAML::NodeType::Sequence);
    vertex.SetStyle(YAML::EmitterStyle::Flow);
    vertex.push_back(pt[0]);
    vertex.push_back(pt[1]);
    n.push_back(vertex);
  }
  return n;
}

YAML::Node emitSolverSystem(const SolverSystemConfig& s) {
  YAML::Node n;
  n["preconditioner"] = s.preconditioner;
  n["mat_type"] = s.matType;
  n["rtol"] = s.rtol;
  n["atol"] = s.atol;
  n["max_iterations"] = s.maxIterations;
  n["reuse_max_solves"] = s.reuseMaxSolves;
  n["reuse_iteration_factor"] = s.reuseIterationFactor;
  return n;
}

}  // namespace

std::string resolvedConfigYaml(const FrehgConfig& cfg) {
  YAML::Node root;

  root["simulation"]["id"] = cfg.simulation.id;
  if (!cfg.simulation.title.empty()) {
    root["simulation"]["title"] = cfg.simulation.title;
  }

  YAML::Node domain;
  domain["nx"] = cfg.domain.nx;
  domain["ny"] = cfg.domain.ny;
  domain["nz"] = cfg.domain.nz;
  domain["dx"] = cfg.domain.dx;
  domain["dy"] = cfg.domain.dy;
  domain["dz"] = cfg.domain.dz;
  domain["dz_stretch"] = cfg.domain.dzStretch;
  domain["bottom_elevation"] = emitFileOrConstant(cfg.domain.bottomElevation);
  domain["follow_terrain"] = cfg.domain.followTerrain;
  domain["periodic"]["x"] = cfg.domain.periodicX;
  domain["periodic"]["y"] = cfg.domain.periodicY;
  if (cfg.domain.followTerrain) {
    domain["terrain_layers"] =
        cfg.domain.terrainLayers == DomainConfig::TerrainLayers::Uniform ? "uniform" : "scaled";
  }
  // 0 means auto (MPI_Dims_create); the schema accepts "auto" or >= 1.
  domain["decomposition"]["mpi_nx"] =
      cfg.domain.decomposition.mpiNx > 0 ? YAML::Node(cfg.domain.decomposition.mpiNx)
                                         : YAML::Node("auto");
  domain["decomposition"]["mpi_ny"] =
      cfg.domain.decomposition.mpiNy > 0 ? YAML::Node(cfg.domain.decomposition.mpiNy)
                                         : YAML::Node("auto");
  root["domain"] = domain;

  YAML::Node time;
  time["dt"] = cfg.time.dt;
  time["t_start"] = cfg.time.tStart;
  time["t_end"] = cfg.time.tEnd;
  time["output_interval"] = cfg.time.outputInterval;
  root["time"] = time;

  root["modules"]["surface_water"] = cfg.modules.surfaceWater;
  root["modules"]["groundwater"] = cfg.modules.groundwater;
  root["modules"]["transport"] = cfg.modules.transport;
  root["modules"]["temperature"] = cfg.modules.temperature;

  if (cfg.modules.surfaceWater) {
    const SurfaceWaterConfig& sw = cfg.surfaceWater;
    YAML::Node node;
    node["gravity"] = sw.gravity;
    node["friction"]["law"] =
        sw.friction.law == FrictionConfig::Law::Chezy ? "chezy" : "manning";
    node["friction"]["coefficient"] = emitFileOrConstant(sw.friction.coefficient);
    node["friction"]["thin_layer_depth"] = sw.friction.thinLayerDepth;
    node["viscosity"]["x"] = sw.viscosityX;
    node["viscosity"]["y"] = sw.viscosityY;
    node["min_depth"] = sw.minDepth;
    node["wetting_face_depth"] = sw.wettingFaceDepth;
    YAML::Node wind;
    wind["enabled"] = sw.wind.enabled;
    wind["law"] = sw.wind.law == WindConfig::DragLaw::Garratt      ? "garratt"
                  : sw.wind.law == WindConfig::DragLaw::SmithBanke ? "smith-banke"
                  : sw.wind.law == WindConfig::DragLaw::Wu         ? "wu"
                  : sw.wind.law == WindConfig::DragLaw::LargePond  ? "large-pond"
                                                                   : "constant";
    // Mirror the schema's cross-field rules (cd belongs to the constant
    // law, cap to the U10 laws) so the resolved config revalidates — the
    // r1 round-trip requirement.
    if (sw.wind.law == WindConfig::DragLaw::Constant) {
      wind["cd"] = sw.wind.cd;
    } else {
      wind["cap"] = sw.wind.cap;
    }
    wind["attenuation_depth"] = sw.wind.attenuationDepth;
    wind["north_angle"] = sw.wind.northAngle;
    if (sw.wind.componentForm) {
      wind["u10"] = emitSeriesOrConstant(sw.wind.u10);
      wind["v10"] = emitSeriesOrConstant(sw.wind.v10);
    } else {
      wind["speed"] = emitSeriesOrConstant(sw.wind.speed);
      wind["direction"] = emitSeriesOrConstant(sw.wind.direction);
    }
    node["wind"] = wind;
    YAML::Node rainfall = emitSeriesOrConstant(sw.rainfall);
    if (!sw.rainfallExcludePolygon.empty()) {
      rainfall["exclude"]["polygon"] = emitPolygon(sw.rainfallExcludePolygon);
    }
    node["rainfall"] = rainfall;
    if (cfg.evapotranspiration.enabled) {
      // Stage-3 T5: the ET block replaces the legacy key (the schema
      // makes them exclusive — the r1 round trip must not re-emit it).
    } else if (sw.evapMode == SurfaceWaterConfig::EvapMode::Bulk) {
      node["evaporation"]["mode"] = "bulk";
    } else {
      YAML::Node evap = emitSeriesOrConstant(sw.evaporation);
      if (!sw.evaporationExcludePolygon.empty()) {
        evap["exclude"]["polygon"] = emitPolygon(sw.evaporationExcludePolygon);
      }
      node["evaporation"] = evap;
    }
    root["surface_water"] = node;
  }

  if (cfg.atmosphere.present) {
    YAML::Node atm;
    atm["air_temperature"] = emitSeriesOrConstant(cfg.atmosphere.airTemperature);
    if (cfg.atmosphere.hasSurfaceTemperature) {
      atm["surface_temperature"] = emitSeriesOrConstant(cfg.atmosphere.surfaceTemperature);
    }
    atm["pressure"] = emitSeriesOrConstant(cfg.atmosphere.pressure);
    if (cfg.atmosphere.humidityIsRelative) {
      atm["relative_humidity"] = emitSeriesOrConstant(cfg.atmosphere.relativeHumidity);
    } else {
      atm["specific_humidity"] = emitSeriesOrConstant(cfg.atmosphere.specificHumidity);
    }
    atm["wind_speed"] = emitSeriesOrConstant(cfg.atmosphere.windSpeed);
    atm["shortwave"] = emitSeriesOrConstant(cfg.atmosphere.shortwave);
    atm["longwave_in"] = emitSeriesOrConstant(cfg.atmosphere.longwaveIn);
    atm["shortwave_incoming"] = emitSeriesOrConstant(cfg.atmosphere.shortwaveIncoming);
    atm["wind_speed_floor"] = cfg.atmosphere.windSpeedFloor;
    YAML::Node aero;
    aero["law"] = cfg.atmosphere.aeroResistance == AtmosphereConfig::AeroResistance::Fao56
                      ? "fao56"
                  : cfg.atmosphere.aeroResistance ==
                          AtmosphereConfig::AeroResistance::LogProfile
                      ? "log_profile"
                      : "liu2006";
    aero["reference_height"] = cfg.atmosphere.aeroRefHeight;
    aero["z0m"] = cfg.atmosphere.aeroZ0m;
    aero["z0h"] = cfg.atmosphere.aeroZ0h;
    aero["displacement"] = cfg.atmosphere.aeroDisplacement;
    atm["aerodynamic_resistance"] = aero;
    root["atmosphere"] = atm;
  }

  if (cfg.evapotranspiration.enabled) {
    const EvapotranspirationConfig& et = cfg.evapotranspiration;
    YAML::Node node;
    node["potential"] =
        et.potential == EvapotranspirationConfig::Potential::PenmanMonteith
            ? "penman_monteith"
        : et.potential == EvapotranspirationConfig::Potential::Prescribed
            ? "prescribed"
            : "bulk";
    if (et.potential == EvapotranspirationConfig::Potential::Prescribed) {
      node["prescribed"] = emitSeriesOrConstant(et.prescribed);
    }
    node["open_water"] = et.openWater;
    if (et.soil) {
      YAML::Node soil;
      if (!et.soilPolygon.empty()) {
        soil["region"]["polygon"] = emitPolygon(et.soilPolygon);
      }
      soil["limiter"] =
          et.limiter == EvapotranspirationConfig::Limiter::BetaLeePielke
              ? "beta_lee_pielke"
          : et.limiter == EvapotranspirationConfig::Limiter::None ? "none"
                                                                  : "alpha_barton";
      node["soil"] = soil;
    }
    YAML::Node energy;
    energy["albedo_water"] = et.energy.albedoWater;
    energy["ground_heat_fraction"]["water"] = et.energy.gFractionWater;
    energy["ground_heat_fraction"]["soil"] = et.energy.gFractionSoil;
    energy["ground_heat_fraction"]["canopy"] = et.energy.gFractionCanopy;
    node["energy"] = energy;
    root["evapotranspiration"] = node;
  }

  if (cfg.vegetation.enabled) {
    YAML::Node node;
    for (const VegetationType& type : cfg.vegetation.types) {
      YAML::Node t;
      t["name"] = type.name;
      t["lai"] = type.lai;
      t["extinction"] = type.extinction;
      t["leaf_resistance"] = type.leafResistance;
      t["root_depth"] = type.rootDepth;
      t["root_profile"] =
          type.profile == VegetationType::RootProfile::Linear ? "linear" : "uniform";
      YAML::Node fed;
      if (type.hasH12) {
        fed["h1"] = type.h1;
        fed["h2"] = type.h2;
      }
      fed["h3"] = type.h3;
      fed["h4"] = type.h4;
      t["feddes"] = fed;
      if (type.hasAlbedo) {
        t["albedo"] = type.albedo;
      }
      node["types"].push_back(t);
    }
    if (cfg.vegetation.mapFromFile) {
      node["map"]["file"] = cfg.vegetation.mapFile;
    } else {
      node["map"]["constant"] = cfg.vegetation.mapConstantName;
    }
    if (cfg.vegetation.hasLaiFile) {
      node["lai"]["file"] = cfg.vegetation.laiFile;
    }
    root["vegetation"] = node;
  }

  if (cfg.modules.groundwater) {
    const GroundwaterConfig& gw = cfg.groundwater;
    YAML::Node node;
    node["scheme"] = gw.scheme;
    node["use_full3d"] = gw.useFull3d;
    node["timestep"]["dt_init"] = gw.timestep.dtInit;
    node["timestep"]["dt_min"] = gw.timestep.dtMin;
    node["timestep"]["dt_max"] = gw.timestep.dtMax;
    node["timestep"]["dq_grow"] = gw.timestep.dqGrow;
    node["timestep"]["dq_shrink"] = gw.timestep.dqShrink;
    node["timestep"]["courant_max"] = gw.timestep.courantMax;
    node["specific_storage"] = gw.specificStorage;
    node["reallocation"] =
        gw.reallocation == GroundwaterConfig::Reallocation::Off ? "disabled" : "legacy";
    node["predictor_corrector"] =
        gw.predictorCorrector == GroundwaterConfig::PredictorCorrector::Consistent
            ? "consistent"
            : "legacy";
    node["face_conductivity"] =
        gw.faceConductivity == GroundwaterConfig::FaceConductivity::Harmonic
            ? "harmonic"
            : "arithmetic";
    node["reallocation_surplus"] =
        gw.reallocationSurplus == GroundwaterConfig::ReallocationSurplus::Redistribute
            ? "redistribute"
            : "drop";
    node["density_coupling"]["enabled"] = gw.densityCoupling.enabled;
    node["density_coupling"]["reallocation_gradient"] =
        gw.densityCoupling.densityReallocationGradient ? "density" : "freshwater";
    // The per-term coefficients require their scalar's module (schema
    // cross-check) — mirror it or the r1 round-trip fails.
    if (cfg.modules.transport) {
      node["density_coupling"]["beta_saline"] = gw.densityCoupling.betaSaline;
      node["density_coupling"]["beta_saline_viscosity"] =
          gw.densityCoupling.betaSalineViscosity;
    }
    if (cfg.modules.temperature) {
      node["density_coupling"]["thermal_expansion"] = gw.densityCoupling.thermalExpansion;
      node["density_coupling"]["reference_temperature"] =
          gw.densityCoupling.referenceTemperature;
    }
    if (gw.evaporation.enabled) {
      node["evaporation"]["mode"] = "bulk";
      node["evaporation"]["region"]["polygon"] = emitPolygon(gw.evaporation.polygon);
    }
    root["groundwater"] = node;

    YAML::Node soil;
    for (const SoilType& type : cfg.soil.types) {
      YAML::Node t;
      t["name"] = type.name;
      t["ksx"] = type.ksx;
      t["ksy"] = type.ksy;
      t["ksz"] = type.ksz;
      t["theta_s"] = type.thetaS;
      t["theta_r"] = type.thetaR;
      t["vg_alpha"] = type.vgAlpha;
      t["vg_n"] = type.vgN;
      t["aev"] = type.aev;
      if (type.hasThetaFc) {
        t["theta_fc"] = type.thetaFc;
      }
      if (type.hasAlbedo) {
        t["albedo"] = type.albedo;
      }
      soil["types"].push_back(t);
    }
    if (cfg.soil.map.fromFile) {
      soil["map"]["file"] = cfg.soil.map.file;
    } else {
      soil["map"]["constant"] = cfg.soil.map.constantName;
    }
    root["soil"] = soil;
  }

  if (cfg.modules.surfaceWater && cfg.modules.groundwater) {
    root["coupling"]["mode"] =
        cfg.coupling.mode == CouplingConfig::Mode::Subcycled ? "subcycled" : "sync";
  }

  YAML::Node ic;
  if (cfg.modules.surfaceWater) {
    ic["surface"]["eta"] = emitFileOrConstant(cfg.initialConditions.surface.eta);
    if (cfg.initialConditions.surface.hasUu) {
      ic["surface"]["uu"] = emitFileOrConstant(cfg.initialConditions.surface.uu);
    }
    if (cfg.initialConditions.surface.hasVv) {
      ic["surface"]["vv"] = emitFileOrConstant(cfg.initialConditions.surface.vv);
    }
  }
  if (cfg.modules.groundwater) {
    const GroundwaterInitialConfig& gwIc = cfg.initialConditions.groundwater;
    const char* key = gwIc.form == GroundwaterInitialConfig::Form::WaterTable ? "water_table"
                      : gwIc.form == GroundwaterInitialConfig::Form::Head    ? "head"
                                                                             : "moisture";
    ic["groundwater"][key] = emitFileOrConstant(gwIc.value);
  }
  if (cfg.modules.transport) {
    if (cfg.modules.surfaceWater) {
      ic["transport"]["surface"] = emitFileOrConstant(cfg.initialConditions.transport.surface);
    }
    if (cfg.modules.groundwater) {
      ic["transport"]["groundwater"] =
          emitFileOrConstant(cfg.initialConditions.transport.groundwater);
    }
  }
  if (cfg.modules.temperature) {
    if (cfg.modules.surfaceWater) {
      ic["temperature"]["surface"] =
          emitFileOrConstant(cfg.initialConditions.temperature.surface);
    }
    if (cfg.modules.groundwater) {
      ic["temperature"]["groundwater"] =
          emitFileOrConstant(cfg.initialConditions.temperature.groundwater);
    }
  }
  if (ic.IsDefined() && ic.IsMap() && ic.size() > 0) {
    root["initial_conditions"] = ic;
  }

  if (!cfg.boundaryConditions.empty()) {
    YAML::Node bcs(YAML::NodeType::Sequence);
    for (const BoundaryConditionConfig& bc : cfg.boundaryConditions) {
      YAML::Node b;
      b["name"] = bc.name;
      b["region"]["polygon"] = emitPolygon(bc.polygon);
      b["target"] = bc.target == BcTarget::Surface             ? "surface"
                    : bc.target == BcTarget::GroundwaterTop    ? "groundwater_top"
                    : bc.target == BcTarget::GroundwaterBottom ? "groundwater_bottom"
                                                               : "groundwater_side";
      b["kind"] = bc.kind == BcKind::Eta          ? "eta"
                  : bc.kind == BcKind::Discharge  ? "discharge"
                  : bc.kind == BcKind::Velocity   ? "velocity"
                  : bc.kind == BcKind::Outflow    ? "outflow"
                  : bc.kind == BcKind::Head       ? "head"
                  : bc.kind == BcKind::Flux       ? "flux"
                  : bc.kind == BcKind::ScalarCauchy ? "scalar_cauchy"
                                                  : "scalar_value";
      if (bc.kind == BcKind::ScalarValue) {
        // Emit the selector only for temperature (salinity is the
        // default; emitting it unconditionally would still revalidate,
        // but the minimal form keeps hand-written configs fixed points).
        if (bc.scalarField == BcScalar::Temperature) {
          b["scalar"] = "temperature";
        }
      }
      // Outflow and scalar_cauchy take no value (schema cross-check:
      // transmissive / homogeneous zero-total-flux).
      if (bc.kind != BcKind::Outflow && bc.kind != BcKind::ScalarCauchy) {
        switch (bc.value.form) {
          case BcValueConfig::Form::Constant:
            b["value"]["constant"] = bc.value.constant;
            break;
          case BcValueConfig::Form::Series:
            b["value"]["series"]["file"] = bc.value.seriesFile;
            break;
          case BcValueConfig::Form::Gravity:
            b["value"]["gravity"] = true;
            break;
          case BcValueConfig::Form::Hydrostatic:
            b["value"]["hydrostatic"]["eta"] = bc.value.hydrostaticEta;
            break;
        }
      }
      // Stage-3 T2: every eta condition resolves its inflow closure (the
      // default materialized, like every other default — gate r1's
      // round-trip contract), with the velocity components only under
      // external and only when non-still (hand-written configs stay fixed
      // points).
      if (bc.kind == BcKind::Eta) {
        b["inflow"]["closure"] =
            bc.inflowClosure == BcInflowClosure::External ? "external" : "zero_gradient";
        const auto emitComponent = [&](const char* key, const BcValueConfig& comp) {
          if (comp.form == BcValueConfig::Form::Series) {
            b["inflow"]["velocity"][key]["series"]["file"] = comp.seriesFile;
          } else if (comp.constant != 0.0) {
            b["inflow"]["velocity"][key]["constant"] = comp.constant;
          }
        };
        if (bc.inflowClosure == BcInflowClosure::External) {
          emitComponent("u", bc.inflowU);
          emitComponent("v", bc.inflowV);
        }
      }
      bcs.push_back(b);
    }
    root["boundary_conditions"] = bcs;
  }

  if (cfg.modules.transport) {
    const TransportConfig& tr = cfg.transport;
    YAML::Node node;
    node["scheme"]["advection"] =
        tr.scheme.advection == TransportSchemeConfig::Advection::Superbee ? "superbee" : "upwind";
    node["surface_diffusivity"]["x"] = tr.surfaceDiffusivityX;
    node["surface_diffusivity"]["y"] = tr.surfaceDiffusivityY;
    node["dispersion"]["longitudinal"] = tr.dispersionLongitudinal;
    node["dispersion"]["transverse"] = tr.dispersionTransverse;
    node["dispersion"]["molecular"] = tr.dispersionMolecular;
    node["bounds"]["min"] = tr.boundMin;
    if (tr.hasBoundMax) {
      node["bounds"]["max"] = tr.boundMax;
    }
    node["legacy_evap_allowance"] = tr.legacyEvapAllowance;
    node["surface_update"] = tr.conservativeSurfaceUpdate ? "conservative"
                             : tr.consistentSurfaceUpdate ? "consistent"
                                                          : "legacy";
    node["subsurface_update"] = tr.consistentSubsurfaceUpdate ? "consistent" : "legacy";
    node["surface_dryout"] =
        tr.surfaceDryout == TransportConfig::SurfaceDryout::Store ? "store" : "discard";
    root["transport"] = node;
  }

  if (cfg.modules.temperature) {
    const TemperatureConfig& tp = cfg.temperature;
    YAML::Node node;
    node["scheme"]["advection"] =
        tp.scheme.advection == TransportSchemeConfig::Advection::Superbee ? "superbee" : "upwind";
    // The thermal parameters carry module dependencies (schema
    // cross-checks) — emit per enabled module so the output revalidates.
    if (cfg.modules.surfaceWater) {
      node["surface_diffusivity"]["x"] = tp.surfaceDiffusivityX;
      node["surface_diffusivity"]["y"] = tp.surfaceDiffusivityY;
    }
    node["heat_capacity_water"] = tp.heatCapacityWater;
    if (cfg.modules.groundwater) {
      node["thermal_conductivity"] = tp.thermalConductivity;
      node["heat_capacity_solid"] = tp.heatCapacitySolid;
      node["dispersivity"]["longitudinal"] = tp.dispersivityLongitudinal;
      node["dispersivity"]["transverse"] = tp.dispersivityTransverse;
    }
    if (tp.hasBoundMin) {
      node["bounds"]["min"] = tp.boundMin;
    }
    if (tp.hasBoundMax) {
      node["bounds"]["max"] = tp.boundMax;
    }
    if (tp.surfaceExchange.mode == SurfaceExchangeConfig::Mode::Equilibrium) {
      node["surface_exchange"]["mode"] = "equilibrium";
      node["surface_exchange"]["equilibrium"]["temperature"] =
          tp.surfaceExchange.equilibriumTemperature;
      node["surface_exchange"]["equilibrium"]["coefficient"] =
          tp.surfaceExchange.equilibriumCoefficient;
    } else if (tp.surfaceExchange.mode == SurfaceExchangeConfig::Mode::Bulk) {
      node["surface_exchange"]["mode"] = "bulk";
    }
    root["temperature"] = node;
  }

  YAML::Node output;
  output["filename"] = cfg.output.filename;
  if (!cfg.output.surfaceVariables.empty()) {
    output["variables"]["surface"] = cfg.output.surfaceVariables;
  }
  if (!cfg.output.groundwaterVariables.empty()) {
    output["variables"]["groundwater"] = cfg.output.groundwaterVariables;
  }
  if (!cfg.output.transportVariables.empty()) {
    output["variables"]["transport"] = cfg.output.transportVariables;
  }
  if (!cfg.output.temperatureVariables.empty()) {
    output["variables"]["temperature"] = cfg.output.temperatureVariables;
  }
  if (!cfg.output.atmosphereVariables.empty()) {
    output["variables"]["atmosphere"] = cfg.output.atmosphereVariables;
  }
  for (const MonitorConfig& mon : cfg.output.monitors) {
    YAML::Node m;
    m["name"] = mon.name;
    m["i"] = mon.i;
    m["j"] = mon.j;
    m["variables"] = mon.variables;
    output["monitors"].push_back(m);
  }
  output["monitor_interval"] = cfg.output.monitorInterval;
  output["checkpoint"]["interval"] = cfg.output.checkpointInterval;
  root["output"] = output;

  YAML::Node restart;
  restart["enabled"] = cfg.restart.enabled;
  if (cfg.restart.enabled) {
    restart["file"] = cfg.restart.file;
    restart["time"] = cfg.restart.time;
  }
  root["restart"] = restart;

  root["solver"]["surface"] = emitSolverSystem(cfg.solver.surface);
  root["solver"]["groundwater"] = emitSolverSystem(cfg.solver.groundwater);
  if (!cfg.solver.petscOptionsFile.empty()) {
    root["solver"]["petsc_options_file"] = cfg.solver.petscOptionsFile;
  }

  root["runtime"]["gpu_aware_mpi"] =
      cfg.runtime.gpuAwareMpi == RuntimeConfig::GpuAwareMpi::On    ? "on"
      : cfg.runtime.gpuAwareMpi == RuntimeConfig::GpuAwareMpi::Off ? "off"
                                                                   : "auto";

  YAML::Emitter emitter;
  emitter.SetDoublePrecision(17);
  emitter.SetFloatPrecision(9);
  emitter << root;
  return std::string(emitter.c_str()) + "\n";
}

}  // namespace frehg
