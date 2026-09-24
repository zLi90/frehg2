/// \file MetForcing.hpp
/// \brief Host-side evaluation of the configured atmosphere block into the
///        plain met sample the bulk-aerodynamic kernels consume (v2 plan
///        §3.2). Each consumer module owns one MetForcing; series files
///        load once at construction, sample(t) is cheap per step.

#ifndef FREHG_ATM_METFORCING_HPP
#define FREHG_ATM_METFORCING_HPP

#include "atm/BulkAerodynamic.hpp"
#include "core/Config.hpp"
#include "core/TimeSeries.hpp"

namespace frehg::atm {

/// The met state at one instant, resolved to the quantities the kernels
/// take. The configured relative humidity (of q_sat at the air
/// temperature) is already folded into airSpecificHumidity.
struct MetSample {
  real_t airTemperatureC = 0.0;      ///< [C]
  /// T_s [C] for the Q4 evaporation consumers (0 when the optional key is
  /// absent — the schema requires it exactly when a Q4 consumer reads it;
  /// the Q5 heat exchange uses the local water temperature instead).
  real_t surfaceTemperatureC = 0.0;
  real_t pressureKpa = 0.0;          ///< [kPa]
  real_t airSpecificHumidity = 0.0;  ///< q_a [-]
  /// U [m/s], floored at atmosphere.wind_speed_floor (the GLM still-air
  /// lower bound on the bulk transfer functions, v2 Q5 — applied here so
  /// every bulk consumer sees it).
  real_t windSpeed = 0.0;
  real_t shortwave = 0.0;            ///< absorbed SW [W/m^2] (Q5)
  real_t longwaveIn = 0.0;           ///< incident LW [W/m^2] (Q5)
};

/// Configured met forcing, evaluated per step on the host.
class MetForcing {
 public:
  MetForcing() = default;
  /// Load the configured series (paths relative to the configuration).
  MetForcing(const AtmosphereConfig& atm, const FrehgConfig& config) : atm_(atm) {
    const auto load = [&config](const SeriesOrConstant& v, TimeSeries& out) {
      if (v.fromSeries) {
        out = TimeSeries::fromFile(config.resolvePath(v.file));
      }
    };
    load(atm_.airTemperature, airTemperature_);
    load(atm_.surfaceTemperature, surfaceTemperature_);
    load(atm_.pressure, pressure_);
    load(atm_.specificHumidity, specificHumidity_);
    load(atm_.relativeHumidity, relativeHumidity_);
    load(atm_.windSpeed, windSpeed_);
    load(atm_.shortwave, shortwave_);
    load(atm_.longwaveIn, longwaveIn_);
  }

  /// True when an atmosphere block was configured (sample() is only
  /// meaningful then).
  bool present() const { return atm_.present; }

  /// The met state at time \p t.
  MetSample sample(real_t t) const {
    const auto value = [t](const SeriesOrConstant& v, const TimeSeries& series) {
      return v.fromSeries ? series.value(t) : v.constant;
    };
    MetSample out;
    out.airTemperatureC = value(atm_.airTemperature, airTemperature_);
    out.surfaceTemperatureC = value(atm_.surfaceTemperature, surfaceTemperature_);
    out.pressureKpa = value(atm_.pressure, pressure_);
    out.windSpeed = value(atm_.windSpeed, windSpeed_);
    if (out.windSpeed < atm_.windSpeedFloor) {
      out.windSpeed = atm_.windSpeedFloor;
    }
    out.shortwave = value(atm_.shortwave, shortwave_);
    out.longwaveIn = value(atm_.longwaveIn, longwaveIn_);
    if (atm_.humidityIsRelative) {
      out.airSpecificHumidity =
          value(atm_.relativeHumidity, relativeHumidity_) *
          saturatedSpecificHumidity(out.airTemperatureC, out.pressureKpa);
    } else {
      out.airSpecificHumidity = value(atm_.specificHumidity, specificHumidity_);
    }
    return out;
  }

 private:
  AtmosphereConfig atm_;
  TimeSeries airTemperature_, surfaceTemperature_, pressure_;
  TimeSeries specificHumidity_, relativeHumidity_, windSpeed_;
  TimeSeries shortwave_, longwaveIn_;
};

}  // namespace frehg::atm

#endif  // FREHG_ATM_METFORCING_HPP
