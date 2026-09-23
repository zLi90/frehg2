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
  real_t surfaceTemperatureC = 0.0;  ///< T_s [C] (prescribed until Q5)
  real_t pressureKpa = 0.0;          ///< [kPa]
  real_t airSpecificHumidity = 0.0;  ///< q_a [-]
  real_t windSpeed = 0.0;            ///< U [m/s]
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
  }

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
};

}  // namespace frehg::atm

#endif  // FREHG_ATM_METFORCING_HPP
