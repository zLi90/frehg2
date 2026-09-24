/// \file WindForcing.hpp
/// \brief Host-side per-step evaluation of the configured wind forcing
///        (v2 Q6, plan §5.2): the selectable Cd(U10) laws with cap, the
///        (u10, v10) component form, and the wrap-safe direction
///        interpolation of the legacy speed/direction form.
///
/// The momentum kernel consumes three per-step scalars (speed, angle from
/// +x in radians, drag coefficient); this class produces them from either
/// input form. Direction series interpolate on the circle: the sampled
/// directions are converted to unit vectors at load, the vectors
/// interpolate linearly in time, and the angle is recovered by atan2 — so
/// 350° → 10° passes through 0°/360°, never through 180° (the legacy
/// linear-in-degrees interpolation took the long way; the direction
/// series path was never exercised by any gate, and the convention is
/// unit-tested across the wrap). Chord interpolation is not constant
/// angular rate: between samples the recovered angle deviates from the
/// uniform sweep by O(delta^2/8) — under 0.05° for a 20° step — and a
/// 180° step between samples is degenerate (sample the series finely
/// enough that successive directions differ by well under 180°). The legacy constant-direction arithmetic
/// — omega = (direction + north_angle) · pi_legacy / 180 with the
/// truncated legacy literal — is preserved bitwise.

#ifndef FREHG_SWE_WINDFORCING_HPP
#define FREHG_SWE_WINDFORCING_HPP

#include "core/Config.hpp"
#include "core/TimeSeries.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace frehg::swe {

/// The legacy pi literal (wind_source, shallowwater.c:263).
inline constexpr real_t kLegacyPi = 3.1415926;

/// Cd(U10) [-] under the configured law (plan §5.2; published formulas).
/// Constant is the legacy Cw and takes NO cap (v1 behavior); every U10
/// law is capped at \p cap (default 3.5e-3, the ADCIRC cap practice).
inline real_t windDragCoefficient(WindConfig::DragLaw law, real_t u10, real_t cap,
                                  real_t constantCd) {
  real_t cd = constantCd;
  switch (law) {
    case WindConfig::DragLaw::Constant:
      return constantCd;
    case WindConfig::DragLaw::Garratt:  // Garratt (1977)
      cd = (0.75 + 0.067 * u10) * 1.0e-3;
      break;
    case WindConfig::DragLaw::SmithBanke:  // Smith & Banke (1975)
      cd = (0.63 + 0.066 * u10) * 1.0e-3;
      break;
    case WindConfig::DragLaw::Wu:  // Wu (1982)
      cd = (0.8 + 0.065 * u10) * 1.0e-3;
      break;
    case WindConfig::DragLaw::LargePond:  // Large & Pond (1981), piecewise
      cd = (u10 < 11.0) ? 1.2e-3 : (0.49 + 0.065 * u10) * 1.0e-3;
      break;
  }
  return cd < cap ? cd : cap;
}

/// Per-step wind sample the momentum kernel consumes.
struct WindSample {
  real_t speed = 0.0;   ///< |U10| [m/s]
  real_t omega = 0.0;   ///< direction, radians from +x (grid frame)
  real_t dragCd = 0.0;  ///< Cd under the configured law at this speed
};

/// Configured wind forcing, evaluated per step on the host.
class WindForcing {
 public:
  WindForcing() = default;

  /// Load the configured series (paths resolved by the caller-provided
  /// resolver, so this header stays free of FrehgConfig plumbing).
  template <typename Resolve>
  WindForcing(const WindConfig& wind, Resolve&& resolvePath) : cfg_(wind) {
    const auto load = [&](const SeriesOrConstant& v, TimeSeries& out) {
      if (v.fromSeries) {
        out = TimeSeries::fromFile(resolvePath(v.file));
      }
    };
    if (cfg_.componentForm) {
      load(cfg_.u10, u10_);
      load(cfg_.v10, v10_);
    } else {
      load(cfg_.speed, speed_);
      if (cfg_.direction.fromSeries) {
        // Circle-safe direction interpolation: sample the file's angles
        // into unit-vector component series (see the file comment).
        const TimeSeries raw = TimeSeries::fromFile(resolvePath(cfg_.direction.file));
        std::vector<real_t> s(raw.size());
        std::vector<real_t> c(raw.size());
        for (std::size_t n = 0; n < raw.size(); ++n) {
          const real_t rad = raw.values()[n] * kLegacyPi / 180.0;
          s[n] = std::sin(rad);
          c[n] = std::cos(rad);
        }
        dirSin_ = TimeSeries(raw.times(), std::move(s));
        dirCos_ = TimeSeries(raw.times(), std::move(c));
      }
    }
  }

  /// The wind state at time \p t.
  WindSample sample(real_t t) const {
    WindSample out;
    if (cfg_.componentForm) {
      const real_t u = cfg_.u10.fromSeries ? u10_.value(t) : cfg_.u10.constant;
      const real_t v = cfg_.v10.fromSeries ? v10_.value(t) : cfg_.v10.constant;
      out.speed = std::hypot(u, v);
      // atan2 is already "radians from +x"; the grid-to-north rotation
      // applies to compass input only, never to grid-frame components.
      out.omega = (out.speed > 0.0) ? std::atan2(v, u) : 0.0;
    } else {
      out.speed = cfg_.speed.fromSeries ? speed_.value(t) : cfg_.speed.constant;
      real_t direction = cfg_.direction.constant;
      if (cfg_.direction.fromSeries) {
        direction = std::atan2(dirSin_.value(t), dirCos_.value(t)) * 180.0 / kLegacyPi;
      }
      // The preserved legacy arithmetic (shallowwater.c:263-265).
      out.omega = (direction + cfg_.northAngle) * kLegacyPi / 180.0;
    }
    out.dragCd = windDragCoefficient(cfg_.law, out.speed, cfg_.cap, cfg_.cd);
    return out;
  }

 private:
  WindConfig cfg_;
  TimeSeries speed_, u10_, v10_;
  TimeSeries dirSin_, dirCos_;
};

}  // namespace frehg::swe

#endif  // FREHG_SWE_WINDFORCING_HPP
