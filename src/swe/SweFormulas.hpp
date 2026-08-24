/// \file SweFormulas.hpp
/// \brief Pointwise SWE closure formulas and their named legacy constants.
///
/// Every formula reproduces a specific legacy Frehg expression (provenance in
/// each function's documentation; plan Appendix B). They are host/device
/// callable so the unit suite exercises them against hand-computed values
/// (plan §8.1) with exactly the code the kernels run.

#ifndef FREHG_SWE_SWEFORMULAS_HPP
#define FREHG_SWE_SWEFORMULAS_HPP

#include "core/Types.hpp"

#include <Kokkos_Core.hpp>

namespace frehg::swe {

/// Lower edge of the advection damping ramp (legacy shallowwater.c:167).
inline constexpr real_t kCflRampLow = 0.5;
/// Upper edge of the advection damping ramp; advection is zeroed above it
/// (legacy shallowwater.c:166).
inline constexpr real_t kCflRampHigh = 0.7;
/// Air density used by the quadratic wind stress [kg/m^3] (legacy input
/// default rhoa; no benchmark varies it, so it is a named constant).
inline constexpr real_t kAirDensity = 1.225;
/// Water density dividing the wind stress [kg/m^3] (legacy rhow).
inline constexpr real_t kWaterDensity = 998.0;

/// Damp an explicit advection term by the local CFL number: zero at rest or
/// above the ramp, linearly reduced across [0.5, 0.7]
/// (legacy shallowwater.c:166-169).
KOKKOS_INLINE_FUNCTION real_t cflDampedAdvection(real_t advection, real_t velocity, real_t cfl) {
  if (velocity == 0.0 || cfl > kCflRampHigh) {
    return 0.0;
  }
  if (cfl > kCflRampLow) {
    return advection * (kCflRampHigh - cfl) / (kCflRampHigh - kCflRampLow);
  }
  return advection;
}

/// Manning drag coefficient with the thin-layer exponent switch: CD =
/// g n^2 / h^(2/3) below the thin-layer depth hD, g n^2 / h^(1/3) above it
/// (legacy update_drag_coef, shallowwater.c:977-994).
KOKKOS_INLINE_FUNCTION real_t manningDrag(real_t gravity, real_t manningN, real_t effectiveDepth,
                                          real_t thinLayerDepth) {
  const real_t expo = (effectiveDepth < thinLayerDepth) ? (2.0 / 3.0) : (1.0 / 3.0);
  return gravity * manningN * manningN / Kokkos::pow(effectiveDepth, expo);
}

/// Chezy drag coefficient CD = g / C^2 — the plan §5.8 additive extension.
/// No thin-layer exponent switch: that switch is Manning-specific.
KOKKOS_INLINE_FUNCTION real_t chezyDrag(real_t gravity, real_t chezyC) {
  return gravity / (chezyC * chezyC);
}

/// Point-implicit drag inversion factor
/// D = 1 / (0.5 dt CD |u| (Asz_face / V_face) + 1)
/// (legacy shallowwater.c:201-202).
KOKKOS_INLINE_FUNCTION real_t pointImplicitFactor(real_t dt, real_t dragCoef, real_t speed,
                                                  real_t faceAreaOverVolume) {
  return 1.0 / (0.5 * dt * dragCoef * speed * faceAreaOverVolume + 1.0);
}

/// Quadratic wind stress magnitude tau = rho_a Cw (W - u cos w - v sin w)^2
/// (legacy wind_source, shallowwater.c:267-269).
KOKKOS_INLINE_FUNCTION real_t windStress(real_t windDrag, real_t windSpeed, real_t u, real_t v,
                                         real_t omega) {
  const real_t relative = windSpeed - u * Kokkos::cos(omega) - v * Kokkos::sin(omega);
  return kAirDensity * windDrag * relative * relative;
}

/// Thin-layer wind-stress attenuation: full stress at face depths >= hD,
/// exponentially attenuated below, zero below hD/2
/// (legacy wind_source, shallowwater.c:271-286).
KOKKOS_INLINE_FUNCTION real_t windStressAttenuated(real_t tau, real_t faceDepth,
                                                   real_t thinLayerDepth, real_t attenuation) {
  if (faceDepth < thinLayerDepth) {
    if (faceDepth < 0.5 * thinLayerDepth) {
      return 0.0;
    }
    return tau * Kokkos::exp(attenuation * (faceDepth - thinLayerDepth) / thinLayerDepth);
  }
  return tau;
}

}  // namespace frehg::swe

#endif  // FREHG_SWE_SWEFORMULAS_HPP
