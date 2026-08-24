/// \file VanGenuchten.hpp
/// \brief Pointwise van Genuchten-Mualem retention and conductivity closures
///        (plan §3.1 item 2; legacy subroutines.c:280-515).
///
/// Every function reproduces the legacy van Genuchten branch (use_vg = 1
/// with the modified-vG switch off) — the only branch any benchmark
/// exercises; the modified-vG and exponential fallbacks are dropped per
/// plan §3.2. The conventions the legacy code
/// bakes in and this header preserves exactly:
///
///  - m = 1 - 1/n (plan §3.1 item 5);
///  - the **aev saturation-cutoff head applies to θ(h) only**
///    (subroutines.c:292 sits outside the modified-vG guard, which makes
///    aev a live parameter in b1/b2/b6): θ = θs for h > aev, while C(h) and
///    K(h) cut off at h > 0;
///  - θ is clamped to [θr, θs] after evaluation;
///  - K is capped at Ks;
///  - h(θ) guards θ - θr ≥ 1e-7 and returns 0 at saturation;
///  - dK/dθ evaluates saturation via the 0.9999·θs limiter
///    (subroutines.c:462-463).
///
/// All functions are scalar and callable from Kokkos kernels; test_vangenuchten
/// pins them against hand-computed values at 1e-12 relative (plan §8.1).

#ifndef FREHG_GW_VANGENUCHTEN_HPP
#define FREHG_GW_VANGENUCHTEN_HPP

#include "core/Types.hpp"

#include <Kokkos_Core.hpp>

namespace frehg::gw {

/// One cell's soil parameters staged for the pointwise closures.
struct VgSoil {
  real_t vga = 0.0;     ///< van Genuchten alpha [1/m] (legacy soil_a)
  real_t vgn = 0.0;     ///< van Genuchten n (legacy soil_n)
  real_t thetaS = 0.0;  ///< saturated water content (legacy wcs)
  real_t thetaR = 0.0;  ///< residual water content (legacy wcr)
  real_t aev = 0.0;     ///< saturation-cutoff (air-entry) head [m], <= 0
};

/// θ(h): water content from pressure head (legacy compute_wch,
/// subroutines.c:280-312). The aev cutoff applies here and only here.
KOKKOS_INLINE_FUNCTION real_t waterContentFromHead(const VgSoil& soil, real_t h) {
  const real_t m = 1.0 - 1.0 / soil.vgn;
  real_t wc;
  if (h > soil.aev) {
    wc = soil.thetaS;
  } else {
    const real_t s = Kokkos::pow(1.0 + Kokkos::pow(Kokkos::fabs(soil.vga * h), soil.vgn), -m);
    wc = soil.thetaR + (soil.thetaS - soil.thetaR) * s;
  }
  if (wc > soil.thetaS) {
    wc = soil.thetaS;
  } else if (wc < soil.thetaR) {
    wc = soil.thetaR;
  }
  return wc;
}

/// h(θ): pressure head from water content (legacy compute_hwc,
/// subroutines.c:314-346). Returns 0 at or above saturation; θ is floored
/// at θr + 1e-7 exactly as the legacy eps guard does.
KOKKOS_INLINE_FUNCTION real_t headFromWaterContent(const VgSoil& soil, real_t wc) {
  constexpr real_t eps = 1.0e-7;
  const real_t m = 1.0 - 1.0 / soil.vgn;
  if (wc - soil.thetaR < eps) {
    wc = soil.thetaR + eps;
  }
  if (wc < soil.thetaS) {
    const real_t ratio = (soil.thetaS - soil.thetaR) / (wc - soil.thetaR);
    return -(1.0 / soil.vga) *
           Kokkos::pow(Kokkos::pow(ratio, 1.0 / m) - 1.0, 1.0 / soil.vgn);
  }
  return 0.0;
}

/// C(h) = dθ/dh: specific moisture capacity (legacy compute_ch,
/// subroutines.c:349-376). Cuts off at h > 0 (not aev — legacy applies the
/// aev cutoff to C only under the dropped modified-vG branch).
KOKKOS_INLINE_FUNCTION real_t capacityFromHead(const VgSoil& soil, real_t h) {
  if (h > 0.0) {
    return 0.0;
  }
  const real_t m = 1.0 - 1.0 / soil.vgn;
  const real_t ah = Kokkos::fabs(soil.vga * h);
  const real_t nume = soil.vga * soil.vgn * m * (soil.thetaS - soil.thetaR) *
                      Kokkos::pow(ah, soil.vgn - 1.0);
  const real_t deno = Kokkos::pow(1.0 + Kokkos::pow(ah, soil.vgn), m + 1.0);
  return nume / deno;
}

/// K(h): unsaturated hydraulic conductivity by Mualem (legacy compute_K,
/// subroutines.c:378-405). \p ks is the saturated conductivity of whichever
/// cell the caller designates (the legacy face rules pass different cells'
/// Ks and vG parameters; see Predictor.cpp / Corrector.cpp). Cuts off at
/// h > 0 and is capped at ks.
KOKKOS_INLINE_FUNCTION real_t conductivityFromHead(const VgSoil& soil, real_t h, real_t ks) {
  const real_t m = 1.0 - 1.0 / soil.vgn;
  const real_t s = Kokkos::pow(1.0 + Kokkos::pow(Kokkos::fabs(soil.vga * h), soil.vgn), -m);
  const real_t bracket = 1.0 - Kokkos::pow(1.0 - Kokkos::pow(s, 1.0 / m), m);
  real_t keff = ks * Kokkos::pow(s, 0.5) * bracket * bracket;
  if (keff > ks) {
    keff = ks;
  }
  if (h > 0.0) {
    keff = ks;
  }
  return keff;
}

/// dK/dθ for the adaptive-step Courant limit (legacy compute_dKdwc,
/// subroutines.c:451-483, standard-vG branch). Callers evaluate this only
/// for unsaturated cells; the 0.9999·θs limiter keeps the expression finite
/// as θ approaches θs.
KOKKOS_INLINE_FUNCTION real_t conductivityDerivative(const VgSoil& soil, real_t wc, real_t ks) {
  if (wc > 0.9999 * soil.thetaS && wc < soil.thetaS) {
    wc = 0.9999 * soil.thetaS;
  }
  const real_t m = 1.0 - 1.0 / soil.vgn;
  const real_t s = (wc - soil.thetaR) / (soil.thetaS - soil.thetaR);
  const real_t term0 = Kokkos::pow(1.0 - Kokkos::pow(s, 1.0 / m), m);
  const real_t term1 = 0.5 * ks * Kokkos::pow(s, -0.5) * (1.0 - term0) * (1.0 - term0);
  const real_t term2 = 2.0 * ks * Kokkos::pow(s, (2.0 - m) / (2.0 * m)) * (1.0 - term0) *
                       Kokkos::pow(1.0 - Kokkos::pow(s, 1.0 / m), m - 1.0);
  return (term1 + term2) / (soil.thetaS - soil.thetaR);
}

}  // namespace frehg::gw

#endif  // FREHG_GW_VANGENUCHTEN_HPP
