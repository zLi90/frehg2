/// \file EtFormulas.hpp
/// \brief The Stage-3 T5 evapotranspiration formulas (plan §9 W5.2):
///        selectable aerodynamic resistance, the bulk potential rate,
///        the Lee & Pielke beta limiter, and the supply cap.
///
/// Provenance: liu2006 is the v2.0 Liu et al. (2006) fit moved here from
/// BulkAerodynamic.hpp (which keeps a forwarding call so the Q4/Q5
/// consumers stay bitwise); fao56 is FAO-56 eq. 47's 208/u2 with the
/// standard measurement-height wind conversion u2 = U·4.87/ln(67.8 z − 5.42)
/// (Allen et al. 1998); log_profile is the neutral two-roughness form.
/// potentialBulk at r_a = liu2006, r_s = 0 is THE SAME EXPRESSION TREE as
/// atm::evaporationRate — asserted bitwise by AtmEt.BulkPotentialMatches…
/// (gate et0). The beta limiter is Lee & Pielke (1992) eq. 16; the supply
/// cap bounds an evaporative demand by the drainable store above θ_r and
/// a condensation flux by the pore room below the legacy 0.9999 θ_s mark
/// (the plan §3.1 named constant).
///
/// All functions are host/device (KOKKOS_INLINE_FUNCTION) and pure; the
/// per-step series evaluation lives in MetForcing.hpp, the per-column
/// application in gw/Evapotranspiration.cpp and swe/SurfaceSources.cpp.

#ifndef FREHG_ATM_ETFORMULAS_HPP
#define FREHG_ATM_ETFORMULAS_HPP

#include "atm/BulkAerodynamic.hpp"
#include "core/Config.hpp"
#include "core/Types.hpp"

namespace frehg::atm {

inline constexpr real_t kVonKarman = 0.41;  ///< von Karman constant [-]

/// Parameters of the log_profile / fao56 resistance laws (Config W5.1).
struct AeroParams {
  real_t refHeight = 2.0;     ///< measurement height z_ref [m]
  real_t z0m = 0.01;          ///< momentum roughness [m]
  real_t z0h = 0.001;         ///< scalar roughness [m]
  real_t displacement = 0.0;  ///< displacement height d [m]
};

/// Aerodynamic resistance r_a [s/m] under the selected law at wind speed
/// \p windMs [m/s] (already wind-floored by MetForcing).
KOKKOS_INLINE_FUNCTION real_t aeroResistance(
    AtmosphereConfig::AeroResistance law, real_t windMs, const AeroParams& p) {
  if (law == AtmosphereConfig::AeroResistance::Fao56) {
    // FAO-56 eq. 47: r_a for the reference crop reduces to 208/u2 with u2
    // the wind at 2 m, converted from z_ref by the log-profile ratio.
    const real_t u2 = windMs * 4.87 / Kokkos::log(67.8 * p.refHeight - 5.42);
    return 208.0 / u2;
  }
  if (law == AtmosphereConfig::AeroResistance::LogProfile) {
    const real_t z = p.refHeight - p.displacement;
    return Kokkos::log(z / p.z0m) * Kokkos::log(z / p.z0h) /
           (kVonKarman * kVonKarman * windMs);
  }
  // liu2006: forward to the v2.0 formula (BulkAerodynamic.hpp) so the two
  // sites cannot drift apart.
  return aerodynamicResistance(windMs);
}

/// The bulk potential vapor flux [m/s of liquid, positive = evaporation]:
///   E = rho_a(T_s, P) (alpha q_sat(T_s) - q_a) / (rho_w (r_a + r_s)).
/// At r_a = liu2006, r_s = 0, this is the v2.0 atm::evaporationRate
/// expression tree term for term (et0 asserts bitwise); r_s enters for
/// the T6 canopy resistance.
KOKKOS_INLINE_FUNCTION real_t potentialBulk(real_t surfaceTempC, real_t pressureKpa,
                                            real_t rA, real_t rS,
                                            real_t airSpecificHumidity,
                                            real_t alpha1 = 1.0) {
  const real_t qGround = alpha1 * saturatedSpecificHumidity(surfaceTempC, pressureKpa);
  return airDensity(surfaceTempC, pressureKpa) / (rA + rS) *
         (qGround - airSpecificHumidity) / kWaterDensity;
}

/// Lee & Pielke (1992) soil-moisture beta limiter [-]:
///   beta = 1 for theta >= theta_fc, else 0.25 (1 - cos(pi theta/theta_fc))^2.
KOKKOS_INLINE_FUNCTION real_t limiterBeta(real_t theta, real_t thetaFc) {
  if (theta >= thetaFc) {
    return 1.0;
  }
  const real_t c = 1.0 - Kokkos::cos(M_PI * theta / thetaFc);
  return 0.25 * c * c;
}

/// Supply cap on the applied top vapor flux [m/s]: an evaporative demand
/// takes at most the top cell's drainable store above theta_r per substep;
/// a condensation flux fills at most the pore room below 0.9999 theta_s
/// (the legacy effective-saturation mark — the same bound the coupled
/// supply-limited mode uses, Predictor.cpp).
KOKKOS_INLINE_FUNCTION real_t supplyCap(real_t demand, real_t theta, real_t thetaR,
                                        real_t thetaS, real_t dz, real_t dtg) {
  if (demand > 0.0) {
    const real_t store = Kokkos::fmax(0.0, theta - thetaR) * dz / dtg;
    return Kokkos::fmin(demand, store);
  }
  const real_t room = Kokkos::fmax(0.0, 0.9999 * thetaS - theta) * dz / dtg;
  return Kokkos::fmax(demand, -room);
}

/// Beer–Lambert canopy transmission tau = exp(-k LAI) (Stage-3 T6):
/// exactly 1.0 at lai = 0 (bare cells leave the ground rates untouched).
KOKKOS_INLINE_FUNCTION real_t beerLambert(real_t extinction, real_t lai) {
  if (lai == 0.0) {
    return 1.0;
  }
  return Kokkos::exp(-extinction * lai);
}

/// Canopy resistance r_c = r_l / (0.5 LAI) [s/m] (big-leaf; T6).
KOKKOS_INLINE_FUNCTION real_t canopyResistance(real_t leafResistance, real_t lai) {
  return leafResistance / (0.5 * lai);
}

/// Feddes (1978) root-uptake stress alpha_F(h) [-], piecewise linear:
/// 0 at/above h1, ramp h1->h2 (0->1), 1 on [h3, h2], ramp h3->h4
/// (1->0), 0 at/below h4; without h1/h2 (wetland species, O7), 1 for
/// every h >= h3.
KOKKOS_INLINE_FUNCTION real_t feddes(real_t h, real_t h1, real_t h2, real_t h3,
                                     real_t h4, bool hasH12) {
  if (h <= h4) {
    return 0.0;
  }
  if (h < h3) {
    return (h - h4) / (h3 - h4);
  }
  if (!hasH12 || h <= h2) {
    return 1.0;
  }
  if (h >= h1) {
    return 0.0;
  }
  return (h1 - h) / (h1 - h2);
}

/// The root-density integral over the depth interval [dA, dB] below the
/// column's top surface, for a root zone of depth L (T6 W6.2):
/// uniform b(d) = 1/L; linear b(d) = 2(L - d)/L^2 (triangular). Closed
/// forms; 0 when the interval starts at or below L.
KOKKOS_INLINE_FUNCTION real_t rootFractionIntegral(bool linearProfile, real_t L,
                                                   real_t dA, real_t dB) {
  if (dA >= L) {
    return 0.0;
  }
  const real_t b = Kokkos::fmin(dB, L);
  if (!linearProfile) {
    return (b - dA) / L;
  }
  const real_t fa = (2.0 * L * dA - dA * dA) / (L * L);
  const real_t fb = (2.0 * L * b - b * b) / (L * L);
  return fb - fa;
}

/// The theta_fc default: theta(h = -3.3 m) from the van Genuchten curve
/// (used when soil.types[].theta_fc is absent; W5.1).
KOKKOS_INLINE_FUNCTION real_t thetaFcFromRetention(real_t thetaR, real_t thetaS,
                                                   real_t vgAlpha, real_t vgN) {
  const real_t m = 1.0 - 1.0 / vgN;
  const real_t se = Kokkos::pow(1.0 + Kokkos::pow(vgAlpha * 3.3, vgN), -m);
  return thetaR + se * (thetaS - thetaR);
}

}  // namespace frehg::atm

#endif  // FREHG_ATM_ETFORMULAS_HPP
