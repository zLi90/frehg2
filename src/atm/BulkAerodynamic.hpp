/// \file BulkAerodynamic.hpp
/// \brief The bulk-aerodynamic (mass-transfer) vapor-flux formulas
///        (v2 plan §3.2; one module, several consumers: open-water and
///        soil evaporation in Q4, the surface latent-heat term in Q5).
///
/// Provenance: Geng & Boufadel (2015) Eqs. (1)-(6) — the paper's own
/// formulation chain, which the V2-A13 derivation showed reproduces its
/// Table-1 values exactly (q_a = 2.896e-3 vs the stated 2.9e-3; E(0) =
/// 1.470e-7 m/s vs the digitized 1.442e-7). The aerodynamic resistance is
/// the Liu et al. (2006) fit named by the plan. The legacy aerodynamic
/// model hardwired these constants (v1 removed it for that, not its
/// physics); every input is configuration (core/Config.hpp).
///
/// All functions are host/device (KOKKOS_INLINE_FUNCTION) and pure; the
/// per-step evaluation of the configured series lives in MetForcing.hpp.

#ifndef FREHG_ATM_BULKAERODYNAMIC_HPP
#define FREHG_ATM_BULKAERODYNAMIC_HPP

#include "core/Types.hpp"

namespace frehg::atm {

inline constexpr real_t kGasConstantDryAir = 287.05;  ///< [J/(kg K)]
inline constexpr real_t kWaterDensity = 1000.0;       ///< [kg/m^3]

/// Saturated vapor pressure [kPa] at \p tempC [C] (Tetens; paper Eq. 5).
KOKKOS_INLINE_FUNCTION real_t tetensSaturationPressure(real_t tempC) {
  return 0.6108 * Kokkos::exp(17.27 * tempC / (tempC + 237.3));
}

/// Saturated specific humidity [-] at \p tempC [C], \p pressureKpa [kPa]
/// (paper Eq. 4).
KOKKOS_INLINE_FUNCTION real_t saturatedSpecificHumidity(real_t tempC,
                                                        real_t pressureKpa) {
  const real_t esat = tetensSaturationPressure(tempC);
  return 0.622 * esat / (pressureKpa - 0.376 * esat);
}

/// Aerodynamic resistance R_air [s/m] at wind speed \p windMs [m/s]
/// (Liu et al. 2006 fit; the plan §3.2 selectable form shipped in v2.0).
KOKKOS_INLINE_FUNCTION real_t aerodynamicResistance(real_t windMs) {
  return 94.909 * Kokkos::pow(windMs, static_cast<real_t>(-0.9036));
}

/// Dry-air density [kg/m^3] from the ideal gas law.
KOKKOS_INLINE_FUNCTION real_t airDensity(real_t tempC, real_t pressureKpa) {
  return pressureKpa * 1000.0 / (kGasConstantDryAir * (tempC + 273.15));
}

/// Soil relative humidity alpha_1(w_g) [-] from the surface volumetric
/// water content \p waterContent (paper Eq. 6, Barton / Lee & Pielke):
/// alpha_1 = min(1, 1.8 w_g / (w_g + 0.30)).
KOKKOS_INLINE_FUNCTION real_t soilRelativeHumidity(real_t waterContent) {
  const real_t a = 1.8 * waterContent / (waterContent + 0.30);
  return a < 1.0 ? a : 1.0;
}

/// Evaporative water flux [m/s of liquid water, positive = evaporation]
/// (Mahfouf & Noilhan form, plan §3.2):
///   E = (rho_a / R_air) (q_g - q_a) / rho_w,  q_g = alpha1 q_sat(T_s).
/// alpha1 = 1 is the open-water / potential rate. Negative results
/// (condensation) are returned as-is — clamping would hide a
/// humidity-gradient sign error from the mass audits (V2-A13 decision).
KOKKOS_INLINE_FUNCTION real_t evaporationRate(real_t surfaceTempC,
                                              real_t pressureKpa, real_t windMs,
                                              real_t airSpecificHumidity,
                                              real_t alpha1 = 1.0) {
  const real_t qGround = alpha1 * saturatedSpecificHumidity(surfaceTempC, pressureKpa);
  return airDensity(surfaceTempC, pressureKpa) / aerodynamicResistance(windMs) *
         (qGround - airSpecificHumidity) / kWaterDensity;
}

}  // namespace frehg::atm

#endif  // FREHG_ATM_BULKAERODYNAMIC_HPP
