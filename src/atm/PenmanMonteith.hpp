/// \file PenmanMonteith.hpp
/// \brief The Penman–Monteith potential-rate formula in specific-humidity
///        form (Stage-3 plan §11 W7.2). It shares the bulk chain's
///        functions — Tetens q_sat, airDensity, latentHeatOfVaporization,
///        kAirHeatCapacity, kWaterEmissivity — so the PM and bulk
///        potentials are commensurable: at the surface temperature that
///        balances the available energy, the two rates agree to second
///        order in (T_s − T_air) (the et0 part-3 unit pins it).
///
/// FAO-56 convention note (theory/evapotranspiration.md): FAO's reference
/// chain is the vapor-PRESSURE form with the psychrometric constant
/// γ = c_pa P / (0.622 λ) and tabulated radiation/γ conventions; this
/// q-form with the model's constants reproduces its Chapter-4 hourly
/// example to the measured δ_conv ≈ 2.6 % (gate et10's tolerance).
///
/// All functions are host/device (KOKKOS_INLINE_FUNCTION) and pure.

#ifndef FREHG_ATM_PENMANMONTEITH_HPP
#define FREHG_ATM_PENMANMONTEITH_HPP

#include "atm/BulkAerodynamic.hpp"

namespace frehg::atm {

/// Analytic d(q_sat)/dT [1/K] of saturatedSpecificHumidity at \p tempC
/// [C], \p pressureKpa [kPa]:
///   q_sat = 0.622 e / (P − 0.376 e),  de/dT = e · 17.27·237.3/(T+237.3)²
///   dq_sat/dT = 0.622 P (de/dT) / (P − 0.376 e)².
KOKKOS_INLINE_FUNCTION real_t dqsatdT(real_t tempC, real_t pressureKpa) {
  const real_t e = tetensSaturationPressure(tempC);
  const real_t t237 = tempC + 237.3;
  const real_t dedT = e * 17.27 * 237.3 / (t237 * t237);
  const real_t denom = pressureKpa - 0.376 * e;
  return 0.622 * dedT * (pressureKpa / (denom * denom));
}

/// Penman–Monteith potential rate E_x(r_s) [m/s of liquid water,
/// positive = evaporation] in specific-humidity form (plan §11):
///   λ ρ_w E_x = [ s·A + λ ρ_a (q_sat(T_air) − q_a)/r_a ] / [ s + 1 + r_s/r_a ]
/// with s = λ (dq_sat/dT)/c_pa at T_air, ρ_a = airDensity(T_air, P), and
/// the available energy
///   A = R_n − G,  R_n = (1 − albedo) R_s + ε (LW_in − σ T_K⁴),
///   G = gFraction · R_n,  T_K = T_air + 273.15,  ε = kWaterEmissivity.
/// \p shortwaveIncoming is INCIDENT R_s [W/m²] (atmosphere
/// .shortwave_incoming — not the Q5 absorbed `shortwave`); \p albedo and
/// \p gFraction are the calling channel's (water / soil type / vegetation
/// type). Negative results (condensation-dominated) are returned as-is,
/// matching the bulk chain's V2-A13 decision.
KOKKOS_INLINE_FUNCTION real_t penmanMonteithRate(
    real_t airTempC, real_t pressureKpa, real_t airSpecificHumidity,
    real_t shortwaveIncoming, real_t longwaveIn, real_t albedo,
    real_t gFraction, real_t aeroResistance, real_t surfaceResistance) {
  const real_t lambda = latentHeatOfVaporization(airTempC);
  const real_t s = lambda * dqsatdT(airTempC, pressureKpa) / kAirHeatCapacity;
  const real_t rhoA = airDensity(airTempC, pressureKpa);
  const real_t tK = airTempC + 273.15;
  const real_t tK2 = tK * tK;
  const real_t rn = (1.0 - albedo) * shortwaveIncoming +
                    kWaterEmissivity * (longwaveIn - kStefanBoltzmann * tK2 * tK2);
  const real_t available = rn - gFraction * rn;
  const real_t qDeficit =
      saturatedSpecificHumidity(airTempC, pressureKpa) - airSpecificHumidity;
  const real_t numerator = s * available + lambda * rhoA * qDeficit / aeroResistance;
  const real_t denominator = s + 1.0 + surfaceResistance / aeroResistance;
  return numerator / denominator / (lambda * kWaterDensity);
}

}  // namespace frehg::atm

#endif  // FREHG_ATM_PENMANMONTEITH_HPP
