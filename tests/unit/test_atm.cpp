/// \file test_atm.cpp
/// \brief Unit tests for the bulk-aerodynamic module (v2 Q4, plan §3.2 —
///        the §9 "atm/ module unit-covered" requirement). Spot values are
///        computed independently (hand + the scripts/test_g45_gates.py
///        python battery, §6.3 rule 1); the Table-1 chain reproduces the
///        paper's own stated q_a and the V2-A13 closed-form E(0).

#include "atm/BulkAerodynamic.hpp"
#include "atm/PenmanMonteith.hpp"
#include "atm/EtFormulas.hpp"
#include "atm/MetForcing.hpp"
#include "core/Config.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <vector>
#include <fstream>

namespace {

using frehg::real_t;
namespace atm = frehg::atm;

TEST(BulkAerodynamic, TetensSaturationPressure) {
  // Hand values: 0.6108 exp(17.27*20/257.3) = 2.3383 kPa;
  //              0.6108 exp(17.27*25/262.3) = 3.1678 kPa.
  EXPECT_NEAR(atm::tetensSaturationPressure(20.0), 2.3383, 5.0e-4);
  EXPECT_NEAR(atm::tetensSaturationPressure(25.0), 3.1678, 5.0e-4);
}

TEST(BulkAerodynamic, SaturatedSpecificHumidityMatchesTable1) {
  // q_sat(20 C, 101.325 kPa) = 0.622*2.3383/(101.325 - 0.376*2.3383)
  //                          = 0.014480; Table 1 states q_a = 20 % of
  // q_sat as 2.9e-3 — the chain reproduces the paper's own number.
  const real_t qsat = atm::saturatedSpecificHumidity(20.0, 101.325);
  EXPECT_NEAR(qsat, 0.014480, 2.0e-6);
  EXPECT_NEAR(0.2 * qsat, 2.9e-3, 1.0e-5);
}

TEST(BulkAerodynamic, AerodynamicResistanceLiuFit) {
  EXPECT_NEAR(atm::aerodynamicResistance(1.0), 94.909, 1.0e-9);
  EXPECT_NEAR(atm::aerodynamicResistance(2.0), 50.734, 5.0e-3);
}

TEST(BulkAerodynamic, PotentialRateMatchesClosedForm) {
  // The V2-A13 anchor: E(0) = 1.4696e-7 m/s at the Table-1 forcing
  // (rho_a = 1.2041 kg/m^3, R_air = 94.909 s/m, q gradient 0.011584).
  const real_t qa = 0.2 * atm::saturatedSpecificHumidity(20.0, 101.325);
  const real_t e = atm::evaporationRate(20.0, 101.325, 1.0, qa);
  EXPECT_NEAR(e, 1.4696e-7, 2.0e-10);
}

TEST(BulkAerodynamic, SoilRelativeHumidityAndEquilibrium) {
  // Eq. (6): alpha_1(0.0375) = 0.2 (the V2-A13 equilibrium inversion),
  // capped at 1 from w_g = 0.375 (the paper's 0.37 note is this cap).
  EXPECT_NEAR(atm::soilRelativeHumidity(0.0375), 0.2, 1.0e-3);
  EXPECT_DOUBLE_EQ(atm::soilRelativeHumidity(0.41), 1.0);
  EXPECT_DOUBLE_EQ(atm::soilRelativeHumidity(0.375), 1.0);
  // At equilibrium the flux vanishes; below it the flux reverses sign
  // (condensation passes through — the V2-A13 design decision).
  const real_t qa = 0.2 * atm::saturatedSpecificHumidity(20.0, 101.325);
  EXPECT_NEAR(atm::evaporationRate(20.0, 101.325, 1.0, qa,
                                   atm::soilRelativeHumidity(0.0375)),
              0.0, 2.0e-10);
  EXPECT_LT(atm::evaporationRate(20.0, 101.325, 1.0, qa,
                                 atm::soilRelativeHumidity(0.02)),
            0.0);
}

TEST(MetForcing, SamplesConstantsAndFoldsRelativeHumidity) {
  frehg::FrehgConfig cfg;
  cfg.configDir = ".";
  cfg.atmosphere.present = true;
  cfg.atmosphere.airTemperature.constant = 20.0;
  cfg.atmosphere.surfaceTemperature.constant = 20.0;
  cfg.atmosphere.pressure.constant = 101.325;
  cfg.atmosphere.humidityIsRelative = true;
  cfg.atmosphere.relativeHumidity.constant = 0.2;
  cfg.atmosphere.windSpeed.constant = 1.0;
  const atm::MetForcing met(cfg.atmosphere, cfg);
  const atm::MetSample sample = met.sample(0.0);
  EXPECT_DOUBLE_EQ(sample.surfaceTemperatureC, 20.0);
  EXPECT_NEAR(sample.airSpecificHumidity, 2.896e-3, 1.0e-5);
  // The specific-humidity form passes the value through untouched.
  cfg.atmosphere.humidityIsRelative = false;
  cfg.atmosphere.specificHumidity.constant = 2.9e-3;
  const atm::MetForcing met2(cfg.atmosphere, cfg);
  EXPECT_DOUBLE_EQ(met2.sample(0.0).airSpecificHumidity, 2.9e-3);
}

TEST(MetForcing, EvaluatesSeries) {
  const std::string path = ::testing::TempDir() + "/atm_wind_series.dat";
  {
    std::ofstream out(path);
    out << "0.0 1.0\n100.0 3.0\n";
  }
  frehg::FrehgConfig cfg;
  cfg.configDir = ::testing::TempDir();
  cfg.atmosphere.present = true;
  cfg.atmosphere.airTemperature.constant = 20.0;
  cfg.atmosphere.surfaceTemperature.constant = 20.0;
  cfg.atmosphere.pressure.constant = 101.325;
  cfg.atmosphere.specificHumidity.constant = 2.9e-3;
  cfg.atmosphere.windSpeed.fromSeries = true;
  cfg.atmosphere.windSpeed.file = "atm_wind_series.dat";
  const atm::MetForcing met(cfg.atmosphere, cfg);
  EXPECT_NEAR(met.sample(50.0).windSpeed, 2.0, 1.0e-12);
  std::remove(path.c_str());
}

TEST(MetForcing, WindSpeedFloorBoundsStillAir) {
  // Q7 §8.3: the configurable still-air floor (default 0.5 m/s, applied to
  // every bulk consumer at MetForcing::sample) had never been set to a
  // non-default value.
  frehg::FrehgConfig cfg;
  cfg.configDir = ".";
  cfg.atmosphere.present = true;
  cfg.atmosphere.airTemperature.constant = 20.0;
  cfg.atmosphere.surfaceTemperature.constant = 20.0;
  cfg.atmosphere.pressure.constant = 101.325;
  cfg.atmosphere.specificHumidity.constant = 2.9e-3;
  cfg.atmosphere.windSpeed.constant = 0.1;  // below every floor in play
  const atm::MetForcing byDefault(cfg.atmosphere, cfg);
  EXPECT_DOUBLE_EQ(byDefault.sample(0.0).windSpeed, 0.5);
  cfg.atmosphere.windSpeedFloor = 1.25;
  const atm::MetForcing raised(cfg.atmosphere, cfg);
  EXPECT_DOUBLE_EQ(raised.sample(0.0).windSpeed, 1.25);
  cfg.atmosphere.windSpeed.constant = 3.0;  // above the floor: untouched
  const atm::MetForcing untouched(cfg.atmosphere, cfg);
  EXPECT_DOUBLE_EQ(untouched.sample(0.0).windSpeed, 3.0);
}


// ---------------------------------------------------------------------------
// Stage-3 T5 et0 (plan §9): the EtFormulas chain. Reference values from
// scripts/test_et_gates.py (two independent evaluations).
// ---------------------------------------------------------------------------

TEST(AtmEt, BulkPotentialMatchesV20EvaporationRateBitwise) {
  // At r_a = liu2006, r_s = 0, potentialBulk IS the v2.0 expression tree.
  for (const real_t u : {0.5, 1.0, 2.0, 8.0}) {
    for (const real_t alpha : {0.2, 0.7, 1.0}) {
      const real_t qa = 0.3 * atm::saturatedSpecificHumidity(25.0, 101.3);
      const real_t rA = atm::aeroResistance(
          frehg::AtmosphereConfig::AeroResistance::Liu2006, u, {});
      EXPECT_EQ(atm::potentialBulk(30.0, 101.3, rA, 0.0, qa, alpha),
                atm::evaporationRate(30.0, 101.3, u, qa, alpha));
    }
  }
}

TEST(AtmEt, AeroResistanceLaws) {
  // liu2006 forwards to the v2.0 formula.
  EXPECT_EQ(atm::aeroResistance(frehg::AtmosphereConfig::AeroResistance::Liu2006, 1.0, {}),
            atm::aerodynamicResistance(1.0));
  EXPECT_EQ(atm::aeroResistance(frehg::AtmosphereConfig::AeroResistance::Liu2006, 8.0, {}),
            atm::aerodynamicResistance(8.0));
  // fao56 at z_ref = 2 m: u2 = U * 4.87 / ln(130.18) (hand value,
  // test_et_gates.py).
  const real_t u2 = 8.0 * 4.87 / std::log(67.8 * 2.0 - 5.42);
  EXPECT_NEAR(atm::aeroResistance(frehg::AtmosphereConfig::AeroResistance::Fao56, 8.0,
                                  {2.0, 0.01, 0.001, 0.0}),
              208.0 / u2, 1.0e-12);
  // log profile hand value: ln(200) ln(2000) / (0.41^2 * 8).
  const real_t hand = std::log(200.0) * std::log(2000.0) / (0.41 * 0.41 * 8.0);
  EXPECT_NEAR(atm::aeroResistance(frehg::AtmosphereConfig::AeroResistance::LogProfile,
                                  8.0, {2.0, 0.01, 0.001, 0.0}),
              hand, 1.0e-12);
}

TEST(AtmEt, LeePielkeBeta) {
  const real_t fc = 0.2694;  // the silt theta_fc default (below)
  EXPECT_DOUBLE_EQ(atm::limiterBeta(fc, fc), 1.0);
  EXPECT_DOUBLE_EQ(atm::limiterBeta(fc + 0.05, fc), 1.0);
  EXPECT_DOUBLE_EQ(atm::limiterBeta(0.0, fc), 0.0);
  const real_t half = atm::limiterBeta(0.5 * fc, fc);
  EXPECT_GT(half, 0.0);
  EXPECT_LT(half, 1.0);
  // The closed form at theta = fc/2: 0.25 (1 - cos(pi/2))^2 = 0.25.
  EXPECT_NEAR(half, 0.25, 1.0e-12);
}

TEST(AtmEt, SupplyCap) {
  // Evaporation: an unfundable demand takes exactly the drainable store.
  EXPECT_DOUBLE_EQ(atm::supplyCap(1.0e-5, 0.0501, 0.05, 0.41, 0.05, 60.0),
                   (0.0501 - 0.05) * 0.05 / 60.0);
  // A fundable demand passes untouched.
  EXPECT_DOUBLE_EQ(atm::supplyCap(1.0e-9, 0.2, 0.05, 0.41, 0.05, 60.0), 1.0e-9);
  // Condensation: bounded by the 0.9999 theta_s pore room.
  EXPECT_DOUBLE_EQ(atm::supplyCap(-1.0e-5, 0.40999, 0.05, 0.41, 0.05, 60.0),
                   -std::max(0.0, 0.9999 * 0.41 - 0.40999) * 0.05 / 60.0);
}

TEST(AtmEt, ThetaFcDefaultIsRetentionAt330cm) {
  // Silt (et4): theta(h = -3.3) = 0.067 + (1 + (1.6*3.3)^1.37)^(-m) * 0.383
  // with m = 1 - 1/1.37 (test_et_gates.py prints 0.2694-class).
  const real_t fc = atm::thetaFcFromRetention(0.067, 0.45, 1.6, 1.37);
  const real_t m = 1.0 - 1.0 / 1.37;
  const real_t se = std::pow(1.0 + std::pow(1.6 * 3.3, 1.37), -m);
  EXPECT_NEAR(fc, 0.067 + se * (0.45 - 0.067), 1.0e-15);
  EXPECT_GT(fc, 0.067);
  EXPECT_LT(fc, 0.45);
}



TEST(AtmEt, BeerLambertLimits) {
  EXPECT_DOUBLE_EQ(atm::beerLambert(0.5, 0.0), 1.0);   // exactly 1 when bare
  EXPECT_NEAR(atm::beerLambert(1.0, 30.0), std::exp(-30.0), 1.0e-25);
  EXPECT_NEAR(atm::beerLambert(0.5, 2.5), std::exp(-1.25), 1.0e-15);
  EXPECT_LT(atm::beerLambert(0.5, 5.0), atm::beerLambert(0.5, 2.0));
}

TEST(AtmEt, FeddesPiecewise) {
  // Full form: h1 = -0.1, h2 = -0.25, h3 = -4, h4 = -80.
  const real_t h1 = -0.1, h2 = -0.25, h3 = -4.0, h4 = -80.0;
  EXPECT_DOUBLE_EQ(atm::feddes(0.5, h1, h2, h3, h4, true), 0.0);    // above h1
  EXPECT_DOUBLE_EQ(atm::feddes(h1, h1, h2, h3, h4, true), 0.0);
  EXPECT_NEAR(atm::feddes(-0.175, h1, h2, h3, h4, true), 0.5, 1.0e-14);  // mid ramp
  EXPECT_DOUBLE_EQ(atm::feddes(h2, h1, h2, h3, h4, true), 1.0);
  EXPECT_DOUBLE_EQ(atm::feddes(-2.0, h1, h2, h3, h4, true), 1.0);   // plateau
  EXPECT_DOUBLE_EQ(atm::feddes(h3, h1, h2, h3, h4, true), 1.0);
  EXPECT_NEAR(atm::feddes(-42.0, h1, h2, h3, h4, true), 0.5, 1.0e-14);  // dry ramp
  EXPECT_DOUBLE_EQ(atm::feddes(h4, h1, h2, h3, h4, true), 0.0);
  EXPECT_DOUBLE_EQ(atm::feddes(-200.0, h1, h2, h3, h4, true), 0.0);
  // Wetland form (no h1/h2, O7): 1 for every h >= h3, incl. inundation.
  EXPECT_DOUBLE_EQ(atm::feddes(5.0, 0.0, 0.0, h3, h4, false), 1.0);
  EXPECT_DOUBLE_EQ(atm::feddes(h3, 0.0, 0.0, h3, h4, false), 1.0);
  EXPECT_NEAR(atm::feddes(-42.0, 0.0, 0.0, h3, h4, false), 0.5, 1.0e-14);
}

TEST(AtmEt, CanopyResistanceAndRootIntegral) {
  EXPECT_DOUBLE_EQ(atm::canopyResistance(100.0, 2.0), 100.0);  // r_l/(0.5 LAI)
  // Uniform: the interval share; zero at/below L; clipped at L.
  EXPECT_DOUBLE_EQ(atm::rootFractionIntegral(false, 0.4, 0.0, 0.1), 0.25);
  EXPECT_DOUBLE_EQ(atm::rootFractionIntegral(false, 0.4, 0.3, 0.6), 0.25);
  EXPECT_DOUBLE_EQ(atm::rootFractionIntegral(false, 0.4, 0.4, 0.6), 0.0);
  // Linear (triangular): integral of 2(L-d)/L^2; full zone = 1.
  EXPECT_NEAR(atm::rootFractionIntegral(true, 0.4, 0.0, 0.4), 1.0, 1.0e-15);
  EXPECT_NEAR(atm::rootFractionIntegral(true, 0.4, 0.0, 0.2) +
                  atm::rootFractionIntegral(true, 0.4, 0.2, 0.4),
              1.0, 1.0e-15);
  EXPECT_GT(atm::rootFractionIntegral(true, 0.4, 0.0, 0.2),
            atm::rootFractionIntegral(true, 0.4, 0.2, 0.4));  // top-heavy
}


// ---------------------------------------------------------------------------
// Stage-3 T7 (plan §11, et0 part 3): the Penman-Monteith q-form.
// ---------------------------------------------------------------------------

TEST(AtmEt, DqsatDtMatchesCentralDifference) {
  for (real_t t : {5.0, 15.0, 25.0, 38.0}) {
    for (real_t p : {90.0, 101.3}) {
      const real_t h = 1.0e-4;
      const real_t fd = (atm::saturatedSpecificHumidity(t + h, p) -
                         atm::saturatedSpecificHumidity(t - h, p)) /
                        (2.0 * h);
      EXPECT_NEAR(atm::dqsatdT(t, p), fd, 1.0e-8 * fd);
    }
  }
}

TEST(AtmEt, PmEquilibriumLimit) {
  // q_sat(T_air) = q_a: the aerodynamic term vanishes exactly and
  // lambda rho_w E = s A / (s + 1) (r_s = 0).
  const real_t t = 25.0, p = 101.3, rs = 500.0, lw = 350.0;
  const real_t albedo = 0.1, g = 0.1, ra = 60.0;
  const real_t qa = atm::saturatedSpecificHumidity(t, p);
  const real_t e = atm::penmanMonteithRate(t, p, qa, rs, lw, albedo, g, ra, 0.0);
  const real_t lambda = atm::latentHeatOfVaporization(t);
  const real_t s = lambda * atm::dqsatdT(t, p) / atm::kAirHeatCapacity;
  const real_t tk = t + 273.15;
  const real_t tk2 = tk * tk;
  const real_t rn = (1.0 - albedo) * rs +
                    atm::kWaterEmissivity * (lw - atm::kStefanBoltzmann * tk2 * tk2);
  const real_t a = rn - g * rn;
  EXPECT_DOUBLE_EQ(e, s * a / (s + 1.0) / (lambda * atm::kWaterDensity));
  // Negative (spec): a perturbed s breaks the identity.
  const real_t sBad = s * (1.0 + 1.0e-6);
  EXPECT_NE(e, sBad * a / (sBad + 1.0) / (lambda * atm::kWaterDensity));
}

TEST(AtmEt, PmAerodynamicLimit) {
  // R_s = 0 and LW_in = sigma T_K^4 make R_n (and so A and G) exactly
  // zero: only the aerodynamic term remains,
  //   lambda rho_w E = lambda rho_a (q_sat - q_a)/r_a / (s + 1 + r_s/r_a),
  // and the G fraction is irrelevant.
  const real_t t = 25.0, p = 101.3, ra = 60.0, rsRes = 70.0;
  const real_t tk = t + 273.15;
  const real_t tk2 = tk * tk;
  const real_t lw = atm::kStefanBoltzmann * tk2 * tk2;
  const real_t qsat = atm::saturatedSpecificHumidity(t, p);
  const real_t qa = 0.5 * qsat;
  const real_t e = atm::penmanMonteithRate(t, p, qa, 0.0, lw, 0.1, 0.1, ra, rsRes);
  const real_t lambda = atm::latentHeatOfVaporization(t);
  const real_t s = lambda * atm::dqsatdT(t, p) / atm::kAirHeatCapacity;
  const real_t expected = lambda * atm::airDensity(t, p) * (qsat - qa) / ra /
                          (s + 1.0 + rsRes / ra) / (lambda * atm::kWaterDensity);
  EXPECT_DOUBLE_EQ(e, expected);
  EXPECT_DOUBLE_EQ(e, atm::penmanMonteithRate(t, p, qa, 0.0, lw, 0.1, 0.9, ra, rsRes));
}

TEST(AtmEt, PmHighResistanceLimit) {
  // The suppression ratio is (s + 1)/(s + 1 + r_s/r_a): the plan's
  // 1e-12 bound at r_s = 1e12 needs r_a <= ~0.25 at 25 C (s ~ 2.9), so
  // this is a pure math-limit check at r_a = 0.1 (S3-A10).
  const real_t t = 25.0, p = 101.3, rs = 500.0, lw = 350.0, ra = 0.1;
  const real_t qa = 0.5 * atm::saturatedSpecificHumidity(t, p);
  const real_t e0 = atm::penmanMonteithRate(t, p, qa, rs, lw, 0.1, 0.1, ra, 0.0);
  const real_t eHi = atm::penmanMonteithRate(t, p, qa, rs, lw, 0.1, 0.1, ra, 1.0e12);
  EXPECT_GT(e0, 0.0);
  EXPECT_GT(eHi, 0.0);
  EXPECT_LT(eHi, 1.0e-12 * e0);
}

TEST(AtmEt, PmCommensurateWithBulkAtEnergyBalanceTemperature) {
  // For each available energy A, bisect the surface temperature T_s*
  // that closes A = lambda rho_a (q_sat(T_s) - q_a)/r_a
  //               + rho_a c_pa (T_s - T_air)/r_a   (rho_a at T_air),
  // and compare the bulk rate at T_s* (same rho_a convention) with the
  // PM rate at T_air. PM is that balance with q_sat linearized at T_air:
  // the ABSOLUTE gap E_PM - E_bulk is second order in (T_s* - T_air)
  // (both rates vanish linearly; the secant-vs-tangent error enters at
  // first order RELATIVE — measured 0.98, recorded in S3-A10, not
  // gated). The fitted log-log slope of the absolute gap over A in
  // {50, 100, 200, 400} W/m^2 must be >= 1.8. (The SHIPPED bulk
  // evaluates rho_a at T_s — a convention difference recorded here,
  // not gated.)
  // q_a = q_sat(T_air) zeroes the aerodynamic baseline flux, so T_s*
  // rises monotonically with A and the gap isolates the linearization
  // error (with q_a < q_sat the baseline can exceed A and put T_s*
  // below T_air — the balance still closes, but the fit's abscissa
  // would cross zero; S3-A10).
  const real_t tA = 25.0, p = 101.3, ra = 60.0;
  const real_t qa = atm::saturatedSpecificHumidity(tA, p);
  const real_t rhoA = atm::airDensity(tA, p);
  const real_t lambda = atm::latentHeatOfVaporization(tA);
  const real_t tk = tA + 273.15;
  const real_t tk2 = tk * tk;
  std::vector<real_t> logDt, logGap;
  for (const real_t aW : {50.0, 100.0, 200.0, 400.0}) {
    // Radiation producing available = aW exactly (r_s = 0, albedo 0,
    // G fraction 0): LW_in = sigma T_K^4 + aW / eps.
    const real_t lw =
        atm::kStefanBoltzmann * tk2 * tk2 + aW / atm::kWaterEmissivity;
    const real_t ePm = atm::penmanMonteithRate(tA, p, qa, 0.0, lw, 0.0, 0.0, ra, 0.0);
    const auto residual = [&](real_t ts) {
      return lambda * rhoA * (atm::saturatedSpecificHumidity(ts, p) - qa) / ra +
             rhoA * atm::kAirHeatCapacity * (ts - tA) / ra - aW;
    };
    real_t lo = tA, hi = tA + 30.0;
    ASSERT_LT(residual(lo), 0.0);
    ASSERT_GT(residual(hi), 0.0);
    for (int it = 0; it < 200; ++it) {
      const real_t mid = 0.5 * (lo + hi);
      (residual(mid) > 0.0 ? hi : lo) = mid;
    }
    const real_t ts = 0.5 * (lo + hi);
    const real_t eBulk =
        rhoA * (atm::saturatedSpecificHumidity(ts, p) - qa) / (ra * atm::kWaterDensity);
    ASSERT_GT(ts, tA);
    logDt.push_back(std::log(ts - tA));
    // Absolute gap in m/s (any fixed normalization shifts the log by a
    // constant and leaves the fitted slope unchanged).
    logGap.push_back(std::log(std::abs(ePm - eBulk)));
  }
  // Least-squares slope of log(gap) against log(dT).
  const auto n = static_cast<real_t>(logDt.size());
  real_t sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
  for (std::size_t k = 0; k < logDt.size(); ++k) {
    sx += logDt[k];
    sy += logGap[k];
    sxx += logDt[k] * logDt[k];
    sxy += logDt[k] * logGap[k];
  }
  const real_t slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  EXPECT_GE(slope, 1.8);
}

}  // namespace
