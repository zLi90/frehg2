/// \file test_atm.cpp
/// \brief Unit tests for the bulk-aerodynamic module (v2 Q4, plan §3.2 —
///        the §9 "atm/ module unit-covered" requirement). Spot values are
///        computed independently (hand + the scripts/test_g45_gates.py
///        python battery, §6.3 rule 1); the Table-1 chain reproduces the
///        paper's own stated q_a and the V2-A13 closed-form E(0).

#include "atm/BulkAerodynamic.hpp"
#include "atm/MetForcing.hpp"
#include "core/Config.hpp"

#include <gtest/gtest.h>

#include <cstdio>
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

}  // namespace
