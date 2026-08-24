/// \file test_vangenuchten.cpp
/// \brief Van Genuchten closure tests (plan §8.1 test_vangenuchten): θ(h),
///        C(h), K(h) against hand-computed values for b2's Warrick soil and
///        b6's sand at the prescribed heads, explicitly covering the
///        aev = -0.02 cutoff and the θr/θs clamps, at 1e-12 relative.
///
/// Reference values were hand-evaluated from the exact legacy formulas
/// (subroutines.c:280-483, use_vg = 1 / use_mvg = 0 branch) in IEEE double
/// precision with the code's operation order; a 40-digit cross-check agrees
/// to each expression's double-precision conditioning (the Mualem bracket
/// amplifies rounding to ~2e-12 for the driest conductivities, so exact
/// arithmetic is not the right reference for a double-precision closure).

#include "gw/VanGenuchten.hpp"

#include <gtest/gtest.h>

namespace {

using frehg::real_t;
using frehg::gw::VgSoil;

constexpr real_t kRelTol = 1.0e-12;

void expectNear(real_t actual, real_t expected) {
  if (expected == 0.0) {
    EXPECT_EQ(actual, 0.0);
    return;
  }
  EXPECT_NEAR(actual / expected, 1.0, kRelTol);
}

/// b2's Warrick (1971) soil (benchmarks/b2-gw): alpha 1.43, n 1.56,
/// theta_s 0.33, theta_r 0, aev -0.02, Ksz 2.89e-6.
const VgSoil kWarrick{1.43, 1.56, 0.33, 0.0, -0.02};
constexpr real_t kWarrickKs = 2.89e-6;

/// b6's Kuan sand (benchmarks/b6-kuan): alpha 5.9, n 2.68, theta_s 0.46,
/// theta_r 0.04, aev -0.02, Ksz 5.26e-3.
const VgSoil kSand{5.9, 2.68, 0.46, 0.04, -0.02};
constexpr real_t kSandKs = 5.26e-3;

TEST(VanGenuchten, WaterContentMatchesHandValuesWarrick) {
  expectNear(frehg::gw::waterContentFromHead(kWarrick, -10.0), 0.0739753893647920896);
  expectNear(frehg::gw::waterContentFromHead(kWarrick, -1.0), 0.229597833180606);
  expectNear(frehg::gw::waterContentFromHead(kWarrick, -0.02), 0.329538302833622731);
  // The aev cutoff: h = -0.01 sits above aev = -0.02, so θ jumps to θs even
  // though h < 0 (the live-aev behavior the plan §3.1 mandates).
  EXPECT_EQ(frehg::gw::waterContentFromHead(kWarrick, -0.01), 0.33);
  EXPECT_EQ(frehg::gw::waterContentFromHead(kWarrick, 0.0), 0.33);
  EXPECT_EQ(frehg::gw::waterContentFromHead(kWarrick, 0.5), 0.33);
}

TEST(VanGenuchten, WaterContentMatchesHandValuesSand) {
  expectNear(frehg::gw::waterContentFromHead(kSand, -10.0), 0.0404448494147192347);
  expectNear(frehg::gw::waterContentFromHead(kSand, -1.0), 0.0611781812182591983);
  expectNear(frehg::gw::waterContentFromHead(kSand, -0.02), 0.459145092411154812);
  EXPECT_EQ(frehg::gw::waterContentFromHead(kSand, -0.01), 0.46);
  EXPECT_EQ(frehg::gw::waterContentFromHead(kSand, 0.0), 0.46);
  EXPECT_EQ(frehg::gw::waterContentFromHead(kSand, 0.5), 0.46);
}

TEST(VanGenuchten, CapacityMatchesHandValues) {
  expectNear(frehg::gw::capacityFromHead(kWarrick, -10.0), 0.00407832963118645968);
  expectNear(frehg::gw::capacityFromHead(kWarrick, -1.0), 0.0817714043754333137);
  expectNear(frehg::gw::capacityFromHead(kWarrick, -0.02), 0.0359170881291010571);
  // C(h) cuts off at h > 0, not at aev (legacy compute_ch:371-374): the
  // capacity stays live between aev and 0.
  expectNear(frehg::gw::capacityFromHead(kWarrick, -0.01), 0.0244480654553699461);
  EXPECT_EQ(frehg::gw::capacityFromHead(kWarrick, 0.0), 0.0);
  EXPECT_EQ(frehg::gw::capacityFromHead(kWarrick, 0.5), 0.0);
  expectNear(frehg::gw::capacityFromHead(kSand, -10.0), 7.47333600479079924e-5);
  expectNear(frehg::gw::capacityFromHead(kSand, -1.0), 0.0352762354400144718);
  expectNear(frehg::gw::capacityFromHead(kSand, -0.02), 0.114255158543752424);
  expectNear(frehg::gw::capacityFromHead(kSand, -0.01), 0.0358164977166531603);
  EXPECT_EQ(frehg::gw::capacityFromHead(kSand, 0.0), 0.0);
}

TEST(VanGenuchten, ConductivityMatchesHandValues) {
  expectNear(frehg::gw::conductivityFromHead(kWarrick, -10.0, kWarrickKs),
             4.2896686667440645e-11);
  expectNear(frehg::gw::conductivityFromHead(kWarrick, -1.0, kWarrickKs),
             5.420575919043507e-8);
  expectNear(frehg::gw::conductivityFromHead(kWarrick, -0.02, kWarrickKs),
             2.153644648328222e-6);
  // K(h) also cuts off at h > 0, not aev: between aev and 0 it stays on the
  // Mualem curve (legacy compute_K:401-403).
  expectNear(frehg::gw::conductivityFromHead(kWarrick, -0.01, kWarrickKs),
             2.378796558153712e-6);
  EXPECT_EQ(frehg::gw::conductivityFromHead(kWarrick, 0.0, kWarrickKs), kWarrickKs);
  EXPECT_EQ(frehg::gw::conductivityFromHead(kWarrick, 0.5, kWarrickKs), kWarrickKs);
  expectNear(frehg::gw::conductivityFromHead(kSand, -10.0, kSandKs), 2.1678893902335157e-14);
  expectNear(frehg::gw::conductivityFromHead(kSand, -1.0, kSandKs), 3.379411314506458e-8);
  expectNear(frehg::gw::conductivityFromHead(kSand, -0.02, kSandKs), 0.004969259950282671);
  expectNear(frehg::gw::conductivityFromHead(kSand, -0.01, kSandKs), 0.005169012561134972);
  EXPECT_EQ(frehg::gw::conductivityFromHead(kSand, 0.5, kSandKs), kSandKs);
}

TEST(VanGenuchten, HeadInversionRoundTripsAndClamps) {
  // b2's initial state: θ = 0.033 inverts to the value its golden embeds.
  expectNear(frehg::gw::headFromWaterContent(kWarrick, 0.033), -42.6502807991778219);
  expectNear(frehg::gw::headFromWaterContent(kSand, 0.25), -0.220394548940179914);
  // At or above θs the inversion returns 0 (saturation).
  EXPECT_EQ(frehg::gw::headFromWaterContent(kWarrick, 0.33), 0.0);
  EXPECT_EQ(frehg::gw::headFromWaterContent(kWarrick, 0.4), 0.0);
  // θ at or below θr floors at θr + 1e-7 (the legacy eps guard) — finite,
  // very dry.
  const real_t floored = frehg::gw::headFromWaterContent(kWarrick, 0.0);
  EXPECT_EQ(floored, frehg::gw::headFromWaterContent(kWarrick, -1.0));
  EXPECT_LT(floored, -1.0e3);
  EXPECT_TRUE(std::isfinite(floored));
  // Round trip through θ(h) within the unsaturated range.
  const real_t h = -2.5;
  const real_t wc = frehg::gw::waterContentFromHead(kSand, h);
  EXPECT_NEAR(frehg::gw::headFromWaterContent(kSand, wc) / h, 1.0, 1.0e-10);
}

TEST(VanGenuchten, ConductivityDerivativeMatchesHandValuesAndLimiter) {
  expectNear(frehg::gw::conductivityDerivative(kWarrick, 0.2, kWarrickKs),
             7.03506545539082082e-7);
  // The 0.9999 θs limiter keeps the derivative finite as θ approaches θs:
  // both evaluations collapse onto the limiter value.
  const real_t nearSat = frehg::gw::conductivityDerivative(kSand, 0.459999, kSandKs);
  expectNear(nearSat, 0.635589131661863617);
  EXPECT_EQ(frehg::gw::conductivityDerivative(kSand, 0.4599999, kSandKs), nearSat);
}

}  // namespace
