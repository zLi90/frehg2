/// \file FaceConductivity.hpp
/// \brief The Stage-3 T11 shared face-conductivity rule (plan §12C P1/P2/P6,
///        owner decision O17; amendment S3-A16).
///
/// Under groundwater.predictor_corrector: consistent the predictor's matrix
/// legs and the corrector's Darcy fluxes evaluate every interior face
/// through ONE of these functions, so the two halves of the PCA step cannot
/// disagree on a face again (the legacy one-sided corrector Ks — x/y: the
/// minus cell's, z: the lower cell's — created water at every saturated Ks
/// contrast through the reallocation restore θ := θ(h): 72 % of the net
/// exfiltration on the v30 H2 beach, S3-A16).
///
/// Two means (groundwater.face_conductivity):
///  - arithmetic (default): 0.5 (K_a(h_a) + K_b(h_b)), each flank's own Ks
///    and vG closure — the predictor's v1 rule, now in both halves;
///  - harmonic: Ks_H · ½ (kr_a + kr_b), with Ks_H the thickness-weighted
///    harmonic mean of the flanks' SATURATED conductivities and kr each
///    flank's own relative permeability K(h)/Ks. Exact for saturated series
///    flow; a Ks = 0 flank seals the face exactly; equal-Ks flanks reduce
///    it to the arithmetic rule in exact arithmetic (Ks_H = Ks and
///    Ks · ½(kr_a + kr_b) = ½(K_a + K_b)) — to rounding, not bitwise, so
///    the K3 bitwise statement is the ARITHMETIC mode's alone.
///
/// The legacy mode never calls these: the pinned corrector kernels keep
/// their literal v1 expression trees behind `if constexpr` (the S3-A8.1
/// fma rule), and the predictor's legacy path is its existing code.

#ifndef FREHG_GW_FACECONDUCTIVITY_HPP
#define FREHG_GW_FACECONDUCTIVITY_HPP

#include "core/Types.hpp"
#include "gw/VanGenuchten.hpp"

namespace frehg::gw {

/// Arithmetic mean of the flanks' own-K closures (P2).
KOKKOS_INLINE_FUNCTION real_t faceKArithmetic(const VgSoil& soilA, real_t hA,
                                              real_t ksA, const VgSoil& soilB,
                                              real_t hB, real_t ksB) {
  const real_t ka = conductivityFromHead(soilA, hA, ksA);
  const real_t kb = conductivityFromHead(soilB, hB, ksB);
  return 0.5 * (ka + kb);
}

/// Harmonic Ks × arithmetic kr (P6), lengths = the flank thicknesses
/// normal to the face (dz for z faces; dx = dx / dy = dy reduce the
/// weighting to the plain harmonic mean).
KOKKOS_INLINE_FUNCTION real_t faceKHarmonic(const VgSoil& soilA, real_t hA,
                                            real_t ksA, real_t lenA,
                                            const VgSoil& soilB, real_t hB,
                                            real_t ksB, real_t lenB) {
  if (ksA <= 0.0 || ksB <= 0.0) {
    return 0.0;  // a zero-Ks flank seals the face exactly (pc4)
  }
  const real_t ksH = (lenA + lenB) / (lenA / ksA + lenB / ksB);
  const real_t krA = conductivityFromHead(soilA, hA, ksA) / ksA;
  const real_t krB = conductivityFromHead(soilB, hB, ksB) / ksB;
  return ksH * 0.5 * (krA + krB);
}

/// The consistent face rule, compile-time selected (the fma rule: no
/// runtime branch inside pinned kernels).
template <bool Harmonic>
KOKKOS_INLINE_FUNCTION real_t faceK(const VgSoil& soilA, real_t hA, real_t ksA,
                                    real_t lenA, const VgSoil& soilB, real_t hB,
                                    real_t ksB, real_t lenB) {
  if constexpr (Harmonic) {
    return faceKHarmonic(soilA, hA, ksA, lenA, soilB, hB, ksB, lenB);
  } else {
    (void)lenA;
    (void)lenB;
    return faceKArithmetic(soilA, hA, ksA, soilB, hB, ksB);
  }
}

}  // namespace frehg::gw

#endif  // FREHG_GW_FACECONDUCTIVITY_HPP
