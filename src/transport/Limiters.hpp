/// \file Limiters.hpp
/// \brief Advection face-value schemes for the scalar transport (plan §3.1).
///
/// Provenance: tvd_superbee (subroutines.c:541-561). The scheme returns the
/// second-order TVD face value s_c + φ(r)/2 (1 − |u| dt/δ)(s_p − s_c) with
/// the superbee limiter φ(r) = max(0, min(2r, 1), min(r, 2)) evaluated from
/// the upwind ratio r = (s_c − s_m)/(s_p − s_c); the first-order upwind
/// alternative is the plain donor value chosen by the caller's flux sign.

#ifndef FREHG_TRANSPORT_LIMITERS_HPP
#define FREHG_TRANSPORT_LIMITERS_HPP

#include "core/Types.hpp"

namespace frehg::transport {

/// Superbee limiter function φ(r) exactly as legacy computes it
/// (subroutines.c:546-559): r1 = min(2r, 1), r2 = min(r, 2),
/// φ = max(r1, r2) when either is positive, else 0.
KOKKOS_INLINE_FUNCTION real_t superbeePhi(real_t r) {
  real_t r1 = 1.0;
  real_t r2 = 2.0;
  if (2.0 * r < 1.0) {
    r1 = 2.0 * r;
  }
  if (r < 2.0) {
    r2 = r;
  }
  if (r1 > 0.0 || r2 > 0.0) {
    return (r1 > r2) ? r1 : r2;
  }
  return 0.0;
}

/// TVD superbee face value (tvd_superbee, subroutines.c:541-561).
/// \param sp scalar of the acceptor cell (downwind of the face)
/// \param sc scalar of the donor cell (upwind of the face)
/// \param sm scalar of the donor's upwind neighbor
/// \param u face velocity used for the Courant factor [m/s]
/// \param delta cell size along the face normal [m]
/// \param dt transport step [s]
KOKKOS_INLINE_FUNCTION real_t tvdSuperbee(real_t sp, real_t sc, real_t sm, real_t u,
                                          real_t delta, real_t dt) {
  const real_t coef = Kokkos::fabs(u) * dt / delta;
  real_t phi = 0.0;
  if (sp != sc) {
    phi = superbeePhi((sc - sm) / (sp - sc));
  }
  return sc + 0.5 * phi * (1.0 - coef) * (sp - sc);
}

}  // namespace frehg::transport

#endif  // FREHG_TRANSPORT_LIMITERS_HPP
