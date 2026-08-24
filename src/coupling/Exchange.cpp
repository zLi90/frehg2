/// \file Exchange.cpp
/// \brief The exchange kernels: window staging, the seepage application
///        (legacy subsurface_source, shallowwater.c:639-683), and the
///        coupled audit reductions.

#include "coupling/Coupler.hpp"

namespace frehg::coupling {

void Coupler::beginWindow() {
  // Seed the window's infiltration budget from the post-solve depth
  // (amendment A9: the limit debits this budget across subcycles; every
  // coupled deposit credits it, mirroring the live legacy dept the limit
  // read) and clear the per-step deposit tally.
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  Field2<real_t> avail = avail_, gain = gain_, dept = surface_.depth();
  Kokkos::parallel_for(
      "cpl_begin_window",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        avail(j, i) = dept(j, i);
        gain(j, i) = 0.0;
      });
}

void Coupler::applySeepage(real_t dt) {
  // subsurface_source (shallowwater.c:639-683) on the volume accumulator:
  // infiltration applies in full with the eta >= bottom clamp (:652-658, the
  // clamp remainder is audited); seepage onto a wet surface applies in full
  // (:677-680); seepage onto a dry surface is withheld until the
  // accumulated depth exceeds min_depth (:668-675 — the legacy reset_seepage
  // hold). qss keeps the applied-rate observable the legacy monitors wrote.
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t minDepth = minDepth_;
  Field2<real_t> seepAccum = seepAccum_, gain = gain_, qss = qss_;
  Field2<real_t> eta = surface_.eta(), dept = surface_.depth(), bottom = surface_.bottom();

  real_t applied = 0.0;
  real_t discarded = 0.0;
  Kokkos::parallel_reduce(
      "cpl_apply_seepage",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumA, real_t& sumD) {
        const real_t a = seepAccum(j, i);
        qss(j, i) = a / dt;
        if (a < 0.0) {
          const real_t before = eta(j, i);
          eta(j, i) += a;
          if (eta(j, i) < bottom(j, i)) {
            eta(j, i) = bottom(j, i);
          }
          const real_t change = eta(j, i) - before;
          gain(j, i) += change;
          sumA += change;
          sumD += change - a;  // >= 0: remainder cut by the clamp
          seepAccum(j, i) = 0.0;
        } else if (dept(j, i) > minDepth || a > minDepth) {
          eta(j, i) += a;
          gain(j, i) += a;
          sumA += a;
          seepAccum(j, i) = 0.0;
        }
        // else: hold (legacy reset_seepage = 0) — the accumulator carries
        // the volume into the following steps.
      },
      applied, discarded);
  const real_t cellArea = grid_.dx() * grid_.dy();
  audit_.applied = applied * cellArea;
  audit_.discarded = discarded * cellArea;
}

void Coupler::reduceAudit() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  Field2<real_t> gain = gain_, seepAccum = seepAccum_;
  real_t gainSum = 0.0;
  real_t heldSum = 0.0;
  Kokkos::parallel_reduce(
      "cpl_reduce_audit",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumG, real_t& sumH) {
        sumG += gain(j, i);
        sumH += seepAccum(j, i);
      },
      gainSum, heldSum);
  const real_t cellArea = grid_.dx() * grid_.dy();
  audit_.surfaceGain = gainSum * cellArea;
  audit_.accumVolume = heldSum * cellArea;
}

}  // namespace frehg::coupling
