/// \file AdaptiveStep.cpp
/// \brief The adaptive subsurface time-step controller (legacy
///        adaptive_time_step, groundwater.c:1651-1727).
///
/// The controller grows dtg by 1.25 when the largest per-cell flux change
/// dq = |q_in - q_out| dtg / dz stays below dq_grow, shrinks it by 0.75 when
/// dq exceeds dq_shrink, caps it with the Courant limit
/// Co_max dz / (dK/dθ) evaluated on unsaturated cells, clamps it to
/// [dt_min, dt_max], and min-reduces across ranks (:1668-1726). Column top
/// cells are excluded from both criteria, exactly as legacy's istop guard
/// does (:1672). The per-cell flux sums divide by the cell's own face areas
/// (not the per-face areas), preserving the legacy expression verbatim
/// (:1674-1675).
///
/// The truncation-error estimator (err/err_max, :1670-1671) and the
/// saturated/unsaturated-zone scan (:1662-1665) compute values legacy never
/// reads; both are dropped as dead code (docs/theory/removed-features.md).
/// The grow/shrink factors are the plan §3.1 item 5 named constants.

#include "gw/RichardsSolver.hpp"
#include "gw/VanGenuchten.hpp"

#include <mpi.h>

namespace frehg::gw {

namespace {

/// Adaptive-step factors (legacy r_red/r_inc, groundwater.c:1657-1658;
/// named constants per plan §3.1 item 5).
constexpr real_t kShrinkFactor = 0.75;
constexpr real_t kGrowFactor = 1.25;
/// Courant-limit search start (legacy dt_Comin = 1e8, :1660).
constexpr real_t kCourantStart = 1.0e8;

template <class V3>
KOKKOS_INLINE_FUNCTION VgSoil soilAt(const V3& vga, const V3& vgn, const V3& wcs, const V3& wcr,
                                     const V3& aev, int j, int i, int k) {
  return VgSoil{vga(j, i, k), vgn(j, i, k), wcs(j, i, k), wcr(j, i, k), aev(j, i, k)};
}

}  // namespace

void RichardsSolver::adaptTimeStep(real_t dtg) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t az = mesh_.areaZ();
  const real_t courantMax = courantMax_;

  Field3<real_t> wc = wc_, ksz = ksz_;
  Field3<real_t> qx = qx_, qy = qy_, qzF = qzF_;
  Field3<real_t> vga = vga_, vgn = vgn_, wcs = wcs_, wcr = wcr_, aev = aev_;
  Field3<real_t> ax = mesh_.areaX(), ay = mesh_.areaY();
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<PetscInt> gid = grid_.gid3();

  real_t dqMax = 0.0;
  Kokkos::parallel_reduce(
      "gw_adaptive_dq",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& value) {
        const bool isTop = (k == 0) || (gid(j, i, k - 1) < 0);
        if (gid(j, i, k) < 0 || isTop) {
          return;
        }
        const real_t qin = qx(j, i, k) / ax(j, i, k) + qy(j, i, k) / ay(j, i, k) +
                           qzF(j, i, k + 1) / az;
        const real_t qou = qx(j, i - 1, k) / ax(j, i, k) + qy(j - 1, i, k) / ay(j, i, k) +
                           qzF(j, i, k) / az;
        const real_t dq = Kokkos::fabs(qin - qou) * dtg / dz3d(j, i, k);
        if (dq > value) {
          value = dq;
        }
      },
      Kokkos::Max<real_t>(dqMax));

  real_t dtCoMin = kCourantStart;
  Kokkos::parallel_reduce(
      "gw_adaptive_courant",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& value) {
        const bool isTop = (k == 0) || (gid(j, i, k - 1) < 0);
        if (gid(j, i, k) < 0 || isTop) {
          return;
        }
        if (wc(j, i, k) < wcs(j, i, k)) {
          const real_t dKdwc = conductivityDerivative(
              soilAt(vga, vgn, wcs, wcr, aev, j, i, k), wc(j, i, k), ksz(j, i, k));
          const real_t dtCo = courantMax * dz3d(j, i, k) / dKdwc;
          if (dtCo < value) {
            value = dtCo;
          }
        }
      },
      Kokkos::Min<real_t>(dtCoMin));

  real_t next = dtg;
  if (dqMax > dqShrink_) {
    next = dtg * kShrinkFactor;
  } else if (dqMax < dqGrow_) {
    next = dtg * kGrowFactor;
  }
  if (next > dtCoMin) {
    next = dtCoMin;
  }
  if (next > dtMax_) {
    next = dtMax_;
  }
  if (next < dtMin_) {
    next = dtMin_;
  }
  real_t global = next;
  MPI_Allreduce(&next, &global, 1, MPI_DOUBLE, MPI_MIN, grid_.comm());
  dtgNext_ = global;
}

}  // namespace frehg::gw
