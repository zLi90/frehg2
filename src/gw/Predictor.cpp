/// \file Predictor.cpp
/// \brief The PCA predictor: face conductivities, the 7-point linear head
///        system, and boundary-ghost head enforcement.
///
/// Provenance: compute_K_face (groundwater.c:202-335), groundwater_mat_coeff
/// (:503-569), groundwater_rhs (:572-682), build/solve_groundwater_system
/// (:686-748, LASPack CG+SSOR at rtol 1e-8 replaced by the plan §5.1 PETSc
/// CG at the same rtol), enforce_head_bc (:751-813).
///
/// Legacy defects fixed on port (§2.1 hazard class; the full table lives in
/// docs/theory/groundwater.md):
///  - interface face conductivities are evaluated two-sided from exchanged
///    neighbor state on every rank (legacy computed one-sided values that
///    disagreed across the interface);
///  - ghost soil parameters are mirrored instead of read uninitialized;
///  - the prescribed-flux bottom RHS term carries dtg/dz (legacy omitted
///    both factors on qbot, groundwater.c:631);
///  - the prescribed-head bottom RHS term enters with the Dirichlet sign
///    (legacy added Gzp*h instead of subtracting, groundwater.c:633);
///  - free drainage acts on the bottom condition (legacy tested
///    bctype_GW[5], the top code, in the bottom branch — darcy_flux
///    subroutines.c:190).
/// None of these paths is exercised by any benchmark configuration.

#include "core/Logger.hpp"
#include "core/Timer.hpp"
#include "gw/RichardsSolver.hpp"
#include "gw/VanGenuchten.hpp"

namespace frehg::gw {

namespace {

/// Soil parameters of one cell, staged for the closures.
template <class V3>
KOKKOS_INLINE_FUNCTION VgSoil soilAt(const V3& vga, const V3& vgn, const V3& wcs, const V3& wcr,
                                     const V3& aev, int j, int i, int k) {
  return VgSoil{vga(j, i, k), vgn(j, i, k), wcs(j, i, k), wcr(j, i, k), aev(j, i, k)};
}

/// Ghost-head rule for one domain-edge cell (enforce_head_bc:757-788).
/// A hydrostatic edge is classified by its configured stage `value` (the
/// water-table elevation) relative to the edge bed. A stage ABOVE the bed is a
/// ponded/tidal boundary: in coupled runs the ghost follows the live local
/// surface exactly as legacy does (:773-788) — the density-scaled hydrostatic
/// column of the edge cell's own bed plus its current surface depth — so a
/// seaward stage tracks the tide. A stage at or below the bed is a subsurface
/// water table (a lateral hillslope edge, or any uncoupled run): the ghost is
/// hydrostatic about the configured stage, (value - zc)*rFace (the P0 form,
/// r-scaled once density coupling is active). Without this split a coupled dry
/// edge would be pinned to the ground surface and flood.
/// Legacy scales the y+ depth term by the face density ratio but leaves the
/// y- depth term unscaled (:777 vs :787) — preserved, with the x edges
/// mirroring the y rules (legacy had no x head conditions; P2
/// generalization, docs/theory/groundwater.md).
/// \param rFace boundary-face density ratio (1 until P4 activates it)
/// \param plusSide true for the x+/y+ edges (the r-scaled depth term)
KOKKOS_INLINE_FUNCTION real_t ghostHead(int code, real_t value, real_t interior, real_t zc,
                                        bool coupled, real_t bed, real_t depth, real_t rFace,
                                        bool plusSide) {
  if (code == static_cast<int>(GwBcCode::Head)) {
    return value;
  }
  if (code == static_cast<int>(GwBcCode::HeadHydrostatic)) {
    // A configured stage ABOVE the edge bed is a ponded/tidal boundary: in a
    // coupled run its ghost follows the live local surface (bed + depth) so a
    // seaward stage tracks the tide. A stage at or below the bed is a
    // subsurface water table (a lateral hillslope edge, or any uncoupled run):
    // the ghost is hydrostatic about the configured stage, the same fixed form
    // the uncoupled path uses. Without this split a coupled dry edge would be
    // pinned to the ground surface and flood, ignoring the configured stage.
    if (coupled && value > bed) {
      const real_t hydro = (bed - zc) * rFace;
      return plusSide ? hydro + depth * rFace : hydro + depth;
    }
    return (value - zc) * rFace;
  }
  return interior;
}

}  // namespace

void RichardsSolver::computeFaceConductivity() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const bool full3d = useFull3d_;

  Field3<real_t> h = h_, wc = wc_, kx = kx_, ky = ky_, kzF = kzF_;
  Field3<real_t> ksx = ksx_, ksy = ksy_, ksz = ksz_;
  Field3<real_t> vga = vga_, vgn = vgn_, wcs = wcs_, wcr = wcr_, aev = aev_;
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> topCode = topCode_, botCode = botCode_;
  Field2<int> sideXm = sideCodeXm_, sideXp = sideCodeXp_;
  Field2<int> sideYm = sideCodeYm_, sideYp = sideCodeYp_;

  // Arithmetic-mean face conductivities from each flank's own state and
  // saturated conductivity (compute_K_face:210-235). Interface faces use
  // the exchanged neighbor state, so every rank sees the same value.
  Kokkos::parallel_for(
      "gw_face_k_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 0, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const real_t km = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                               h(j, i, k), ksx(j, i, k));
        const real_t kp = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i + 1, k),
                                               h(j, i + 1, k), ksx(j, i + 1, k));
        kx(j, i, k) = 0.5 * (km + kp);
      });
  Kokkos::parallel_for(
      "gw_face_k_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {0, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const real_t km = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                               h(j, i, k), ksy(j, i, k));
        const real_t kp = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j + 1, i, k),
                                               h(j + 1, i, k), ksy(j + 1, i, k));
        ky(j, i, k) = 0.5 * (km + kp);
      });
  Kokkos::parallel_for(
      "gw_face_k_z",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz + 1}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        if (k == 0) {
          // Face above the k = 0 cell: zero, then the top-face rule below
          // (legacy Kz[kMou] = 0, groundwater.c:266).
          kzF(j, i, k) = 0.0;
          return;
        }
        if (k == nz) {
          // Bottom boundary face: the cell's own conductivity (legacy
          // icjckP > n3ci branch, :234); no-flux bottoms are sealed here
          // (:267).
          if (botCode(j, i) == static_cast<int>(GwBcCode::NoFlux)) {
            kzF(j, i, k) = 0.0;
          } else {
            kzF(j, i, k) =
                conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, nz - 1),
                                     h(j, i, nz - 1), ksz(j, i, nz - 1));
          }
          return;
        }
        // Interior face between cells k-1 (above) and k (below).
        const bool belowIsTop = (gid(j, i, k) >= 0) && (gid(j, i, k - 1) < 0);
        if (belowIsTop) {
          // The cell below is a column's top cell: one-sided from below
          // (legacy istop[icjckP] branch, :232). The seal pass zeroes this
          // again for masked interior columns, exactly as legacy does.
          kzF(j, i, k) = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                              h(j, i, k), ksz(j, i, k));
          return;
        }
        if (gid(j, i, k - 1) < 0) {
          kzF(j, i, k) = 0.0;
          return;
        }
        const real_t km = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k - 1),
                                               h(j, i, k - 1), ksz(j, i, k - 1));
        const real_t kp = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                               h(j, i, k), ksz(j, i, k));
        kzF(j, i, k) = 0.5 * (km + kp);
      });

  // Domain-edge side overrides (compute_K_face:222-262): no-flux edges are
  // sealed; prescribed-head edges take the interior cell's conductivity.
  Kokkos::parallel_for(
      "gw_face_k_edges",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        const auto applyEdge = [&](const Field3<real_t>& faces, int fj, int fi, int code,
                                   const Field3<real_t>& ksat) {
          if (code == static_cast<int>(GwBcCode::NoFlux)) {
            faces(fj, fi, k) = 0.0;
          } else if (code == static_cast<int>(GwBcCode::Head) ||
                     code == static_cast<int>(GwBcCode::HeadHydrostatic)) {
            faces(fj, fi, k) = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                                    h(j, i, k), ksat(j, i, k));
          }
        };
        if (gi == 0) {
          applyEdge(kx, j, i - 1, sideXm(j, i), ksx);
        }
        if (gi == nxGlobal - 1) {
          applyEdge(kx, j, i, sideXp(j, i), ksx);
        }
        if (gj == 0) {
          applyEdge(ky, j - 1, i, sideYm(j, i), ksy);
        }
        if (gj == nyGlobal - 1) {
          applyEdge(ky, j, i, sideYp(j, i), ksy);
        }
      });

  // The top-face rule (compute_K_face:269-292): prescribed head uses the
  // saturated conductivity, prescribed flux the cell's own K(h), no-flux
  // seals the face. Coupled runs replace the configured code (:276-284):
  // a wet surface averages the saturated and actual conductivities
  // (0.5 (Ksz + K(h)); the corrector uses pure Ksz — the preserved
  // predictor/corrector asymmetry), a dry surface is a seepage face with
  // the cell's own K(h). Legacy's scalar param->Ksz generalizes to the top
  // cell's own saturated conductivity (docs/theory/exchange-flux.md).
  const bool coupled = cpl_.active;
  Field2<real_t> cplDept = coupled ? cpl_.dept : Field2<real_t>();
  Kokkos::parallel_for(
      "gw_face_k_top",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const bool isTop = (gid(j, i, k) >= 0) && (k == 0 || gid(j, i, k - 1) < 0);
        if (!isTop) {
          return;
        }
        const real_t own = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                                h(j, i, k), ksz(j, i, k));
        if (coupled) {
          kzF(j, i, k) = (cplDept(j, i) > 0.0) ? 0.5 * (ksz(j, i, k) + own) : own;
          return;
        }
        const int code = topCode(j, i);
        if (code == static_cast<int>(GwBcCode::Head)) {
          kzF(j, i, k) = ksz(j, i, k);
        } else if (code == static_cast<int>(GwBcCode::NoFlux)) {
          kzF(j, i, k) = 0.0;
        } else {
          kzF(j, i, k) = own;
        }
      });

  // use_full3d = false: lateral conductivity vanishes at faces flanked by an
  // unsaturated active cell (compute_K_face:295-307; single-rank-equivalent
  // face-centric form, so every rank agrees on interface faces).
  if (!full3d) {
    Kokkos::parallel_for(
        "gw_face_k_full3d_x",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
            {1, 0, 0}, {nyl + 1, nxl + 1, nz}),
        KOKKOS_LAMBDA(const int j, const int i, const int k) {
          const auto unsat = [&](int jj, int ii) {
            return gid(jj, ii, k) >= 0 && wc(jj, ii, k) < wcs(jj, ii, k);
          };
          if (unsat(j, i) || unsat(j, i + 1)) {
            kx(j, i, k) = 0.0;
          }
        });
    Kokkos::parallel_for(
        "gw_face_k_full3d_y",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
            {0, 1, 0}, {nyl + 1, nxl + 1, nz}),
        KOKKOS_LAMBDA(const int j, const int i, const int k) {
          const auto unsat = [&](int jj, int ii) {
            return gid(jj, ii, k) >= 0 && wc(jj, ii, k) < wcs(jj, ii, k);
          };
          if (unsat(j, i) || unsat(j + 1, i)) {
            ky(j, i, k) = 0.0;
          }
        });
  }

  // Seal every face flanked by an inactive interior cell — LAST, exactly as
  // legacy orders it (compute_K_face:310-316). Face-centric form of the
  // per-cell legacy loop (identical single-rank semantics, race-free, and
  // rank-invariant because halo cells carry the neighbor's activity).
  Kokkos::parallel_for(
      "gw_face_k_seal_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 0, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const auto inactive = [&](int ii) {
          const int gi = i0 + ii - 1;
          return gi >= 0 && gi < nxGlobal && gid(j, ii, k) < 0;
        };
        if (inactive(i) || inactive(i + 1)) {
          kx(j, i, k) = 0.0;
        }
      });
  Kokkos::parallel_for(
      "gw_face_k_seal_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {0, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const auto inactive = [&](int jj) {
          const int gj = j0 + jj - 1;
          return gj >= 0 && gj < nyGlobal && gid(jj, i, k) < 0;
        };
        if (inactive(j) || inactive(j + 1)) {
          ky(j, i, k) = 0.0;
        }
      });
  Kokkos::parallel_for(
      "gw_face_k_seal_z",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz + 1}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const bool aboveInactive = (k >= 1) && (gid(j, i, k - 1) < 0);
        const bool belowInactive = (k < nz) && (gid(j, i, k) < 0);
        if (aboveInactive || belowInactive) {
          kzF(j, i, k) = 0.0;
        }
      });
}

void RichardsSolver::assembleSystem(real_t dtg) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t az = mesh_.areaZ();
  const real_t ss = ss_;
  const bool terrain = mesh_.followTerrain();
  const PetscInt offset = grid_.offset3();

  Field3<real_t> hn = hn_, wcn = wcn_, ch = ch_;
  Field3<real_t> kx = kx_, ky = ky_, kzF = kzF_;
  Field3<real_t> vga = vga_, vgn = vgn_, wcs = wcs_, wcr = wcr_, aev = aev_;
  Field3<real_t> rRho = rRho_, rXp = rRhoXp_, rYp = rRhoYp_, rZp = rRhoZp_;
  Field3<real_t> vXp = rViscXp_, vYp = rViscYp_, vZp = rViscZp_;
  Field3<real_t> dz3d = mesh_.dz3d(), bot3d = mesh_.bot3d();
  Field3<real_t> ax = mesh_.areaX(), ay = mesh_.areaY();
  Field3<real_t> sinx = mesh_.sinX(), cosx = mesh_.cosX();
  Field3<real_t> siny = mesh_.sinY(), cosy = mesh_.cosY();
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> topCode = topCode_, botCode = botCode_;
  Field2<int> sideXm = sideCodeXm_, sideXp = sideCodeXp_;
  Field2<int> sideYm = sideCodeYm_, sideYp = sideCodeYp_;
  Field2<real_t> topValue = topValue_, botValue = botValue_;
  Field2<real_t> sideVXm = sideValueXm_, sideVXp = sideValueXp_;
  Field2<real_t> sideVYm = sideValueYm_, sideVYp = sideValueYp_;
  Kokkos::View<real_t*, MemSpace> values = cooValues_;
  Kokkos::View<real_t*, MemSpace> rhs = rhsVec_;
  const bool coupled = cpl_.active;
  Field2<real_t> cplDept = coupled ? cpl_.dept : Field2<real_t>();
  Field2<real_t> cplAvail = coupled ? cpl_.avail : Field2<real_t>();
  Field2<int> cplMode = coupled ? cplMode_ : Field2<int>();

  // Specific capacity from the step-start head (legacy computes ch before
  // the predictor, groundwater.c:74).
  Kokkos::parallel_for(
      "gw_capacity",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        ch(j, i, k) =
            capacityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k), hn(j, i, k));
      });

  Kokkos::parallel_for(
      "gw_assemble",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const PetscInt row = gid(j, i, k);
        if (row < 0) {
          return;
        }
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        const bool isTop = (k == 0) || (gid(j, i, k - 1) < 0);
        const real_t volume = az * dz3d(j, i, k);

        // Off-diagonal coefficients (groundwater_mat_coeff:510-532).
        real_t gxp = -kx(j, i, k) * dtg * rXp(j, i, k) * vXp(j, i, k) * cosx(j, i, k) /
                     (dx * dx) * (ax(j, i, k) * dx);
        real_t gxm = -kx(j, i - 1, k) * dtg * rXp(j, i - 1, k) * vXp(j, i - 1, k) *
                     cosx(j, i - 1, k) / (dx * dx) * (ax(j, i - 1, k) * dx);
        real_t gyp = -ky(j, i, k) * dtg * rYp(j, i, k) * vYp(j, i, k) * cosy(j, i, k) /
                     (dy * dy) * (ay(j, i, k) * dy);
        real_t gym = -ky(j - 1, i, k) * dtg * rYp(j - 1, i, k) * vYp(j - 1, i, k) *
                     cosy(j - 1, i, k) / (dy * dy) * (ay(j - 1, i, k) * dy);
        // Prescribed-head edges sit at half spacing: the leg doubles
        // (groundwater_mat_coeff:519, :523; generalized to x edges).
        const auto isHead = [](int code) {
          return code == static_cast<int>(GwBcCode::Head) ||
                 code == static_cast<int>(GwBcCode::HeadHydrostatic);
        };
        if (gi == nxGlobal - 1 && isHead(sideXp(j, i))) {
          gxp *= 2.0;
        }
        if (gi == 0 && isHead(sideXm(j, i))) {
          gxm *= 2.0;
        }
        if (gj == nyGlobal - 1 && isHead(sideYp(j, i))) {
          gyp *= 2.0;
        }
        if (gj == 0 && isHead(sideYm(j, i))) {
          gym *= 2.0;
        }

        real_t dzf = (k == nz - 1) ? 0.5 * dz3d(j, i, k)
                                   : 0.5 * (dz3d(j, i, k) + dz3d(j, i, k + 1));
        const real_t gzp =
            -kzF(j, i, k + 1) * dtg * rZp(j, i, k + 1) * vZp(j, i, k + 1) /
            (dz3d(j, i, k) * dzf) * volume;
        dzf = isTop ? 0.5 * dz3d(j, i, k) : 0.5 * (dz3d(j, i, k) + dz3d(j, i, k - 1));
        const real_t gzm = -kzF(j, i, k) * dtg * rZp(j, i, k) * vZp(j, i, k) /
                           (dz3d(j, i, k) * dzf) * volume;

        // Diagonal (groundwater_mat_coeff:534-559).
        const real_t storage = (ch(j, i, k) + ss * wcn(j, i, k) / wcs(j, i, k)) * rRho(j, i, k);
        real_t gct = storage * volume - (gxp + gxm + gyp + gym);
        if (gi == nxGlobal - 1 && sideXp(j, i) == static_cast<int>(GwBcCode::Flux)) {
          gct += gxp;
        }
        if (gi == 0 && sideXm(j, i) == static_cast<int>(GwBcCode::Flux)) {
          gct += gxm;
        }
        if (gj == nyGlobal - 1 && sideYp(j, i) == static_cast<int>(GwBcCode::Flux)) {
          gct += gyp;
        }
        if (gj == 0 && sideYm(j, i) == static_cast<int>(GwBcCode::Flux)) {
          gct += gym;
        }
        if (!(k == nz - 1 && botCode(j, i) != static_cast<int>(GwBcCode::Head))) {
          gct -= gzp;
        }
        if (isTop) {
          // The top face couples into the diagonal only under a Dirichlet
          // head: the configured head condition, or a capacity-limited wet
          // surface in coupled runs (groundwater_mat_coeff:550-557;
          // supply-limited columns take a flux instead — amendment A13).
          const bool dirichletTop = coupled ? (cplMode(j, i) == 1)
                                            : (topCode(j, i) == static_cast<int>(GwBcCode::Head));
          if (dirichletTop) {
            gct -= gzm;
          }
        } else {
          gct -= gzm;
        }

        // Right-hand side (groundwater_rhs:578-675).
        real_t b = storage * hn(j, i, k) * volume;
        b -= dtg * volume * kzF(j, i, k + 1) * rZp(j, i, k + 1) * rZp(j, i, k + 1) *
             vZp(j, i, k + 1) / dz3d(j, i, k);
        b += dtg * volume * kzF(j, i, k) * rZp(j, i, k) * rZp(j, i, k) * vZp(j, i, k) /
             dz3d(j, i, k);
        if (terrain) {
          real_t sign = (bot3d(j, i + 1, k) > bot3d(j, i, k)) ? 1.0 : -1.0;
          b += sign * dtg * ax(j, i, k) * kx(j, i, k) * sinx(j, i, k);
          sign = (bot3d(j, i, k) > bot3d(j, i - 1, k)) ? -1.0 : 1.0;
          b += sign * dtg * ax(j, i - 1, k) * kx(j, i - 1, k) * sinx(j, i - 1, k);
          sign = (bot3d(j + 1, i, k) > bot3d(j, i, k)) ? 1.0 : -1.0;
          b += sign * dtg * ay(j, i, k) * ky(j, i, k) * siny(j, i, k);
          sign = (bot3d(j, i, k) > bot3d(j - 1, i, k)) ? -1.0 : 1.0;
          b += sign * dtg * ay(j - 1, i, k) * ky(j - 1, i, k) * siny(j - 1, i, k);
        }
        // Domain-edge legs: prescribed flux replaces the leg; anything else
        // folds the ghost head (zero-gradient copies for open/no-flux edges,
        // boundary values for head edges) — groundwater_rhs:601-626. The
        // side-flux sign convention is "positive into the domain".
        if (gi == nxGlobal - 1) {
          if (sideXp(j, i) == static_cast<int>(GwBcCode::Flux)) {
            b += rXp(j, i, k) * ax(j, i, k) * sideVXp(j, i) * dtg;
          } else {
            b -= gxp * hn(j, i + 1, k);
          }
        }
        if (gi == 0) {
          if (sideXm(j, i) == static_cast<int>(GwBcCode::Flux)) {
            b += rXp(j, i - 1, k) * ax(j, i - 1, k) * sideVXm(j, i) * dtg;
          } else {
            b -= gxm * hn(j, i - 1, k);
          }
        }
        if (gj == nyGlobal - 1) {
          if (sideYp(j, i) == static_cast<int>(GwBcCode::Flux)) {
            b += rYp(j, i, k) * ay(j, i, k) * sideVYp(j, i) * dtg;
          } else {
            b -= gyp * hn(j + 1, i, k);
          }
        }
        if (gj == 0) {
          if (sideYm(j, i) == static_cast<int>(GwBcCode::Flux)) {
            b += rYp(j - 1, i, k) * ay(j - 1, i, k) * sideVYm(j, i) * dtg;
          } else {
            b -= gym * hn(j - 1, i, k);
          }
        }
        // Bottom condition (groundwater_rhs:628-634 with the dtg/dz and
        // Dirichlet-sign fixes noted in the file header).
        if (k == nz - 1) {
          if (botCode(j, i) == static_cast<int>(GwBcCode::Flux)) {
            b += dtg * rZp(j, i, nz) * volume *
                 (botValue(j, i) + kzF(j, i, nz) * rZp(j, i, nz) * vZp(j, i, nz)) /
                 dz3d(j, i, k);
          } else if (botCode(j, i) == static_cast<int>(GwBcCode::Head)) {
            b -= gzp * botValue(j, i);
          }
        }
        // Top condition (groundwater_rhs:635-675). Coupled runs: a
        // capacity-limited wet surface is a Dirichlet head equal to the
        // local depth (:639-640); a supply-limited wet surface feeds the
        // remaining depth as a top flux (amendment A13; SERGHEI
        // GwBC.h:605-609); a dry surface takes the seepage-face form — the
        // flux expression with the legacy moisture guards on any configured
        // top flux (:641-659; qtop = 0 without one, leaving the pure
        // cancellation of the base gravity term).
        if (isTop) {
          if (coupled) {
            if (cplMode(j, i) == 1) {
              b -= gzm * cplDept(j, i);
            } else if (cplMode(j, i) == 2) {
              // Supply-limited flux (amendment A13): the exchange takes the
              // remaining surface water, bounded by the top cell's pore
              // room up to the legacy effective-saturation mark (0.9999
              // theta_s, the plan §3.1 named constant) so one large
              // adaptive substep can neither over-pressurize a filling
              // column nor land it exactly on the singular C(h) = 0 state
              // (SERGHEI's CFL-sized surface steps keep the same formula
              // bounded implicitly; our dtg does not).
              const real_t roomTop = Kokkos::fmax(
                  0.0, (0.9999 * wcs(j, i, k) - wcn(j, i, k)) * dz3d(j, i, k));
              const real_t take = Kokkos::fmin(cplAvail(j, i), roomTop);
              const real_t qSupply = -take / dtg;
              b += -dtg * rZp(j, i, k) * volume *
                   (qSupply + kzF(j, i, k) * rZp(j, i, k) * vZp(j, i, k)) / dz3d(j, i, k);
            } else {
              real_t qtop = 0.0;
              if (topCode(j, i) == static_cast<int>(GwBcCode::Flux)) {
                const real_t v = topValue(j, i);
                if ((v > 0.0 && wcn(j, i, k) < wcs(j, i, k)) ||
                    (v < 0.0 && wcn(j, i, k) > wcr(j, i, k))) {
                  qtop = v;
                }
              }
              b += -dtg * rZp(j, i, k) * volume *
                   (qtop + kzF(j, i, k) * rZp(j, i, k) * vZp(j, i, k)) / dz3d(j, i, k);
            }
          } else if (topCode(j, i) == static_cast<int>(GwBcCode::Head)) {
            b -= gzm * topValue(j, i);
          } else if (topCode(j, i) == static_cast<int>(GwBcCode::Flux)) {
            b += -dtg * rZp(j, i, k) * volume *
                 (topValue(j, i) + kzF(j, i, k) * rZp(j, i, k) * vZp(j, i, k)) /
                 dz3d(j, i, k);
          }
        }

        const std::size_t base = 7 * static_cast<std::size_t>(row - offset);
        values(base + 0) = gct;
        values(base + 1) = gxm;
        values(base + 2) = gxp;
        values(base + 3) = gym;
        values(base + 4) = gyp;
        values(base + 5) = gzm;
        values(base + 6) = gzp;
        rhs(static_cast<std::size_t>(row - offset)) = b;
      });
}

void RichardsSolver::fillAndSolve() {
  {
    Timer::Scoped timer("solve");
    system_->setValues(cooValues_);
    lastSolve_ = system_->solve(rhsVec_, solVec_);
  }
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const PetscInt offset = grid_.offset3();
  Field3<real_t> h = h_;
  Field3<PetscInt> gid = grid_.gid3();
  Kokkos::View<real_t*, MemSpace> sol = solVec_;
  Kokkos::parallel_for(
      "gw_scatter_head",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const PetscInt row = gid(j, i, k);
        if (row >= 0) {
          h(j, i, k) = sol(static_cast<std::size_t>(row - offset));
        }
      });
}

void RichardsSolver::enforceHeadBc() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const bool west = (grid_.rankWest() == MPI_PROC_NULL);
  const bool east = (grid_.rankEast() == MPI_PROC_NULL);
  const bool south = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool north = (grid_.rankNorth() == MPI_PROC_NULL);

  Field3<real_t> h = h_;
  Field3<real_t> dz3d = mesh_.dz3d(), bot3d = mesh_.bot3d();
  Field2<real_t> bath = mesh_.bath();
  Field3<real_t> rXp = rRhoXp_, rYp = rRhoYp_;
  const bool coupled = cpl_.active;
  Field2<real_t> dept = coupled ? cpl_.dept : Field2<real_t>();
  Field2<int> sideXm = sideCodeXm_, sideXp = sideCodeXp_;
  Field2<int> sideYm = sideCodeYm_, sideYp = sideCodeYp_;
  Field2<real_t> sideVXm = sideValueXm_, sideVXp = sideValueXp_;
  Field2<real_t> sideVYm = sideValueYm_, sideVYp = sideValueYp_;

  // Domain-edge ghost heads (enforce_head_bc:757-788): zero-gradient copies
  // by default; prescribed-head edges take the boundary value; hydrostatic
  // edges take the ghostHead hydrostatic rule (coupled runs follow the live
  // local surface, r-scaled — the legacy :769-788 branches — recomputed
  // from live state at every call so restarts reproduce them bitwise).
  // Corner halos are never read by the 7-point stencils and stay untouched.
  Kokkos::parallel_for(
      "gw_head_ghost_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nyl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int k) {
        if (west) {
          const real_t zc = bot3d(j, 1, k) + 0.5 * dz3d(j, 1, k);
          h(j, 0, k) = ghostHead(sideXm(j, 1), sideVXm(j, 1), h(j, 1, k), zc, coupled,
                                 bath(j, 1), coupled ? dept(j, 1) : 0.0, rXp(j, 0, k), false);
        }
        if (east) {
          const real_t zc = bot3d(j, nxl, k) + 0.5 * dz3d(j, nxl, k);
          h(j, nxl + 1, k) =
              ghostHead(sideXp(j, nxl), sideVXp(j, nxl), h(j, nxl, k), zc, coupled,
                        bath(j, nxl), coupled ? dept(j, nxl) : 0.0, rXp(j, nxl, k), true);
        }
      });
  Kokkos::parallel_for(
      "gw_head_ghost_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nxl + 1, nz}),
      KOKKOS_LAMBDA(const int i, const int k) {
        if (south) {
          const real_t zc = bot3d(1, i, k) + 0.5 * dz3d(1, i, k);
          h(0, i, k) = ghostHead(sideYm(1, i), sideVYm(1, i), h(1, i, k), zc, coupled,
                                 bath(1, i), coupled ? dept(1, i) : 0.0, rYp(0, i, k), false);
        }
        if (north) {
          const real_t zc = bot3d(nyl, i, k) + 0.5 * dz3d(nyl, i, k);
          h(nyl + 1, i, k) =
              ghostHead(sideYp(nyl, i), sideVYp(nyl, i), h(nyl, i, k), zc, coupled,
                        bath(nyl, i), coupled ? dept(nyl, i) : 0.0, rYp(nyl, i, k), true);
        }
      });
}

}  // namespace frehg::gw
