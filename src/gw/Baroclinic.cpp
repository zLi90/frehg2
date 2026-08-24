/// \file Baroclinic.cpp
/// \brief Density/viscosity ratios from the transport scalar and their face
///        means (plan §10 P4).
///
/// Provenance: update_rhovisc (scalar.c:921-955) and baroclinic_face
/// (groundwater.c:338-415). The cell ratios are r_rho = 1 + kBaroclinicBetaRho
/// s and r_visc = 1/(1 + kBaroclinicBetaVisc s) over every cell including
/// ghosts (the ghost scalar carries the side Dirichlet values the transport
/// module enforces, so a saline head boundary densifies its boundary face).
/// Face ratios are arithmetic means between active flanks with the legacy
/// edge rules preserved:
///  - minus-side domain-edge faces average the interior cell with the ghost
///    (legacy iMou/jMou loops, groundwater.c:395-411);
///  - plus-side domain-edge faces stay at 1 unless the edge holds a head
///    condition, in which case both edge kinds take the ghost's own ratio
///    one-sided (legacy :361-364 for y+, :407-410 for y-);
///  - the face above a column's top cell averages the top cell with the
///    ghost above it — the surface scalar in coupled runs, the top cell
///    itself otherwise (legacy :387-389 with the enforce_scalar_bc ghost,
///    scalar.c:897-909);
///  - the bottom boundary face stays at 1 (legacy never writes it).
/// Legacy defined the head-BC one-sided rule only for the y faces because
/// its x sides had no head conditions at all; Frehg2's x rules mirror the y
/// rules exactly, matching the P2 x-side generalization
/// (docs/theory/groundwater.md).

#include "gw/RichardsSolver.hpp"

namespace frehg::gw {

namespace {

KOKKOS_INLINE_FUNCTION real_t rhoOf(real_t s) { return 1.0 + kBaroclinicBetaRho * s; }
KOKKOS_INLINE_FUNCTION real_t viscOf(real_t s) { return 1.0 / (1.0 + kBaroclinicBetaVisc * s); }

KOKKOS_INLINE_FUNCTION bool isHeadCode(int code) {
  return code == static_cast<int>(GwBcCode::Head) ||
         code == static_cast<int>(GwBcCode::HeadHydrostatic);
}

}  // namespace

void RichardsSolver::updateBaroclinicFaces() {
  if (!baroclinic_ || sSubs_.size() == 0) {
    return;
  }
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const bool coupled = cpl_.active;

  Field3<real_t> s = sSubs_;
  Field2<real_t> sSurf = sSurfGhost_;
  Field3<real_t> rRho = rRho_, rXp = rRhoXp_, rYp = rRhoYp_, rZp = rRhoZp_;
  Field3<real_t> rVisc = rVisc_, vXp = rViscXp_, vYp = rViscYp_, vZp = rViscZp_;
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> sideXm = sideCodeXm_, sideXp = sideCodeXp_;
  Field2<int> sideYm = sideCodeYm_, sideYp = sideCodeYp_;

  // Cell ratios over the extended box (update_rhovisc runs over n3ct).
  Kokkos::parallel_for(
      "gw_baroclinic_cells",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {0, 0, 0}, {nyl + 2, nxl + 2, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        rRho(j, i, k) = rhoOf(s(j, i, k));
        rVisc(j, i, k) = viscOf(s(j, i, k));
      });

  // x faces: slot (j, i, k) is the x+ face of cell i; slot i = 0 the west
  // boundary/interface face (the corrector's flux layout).
  Kokkos::parallel_for(
      "gw_baroclinic_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 0, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const int giMinus = i0 + i - 1;
        real_t rr = 1.0;
        real_t rv = 1.0;
        if (giMinus < 0) {
          // West domain edge (legacy iMou rule; head override mirrors y-).
          if (gid(j, 1, k) >= 0) {
            if (isHeadCode(sideXm(j, 1))) {
              rr = rRho(j, 0, k);
              rv = rVisc(j, 0, k);
            } else {
              rr = 0.5 * (rRho(j, 1, k) + rRho(j, 0, k));
              rv = 0.5 * (rVisc(j, 1, k) + rVisc(j, 0, k));
            }
          }
        } else if (giMinus == nxGlobal - 1) {
          // East domain edge (legacy y+ rule): ghost value under a head
          // condition, otherwise 1 (the actv pair guard fails there).
          if (gid(j, i, k) >= 0 && isHeadCode(sideXp(j, i))) {
            rr = rRho(j, i + 1, k);
            rv = rVisc(j, i + 1, k);
          }
        } else if (gid(j, i, k) >= 0 && gid(j, i + 1, k) >= 0) {
          rr = 0.5 * (rRho(j, i, k) + rRho(j, i + 1, k));
          rv = 0.5 * (rVisc(j, i, k) + rVisc(j, i + 1, k));
        }
        rXp(j, i, k) = rr;
        vXp(j, i, k) = rv;
      });

  // y faces (legacy :361-375 and the jMou loop :401-411).
  Kokkos::parallel_for(
      "gw_baroclinic_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {0, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const int gjMinus = j0 + j - 1;
        real_t rr = 1.0;
        real_t rv = 1.0;
        if (gjMinus < 0) {
          if (gid(1, i, k) >= 0) {
            if (isHeadCode(sideYm(1, i))) {
              rr = rRho(0, i, k);
              rv = rVisc(0, i, k);
            } else {
              rr = 0.5 * (rRho(1, i, k) + rRho(0, i, k));
              rv = 0.5 * (rVisc(1, i, k) + rVisc(0, i, k));
            }
          }
        } else if (gjMinus == nyGlobal - 1) {
          if (gid(j, i, k) >= 0 && isHeadCode(sideYp(j, i))) {
            rr = rRho(j + 1, i, k);
            rv = rVisc(j + 1, i, k);
          }
        } else if (gid(j, i, k) >= 0 && gid(j + 1, i, k) >= 0) {
          rr = 0.5 * (rRho(j, i, k) + rRho(j + 1, i, k));
          rv = 0.5 * (rVisc(j, i, k) + rVisc(j + 1, i, k));
        }
        rYp(j, i, k) = rr;
        vYp(j, i, k) = rv;
      });

  // z faces: plane k is the face above cell k, plane nz the bottom boundary
  // face (stays 1, legacy never writes it). The plane above each column's
  // top cell averages against the surface scalar in coupled runs (the
  // legacy kM-ghost scalar, scalar.c:901) and collapses one-sided otherwise
  // (zero-gradient ghost).
  Kokkos::parallel_for(
      "gw_baroclinic_z",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz + 1}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const bool belowActive = (k < nz) && (gid(j, i, k) >= 0);
        const bool aboveActive = (k >= 1) && (gid(j, i, k - 1) >= 0);
        real_t rr = 1.0;
        real_t rv = 1.0;
        if (belowActive && aboveActive) {
          rr = 0.5 * (rRho(j, i, k - 1) + rRho(j, i, k));
          rv = 0.5 * (rVisc(j, i, k - 1) + rVisc(j, i, k));
        } else if (belowActive) {
          // Top face of the column (legacy istop branch, :387-389).
          const real_t sGhost = coupled ? sSurf(j, i) : s(j, i, k);
          rr = 0.5 * (rRho(j, i, k) + rhoOf(sGhost));
          rv = 0.5 * (rVisc(j, i, k) + viscOf(sGhost));
        }
        rZp(j, i, k) = rr;
        vZp(j, i, k) = rv;
      });
}

}  // namespace frehg::gw
