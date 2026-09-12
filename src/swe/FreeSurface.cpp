/// \file FreeSurface.cpp
/// \brief The implicit 5-point free-surface system: boundary enforcement,
///        right-hand side, matrix coefficients, COO fill, and solve.
///
/// Provenance: enforce_surf_bc (shallowwater.c:502-540), shallowwater_rhs
/// (shallowwater.c:296-317), shallowwater_mat_coeff (shallowwater.c:320-413),
/// build/solve_shallowwater_system (shallowwater.c:416-499). Deviations
/// (documented in docs/theory/surface-water.md): the system is one global
/// PETSc solve — legacy solved per-rank blocks with the neighbor's previous
/// eta as Dirichlet data, which cannot be rank-invariant (plan §8.2); eta-
/// condition columns are eliminated into the right-hand side so the matrix
/// stays SPD (the solution is identical); a discharge condition divides by
/// its global member count — legacy divided by the per-rank count
/// (shallowwater.c:313), which multiplies the inflow when the region spans
/// ranks.

#include "swe/SurfaceSolver.hpp"

#include "core/Logger.hpp"
#include "core/Timer.hpp"

namespace frehg::swe {

void SurfaceSolver::enforceSurfBc(bool prescribeStage) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t cellArea = grid_.dx() * grid_.dy();
  Field2<real_t> eta = eta_, bottom = bottom_;

  // Surface below the bed is lifted onto it (shallowwater.c:507-511); the
  // lift creates volume with no compensating flux, so it is measured into
  // the step audit (the legacy defect stays, its size becomes data).
  real_t clamped = 0.0;
  Kokkos::parallel_reduce(
      "swe_eta_clamp",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& acc) {
        if (eta(j, i) < bottom(j, i)) {
          acc += (bottom(j, i) - eta(j, i)) * cellArea;
          eta(j, i) = bottom(j, i);
        }
      },
      clamped);
  audit_.clampVolume += clamped;

  // Prescribed stage on eta-condition cells, after the clamp so a stage
  // below the bed survives (shallowwater.c:513-525). Skipped on restart
  // refreshes (see the header note).
  if (prescribeStage) {
    for (const DeviceBcList& list : etaBcs_) {
      auto lj = list.j;
      auto li = list.i;
      const real_t current = list.current;
      Kokkos::parallel_for(
          "swe_eta_bc_enforce",
          Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
          KOKKOS_LAMBDA(const std::size_t m) { eta(lj(m), li(m)) = current; });
    }
  }

  // Zero-gradient ghosts along the physical domain edges
  // (shallowwater.c:527-539).
  const bool westEdge = (grid_.rankWest() == MPI_PROC_NULL);
  const bool eastEdge = (grid_.rankEast() == MPI_PROC_NULL);
  const bool southEdge = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool northEdge = (grid_.rankNorth() == MPI_PROC_NULL);
  Kokkos::parallel_for(
      "swe_eta_edge_ghosts_i",
      Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nyl + 1),
      KOKKOS_LAMBDA(const int j) {
        if (westEdge) {
          eta(j, 0) = eta(j, 1);
        }
        if (eastEdge) {
          eta(j, nxl + 1) = eta(j, nxl);
        }
      });
  Kokkos::parallel_for(
      "swe_eta_edge_ghosts_j",
      Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nxl + 1),
      KOKKOS_LAMBDA(const int i) {
        if (southEdge) {
          eta(0, i) = eta(1, i);
        }
        if (northEdge) {
          eta(nyl + 1, i) = eta(nyl, i);
        }
      });

  // Free-outflow faces extrapolate the stage down the continued bed slope:
  // the ghost surface sits one bed drop below the cell, so the boundary
  // face keeps the interior momentum balance instead of the closed-edge
  // zero gradient (kind outflow, plan §5.6 amendment; decided by the P1 b4
  // gate).
  for (const DeviceBcList& list : outflowBcs_) {
    auto lj = list.j;
    auto li = list.i;
    auto lface = list.face;
    auto drop = list.bedDrop;
    Kokkos::parallel_for(
        "swe_outflow_ghost_eta",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m) {
          const int j = lj(m);
          const int i = li(m);
          const real_t ghost = eta(j, i) - drop(m);
          switch (static_cast<BcFace>(lface(m))) {
            case BcFace::XPlus:
              eta(j, i + 1) = ghost;
              break;
            case BcFace::XMinus:
              eta(j, i - 1) = ghost;
              break;
            case BcFace::YPlus:
              eta(j + 1, i) = ghost;
              break;
            case BcFace::YMinus:
              eta(j - 1, i) = ghost;
              break;
            default:
              break;
          }
        });
  }
}

void SurfaceSolver::assembleRhs() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t dt = dt_;
  Field2<real_t> Srhs = Srhs_, eta = eta_, Asz = Asz_, Asx = Asx_, Asy = Asy_, Ex = Ex_, Ey = Ey_;

  Kokkos::parallel_for(
      "swe_rhs",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        Srhs(j, i) = eta(j, i) * Asz(j, i) -
                     dt * (Asx(j, i) * Ex(j, i) - Asx(j, i - 1) * Ex(j, i - 1) +
                           Asy(j, i) * Ey(j, i) - Asy(j - 1, i) * Ey(j - 1, i));
      });

  // Discharge conditions distribute the inflow volume over the region's
  // global member count (shallowwater.c:306-316; global-count fix above).
  for (const DeviceBcList& list : dischargeBcs_) {
    auto lj = list.j;
    auto li = list.i;
    const real_t perCell = list.current * dt / static_cast<real_t>(list.bc->globalCellCount());
    real_t added = 0.0;
    Kokkos::parallel_reduce(
        "swe_discharge_bc",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m, real_t& sum) {
          Srhs(lj(m), li(m)) += perCell;
          sum += perCell;
        },
        added);
    audit_.bcInflow += added;
  }

  // Velocity conditions prescribe the normal velocity of a domain-edge face
  // (sign convention of uu/vv). The face flux replaces the default edge
  // treatment: open east/north faces carry dt*As*E, west/south faces carry
  // nothing.
  for (const DeviceBcList& list : velocityBcs_) {
    auto lj = list.j;
    auto li = list.i;
    auto lface = list.face;
    const real_t u = list.current;
    real_t added = 0.0;
    Kokkos::parallel_reduce(
        "swe_velocity_bc_rhs",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m, real_t& sum) {
          const int j = lj(m);
          const int i = li(m);
          real_t add = 0.0;
          switch (static_cast<BcFace>(lface(m))) {
            case BcFace::XMinus:
              add = dt * Asx(j, 0) * u;
              break;
            case BcFace::XPlus:
              add = dt * Asx(j, i) * (Ex(j, i) - u);
              break;
            case BcFace::YMinus:
              add = dt * Asy(0, i) * u;
              break;
            case BcFace::YPlus:
              add = dt * Asy(j, i) * (Ey(j, i) - u);
              break;
            default:
              break;
          }
          Srhs(j, i) += add;
          sum += add;
        },
        added);
    audit_.bcInflow += added;
  }
}

void SurfaceSolver::assembleCoefficients() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const real_t coef = gravity_ * dt_ * dt_;
  const real_t cellArea = grid_.dx() * grid_.dy();
  Field2<real_t> Sxp = Sxp_, Sxm = Sxm_, Syp = Syp_, Sym = Sym_, Sct = Sct_, Srhs = Srhs_;
  Field2<real_t> Asx = Asx_, Asy = Asy_, Asz = Asz_, Vsx = Vsx_, Vsy = Vsy_;
  Field2<real_t> Dx = Dx_, Dy = Dy_;
  Field2<real_t> dept = dept_, eta = eta_;

  Kokkos::parallel_for(
      "swe_mat_coeff",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        real_t sxp = 0.0, sxm = 0.0, syp = 0.0, sym = 0.0;
        if (Vsx(j, i) > 0.0) {
          sxp = coef * Asx(j, i) * Asx(j, i) * Dx(j, i) / Vsx(j, i);
        }
        if (Vsx(j, i - 1) > 0.0) {
          sxm = coef * Asx(j, i - 1) * Asx(j, i - 1) * Dx(j, i - 1) / Vsx(j, i - 1);
        }
        if (Vsy(j, i) > 0.0) {
          syp = coef * Asy(j, i) * Asy(j, i) * Dy(j, i) / Vsy(j, i);
        }
        if (Vsy(j - 1, i) > 0.0) {
          sym = coef * Asy(j - 1, i) * Asy(j - 1, i) * Dy(j - 1, i) / Vsy(j - 1, i);
        }
        real_t sct = Asz(j, i) + sxp + sxm + syp + sym;

        // Dry-cell closure (amendment A13; legacy shallowwater.c:349-361
        // replaced). Legacy held dry rows at eta with the off-diagonal legs
        // kept whenever a face velocity was nonzero and the right-hand side
        // overwritten to eta dx dy: the held row then reads
        // eta = eta_old + (S/A) eta_neighbor — a frame-dependent pull on
        // the neighbors' *absolute* stage that both mints volume out of
        // the datum and (velocity-zero case) breaks the symmetry CG
        // requires. Below the b1/b4 tolerance floor there (legs / area
        // ~ 1e-4), it floods b5's wall-adjacent column within steps and
        // stalls the solver. The closure here treats a dry cell as a
        // zero-depth continuity row: its free-surface area is the cell
        // area (dV/d eta at the bed) and its legs stay live. The face
        // depths built from the higher-of-two-bottoms rule vanish whenever
        // the neighbors' stage sits below this cell's bed, so live legs
        // can only carry flux *into* the dry cell — wetting is implicit,
        // conservative, symmetric, and translation-invariant.
        if (dept(j, i) == 0.0) {
          sct = cellArea + sxp + sxm + syp + sym;
          Srhs(j, i) += eta(j, i) * cellArea;
        }

        // Fold the closed-boundary legs into the diagonal at the global
        // domain edges (shallowwater.c:367-393): the ghost surface is a
        // zero-gradient copy, so the implicit face flux vanishes.
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        if (gi == 0) {
          sct -= sxm;
        }
        if (gi == nxGlobal - 1) {
          sct -= sxp;
        }
        if (gj == 0) {
          sct -= sym;
        }
        if (gj == nyGlobal - 1) {
          sct -= syp;
        }

        Sxp(j, i) = sxp;
        Sxm(j, i) = sxm;
        Syp(j, i) = syp;
        Sym(j, i) = sym;
        Sct(j, i) = sct;
      });
}

void SurfaceSolver::applyOutflowCorrections() {
  // The transmissive ghost stage tracks the cell stage (eta_ghost' =
  // eta_cell' - drop), so the boundary face's implicit coupling cancels and
  // only the constant S_face * drop survives as an explicit outflow volume.
  // The face coefficient mirrors the interior form with the cell's drag
  // factor; west/south faces use the halo-column face geometry.
  const real_t coef = gravity_ * dt_ * dt_;
  Field2<real_t> Srhs = Srhs_, Sxp = Sxp_, Syp = Syp_;
  Field2<real_t> Asx = Asx_, Asy = Asy_, Vsx = Vsx_, Vsy = Vsy_, Dx = Dx_, Dy = Dy_;
  for (const DeviceBcList& list : outflowBcs_) {
    auto lj = list.j;
    auto li = list.i;
    auto lface = list.face;
    auto drop = list.bedDrop;
    real_t released = 0.0;
    Kokkos::parallel_reduce(
        "swe_outflow_rhs",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m, real_t& sum) {
          const int j = lj(m);
          const int i = li(m);
          real_t s = 0.0;
          switch (static_cast<BcFace>(lface(m))) {
            case BcFace::XPlus:
              s = Sxp(j, i);
              break;
            case BcFace::XMinus:
              s = (Vsx(j, 0) > 0.0) ? coef * Asx(j, 0) * Asx(j, 0) * Dx(j, i) / Vsx(j, 0) : 0.0;
              break;
            case BcFace::YPlus:
              s = Syp(j, i);
              break;
            case BcFace::YMinus:
              s = (Vsy(0, i) > 0.0) ? coef * Asy(0, i) * Asy(0, i) * Dy(j, i) / Vsy(0, i) : 0.0;
              break;
            default:
              break;
          }
          const real_t volume = s * drop(m);
          Srhs(j, i) -= volume;
          sum += volume;
        },
        released);
    audit_.boundaryOutflow += released;
  }
}

void SurfaceSolver::fillAndSolve() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const PetscInt offset = grid_.offset2();
  Field2<PetscInt> gid = grid_.gid2();
  Field2<real_t> Sxp = Sxp_, Sxm = Sxm_, Syp = Syp_, Sym = Sym_, Sct = Sct_, Srhs = Srhs_;
  Field2<real_t> isEtaBc = isEtaBc_, etaBcValue = etaBcValue_, eta = eta_;
  Kokkos::View<real_t*, MemSpace> values = cooValues_;
  Kokkos::View<real_t*, MemSpace> rhs = rhsVec_;
  Kokkos::View<real_t*, MemSpace> sol = solVec_;

  // Wet/dry-mask rebuild trigger (v2 plan §2.2.4): a mask flip is a
  // structural change of the operator (the A13 dry-cell closure swaps a
  // row's diagonal between Asz and the cell area), which a frozen AMG
  // hierarchy does not represent. The hash is Allreduced so every rank
  // reaches the same rebuild decision (PCSetUp is collective).
  if (system_->reusesHierarchy()) {
    Field2<real_t> dept = dept_;
    long long localHash = 0;
    Kokkos::parallel_reduce(
        "swe_wet_mask_hash",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i, long long& sum) {
          if (dept(j, i) > 0.0) {
            sum += static_cast<long long>(gid(j, i)) + 1;
          }
        },
        localHash);
    long long hash = 0;
    MPI_Allreduce(&localHash, &hash, 1, MPI_LONG_LONG, MPI_SUM, grid_.comm());
    if (hash != wetMaskHash_) {
      wetMaskHash_ = hash;
      system_->forceRebuild();
    }
  }

  const real_t cellArea = grid_.dx() * grid_.dy();
  Kokkos::parallel_for(
      "swe_coo_fill",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const std::size_t row = static_cast<std::size_t>(gid(j, i) - offset);
        const std::size_t base = 5 * row;
        if (isEtaBc(j, i) != 0.0) {
          // Prescribed-stage row: decoupled identity (shallowwater.c:396-411),
          // scaled to the cell area so its diagonal matches the wet rows'
          // magnitude — a literal 1.0 among ~dx*dy diagonals is the poorly
          // scaled identity row AMG coarsening mishandles (v2 plan §2.2.5).
          // The row solution eta = value is unchanged by the scaling.
          values(base + 0) = cellArea;
          values(base + 1) = 0.0;
          values(base + 2) = 0.0;
          values(base + 3) = 0.0;
          values(base + 4) = 0.0;
          rhs(row) = etaBcValue(j, i) * cellArea;
          return;
        }
        real_t b = Srhs(j, i);
        real_t vxm = -Sxm(j, i);
        real_t vxp = -Sxp(j, i);
        real_t vym = -Sym(j, i);
        real_t vyp = -Syp(j, i);
        // Eliminate prescribed-stage columns into the right-hand side so the
        // matrix stays SPD; the linear solution is unchanged.
        if (isEtaBc(j, i - 1) != 0.0) {
          b += Sxm(j, i) * etaBcValue(j, i - 1);
          vxm = 0.0;
        }
        if (isEtaBc(j, i + 1) != 0.0) {
          b += Sxp(j, i) * etaBcValue(j, i + 1);
          vxp = 0.0;
        }
        if (isEtaBc(j - 1, i) != 0.0) {
          b += Sym(j, i) * etaBcValue(j - 1, i);
          vym = 0.0;
        }
        if (isEtaBc(j + 1, i) != 0.0) {
          b += Syp(j, i) * etaBcValue(j + 1, i);
          vyp = 0.0;
        }
        values(base + 0) = Sct(j, i);
        values(base + 1) = vxm;
        values(base + 2) = vxp;
        values(base + 3) = vym;
        values(base + 4) = vyp;
        rhs(row) = b;
      });

  {
    Timer::Scoped timer("solve");
    system_->setValues(cooValues_);
    lastSolve_ = system_->solve(rhsVec_, solVec_);
  }

  Kokkos::parallel_for(
      "swe_scatter_eta",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        eta(j, i) = sol(static_cast<std::size_t>(gid(j, i) - offset));
      });
}

void SurfaceSolver::accumulateBoundaryFluxes() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const real_t dt = dt_;
  Field2<real_t> Asx = Asx_, Asy = Asy_, Ex = Ex_, Ey = Ey_;
  Field2<real_t> Sxp = Sxp_, Sxm = Sxm_, Syp = Syp_, Sym = Sym_;
  Field2<real_t> isEtaBc = isEtaBc_, etaBcValue = etaBcValue_, eta = eta_;

  // Volume leaving the domain this step through open edges (the explicit
  // face flux is the only transport there once the implicit leg is folded)
  // and across faces into prescribed-stage cells, which act as external
  // reservoirs. Face fluxes mirror the assembled continuity equation, so
  // interior faces cancel exactly. All four domain edges are counted for
  // completeness (P3): under the preserved legacy asymmetry the ghost-side
  // E at west/south faces is identically zero (docs/theory/
  // surface-water.md, "West/south closed edges carry no explicit flux"),
  // so those terms measure zero today and become live only if a later
  // boundary kind fills the ghost momentum.
  real_t outflow = 0.0;
  Kokkos::parallel_reduce(
      "swe_boundary_outflow",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sum) {
        if (isEtaBc(j, i) != 0.0) {
          return;
        }
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        if (gi == 0) {
          sum -= dt * Asx(j, i - 1) * Ex(j, i - 1);
        }
        if (gi == nxGlobal - 1) {
          sum += dt * Asx(j, i) * Ex(j, i);
        }
        if (gj == 0) {
          sum -= dt * Asy(j - 1, i) * Ey(j - 1, i);
        }
        if (gj == nyGlobal - 1) {
          sum += dt * Asy(j, i) * Ey(j, i);
        }
        if (isEtaBc(j, i + 1) != 0.0) {
          sum += dt * Asx(j, i) * Ex(j, i) - Sxp(j, i) * (etaBcValue(j, i + 1) - eta(j, i));
        }
        if (isEtaBc(j, i - 1) != 0.0) {
          sum += -dt * Asx(j, i - 1) * Ex(j, i - 1) -
                 Sxm(j, i) * (etaBcValue(j, i - 1) - eta(j, i));
        }
        if (isEtaBc(j + 1, i) != 0.0) {
          sum += dt * Asy(j, i) * Ey(j, i) - Syp(j, i) * (etaBcValue(j + 1, i) - eta(j, i));
        }
        if (isEtaBc(j - 1, i) != 0.0) {
          sum += -dt * Asy(j - 1, i) * Ey(j - 1, i) -
                 Sym(j, i) * (etaBcValue(j - 1, i) - eta(j, i));
        }
      },
      outflow);
  audit_.boundaryOutflow += outflow;
}

}  // namespace frehg::swe
