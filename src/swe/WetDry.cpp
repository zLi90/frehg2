/// \file WetDry.cpp
/// \brief Depth and face-geometry updates, wetting/drying limiters, and
///        velocity boundary handling.
///
/// Provenance: update_depth (initialize.c:918-995), the plain-geometry
/// branch of the legacy cell-volume/face-area update
/// (shallowwater.c:1046-1099) and its ic_surface twin
/// (initialize.c:584-611), cfl_limiter (shallowwater.c:543-574), the
/// velocity limiters of update_velocity (shallowwater.c:796-833), and
/// enforce_velo_bc (shallowwater.c:891-959).
///
/// The velocity limiters are face-centric here: legacy iterated cells and
/// wrote to the west/south neighbor's face; because every write stores zero
/// and the trigger conditions read strict signs, the per-face formulation
/// reproduces the sequential result exactly while staying deterministic in
/// parallel. Legacy's zeroing of physical-edge ghost faces was dead code
/// (enforce_velo_bc overwrites those ghosts with interior copies right
/// after).

#include "swe/SurfaceSolver.hpp"

namespace frehg::swe {

void SurfaceSolver::updateDepth() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t minDepth = minDepth_;
  Field2<real_t> eta = eta_, bottom = bottom_, dept = dept_, deptx = deptx_, depty = depty_;

  // Legacy update_depth first repeats the small-depth removal; every caller
  // in the P1 flow reaches it with eta already clamped (evaprain's final
  // clamp, initialize.c:627-635), so only the depth evaluation remains.
  Kokkos::parallel_for(
      "swe_depth_center",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 0}, {nyl + 2, nxl + 2}),
      KOKKOS_LAMBDA(const int j, const int i) {
        real_t d = eta(j, i) - bottom(j, i);
        if (d <= minDepth) {
          d = 0.0;
        }
        dept(j, i) = d;
      });

  // Face depths from the higher of the two surfaces over the higher of the
  // two bottoms (initialize.c:937-953), clamped at zero. The x-face column
  // i = 0 and y-face row j = 0 cover the halo-adjacent faces; they are
  // consistent across ranks because eta halos were exchanged.
  Kokkos::parallel_for(
      "swe_depth_faces_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const real_t etaHi = Kokkos::fmax(eta(j, i), eta(j, i + 1));
        const real_t botHi = Kokkos::fmax(bottom(j, i), bottom(j, i + 1));
        deptx(j, i) = Kokkos::fmax(etaHi - botHi, 0.0);
      });
  Kokkos::parallel_for(
      "swe_depth_faces_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const real_t etaHi = Kokkos::fmax(eta(j, i), eta(j + 1, i));
        const real_t botHi = Kokkos::fmax(bottom(j, i), bottom(j + 1, i));
        depty(j, i) = Kokkos::fmax(etaHi - botHi, 0.0);
      });
}

void SurfaceSolver::updateGeometry() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const bool oneColumn = (grid_.nx() == 1);
  const bool oneRow = (grid_.ny() == 1);
  Field2<real_t> dept = dept_, deptx = deptx_, depty = depty_;
  Field2<real_t> Vs = Vs_, Vsx = Vsx_, Vsy = Vsy_;
  Field2<real_t> Asx = Asx_, Asy = Asy_, Asz = Asz_, Aszx = Aszx_, Aszy = Aszy_;

  // Cell volumes and face areas over the full extended box; every value is
  // a pointwise function of exchanged data, so interface halos need no
  // second exchange (shallowwater.c:1048-1058).
  Kokkos::parallel_for(
      "swe_geometry_cells",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 0}, {nyl + 2, nxl + 2}),
      KOKKOS_LAMBDA(const int j, const int i) {
        Vs(j, i) = dept(j, i) * dx * dy;
        Asz(j, i) = (dept(j, i) > 0.0) ? dx * dy : 0.0;
        Asx(j, i) = oneColumn ? 0.0 : deptx(j, i) * dy;
        Asy(j, i) = oneRow ? 0.0 : depty(j, i) * dx;
      });

  // Face volumes (shallowwater.c:1059-1065), extended one halo column/row
  // so the west/south matrix legs have their face quantities locally.
  Kokkos::parallel_for(
      "swe_geometry_faces_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        Vsx(j, i) = 0.5 * (Vs(j, i) + Vs(j, i + 1));
        Aszx(j, i) = 0.5 * (Asz(j, i) + Asz(j, i + 1));
      });
  Kokkos::parallel_for(
      "swe_geometry_faces_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        Vsy(j, i) = 0.5 * (Vs(j, i) + Vs(j + 1, i));
        Aszy(j, i) = 0.5 * (Asz(j, i) + Asz(j + 1, i));
      });

  // Physical-edge ghost rules (shallowwater.c:1067-1099): interior copies
  // for the areas, and the legacy full-volume (not averaged) assignments
  // for the west/south face volumes.
  const bool westEdge = (grid_.rankWest() == MPI_PROC_NULL);
  const bool southEdge = (grid_.rankSouth() == MPI_PROC_NULL);
  if (westEdge) {
    Kokkos::parallel_for(
        "swe_geometry_west_edge",
        Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nyl + 1),
        KOKKOS_LAMBDA(const int j) {
          Vs(j, 0) = Vs(j, 1);
          Asx(j, 0) = Asx(j, 1);
          Asy(j, 0) = Asy(j, 1);
          Asz(j, 0) = Asz(j, 1);
          Vsx(j, 0) = Vs(j, 1);
          Aszx(j, 0) = Asz(j, 1);
        });
  }
  if (southEdge) {
    Kokkos::parallel_for(
        "swe_geometry_south_edge",
        Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nxl + 1),
        KOKKOS_LAMBDA(const int i) {
          Vs(0, i) = Vs(1, i);
          Asx(0, i) = Asx(1, i);
          Asy(0, i) = Asy(1, i);
          Asz(0, i) = Asz(1, i);
          Vsy(0, i) = Vs(1, i);
          Aszy(0, i) = Asz(1, i);
        });
  }
}

void SurfaceSolver::cflLimiter() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t minDepth = minDepth_;
  Field2<real_t> eta = eta_, bottom = bottom_, dept = dept_, cflActive = cflActive_;

  // Wetting is restricted to cells adjacent to water: an isolated cell that
  // the solve wetted is forced dry and its faces flagged for zeroing
  // (shallowwater.c:548-567). dept still holds the start-of-step state.
  Kokkos::parallel_for(
      "swe_cfl_limiter_isolated",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const real_t diff = eta(j, i) - bottom(j, i);
        if (dept(j, i) <= 0.0 && diff > 0.0) {
          const bool wet = dept(j, i + 1) > 0.0 || dept(j, i - 1) > 0.0 ||
                           dept(j + 1, i) > 0.0 || dept(j - 1, i) > 0.0;
          if (!wet) {
            eta(j, i) = bottom(j, i);
            cflActive(j, i) = 1.0;
          }
        }
      });

  // Depths thinner than min_depth are dried (shallowwater.c:569-573).
  Kokkos::parallel_for(
      "swe_cfl_limiter_small",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const real_t diff = eta(j, i) - bottom(j, i);
        if (diff > 0.0 && diff < minDepth) {
          eta(j, i) = bottom(j, i);
        }
      });
}

void SurfaceSolver::applyVelocityLimiters() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t dt = dt_;
  const real_t wtfh = wettingFaceDepth_;
  Field2<real_t> uu = uu_, vv = vv_, dept = dept_, cflActive = cflActive_;
  Field2<real_t> Asx = Asx_, Asy = Asy_, Fu = Fu_, Fv = Fv_, cflx = cflx_, cfly = cfly_;

  Kokkos::parallel_for(
      "swe_velocity_limiters",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        // Zero velocity across vanished faces (shallowwater.c:799-800).
        if (Asx(j, i) < wtfh * dy) {
          uu(j, i) = 0.0;
        }
        if (Asy(j, i) < wtfh * dx) {
          vv(j, i) = 0.0;
        }
        // Zero outflow from effectively dry cells (shallowwater.c:802-808).
        if (dept(j, i) < wtfh && uu(j, i) > 0.0) {
          uu(j, i) = 0.0;
        }
        if (dept(j, i + 1) < wtfh && uu(j, i) < 0.0) {
          uu(j, i) = 0.0;
        }
        if (dept(j, i) < wtfh && vv(j, i) > 0.0) {
          vv(j, i) = 0.0;
        }
        if (dept(j + 1, i) < wtfh && vv(j, i) < 0.0) {
          vv(j, i) = 0.0;
        }
        // Faces of wetting-limited cells (shallowwater.c:810-817).
        if (cflActive(j, i) == 1.0 || cflActive(j, i + 1) == 1.0) {
          uu(j, i) = 0.0;
        }
        if (cflActive(j, i) == 1.0 || cflActive(j + 1, i) == 1.0) {
          vv(j, i) = 0.0;
        }
      });

  // Flow rates and CFL numbers from the limited velocities
  // (shallowwater.c:820-826); the wetting flags are consumed and reset.
  Kokkos::parallel_for(
      "swe_flow_rates",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        Fu(j, i) = uu(j, i) * Asx(j, i);
        Fv(j, i) = vv(j, i) * Asy(j, i);
        cflx(j, i) = Kokkos::fabs(uu(j, i) * dt / dx);
        cfly(j, i) = Kokkos::fabs(vv(j, i) * dt / dy);
        cflActive(j, i) = 0.0;
      });
}

void SurfaceSolver::enforceVeloBc(VeloBcApply apply) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const bool westEdge = (grid_.rankWest() == MPI_PROC_NULL);
  const bool eastEdge = (grid_.rankEast() == MPI_PROC_NULL);
  const bool southEdge = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool northEdge = (grid_.rankNorth() == MPI_PROC_NULL);
  Field2<real_t> uu = uu_, vv = vv_;

  // Zero-gradient velocity ghosts along physical edges
  // (shallowwater.c:896-920).
  Kokkos::parallel_for(
      "swe_velo_edge_ghosts_i",
      Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nyl + 1),
      KOKKOS_LAMBDA(const int j) {
        if (westEdge) {
          uu(j, 0) = uu(j, 1);
          vv(j, 0) = vv(j, 1);
        }
        if (eastEdge) {
          uu(j, nxl + 1) = uu(j, nxl);
          vv(j, nxl + 1) = vv(j, nxl);
        }
      });
  Kokkos::parallel_for(
      "swe_velo_edge_ghosts_j",
      Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nxl + 1),
      KOKKOS_LAMBDA(const int i) {
        if (southEdge) {
          uu(0, i) = uu(1, i);
          vv(0, i) = vv(1, i);
        }
        if (northEdge) {
          uu(nyl + 1, i) = uu(nyl, i);
          vv(nyl + 1, i) = vv(nyl, i);
        }
      });

  // Prescribed-face velocities of velocity conditions, applied after the
  // limiters like the legacy stage-boundary correction.
  for (const DeviceBcList& list : velocityBcs_) {
    auto lj = list.j;
    auto li = list.i;
    auto lface = list.face;
    const real_t u = list.current;
    Kokkos::parallel_for(
        "swe_velocity_bc_faces",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m) {
          const int j = lj(m);
          const int i = li(m);
          switch (static_cast<BcFace>(lface(m))) {
            case BcFace::XMinus:
              uu(j, 0) = u;
              break;
            case BcFace::XPlus:
              uu(j, i) = u;
              break;
            case BcFace::YMinus:
              vv(0, i) = u;
              break;
            case BcFace::YPlus:
              vv(j, i) = u;
              break;
            default:
              break;
          }
        });
  }

  // Mass-consistent boundary-face velocity at prescribed-stage cells on the
  // global domain edges (shallowwater.c:922-958). Legacy's north-edge rank
  // test used ">" where ">=" was meant, which silently skipped the first
  // rank of the top row — including single-rank runs; the corrected global
  // test is recorded in docs/theory/surface-water.md. In Refresh mode only
  // the west/south edge-slot writes run (see VeloBcApply); the inputs are
  // the reconstructed flow rates, the restored eta^n, and the completed
  // step's dt, so the slots restore bitwise.
  if (!etaBcs_.empty()) {
    const bool edgeOnly = (apply == VeloBcApply::Refresh);
    halo_.exchange({"swe_Fu", "swe_Fv"});
    const real_t dt = dt_;
    const int nxGlobal = grid_.nx();
    const int nyGlobal = grid_.ny();
    const int i0 = grid_.i0();
    const int j0 = grid_.j0();
    Field2<real_t> eta = eta_, etan = etan_, Asz = Asz_, Asx = Asx_, Asy = Asy_;
    Field2<real_t> Fu = Fu_, Fv = Fv_;
    for (const DeviceBcList& list : etaBcs_) {
      auto lj = list.j;
      auto li = list.i;
      Kokkos::parallel_for(
          "swe_eta_bc_velocity",
          Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
          KOKKOS_LAMBDA(const std::size_t m) {
            const int j = lj(m);
            const int i = li(m);
            const int gi = i0 + i - 1;
            const int gj = j0 + j - 1;
            const real_t Fw = (eta(j, i) - etan(j, i)) * Asz(j, i) / dt;
            if (gi == 0) {
              const real_t Ftot = -Fu(j, i) - Fv(j, i) + Fv(j - 1, i) - Fw;
              uu(j, 0) = (Asx(j, i) > 0.0) ? -Ftot / Asx(j, i) : 0.0;
            }
            if (gi == nxGlobal - 1 && !edgeOnly) {
              const real_t Ftot = Fu(j, i - 1) - Fv(j, i) + Fv(j - 1, i) - Fw;
              uu(j, i) = (Asx(j, i) > 0.0) ? Ftot / Asx(j, i) : 0.0;
            }
            if (gj == 0) {
              const real_t Ftot = -Fv(j, i) - Fu(j, i) + Fu(j, i - 1) - Fw;
              vv(0, i) = (Asy(j, i) > 0.0) ? -Ftot / Asy(j, i) : 0.0;
            }
            if (gj == nyGlobal - 1 && !edgeOnly) {
              const real_t Ftot = Fv(j - 1, i) - Fu(j, i) + Fu(j, i - 1) - Fw;
              vv(j, i) = (Asy(j, i) > 0.0) ? Ftot / Asy(j, i) : 0.0;
            }
          });
    }
  }
}

void SurfaceSolver::fillVelocityGhostCorners() {
  // After the corner exchange, the ghost-row halo columns along physical
  // j-edges are the only unfilled diagonal-stencil cells: extend the
  // zero-gradient copy across the full row width. Columns settled by the
  // exchange (or by the physical-edge fills) are already consistent.
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const bool southEdge = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool northEdge = (grid_.rankNorth() == MPI_PROC_NULL);
  if (!southEdge && !northEdge) {
    return;
  }
  Field2<real_t> uu = uu_, vv = vv_;
  Kokkos::parallel_for(
      "swe_velo_ghost_corners", Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(0, 2),
      KOKKOS_LAMBDA(const int side) {
        const int ii = (side == 0) ? 0 : nxl + 1;
        if (southEdge) {
          uu(0, ii) = uu(1, ii);
          vv(0, ii) = vv(1, ii);
        }
        if (northEdge) {
          uu(nyl + 1, ii) = uu(nyl, ii);
          vv(nyl + 1, ii) = vv(nyl, ii);
        }
      });
}

}  // namespace frehg::swe
