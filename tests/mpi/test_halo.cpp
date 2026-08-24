/// \file test_halo.cpp
/// \brief MPI halo-exchange driver: analytic f(i, j, k) fields must arrive
///        in the ghost cells bitwise-exact at 1/2/4 ranks (plan §8.2), in
///        both GPU-aware and host-staged modes, for full and targeted
///        exchanges. Domain-edge halos must stay untouched.

#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Logger.hpp"
#include "core/PetscSession.hpp"
#include "core/Types.hpp"

#include <mpi.h>

#include <cmath>
#include <cstring>
#include <string>

namespace {

using frehg::real_t;

int failures = 0;

void check(bool ok, const std::string& what) {
  if (!ok) {
    ++failures;
    frehg::log::error("FAIL: " + what);
  }
}

/// Analytic global cell values: integers, exactly representable, so
/// "bitwise-exact" is meaningful.
real_t f2(int iGlob, int jGlob) { return 1000.0 * iGlob + jGlob; }
real_t f3(int iGlob, int jGlob, int k) { return 100000.0 * iGlob + 100.0 * jGlob + k; }

constexpr real_t kSentinel = -777777.0;

bool bitEqual(real_t a, real_t b) { return std::memcmp(&a, &b, sizeof(real_t)) == 0; }

void runCase(const frehg::Grid& grid, bool gpuAware) {
  const std::string tag = gpuAware ? "[gpu-aware] " : "[host-staged] ";
  const std::size_t nyl = static_cast<std::size_t>(grid.nyLocal());
  const std::size_t nxl = static_cast<std::size_t>(grid.nxLocal());
  const std::size_t nzg = static_cast<std::size_t>(grid.nz());

  frehg::Field2<real_t> eta("eta", nyl + 2, nxl + 2);
  frehg::Field2<real_t> depth("depth", nyl + 2, nxl + 2);
  frehg::Field3<real_t> head("head", nyl + 2, nxl + 2, nzg);

  auto fill = [&](real_t offset) {
    auto etaH = Kokkos::create_mirror_view(eta);
    auto depthH = Kokkos::create_mirror_view(depth);
    auto headH = Kokkos::create_mirror_view(head);
    Kokkos::deep_copy(etaH, kSentinel);
    Kokkos::deep_copy(depthH, kSentinel);
    Kokkos::deep_copy(headH, kSentinel);
    for (std::size_t j = 1; j <= nyl; ++j) {
      for (std::size_t i = 1; i <= nxl; ++i) {
        const int iGlob = grid.i0() + static_cast<int>(i) - 1;
        const int jGlob = grid.j0() + static_cast<int>(j) - 1;
        etaH(j, i) = f2(iGlob, jGlob) + offset;
        depthH(j, i) = 2.0 * f2(iGlob, jGlob) + offset;
        for (std::size_t k = 0; k < nzg; ++k) {
          headH(j, i, k) = f3(iGlob, jGlob, static_cast<int>(k)) + offset;
        }
      }
    }
    Kokkos::deep_copy(eta, etaH);
    Kokkos::deep_copy(depth, depthH);
    Kokkos::deep_copy(head, headH);
  };

  frehg::HaloExchanger halo(grid, gpuAware);
  halo.add("eta", eta);
  halo.add("depth", depth);
  halo.add("head", head);
  check(halo.fieldCount() == 3, tag + "field count");
  check(halo.gpuAware() == gpuAware, tag + "gpu-aware flag");

  // ---- exchangeAll: every ghost matches the analytic neighbor value. ----
  fill(0.0);
  halo.exchangeAll();

  auto etaH = Kokkos::create_mirror_view(eta);
  auto depthH = Kokkos::create_mirror_view(depth);
  auto headH = Kokkos::create_mirror_view(head);
  Kokkos::deep_copy(etaH, eta);
  Kokkos::deep_copy(depthH, depth);
  Kokkos::deep_copy(headH, head);

  auto checkGhostColumn = [&](std::size_t ih, int iGlob, bool hasNeighbor,
                              const std::string& side) {
    for (std::size_t j = 1; j <= nyl; ++j) {
      const int jGlob = grid.j0() + static_cast<int>(j) - 1;
      if (hasNeighbor) {
        check(bitEqual(etaH(j, ih), f2(iGlob, jGlob)), tag + side + " eta ghost bitwise");
        check(bitEqual(depthH(j, ih), 2.0 * f2(iGlob, jGlob)),
              tag + side + " depth ghost bitwise");
        for (std::size_t k = 0; k < nzg; ++k) {
          check(bitEqual(headH(j, ih, k), f3(iGlob, jGlob, static_cast<int>(k))),
                tag + side + " head ghost bitwise");
        }
      } else {
        check(bitEqual(etaH(j, ih), kSentinel), tag + side + " domain-edge halo untouched");
      }
    }
  };
  checkGhostColumn(0, grid.i0() - 1, grid.rankWest() != MPI_PROC_NULL, "west");
  checkGhostColumn(nxl + 1, grid.i0() + static_cast<int>(nxl),
                   grid.rankEast() != MPI_PROC_NULL, "east");

  auto checkGhostRow = [&](std::size_t jh, int jGlob, bool hasNeighbor,
                           const std::string& side) {
    for (std::size_t i = 1; i <= nxl; ++i) {
      const int iGlob = grid.i0() + static_cast<int>(i) - 1;
      if (hasNeighbor) {
        check(bitEqual(etaH(jh, i), f2(iGlob, jGlob)), tag + side + " eta ghost bitwise");
        for (std::size_t k = 0; k < nzg; ++k) {
          check(bitEqual(headH(jh, i, k), f3(iGlob, jGlob, static_cast<int>(k))),
                tag + side + " head ghost bitwise");
        }
      } else {
        check(bitEqual(etaH(jh, i), kSentinel), tag + side + " domain-edge halo untouched");
      }
    }
  };
  checkGhostRow(0, grid.j0() - 1, grid.rankSouth() != MPI_PROC_NULL, "south");
  checkGhostRow(nyl + 1, grid.j0() + static_cast<int>(nyl),
                grid.rankNorth() != MPI_PROC_NULL, "north");

  // ---- Targeted exchange: only the named field updates. ----
  fill(0.5);
  halo.exchange({"eta"});
  Kokkos::deep_copy(etaH, eta);
  Kokkos::deep_copy(depthH, depth);
  if (grid.rankWest() != MPI_PROC_NULL) {
    const int iGlob = grid.i0() - 1;
    for (std::size_t j = 1; j <= nyl; ++j) {
      const int jGlob = grid.j0() + static_cast<int>(j) - 1;
      check(bitEqual(etaH(j, 0), f2(iGlob, jGlob) + 0.5), tag + "targeted eta updated");
      check(bitEqual(depthH(j, 0), kSentinel), tag + "targeted depth untouched");
    }
  }

  // ---- Corner exchange: the two-phase mode also fills the four halo
  // corners of 2D fields wherever an interior (or boundary-filled) source
  // exists (plan §5.3 amendment for the uy/vx interpolation stencil). ----
  fill(0.25);
  // Emulate the boundary-condition ghost fills that must precede the
  // exchange: analytic values in the physical-edge halo columns, so the
  // wide rows carry them to the j-neighbors.
  {
    Kokkos::deep_copy(etaH, eta);
    if (grid.rankWest() == MPI_PROC_NULL) {
      for (std::size_t j = 1; j <= nyl; ++j) {
        etaH(j, 0) = f2(grid.i0() - 1, grid.j0() + static_cast<int>(j) - 1) + 0.25;
      }
    }
    if (grid.rankEast() == MPI_PROC_NULL) {
      for (std::size_t j = 1; j <= nyl; ++j) {
        etaH(j, nxl + 1) =
            f2(grid.i0() + static_cast<int>(nxl), grid.j0() + static_cast<int>(j) - 1) + 0.25;
      }
    }
    Kokkos::deep_copy(eta, etaH);
  }
  halo.exchangeWithCorners({"eta"});
  Kokkos::deep_copy(etaH, eta);
  // Edge halos behave exactly like the plain exchange.
  if (grid.rankWest() != MPI_PROC_NULL) {
    for (std::size_t j = 1; j <= nyl; ++j) {
      check(bitEqual(etaH(j, 0), f2(grid.i0() - 1, grid.j0() + static_cast<int>(j) - 1) + 0.25),
            tag + "corner-mode west edge bitwise");
    }
  }
  // Corner halos arrive whenever the j-neighbor exists (its wide row
  // carries the i-halo columns it received in phase 1 or filled at its own
  // physical edge); without a j-neighbor the corners stay untouched.
  auto checkCorner = [&](std::size_t jh, std::size_t ih, int jNbr, int iGlobCorner,
                         int jGlobCorner) {
    if (jNbr == MPI_PROC_NULL) {
      check(bitEqual(etaH(jh, ih), kSentinel), tag + "corner without j-neighbor untouched");
      return;
    }
    check(bitEqual(etaH(jh, ih), f2(iGlobCorner, jGlobCorner) + 0.25),
          tag + "corner ghost bitwise");
  };
  checkCorner(0, 0, grid.rankSouth(), grid.i0() - 1, grid.j0() - 1);
  checkCorner(0, nxl + 1, grid.rankSouth(), grid.i0() + static_cast<int>(nxl), grid.j0() - 1);
  checkCorner(nyl + 1, 0, grid.rankNorth(), grid.i0() - 1, grid.j0() + static_cast<int>(nyl));
  checkCorner(nyl + 1, nxl + 1, grid.rankNorth(), grid.i0() + static_cast<int>(nxl),
              grid.j0() + static_cast<int>(nyl));
}

}  // namespace

int main(int argc, char** argv) {
  frehg::PetscSession session(argc, argv);
  {
    frehg::DomainConfig dom;
    dom.nx = 7;   // deliberately non-divisible over 2 and 4 ranks
    dom.ny = 5;
    dom.nz = 3;
    dom.dx = 1.0;
    dom.dy = 1.0;
    dom.dz = 0.5;
    const frehg::Grid grid(session.comm(), dom);

    runCase(grid, false);
    runCase(grid, true);

    // Unknown-name discipline (single-rank check to keep ranks collective).
    if (session.size() == 1) {
      frehg::HaloExchanger halo(grid, false);
      bool threw = false;
      try {
        halo.exchange({"absent"});
      } catch (const frehg::FatalError&) {
        threw = true;
      }
      check(threw, "unknown field name is fatal");
    }
  }

  int globalFailures = 0;
  MPI_Allreduce(&failures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (session.rank() == 0) {
    frehg::log::info(globalFailures == 0 ? "test_halo: PASS"
                                         : "test_halo: FAILED (" +
                                               std::to_string(globalFailures) + " checks)");
  }
  return globalFailures == 0 ? 0 : 1;
}
