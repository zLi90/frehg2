/// \file test_parallel_core.cpp
/// \brief Rank-count-portable checks of the parallel foundation (plan §8.2):
///        gid compression across ranks (incl. the NX = 101 / 4-rank case),
///        polygon rasterization spanning rank boundaries with sub-
///        communicators, parallel HDF5 layout, monitors, and checkpoint
///        round-trip. Runs identically at 1, 2, and 4 ranks.

#include "bc/BoundarySet.hpp"
#include "core/Grid.hpp"
#include "core/Logger.hpp"
#include "core/PetscSession.hpp"
#include "core/Types.hpp"
#include "io/Checkpoint.hpp"
#include "io/Hdf5Output.hpp"
#include "io/Monitor.hpp"

#include <mpi.h>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

namespace {

using frehg::real_t;

int failures = 0;

void check(bool ok, const std::string& what) {
  if (!ok) {
    ++failures;
    frehg::log::error("FAIL: " + what);
  }
}

bool bitEqual(real_t a, real_t b) { return std::memcmp(&a, &b, sizeof(real_t)) == 0; }

/// gid coverage: gathered owned gids must be exactly {0 .. N-1}, and halo
/// gids must match what the neighbor rank assigned.
void testGridGids(const frehg::PetscSession& session) {
  frehg::DomainConfig dom;
  dom.nx = 101;  // b5's non-divisible extent (plan §5.2)
  dom.ny = 6;
  dom.nz = 3;
  dom.dx = 1.0;
  dom.dy = 1.0;
  dom.dz = 0.2;
  // Force the whole rank count onto i so the non-divisible 101-block math is
  // exercised (26/25/25/25 at 4 ranks).
  dom.decomposition.mpiNx = session.size();
  dom.decomposition.mpiNy = 1;
  frehg::Grid grid(session.comm(), dom);

  // Mask: a diagonal band of partially inactive columns, one fully dead
  // column per rank block.
  frehg::HostField2<int> ktop("ktop", static_cast<std::size_t>(grid.nyLocal()),
                              static_cast<std::size_t>(grid.nxLocal()));
  for (int j = 0; j < grid.nyLocal(); ++j) {
    for (int i = 0; i < grid.nxLocal(); ++i) {
      const int iGlob = grid.i0() + i;
      const int jGlob = grid.j0() + j;
      int kt = (iGlob + jGlob) % 3 == 0 ? 1 : 0;
      if (iGlob % 17 == 3 && jGlob == 2) {
        kt = grid.nz();  // fully inactive
      }
      ktop(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) = kt;
    }
  }
  grid.buildGlobalIds(ktop);

  // Collect owned 3D gids and verify global coverage 0 .. N-1.
  std::vector<PetscInt> mine;
  const auto& gid3 = grid.gid3Host();
  for (int j = 1; j <= grid.nyLocal(); ++j) {
    for (int i = 1; i <= grid.nxLocal(); ++i) {
      for (int k = 0; k < grid.nz(); ++k) {
        const PetscInt g = gid3(static_cast<std::size_t>(j), static_cast<std::size_t>(i),
                                static_cast<std::size_t>(k));
        if (g >= 0) {
          mine.push_back(g);
        }
      }
    }
  }
  check(static_cast<PetscInt>(mine.size()) == grid.activeCount3Local(),
        "gid3 owned count matches activeCount3Local");

  const int localCount = static_cast<int>(mine.size());
  std::vector<int> counts(static_cast<std::size_t>(session.size()), 0);
  MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, session.comm());
  std::vector<int> displs(counts.size(), 0);
  int total = 0;
  for (std::size_t r = 0; r < counts.size(); ++r) {
    displs[r] = total;
    total += counts[r];
  }
  std::vector<PetscInt> all(static_cast<std::size_t>(total));
  const MPI_Datatype dtype = sizeof(PetscInt) == 8 ? MPI_LONG_LONG : MPI_INT;
  MPI_Allgatherv(mine.data(), localCount, dtype, all.data(), counts.data(), displs.data(), dtype,
                 session.comm());
  check(static_cast<PetscInt>(all.size()) == grid.activeCount3Global(),
        "global gid3 count matches activeCount3Global");
  std::vector<char> seen(all.size(), 0);
  bool coverage = true;
  for (const PetscInt g : all) {
    if (g < 0 || g >= static_cast<PetscInt>(all.size()) ||
        seen[static_cast<std::size_t>(g)] != 0) {
      coverage = false;
      break;
    }
    seen[static_cast<std::size_t>(g)] = 1;
  }
  check(coverage, "gid3 ids are a permutation of 0..N-1 across ranks");

  // Halo gid values match the neighbor's analytic re-computation: rebuild
  // the expected gid for a global column by replaying the deterministic
  // mask, on the rank that owns the halo cell's column block.
  // (Cheap cross-check: halo gid2 of my east halo equals the first owned
  // gid2 of the east neighbor's block column - verified via sendrecv.)
  const auto& gid2 = grid.gid2Host();
  PetscInt sendWest = gid2(1, 1);
  PetscInt recvEast = -12345;
  MPI_Sendrecv(&sendWest, 1, dtype, grid.rankWest(), 7, &recvEast, 1, dtype, grid.rankEast(), 7,
               session.comm(), MPI_STATUS_IGNORE);
  if (grid.rankEast() != MPI_PROC_NULL) {
    check(gid2(1, static_cast<std::size_t>(grid.nxLocal()) + 1) == recvEast,
          "east halo gid2 equals neighbor's owned gid2");
  }

  // 4-rank b5 case from the plan: block sizes 26/25/25/25.
  if (session.size() == 4 && grid.px() == 4) {
    check(grid.nxLocal() == (grid.pi() == 0 ? 26 : 25), "NX=101 four-rank block sizes");
  }
}

/// b5's outlet polygon rasterized across rank boundaries: the global count
/// must be 5 regardless of rank count, and the sub-communicator spans
/// exactly the owning ranks.
/// The automatic decomposition must respect the grid extents: 4 ranks on a
/// 1 x 10 domain (b1's shape) must choose 1 x 4, not MPI_Dims_create's
/// 2 x 2, which would leave empty blocks.
void testNarrowDomainDecomposition(const frehg::PetscSession& session) {
  frehg::DomainConfig dom;
  dom.nx = 1;
  dom.ny = 10;
  dom.nz = 1;
  dom.dx = 1.0;
  dom.dy = 1.0;
  dom.dz = 1.0;
  const frehg::Grid grid(session.comm(), dom);
  check(grid.px() == 1, "narrow domain keeps a single rank column");
  check(grid.py() == session.size(), "narrow domain splits all ranks along j");
  check(grid.nyLocal() >= 1 && grid.nxLocal() == 1, "narrow domain block is non-empty");
}

void testBoundarySpanningRanks(const frehg::PetscSession& session) {  frehg::DomainConfig dom;
  dom.nx = 101;
  dom.ny = 55;
  dom.nz = 2;
  dom.dx = 1.0;
  dom.dy = 1.0;
  dom.dz = 0.2;
  const frehg::Grid grid(session.comm(), dom);

  frehg::BoundaryConditionConfig outlet;
  outlet.name = "outlet";
  outlet.polygon = {{-0.1, 50.1}, {0.9, 50.1}, {0.9, 55.1}, {-0.1, 55.1}};
  outlet.target = frehg::BcTarget::Surface;
  outlet.kind = frehg::BcKind::Eta;
  outlet.value.form = frehg::BcValueConfig::Form::Constant;
  outlet.value.constant = 1.5;

  // A band polygon crossing every rank boundary in i.
  frehg::BoundaryConditionConfig band;
  band.name = "band";
  band.polygon = {{-0.5, 9.9}, {101.5, 9.9}, {101.5, 12.1}, {-0.5, 12.1}};
  band.target = frehg::BcTarget::Surface;
  band.kind = frehg::BcKind::Discharge;
  band.value.form = frehg::BcValueConfig::Form::Constant;
  band.value.constant = 2.0;

  const frehg::BoundarySet set(grid, {outlet, band}, ".");

  check(set.byName("outlet").globalCellCount() == 5, "outlet selects 5 cells at any rank count");
  // Band: full x extent (101 cells) over global rows j = 10, 11.
  check(set.byName("band").globalCellCount() == 202, "band selects 202 cells at any rank count");

  // Sub-communicator size equals the number of ranks owning members.
  const frehg::BoundaryCondition& bandBc = set.byName("band");
  const int owns = bandBc.active() ? 1 : 0;
  int owningRanks = 0;
  MPI_Allreduce(&owns, &owningRanks, 1, MPI_INT, MPI_SUM, session.comm());
  if (bandBc.active()) {
    int subSize = 0;
    MPI_Comm_size(bandBc.subComm(), &subSize);
    check(subSize == owningRanks, "band sub-communicator spans exactly the owning ranks");
    // The band crosses every i-block: all ranks whose j-range covers 10-11.
    long localMembers = static_cast<long>(bandBc.cells().size());
    long subTotal = 0;
    MPI_Allreduce(&localMembers, &subTotal, 1, MPI_LONG, MPI_SUM, bandBc.subComm());
    check(subTotal == 202, "band member total over the sub-communicator");
  }
}

/// Parallel HDF5: a globally-defined analytic field written by N ranks must
/// read back bit-exact in the §7 flattening on every rank.
void testParallelHdf5(const frehg::PetscSession& session, const std::string& scratchDir) {
  frehg::DomainConfig dom;
  dom.nx = 9;  // non-divisible over 2 and 4
  dom.ny = 5;
  dom.nz = 3;
  dom.dx = 1.0;
  dom.dy = 1.0;
  dom.dz = 0.5;
  const frehg::Grid grid(session.comm(), dom);

  const std::size_t nyl = static_cast<std::size_t>(grid.nyLocal());
  const std::size_t nxl = static_cast<std::size_t>(grid.nxLocal());
  const std::size_t nzg = static_cast<std::size_t>(grid.nz());

  auto value2 = [](int iGlob, int jGlob) -> real_t {
    return (1.0 / 3.0) + 7.0 * iGlob + 1000.0 * jGlob;
  };
  auto value3 = [&value2](int iGlob, int jGlob, int k) -> real_t {
    return value2(iGlob, jGlob) + (1.0 / 7.0) * (k + 1);
  };

  frehg::Field2<real_t> eta("eta", nyl + 2, nxl + 2);
  frehg::Field3<real_t> head("head", nyl + 2, nxl + 2, nzg);
  auto etaH = Kokkos::create_mirror_view(eta);
  auto headH = Kokkos::create_mirror_view(head);
  for (std::size_t j = 1; j <= nyl; ++j) {
    for (std::size_t i = 1; i <= nxl; ++i) {
      const int iGlob = grid.i0() + static_cast<int>(i) - 1;
      const int jGlob = grid.j0() + static_cast<int>(j) - 1;
      etaH(j, i) = value2(iGlob, jGlob);
      for (std::size_t k = 0; k < nzg; ++k) {
        headH(j, i, k) = value3(iGlob, jGlob, static_cast<int>(k));
      }
    }
  }
  Kokkos::deep_copy(eta, etaH);
  Kokkos::deep_copy(head, headH);

  const std::string filename = scratchDir + "/parallel_core_n" +
                               std::to_string(session.size()) + ".h5";
  {
    frehg::io::Hdf5Output out(grid, filename, "mpi-test-config");
    out.writeField2("surface", "eta", 60.0, eta, {"m", "eta"});
    out.writeField3("groundwater", "hydraulic_head", 60.0, head, {"m", "head"});

    const std::vector<real_t> flat2 = out.readAll1D("/surface/eta/60");
    check(flat2.size() == static_cast<std::size_t>(dom.nx * dom.ny), "eta dataset size");
    bool ok2 = true;
    for (int j = 0; j < dom.ny; ++j) {
      for (int i = 0; i < dom.nx; ++i) {
        ok2 = ok2 && bitEqual(flat2[static_cast<std::size_t>(j * dom.nx + i)], value2(i, j));
      }
    }
    check(ok2, "parallel 2D layout j*NX+i bit-exact");

    const std::vector<real_t> flat3 = out.readAll1D("/groundwater/hydraulic_head/60");
    bool ok3 = true;
    for (int j = 0; j < dom.ny; ++j) {
      for (int i = 0; i < dom.nx; ++i) {
        for (int k = 0; k < dom.nz; ++k) {
          ok3 = ok3 && bitEqual(flat3[static_cast<std::size_t>((j * dom.nx + i) * dom.nz + k)],
                                value3(i, j, k));
        }
      }
    }
    check(ok3, "parallel 3D layout (j*NX+i)*NZ+k bit-exact");

    // Monitor at a fixed global point: exactly one rank owns it.
    frehg::MonitorConfig monCfg;
    monCfg.name = "probe";
    monCfg.i = 4;
    monCfg.j = 2;
    monCfg.variables = {"eta"};
    frehg::io::Monitor monitor(out, monCfg);
    const int owns = monitor.ownsPoint() ? 1 : 0;
    int owners = 0;
    MPI_Allreduce(&owns, &owners, 1, MPI_INT, MPI_SUM, session.comm());
    check(owners == 1, "exactly one rank owns the monitor point");
    monitor.record(0.0, {value2(4, 2)});
    monitor.record(60.0, {value2(4, 2) + 1.0});
    monitor.flush();
    check(monitor.rowsWritten() == 2, "monitor rows written on all ranks agree");

    // Checkpoint round-trip across ranks.
    frehg::io::Checkpoint checkpoint(out);
    checkpoint.write(60.0, 12, {{"dtg", 0.75}}, {{"eta", eta}}, {{"head", head}});
    Kokkos::deep_copy(eta, 0.0);
    Kokkos::deep_copy(head, 0.0);
    const auto header = checkpoint.read(60.0, {{"eta", eta}}, {{"head", head}});
    check(header.step == 12 && bitEqual(header.scalars.at("dtg"), 0.75), "checkpoint header");
    Kokkos::deep_copy(etaH, eta);
    Kokkos::deep_copy(headH, head);
    bool okCkpt = true;
    for (std::size_t j = 1; j <= nyl; ++j) {
      for (std::size_t i = 1; i <= nxl; ++i) {
        const int iGlob = grid.i0() + static_cast<int>(i) - 1;
        const int jGlob = grid.j0() + static_cast<int>(j) - 1;
        okCkpt = okCkpt && bitEqual(etaH(j, i), value2(iGlob, jGlob));
        for (std::size_t k = 0; k < nzg; ++k) {
          okCkpt = okCkpt &&
                   bitEqual(headH(j, i, k), value3(iGlob, jGlob, static_cast<int>(k)));
        }
      }
    }
    check(okCkpt, "checkpoint fields bit-exact after parallel round-trip");
  }
  MPI_Barrier(session.comm());
  if (session.rank() == 0) {
    std::filesystem::remove(filename);
  }
}

}  // namespace

int main(int argc, char** argv) {
  frehg::PetscSession session(argc, argv);
  {
    std::string scratchDir;
    if (session.rank() == 0) {
      scratchDir = (std::filesystem::temp_directory_path() / "frehg_mpi_core").string();
      std::filesystem::create_directories(scratchDir);
    }
    // Share the scratch dir so all ranks open the same file.
    int len = static_cast<int>(scratchDir.size());
    MPI_Bcast(&len, 1, MPI_INT, 0, session.comm());
    scratchDir.resize(static_cast<std::size_t>(len));
    MPI_Bcast(scratchDir.data(), len, MPI_CHAR, 0, session.comm());

    testGridGids(session);
    testNarrowDomainDecomposition(session);
    testBoundarySpanningRanks(session);
    testParallelHdf5(session, scratchDir);
  }

  int globalFailures = 0;
  MPI_Allreduce(&failures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (session.rank() == 0) {
    frehg::log::info(globalFailures == 0
                         ? "test_parallel_core: PASS"
                         : "test_parallel_core: FAILED (" + std::to_string(globalFailures) +
                               " checks)");
  }
  return globalFailures == 0 ? 0 : 1;
}
