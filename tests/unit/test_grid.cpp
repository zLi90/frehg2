/// \file test_grid.cpp
/// \brief Grid tests: non-divisible block math (NX = 101 over 4 blocks),
///        gid compression with masked columns, vertical stretch (plan §8.1).

#include "core/Grid.hpp"
#include "core/Types.hpp"

#include <gtest/gtest.h>
#include <mpi.h>

#include <set>

namespace {

using frehg::DomainConfig;
using frehg::Grid;

DomainConfig smallDomain() {
  DomainConfig dom;
  dom.nx = 5;
  dom.ny = 4;
  dom.nz = 3;
  dom.dx = 2.0;
  dom.dy = 0.5;
  dom.dz = 0.1;
  dom.dzStretch = 1.0;
  return dom;
}

TEST(GridBlocks, NonDivisible101Over4) {
  // b5's extent: NX = 101 split over 4 ranks (plan §5.2).
  EXPECT_EQ(Grid::blockSize(101, 4, 0), 26);
  EXPECT_EQ(Grid::blockSize(101, 4, 1), 25);
  EXPECT_EQ(Grid::blockSize(101, 4, 2), 25);
  EXPECT_EQ(Grid::blockSize(101, 4, 3), 25);
  EXPECT_EQ(Grid::blockStart(101, 4, 0), 0);
  EXPECT_EQ(Grid::blockStart(101, 4, 1), 26);
  EXPECT_EQ(Grid::blockStart(101, 4, 2), 51);
  EXPECT_EQ(Grid::blockStart(101, 4, 3), 76);
}

TEST(GridBlocks, PartitionCoversExactly) {
  for (const int N : {1, 7, 55, 101, 200}) {
    for (const int P : {1, 2, 3, 4, 7}) {
      int total = 0;
      for (int p = 0; p < P; ++p) {
        EXPECT_EQ(Grid::blockStart(N, P, p), total) << "N=" << N << " P=" << P << " p=" << p;
        total += Grid::blockSize(N, P, p);
      }
      EXPECT_EQ(total, N) << "N=" << N << " P=" << P;
    }
  }
}

TEST(GridSerial, ExtentsAndCoordinates) {
  const Grid grid(MPI_COMM_SELF, smallDomain());
  EXPECT_EQ(grid.nx(), 5);
  EXPECT_EQ(grid.nyLocal(), 4);
  EXPECT_EQ(grid.nxLocal(), 5);
  EXPECT_EQ(grid.i0(), 0);
  EXPECT_EQ(grid.j0(), 0);
  EXPECT_EQ(grid.px(), 1);
  EXPECT_EQ(grid.py(), 1);
  EXPECT_EQ(grid.rankWest(), MPI_PROC_NULL);
  EXPECT_EQ(grid.rankEast(), MPI_PROC_NULL);
  EXPECT_DOUBLE_EQ(grid.xCenter(0), 1.0);
  EXPECT_DOUBLE_EQ(grid.xCenter(4), 9.0);
  EXPECT_DOUBLE_EQ(grid.yCenter(3), 1.75);
}

TEST(GridSerial, VerticalStretch) {
  DomainConfig dom = smallDomain();
  dom.dzStretch = 2.0;
  const Grid grid(MPI_COMM_SELF, dom);
  EXPECT_DOUBLE_EQ(grid.dzK(0), 0.1);
  EXPECT_DOUBLE_EQ(grid.dzK(1), 0.2);
  EXPECT_DOUBLE_EQ(grid.dzK(2), 0.4);
  EXPECT_DOUBLE_EQ(grid.totalDepth(), 0.7);
  EXPECT_DOUBLE_EQ(grid.zCenter(0), -0.05);
  EXPECT_DOUBLE_EQ(grid.zCenter(1), -0.2);
  EXPECT_DOUBLE_EQ(grid.zCenter(2), -0.5);
}

TEST(GridSerial, DefaultAllActive) {
  const Grid grid(MPI_COMM_SELF, smallDomain());
  EXPECT_EQ(grid.activeCount2Local(), 20);
  EXPECT_EQ(grid.activeCount2Global(), 20);
  EXPECT_EQ(grid.activeCount3Global(), 60);
  EXPECT_EQ(grid.offset2(), 0);
  EXPECT_EQ(grid.offset3(), 0);
}

TEST(GridSerial, GidCompressionWithMaskedColumns) {
  Grid grid(MPI_COMM_SELF, smallDomain());
  // ktop: column (0,0) fully inactive (nz), column (1,2) starts at k=2,
  // column (2,3) starts at k=1, everything else fully active.
  frehg::HostField2<int> ktop("ktop", 4, 5);
  Kokkos::deep_copy(ktop, 0);
  ktop(0, 0) = 3;
  ktop(1, 2) = 2;
  ktop(2, 3) = 1;
  grid.buildGlobalIds(ktop);

  // 2D: 19 active surface cells (one fully masked column).
  EXPECT_EQ(grid.activeCount2Local(), 19);
  // 3D: 20 columns * 3 layers - 3 (full col) - 2 - 1 = 54.
  EXPECT_EQ(grid.activeCount3Local(), 54);

  const auto& gid2 = grid.gid2Host();
  const auto& gid3 = grid.gid3Host();

  // Masked column: -1 everywhere.
  EXPECT_EQ(gid2(1, 1), -1);
  EXPECT_EQ(gid3(1, 1, 0), -1);
  EXPECT_EQ(gid3(1, 1, 2), -1);

  // Partially masked columns: -1 above ktop, ids below.
  EXPECT_EQ(gid3(2, 3, 0), -1);
  EXPECT_EQ(gid3(2, 3, 1), -1);
  EXPECT_GE(gid3(2, 3, 2), 0);
  EXPECT_EQ(gid3(3, 4, 0), -1);
  EXPECT_GE(gid3(3, 4, 1), 0);

  // Ids are a compressed permutation: consecutive, unique, covering
  // [0, activeCount).
  std::set<PetscInt> seen2;
  std::set<PetscInt> seen3;
  for (std::size_t j = 1; j <= 4; ++j) {
    for (std::size_t i = 1; i <= 5; ++i) {
      if (gid2(j, i) >= 0) {
        seen2.insert(gid2(j, i));
      }
      for (std::size_t k = 0; k < 3; ++k) {
        if (gid3(j, i, k) >= 0) {
          seen3.insert(gid3(j, i, k));
        }
      }
    }
  }
  EXPECT_EQ(seen2.size(), 19u);
  EXPECT_EQ(*seen2.begin(), 0);
  EXPECT_EQ(*seen2.rbegin(), 18);
  EXPECT_EQ(seen3.size(), 54u);
  EXPECT_EQ(*seen3.begin(), 0);
  EXPECT_EQ(*seen3.rbegin(), 53);

  // k is innermost: within a column, ids are consecutive.
  EXPECT_EQ(gid3(1, 2, 1), gid3(1, 2, 0) + 1);
  EXPECT_EQ(gid3(1, 2, 2), gid3(1, 2, 1) + 1);

  // Domain-edge halos carry no ids.
  EXPECT_EQ(gid2(0, 1), -1);
  EXPECT_EQ(gid2(1, 0), -1);
  EXPECT_EQ(gid3(0, 1, 0), -1);
}

TEST(GridSerial, KtopValidation) {
  Grid grid(MPI_COMM_SELF, smallDomain());
  frehg::HostField2<int> wrongShape("ktop", 3, 5);
  EXPECT_THROW(grid.buildGlobalIds(wrongShape), frehg::FatalError);

  frehg::HostField2<int> outOfRange("ktop", 4, 5);
  Kokkos::deep_copy(outOfRange, 0);
  outOfRange(0, 0) = 4;  // > nz
  EXPECT_THROW(grid.buildGlobalIds(outOfRange), frehg::FatalError);
}

TEST(GridSerial, ExplicitDecompositionMismatchIsFatal) {
  DomainConfig dom = smallDomain();
  dom.decomposition.mpiNx = 2;
  dom.decomposition.mpiNy = 2;  // needs 4 ranks, COMM_SELF has 1
  EXPECT_THROW(Grid(MPI_COMM_SELF, dom), frehg::FatalError);
}

}  // namespace
