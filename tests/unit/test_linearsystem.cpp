/// \file test_linearsystem.cpp
/// \brief LinearSystem tests: 1D/2D Poisson to 1e-9, symmetry check, COO
///        boundary folding and duplicate summation (plan §8.1).

#include "core/LinearSystem.hpp"
#include "core/Types.hpp"

#include <gtest/gtest.h>
#include <mpi.h>

#include <cmath>
#include <numbers>

namespace {

using frehg::LinearSystem;
using frehg::real_t;
using frehg::SolverSettings;

using IdxView = Kokkos::View<PetscInt*, frehg::MemSpace>;
using ValView = Kokkos::View<real_t*, frehg::MemSpace>;

SolverSettings tightSettings() {
  SolverSettings s;
  s.rtol = 1.0e-12;
  s.atol = 1.0e-30;
  s.maxIterations = 5000;
  return s;
}

/// Assemble the 1D Poisson matrix -u'' with Dirichlet values folded through
/// COO index -1 legs. Returns the infinity-norm error against the exact
/// discrete solution u(x) = 1 + x (which second differences reproduce
/// exactly, so the linear-system error is isolated from discretization
/// error).
real_t solve1DPoisson(bool splitDuplicates, SolverSettings settings = tightSettings(),
                      const std::string& prefix = "t1d_") {
  const int n = 101;
  const real_t h = 1.0 / (n + 1);
  const real_t invH2 = 1.0 / (h * h);
  const real_t uLeft = 1.0;   // u(0)
  const real_t uRight = 2.0;  // u(1)

  // Stencil: 3 legs per row; optionally split each off-diagonal entry into
  // two half-weight duplicates to verify COO duplicate summation.
  const int legsPerRow = splitDuplicates ? 5 : 3;
  const std::size_t nCoo = static_cast<std::size_t>(n * legsPerRow);
  IdxView rows("rows", nCoo);
  IdxView cols("cols", nCoo);
  ValView vals("vals", nCoo);
  auto rowsH = Kokkos::create_mirror_view(rows);
  auto colsH = Kokkos::create_mirror_view(cols);
  auto valsH = Kokkos::create_mirror_view(vals);

  std::size_t e = 0;
  for (int r = 0; r < n; ++r) {
    rowsH(e) = r;
    colsH(e) = r;
    valsH(e) = 2.0 * invH2;
    ++e;
    const PetscInt west = (r > 0) ? r - 1 : -1;      // -1: folded into rhs
    const PetscInt east = (r < n - 1) ? r + 1 : -1;  // -1: folded into rhs
    if (splitDuplicates) {
      rowsH(e) = r; colsH(e) = west; valsH(e) = -0.5 * invH2; ++e;
      rowsH(e) = r; colsH(e) = west; valsH(e) = -0.5 * invH2; ++e;
      rowsH(e) = r; colsH(e) = east; valsH(e) = -0.5 * invH2; ++e;
      rowsH(e) = r; colsH(e) = east; valsH(e) = -0.5 * invH2; ++e;
    } else {
      rowsH(e) = r; colsH(e) = west; valsH(e) = -invH2; ++e;
      rowsH(e) = r; colsH(e) = east; valsH(e) = -invH2; ++e;
    }
  }
  Kokkos::deep_copy(rows, rowsH);
  Kokkos::deep_copy(cols, colsH);
  Kokkos::deep_copy(vals, valsH);

  LinearSystem system(MPI_COMM_SELF, prefix, n, n, settings);
  system.setPattern(rows, cols);
  system.setValues(vals);
  EXPECT_TRUE(system.isSymmetric(1.0e-12));

  ValView rhs("rhs", static_cast<std::size_t>(n));
  ValView sol("sol", static_cast<std::size_t>(n));
  auto rhsH = Kokkos::create_mirror_view(rhs);
  for (int r = 0; r < n; ++r) {
    real_t b = 0.0;  // -u'' = 0 for the linear exact solution
    if (r == 0) {
      b += uLeft * invH2;  // folded Dirichlet contribution
    }
    if (r == n - 1) {
      b += uRight * invH2;
    }
    rhsH(static_cast<std::size_t>(r)) = b;
  }
  Kokkos::deep_copy(rhs, rhsH);

  const frehg::SolveStats stats = system.solve(rhs, sol);
  EXPECT_GT(stats.iterations, 0);
  EXPECT_FALSE(stats.reason.empty());

  auto solH = Kokkos::create_mirror_view(sol);
  Kokkos::deep_copy(solH, sol);
  real_t errInf = 0.0;
  for (int r = 0; r < n; ++r) {
    const real_t x = (r + 1) * h;
    const real_t exact = 1.0 + x;
    errInf = std::max(errInf, std::abs(solH(static_cast<std::size_t>(r)) - exact));
  }
  return errInf;
}

TEST(LinearSystem, Poisson1DWithBoundaryFold) {
  EXPECT_LT(solve1DPoisson(false), 1.0e-9);
}

TEST(LinearSystem, CooDuplicateEntriesAreSummed) {
  EXPECT_LT(solve1DPoisson(true), 1.0e-9);
}

// v2 Q1 (plan §2.2): the same Poisson solve must land at the same accuracy
// under each selectable preconditioner — the unit-level face of the g1
// solver-invariance gate.
TEST(LinearSystem, PoissonSolvesUnderGamg) {
  SolverSettings s = tightSettings();
  s.preconditioner = "gamg";
  EXPECT_LT(solve1DPoisson(false, s, "tgamg_"), 1.0e-9);
}

TEST(LinearSystem, PoissonSolvesUnderHypreAmg) {
  PetscBool hasHypre = PETSC_FALSE;
  ASSERT_EQ(PetscHasExternalPackage("hypre", &hasHypre), PETSC_SUCCESS);
  if (hasHypre == PETSC_FALSE) {
    GTEST_SKIP() << "PETSc built without hypre (configure with --download-hypre)";
  }
  SolverSettings s = tightSettings();
  s.preconditioner = "amg";
  EXPECT_LT(solve1DPoisson(false, s, "tamg_"), 1.0e-9);
}

TEST(LinearSystem, UnknownPreconditionerIsFatal) {
  SolverSettings s;
  s.preconditioner = "ilu";
  EXPECT_THROW(LinearSystem(MPI_COMM_SELF, "tbadpc_", 4, 4, s), frehg::FatalError);
}

// v2 Q1 (plan §2.2.4): the hierarchy-reuse policy rebuilds on the cadence,
// counts every rebuild, and honors forceRebuild().
TEST(LinearSystem, HierarchyReuseRebuildPolicy) {
  const int n = 50;
  const std::size_t nCoo = static_cast<std::size_t>(n) * 3;
  IdxView rows("rows", nCoo);
  IdxView cols("cols", nCoo);
  ValView vals("vals", nCoo);
  auto rowsH = Kokkos::create_mirror_view(rows);
  auto colsH = Kokkos::create_mirror_view(cols);
  auto valsH = Kokkos::create_mirror_view(vals);
  std::size_t e = 0;
  for (int r = 0; r < n; ++r) {
    rowsH(e) = r; colsH(e) = r; valsH(e) = 2.0; ++e;
    rowsH(e) = r; colsH(e) = (r > 0) ? r - 1 : -1; valsH(e) = -1.0; ++e;
    rowsH(e) = r; colsH(e) = (r < n - 1) ? r + 1 : -1; valsH(e) = -1.0; ++e;
  }
  Kokkos::deep_copy(rows, rowsH);
  Kokkos::deep_copy(cols, colsH);
  Kokkos::deep_copy(vals, valsH);

  SolverSettings s = tightSettings();
  s.preconditioner = "gamg";
  s.reuseMaxSolves = 3;
  s.reuseIterationFactor = 100.0;  // cadence is the only live trigger here
  LinearSystem system(MPI_COMM_SELF, "treuse_", n, n, s);
  EXPECT_TRUE(system.reusesHierarchy());
  system.setPattern(rows, cols);

  ValView rhs("rhs", static_cast<std::size_t>(n));
  ValView sol("sol", static_cast<std::size_t>(n));
  Kokkos::deep_copy(rhs, 1.0);
  // Solves 1..7 with the cadence at 3: rebuilds at solves 1, 4, and 7.
  for (int k = 0; k < 7; ++k) {
    system.setValues(vals);  // marks the operator changed, as the model does
    system.solve(rhs, sol);
  }
  EXPECT_EQ(system.telemetry().solves, 7);
  EXPECT_EQ(system.telemetry().rebuilds, 3);
  EXPECT_EQ(system.telemetry().retries, 0);

  system.forceRebuild();  // the wet/dry-mask trigger path
  system.setValues(vals);
  system.solve(rhs, sol);
  EXPECT_EQ(system.telemetry().solves, 8);
  EXPECT_EQ(system.telemetry().rebuilds, 4);
  EXPECT_GT(system.telemetry().totalIterations, 0);
}

// Without reuse (the bjacobi-icc default) every solve refreshes the
// preconditioner — the golden-pinned v1 behavior.
TEST(LinearSystem, DefaultPreconditionerRebuildsEverySolve) {
  SolverSettings s = tightSettings();
  ASSERT_EQ(s.preconditioner, "bjacobi-icc");
  const int n = 8;
  IdxView rows("rows", static_cast<std::size_t>(n));
  IdxView cols("cols", static_cast<std::size_t>(n));
  ValView vals("vals", static_cast<std::size_t>(n));
  auto rowsH = Kokkos::create_mirror_view(rows);
  auto colsH = Kokkos::create_mirror_view(cols);
  for (int r = 0; r < n; ++r) {
    rowsH(static_cast<std::size_t>(r)) = r;
    colsH(static_cast<std::size_t>(r)) = r;
  }
  Kokkos::deep_copy(rows, rowsH);
  Kokkos::deep_copy(cols, colsH);
  Kokkos::deep_copy(vals, 3.0);

  LinearSystem system(MPI_COMM_SELF, "tnoreuse_", n, n, s);
  EXPECT_FALSE(system.reusesHierarchy());
  system.setPattern(rows, cols);
  ValView rhs("rhs", static_cast<std::size_t>(n));
  ValView sol("sol", static_cast<std::size_t>(n));
  Kokkos::deep_copy(rhs, 1.0);
  for (int k = 0; k < 3; ++k) {
    system.setValues(vals);
    system.solve(rhs, sol);
  }
  EXPECT_EQ(system.telemetry().solves, 3);
  EXPECT_EQ(system.telemetry().rebuilds, 3);
}

TEST(LinearSystem, Poisson2DManufacturedSolution) {
  // 2D 5-point Poisson on a 24 x 24 unit-square interior grid; the target
  // x* is imposed by manufacturing b = A x*, isolating solver accuracy from
  // discretization error.
  const int m = 24;
  const int n = m * m;
  const real_t h = 1.0 / (m + 1);
  const real_t invH2 = 1.0 / (h * h);

  const std::size_t nCoo = static_cast<std::size_t>(n) * 5;
  IdxView rows("rows", nCoo);
  IdxView cols("cols", nCoo);
  ValView vals("vals", nCoo);
  auto rowsH = Kokkos::create_mirror_view(rows);
  auto colsH = Kokkos::create_mirror_view(cols);
  auto valsH = Kokkos::create_mirror_view(vals);

  // m is a constant expression, so the lambda reads it without capturing.
  auto gid = [](int ix, int iy) -> PetscInt {
    if (ix < 0 || ix >= m || iy < 0 || iy >= m) {
      return -1;
    }
    return iy * m + ix;
  };

  std::size_t e = 0;
  for (int iy = 0; iy < m; ++iy) {
    for (int ix = 0; ix < m; ++ix) {
      const PetscInt r = gid(ix, iy);
      rowsH(e) = r; colsH(e) = r;               valsH(e) = 4.0 * invH2; ++e;
      rowsH(e) = r; colsH(e) = gid(ix - 1, iy); valsH(e) = -invH2; ++e;
      rowsH(e) = r; colsH(e) = gid(ix + 1, iy); valsH(e) = -invH2; ++e;
      rowsH(e) = r; colsH(e) = gid(ix, iy - 1); valsH(e) = -invH2; ++e;
      rowsH(e) = r; colsH(e) = gid(ix, iy + 1); valsH(e) = -invH2; ++e;
    }
  }
  Kokkos::deep_copy(rows, rowsH);
  Kokkos::deep_copy(cols, colsH);
  Kokkos::deep_copy(vals, valsH);

  LinearSystem system(MPI_COMM_SELF, "t2d_", n, n, tightSettings());
  system.setPattern(rows, cols);
  system.setValues(vals);
  EXPECT_TRUE(system.isSymmetric(1.0e-12));

  // Manufactured target and its image b = A x* (Dirichlet 0 outside).
  std::vector<real_t> xStar(static_cast<std::size_t>(n));
  for (int iy = 0; iy < m; ++iy) {
    for (int ix = 0; ix < m; ++ix) {
      const real_t x = (ix + 1) * h;
      const real_t y = (iy + 1) * h;
      xStar[static_cast<std::size_t>(gid(ix, iy))] =
          std::sin(std::numbers::pi * x) * std::sin(std::numbers::pi * y);
    }
  }
  auto starAt = [&xStar, &gid](int ix, int iy) -> real_t {
    const PetscInt g = gid(ix, iy);
    return g < 0 ? 0.0 : xStar[static_cast<std::size_t>(g)];
  };

  ValView rhs("rhs", static_cast<std::size_t>(n));
  ValView sol("sol", static_cast<std::size_t>(n));
  auto rhsH = Kokkos::create_mirror_view(rhs);
  for (int iy = 0; iy < m; ++iy) {
    for (int ix = 0; ix < m; ++ix) {
      rhsH(static_cast<std::size_t>(gid(ix, iy))) =
          invH2 * (4.0 * starAt(ix, iy) - starAt(ix - 1, iy) - starAt(ix + 1, iy) -
                   starAt(ix, iy - 1) - starAt(ix, iy + 1));
    }
  }
  Kokkos::deep_copy(rhs, rhsH);

  system.solve(rhs, sol);

  auto solH = Kokkos::create_mirror_view(sol);
  Kokkos::deep_copy(solH, sol);
  real_t errInf = 0.0;
  for (std::size_t r = 0; r < static_cast<std::size_t>(n); ++r) {
    errInf = std::max(errInf, std::abs(solH(r) - xStar[r]));
  }
  EXPECT_LT(errInf, 1.0e-9);
}

TEST(LinearSystem, AsymmetricMatrixDetected) {
  const int n = 4;
  const std::size_t nCoo = 6;
  IdxView rows("rows", nCoo);
  IdxView cols("cols", nCoo);
  ValView vals("vals", nCoo);
  auto rowsH = Kokkos::create_mirror_view(rows);
  auto colsH = Kokkos::create_mirror_view(cols);
  auto valsH = Kokkos::create_mirror_view(vals);
  for (int r = 0; r < n; ++r) {
    rowsH(static_cast<std::size_t>(r)) = r;
    colsH(static_cast<std::size_t>(r)) = r;
    valsH(static_cast<std::size_t>(r)) = 2.0;
  }
  rowsH(4) = 0; colsH(4) = 1; valsH(4) = -1.0;
  rowsH(5) = 1; colsH(5) = 0; valsH(5) = -0.5;  // breaks symmetry
  Kokkos::deep_copy(rows, rowsH);
  Kokkos::deep_copy(cols, colsH);
  Kokkos::deep_copy(vals, valsH);

  LinearSystem system(MPI_COMM_SELF, "tasym_", n, n);
  system.setPattern(rows, cols);
  system.setValues(vals);
  EXPECT_FALSE(system.isSymmetric(1.0e-12));
}

TEST(LinearSystem, ApiMisuseIsFatal) {
  LinearSystem system(MPI_COMM_SELF, "tmis_", 4, 4);
  ValView vals("vals", 4);
  EXPECT_THROW(system.setValues(vals), frehg::FatalError);  // pattern unset

  ValView rhs("rhs", 3);  // wrong size
  ValView sol("sol", 4);
  EXPECT_THROW(system.solve(rhs, sol), frehg::FatalError);

  IdxView rows("rows", 4);
  IdxView cols("cols", 5);  // mismatched lengths
  EXPECT_THROW(system.setPattern(rows, cols), frehg::FatalError);
}

}  // namespace
