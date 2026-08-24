/// \file LinearSystem.cpp
/// \brief Implementation of the PETSc COO-assembled linear system.

#include "core/LinearSystem.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <vector>

/// Every PETSc return code is checked; failures are fatal (plan §11.1
/// rule 5).
#define FREHG_PETSC_CHECK(call)                                                        \
  do {                                                                                 \
    const PetscErrorCode frehgPetscErr = (call);                                       \
    if (frehgPetscErr != PETSC_SUCCESS) {                                              \
      ::frehg::log::fatal(::frehg::log::msg()                                          \
                          << "PETSc error " << static_cast<int>(frehgPetscErr)         \
                          << " from " << #call);                                       \
    }                                                                                  \
  } while (0)

namespace frehg {

LinearSystem::LinearSystem(MPI_Comm comm, const std::string& prefix, PetscInt nLocalRows,
                           PetscInt nGlobalRows, const SolverSettings& settings)
    : comm_(comm), prefix_(prefix), nLocalRows_(nLocalRows), nGlobalRows_(nGlobalRows) {
  if (nLocalRows_ < 0 || nGlobalRows_ < 1) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): invalid sizes (local "
                          << nLocalRows_ << ", global " << nGlobalRows_ << ")");
  }

  FREHG_PETSC_CHECK(MatCreate(comm_, &mat_));
  FREHG_PETSC_CHECK(MatSetOptionsPrefix(mat_, prefix_.c_str()));
  FREHG_PETSC_CHECK(MatSetSizes(mat_, nLocalRows_, nLocalRows_, nGlobalRows_, nGlobalRows_));
  FREHG_PETSC_CHECK(MatSetType(mat_, MATAIJ));
  FREHG_PETSC_CHECK(MatSetFromOptions(mat_));

  FREHG_PETSC_CHECK(VecCreate(comm_, &rhs_));
  FREHG_PETSC_CHECK(VecSetOptionsPrefix(rhs_, prefix_.c_str()));
  FREHG_PETSC_CHECK(VecSetSizes(rhs_, nLocalRows_, nGlobalRows_));
  FREHG_PETSC_CHECK(VecSetFromOptions(rhs_));
  FREHG_PETSC_CHECK(VecDuplicate(rhs_, &sol_));

  FREHG_PETSC_CHECK(KSPCreate(comm_, &ksp_));
  FREHG_PETSC_CHECK(KSPSetOptionsPrefix(ksp_, prefix_.c_str()));
  FREHG_PETSC_CHECK(KSPSetType(ksp_, KSPCG));
  FREHG_PETSC_CHECK(KSPSetTolerances(ksp_, settings.rtol, settings.atol, PETSC_DEFAULT,
                                     settings.maxIterations));
  PC pc = nullptr;
  FREHG_PETSC_CHECK(KSPGetPC(ksp_, &pc));
  FREHG_PETSC_CHECK(PCSetType(pc, PCBJACOBI));
  // Default the block solver to icc(0) unless the user overrode it; this is
  // the plan §5.1 default preconditioner (hypre/gamg selectable through the
  // options file).
  {
    const std::string subKey = "-" + prefix_ + "sub_pc_type";
    PetscBool present = PETSC_FALSE;
    FREHG_PETSC_CHECK(PetscOptionsHasName(nullptr, nullptr, subKey.c_str(), &present));
    if (present == PETSC_FALSE) {
      FREHG_PETSC_CHECK(PetscOptionsSetValue(nullptr, subKey.c_str(), "icc"));
    }
  }
  FREHG_PETSC_CHECK(KSPSetFromOptions(ksp_));
}

LinearSystem::~LinearSystem() {
  PetscBool finalized = PETSC_FALSE;
  PetscFinalized(&finalized);
  if (finalized == PETSC_TRUE) {
    return;  // teardown after PetscFinalize: handles are already invalid
  }
  if (ksp_ != nullptr) {
    KSPDestroy(&ksp_);
  }
  if (sol_ != nullptr) {
    VecDestroy(&sol_);
  }
  if (rhs_ != nullptr) {
    VecDestroy(&rhs_);
  }
  if (mat_ != nullptr) {
    MatDestroy(&mat_);
  }
}

void LinearSystem::setPattern(const Kokkos::View<PetscInt*, MemSpace>& cooRows,
                              const Kokkos::View<PetscInt*, MemSpace>& cooCols) {
  if (cooRows.extent(0) != cooCols.extent(0)) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): COO row/col arrays differ ("
                          << cooRows.extent(0) << " vs " << cooCols.extent(0) << ")");
  }
  // PETSc may reorder the index arrays internally, so hand it scratch copies.
  Kokkos::View<PetscInt*, MemSpace> rows("coo_rows_scratch", cooRows.extent(0));
  Kokkos::View<PetscInt*, MemSpace> cols("coo_cols_scratch", cooCols.extent(0));
  Kokkos::deep_copy(rows, cooRows);
  Kokkos::deep_copy(cols, cooCols);
  Kokkos::fence();

  nCoo_ = static_cast<PetscCount>(cooRows.extent(0));
  FREHG_PETSC_CHECK(MatSetPreallocationCOO(mat_, nCoo_, rows.data(), cols.data()));
  patternSet_ = true;
}

void LinearSystem::setValues(const Kokkos::View<real_t*, MemSpace>& values) {
  if (!patternSet_) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_
                          << "'): setValues before setPattern");
  }
  if (static_cast<PetscCount>(values.extent(0)) != nCoo_) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): " << values.extent(0)
                          << " values for " << nCoo_ << " COO entries");
  }
  Kokkos::fence();
  FREHG_PETSC_CHECK(MatSetValuesCOO(mat_, values.data(), INSERT_VALUES));
}

SolveStats LinearSystem::solve(const Kokkos::View<real_t*, MemSpace>& rhs,
                               const Kokkos::View<real_t*, MemSpace>& solution) {
  if (rhs.extent(0) != static_cast<std::size_t>(nLocalRows_) ||
      solution.extent(0) != static_cast<std::size_t>(nLocalRows_)) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): rhs/solution sized "
                          << rhs.extent(0) << "/" << solution.extent(0) << ", expected "
                          << nLocalRows_);
  }

#ifndef NDEBUG
  if (!firstSolveChecked_) {
    firstSolveChecked_ = true;
    if (!isSymmetric(1.0e-10)) {
      log::fatal(log::msg() << "LinearSystem('" << prefix_
                            << "'): matrix is not symmetric on the first solve; both Frehg2 "
                               "systems must be SPD by construction (plan §5.1)");
    }
  }
#endif

  // Stage the right-hand side through a host mirror into the PETSc vector.
  auto rhsHost = Kokkos::create_mirror_view(rhs);
  Kokkos::deep_copy(rhsHost, rhs);
  {
    PetscScalar* ptr = nullptr;
    FREHG_PETSC_CHECK(VecGetArrayWrite(rhs_, &ptr));
    for (PetscInt n = 0; n < nLocalRows_; ++n) {
      ptr[n] = rhsHost(static_cast<std::size_t>(n));
    }
    FREHG_PETSC_CHECK(VecRestoreArrayWrite(rhs_, &ptr));
  }

  FREHG_PETSC_CHECK(KSPSetOperators(ksp_, mat_, mat_));
  FREHG_PETSC_CHECK(KSPSolve(ksp_, rhs_, sol_));

  KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
  FREHG_PETSC_CHECK(KSPGetConvergedReason(ksp_, &reason));
  PetscInt iterations = 0;
  FREHG_PETSC_CHECK(KSPGetIterationNumber(ksp_, &iterations));
  PetscReal residual = 0.0;
  FREHG_PETSC_CHECK(KSPGetResidualNorm(ksp_, &residual));

  if (reason < 0) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): solver diverged after "
                          << iterations << " iterations, residual " << residual << ", reason "
                          << KSPConvergedReasons[reason]
                          << "; aborting (silent continuation was a legacy bug, plan §5.1)");
  }

  auto solHost = Kokkos::create_mirror_view(solution);
  {
    const PetscScalar* ptr = nullptr;
    FREHG_PETSC_CHECK(VecGetArrayRead(sol_, &ptr));
    for (PetscInt n = 0; n < nLocalRows_; ++n) {
      solHost(static_cast<std::size_t>(n)) = ptr[n];
    }
    FREHG_PETSC_CHECK(VecRestoreArrayRead(sol_, &ptr));
  }
  Kokkos::deep_copy(solution, solHost);
  Kokkos::fence();

  SolveStats stats;
  stats.iterations = static_cast<int>(iterations);
  stats.residualNorm = static_cast<real_t>(residual);
  stats.reason = KSPConvergedReasons[reason];
  return stats;
}

bool LinearSystem::isSymmetric(real_t tolerance) const {
  // PETSc's MatIsSymmetric rejects some parallel matrix types, so measure
  // ||A - A^T||_inf directly.
  Mat transposed = nullptr;
  FREHG_PETSC_CHECK(MatTranspose(mat_, MAT_INITIAL_MATRIX, &transposed));
  FREHG_PETSC_CHECK(MatAXPY(transposed, -1.0, mat_, DIFFERENT_NONZERO_PATTERN));
  PetscReal diffNorm = 0.0;
  FREHG_PETSC_CHECK(MatNorm(transposed, NORM_INFINITY, &diffNorm));
  PetscReal matNorm = 0.0;
  FREHG_PETSC_CHECK(MatNorm(mat_, NORM_INFINITY, &matNorm));
  FREHG_PETSC_CHECK(MatDestroy(&transposed));
  return diffNorm <= tolerance * (1.0 + matNorm);
}

}  // namespace frehg
