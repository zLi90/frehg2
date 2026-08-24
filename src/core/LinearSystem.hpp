/// \file LinearSystem.hpp
/// \brief PETSc KSP wrapper with device-capable COO assembly (plan §5.1).
///
/// One LinearSystem instance wraps a PETSc Mat/Vec/KSP triplet. The sparsity
/// pattern is fixed once with MatSetPreallocationCOO (both Frehg2 stencils
/// are static: 5-point surface, 7-point subsurface); each solve updates the
/// values with MatSetValuesCOO from a Kokkos view. Out-of-domain stencil legs
/// use COO index -1 (ignored by PETSc) with their coefficient folded into the
/// diagonal by the assembly kernel; duplicate (i, j) entries are summed,
/// which handles symmetric scatter.
///
/// With `-<prefix>mat_type aijkokkos` and a Kokkos-enabled PETSc, the same
/// calls accept device pointers, so the CPU-to-GPU switch is a runtime
/// option, not a code path.

#ifndef FREHG_CORE_LINEARSYSTEM_HPP
#define FREHG_CORE_LINEARSYSTEM_HPP

#include "core/Types.hpp"

#include <mpi.h>
#include <petscksp.h>

#include <string>

namespace frehg {

/// Convergence result of one linear solve.
struct SolveStats {
  int iterations = 0;        ///< KSP iterations performed
  real_t residualNorm = 0.0; ///< final (preconditioned) residual norm
  std::string reason;        ///< PETSc KSPConvergedReason string
};

/// Solver tolerances applied at construction (overridable per prefix through
/// the PETSc options database, plan §5.1 table).
struct SolverSettings {
  real_t rtol = 1.0e-8;     ///< relative tolerance (matches legacy SetRTCAccuracy)
  real_t atol = 1.0e-14;    ///< absolute tolerance
  int maxIterations = 500;  ///< iteration cap
};

/// A distributed sparse linear system with COO assembly.
class LinearSystem {
 public:
  /// Create the Mat/Vec/KSP set.
  /// \param comm communicator matching the grid decomposition.
  /// \param prefix PETSc options prefix, e.g. "fs_" or "gw_".
  /// \param nLocalRows rows owned by this rank (grid activeCount*Local()).
  /// \param nGlobalRows total rows (grid activeCount*Global()).
  /// \param settings default tolerances (KSP CG, bjacobi + icc(0) unless
  ///        overridden via options).
  LinearSystem(MPI_Comm comm, const std::string& prefix, PetscInt nLocalRows,
               PetscInt nGlobalRows, const SolverSettings& settings = SolverSettings{});

  ~LinearSystem();
  LinearSystem(const LinearSystem&) = delete;
  LinearSystem& operator=(const LinearSystem&) = delete;

  /// Fix the sparsity pattern from COO index arrays (device or host views).
  /// Entries with row or column -1 are ignored by PETSc. Called once; the
  /// number of values passed to setValues() must equal cooRows.extent(0).
  void setPattern(const Kokkos::View<PetscInt*, MemSpace>& cooRows,
                  const Kokkos::View<PetscInt*, MemSpace>& cooCols);

  /// Update the matrix values (same ordering as the pattern arrays;
  /// duplicate indices are summed).
  void setValues(const Kokkos::View<real_t*, MemSpace>& values);

  /// Solve A x = b.
  /// \param rhs local right-hand-side values (nLocalRows entries).
  /// \param solution local solution values, overwritten (nLocalRows).
  /// \return iteration count, residual norm, and converged reason.
  /// Divergence is fatal: the plan treats silent solver failure as a legacy
  /// bug, not fidelity (§5.1).
  SolveStats solve(const Kokkos::View<real_t*, MemSpace>& rhs,
                   const Kokkos::View<real_t*, MemSpace>& solution);

  /// \return true if the assembled matrix is symmetric within \p tolerance.
  /// Runs PETSc's MatIsSymmetric; used by the debug-build step-1 check and
  /// by tests (both Frehg2 systems are SPD by construction).
  bool isSymmetric(real_t tolerance) const;

  /// \return rows owned by this rank.
  PetscInt localRows() const { return nLocalRows_; }

  /// \return total rows.
  PetscInt globalRows() const { return nGlobalRows_; }

  /// \return the options prefix.
  const std::string& prefix() const { return prefix_; }

 private:
  MPI_Comm comm_ = MPI_COMM_NULL;
  std::string prefix_;
  PetscInt nLocalRows_ = 0;
  PetscInt nGlobalRows_ = 0;
  PetscCount nCoo_ = 0;
  bool patternSet_ = false;
  /// Only consulted by the debug-build first-solve symmetry check.
  [[maybe_unused]] bool firstSolveChecked_ = false;

  Mat mat_ = nullptr;
  Vec rhs_ = nullptr;
  Vec sol_ = nullptr;
  KSP ksp_ = nullptr;
};

}  // namespace frehg

#endif  // FREHG_CORE_LINEARSYSTEM_HPP
