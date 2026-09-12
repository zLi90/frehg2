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
#include <type_traits>

namespace frehg {

/// Compile-time backend-consistency invariant (v2 plan §2B.2 B1, gate p5):
/// a non-host MemSpace (CUDA/HIP build) requires a Kokkos-enabled PETSc,
/// because the COO assembly path hands MemSpace pointers to MatSetValuesCOO
/// and a host MATAIJ would read them as host memory. LinearSystem.cpp
/// instantiates this with the build's actual (PETSC_HAVE_KOKKOS, MemSpace)
/// pair; the p5 negative test instantiates the forbidden combination and
/// asserts the compile fails.
template <bool PetscHasKokkos, class Mem>
struct PetscBackendConsistent {
  static constexpr bool value = PetscHasKokkos || std::is_same_v<Mem, Kokkos::HostSpace>;
  static_assert(PetscHasKokkos || std::is_same_v<Mem, Kokkos::HostSpace>,
                "Frehg2: a device MemSpace requires a Kokkos-enabled PETSc "
                "(PETSC_HAVE_KOKKOS). Rebuild PETSc with --with-kokkos-dir (or "
                "--download-kokkos) --download-kokkos-kernels and the matching "
                "device backend (v2 plan section 2B.2).");
};

/// Convergence result of one linear solve.
struct SolveStats {
  int iterations = 0;        ///< KSP iterations performed
  real_t residualNorm = 0.0; ///< final (preconditioned) residual norm
  std::string reason;        ///< PETSc KSPConvergedReason string
};

/// Cumulative per-system solver telemetry (v2 plan §2.2.6): the g2/g3
/// gates and the end-of-run "solver summary" lines read these. Iteration
/// counts are identical on every rank (KSP is collective); times are the
/// local rank's.
struct SolverTelemetry {
  long solves = 0;           ///< linear solves performed
  long totalIterations = 0;  ///< sum of iteration counts
  int maxIterations = 0;     ///< largest single-solve iteration count
  long rebuilds = 0;         ///< preconditioner setups (== solves without reuse)
  int retries = 0;           ///< stale-hierarchy divergences retried after rebuild
  double setupSeconds = 0.0; ///< time in KSPSetUp (preconditioner setup)
  double solveSeconds = 0.0; ///< time in KSPSolve
};

/// Solver tolerances and preconditioner selection applied at construction
/// (v2 plan §2.2; every PETSc detail stays overridable per prefix through
/// the options database — KSPSetFromOptions runs last).
struct SolverSettings {
  real_t rtol = 1.0e-8;     ///< relative tolerance (matches legacy SetRTCAccuracy)
  real_t atol = 1.0e-14;    ///< absolute tolerance
  int maxIterations = 500;  ///< iteration cap
  /// "bjacobi-icc" (the v1 default: block-Jacobi with ICC(0) blocks),
  /// "amg" (hypre BoomerAMG; requires a hypre-enabled PETSc), or "gamg"
  /// (PETSc's built-in smoothed aggregation).
  std::string preconditioner = "bjacobi-icc";
  /// PETSc matrix/vector backend (v2 plan §2B.2 B1): "aij" (host default)
  /// or "aijkokkos" (Kokkos Kernels on the build's execution space;
  /// requires PETSC_HAVE_KOKKOS). Device builds force "aijkokkos".
  std::string matType = "aij";
  /// BoomerAMG strength-of-connection threshold default: 0.25 suits the 2D
  /// 5-point system; the 3D thin-layer Richards system passes 0.5 (hypre's
  /// 3D guidance; below 0.5 in 3D risks complexity blowup — v2 plan §2.2.3).
  real_t amgStrongThreshold = 0.25;
  /// BoomerAMG aggressive-coarsening levels (0 = none; the 3D system uses 1).
  int amgAggressiveLevels = 0;
  /// AMG hierarchy reuse cadence: rebuild at least every this many solves
  /// (0 = rebuild every solve). Ignored for bjacobi-icc, whose per-solve
  /// refresh is cheap and is the golden-pinned v1 behavior.
  int reuseMaxSolves = 50;
  /// Early rebuild when iterations exceed this factor times the count
  /// measured right after the last rebuild (plus a +2 absolute allowance).
  real_t reuseIterationFactor = 1.5;
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

  /// \return true when the preconditioner hierarchy is reused across solves
  /// (amg/gamg with reuse_max_solves > 0).
  bool reusesHierarchy() const { return reusable_; }

  /// Request a preconditioner rebuild before the next solve. Collective in
  /// effect: callers must reach the same decision on every rank of the
  /// communicator (PCSetUp is collective), e.g. from an Allreduced
  /// wet/dry-mask hash (v2 plan §2.2.4). No-op without hierarchy reuse.
  void forceRebuild() { rebuildNext_ = true; }

  /// \return the cumulative solver telemetry (v2 plan §2.2.6).
  const SolverTelemetry& telemetry() const { return telemetry_; }

  /// \return rows owned by this rank.
  PetscInt localRows() const { return nLocalRows_; }

  /// \return total rows.
  PetscInt globalRows() const { return nGlobalRows_; }

  /// \return the options prefix.
  const std::string& prefix() const { return prefix_; }

  /// \return the resolved matrix/vector backend actually in use ("aij" or
  /// "aijkokkos"; device builds force the latter — v2 plan §2B.2 B1). Feeds
  /// the run record's solver section.
  const std::string& matType() const { return matType_; }

  /// \return the resolved BoomerAMG coarsening / smoother choice (empty
  /// unless preconditioner == "amg"; "petsc-default" when nothing set it).
  /// Feeds the run record's solver section (v2 plan §2B.2 B2).
  const std::string& amgCoarsenType() const { return amgCoarsenType_; }
  const std::string& amgRelaxType() const { return amgRelaxType_; }

 private:
  /// One timed KSPSetUp + KSPSolve pass; returns the converged reason and
  /// fills iterations/residual. Factored out so the stale-hierarchy retry
  /// path (v2 plan §2.2.4) can rerun it after a forced rebuild.
  PetscInt solveOnce(KSPConvergedReason& reason, PetscReal& residual);

  MPI_Comm comm_ = MPI_COMM_NULL;
  std::string prefix_;
  /// Resolved matrix/vector backend ("aij" | "aijkokkos"); when
  /// "aijkokkos", solve() moves the RHS/solution through PETSc's Kokkos
  /// view API instead of host mirrors (v2 plan §2B.2 B1).
  std::string matType_ = "aij";
  /// Resolved BoomerAMG coarsening/smoother (amg only; run-record fields).
  std::string amgCoarsenType_;
  std::string amgRelaxType_;
  bool kokkosVec_ = false;  ///< rhs_/sol_ are VECKOKKOS (read back at construction)
  PetscInt nLocalRows_ = 0;
  PetscInt nGlobalRows_ = 0;
  PetscCount nCoo_ = 0;
  bool patternSet_ = false;
  /// Only consulted by the debug-build first-solve symmetry check.
  [[maybe_unused]] bool firstSolveChecked_ = false;

  // Hierarchy-reuse state (v2 plan §2.2.4). rebuildNext_ starts true so the
  // first solve always builds.
  bool reusable_ = false;
  bool rebuildNext_ = true;
  int solvesSinceRebuild_ = 0;
  int baselineIterations_ = 0;  ///< iteration count right after the last rebuild
  int reuseMaxSolves_ = 0;
  real_t reuseIterationFactor_ = 1.5;

  SolverTelemetry telemetry_;

  Mat mat_ = nullptr;
  Vec rhs_ = nullptr;
  Vec sol_ = nullptr;
  KSP ksp_ = nullptr;
};

}  // namespace frehg

#endif  // FREHG_CORE_LINEARSYSTEM_HPP
