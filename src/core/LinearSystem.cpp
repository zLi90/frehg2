/// \file LinearSystem.cpp
/// \brief Implementation of the PETSc COO-assembled linear system.

// PETSc requires petscvec_kokkos.hpp to precede every other PETSc header in
// the translation unit, so probe petscconf.h (macro-only, no declarations)
// for Kokkos support before the class header pulls in petscksp.h.
#include <petscconf.h>
#if defined(PETSC_HAVE_KOKKOS)
#include <petscvec_kokkos.hpp>
#endif

#include "core/LinearSystem.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
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

namespace {

/// True when simulation fields live in device memory (CUDA/HIP builds).
/// Drives forced Kokkos matrix/vector types (B1) and the backend-aware
/// BoomerAMG defaults (B2) — v2 plan §2B.2.
constexpr bool kOnDevice = !std::is_same_v<MemSpace, Kokkos::HostSpace>;

/// Whether this PETSc was built with Kokkos support.
#if defined(PETSC_HAVE_KOKKOS)
constexpr bool kPetscHasKokkos = true;
#else
constexpr bool kPetscHasKokkos = false;
#endif

/// The p5 compile-time invariant (v2 plan §2B.2 B1): a device build must
/// link a Kokkos-enabled PETSc. Fails the build, not the first solve.
static_assert(PetscBackendConsistent<kPetscHasKokkos, MemSpace>::value);

/// Seconds since an arbitrary epoch, for the setup/solve time split.
double nowSeconds() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

/// Set a prefixed PETSc option only when the user has not already set it,
/// so defaults injected here never override an options file or the command
/// line (the v1 sub_pc_type pattern, generalized).
void setOptionDefault(const std::string& prefix, const std::string& key,
                      const std::string& value) {
  const std::string full = "-" + prefix + key;
  PetscBool present = PETSC_FALSE;
  FREHG_PETSC_CHECK(PetscOptionsHasName(nullptr, nullptr, full.c_str(), &present));
  if (present == PETSC_FALSE) {
    FREHG_PETSC_CHECK(PetscOptionsSetValue(nullptr, full.c_str(), value.c_str()));
  }
}

}  // namespace

LinearSystem::LinearSystem(MPI_Comm comm, const std::string& prefix, PetscInt nLocalRows,
                           PetscInt nGlobalRows, const SolverSettings& settings)
    : comm_(comm), prefix_(prefix), nLocalRows_(nLocalRows), nGlobalRows_(nGlobalRows) {
  if (nLocalRows_ < 0 || nGlobalRows_ < 1) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): invalid sizes (local "
                          << nLocalRows_ << ", global " << nGlobalRows_ << ")");
  }

  // Matrix/vector backend selection (v2 plan §2B.2 B1). Device builds force
  // the Kokkos types: the COO assembly path (setValues) hands PETSc a
  // MemSpace pointer, which a host MATAIJ would misread as host memory. On
  // host builds the type follows settings.matType ("aij" default;
  // "aijkokkos" runs the solve through Kokkos Kernels on the host execution
  // space — the threaded-CPU lane). Injected via setOptionDefault so
  // -<prefix>mat_type / -<prefix>vec_type on the command line or in an
  // options file still win; the resolved truth is read back below.
  std::string requestedMatType = settings.matType;
  if (kOnDevice) {
    requestedMatType = "aijkokkos";
  }
  if (requestedMatType == "aijkokkos" && !kPetscHasKokkos) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): mat_type 'aijkokkos' "
                          << "needs a Kokkos-enabled PETSc (configure with "
                             "--with-kokkos-dir or --download-kokkos plus "
                             "--download-kokkos-kernels — v2 plan §2B.2)");
  } else if (requestedMatType == "aijkokkos") {
    setOptionDefault(prefix_, "mat_type", "aijkokkos");
    setOptionDefault(prefix_, "vec_type", "kokkos");
  } else if (requestedMatType != "aij") {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): unknown mat_type '"
                          << requestedMatType << "' (aij | aijkokkos)");
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

  // Read back what PETSc actually instantiated (an options-database
  // override may have changed it): the run record reports this resolved
  // type, and the solve() staging path keys off the vector type.
  {
    MatType matTypeActual = nullptr;
    FREHG_PETSC_CHECK(MatGetType(mat_, &matTypeActual));
    matType_ = matTypeActual;
    VecType vecTypeActual = nullptr;
    FREHG_PETSC_CHECK(VecGetType(rhs_, &vecTypeActual));
    kokkosVec_ = std::strstr(vecTypeActual, "kokkos") != nullptr;
    if (kOnDevice && !kokkosVec_) {
      log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): device build with a "
                            << "non-Kokkos vector type '" << vecTypeActual
                            << "' — a -" << prefix_ << "vec_type override cannot "
                               "select host types on a device build (v2 plan §2B.2)");
    }
  }

  FREHG_PETSC_CHECK(KSPCreate(comm_, &ksp_));
  FREHG_PETSC_CHECK(KSPSetOptionsPrefix(ksp_, prefix_.c_str()));
  FREHG_PETSC_CHECK(KSPSetType(ksp_, KSPCG));
  FREHG_PETSC_CHECK(KSPSetTolerances(ksp_, settings.rtol, settings.atol, PETSC_DEFAULT,
                                     settings.maxIterations));
  PC pc = nullptr;
  FREHG_PETSC_CHECK(KSPGetPC(ksp_, &pc));

  // Preconditioner selection (v2 plan §2.2). Every default injected below
  // uses setOptionDefault, so a PETSc options file or command line always
  // wins; KSPSetFromOptions runs last for the same reason.
  if (settings.preconditioner == "amg") {
    // hypre BoomerAMG: classical AMG, the primary scalable path. Fail with
    // a plain message when PETSc lacks hypre rather than deep inside
    // PCSetUp (v2 plan §2.2.2).
    PetscBool hasHypre = PETSC_FALSE;
    FREHG_PETSC_CHECK(PetscHasExternalPackage("hypre", &hasHypre));
    if (hasHypre == PETSC_FALSE) {
      log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): preconditioner 'amg' "
                            << "needs a hypre-enabled PETSc (configure with "
                               "--download-hypre, scripts/ci_install_deps.sh) — or "
                               "select 'gamg', which is built into PETSc");
    }
    FREHG_PETSC_CHECK(PCSetType(pc, PCHYPRE));
    FREHG_PETSC_CHECK(PCHYPRESetType(pc, "boomeramg"));
    // The v2 plan §2.2.3 baseline, made backend-aware in §2B.2 B2: HMIS
    // coarsening on host builds; on device builds PMIS, the only coarsening
    // hypre implements on GPUs (HMIS would silently fall back or error).
    // ext+i distance-two interpolation is device-supported and stays for
    // both. Aggressive-coarsening interpolation needs no forcing here:
    // PETSc >= 3.25 selects a device-capable two-stage MM operator itself
    // when the matrix lives on device (agg_interptype 7 vs multipass 4).
    setOptionDefault(prefix_, "pc_hypre_boomeramg_coarsen_type",
                     kOnDevice ? "PMIS" : "HMIS");
    setOptionDefault(prefix_, "pc_hypre_boomeramg_interp_type", "ext+i");
    setOptionDefault(prefix_, "pc_hypre_boomeramg_P_max", "4");
    setOptionDefault(prefix_, "pc_hypre_boomeramg_strong_threshold",
                     std::to_string(settings.amgStrongThreshold));
    setOptionDefault(prefix_, "pc_hypre_boomeramg_agg_nl",
                     std::to_string(settings.amgAggressiveLevels));
    // Smoother (v2 plan §2B.2 B2): hypre's host default (hybrid symmetric
    // SOR/Jacobi) is thread-count-DEPENDENT by construction — Jacobi between
    // OpenMP threads, Gauss-Seidel within — so iteration counts drift with
    // OMP_NUM_THREADS, which breaks p2's thread-invariance assertion. On
    // the Kokkos lanes (host-threaded or device) default to l1-scaled
    // Jacobi: symmetric (CG-safe), thread-invariant, and the same smoother
    // PETSc selects on device — the OpenMP lane then rehearses the exact
    // GPU algebra. The host "aij" lane keeps hypre's default (golden-pinned
    // Q1 g2 behavior).
    if (kokkosVec_) {
      setOptionDefault(prefix_, "pc_hypre_boomeramg_relax_type_all", "l1scaled-Jacobi");
    }
  } else if (settings.preconditioner == "gamg") {
    // PETSc smoothed aggregation: no external package, and the native
    // device path once a Kokkos-enabled PETSc exists (v2 plan §2.2.2).
    FREHG_PETSC_CHECK(PCSetType(pc, PCGAMG));
    setOptionDefault(prefix_, "pc_gamg_type", "agg");
    setOptionDefault(prefix_, "pc_gamg_agg_nsmooths", "1");
    setOptionDefault(prefix_, "pc_gamg_threshold", "0.02");
  } else if (settings.preconditioner == "bjacobi-icc") {
    FREHG_PETSC_CHECK(PCSetType(pc, PCBJACOBI));
    // Default the block solver to icc(0) unless the user overrode it; this
    // is the plan §5.1 default preconditioner.
    setOptionDefault(prefix_, "sub_pc_type", "icc");
    // On the Kokkos lanes, keep the IC(0) factorization on PETSc's host
    // path: Kokkos Kernels' IC(0) produces a measurably weaker factor on
    // ill-scaled systems — CG stalls an order short of rtol on the b6 Kuan
    // matrix (dy = 0.05, dz = 0.04 cells), caught by gate p1. Host factors
    // make the algebra identical to the aij lane by construction; SpMV and
    // vector work still run through Kokkos. The device-intended
    // preconditioner is amg (v2 plan §2B.2 B2) — sequential triangular
    // solves are not a device path either way.
    if (kokkosVec_) {
      setOptionDefault(prefix_, "sub_pc_factor_mat_solver_type", "petsc");
    }
  } else {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): unknown preconditioner '"
                          << settings.preconditioner
                          << "' (bjacobi-icc | amg | gamg)");
  }
  FREHG_PETSC_CHECK(KSPSetFromOptions(ksp_));

  // Read back the resolved BoomerAMG choices for the run record (v2 plan
  // §2B.2 B2): the options database at this point holds whichever value won
  // (our setOptionDefault, an options file, or the command line), so a run
  // is self-describing about which AMG it actually got.
  if (settings.preconditioner == "amg") {
    const auto resolvedOption = [this](const char* key) {
      char buf[128] = {0};
      PetscBool set = PETSC_FALSE;
      const std::string full = "-" + prefix_ + key;
      PetscOptionsGetString(nullptr, nullptr, full.c_str(), buf, sizeof(buf), &set);
      return std::string(set == PETSC_TRUE ? buf : "petsc-default");
    };
    amgCoarsenType_ = resolvedOption("pc_hypre_boomeramg_coarsen_type");
    amgRelaxType_ = resolvedOption("pc_hypre_boomeramg_relax_type_all");
  }

  // Hierarchy reuse applies to the AMG-class preconditioners only: their
  // setup costs several solves, while bjacobi-icc's per-solve refresh is
  // cheap and is the golden-pinned v1 behavior (v2 plan §2.2.4).
  reusable_ = (settings.preconditioner != "bjacobi-icc") && settings.reuseMaxSolves > 0;
  reuseMaxSolves_ = settings.reuseMaxSolves;
  reuseIterationFactor_ = settings.reuseIterationFactor;
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

PetscInt LinearSystem::solveOnce(KSPConvergedReason& reason, PetscReal& residual) {
  // KSPSetUp performs the preconditioner setup when the operator changed
  // (or when a rebuild was forced); timing it separately from KSPSolve
  // gives the g3 setup-fraction telemetry (v2 plan §2.3).
  const double setupStart = nowSeconds();
  FREHG_PETSC_CHECK(KSPSetUp(ksp_));
  telemetry_.setupSeconds += nowSeconds() - setupStart;

  const double solveStart = nowSeconds();
  FREHG_PETSC_CHECK(KSPSolve(ksp_, rhs_, sol_));
  telemetry_.solveSeconds += nowSeconds() - solveStart;

  reason = KSP_CONVERGED_ITERATING;
  FREHG_PETSC_CHECK(KSPGetConvergedReason(ksp_, &reason));
  PetscInt iterations = 0;
  FREHG_PETSC_CHECK(KSPGetIterationNumber(ksp_, &iterations));
  residual = 0.0;
  FREHG_PETSC_CHECK(KSPGetResidualNorm(ksp_, &residual));
  return iterations;
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
    // p5 memtype discipline (v2 plan §2B.3): on a device build the PETSc
    // vector storage must be device memory — a host vector here means the
    // solve would silently stage through the wrong memory space.
    if (kOnDevice) {
      const PetscScalar* probe = nullptr;
      PetscMemType memType = PETSC_MEMTYPE_HOST;
      FREHG_PETSC_CHECK(VecGetArrayReadAndMemType(rhs_, &probe, &memType));
      FREHG_PETSC_CHECK(VecRestoreArrayReadAndMemType(rhs_, &probe));
      if (!PetscMemTypeDevice(memType)) {
        log::fatal(log::msg() << "LinearSystem('" << prefix_
                              << "'): device build but the PETSc vector reports host "
                                 "memory — matrix/vector type resolution is broken "
                                 "(v2 plan §2B.2 B1)");
      }
    }
    if (!isSymmetric(1.0e-10)) {
      log::fatal(log::msg() << "LinearSystem('" << prefix_
                            << "'): matrix is not symmetric on the first solve; both Frehg2 "
                               "systems must be SPD by construction (plan §5.1)");
    }
  }
#endif

  // Move the right-hand side into the PETSc vector. With VECKOKKOS the
  // handoff stays in MemSpace through PETSc's Kokkos view API — on a GPU
  // build this is the difference between a device-side copy and two
  // device<->host round trips per solve (v2 plan §2B.2 B1). With host
  // vectors, stage through a host mirror (a no-op mirror on CPU builds).
#if defined(PETSC_HAVE_KOKKOS)
  if (kokkosVec_) {
    Kokkos::View<PetscScalar*, MemSpace> view;
    FREHG_PETSC_CHECK(VecGetKokkosViewWrite(rhs_, &view));
    Kokkos::deep_copy(view, rhs);
    FREHG_PETSC_CHECK(VecRestoreKokkosViewWrite(rhs_, &view));
  } else
#endif
  {
    auto rhsHost = Kokkos::create_mirror_view(rhs);
    Kokkos::deep_copy(rhsHost, rhs);
    PetscScalar* ptr = nullptr;
    FREHG_PETSC_CHECK(VecGetArrayWrite(rhs_, &ptr));
    for (PetscInt n = 0; n < nLocalRows_; ++n) {
      ptr[n] = rhsHost(static_cast<std::size_t>(n));
    }
    FREHG_PETSC_CHECK(VecRestoreArrayWrite(rhs_, &ptr));
  }

  FREHG_PETSC_CHECK(KSPSetOperators(ksp_, mat_, mat_));

  // Hierarchy-reuse policy (v2 plan §2.2.4). Without reuse (bjacobi-icc, or
  // reuse_max_solves 0) PETSc rebuilds the PC whenever the operator changed
  // — the v1 behavior. With reuse, the hierarchy is frozen after each build
  // and refreshed on the schedule below. Every branch here is driven by
  // KSP-collective quantities (iteration counts, solve counts, the caller's
  // Allreduced forceRebuild), so all ranks take the same path.
  bool rebuilt = false;
  if (reusable_) {
    if (rebuildNext_) {
      FREHG_PETSC_CHECK(KSPSetReusePreconditioner(ksp_, PETSC_FALSE));
      rebuilt = true;
    }
  }

  KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
  PetscReal residual = 0.0;
  PetscInt iterations = solveOnce(reason, residual);

  if (reason < 0 && reusable_ && !rebuilt) {
    // A reused hierarchy can go stale (the operator drifted past what the
    // frozen coarse spaces represent). Divergence under a *fresh*
    // preconditioner stays fatal below; under a reused one, rebuild and
    // retry once before giving up.
    log::warn(log::msg() << "LinearSystem('" << prefix_ << "'): diverged ("
                         << KSPConvergedReasons[reason]
                         << ") with a reused preconditioner; rebuilding and retrying");
    FREHG_PETSC_CHECK(KSPSetReusePreconditioner(ksp_, PETSC_FALSE));
    FREHG_PETSC_CHECK(KSPSetOperators(ksp_, mat_, mat_));
    rebuilt = true;
    ++telemetry_.retries;
    iterations = solveOnce(reason, residual);
  }

  if (reason < 0) {
    log::fatal(log::msg() << "LinearSystem('" << prefix_ << "'): solver diverged after "
                          << iterations << " iterations, residual " << residual << ", reason "
                          << KSPConvergedReasons[reason]
                          << "; aborting (silent continuation was a legacy bug, plan §5.1)");
  }

  if (reusable_) {
    if (rebuilt) {
      // Freeze the fresh hierarchy and record its iteration baseline.
      FREHG_PETSC_CHECK(KSPSetReusePreconditioner(ksp_, PETSC_TRUE));
      baselineIterations_ = static_cast<int>(iterations);
      solvesSinceRebuild_ = 0;
      rebuildNext_ = false;
      ++telemetry_.rebuilds;
    }
    ++solvesSinceRebuild_;
    const real_t trigger =
        reuseIterationFactor_ * static_cast<real_t>(baselineIterations_) + 2.0;
    if (static_cast<real_t>(iterations) > trigger ||
        solvesSinceRebuild_ >= reuseMaxSolves_) {
      rebuildNext_ = true;
    }
  } else {
    // Without reuse every solve that saw a changed operator set up the PC.
    ++telemetry_.rebuilds;
  }

  ++telemetry_.solves;
  telemetry_.totalIterations += static_cast<long>(iterations);
  telemetry_.maxIterations = std::max(telemetry_.maxIterations, static_cast<int>(iterations));

  // Read the solution back through the same memtype-correct path.
#if defined(PETSC_HAVE_KOKKOS)
  if (kokkosVec_) {
    Kokkos::View<const PetscScalar*, MemSpace> view;
    FREHG_PETSC_CHECK(VecGetKokkosView(sol_, &view));
    Kokkos::deep_copy(solution, view);
    FREHG_PETSC_CHECK(VecRestoreKokkosView(sol_, &view));
  } else
#endif
  {
    auto solHost = Kokkos::create_mirror_view(solution);
    const PetscScalar* ptr = nullptr;
    FREHG_PETSC_CHECK(VecGetArrayRead(sol_, &ptr));
    for (PetscInt n = 0; n < nLocalRows_; ++n) {
      solHost(static_cast<std::size_t>(n)) = ptr[n];
    }
    FREHG_PETSC_CHECK(VecRestoreArrayRead(sol_, &ptr));
    Kokkos::deep_copy(solution, solHost);
  }
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
