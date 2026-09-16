#!/usr/bin/env bash
# ci_build_and_test.sh — the per-PR machine-checkable gate (plan §11.2):
# strict build, unit + mpi labels, fast regressions (P1+), forbidden scan,
# parameter-docs lockstep, and Doxygen with warnings as errors.
#
# The compiler comes from the environment. Dependencies must be consumed
# with the compiler that built them; e.g. with gcc-built dependencies at
# $HOME/frehg-deps, run with
#   CXX=g++-15 CC=gcc-15 PATH=$HOME/frehg-deps/bin:$PATH \
#   CMAKE_PREFIX_PATH=$HOME/frehg-deps scripts/ci_build_and_test.sh
# (report-P0.md §4 explains why CXX must not be the mpicxx wrapper).

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-$ROOT/build-ci}"

# Runtime loader path for the shared dependency libraries. Kokkos is built
# shared (A8 invariant 7: one libkokkoscore runtime in the process, or
# VecKokkos double-inits the singleton and segfaults), and frehg pulls
# libkokkoscontainers.so et al. in transitively through PETSc's pkg-config -L
# flags. CMake does not turn those -L flags into an rpath -- unlike PETSc's own
# libpetsc.so, which carries the rpath PETSc bakes in -- so the ctest runs
# below cannot find the Kokkos .so's without the dep lib dirs on
# LD_LIBRARY_PATH. Derived from CMAKE_PREFIX_PATH (the documented invocation
# above sets it); mirrors build_frehg2_slurm.sh, which exports the same.
_pfx="${CMAKE_PREFIX_PATH:-}"
for _p in ${_pfx//:/ }; do
  [ -n "$_p" ] && LD_LIBRARY_PATH="$_p/lib:$_p/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
done
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Prepend the system MPICH runtime dir so the frehg binary resolves libmpi to
# the same MPICH that mpiexec.mpich launches under. The dependency prefix is on
# LD_LIBRARY_PATH above (for the Kokkos .so's); if a stray libmpi ever lands
# there it would shadow the system one, and the binary would MPI_Init under a
# different MPICH than the launcher's PMI -- every rank then degrades to a
# size-1 MPI_COMM_WORLD ("1 MPI rank" on all ranks) and the parallel-HDF5 tests
# race on a colliding filename (mpi.core.n4). Derived from the MPICH wrapper's
# own -L flags; a no-op when the prefix is clean (same dir the loader would pick
# anyway). DT_RPATH in libpetsc still wins over this, so it does not mask a
# downloaded-and-rpath'd MPI -- the diagnostics below surface that case.
_mpich_libdir="$( { mpicxx.mpich -show 2>/dev/null || mpicxx -show 2>/dev/null || true; } \
  | tr ' ' '\n' | sed -n 's/^-L//p' | grep -i mpich | head -1 )"
[ -n "$_mpich_libdir" ] && export LD_LIBRARY_PATH="$_mpich_libdir:${LD_LIBRARY_PATH}"

# Constrain UCX to shared-memory/self/tcp for every test below. Ubuntu's apt
# MPICH is built on the ch4:ucx netmod, and UCX otherwise probes InfiniBand
# verbs at MPI_Init and aborts on a runner with no RDMA hardware
# (ibv_create_srq: Operation not supported), before any test can run. The
# unit/validate tests carry no per-test ENVIRONMENT, so the pin is set here in
# the ambient environment; mpiexec forwards it to the ranks. This is the
# ch4:ucx analogue of the FI_PROVIDER=tcp pin the regression tests already set
# in tests/CMakeLists.txt (harmless on a ch4:ofi netmod).
export UCX_TLS="${UCX_TLS:-tcp,self,sm}"

# Pin single-threaded execution for this per-PR gate. The unit and validate
# ctest entries carry no per-test OMP pin (only the mpi and regression entries
# do, via ${_mpi_test_env}/${_regress_env} in tests/CMakeLists.txt); with
# OMP_NUM_THREADS unset, Kokkos grabs every core and the strict bit-exact /
# conservation unit tests run threaded. Those have no fixed reduction or
# wet/dry-threshold evaluation order across a threaded partition, so they drift
# and fail intermittently (SweModule.OutflowConditionDrainsASlopedChannel...).
# Every per-PR lane runs at one thread by design (report-P1.md; the invariant
# stated in run_openmp_lane.sh) — threaded coverage is the openmp lane's
# tolerance-based gates, not the bit-exact unit suite.
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false

cmake -B "$BUILD" -S "$ROOT" -DFREHG_WERROR=ON
cmake --build "$BUILD" -j "$(getconf _NPROCESSORS_ONLN)"

# --- MPI launch diagnostics (singleton triage for mpi.core.n4) --------------
# The mpi tests intermittently degrade to size-1 MPI_COMM_WORLD on every rank
# ("1 MPI rank" banner, [frehg:0] prefix on all output) -- four singletons then
# race on the parallel-HDF5 filename and the bit-exact layout check flakes. That
# happens only when the frehg binary's runtime libmpi does not match the MPICH
# that mpiexec launches under. These lines make the cause readable from the log
# instead of guessed: the launcher identity, the exact libmpi the binary
# resolves (DT_RPATH first, then LD_LIBRARY_PATH), and whether the dependency
# prefix ever grew its own MPI. All are read-only and cannot fail the gate.
echo "==== MPI launch diagnostics ===="
echo "-- launcher:"; command -v mpiexec mpiexec.mpich 2>/dev/null || true
mpiexec --version 2>&1 | head -3 || true
echo "-- frehg runtime MPI/PETSc linkage (the libmpi actually loaded):"
ldd "$BUILD/src/frehg" 2>/dev/null | grep -iE 'mpi|petsc' || echo "  (none reported)"
echo "-- dependency-prefix MPI libraries (MUST be empty; a hit here is the bug):"
for _p in ${CMAKE_PREFIX_PATH//:/ }; do
  ls -1 "$_p"/lib*/libmpi* "$_p"/lib*/libpmi* "$_p"/bin/mpiexec 2>/dev/null || true
done
echo "================================"

ctest --test-dir "$BUILD" -L unit --output-on-failure
# -LE regression: the rank-invariance regressions also carry the mpi label
# and are picked up (once) by the regression step below.
ctest --test-dir "$BUILD" -L mpi -LE regression --output-on-failure
# Fast regression cases gate b1-b4 from P1 on, plus the per-PR-sized
# coupled gates from P3 (restart determinism and the amendment-A14
# rank-invariance lanes) and the P4 transport restart gate. The four b5
# envelope runs and the two b6 variants carry the regression_nightly label
# (plan §11.2 keeps b5/b6 out of the per-PR set), as does the
# recorded-not-gated b6_gw_smoke. All of these compare against the legacy
# goldens, which are a development-side archive not distributed with the
# repository — skip them loudly when the archive is absent.
LEGACY="$(sed -n 's/^FREHG_LEGACY_BENCHMARKS:[^=]*=//p' "$BUILD/CMakeCache.txt")"
if [ -d "$LEGACY" ]; then
  ctest --test-dir "$BUILD" -L regression \
    -R "b1|b2|b3|b4|b5_restart|b5_rank_invariance|b6_restart" --output-on-failure \
    --no-tests=ignore
else
  echo "ci_build_and_test.sh: legacy goldens not found at '$LEGACY'"
  echo "  (FREHG_LEGACY_BENCHMARKS; development-side archive, not distributed)."
  echo "  Running only the self-contained regression gates (restart determinism"
  echo "  and rank invariance); the golden-comparison gates are SKIPPED."
  ctest --test-dir "$BUILD" -L '^regression$' \
    -R "restart|rank_invariance" --output-on-failure --no-tests=ignore
fi

"$ROOT/scripts/check_forbidden.sh"
python3 "$ROOT/scripts/check_parameter_docs.py"
(cd "$ROOT" && doxygen docs/Doxyfile)
# The documentation site must build warning-free (plan §10 P5); mkdocs and
# mkdocs-material come from pip (see .github/workflows/build.yml).
(cd "$ROOT" && python3 -m mkdocs build --strict --site-dir "$BUILD/site")
echo "ci_build_and_test.sh: all gates green"
