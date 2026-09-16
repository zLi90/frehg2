#!/usr/bin/env bash
# run_openmp_lane.sh — the Kokkos OpenMP lane (plan §10 P5): run the suite
# with the OpenMP host backend actively threaded. Every per-PR lane pins
# OMP_NUM_THREADS=1 (benchmark-scale grids are launch-latency-bound and the
# gates stay deterministic — report-P1.md); this lane is what exercises the
# threaded kernels: parallel_for races, reduction identities, and
# thread-count-dependent rounding all surface here and nowhere else.
#
# Usage: scripts/run_openmp_lane.sh [build-dir] [threads]
#
# Runs: the unit label and the plain mpi drivers at OMP_NUM_THREADS=<threads>,
# then three threaded regression gates invoked directly (the ctest entries
# pin OMP_NUM_THREADS=1 in their ENVIRONMENT, so they cannot be reused):
#   b1          surface-water gate  (tolerance-based; absorbs threaded
#                                    reduction-order rounding)
#   b2          groundwater gate    (same)
#   b6-restart  coupled + transport bitwise restart — both runs execute at
#               the same thread count, so determinism must hold exactly
# The strict rank-invariance lanes stay OMP_NUM_THREADS=1 by design: their
# 1e-12 cross-decomposition bound assumes a fixed summation order, which a
# threaded partition of rank-dependent extents does not provide.

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-$ROOT/build}"
THREADS="${2:-4}"

# Runtime loader path for the shared dependency libraries: Kokkos is built
# shared (A8 invariant 7), and frehg pulls libkokkoscontainers.so et al. in
# transitively through PETSc's pkg-config -L flags, which CMake does not turn
# into an rpath. The ctest and regression runs below launch the frehg binary,
# so the dep lib dirs must be on LD_LIBRARY_PATH. Derived from CMAKE_PREFIX_PATH
# (the workflow step exports it); mirrors build_frehg2_slurm.sh.
_pfx="${CMAKE_PREFIX_PATH:-}"
for _p in ${_pfx//:/ }; do
  [ -n "$_p" ] && LD_LIBRARY_PATH="$_p/lib:$_p/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
done
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Prepend the system MPICH runtime dir so the frehg binary resolves libmpi to
# the same MPICH that mpiexec.mpich launches under. If a stray libmpi ever lands
# in the dependency prefix (on LD_LIBRARY_PATH above), it would shadow the
# system one and the binary would MPI_Init under a different MPICH than the
# launcher's PMI -- every rank then degrades to a size-1 MPI_COMM_WORLD and the
# parallel-HDF5 mpi tests race on a colliding filename. A no-op when the prefix
# is clean. (DT_RPATH in libpetsc still wins over LD_LIBRARY_PATH; the ci build
# lane's diagnostics surface a downloaded-and-rpath'd MPI, if any.)
_mpich_libdir="$( { mpicxx.mpich -show 2>/dev/null || mpicxx -show 2>/dev/null || true; } \
  | tr ' ' '\n' | sed -n 's/^-L//p' | grep -i mpich | head -1 )"
[ -n "$_mpich_libdir" ] && export LD_LIBRARY_PATH="$_mpich_libdir:${LD_LIBRARY_PATH}"

# Constrain UCX to shared-memory/self/tcp for the unit/mpi ctest runs below
# (they inherit this ambient env; the ctest entries carry no UCX pin). Ubuntu's
# apt MPICH uses the ch4:ucx netmod, which otherwise probes InfiniBand verbs at
# MPI_Init and aborts on a runner with no RDMA hardware. ch4:ucx analogue of the
# FI_PROVIDER=tcp pin set further down for the regression phase.
export UCX_TLS="${UCX_TLS:-tcp,self,sm}"

export OMP_NUM_THREADS="$THREADS"
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# MPI launch diagnostics (singleton triage for mpi.core.n4): print the launcher
# identity and the exact libmpi the frehg binary resolves, so a degraded size-1
# world is readable from this lane's log too. Read-only; cannot fail the lane.
# (unit.all pins OMP_NUM_THREADS=1 in its own ctest ENVIRONMENT, overriding the
# OMP_NUM_THREADS=$THREADS above -- the strict bit-exact unit tests are
# single-threaded by design; this lane's threaded coverage is its b1/b2/
# b6-restart tolerance gates below, which run at $THREADS.)
echo "==== MPI launch diagnostics (OMP_NUM_THREADS=$THREADS) ===="
mpiexec --version 2>&1 | head -3 || true
ldd "$BUILD/src/frehg" 2>/dev/null | grep -iE 'mpi|petsc' || echo "  (none reported)"
echo "==========================================================="

ctest --test-dir "$BUILD" -L unit --output-on-failure
ctest --test-dir "$BUILD" -L mpi -LE regression --output-on-failure

export FI_PROVIDER=tcp
MPIEXEC="$(sed -n 's/^MPIEXEC_EXECUTABLE:[^=]*=//p' "$BUILD/CMakeCache.txt")"
LEGACY="$(sed -n 's/^FREHG_LEGACY_BENCHMARKS:[^=]*=//p' "$BUILD/CMakeCache.txt")"
if [ ! -d "$LEGACY" ]; then
  echo "run_openmp_lane.sh: SKIPPED the threaded b1/b2/b6-restart gates —"
  echo "  legacy goldens not found at '$LEGACY' (FREHG_LEGACY_BENCHMARKS)."
  echo "run_openmp_lane.sh: unit+mpi labels clean at OMP_NUM_THREADS=$THREADS"
  exit 0
fi
REGRESS=(python3 "$ROOT/tests/regression/run_regression.py")
REGRESS_ARGS=(--frehg "$BUILD/src/frehg" --mpiexec "$MPIEXEC" --repo "$ROOT" --legacy "$LEGACY")

for gate in b1 b2 b6-restart; do
  "${REGRESS[@]}" "$gate" "${REGRESS_ARGS[@]}" \
    --work "$BUILD/openmp-lane/$gate"
done

echo "run_openmp_lane.sh: clean at OMP_NUM_THREADS=$THREADS"
