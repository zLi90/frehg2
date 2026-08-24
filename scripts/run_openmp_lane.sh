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

export OMP_NUM_THREADS="$THREADS"
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

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
