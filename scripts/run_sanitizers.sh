#!/usr/bin/env bash
# run_sanitizers.sh — configure, build, and test Frehg2 under ASan + UBSan
# (plan §8.4, §10 P5).
#
# Usage: scripts/run_sanitizers.sh [--full] [build-dir] [prefix-path]
#
# Default (per-PR lane): the unit and mpi ctest labels — which include the
# rank-invariance regressions through their mpi label.
# --full (the P5 sanitizer matrix, amendment A21): every per-PR label in
# full (unit, mpi, regression) plus sanitized path-coverage runs of the
# nightly-class cases on shortened horizons (the four b5 scenario x
# coupling combinations at 4 ranks and both b6 variants serial, via the
# regression driver's smoke-b5/smoke-b6 subcommands). The full-horizon b5
# envelope runs are hours long uninstrumented and are not sanitizable in
# practice; the shortened runs execute the same code paths, which is what
# the sanitizers check. Timeouts are scaled at configure time
# (FREHG_TEST_TIMEOUT_SCALE) because instrumented runs are 3-10x slower.

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

FULL=0
if [ "${1:-}" = "--full" ]; then
  FULL=1
  shift
fi
BUILD="${1:-$ROOT/build-asan}"
PREFIX="${2:-${CMAKE_PREFIX_PATH:-}}"
if [ -z "$PREFIX" ]; then
  echo "run_sanitizers.sh: pass the dependency prefix as the second argument" >&2
  echo "  or set CMAKE_PREFIX_PATH (e.g. \$HOME/frehg-deps)." >&2
  exit 2
fi
# The smoke runs below execute from scratch directories; the build path
# must survive that cwd change.
case "$BUILD" in
  /*) ;;
  *) BUILD="$PWD/$BUILD" ;;
esac

cmake -B "$BUILD" -S "$ROOT" \
  -DCMAKE_PREFIX_PATH="$PREFIX" \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DFREHG_SANITIZE=ON \
  -DFREHG_WERROR=ON \
  -DFREHG_TEST_TIMEOUT_SCALE=5 \
  "${FREHG_CMAKE_EXTRA[@]:-}"
cmake --build "$BUILD" -j "$(getconf _NPROCESSORS_ONLN)"

# halt_on_error: any finding fails the lane. Leak checking is enabled where
# the platform ASan supports it (Linux CI); macOS ASan has no leak detector.
export ASAN_OPTIONS="halt_on_error=1:${ASAN_OPTIONS:-}"
export UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1:${UBSAN_OPTIONS:-}"
if [ "$(uname -s)" = "Linux" ]; then
  export ASAN_OPTIONS="detect_leaks=1:$ASAN_OPTIONS"
fi

if [ "$FULL" -eq 0 ]; then
  ctest --test-dir "$BUILD" -L 'unit|mpi' --output-on-failure
  echo "run_sanitizers.sh: clean"
  exit 0
fi

# ---- P5 full matrix (amendment A21) ---------------------------------------
# Anchored label regex: plain 'regression' would also match the
# regression_nightly label by substring. The golden-comparison gates need
# the legacy archive (not distributed); without it, run the self-contained
# subset (restart determinism, rank invariance) and say so.
LEGACY="$(sed -n 's/^FREHG_LEGACY_BENCHMARKS:[^=]*=//p' "$BUILD/CMakeCache.txt")"
if [ -d "$LEGACY" ]; then
  ctest --test-dir "$BUILD" -L '^(unit|mpi|regression)$' --output-on-failure \
    --no-tests=error
else
  echo "run_sanitizers.sh: legacy goldens not found at '$LEGACY' —"
  echo "  golden-comparison gates SKIPPED; running the self-contained set."
  ctest --test-dir "$BUILD" -L '^(unit|mpi)$' --output-on-failure --no-tests=error
  ctest --test-dir "$BUILD" -L '^regression$' -R "restart|rank_invariance" \
    --output-on-failure --no-tests=ignore
fi

# Nightly-class path coverage on shortened horizons. Same environment pins
# as the ctest regression entries (tests/CMakeLists.txt).
export OMP_NUM_THREADS=1 OMP_PROC_BIND=false FI_PROVIDER=tcp
MPIEXEC="$(sed -n 's/^MPIEXEC_EXECUTABLE:[^=]*=//p' "$BUILD/CMakeCache.txt")"
LEGACY="$(sed -n 's/^FREHG_LEGACY_BENCHMARKS:[^=]*=//p' "$BUILD/CMakeCache.txt")"
SMOKE=(python3 "$ROOT/tests/regression/run_regression.py")
SMOKE_ARGS=(--frehg "$BUILD/src/frehg" --mpiexec "$MPIEXEC" --repo "$ROOT" --legacy "$LEGACY")

for scenario in rain norain; do
  for coupling in sync subcycled; do
    "${SMOKE[@]}" smoke-b5 --scenario "$scenario" --coupling "$coupling" \
      --ranks 4 "${SMOKE_ARGS[@]}" \
      --work "$BUILD/sanitize-smoke/b5-$scenario-$coupling"
  done
done
for variant in ss td; do
  "${SMOKE[@]}" smoke-b6 --variant "$variant" "${SMOKE_ARGS[@]}" \
    --work "$BUILD/sanitize-smoke/b6-$variant"
done

echo "run_sanitizers.sh --full: clean"
