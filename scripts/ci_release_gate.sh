#!/usr/bin/env bash
# ci_release_gate.sh — the single pipeline verifying every release-blocking
# benchmark gate (plan §10 P5: "all six benchmark gates green in a single
# CI pipeline"), sized for an 8-core workstation:
#
#   1. scripts/ci_build_and_test.sh — strict zero-warning build, unit + mpi
#      labels, the b1/b2/b3/b4 gates, the restart-determinism and
#      rank-invariance regressions, forbidden-pattern scan, parameter-docs
#      lockstep, Doxygen (undocumented public API = error), and the mkdocs
#      site with --strict.
#   2. the b5 envelope gate, rain/sync on the amendment-A15/A17 shortened
#      24 h horizon at 4 ranks — the gated record P3 shipped with (the four
#      full-horizon runs remain the nightly definition, ctest label
#      regression_nightly).
#   3. both b6 variants in full (the plan §9 b6 gate, ~6 + ~13 min serial).
#
# Same environment expectations as ci_build_and_test.sh (CXX must be the
# real compiler, CMAKE_PREFIX_PATH the dependency prefix).

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-$ROOT/build-ci}"

"$ROOT/scripts/ci_build_and_test.sh" "$BUILD"

export OMP_NUM_THREADS=1 OMP_PROC_BIND=false FI_PROVIDER=tcp
MPIEXEC="$(sed -n 's/^MPIEXEC_EXECUTABLE:[^=]*=//p' "$BUILD/CMakeCache.txt")"
LEGACY="$(sed -n 's/^FREHG_LEGACY_BENCHMARKS:[^=]*=//p' "$BUILD/CMakeCache.txt")"

python3 "$ROOT/tests/regression/run_regression.py" b5 \
  --scenario rain --coupling sync --t-end 86400 --ranks 4 \
  --frehg "$BUILD/src/frehg" --mpiexec "$MPIEXEC" \
  --repo "$ROOT" --legacy "$LEGACY" \
  --work "$BUILD/release-gate/b5-rain-sync"

ctest --test-dir "$BUILD" -R 'regression\.b6\.(ss|td)' --output-on-failure

echo "ci_release_gate.sh: all six benchmark gates green"
