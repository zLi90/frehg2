#!/usr/bin/env bash
# ci_release_gate.sh — the single pipeline verifying every release-blocking
# gate. v1.0.0 scope was the six benchmark gates (plan §10 P5); the v2.0.0
# scope (v2 plan §9 Q7) is **all g + all b + s1–s4 + p1–p5 in one pipeline**:
#
#   1. scripts/ci_build_and_test.sh — strict zero-warning build, unit + mpi
#      labels (the unit tier carries the §8.2 BC-matrix and §8.3
#      feature-coverage lockstep checkers, the p5 negative compile test, and
#      the parameter-docs lockstep), the b1/b2/b3/b4 gates, the
#      restart-determinism and rank-invariance regressions, forbidden-pattern
#      scan, Doxygen (undocumented public API = error), and the mkdocs site
#      with --strict.
#   2. the FULL per-PR regression label: g1's fast lanes (amg/gamg on b1–b4),
#      g4, g6, g7, g9, g10, the p1 backend-invariance lanes, r2, the §8.1
#      dihedral batteries (swe/gw/transport/wind/heat), the §8.2 side sweeps
#      (outflow staircase ×4 edges, surface_bc_sides), the §8.3 pair gates
#      and §8.4 restart-mid-transient rows, and the perf baseline.
#   3. the b5 envelope gate, rain/sync on the amendment-A15/A17 shortened
#      24 h horizon at 4 ranks — the adjudicated record P3 shipped with
#      (V2-A4: the full horizon was never adjudicated as a gate).
#   4. the nightly regression label's adjudicated gates: both b6 variants,
#      g5 (Geng 2015 + its transposed slice), g8, the g1 nightly amg lanes
#      (b5@86400 + b6 + rank invariance), and the aijkokkos full-physics
#      pair lane. The four full-horizon b5 envelope runs are excluded —
#      never adjudicated as gates (A15/A17, V2-A4, V2-A20); they stay in
#      the regression_nightly label as on-demand local measurements. This
#      step is still hours (g5, b6 full).
#   5. the nightly scaling label: s1–s4 (s2 doubles as g3 weak scaling), the
#      g2 iteration-flatness pair (amg + the bjacobi record), p2 solver
#      thread scaling, p3 hybrid placements. Run on an otherwise idle
#      machine — the s-gates assert timing at the largest rank count that
#      leaves a spare performance core (V2-A12).
#   6. p4 (CUDA compile+link): built here when nvcc is available; otherwise
#      the evidence is the cuda-compile.yml CI run on the release SHA, named
#      via FREHG_P4_EVIDENCE (recorded into this log; the run happens on the
#      official remote after the phase-boundary sync, before tagging).
#
# The sanitizer matrix (scripts/run_sanitizers.sh --full) is a separate
# release requirement with its own toolchain (plan §9 Q7 "also required");
# run it independently and record it in the DoD.
#
# Same environment expectations as ci_build_and_test.sh (CXX must be the
# real compiler, CMAKE_PREFIX_PATH the dependency prefix). The legacy
# goldens archive is REQUIRED here — a release cannot gate without it.

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:-$ROOT/build-ci}"
# The b5 step below launches from a staged scratch directory, so the build
# path must survive that cwd change (same rule as run_sanitizers.sh).
case "$BUILD" in
  /*) ;;
  *) BUILD="$PWD/$BUILD" ;;
esac

"$ROOT/scripts/ci_build_and_test.sh" "$BUILD"

export OMP_NUM_THREADS=1 OMP_PROC_BIND=false FI_PROVIDER=tcp UCX_TLS=tcp,self,sm
MPIEXEC="$(sed -n 's/^MPIEXEC_EXECUTABLE:[^=]*=//p' "$BUILD/CMakeCache.txt")"
LEGACY="$(sed -n 's/^FREHG_LEGACY_BENCHMARKS:[^=]*=//p' "$BUILD/CMakeCache.txt")"
if [ ! -d "$LEGACY" ]; then
  echo "ci_release_gate.sh: legacy goldens REQUIRED for a release gate" >&2
  echo "  (FREHG_LEGACY_BENCHMARKS='$LEGACY' not found)" >&2
  exit 2
fi

# 2. The full per-PR regression tier (anchored: the plain regex 'regression'
# would also select regression_nightly by substring).
ctest --test-dir "$BUILD" -L '^regression$' --output-on-failure

# 3. The adjudicated b5 record (rain/sync, 24 h horizon, 4 ranks).
python3 "$ROOT/tests/regression/run_regression.py" b5 \
  --scenario rain --coupling sync --t-end 86400 --ranks 4 \
  --frehg "$BUILD/src/frehg" --mpiexec "$MPIEXEC" \
  --repo "$ROOT" --legacy "$LEGACY" \
  --work "$BUILD/release-gate/b5-rain-sync"

# 4. The nightly regression tier's ADJUDICATED gates: b6 both variants, g5
# (Geng 2015 + transposed slice), g8, the g1 amg nightly lanes, the
# aijkokkos full-physics pair lane. The four FULL-HORIZON b5 envelope runs
# (regression.b5.{rain,norain}.{sync,subcycled}) are deliberately excluded:
# their full-horizon bounds were never adjudicated as gates (A15/A17 — the
# owner concluded the b5 battery early with rain/sync adjudicated at
# --t-end 86400, which step 3 above asserts; V2-A4 reaffirmed that record;
# norain was stopped mid-flight at P3 and first completed 2026-09-26;
# rain/subcycled is documented-failing at the committed dt, A17). They
# stay in the regression_nightly label as on-demand local measurements:
# ~8 h each at 4 ranks, beyond a hosted runner's 6 h job cap, and the
# hosted nightly has never executed them (goldens not provisioned there —
# V2-A20, which also records the release-time measurements). Asserting
# never-adjudicated bounds in the release pipeline would gate on noise;
# hiding the exclusion would be the vacuous-gate pattern — hence this
# block and the amendment.
ctest --test-dir "$BUILD" -L regression_nightly \
  -E '^regression\.b5\.(rain|norain)\.' --output-on-failure
echo "note: the four full-horizon b5 envelope runs are excluded here by"
echo "      design (never-adjudicated bounds; see V2-A20 and the comment"
echo "      above). The adjudicated b5 record is step 3; the full-horizon"
echo "      lanes stay runnable on demand: ctest -R '^regression\\.b5\\.'"

# 5. The nightly scaling tier (s1-s4/g2/g3, p2, p3).
ctest --test-dir "$BUILD" -L scaling_nightly --output-on-failure

# 6. p4 — CUDA compile+link.
if command -v nvcc >/dev/null 2>&1; then
  cmake -B "$BUILD-cuda" -S "$ROOT" -DFREHG_WERROR=ON \
    -DKokkos_ENABLE_CUDA=ON
  cmake --build "$BUILD-cuda" -j "$(getconf _NPROCESSORS_ONLN)"
  echo "p4: CUDA lane compiled and linked locally"
elif [ -n "${FREHG_P4_EVIDENCE:-}" ]; then
  echo "p4: no nvcc on this machine; evidence recorded:"
  echo "    $FREHG_P4_EVIDENCE"
  echo "    (cuda-compile.yml must be green on the pushed release SHA"
  echo "     before the official tag — verify after the §1.2 sync.)"
else
  echo "ci_release_gate.sh: p4 requires nvcc or FREHG_P4_EVIDENCE" >&2
  echo "  (set FREHG_P4_EVIDENCE to the cuda-compile.yml run that covers" >&2
  echo "   this tree, e.g. a run id + the statement that the workflow and" >&2
  echo "   device-facing sources are unchanged since.)" >&2
  exit 2
fi

echo "ci_release_gate.sh: every release-blocking gate green"
echo "  (b1-b6, g1-g10, s1-s4, p1-p3+p5 executed here; p4 per the evidence"
echo "   above; sanitizer matrix and the owner p6 bundle tracked in dod-Q7)"
