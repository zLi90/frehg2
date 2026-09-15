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

cmake -B "$BUILD" -S "$ROOT" -DFREHG_WERROR=ON
cmake --build "$BUILD" -j "$(getconf _NPROCESSORS_ONLN)"

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
