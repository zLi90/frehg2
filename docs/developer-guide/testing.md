# Testing and CI

The test architecture is plan §8; this page is the practical guide.

## Test labels

| Label | Contents | When it runs |
|---|---|---|
| `unit` | the GoogleTest suite (`frehg_unit_tests`), the migration-tool test, the schema↔docs lockstep check, `frehg --validate` on all benchmark configs | every PR |
| `mpi` | the halo and parallel-core drivers at 1/2/4 ranks, plus the rank-invariance regressions (dual-labeled `regression;mpi`) | every PR |
| `regression` | the b1–b4 gates, the restart-determinism gates (b1/b2/b5/b6), the rank-invariance lanes (b1/b2/b5 × strict/default), the recorded b6 groundwater smoke | every PR (fast subset via `ci_build_and_test.sh`), fully in the nightly |
| `regression_nightly` | the four full-horizon b5 envelope runs (wall-hours each) and the two full b6 variants | nightly / deliberately |

```bash
ctest --test-dir build -L unit --output-on-failure
ctest --test-dir build -L mpi -LE regression
ctest --test-dir build -L '^regression$'          # anchored: excludes _nightly
```

Regression gates need the legacy goldens (`FREHG_LEGACY_BENCHMARKS` cache
variable, default `../legacy/benchmarks`); goldens are never committed.

## The regression harness

`tests/regression/run_regression.py` stages a benchmark into a scratch
directory, runs it, and applies the plan §9 gate, always printing
achieved-vs-allowed for every metric (so tolerances can be revisited from
data). Tolerances live in `tests/regression/tolerances/<case>.yaml`;
element-wise comparison logic in `tests/regression/compare_h5.py`;
`tools/ascii_golden_to_h5.py` converts legacy ASCII goldens on the fly.

Environment pins (set on every ctest regression entry): `OMP_NUM_THREADS=1`
(benchmark grids are launch-latency-bound and the gates stay
deterministic), `FI_PROVIDER=tcp` (the libfabric sockets provider can wedge
MPICH's `MPI_Finalize` — post-completion only, but it hangs harness runs).

## Rank invariance

Two lanes per case (plan §8.2, amendments A5/A8/A14):

- **strict** — rank-invariant `jacobi` preconditioning at machine-precision
  tolerances proves the assembly, physics, and exchange operators
  rank-invariant (bounds at 1e-12, measured down to ~1e-14). Horizons are
  limited to windows before discrete threshold crossings (wet/dry fronts,
  saturation-front cell crossings) amplify rounding differences.
- **default** — production solver settings over longer windows, gated at
  data-derived bounds (the harness prints achieved values on every run).

## Per-PR gate

```bash
scripts/ci_build_and_test.sh
```

Strict zero-warning build, unit + mpi labels, the fast regression set,
`scripts/check_forbidden.sh`, `scripts/check_parameter_docs.py`, and
Doxygen with warnings as errors (undocumented public API fails the build).

## Sanitizers

```bash
scripts/run_sanitizers.sh              # per-PR lane: unit + mpi labels
scripts/run_sanitizers.sh --full       # the P5 matrix: all per-PR labels +
                                       # shortened nightly-class smokes
```

ASan + UBSan with `halt_on_error=1` (leak detection on Linux; macOS ASan
has no leak detector — the gcc toolchain on macOS/arm64 also ships no
sanitizer runtimes, so the local lane builds with Apple clang against
clang-built dependencies). `--full` runs the anchored
`unit|mpi|regression` labels plus sanitized shortened runs of every
nightly-class configuration (the four b5 scenario×coupling combinations
and both b6 variants, `smoke-b5`/`smoke-b6` subcommands) — the full-horizon
envelope runs are not sanitizable in practice (hours uninstrumented ×
3–10 slowdown), and path coverage, not gate metrics, is what sanitizers
check (amendment A21). Test timeouts are scaled at configure time
(`FREHG_TEST_TIMEOUT_SCALE=5`).

## The OpenMP lane

```bash
scripts/run_openmp_lane.sh build 4
```

Everything else pins `OMP_NUM_THREADS=1`, so this lane is what exercises
threaded Kokkos kernels: the unit and plain-mpi labels at 4 threads plus
threaded b1/b2/b6-restart gates (the restart gate must stay **bitwise**
under threads — both runs use the same thread count, so any
scheduling-dependent nondeterminism fails it). The strict rank-invariance
lanes stay single-threaded by design: their cross-decomposition bound
assumes a fixed summation order.

## CI workflows (`.github/workflows/`)

| Workflow | What it runs |
|---|---|
| `build.yml` | gcc + clang matrix: the per-PR gate |
| `sanitize.yml` | the ASan/UBSan lane |
| `openmp.yml` | the threaded lane |
| `cuda-compile.yml` | Kokkos CUDA backend, compile-only, zero warnings (no GPU executes) |
| `regression-nightly.yml` | full regression + nightly labels on a schedule |

The workflows were authored and version-pinned before the project had CI
history (report-P5.md records that context) — **watch their first runs**.
The golden-comparison regression steps need the legacy-goldens archive (a
development-side artifact, not distributed); every entry point
(`ci_build_and_test.sh`, the nightly workflow, the full sanitizer lane)
detects its absence and falls back loudly to the self-contained gates
(restart determinism, rank invariance). To arm the golden gates on a
runner, provision the archive and point `FREHG_LEGACY_BENCHMARKS` at it.

## Forbidden patterns

`scripts/check_forbidden.sh` fails on any hit in `src/`: TODO/FIXME-class
markers, warning suppression, `file(GLOB`, raw allocation, `printf`/`cout`
outside the logger, `catch (...)`, UVM/managed-memory types, backend
`#ifdef`s in physics directories, and the dropped-feature keyword tripwire.
Run it early and often — it is cheap and it is a hard CI gate.
