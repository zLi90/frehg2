# Testing and CI

The test architecture is plan §8; this page is the practical guide.

## Test labels

| Label | Contents | When it runs |
|---|---|---|
| `unit` | the GoogleTest suite (`frehg_unit_tests`), the migration-tool test, the schema↔docs lockstep check, the §8.2 BC-matrix and §8.3 feature-coverage lockstep checkers, the run-record checker self-test, `frehg --validate` on all benchmark configs (incl. the wind and x-orient bases) | every PR |
| `mpi` | the halo and parallel-core drivers at 1/2/3/4 ranks (n=3 is the §8.4 non-divisible split; `FREHG_EXPECT_RANKS` asserts the launched world size), plus the rank-invariance regressions (dual-labeled `regression;mpi`) | every PR |
| `regression` | the b1–b4 gates; the restart-determinism gates (b1/b2/b5/b6); the rank-invariance lanes (b1/b2/b5/tracer/heat × strict/default); the recorded b6 groundwater smoke; the fast v2 gates (g4, g6, g7, g9, g10, the g1 amg/gamg b1–b4 lanes, the p1 backend-invariance lanes, r2); the §8.1 dihedral batteries (swe/gw/transport/wind/heat); the §8.2 side sweeps (`outflow_staircase` ×4 edges, `surface_bc_sides`); the §8.3 pair gates and §8.4 restart-mid-transient rows; the perf baseline | every PR (fast subset via `ci_build_and_test.sh`), fully in the nightly |
| `regression_nightly` | the four full-horizon b5 envelope runs (wall-hours each), the two full b6 variants, g5 (Geng 2015 + transposed slice), g8, the g1 nightly amg lanes (b5/b6 + rank invariance), the aijkokkos full-physics pair lane | nightly / deliberately |
| `scaling_nightly` | s1–s4 (s2 doubles as g3 weak scaling), the g2 iteration-flatness pair (amg + the bjacobi record), p2 solver-thread scaling, p3 hybrid placements | nightly, on an otherwise idle machine (V2-A12: timing gates assert at the largest rank count leaving a spare performance core) |

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

Two lanes per case (plan §8.2, amendments A5/A8/A14, V2-A9):

- **strict** — rank-invariant `jacobi` preconditioning at machine-precision
  tolerances proves the assembly, physics, and exchange operators
  rank-invariant (b1/b2 bounds at 1e-12, measured down to ~1e-14). Horizons
  are limited to windows before discrete threshold crossings (wet/dry fronts,
  saturation-front cell crossings) amplify rounding differences — but *which
  step* the first crossing lands on is platform rounding, so the coupled b5
  one-step bound is 1e-9 (V2-A9: one threshold-class branch flips inside
  step 1 on the x86-64 CI platform, isolated footprint 2.2e-10; the harness
  prints the per-field breakdown for the record).
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
sanitizer runtimes, so the local lane builds with Apple clang against a
clang-built Kokkos, while PETSc comes from the main gcc-built prefix:
PETSc's C API crosses that C++-runtime boundary safely, its Kokkos
objects do not). `--full` runs the anchored
`unit|mpi|regression` labels minus two exclusions it prints by name —
`regression.perf_baseline` (a wall-clock gate, which instrumented runs
fail by construction) and, when the sanitized binary and PETSc link
different C++ runtimes, the `aijkokkos` lanes (two Kokkos builds in one
process; V2-A20) — plus sanitized shortened runs of every nightly-class
configuration (the four b5 scenario×coupling combinations and both b6
variants, `smoke-b5`/`smoke-b6` subcommands, and the v2 `smoke-g5`/
`smoke-g8` runs) — the full-horizon
envelope runs are not sanitizable in practice (hours uninstrumented ×
3–10 slowdown), and path coverage, not gate metrics, is what sanitizers
check (amendment A21). Test timeouts are scaled at configure time
(`FREHG_TEST_TIMEOUT_SCALE=5`).

Linux leak checking runs with the third-party suppressions in
`scripts/lsan.supp`: MPICH's singleton-mode hwloc topology and PETSc's
registry state are process-lifetime allocations their finalizers never
free, and without the suppressions every direct (non-`mpiexec`) binary
invocation — `frehg --validate`, the r1 run-record revalidation, the unit
gtest binary — exits nonzero on a constant ~22 KB report. The file
documents the per-module rationale and the one accepted blind spot
(leaks of PETSc-held objects match the libpetsc rule; frehg's own
allocations stay fully covered).

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

## Generality batteries and coverage lockstep (v2 plan §8, Q7)

Accuracy gates fix one orientation, one BC side, one decomposition, one
feature combination; the §8 instruments close what they cannot see:

- **§8.1 dihedral batteries** (`regression.{swe,gw,transport}_orient`,
  plus the Q5/Q6 `heat_orient`/`wind_orient`): each base case in
  `benchmarks/x-orient/` runs through the 8 dihedral transforms
  (bathymetry, ICs, BC polygons, rasters, and forcing directions
  transformed together) and the outputs must inverse-map onto the identity
  run. Two tolerance classes per battery
  (`tests/regression/tolerances/x-orient.yaml`): `transpose` preserves
  each edge's minus/plus class and gates at the strict floor; the six
  mixed transforms carry the measured legacy plus/minus edge-arithmetic
  floor. The waiver mechanism is the
  [symmetry-exemption table](../theory/symmetry-exemptions.md) — an
  undocumented over-tolerance deviation is a bug, and new v2 physics gets
  no exemptions.
- **§8.2 BC matrix** (`unit.bc_matrix`,
  `tests/coverage/bc_matrix.csv`, `scripts/check_bc_matrix.py`): the full
  kind × target × coupling enumeration; every cell is `tested` (ids
  cross-referenced against the gtest/ctest/validation universes) or
  `schema-rejected` (with a test pinning the message); `limitation` rows
  fail at release. The kind/target vocabularies are re-read from
  `ConfigSchema.cpp`, so a new BC kind cannot land without extending the
  matrix. The side sweeps behind it: `regression.outflow_staircase` (the
  V2-A11 cells — a transmissive outlet on **every** edge, both feeding
  directions, absolute Manning/volume criteria) and
  `regression.surface_bc_sides` (the conveyor case in 4 rotations ×
  {uncoupled, coupled}).
- **§8.3 feature-interaction coverage** (`unit.feature_coverage`,
  [feature-coverage.md](feature-coverage.md),
  `scripts/check_feature_coverage.py`): zero empty cells; every cited
  test exists; GPU rows carry an explicit `experimental`/`supported`
  status. The pair gates it added: `regression.{wind,evap,heat}_restart`
  (restart mid-transient), `regression.{b6,heat}_subcycled`,
  `regression.{wind_drydown,terrain_heat,two_scalar_restart,evap_heat}`,
  and the AMG/aijkokkos + restart lanes.
- **§8.4 decomposition/edge cases**: the n=3 mpi lanes, the transposed g5
  slice (`regression.g5.transposed`), single-cell and corner BC polygons
  and the masked-column-per-edge sweep (unit tier), and the coupled ×
  masked-column loud refusal (`unit.masked_coupled_fatal`).

## The release pipeline (v2.0.0)

```bash
scripts/ci_release_gate.sh build-ci
```

The plan §9 Q7 blocking criterion — **all g + all b + s1–s4 + p1–p5 in
one pipeline**: the per-PR gate (step 1), the full anchored regression
label (step 2), the adjudicated b5 record (rain/sync, 24 h horizon, step
3), the nightly regression label's adjudicated gates (step 4 — the four
full-horizon b5 envelope runs are excluded and the exclusion is printed:
their bounds were never adjudicated as gates, V2-A20), the nightly
scaling label (step 5), and p4 (CUDA compile+link where nvcc exists;
otherwise the `FREHG_P4_EVIDENCE` string names the CI run covering the
tree). The sanitizer matrix (`run_sanitizers.sh --full`) is a separate
release requirement with its own toolchain. Runs wall-hours; use an
otherwise idle machine so the scaling gates see their calibrated
headroom.

## CI workflows (`.github/workflows/`)

| Workflow | What it runs |
|---|---|
| `build.yml` | gcc + clang matrix: the per-PR gate |
| `sanitize.yml` | the ASan/UBSan lane |
| `openmp.yml` | the threaded lane |
| `cuda-compile.yml` | Kokkos CUDA backend, compile+link, zero warnings (no GPU executes; the p4 gate) |
| `regression-nightly.yml` | full regression + nightly regression + nightly scaling labels on a schedule |

The workflows were authored and version-pinned before the project had CI
history (report-P5.md records that context). The golden-comparison
regression steps need the legacy-goldens archive (a development-side
artifact, not distributed); every entry point (`ci_build_and_test.sh`,
the nightly workflow, the full sanitizer lane) detects its absence and
falls back to the self-contained gates (restart determinism, rank
invariance), printing a notice — **but the job still reports success**.
The archive has never been provisioned on the hosted runner, so every
scheduled nightly to date has validated build, unit, mpi,
restart/rank-invariance, and the s-gates only; all golden and envelope
evidence is local (V2-A20). To arm the golden gates on a runner,
provision the archive and point `FREHG_LEGACY_BENCHMARKS` at it; the four
full-horizon b5 envelope runs (~8 h each at 4 ranks) additionally exceed
a hosted runner's 6 h job cap and need a self-hosted runner or their own
jobs.

## Forbidden patterns

`scripts/check_forbidden.sh` fails on any hit in `src/`: TODO/FIXME-class
markers, warning suppression, `file(GLOB`, raw allocation, `printf`/`cout`
outside the logger, `catch (...)`, UVM/managed-memory types, backend
`#ifdef`s in physics directories, and the dropped-feature keyword tripwire.
Run it early and often — it is cheap and it is a hard CI gate.
