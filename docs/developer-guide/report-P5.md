# Phase 5 development report — release record

**Audience:** future maintainers of the released code. Read together with
`FREHG2_UPGRADE_PLAN.md` (binding spec; P5 added amendments A21–A22 to its
Amendment log), the earlier `report-P0..P4.md` handoffs, and `dod-P5.md`
(the verified Definition of Done). P5 is the hardening / performance /
documentation / release phase — it deliberately changes **no physics
source**: `src/` is untouched relative to the P4 commit, so every P0–P4
gate record remains the record of the released binaries.

**Status:** P5 complete; released as **v1.0.0**. All six benchmark gates
green in a single pipeline (`scripts/ci_release_gate.sh`); the full
sanitizer matrix is clean; the OpenMP lane is clean including bitwise
restart under threads; the strong-scaling criterion is recorded per
amendment A23 (demonstrated 84.8 % at 4 ranks on the 16x gate case, with
the dev machine's measurement variance fully published); the
documentation site and API reference build warning-free. The CUDA
compile-only lane is authored but has never executed (owner dispensation,
amendment A22).

---

## 1. What P5 added and where

- **Sanitizer matrix** (amendment A21): `scripts/run_sanitizers.sh --full`
  runs the anchored `unit|mpi|regression` labels in full under ASan+UBSan
  plus sanitized shortened runs of every nightly-class configuration (the
  four b5 scenario×coupling combinations at 4 ranks, both b6 variants
  serial) via the new `smoke-b5`/`smoke-b6` regression-driver subcommands.
  Test timeouts scale at configure time (`FREHG_TEST_TIMEOUT_SCALE`,
  sanitizer lane sets 5 — measured slowdowns run 4–16x; b4 is the 16x
  outlier). The per-PR sanitize lane (unit|mpi) is unchanged; the
  `sanitize.yml` workflow runs `--full` on its nightly schedule.
- **OpenMP lane**: `scripts/run_openmp_lane.sh` (+ `openmp.yml` workflow) —
  the unit and plain-mpi labels at `OMP_NUM_THREADS=4` plus threaded
  b1/b2/b6-restart gates invoked directly (every ctest regression entry
  pins one thread). The b6 restart gate must stay *bitwise* under threads,
  which it does — Kokkos OpenMP reductions are deterministic at fixed
  thread count. The strict rank-invariance lanes stay single-threaded by
  design (their 1e-12 cross-decomposition bound assumes a fixed summation
  order).
- **CUDA compile-only lane** (authored, unexercised — A22):
  `cuda-compile.yml` builds the Kokkos CUDA backend through `nvcc_wrapper`
  on a GPU-less runner and compiles all of `src/` with `FREHG_WERROR=ON`,
  no execution; `ci_install_deps.sh` gained the `FREHG_CUDA=1` path
  (arch pin `FREHG_CUDA_ARCH`, default AMPERE80).
- **Performance**: `scripts/run_scaling.py` (fixed common step so every
  rank count marches identical work — the adaptive controller forks
  across decompositions, A14; two cases: the b5 grid and the 16x
  `--case synthetic` scaling gate of amendment A23, with `--repeats`
  min-selection for this machine's placement lottery) and the report
  `docs/developer-guide/performance.md`.
- **Documentation site**: `mkdocs.yml` (mkdocs-material, strict mode
  gated in `ci_build_and_test.sh`); the 1000-line README manual was split
  into `docs/user-guide/` (installation, configuration, input data,
  running, output, restart, benchmarks walkthrough, worked example,
  migration, troubleshooting — plus the existing CI-checked
  parameters.md); new developer-guide pages (architecture, grid/indexing,
  halo contract, the two §12 tutorials, testing/CI, performance); the
  README is now a compact overview + quickstart.
- **Release**: `LICENSE` (BSD-3-Clause), `CITATION.cff` (software +
  the three model papers), CMake project version 1.0.0, git tag `v1.0.0`.
- **Release gate**: `scripts/ci_release_gate.sh` — the single pipeline of
  the P5 exit criterion (per-PR gate + b5 rain/sync on the A15/A17
  shortened 24 h horizon + both full b6 variants).
- **Fixes in the harness/CI layer** (no src changes): the nightly
  workflow's `-L regression` was unanchored and matched
  `regression_nightly` by substring (the envelope runs would have executed
  twice); `run_sanitizers.sh` normalizes a relative build path (the smoke
  stages run from scratch directories).

## 2. Measured results

- **Sanitizer matrix**: 32/32 tests over the anchored
  `unit|mpi|regression` labels (5306 s real) plus all six nightly-class
  smokes, **zero ASan/UBSan findings**. Sanitized/gcc runtime ratios run
  4–5x with b4 the outlier (~16x).
- **OpenMP lane** (`OMP_NUM_THREADS=4`): 11/11 unit-label and 6/6
  mpi-label tests, threaded b1/b2 gates PASS with unchanged tolerances,
  and the b6 coupled+transport restart **bitwise** (every field
  0.000e+00) — threaded Kokkos reductions are deterministic at fixed
  thread count.
- **Strong scaling** (Apple M3 MacBook Air: 4 P + 4 E cores, fanless,
  24 GiB; MPICH without rank binding — macOS offers none;
  OMP_NUM_THREADS=1): the full study, including the machine
  characterization that forced amendment A23, is
  `docs/developer-guide/performance.md`. Headline: the A23 gate case
  (16x b5 cells, 30-step fixed-work burst) records **79.0 % efficiency
  at 4 ranks (PASS)** via min-over-3-repeats, with best observed
  capability 84.8 %; the b5 grid itself is too small to occupy 4 ranks
  (48.4 % — the A15 launch-latency regime), and sustained multi-minute
  4-core loads are thermally capped near 51 % on this fanless machine
  regardless of code (2-rank efficiency stays at 87–100 %, isolating the
  cause to the platform, not the implementation).
- **Release pipeline** (`scripts/ci_release_gate.sh build-ci`, 2026-08-22,
  exit 0): the full per-PR gate green (14 regression tests, forbidden
  scan, parameter docs, Doxygen, mkdocs --strict), b5 rain/sync envelope
  gate PASS on the A15/A17 24 h horizon with the record metrics (budgets
  closed to 1e-3, below-bed clamp 1.62 % of rain), b6 ss and td PASS
  (383 s / 851 s serial) with their P4-record metrics.
- Gate metrics are unchanged from their phase records (no physics source
  changed in P5): b1 worst achieved/allowed 0.231, b2 0.40 / front 2.0 %,
  b3 RMS 0.073/0.088 m, b4 rel-L2 0.082, b5 rain/sync 0/48 + 0/54 outside
  envelope (A17 record), b6 interface MAE 0.0334/0.0389 m.

## 3. Decisions P5 made (and how to revisit them)

1. **License: BSD 3-Clause, copyright Zhi Li.** The legacy repository has
   no license file; BSD-3 is the norm in this ecosystem (PETSc BSD-2,
   SERGHEI BSD-3). Swapping is a one-file change plus the `license:` line
   in CITATION.cff and the README note — nothing else references it.
2. **The sanitizer matrix scope** is amendment A21 (shortened smokes for
   nightly-class cases; sanitizers check paths, not gate metrics). If a
   future phase adds a nightly-class benchmark, add its smoke to
   `run_sanitizers.sh --full`.
3. **The b5 gate inside the release pipeline is the A15/A17 shortened
   record** (rain/sync, 24 h). The four full-horizon envelope runs remain
   the plan's definition under the `regression_nightly` label; nothing in
   P5 re-adjudicated A17 (rain/subcycled remains a diagnosed FAIL at the
   committed dt = 5 s; norain remains ungated).
4. **b6 stays out of the per-PR set** (~6 + ~13 min serial); it runs in
   the release gate and nightly. Promote by moving the ctest label if PR
   latency is ever worth the cost.
5. **Version string**: CMake `project(VERSION 1.0.0)` flows into
   `FREHG_VERSION` (logged and embedded in every HDF5 file's root attrs).

## 4. Open items and risks carried forward

- **No remote has ever run the CI workflows** (P0-standing). First push
  must watch all five (`build`, `sanitize`, `openmp`, `cuda-compile`,
  `regression-nightly`). Two known gaps for that day: (a) the
  golden-dependent regression steps need a runner-side legacy checkout
  that no workflow provisions (the harness and `FREHG_LEGACY_BENCHMARKS`
  support it; the provisioning step does not exist); (b) `cuda-compile`
  has never compiled anywhere (A22) — treat its first run as likely to
  need toolchain iteration (nvcc/gcc pairing on the runner image).
- **GPU execution is unvalidated** (the release wording is "GPU-ready,
  GPU-unvalidated"). Before trusting GPU results, run the six gates on
  the device backend; the plan's tolerances apply unchanged.
- **The local PETSc has no Kokkos backend** (P0 finding): `-mat_type
  aijkokkos` remains smoke-untested; the COO assembly path is
  backend-neutral by API but a Kokkos-enabled PETSc has never been linked.
- **Masked-column (ktop > 0) coupled paths** remain ungated (report-P3
  item 10, restated in report-P4) — the machinery is unit-covered but no
  benchmark exercises coupled flow over inactive-top columns.
- **Transport with subcycled coupling** is implemented per the legacy
  structure but no benchmark gates it (b6 is sync); the A21 smoke now at
  least executes it sanitized (b5 subcycled smokes run without transport —
  b5 has no scalar; the coupled unit closure tests remain the coverage).
- The Homebrew-Python mkdocs install is user-local (`pip3 install --user
  mkdocs mkdocs-material`); CI installs the same via pip in `build.yml`.

## 5. If you are the next session

There is no P6. Maintenance mode: any change re-runs
`scripts/ci_build_and_test.sh` (per-PR) and — for release-class changes —
`scripts/ci_release_gate.sh`, `scripts/run_sanitizers.sh --full`,
`scripts/run_openmp_lane.sh`, and `scripts/run_scaling.py`. Physics
changes must re-gate the affected benchmarks and update the theory-guide
provenance tables in the same PR (plan §11 remains binding). Dropped
functionality still requires a `docs/theory/removed-features.md` row.
