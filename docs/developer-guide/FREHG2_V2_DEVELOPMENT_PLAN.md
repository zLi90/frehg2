# FREHG2 v2.0 Development Plan — Solver Scalability, Run Provenance, Performance Portability, Evaporation, Temperature, Wind

**Status:** DRAFT for owner review
**Prepared:** 2026-09-09
**Baseline:** frehg2 v1.0.0 (all six b1–b6 gates green; see `docs/developer-guide/completion-report.md`)
**Binding conventions:** this plan inherits the process rules of `FREHG2_UPGRADE_PLAN.md`
§11 verbatim — no deferred functionality, amendment log for every deviation, DoD file per
phase, forbidden-pattern scan, zero-warning builds, per-PR re-verification of all previous
gates (b1–b6 stay blocking throughout).

---

## 1. Scope

Six capabilities, in dependency order:

| Phase | Capability | Blocking gates |
|---|---|---|
| **Q1** | Scalable AMG preconditioner (hypre BoomerAMG primary, PETSc GAMG fallback), YAML-selectable, with setup-reuse policy | g1 (solver invariance), g2 (iteration-vs-ranks), g3 (weak-scaling, nightly) |
| **Q2** | Run provenance & timing log: every simulation writes a structured run record (resolved config, modules, BCs, hierarchical timers, solver telemetry, mass-audit closure) | r1 (record present + schema-valid + round-trips, every gate case), r2 (timer coverage ≥ 90 % of the time loop) |
| **Q3** | Performance-portable parallelization: one source tree over serial / OpenMP / MPI / MPI+OpenMP / single-GPU / multi-GPU, with a device-consistent linear-algebra path (Kokkos-aware PETSc + hypre) and backend-aware AMG defaults | p1 (backend invariance), p2 (OpenMP solver thread scaling), p3 (hybrid MPI+OpenMP scaling), p4 (device build integrity), p5 (host/device pointer discipline), p6 (owner-executed GPU acceptance bundle) |
| **Q4** | Robust open-water + soil evaporation: bulk-aerodynamic (mass-transfer) formulation, potential→actual moisture limiting, per-cell forcing, salinity coupling | g4 (analytic drawdown + evaporative concentration), g5 (Geng & Boufadel 2015 bare-soil salinization) |
| **Q5** | Temperature transport (surface + subsurface) with surface heat exchange and thermal density coupling | g6 (subsurface heat analytic), g7 (surface equilibrium-temperature relaxation), g8 (optional: thermal convection onset) |
| **Q6** | Robust wind stress: selectable Cd(U10) laws with cap, per-component wind fields | g9 (wind setup analytic), g10 (Merian seiche relaxation), plus Cd unit-test battery |
| **Q7** | Hardening, cross-coupling gates, docs, v2.0.0 release | all g-gates + b1–b6 in one pipeline |

Dependency graph: **Q1 → Q2 → Q3 → {Q4 ∥ Q6} → Q5 → Q7.** Q1 first because every later gate
re-runs the solve path and benefits from the solver telemetry Q1 adds. Q5 after Q4 because
the surface heat-exchange latent term reuses Q4's bulk-aerodynamic machinery (one vapor
pressure/humidity module, two consumers). Q6 is independent of Q4/Q5 physics but its g9/g10
cases become the base configurations for the Q5 surface heat-exchange tests, so it may run
in parallel with Q4. Q2 (run log) is a one-week infrastructure sprint inserted
directly after Q1 (V2-A5): it has no physics, and doing it before the physics
phases means every later gate run is self-recording from its first execution.
Q3 (performance portability) is likewise infrastructure and is placed **before** the
physics phases for the same reason, sharpened (V2-A6): every kernel written in Q4–Q6
must be device-correct from its first PR, and retro-fitting memory-space discipline
onto three shipped physics modules is strictly more expensive than gating it up front.
Q3 also depends on Q1 — the backend-aware AMG defaults it adds are edits to Q1's
preconditioner dispatch — and on Q2, whose run record is where the per-backend timing
evidence lands.

Explicitly **out of scope** (unchanged from v1 dispositions): multiple simultaneous
scalars beyond {salinity, temperature}, vegetation transpiration/root uptake (recorded
below as the follow-on), Coriolis, radiation BCs, sediment, NetCDF/VTK. Any of these
found necessary mid-implementation requires an amendment, not silent scope creep.

### 1.1 Precondition: clear the v1 liabilities first (Q0, one week)

These are cheap and de-risk everything after:

- **Q0.1** Reconcile the dev repo with the official public repo. The official
  `~/Codes/frehg2` (github.com/zLi90/frehg2) **already contains** the local working-tree
  changes as commits ("fixed bugs for building on gpu" = the nvcc extended-lambda
  public-section refactor; "fixed bugs on lateral groundwater bc" = the hydrostatic-edge
  `value > bed` classification) — its `src/` is byte-identical to the local uncommitted
  tree. So Q0.1 is: commit the same changes in the dev repo (with the b-gate suite
  green), then assert `diff -r` cleanliness of `src/` between the two repos.
- **Q0.2** Push the dev repo to a real remote and watch the five authored CI workflows
  execute once; fix runner-side issues (the golden-dependent steps need a documented
  `FREHG_LEGACY_BENCHMARKS` provisioning step or a vendored golden subset).
- **Q0.3** Resolve the swere-superslab broken run flagged 2026-08-30 (15 m³ rain
  injected, zero ponding/outflow/seepage) — either a config defect or a genuine
  mass-balance bug; must be diagnosed before Q4 touches the rain/evap source path.
- **Q0.4** Workspace cleanup — **done 2026-09-09** (recorded here so it isn't
  re-litigated). Deleted as verified-duplicated or regenerable: all five `frehg2/build*`
  trees; `serghei-validation/` and `frehg2-validation/` (pre-merge snapshots — every
  case was merged into `frehg2/validation/` at commit 52c29ea, and the in-repo copies
  are strictly newer); `frehg2-release/` (superseded by the `v1.0.0` git tag);
  `testcases/serghei-*` (raw upstream SERGHEI downloads, re-obtainable);
  `legacy/serghei/.git` (349 MB clone history; the source stays for provenance);
  the root duplicate of `FREHG2_UPGRADE_PLAN.md` (the tracked in-repo copy under
  `docs/developer-guide/` is canonical); all `.DS_Store`. Reference PDFs consolidated
  into `references/`. **Kept deliberately:** `legacy/benchmarks/` (the b1–b6 goldens
  the regression gates require), `legacy/frehg1.0/` and `legacy/serghei/` sources
  (line-level provenance citations throughout `src/`), `legacy/writing-samples/`,
  `deps-clang/` (the clang sanitizer lane's Kokkos/yaml-cpp builds), and
  `frehg2/validation/*/out*` (~1.9 GB of *active* run outputs feeding the manuscript
  plots — regenerable but expensive; prune per-case only after the paper's figures
  are final).

### 1.2 Repository topology (binding)

Development happens in `~/frehg2-upgrade/frehg2` (the full-history dev repo:
P0–P5 commits, tags, goldens machinery beside it). The **official public repo** is
`~/Codes/frehg2` → `github.com/zLi90/frehg2` (squashed history). Release flow: a
version is cut in the dev repo (tag), then synced to the official repo as a single
release commit; after every sync, `diff -r --exclude=.git` between the two trees
must be empty. The other `~/Codes/frehg2-*` directories (dev/legacy/testing) are
outside this plan's scope and are not written to.

### 1.3 Gate taxonomy: three axes (new in v2)

The v1 gates (b1–b6) tested **accuracy only** — which is precisely why the
non-scalable preconditioner shipped undetected. v2 gates every capability on three
axes:

- **a-gates (accuracy):** the b1–b6 legacy set plus the new g4–g10 physics gates.
  "Does it compute the right answer on the reference problem?"
- **s-gates (scalability):** MPI-rank scaling, OpenMP-thread scaling, and hybrid —
  §7. "Does the right answer arrive at scale, and does performance regress?"
  s-gates are *standing* infrastructure: every later capability runs under them
  automatically, so a Q4/Q5/Q6 feature that serializes a loop or breaks halo
  batching fails CI the same day. The **p-gates** (§2B.3) are the portability
  family within this axis: they extend the question from "does it scale across
  ranks and threads" to "does it produce the same answer, and scale, across
  *execution backends*" (serial / OpenMP / CUDA / HIP). They are listed
  separately only because Q3 authors them; once live they are standing gates
  like s1–s4.
- **x-gates (generality):** orientation/symmetry sweeps, BC-side coverage, and
  feature-interaction coverage — §8. "Does it work only in the configuration the
  benchmark happens to use?" This axis exists because the failure mode is
  documented in this codebase's own history: legacy injected subsurface scalars on
  y+ faces only, had no x-edge head conditions at all (P2 added them by mirroring
  the y rules), and scaled the y+ hydrostatic depth term by density but not y−.
  A gate suite that only ever applies a tide on the south edge would never catch
  the north-edge equivalent being broken.

A phase closes only when its capability is green on **all three axes**.

---

## 2. Q1 — Scalable preconditioner

### 2.1 Current state (measured)

Both implicit systems hardcode **CG + block-Jacobi/ICC(0)** at
`src/core/LinearSystem.cpp:47-61`; AMG is only reachable through a raw PETSc options
file, hypre is **not linked** (neither `scripts/ci_install_deps.sh` nor
`build_frehg2_hpc.sh` configures `--download-hypre`), the matrix is host `MATAIJ`
staged element-wise, and PETSc runs pure-MPI (`MPI_THREAD_FUNNELED`). Block-Jacobi's
iteration count grows with rank count — the scalability ceiling is the preconditioner,
not the assembly kernels.

### 2.2 Design

**Precedent.** ParFlow demonstrates multigrid-preconditioned Richards at up to 16k cores
with near-flat iterations (Jones & Woodward 2001, *AWR* 24:763; Kollet et al. 2010,
*WRR* 46:W04201; Maxwell 2013, *AWR* 53:109). ParFlow uses hypre's structured PFMG/SMG;
Frehg2's MATAIJ path maps to **BoomerAMG** (classical AMG — better than smoothed
aggregation on anisotropic problems) with **GAMG** as the no-external-dependency and
future-GPU fallback (GAMG runs natively on MATAIJKOKKOS).

**Deliverables:**

1. `solver:` YAML block (new, schema-validated per §5.4 conventions):
   ```yaml
   solver:
     surface:      { preconditioner: bjacobi-icc | amg | gamg, rtol: ..., atol: ..., max_iterations: ... }
     groundwater:  { preconditioner: bjacobi-icc | amg | gamg, rtol: ..., atol: ..., max_iterations: ... }
     petsc_options_file: ...   # kept; raw options still override (KSPSetFromOptions stays last)
   ```
   Default remains `bjacobi-icc` in v2.0 (b1–b6 goldens were gated under it; flipping the
   default is a v2.1 decision after g3 evidence).
2. PETSc built with `--download-hypre` in `ci_install_deps.sh`, `build_frehg2_hpc.sh`,
   and the docs' install instructions. hypre stays flat-MPI (LLNL's own evaluation finds
   MPI-only ≥ hybrid on most machines; benchmarks already pin `OMP_NUM_THREADS=1`).
3. Baseline BoomerAMG parameter sets wired behind `preconditioner: amg` (overridable via
   options file), from hypre/MOOSE guidance for anisotropic 3D and 2D respectively:
   - `gw_` (3D, thin-layer dz≪dx anisotropy): strong_threshold **0.5** (sweep
     {0.5, 0.6, 0.7} at bring-up — MOOSE warns <0.5 in 3D risks complexity blowup),
     HMIS coarsening, ext+i interpolation, P_max 4, agg_nl 1, symmetric relaxation
     (CG requires a symmetric smoother).
   - `fs_` (2D, 5-point): strong_threshold 0.25, HMIS, ext+i, P_max 4.
4. **Setup-reuse policy** (the main risk: AMG setup can cost several solves and the
   matrix is re-assembled every step):
   - Reuse the hierarchy across steps via `KSPSetReusePreconditioner`; rebuild when the
     CG iteration count exceeds **1.5×** the post-rebuild count, or every N steps
     (configurable, default 50), whichever first.
   - The SWE system additionally forces a rebuild whenever the wet/dry mask changes
     (mask hash comparison — a wet/dry flip is a structural change; the Richards
     active-column mask is static so gw never needs this trigger).
5. **Dry-row hygiene** (required for AMG correctness, review before enabling):
   - dry/inactive rows must be symmetric identity rows (zeroed row *and* column
     contributions) — asymmetric elimination stalls CG under AMG;
   - the placed diagonal must be scaled to the local wet-cell diagonal magnitude,
     not literal 1.0 (poorly scaled identity rows degrade AMG coarsening).
   Audit the current assembly for both properties first; if v1 already satisfies them,
   record that in the Q1 report with the code references.
6. Solver telemetry: per-solve iteration counts, setup/solve time split (`-log_view`
   summarized), and BoomerAMG grid+operator complexity logged to the monitor stream —
   g2/g3 consume these, and users get regression visibility for free.

### 2.3 Gates

| Gate | Basis | Pass criteria |
|---|---|---|
| **g1 — solver invariance** (per-PR) | b1–b6 goldens | All six b-gates pass with `preconditioner: amg` for both systems at 1/2/4 ranks, under the *same tolerances* as the v1 default-mode rank-invariance gates (A5/A8/A14). `-ksp_error_if_not_converged` armed: any KSP divergence is a hard fail. Same for `gamg`. |
| **g2 — iteration flatness vs ranks** (per-PR) | PETSc-style recorded counts | On the 16×-b5 scaling case (the A23 artifact grid): mean CG iterations per system recorded at {1, 2, 4, 8} ranks. Pass: iters(8)/iters(1) ≤ **1.10** for AMG, and max iters ≤ golden + max(3, 10 %). Run the same measurement for bjacobi-icc and *record* (not gate) it — the widening gap is the documented motivation. |
| **g3 — weak scaling** (nightly) | Kollet-2010 protocol, min-over-repeats (A23) | Fixed per-rank subdomain (unit problem replicated per rank), {1, 2, 4, 8} ranks: iteration growth smallest→largest ≤ **15 %**; solver-time weak-scaling efficiency ≥ **0.8** target / **0.5** hard floor on the CI machine (the A23 variance dispensation applies); BoomerAMG operator complexity ≤ **2.0** (memory-blowup guard); setup-time fraction reported so reuse-policy regressions are visible. |

Machine-calibration note: g3 thresholds are provisional until Q0.2 lands real runners;
the first CI execution recalibrates them by amendment, exactly as A23 did for M3 variance.

---

## 2A. Q2 — Run provenance and timing log (inserted 2026-09-10, V2-A5)

*Section numbered 2A, not renumbered into the sequence: source-code comments and
the Q1 records cite this document's section numbers (`v2 plan §2.2.4`, `§7.1.2`,
…), so §3–§12 keep their numbers permanently; only the phase labels moved.*

### 2A.1 Current state (measured)

The timing machinery half-exists and persists nothing:

- A hierarchical `Timer` (`src/core/Timer.hpp`) already wraps the major segments
  (`simulation`, `swe/free_surface`, its `solve`, `swe/velocity`, the gw
  predictor/corrector and its `solve`, `transport`, `coupled_step`) and prints a
  per-section table (count, min/mean/max over ranks) — but only to stdout at the
  end of a successful run (`main.cpp`, `Timer::report`). Q1 added the per-system
  `solver summary fs|gw` lines (solves, iteration stats, rebuilds, setup/solve
  seconds), also stdout-only.
- Nothing records *what* was run: no resolved-parameter echo, no module/BC
  summary, no code version, no rank/thread count, no wall-clock bracket. A
  finished simulation is not self-describing; the timing table dies with the
  terminal scrollback. `run_scaling.py` scrapes stdout — workable for the
  harness, useless for a user's own runs (the owner's Kuan runs prompted this
  phase).
- Not everything is timed: initialization, HDF5 output/checkpoint writes, and
  halo exchange (message passing) have no dedicated sections, so "where did the
  time go" cannot separate I/O and communication from compute.

### 2A.2 Design

**One artifact:** every run writes `run-record.yaml` next to the HDF5 output
(same directory as `output.filename`), rewritten at every output flush and
finalized at exit — so an aborted or wall-clock-killed run still leaves a
truthful partial record with the last completed step. Format YAML (the config
culture of the code); one file per run, overwritten only by the same output
path.

**Contents (top-level keys, schema-checked):**

1. `provenance`: frehg2 version + git SHA (already embedded in the binary),
   build type/compiler, hostname, MPI ranks and decomposition (px × py), OpenMP
   threads, Kokkos backend, start/end wall-clock timestamps (ISO 8601), total
   wall seconds, the input YAML path and its SHA-256, restart parentage
   (checkpoint file + resume time) when restarted.
2. `configuration`: the fully **resolved** config — every parameter with
   defaults materialized (the `FrehgConfig` struct serialized back to YAML),
   not a copy of the input file. What the model actually used, not what the
   user happened to write.
3. `modules`: surface_water / groundwater / transport flags, coupling mode,
   baroclinic activation.
4. `boundary_conditions`: one entry per configured BC — kind, target,
   polygon vertex count and bounding box, value source (constant, or series
   file path + its SHA-256), plus the resolved member-cell count (global).
5. `timers`: the full hierarchical tree exactly as `Timer::report` computes it
   (per section: count, min/mean/max seconds over ranks), **plus the new
   sections**: `init` (everything before the time loop), `io/output`,
   `io/checkpoint`, `halo` (pack/exchange/unpack in HaloExchanger). With these,
   the tree separates shallow-water solver, groundwater solver, transport,
   linear-solve setup/solve (Q1 telemetry), message passing, initialization,
   and I/O — the owner's requested breakdown.
6. `solver`: the Q1 per-system telemetry (solves, iteration mean/max, rebuilds,
   retries, setup/solve seconds) plus preconditioner names.
7. `closure`: the final mass-audit row(s) (surface and gw cumulative budgets)
   and the closure residuals — a one-glance sanity check per run.

**Implementation shape:** a small `io::RunRecord` class (core/io layer; rank 0
writes, values reduced the same way `Timer::report` already reduces);
`Timer` gains a structured accessor (the tree it already builds for the text
report) so the record and the stdout table cannot disagree. The stdout report
stays — the file is additive. Writing the record is excluded from its own
timers' scope and measured by r2's overhead check instead.

### 2A.3 Gates

- **r1 — record correctness (per-PR, a-axis + x-axis).** The regression
  harness gains a post-run hook: for **every** gate case it runs (b1–b6, g-, s-,
  x-lanes), assert the record exists, validates against
  `tools/check_run_record.py` (schema + required keys), and **round-trips**:
  the embedded resolved configuration re-validates against the schema, and
  re-resolving the original input YAML reproduces the embedded configuration
  byte-for-byte. Ranks/threads/decomposition in the record must match the
  launch. Covered in all four regimes (surface-only, gw-only, coupled,
  coupled+transport) and under restart (restart parentage present) — the
  regimes fall out automatically because the existing gate set spans them.
  Negative test (§6.3): a record with a deleted key and a record from a
  mismatched config both fail the checker.
- **r2 — timer coverage and overhead (per-PR, s-axis).** On the perf-baseline
  case: the sum of the top-level timed sections (init + swe + gw + transport +
  halo + io + coupled bookkeeping) ≥ **90 %** of the `simulation` timer at
  1 rank (no major untimed segment), `halo` and `io/*` nonzero in a 4-rank run;
  and total runtime with the record enabled within **1 %** of the Q1 baseline
  (warn) — the record must be free.
- **Physics untouched (s-axis, standing):** the bitwise restart and
  rank-invariance lanes already run per-PR; they are the proof that the new
  timer sections and the record writer perturb nothing. The perf-baseline gate
  guards the overhead.

### 2A.4 Exit

r1 + r2 green; all standing per-PR lanes green; new timer sections documented
in `docs/user-guide/output.md` (run-record chapter) and
`docs/developer-guide/performance.md`; `run_scaling.py` switched to read the
record file instead of scraping stdout (one parser, shared with users);
dod-Q2 / report-Q2 committed.

---

## 2B. Q3 — Performance-portable parallelization (inserted 2026-09-11, V2-A6)

*Section numbered 2B for the same reason 2A is: §3–§12 are cited by source
comments and the Q1/Q2 records and are permanently frozen. Only phase labels
moved.*

**Goal, in the owner's terms:** one source tree, six deployment lanes —
laptop serial, laptop OpenMP, laptop MPI, HPC MPI, HPC MPI+OpenMP, single GPU,
multi-GPU — with the linear solver scalable and backend-correct in every one.

### 2B.1 Current state (measured 2026-09-11)

The abstraction groundwork is genuinely in place; the *linear-algebra boundary*
is where portability stops.

**Already portable (verified by inspection):**

- `src/core/Types.hpp:27-33` defines `ExecSpace = Kokkos::DefaultExecutionSpace`
  and `MemSpace = ExecSpace::memory_space`; every `Field2`/`Field3` is templated
  on `MemSpace`. Physics code never names a backend.
- `scripts/check_forbidden.sh:75-84` already *enforces* that: a CI rule rejects
  `#ifdef KOKKOS_ENABLE_CUDA|HIP` anywhere in `swe/gw/transport/coupling`, and a
  separate rule bans `SharedSpace`/`CudaUVMSpace` so a missing transfer cannot be
  papered over with managed memory.
- Every `create_mirror_view` in the physics modules is **initialization-only**
  (raster ingest, soil-property fill, terrain metrics) — checked file by file.
  The time loop is mirror-free except for the solve, below.
- `HaloExchanger` has a real GPU-aware path with persistent device buffers and a
  single toggle (`runtime.gpu_aware_mpi`, `HaloExchanger.cpp:228-270`).
- Matrix assembly already uses the device-friendly
  `MatSetPreallocationCOO`/`MatSetValuesCOO` pair (`LinearSystem.cpp:169-183`).

**The four gaps:**

1. **`MATAIJ` is hardcoded** (`LinearSystem.cpp:61`) and `MatSetValuesCOO` is
   handed `values.data()` — a `MemSpace` pointer. On a CPU build
   `MemSpace == HostSpace` so this is correct today; on a device build it passes
   a device pointer to a host-memory matrix. **This is a latent defect, not a
   missing feature**, and the existing CUDA lane cannot catch it because that
   lane is compile-only (`scripts/ci_install_deps.sh:28-32`, `FREHG_CUDA=1`
   builds Kokkos with CUDA through `nvcc_wrapper`, no GPU, never linked or run).
2. **The solve stages through host mirrors.** `LinearSystem::solve`
   (`LinearSystem.cpp:227-240`, `305-318`) does `create_mirror_view` + a scalar
   copy loop into `VecGetArrayWrite`, and the mirror image on the way out — two
   device↔host round trips per solve on a device build, which would dominate.
3. **Q1's AMG defaults are CPU-only settings.** §2.2 injects
   `coarsen_type = HMIS`, `interp_type = ext+i`, and exposes `agg_nl`. Per the
   hypre GPU documentation **PMIS is the only device-supported coarsening**
   (HMIS/Falgout/CLJP/Ruge-Stüben are host-only, and hypre's own default silently
   differs by backend); device interpolation is limited to direct (3),
   extended (14), extended+i (6), extended+e (18), BAMG-direct (15); and
   **multipass interpolation — the CPU default at aggressive-coarsening levels —
   has no GPU implementation**, so `agg_nl > 0` requires the two-stage MM-ext+e
   operators instead. `ext+i` is already the right choice (it is the long-range
   interpolation hypre recommends with PMIS/HMIS); the defect is HMIS and the
   unguarded `agg_nl`.
4. **Neither build script configures PETSc for Kokkos.** `build_frehg2_local.sh:218`
   and the CI installer build PETSc with `--download-hypre` only — no
   `--download-kokkos --download-kokkos-kernels`, no `--with-cuda`, no
   `--with-openmp`. So even on CPU there is currently no way to select a threaded
   PETSc backend.

**The consequence for OpenMP, stated plainly:** frehg2's OpenMP thread scaling is
capped today not by AMG but by PETSc's matrix type. PETSc's native `MATAIJ` is not
OpenMP-threaded (`--with-openmp-kernels` threads only *some* routines), which is
exactly what s3 measured as "solve time flat vs thread count". **AMG does not and
cannot fix this**: AMG buys *algorithmic* scalability (iteration count flat vs
problem size and rank count — Q1 delivered it, 23→25 down to 8.5→9.0), which is
orthogonal to whether the remaining work is thread-parallel.

### 2B.2 Design

**One principle:** the execution backend is a *build* choice; the linear-algebra
backend follows from it automatically and is never independently configurable
into an inconsistent state.

| Lane | Kokkos backend | PETSc Mat/Vec | Launch | Notes |
|---|---|---|---|---|
| laptop serial | Serial | `aij` | 1 rank | the debug/validation lane; stays the default |
| laptop OpenMP | OpenMP | `aijkokkos` | 1 rank × N threads | needs ≥ ~1e5 cells to pay off |
| laptop/HPC MPI | Serial or OpenMP | `aij` | N ranks | best CPU lane at benchmark scale (measured) |
| HPC hybrid | OpenMP | `aijkokkos` | ranks × threads, 1 rank/NUMA domain | §2B.3 p3 picks the ratio |
| single GPU | CUDA/HIP | `aijkokkos` | 1 rank | ≥ ~1e6 cells to pay off |
| multi-GPU | CUDA/HIP | `aijkokkos` | 1 rank per GPU + GPU-aware MPI | §2B.3 p6 |

**B1 — device-consistent linear algebra.** `LinearSystem` stops hardcoding
`MATAIJ`. When `MemSpace` is not `HostSpace` the Kokkos types are *forced*
(`MATAIJKOKKOS`/`VECKOKKOS`) rather than defaulted, with a `static_assert`-backed
compile-time invariant; on host builds the type stays selectable (`aij` default,
`aijkokkos` opt-in for the threaded lane) through the existing
`MatSetFromOptions`/`VecSetFromOptions` calls, which already carry per-system
prefixes — so `-gw_mat_type aijkokkos` needs no new plumbing. The host-mirror
staging in `solve()` is replaced by PETSc's Kokkos view accessors
(`VecGetKokkosViewWrite`/`VecGetKokkosView`), making the RHS/solution handoff a
no-copy view alias in every lane. Debug builds assert the returned `PetscMemType`
matches `MemSpace`.

**B2 — backend-aware AMG defaults.** §2.2's `setOptionDefault` block gains a
device branch driven by `MemSpace`: `coarsen_type` PMIS on device / HMIS on host;
`interp_type` `ext+i` in both (already device-supported); `agg_nl` defaulted to 0
on device with `agg_interp_type` forced to a two-stage MM operator whenever a user
raises it. Because these remain `setOptionDefault`, an options file or command
line still wins — the Q1 override contract is preserved exactly. A startup
capability probe logs the resolved hypre configuration into the Q2 run record so a
run is self-describing about which AMG it actually got.

**B3 — the threaded-CPU solver lane.** Both build scripts gain
`--download-kokkos --download-kokkos-kernels` plus `--with-openmp`, which is the
configuration PETSc documents for a CPU-threaded Kokkos backend and which also
propagates OpenMP to downloaded hypre. This is what makes `aijkokkos` meaningful
on a laptop: SpMV, vector ops and the Kokkos-Kernels smoother work run on the
OpenMP execution space instead of one thread. **Expectation management is part of
the design, not an afterthought:** PETSc's own documentation notes the
OpenMP-backed Kokkos path is far less exercised than the GPU backends and that
MPI-only `MATAIJ` is usually the better-tested and faster route on pure CPU, and
hypre's AMG *setup* phase is markedly less thread-scalable than its solve. p2
therefore gates the honest quantity (see below) and records the rest.

**B4 — multi-GPU.** One rank per GPU with explicit rank→device binding from the
node-local rank; PETSc's `-use_gpu_aware_mpi` wired to the same
`runtime.gpu_aware_mpi` toggle the halo exchanger already uses, so there is one
switch for the whole program rather than two that can disagree.

**B5 — build system.** A single `FREHG_BACKEND={serial,openmp,cuda,hip}` cache
variable drives Kokkos backend selection, the PETSc consistency check (CMake
fails at configure time if a device backend is paired with a PETSc lacking
Kokkos support, rather than failing at the first solve), and the CI lane matrix.
`build_frehg2_local.sh` grows the corresponding lanes.

**B6 — guardrails that stand in for absent hardware.** See §2B.4.

### 2B.3 Gates

- **p1 — backend invariance (per-PR, CPU).** b1, b4 and b6 run under
  {Serial + `aij`}, {OpenMP-4 + `aij`}, {OpenMP-4 + `aijkokkos`} and must agree
  within the established rank-invariance tolerances (V2-A4's per-solver bounds
  apply unchanged; `aijkokkos` is registered as its own solver key). This is the
  gate that proves swapping the linear-algebra backend is physics-neutral — the
  single most important property for trusting an untested GPU lane, because it is
  the same code path the GPU will take.
- **p2 — OpenMP solver thread scaling (nightly, CPU).** On the ≥ 1e6 cells/rank
  s3 case with `-mat_type aijkokkos`, at 1 rank × {1, 2, 4, 8} threads.
  **Hard assertions:** iteration counts constant within ±2 % across thread counts
  (thread parallelism must not perturb the algebra), and solve wall time
  **strictly decreasing** from 1→4 threads (the minimum bar: the solve must
  actually thread, which is precisely what today's build fails). **Recorded, not
  gated:** solve-phase thread efficiency, setup-vs-apply split (`PCSetUp` vs
  `PCApply`), and total runtime efficiency — because the ecosystem's threaded-AMG
  performance is not ours to guarantee. p2 supersedes s3's "solve time asserted
  flat" assertion, which encoded the *old* MPI-only-solve design; s3 is amended
  to apply only to the `aij` lane.
- **p3 — hybrid MPI+OpenMP scaling (nightly, CPU).** Extends s4 to the
  `aijkokkos` path at fixed total cores: {8×1, 4×2, 2×4, 1×8}. Results within
  tolerance across all placements (hard); iteration counts within 10 % across
  placements (hard — a placement that changes the algebra indicates a
  decomposition bug, not a tuning result); timings and the winning ratio recorded
  to `scaling-history.md`. The literature's 1-rank-per-NUMA-domain heuristic is
  the expected answer; the gate exists to find out rather than assume.
- **p4 — device build integrity (per-PR, no GPU required).** The CUDA lane is
  promoted from compile-only to **compile + link**, including the unit-test
  binaries, with `FREHG_BACKEND=cuda` and a Kokkos+CUDA PETSc. A HIP
  compile+link lane is added if a ROCm toolchain is available to CI, optional
  otherwise. This cannot execute kernels, and the gate says so explicitly — its
  job is to catch signature, lambda-capture, and linkage breakage at the PR that
  introduces it, which is the failure mode that historically bit this codebase
  ("fixed bugs for building on gpu", Q0.1).
- **p5 — host/device pointer discipline (per-PR, static + runtime).** Three
  mechanisms, because this is the bug class the CPU lanes structurally cannot
  see: (a) a `check_forbidden.sh` rule rejecting raw `.data()` on a `MemSpace`
  view passed to any PETSc API outside the sanctioned `LinearSystem` wrappers;
  (b) the B1 compile-time invariant, with a **negative test** (§6.3) that a
  deliberately inconsistent backend/matrix-type pairing fails to compile;
  (c) debug-build `PetscMemType` assertions on every Vec/Mat handoff. Gap 1 of
  §2B.1 is the worked example this gate is designed around.
- **p6 — GPU acceptance bundle (authored here, executed by the owner).** Not a CI
  gate — a *deliverable*: `scripts/gpu_acceptance.sh` plus
  `docs/developer-guide/gpu-acceptance.md`, which runs, on the owner's GPU HPC,
  (i) single-GPU b1/b4/b6 against the same goldens the CPU lanes use, with the
  p1 tolerances, (ii) a single-GPU vs 1-rank-CPU timing comparison on a ≥ 1e6-cell
  case, (iii) a 2- and 4-GPU weak-scaling run with GPU-aware MPI on and off, and
  (iv) the Q2 run record from each, which is the artifact to send back. The script
  self-checks and prints PASS/FAIL per item against written criteria, so the
  result is unambiguous without interpretation. It is dry-run on the CPU lanes
  before Q3 closes, so the only untested variable on the HPC is the hardware.

### 2B.4 The no-GPU constraint, and one correction to the premise

The owner's plan — gate OpenMP and MPI+OpenMP here, verify GPU on the HPC — is
the right division of labour and is what §2B.3 implements. One premise needs
correcting, and it shapes p4/p5:

> *"Since both OpenMP and GPU are shared-memory, with Kokkos, if the code works
> well with OpenMP, it should work on GPU."*

True for the **kernels** — a Kokkos `parallel_for` that is correct and
thread-invariant under OpenMP is very likely correct under CUDA, and that is
exactly why p1/p2/p3 are worth gating here. It is **not** true for the
host/device **boundary**, for a structural reason: under the OpenMP backend
`MemSpace == HostSpace`, so every mirror, raw pointer and MPI buffer aliases host
memory and every memory-space mistake is invisible. Gap 1 of §2B.1 is a live
example — `MatSetValuesCOO(values.data())` is correct under OpenMP and wrong
under CUDA, and no amount of OpenMP testing would reveal it. p4 (compile+link)
and p5 (static invariant + negative test + memtype assertions) exist specifically
to cover the blind spot OpenMP leaves, and they are the reason Q3 can close
honestly without hardware.

**Release-status consequence (must not be finessed).** §8.3 forbids
"authored-unexercised" features in a release, and P5's CUDA lane is precisely the
authored-unexercised state that rule was written against. Shipping a GPU lane
that has never executed would violate the plan's own rule. Resolution: the GPU
lanes ship as **`experimental`**, a new documented status meaning *compile-verified,
statically invariant-checked, and physics-verified on the equivalent CPU backend,
but not executed on the target hardware by CI*. The status is recorded in the
§8.3 feature-interaction table, printed by the binary at startup on a device
build, and stated in the release notes. It is cleared — feature promoted to
supported — when the owner returns a passing p6 bundle, which is a one-line table
edit plus the artifact in `docs/developer-guide/`. This keeps Q7's release
criteria meaningful instead of quietly weakening them.

### 2B.5 Exit

p1, p2, p3, p5 green on this machine; p4 compile+link green in CI; p6 bundle
authored, documented, and dry-run through the CPU lanes; `aijkokkos` and
`FREHG_BACKEND` documented in `docs/user-guide/installation.md`,
`parameters.md` and `developer-guide/performance.md`, with the lane-selection
table of §2B.2 reproduced for users; the backend and resolved AMG configuration
recorded in every run record (Q2 lockstep: a schema addition requires the
matching `resolvedConfigYaml` change in the same PR); s3 amended per p2;
§8.3 table carries the `experimental` GPU rows; dod-Q3 / report-Q3 committed.

**Q3 closes with GPU performance formally unverified.** That is a tracked,
named open item with an owner action attached, not a silent gap.

---

## 3. Q4 — Evaporation (open water + soil) and salinity coupling

### 3.1 Current state (measured)

- Open water: single spatially-uniform prescribed rate subtracted unconditionally from
  `eta` (`src/swe/SurfaceSources.cpp:43-55`); no per-cell mask (rain has one), no
  atmospheric physics. The legacy aerodynamic model (`evap_model=1`) was dropped for
  hardcoding 20 °C air/pressure/humidity — the *formulation* was sound, the constants
  were not.
- Soil: prescribed positive top-face flux, cut off at residual moisture
  (`src/gw/Corrector.cpp:332-337`, 447-451). No potential→actual reduction.
- Salinity: surface evaporative concentration is already correct and salt-conserving
  (`src/transport/SurfaceTransport.cpp:349-360`); subsurface concentration is throttled
  by the monotonicity limiter with a hardcoded `hi += 0.01` allowance
  (`src/transport/SubsurfaceTransport.cpp:556-560`).

### 3.2 Design

One shared **bulk-aerodynamic module** (`src/atm/`, header + kernels) providing:
saturation vapor pressure e_sat(T) (Tetens/Magnus — the Geng & Boufadel paper uses
Tetens Eq. (5): e_sat = 0.6108·exp(17.27·T/(T+237.3)) kPa), specific humidity
conversions, aerodynamic resistance R_air(U) — with Liu et al. (2006)
R_air = 94.909·U^(−0.9036) as one selectable fit — and the evaporative flux

  E_g = (ρ_a / R_air) · (q_g − q_a)      [Mahfouf & Noilhan 1991 form]

New `atmosphere:` YAML block: air temperature, relative humidity (or specific humidity),
pressure, and wind speed as SeriesOrConstant (wind speed shared with the Q6 wind block);
all previously hardcoded legacy constants become required inputs — that is what makes
the dropped `evap_model=1` general enough to reinstate.

**Open-water evaporation** (three selectable modes on `surface_water.evaporation`):
1. `prescribed` — current behavior, unchanged (b1 golden fidelity).
2. `bulk` — E from the bulk-aerodynamic module with water-surface q_g = q_sat(T_s)
   (T_s prescribed until Q5, then the transported temperature); applied per cell, wet
   cells only, clamped by available depth **with the shortfall audited** (extend the
   existing clampVolume audit rather than silently creating volume).
3. `series` per-cell-masked variant of prescribed (fixes the no-mask asymmetry vs rain).

**Soil evaporation** (potential→actual):
- Potential rate E_p from the same module; **actual** rate limited by the surface-layer
  moisture through the soil relative humidity α₁ (Geng & Boufadel Eq. (6),
  Barton/Lee-Pielke: α₁ = min(1, 1.8·w_g/(w_g + 0.30)) with w_g = φ·S the surface
  water content), so E → 0 smoothly as the surface dries — replacing the hard wc > wcr
  cutoff for the `bulk` mode (the prescribed-flux BC path keeps legacy behavior).
- Applied as the groundwater top-face flux through the existing coupled bookkeeping
  (`cplEvap` audit), so mass accounting is unchanged in structure.

**Salinity coupling:**
- Subsurface: replace the magic `hi += 0.01` limiter allowance with the exact
  evaporative-concentration factor for the top cell (the concentration multiplier
  V/(V − E·A·dt) is known in-step, so the limiter ceiling can be raised by precisely
  the physical amount; behavior with the legacy 0.01 preserved under a
  `transport.legacy_evap_allowance: true` toggle for b6 golden fidelity).
- Add a Cauchy (zero-total-flux) top salt condition under evaporation as in Geng &
  Boufadel Eq. (7): water leaves, salt stays — this is the g5 physics.

### 3.3 Gates

| Gate | Basis | Pass criteria |
|---|---|---|
| **g4 — analytic drawdown + evaporative concentration** (per-PR) | Closed-form conservation | (a) Closed flat basin, uniform prescribed E: η(t) = η₀ − ∫E dt exact to ≤ 1e-6 of total drawdown until the dry threshold. (b) Same run with salinity: s(t) = s₀·V₀/V(t) to relative error ≤ 1e-4, global salt mass constant to 1e-10 relative. (c) bulk mode with *constant* met forcing: E constant computable offline, same criteria — gates the formula end-to-end. (d) drain-to-dry variant: positivity + audited shortfall only. |
| **g5 — Geng & Boufadel (2015) bare-soil salinization** (nightly-class) | *J. Hydrology* 524:427-438, code-to-code (MARUN) | Configuration §5 below. Metrics: (i) evaporation-rate time series shape — initial rate ~1.5e-7 m/s decaying one order of magnitude within the first ~7 h and toward ~1e-10 m/s by 40-50 h (two-regime check: rate at 10 h within a factor of 2 of MARUN's digitized Fig. 3 curve; rate ratio E(50h)/E(0) ≤ 1e-2); (ii) horizontally averaged moisture-ratio and salinity profiles vs digitized Fig. 4/Fig. 9 at 20 h and 50 h — RMS ≤ 0.05 in moisture ratio, salinity profile RMS ≤ 10 g/L with the near-surface peak (>60 g/L above z = 1.9 m at 50 h) reproduced; (iii) salt mass conserved to ≤ 4 % (the paper's own MARUN budget bound — do not gate tighter than the reference's self-consistency); (iv) qualitative density gate: with β = 7.44e-4 the upper saline plume must remain, and the β = 0 control must show measurably deeper plume spreading (Fig. 7 contrast) — an inequality assertion, not a curve match. |

### 3.4 Why Geng & Boufadel 2015 is the right g5 (assessment of the owner's suggestion)

The paper (`Geng_and_Boufadel_JH2015.pdf`) is well suited, with caveats:

**For:** it exercises *exactly* the new capability stack in one case — bulk-aerodynamic
evaporation with humidity feedback from pore moisture (Eqs. 1-6), moisture-limited
actual-vs-potential decay (the two-regime rate curve is a discriminating fingerprint —
a wrong α₁ or missing limiter cannot fake it), evaporative salt concentration under a
Cauchy top BC (Eq. 7), and density feedback on variably saturated flow — with the same
van Genuchten retention and the **same β = 7.44e-4 density coefficient Frehg2 already
hardwires** (Table 1). Domain 50 × 2 m, medium mesh 501 × 27 (~13.5k nodes), 50 h
horizon, dt at Courant < 0.3 — comfortably nightly-class cost, same tier as b6. All
parameters are in Table 1 (α = 4.75 m⁻¹, n = 8.5, K₀ = 9e-4 m/s, φ = 0.41, Sr = 0.02,
S₀ = 1e-5, αL = 0.1 m, αT = 0.01 m, τDm = 1e-10 m²/s; forcing T_s = 20 °C, U = 1 m/s,
q_a = 20 % of q_sat, P₀ = 101.325 kPa; IC saturated, 25 g/L).

**Caveats to encode in the gate design:**
1. It is a **numerical study** (MARUN), not a lab experiment — g5 is a code-to-code
   gate like b5's envelope, not a physical-truth gate like b6. References must be
   digitized from Figs. 3, 4, 9 (and 7 for the density contrast); digitization carries
   the same ~1-cell-class error b3 already budgets for. The horizontally averaged
   profiles (Figs. 4, 9) are the robust targets; individual salt fingers (Fig. 5) are
   instability-set and **must not** be gated pointwise — position of fingers is
   sensitive to perturbations (the paper itself notes onset is Rayleigh-criterion
   physics). Gate finger *existence* qualitatively at most.
2. MARUN is 2D vertical-plane finite element; Frehg2 runs it as an x-z slice (ny = 1
   plus ghost, the b2/b3 pattern). The evaporation zone x ∈ [1, 49] m maps to the
   groundwater top-flux polygon; side/bottom no-flow.
3. The paper's evaporation acts on the *subsurface* top boundary (no ponded water) —
   it gates soil evaporation + salinity only. Open-water evaporation is gated by g4(c)
   and, with temperature, g7. That split is fine: no single benchmark covers both.
4. Frehg2's subsurface transport currently has the y+-only scalar-injection convention
   and dispersion-on-face-flux quirk (`frehg2-validation` findings); the case setup
   must route salt through the top Cauchy condition — a new BC kind, so unit tests for
   it come first (§11.1-style: the BC lands complete, not stubbed for the benchmark).

**Verdict: adopt as g5.** It is the strongest available soil-evaporation benchmark short
of a lab dataset, and no comparable lab dataset with full parameter disclosure exists at
this cost tier (the Gran et al. 2011 column experiment is the lab alternative — reserve
it as extended validation, not a gate, since its mineral precipitation is out of scope).

---

## 4. Q5 — Temperature transport (surface + subsurface)

### 4.1 Design

- Generalize `ScalarSolver` to two registered scalars (salinity, temperature) rather
  than a general N-scalar rewrite: the v1 schema already reserved the list form; state,
  halo entries, BC lists, and ledgers become per-scalar; the advection/dispersion
  kernels are already scalar-agnostic. Temperature advects with the same schemes;
  subsurface effective thermal conduction/dispersion enters through the existing
  dispersion-tensor slot with thermal diffusivity λ_eff/(ρc)_w and the retardation
  factor from the sediment heat capacity ((ρc)_bulk/(ρc)_w — the "thermal front
  velocity" scaling); surface horizontal diffusion likewise.
- Surface heat exchange: net flux Q_net = Q_sw + Q_lw,in − Q_lw,out − Q_lat − Q_sens
  with the latent term from the Q4 bulk module (one shared implementation) and an
  **equilibrium-temperature mode** (Edinger form: dT/dt = −K_e(T − T_e)/(ρ c_p h)) for
  the analytic gate; GLM's still-air lower bound on the transfer functions adopted for
  low-wind robustness.
- Density: extend r_rho/r_visc to ρ(s, T) with a linear (or UNESCO-polynomial,
  selectable) thermal expansion term; the compile-time β constants move to
  schema-validated config with the current values as defaults (fixes the A-noted
  compile-time-constant limitation without changing gated results).

### 4.2 Gates

| Gate | Basis | Pass criteria |
|---|---|---|
| **g6 — subsurface heat, analytic** (per-PR, three parts) | Bredehoeft & Papadopulos (1965) *WRR* 1:325; Ogata–Banks (1961) heat form as specified by the OpenGeoSys benchmark; Stallman (1965) *JGR* 70:2821 | (a) **B&P steady Péclet sweep** (smoke tier): 1D column L = 10 m between fixed temperatures 20/10 °C, k_fs = 2.0 W/m/°C, Pe ∈ {−5, −1, 0, 1, 5}: time-marched steady profile vs the exponential closed form, max-norm ≤ 1 % of ΔT. (b) **Ogata–Banks heat step** with the OGS parameter set (ρc_eff = 2e6 J/m³/K, λ = 2.2 W/m/K → α = 1.1e-6 m²/s, v = 1.5e-6 m/s, T 300→330 K, observation points x = 1/5/10/20/50 m over 500 d): breakthrough vs the erfc closed form, relative error ≤ 1 % on the superbee scheme at the graded mesh (OGS asserts 0.25 % with FEM; the upwind option is recorded, not gated). (c) **Stallman damped sinusoid**: diel forcing (A = 5 °C, P = 86400 s) on a 2–5 m column, sub-cases q_z ∈ {0, +5e-6, −5e-6} m/s with (ρc)_w = 4.184e6, (ρc)_bulk = 2.96e6 J/m³/K, λ = 2.0 W/m/K; after 3 spin-up cycles, amplitude ratio at 2–3 depths within **2–5 %** and phase lag within **10 min** of the closed form. Amplitude and phase are gated separately — amplitude catches numerical dissipation, phase catches dispersion error. |
| **g7 — surface heat exchange, analytic** (per-PR, two parts) | Edinger et al. (1968) *WRR* 4:1137 equilibrium temperature (the CE-QUAL-W2/Delft3D "excess temperature" linearization); van Genuchten et al. (2013) *J. Hydrol. Hydromech.* 61:146/250 exact river ADE+decay solutions (Neilson et al. 2012 heat-native form) | (a) **Equilibrium relaxation**: still closed basin, h = 1 m, T₀ = 30 °C, T_e = 20 °C, K_e = 30 W/m²/K (e-fold ≈ 1.6 d, run 10 d): T(t) vs the exact exponential within 0.5 %; energy ledger closes to 1e-8/step (mirror of the scalar-mass audit). Second stage: full bulk-formula mode under constant met forcing must reach the T_e computed offline by root-finding the net-flux zero, within 0.05 °C — the consistency check on the nonlinear flux terms. (b) **Channel thermal plume**: uniform flow u = 0.5 m/s, D_L = 5 m²/s, h = 1 m, K_e = 25 W/m²/K, upstream step 25 °C into 15 °C ambient = T_e: steady profile L2 ≤ 1 %, mid-channel transient breakthrough ≤ 2 %; the K_e = 0 limit must reduce to the g6(b) machinery (one harness, two configs). |
| **g8 — thermal convection onset** (nightly, *recommended*) | Horton–Rogers–Lapwood: critical Rayleigh number Ra_c = 4π² ≈ 39.48 (Horton & Rogers 1945; Lapwood 1948) | Saturated porous slab heated from below, linear ρ(T), seeded sinusoidal perturbation (wavenumber π/H): at Ra = 30 the perturbation decays to the conduction profile; at Ra = 50–60 the Nusselt number exceeds 1.05 persistently. A deterministic analytic threshold — the rare pass/fail gate available for density-coupled flow; it gates the ρ(T) Darcy coupling that g6/g7 never touch. |

Notes: g6/g7 are analytic and cheap — per-PR tier, unlike b-class gates; the Kuan
td tidal-phase-averaging harness is reusable for the Stallman amplitude/phase
extraction. The **thermal Elder problem is explicitly rejected as a gate** (≥ 11
grid-dependent steady states — central upwelling flips with refinement; *Fluids*
2017 retrospective); HRL replaces it. A coupled surface→subsurface thermal gate is
composed rather than sourced: no fully specified analytic coupled thermal benchmark
exists in the literature (the HGS reference, Brookfield et al. 2009, validates on a
field site) — so the coupled variant runs g6(c) with the sinusoidal temperature
supplied by the *surface module* through the infiltration exchange instead of a
prescribed BC; the subsurface response must still match Stallman. That composition
gates the heat-exchange coupling path with zero new reference data.

### 4.3 Follow-on cross-check (not a gate)

After Q5, rerun **g5 with the full T-dependent module active** at the paper's constant
20 °C: results must be statistically identical to the Q4 run (temperature constant ⇒
no thermal feedback). This is the cross-coupling regression that catches accidental
activation, in the spirit of the v1 "b1 stays golden with modules off" discipline.

---

## 5. Q6 — Wind stress

### 5.1 Current state (measured)

Quadratic stress with constant Cd = 0.0013 and thin-layer attenuation
(`src/swe/SweFormulas.hpp:65-83`, `WindConfig` in `src/core/Config.hpp:104-112`);
speed/direction time series with the direction linearly interpolated (legacy held it
piecewise-constant); **no benchmark has ever exercised the wind path** (implemented,
ungated — completion report §6).

### 5.2 Design

- Selectable Cd(U10): `constant` (current, default), `garratt` (0.75 + 0.067·U10 in
  10⁻³, capped — ADCIRC default and cap practice), `smith-banke`, `wu`, `large-pond`,
  each with a configurable cap (default 3.5e-3). Relative-velocity form (current) kept.
- Wind forcing generalized to (u₁₀, v₁₀) component series (avoids the direction-
  interpolation-through-north ambiguity; speed/direction input retained and converted,
  with the interpolation convention documented and unit-tested across the 360°→0° wrap).
- The thin-layer attenuation stays as-is (legacy-faithful, benchmark-neutral: both gate
  basins are deep relative to hD).

### 5.3 Gates

| Gate | Basis | Pass criteria |
|---|---|---|
| **g9 — steady wind setup** (per-PR) | TELEMAC-2D validation case 10 config + closed-form; Dean & Dalrymple Eq. 5.96 | (a) Nonlinear variant replicating TELEMAC's numbers: 500 × 100 m basin, 2 m rest depth, no friction/viscosity, effective τ/ρ_w = c·U² with c = 1.2615e-3, U = 5 m/s: steady end-to-end setup vs analytic h(x)² = C + 2τx/(ρ_w g) (volume-conserving C) within **0.5 %** when spun up from rest with mild relaxation (TELEMAC achieves 0.002 % initialized *at* the solution — spin-up from rest is the honest harder test); (b) linear-regime variant (physical Garratt Cd, U = 15 m/s, 10 m depth): dη/dx = τ/(ρ g h) within 1 %; (c) sloping-bottom variant vs quadrature reference within 1 % — exercises depth-varying setup and grazes (not crosses) the wet/dry threshold. |
| **g10 — seiche relaxation** (per-PR) | Merian formula T = 2L/√(gH) | From the g9(a) steady state, wind switched off (time-series step to zero — also gates the forcing-series edge handling): fundamental period from zero-crossings of η at the basin end within **2 %** of Merian at Cr ≲ 1 (θ-scheme damping documented at the run's Cr); friction-free amplitude decay ≤ a few %/period recorded (regression-tracked, not hard-gated, since it is scheme dissipation). |
| **unit battery** (per-PR) | Hand-computed | Cd(U10) tables for every law at U10 ∈ {0, 5, 10, 20, 30, 40, 60} m/s including caps and breakpoints, exact to 1e-12 relative; direction-wrap interpolation cases. |

Rejected for gating (recorded so nobody re-litigates): Lake Okeechobee/Erie hindcasts
(data in scanned figures, wind-field construction is itself research — fine for a paper
figure, wrong for CI); Csanady topographic gyres (needs linear-friction mode, published
checks are qualitative). Optional stretch: Kraus & Militello (1999) 1D sea-breeze-forced
bay closed form, if a time-varying-wind gate is later wanted.

---

## 6. Gate authoring: who builds what, and in what order

**None of the g/s/x gates exist at plan time** — no case YAMLs, no reference data, no
tolerance files. Authoring them is part of the plan's execution, not a precondition.
The division of labor is explicit so nothing is silently assumed of the owner.

### 6.1 Gate-first rule (binding)

Within each phase, the gate is authored and merged **before** the physics it gates:
case YAML(s) + reference generator (or digitized data) + tolerance file + harness
entry, demonstrated *failing* (or skipping with "capability absent") against the
pre-phase code, then turned green by the capability PRs. This is the executable
analogue of test-first, and it prevents the gate from being quietly fitted to the
implementation. A gate's own correctness is established at authoring time by the
self-checks in §6.3.

### 6.2 Provenance table: every gate's inputs and reference source

| Gate | Case inputs | Reference | Produced by | Owner action |
|---|---|---|---|---|
| g1 | existing b1–b6 YAMLs + `solver:` override | existing goldens/criteria | AI (config only) | none |
| g2/g3, s1–s4 | existing 16×-b5 + one synthetic replicated-unit case | none (assertions on measured iterations/timings) | AI | approve thresholds after first real-runner calibration (amendment) |
| g4 | new tiny closed-basin YAMLs (§3.3) | closed form, computed in-script | AI | none |
| **g5** | new x-z slice YAML from Geng & Boufadel Table 1 (all parameters published) | **digitized CSVs from Figs. 3, 4, 7, 9** | AI digitizes; committed with `DIGITIZATION.md` (axis calibration, extraction method, point count) | **verify the digitization** — §6.4 |
| g6 (a/b/c) | new 1D-column YAMLs (§4.2 parameters) | closed forms (B&P exponential, Ogata–Banks erfc, Stallman damped sinusoid), computed in-script | AI | none |
| g7 (a/b) | new basin/channel YAMLs | closed forms (Edinger exponential; van Genuchten et al. river ADE solution; offline root-find for the nonlinear T_e) | AI | none |
| g8 | new porous-slab YAML pair (Ra = 30 / 50–60) | analytic threshold Ra_c = 4π² (bracket assertion, no curve) | AI | none |
| g9 (a/b/c) | new basin YAMLs (§5.3 configs incl. the TELEMAC case-10 replica) | closed forms (setup profile; quadrature for the sloping bottom) | AI | none |
| g10 | g9(a) + wind-cutoff series | Merian formula | AI | none |
| Cd battery | none (unit tests) | hand-computed tables from the published formulas | AI | none |
| x-gates §8.1/§8.4 | transforms of existing + new cases, auto-generated by the harness | the case's own identity run | AI | review the **asymmetry-exemption table** (each entry is an owner decision: keep legacy quirk or fix) |
| x-gate §8.2 | none (matrix backfill = unit tests + schema assertions) | n/a | AI | decide fix-vs-reject for each silent-gap row found (amendment each) |

So: **one gate (g5) has externally sourced reference data; everything else is
self-contained.** This is deliberate — the benchmark selection in §§3–5 preferred
closed-form references precisely so that AI-automated implementation does not
bottleneck on data acquisition. (Optional later additions that would need owner-side
data work are already quarantined in §11: Lake Hefner PP-270 Tier-2, Gran et al. 2011.)

### 6.3 Self-checks a new gate must pass at authoring time

1. **Reference sanity:** the in-script closed form is checked against at least two
   independently computed spot values (documented in the gate script's header — e.g.
   a hand calculation and a `scipy` evaluation) before any model output touches it.
2. **Sensitivity:** the gate demonstrably *fails* when fed a deliberately corrupted
   run (e.g. the capability disabled, or a parameter perturbed 10 %) — committed as a
   negative test, so a gate that can't fail can't pass.
3. **Convergence:** for discretization-sensitive gates (g6b, g7b, g9), the metric is
   shown to improve under mesh refinement in a one-off study archived with the gate
   (tolerances are then set at the gated resolution with headroom stated).

### 6.4 The g5 digitization protocol (the one owner touchpoint)

AI extracts the curves from the PDF figures (axis-calibrated, ~20–40 points per
curve), commits CSVs + `DIGITIZATION.md` + an **overlay plot** (digitized points
re-plotted over the original figure image). Owner action is a visual check of the
overlay plots — minutes, not hours — plus sign-off recorded in the DoD. If any
figure proves too low-resolution to digitize within the gate's tolerance budget,
the fallback (by amendment) is to gate on the robust scalar fingerprints only
(§3.3 g5(i) rate-decay metrics + salt-mass bound), which need no digitization.

### 6.5 Harness wiring

- New gate scripts extend `tests/regression/run_regression.py` + `tolerances/*.yaml`
  exactly as b1–b6 (one YAML per gate; the §11.2 machine-checkable-gate harness picks
  them up unmodified).
- g2/g3/s-gates read the Q1 solver telemetry; thresholds live in `tolerances/*.yaml`
  so recalibration-by-amendment is a one-line diff with provenance.

## 7. s-gates — parallel scalability as a standing gate axis (MPI **and** OpenMP)

Lesson encoded: v1 measured performance once (P5 report) but *gated* nothing on it,
so CG+bjacobi's rank-degradation was invisible to CI. v2 makes scalability a standing,
per-PR/nightly gate family that every capability inherits.

### 7.1 The gate matrix

Every a-gate case (b1–b6, g4–g10) declares a **parallel execution matrix** in its
tolerance YAML; the harness runs the case over the matrix and applies two kinds of
assertions:

1. **Parallel correctness (per-PR):**
   - MPI: results at {1, 2, 4} ranks within the established rank-invariance
     tolerances (A5/A8/A14 machinery, reused unchanged). Decompositions must include
     a non-square split (1×4 *and* 4×1 *and* 2×2) — not just the default — because
     row-major vs column-major splits stress different halo paths.
   - OpenMP: results at {1, 4} threads (Kokkos OpenMP backend) bitwise-identical
     where the v1 OpenMP lane already proves it (restart determinism at 4 threads
     exists today) or within rank-invariance-class tolerance where reductions
     reorder; every *new* physics kernel (evaporation, heat, wind) must be
     thread-invariant from its first PR — atomics/reduction discipline is reviewed
     against this gate, not against hope.
   - Hybrid smoke: one 2-rank × 2-thread run per gate label, same tolerance.
2. **Performance regression (per-PR, cheap):** per-module time-per-cell-per-step
   (min over 3 repeats, the A23 protocol) recorded to a tracked baseline file;
   fail if any module regresses > 25 % vs baseline at 1 rank/1 thread, warn > 10 %.
   Baseline updates are explicit commits with justification (same discipline as
   tolerance changes).

### 7.2 Scaling gates proper (extend Q1's g2/g3 beyond the solver)

- **s1 — MPI strong scaling (nightly):** the 16×-b5 case at {1, 2, 4, 8} ranks:
  end-to-end efficiency ≥ 70 % at 4 ranks (the v1 P5 bar, now permanent), plus the
  g2 iteration-flatness assertion. Reported per module so a regression names its
  culprit.
- **s2 — MPI weak scaling (nightly):** g3 as specified in §2.3, permanent.
- **s3 — OpenMP thread scaling (nightly):** same case, 1 rank × {1, 2, 4, 8}
  threads on a grid sized ≥ 1e6 cells/rank (the v1 finding that benchmark-scale
  grids are kernel-launch-bound at OMP>1 means *small* grids must not be used to
  judge threading). Assertions are honest about Amdahl: **Kokkos kernel time**
  efficiency ≥ 60 % at 4 threads, while the PETSc solve time is asserted *flat*
  (the solve is MPI-only; a thread count that slows the solve indicates oversubscription
  misconfig). Total-runtime thread efficiency is recorded, not gated, with the solve
  fraction printed next to it.
  **Amended by V2-A6:** the flat-solve assertion holds **only for the `aij`
  lane**, where it is a true statement about PETSc's non-threaded default matrix
  type. Under `-mat_type aijkokkos` the solve is *expected* to speed up and gate
  p2 (§2B.3) asserts that it does. Applying the flat assertion to the Kokkos lane
  would gate in the very limitation Q3 removes.
- **s4 — hybrid placement sanity (nightly):** {4 ranks × 1 thread} vs
  {2 ranks × 2 threads} vs {1 rank × 4 threads} on the same grid: results within
  tolerance, timings recorded. This is the configuration users will actually run on
  HPC nodes and it exercises rank×thread interaction (pinning, halo/pack threading)
  that neither s1 nor s3 sees alone.

### 7.3 Reporting

`scripts/run_scaling.py` grows a `--matrix` mode emitting one JSON artifact per
nightly run (per-module times, iterations, efficiencies); a tracked
`docs/developer-guide/scaling-history.md` table is appended per release so the
trend across v2.0 → v2.x is inspectable. Threshold recalibration after Q0.2's first
real-runner execution follows the amendment protocol (A23 precedent).

## 8. x-gates — generality beyond the benchmark configurations

Lesson encoded: a gate case fixes one orientation, one BC side, one decomposition,
one feature combination. Code can pass every gate while being broken in every
configuration the gates don't visit. Three instruments close this gap.

### 8.1 Symmetry/orientation battery (per-PR)

For a designated set of small cases (one per module: a b1-derived SWE case, a
b3-derived gw case, a tracer-advection transport case, plus one new case per Q4/Q5/Q6
capability), the harness auto-generates the **8 dihedral transforms** of the
configuration (identity, 90/180/270° rotations, x-mirror, y-mirror, both diagonals —
rotations only where nx=ny or the case is transposable): bathymetry/DEM, ICs, BC
polygons, and forcing directions (wind!) are all transformed together, the case is
run, and the output fields inverse-transformed and compared to the identity run.

- Pass: max-norm agreement within the strict rank-invariance tolerance (these runs
  use the strict solver mode, where v1 proved machine-precision invariance is
  attainable, so the gate can be tight: 1e-12 relative).
- **Documented-asymmetry exemption table:** legacy-faithful asymmetries are *waived
  explicitly, per term, with provenance* — e.g. the preserved legacy y+/y− density
  scaling split (`ghostHead`, legacy :777 vs :787). The table lives in
  `docs/theory/` next to the removed-features registry and is the only mechanism for
  excusing a symmetry failure; an undocumented failure is a bug. (Writing this table
  is itself an audit: every entry is a candidate v2 physics fix, decided by the owner.)
- New v2 physics (evaporation, heat exchange, wind stress, Cd laws) gets **no
  exemptions**: it must be exactly symmetric by construction.

### 8.2 BC-kind × side coverage matrix (per-PR unit tier + schema)

An enumerated matrix: every BC kind (surface: Eta/Discharge/Velocity/Outflow +
scalar/temperature values; groundwater: Head/HeadHydrostatic/Flux/Gravity; new Q4/Q5
kinds) × every applicable target (N/S/E/W edge via polygon; gw top/bottom) ×
{uncoupled, coupled}. For each cell of the matrix, exactly one of:

1. a unit or mini-regression test exercises it (the test id is recorded in the
   matrix), or
2. the **schema rejects it loudly** (with a test asserting the rejection and its
   message), or
3. a documented-limitation row explains why it is accepted but unverified
   (this category must be empty at Q7 release — it exists only as a mid-phase state).

The matrix is a checked-in CSV (`tests/coverage/bc_matrix.csv`) validated by a CI
script that cross-references test names — the same lockstep mechanism as the v1
schema↔docs check. First deliverable of the x-gate work is **backfilling the matrix
for v1 features**, which will surface today's silent gaps (the y+-only subsurface
scalar injection is the known example: v2 either generalizes it to all sides or makes
the schema reject the other sides — silent no-op is the one forbidden outcome).

### 8.3 Feature-interaction coverage (release tier)

A tracked table (`docs/developer-guide/feature-coverage.md`): rows = features
(wind, rain, evaporation modes, transport schemes, terrain-following, subcycled
coupling, masked columns, each preconditioner, restart, each new v2 capability);
columns = the gate/test that exercises each feature *and* each shipped pairwise
interaction that plausibly couples (wind+wetting-drying, evaporation+transport,
temperature+density, AMG+restart, subcycling+transport, terrain-following+heat...).
Rule inherited from v1's hard lesson (wind: implemented P1, never exercised by any
gate, still unvalidated at release): **a feature or interaction with an empty cell
either gets at least a conservation/invariance test before Q7, or its config
combination is schema-rejected**. "Authored-unexercised" is allowed mid-phase and
forbidden in a release.

### 8.4 Decomposition/edge-case battery (per-PR, extends the existing MPI tier)

Standing runs already partially present in v1, completed and made explicit: nx=1 and
ny=1 slice domains **in both orientations** (the Geng g5 case provides one; its
transpose is auto-generated), non-divisible rank splits, a single-cell-wide BC
polygon on each edge, a fully masked column adjacent to each edge, and restart in
the middle of every new capability's forcing transient (wind step-down, evaporation
under drying, heat relaxation).

## 9. Phase exit criteria summary

Every phase additionally obeys the §6.1 gate-first rule: its gates are authored,
self-checked (§6.3), and merged *failing* before the capability code lands.

| Phase | Blocking | Also required |
|---|---|---|
| Q0 | b1–b6 green on real CI runners | repos reconciled (§1.1 Q0.1); superslab diagnosis filed; workspace cleanup done (§1.1 Q0.4) |
| Q1 | g1, g2 (g3 nightly-armed); s-gate harness live (§7.1 matrix + s1–s4 armed) | hypre in all build paths; solver YAML documented; dry-row audit note; performance-baseline file seeded |
| Q2 | r1, r2; restart/rank-invariance lanes unchanged (log must not perturb physics) | run-record schema checker in tools/; new timer sections (init, io, halo); overhead measured < 1 % on perf-baseline; docs (output.md run-record chapter) |
| Q3 | p1, p2, p3, p5; p4 compile+link in CI | p6 GPU acceptance bundle authored + dry-run; PETSc Kokkos in all build paths; `FREHG_BACKEND` + `aijkokkos` documented; backend/AMG config in the run record; s3 amended per p2; §8.3 GPU rows marked `experimental` |
| Q4 | g4, g5; s-matrix green for new kernels; x-battery for evap (symmetry, BC-side rows, transposed g5 slice) | b1/b6 unchanged under `prescribed`/legacy toggles; `atm/` module unit-covered; mass audits extended (evap shortfall); **p1/p5 green for every new kernel** (V2-A6) |
| Q5 | g6, g7 (g8 recommended); s-matrix green; x-battery for heat (symmetry incl. thermal BC on all sides) | b6 unchanged with temperature off; g5 cross-check (§4.3); density-coefficient config migration is golden-neutral |
| Q6 | g9, g10, Cd battery; s-matrix green; x-battery for wind (8-orientation setup case is itself the symmetry gate) | b1–b6 unchanged (wind off); direction-wrap unit tests |
| Q7 | all g + all b + s1–s4 + p1–p5 in one pipeline; §8.2 BC matrix has zero "unverified" rows; §8.3 table has zero empty cells and every GPU row is either `supported` (p6 returned) or explicitly `experimental` | sanitizer matrix over new labels; docs (theory chapters for atm/heat/wind; symmetry-exemption table; parameter-table lockstep); v2.0.0 |

## 10. Risk register (additions to the v1 register)

| Risk | Mitigation |
|---|---|
| AMG setup cost erases iteration savings at benchmark scale | Reuse policy (§2.2.4) + g3 setup-fraction telemetry; keep bjacobi-icc default in v2.0 |
| BoomerAMG complexity blowup on the anisotropic gw system | strong-threshold sweep at bring-up; complexity ≤ 2.0 asserted in g3 |
| g5 digitization error dominates tolerance | Profile-RMS metrics on horizontally averaged fields only; finger positions never gated pointwise |
| Two-scalar refactor perturbs b6 bitwise restart | Restart-determinism lanes extended to the two-scalar checkpoint layout before physics lands (P0-style: infrastructure first) |
| Salinity-limiter change shifts b6 | `legacy_evap_allowance` toggle pins b6; new exact factor gated by g4/g5 only |
| Met-forcing module scope creep (radiation schemes, cloud models) | Q5 ships bulk formulas + prescribed radiation series only; anything more is an amendment |
| s-gate thresholds unreliable on the fanless dev machine (the M3 lottery, A23) | min-over-repeats everywhere; hard gates only on iteration counts and correctness, soft (warn) on timing until real runners calibrate |
| Symmetry battery flags legacy-faithful asymmetries as failures | documented-asymmetry exemption table (§8.1) with per-term provenance; building it doubles as an owner-reviewed audit |
| BC-matrix backfill (§8.2) reveals v1 silent no-ops mid-Q1 | expected (y+-only scalar injection is known); fix-or-reject per row by amendment, never silently accept |
| Coverage machinery (x-gates) grows slower than features | matrix/table checks are CI-enforced lockstep scripts (like schema↔docs), so an uncovered feature fails the build rather than relying on review vigilance |
| **GPU lane ships without ever executing** (Q3) | p4 compile+link + p5 static/compile-time/runtime pointer discipline + p1 proving backend-swap neutrality on CPU; `experimental` status printed at startup and in release notes; p6 bundle converts it to `supported` with one owner run. The honest residual: GPU *performance* is unverified at v2.0 |
| Threaded hypre/PETSc underdelivers on CPU (ecosystem, not ours) | p2 hard-gates only the falsifiable minimum (iterations invariant, solve time decreasing 1→4 threads) and *records* efficiency; MPI stays the recommended and default CPU lane, so a weak threaded AMG degrades an option rather than the product |
| `aijkokkos` on CPU is slower than `aij` at benchmark scale | Expected and documented (PETSc's own guidance); the lane-selection table §2B.2 tells users which lane to pick, and `aij` remains the default — p1 guarantees only that switching is *correct*, not that it is always faster |
| Q3 slips and delays the physics phases | Q3 is strictly infrastructure with no physics coupling; if it slips, Q4/Q6 can start in parallel at the cost of retro-fitting p1/p5 to their kernels (the cost this ordering exists to avoid) — an explicit amendment, not a silent reordering |

## 11. Deferred (recorded now, decided later)

Vegetation transpiration / root-zone uptake (natural Q4 successor — Feddes-type sink in
the Richards corrector, gated on the Simunek/HYDRUS root-uptake verification set);
Coriolis + radiation open boundaries (the "coastal credibility" pair); NetCDF/CF or
XDMF output; flipping the default preconditioner to AMG (needs g3 history on real
hardware); Gran et al. (2011, *HESS*) saline-soil column drying experiment as
*extended validation* for Q4 (lab data, but its mineral precipitation is out of scope
for a gate).

*Removed from this list by V2-A6:* "GPU execution of the AMG path" and the
standing GPU-validation debt are no longer deferred — they are Q3 (§2B), with the
single genuinely-deferred residual being **on-hardware GPU performance
verification**, which is the p6 owner action.

## 12. Key references by gate

**Q1 (solver):** Jones & Woodward 2001, *Adv. Water Resour.* 24:763 (multigrid-preconditioned
Richards); Kollet et al. 2010, *WRR* 46:W04201 (weak scaling to 16k cores — the g3 protocol);
Maxwell 2013, *Adv. Water Resour.* 53:109; hypre BoomerAMG docs + MOOSE hypre tuning guide
(3D strong-threshold ≥ 0.5 warning); PETSc KSP manual (AMG-on-anisotropy caveats);
`KSPSetReusePreconditioner` manual page; Demidov 2021, arXiv:2108.02054 (setup-reuse
amortization); MOOSE issue #8681 (asymmetric Dirichlet rows stall CG+AMG).

**Q3 (performance portability):** hypre GPU wiki, github.com/hypre-space/hypre/wiki/GPUs
(**authoritative device-support list**: PMIS the only GPU coarsening; device interpolation
= direct/BAMG-direct/ext/ext+i/ext+e; multipass has no GPU implementation; device smoothers
Jacobi/two-stage GS/l1-Jacobi/Chebyshev); Falgout et al., *Porting hypre to Heterogeneous
Computer Architectures*, OSTI 1860740 (two-stage MM-ext/MM-ext+e interpolation; aggressive
coarsening less effective on GPU); hypre BoomerAMG manual (backend-dependent coarsening
defaults: HMIS on CPU, PMIS on GPU); PETSc installation guide §Kokkos
(`--download-kokkos --download-kokkos-kernels` + one of `--with-cuda|hip|sycl|openmp`;
`--with-openmp` propagates to downloaded packages; `--with-openmp-kernels` threads only
*some* PETSc routines; oversubscription warning for hybrid runs); PETSc 3.14 release notes
(`-vec_type kokkos -mat_type aijkokkos` introduction and initial MATAIJKOKKOS operation
coverage); PETSc `VecGetKokkosView`/`PetscGetKokkosExecutionSpace` manual pages and
`src/snes/tutorials/ex3k.kokkos.cxx`; Mills et al. 2024, arXiv:2406.08646 (*PETSc/TAO
Developments for GPU-Based Early Exascale Systems*); Lange et al. 2013, arXiv:1303.5275
(*Achieving Efficient Strong Scaling with PETSc using Hybrid MPI/OpenMP Optimisation* —
the 4-ranks×8-threads/one-rank-per-NUMA-domain heuristic p3 tests); petsc-users 2015-05
thread on verifying an OpenMP-linked hypre (`otool -L`/`ldd` check — the silent
serial-hypre trap B3 avoids by configuring in one shot).

**Q4 (evaporation):** Geng & Boufadel 2015, *J. Hydrology* 524:427-438 (g5 — the PDF is at
repo root); Mahfouf & Noilhan 1991 (bulk aerodynamic form); Lee & Pielke 1992 / Barton 1979
(α₁ soil humidity); Tetens 1930 (e_sat); Liu et al. 2006 (R_air fit); USGS Prof. Papers
269/270 (Lake Hefner mass-transfer coefficients — Tier-1 formula unit tests; optional
Tier-2 monthly total after one-time digitization); GLM 3.0 (Hipsey et al. 2019, *GMD*
12:473) for the still-air transfer-function lower bound.

**Q5 (temperature):** Bredehoeft & Papadopulos 1965, *WRR* 1:325; Ogata & Banks 1961,
USGS PP 411-A with the OpenGeoSys heat-benchmark parameterization; Stallman 1965, *JGR*
70:2821 (closed form as transcribed in EPA/600/R-15/454 Eq. 6, with Hatch 2006 amplitude/
phase inversions); Edinger, Duttweiler & Geyer 1968, *WRR* 4:1137 (equilibrium temperature;
CE-QUAL-W2 App. A operationalization); van Genuchten et al. 2013, *J. Hydrol. Hydromech.*
61(2):146 & 61(3):250 (exact river ADE+exchange solutions; open-access USDA PDFs) /
Neilson et al. 2012, *Adv. Water Resour.* (heat-native form); Horton & Rogers 1945 +
Lapwood 1948 (Ra_c = 4π²); SEAWAT v4 TM 6-A22 (thermal-retardation-as-scalar precedent);
*Fluids* 2017 2(1):11 (Elder multiple steady states — why Elder is not a gate).

**Q6 (wind):** TELEMAC-2D Validation Document v5.0 §3.10 (wind set-up in a closed basin —
config and 81.52 cm analytic setup); Dean & Dalrymple 1991 Eq. 5.96; Merian formula per
Rabinovich 2010 (*Handbook of Coastal and Ocean Engineering* ch. 9); Garratt 1977 (ADCIRC
default Cd with 3.5e-3 cap), Smith & Banke 1975 (Delft3D default breakpoints), Wu 1982,
Large & Pond 1981, Powell et al. 2003 (high-wind saturation); Kraus & Militello 1999
(DTIC ADA482229, optional periodic-wind analytic).

---

## Amendment log (v2)

### V2-A1 (Q1, 2026-09-09) — §2.3 g2: the synthetic case's fs system is trivial

Measurement: on the A23 synthetic scaling case the free-surface solves
converge in 0 iterations at every rank count (the draining-box surface
leaves the eta system's residual below atol), so a recorded-count gate on
fs would pin 0 and measure nothing. g2's recorded-count assertion therefore
binds the gw system only (`tolerances/g2.yaml` records `gw_mean`); the
ratio assertion still evaluates both systems, and the fs system's AMG
behavior is exercised by the g1 lanes (b1/b4 run 3600 fs solves per case
with mean ~2 iterations under amg).

### V2-A2 (Q1, 2026-09-09) — §2.2.4: divergence handling under hierarchy reuse

The plan specified the rebuild triggers but not divergence behavior. As
implemented: a solve that diverges under a *reused* hierarchy rebuilds the
preconditioner and retries once (counted in telemetry as a retry) before
the fatal-divergence rule applies; divergence under a freshly built
preconditioner stays immediately fatal (the v1 rule). Rationale: a stale
frozen hierarchy is an artifact of the reuse optimization, not of the
physics, so one refresh attempt is owed before aborting; every branch is
KSP-collective so all ranks agree.

### V2-A3 (Q1, 2026-09-09) — §2.3 g2: E-core ranks are valid for iteration gates

g2 runs at 8 ranks on the 4P+4E reference machine. A23's caveat (ranks
beyond the P-core count are a topology artifact) applies to *timing*
gates; iteration counts are decomposition-determined and
scheduling-independent, so the 8-rank iteration measurement is valid. The
timing-bearing s-gates keep the A23 discipline (min-over-repeats, 4-rank
efficiency criteria).

### V2-A4 (Q1, 2026-09-10) — §2.3 g1: the b5 amg lanes are held to the parent's approved record

Two adjudications from the first nightly g1 execution, both preserving the
"solver choice must not change the physics *criteria*" principle while
correcting what those criteria actually are:

1. **Horizon.** The g1 b5 rain/sync amg lane initially ran the full
   ~120 h horizon and failed one metric: ponding integral −11.7 % vs the
   10 % allowance over [2.3, 69.5] h (discharge 0/84 outside, ponding
   0/81 outside, peaks within bounds — everything else green). But the
   parent gate's *approved record* is the A15/A17 shortened 24 h horizon
   (owner adjudication at P3; `scripts/ci_release_gate.sh` gates b5 with
   `--t-end 86400`), and the full-horizon envelope was never adjudicated
   green under any solver. The g1 lane now runs `--t-end 86400`, matching
   the record it is invariance-testing against. The full-horizon amg
   metrics are archived here for the record.
2. **Default-mode rank-invariance bounds are per-solver.** The A14 bounds
   (subsurface 3e-5, exchanged/surface 2e-3) were *data-derived under
   bjacobi-icc*. Measured at Q1: switching solver at a fixed 1-rank
   decomposition already moves the subsurface volume 5.9e-5 and the
   exchange 4.4e-3 — past those bounds with no rank change at all — and
   AMG coarsening is decomposition-dependent by construction, so amg
   rank-to-rank drift (measured 1.8e-4 / 8.1e-3 / 1.4e-3 at 600 s) is
   the A14 chaos-amplification class on a different solve trajectory,
   not a defect. Per the A5/A8/A14 data-derived-bound discipline, the
   amg/gamg default-mode bounds are set at ~2x the Q1 measurements:
   subsurface 4e-4, exchanged 1.6e-2, surface 3e-3
   (`run_regression.py`, solver-keyed). The strict lane (jacobi CLI
   override, 1e-12) is unchanged and remains the machine-precision
   invariance proof.

### V2-A5 (2026-09-10) — §1: phase inserted; labels renumbered, section numbers frozen

Owner-requested capability (from the first hands-on Kuan runs in frehg2-dev):
every simulation must leave a persistent record of its settings, active
modules, boundary conditions, and per-segment run times. Inserted as the new
**Q2** (§2A) — pure infrastructure, one week, placed before the physics phases
so their gate runs are self-recording from the start. The former Q2–Q5
(evaporation, temperature, wind, release) are relabeled **Q3–Q6**; gate names
(g4–g10) are unchanged. Document section numbers §2–§12 are frozen because
source code and the Q1 records cite them; the new phase is §2A.

### V2-A6 (2026-09-11) — §1: performance-portability phase inserted as Q3 (§2B)

Owner-requested capability, from the Q2 run-record analysis of the Kuan runs
(omp2 slower than serial) and the follow-up question of whether one source tree
can serve laptop CPU, HPC CPU, single GPU and multi-GPU. Diagnosis: it can, and
most of the abstraction work is already done (`MemSpace`-templated fields, the
`check_forbidden.sh` backend-isolation rules, COO assembly, a GPU-aware halo
exchanger) — but the **linear-algebra boundary** is neither threaded nor
device-correct: `MATAIJ` is hardcoded, the solve stages RHS/solution through host
mirrors, Q1's AMG defaults use host-only hypre options (HMIS coarsening, unguarded
`agg_nl`), and no build path configures PETSc with Kokkos. The `MatSetValuesCOO`
call is a **latent device defect** today, invisible to every CPU lane because
`MemSpace == HostSpace` under OpenMP.

Inserted as the new **Q3** (§2B) with gates **p1–p6**; the former Q3–Q6
(evaporation, temperature, wind, release) are relabeled **Q4–Q7**. Gate names
g4–g10 and r1/r2 are unchanged. Document section numbers §3–§12 remain frozen;
the new phase is §2B (precedent: §2A).

Three consequences recorded here because they change earlier text:

1. **§7.2 s3 amended.** Its "PETSc solve time asserted flat under threads"
   assertion was a true statement about the MPI-only `aij` design; it now applies
   to the `aij` lane only. Under `aijkokkos` the solve must *speed up*, and p2
   asserts it. Left unamended, s3 would gate in the limitation Q3 removes.
2. **§8.3 gains an `experimental` status.** A GPU lane that has never executed is
   exactly the "authored-unexercised" state §8.3 forbids in a release (P5's
   compile-only CUDA lane is the precedent). Rather than weaken the rule, device
   lanes ship as `experimental` — compile+link verified, statically
   invariant-checked, physics-verified on the equivalent CPU backend, not executed
   on target hardware — printed at startup and in the release notes, cleared to
   `supported` when the owner returns a passing p6 bundle.
3. **§11 deferred list.** "GPU execution of the AMG path" and the GPU-validation
   debt move out of deferred into Q3. The residual deferred item is on-hardware
   GPU *performance*, owned by p6.

Ordering rationale: Q3 precedes the physics phases for the same reason Q2 did.
Every kernel authored in Q4–Q6 must be device-correct from its first PR, and
p1/p5 are added to Q4–Q6's exit criteria; retro-fitting memory-space discipline
onto three shipped physics modules is strictly more expensive.

Scope boundary: Q3 is infrastructure only — no physics, no numerics changes. Any
change to a discretization or coefficient discovered necessary during Q3 is a
separate amendment.

### V2-A7 (Q3, 2026-09-11) — §2B.2 B2: l1-scaled Jacobi smoother on the Kokkos lanes

**What changed.** §2B.2 B2 as drafted said BoomerAMG relaxation "stays PETSc's
backend default — symmetric SOR/Jacobi on host, l1-Jacobi on device." At
implementation the host default proved incompatible with p2's thread-invariance
assertion: hypre's *hybrid* symmetric SOR/Jacobi is thread-count-DEPENDENT by
construction (Jacobi between OpenMP threads, Gauss-Seidel within each thread's
rows), so gw iteration counts drift with OMP_NUM_THREADS — measured 8 → 9
between smoother partitionings on the Kuan-scale smoke case, and any drift
fails p2's ±2 % bound by design, not by defect.

**Resolution.** When the resolved vector type is Kokkos (`kokkosVec_`), the amg
path defaults `pc_hypre_boomeramg_relax_type_all = l1scaled-Jacobi`: symmetric
(CG-safe), thread-count-invariant, and the same smoother family PETSc selects
on device — so the OpenMP lane rehearses the exact GPU algebra, which is p1's
purpose. Measured cost on the 2.2M-cell p2 case: gw 15 iterations/solve vs 8-9
under the host default (~1.7×), repaid by the solve now threading (34.6 s →
20.4 s at 4 threads) and identical at 1/2/4/8 threads. The host `aij` lane
keeps hypre's default — the Q1 g2 record and the golden-pinned iteration
history are untouched. Still a `setOptionDefault`: any options file or
command line wins.

**Scope-boundary note.** This is a preconditioner-lane default, not a
discretization or coefficient change; solutions remain within the established
rank-invariance tolerances (p1 asserts it). Logged as an amendment because the
B2 text promised backend defaults and the implementation deliberately does not
follow it on the Kokkos host lane.

### V2-A8 (Q3, 2026-09-11) — §2B.2 B3: Kokkos must be built shared, and a build-provenance gate to prove it

**What changed.** §2B.2 B3 required the self-contained build scripts
(`build_frehg2_local.sh`, `build_frehg2_hpc.sh`) to install a Kokkos-aware
PETSc; it did not specify Kokkos's link type. Both scripts built Kokkos with
CMake's default (`BUILD_SHARED_LIBS=OFF` → static `libkokkoscore.a`). This is
silently wrong when the Kokkos runtime is shared with PETSc: frehg links
Kokkos core, and so does PETSc's `--download-kokkos-kernels`
(`libkokkoskernels.dylib`). A static core is absorbed *independently* into the
frehg executable and into `libkokkoskernels.dylib`, so the process holds two
copies of Kokkos's global singleton — `Kokkos::OpenMP::initialize` runs twice
(and prints its banner twice), and the first `VecKokkos` view access
dereferences state owned by the other runtime and **segfaults**.

**How it escaped the p-gates.** p1/p2/p3/p6 all validated the *code*, but every
one of them ran against a Kokkos that predated the self-contained scripts — the
hand-built `local/` prefix, whose Kokkos happened to be **shared**
(`libkokkoscore.dylib`), so it carried a single runtime and never reproduced
the fault. No gate ever built the dependency stack *from* `build_frehg2_local.sh`
and ran a solve. The defect surfaced only when a fresh self-contained build
(the `frehg2-dev` deployment copy) segfaulted on the b6-kuan-td case. The bug
was in packaging, not in any gated artifact — which is exactly the blind spot.

**Resolution.** Add `-DBUILD_SHARED_LIBS=ON` to the Kokkos install in both
scripts. This is precisely the configuration of the `local/` prefix that passed
every p-gate, so the fix is validated by the known-good build, not merely
asserted. A rebuild against shared Kokkos must also rebuild PETSc, because its
`kokkos-kernels` embedded the static core; the scripts' reuse markers make this
a delete-Kokkos-and-PETSc-then-rerun, keeping MPICH/HDF5/yaml-cpp. Because
Kokkos and PETSc have *independent* reuse markers, deleting only Kokkos leaves a
stale PETSc whose `libkokkoskernels` was removed with the old Kokkos — the
script now guards this (when reusing PETSc it asserts the Kokkos Kernels runtime
is present in the prefix, since `petscconf.h` still advertises it after the
library is gone), so the mismatch fails with an actionable message instead of a
dyld abort at the first-run validate step.

**New gate (p4 extension — build provenance).** p4 previously covered only the
device compile+link in CI. It is extended with a **host self-contained-build
smoke gate**: build the full stack from `build_frehg2_local.sh` into a scratch
prefix, then assert (a) `otool -L`/`ldd` on the frehg binary resolves a *shared*
`libkokkoscore` (no statically-absorbed core), (b) a one-step `aijkokkos`+amg
solve runs to completion and emits the Kokkos init banner exactly once, and
(c) the result matches the `aij` lane to rank-invariance tolerance. This closes
the "gates ran against a different dependency build than ships" gap by making
the shipped build script itself a gated artifact.

**Scope-boundary note.** No source or discretization change; this is a
dependency-packaging correction plus the gate that would have caught it.
Physics is untouched — the same binary, built correctly, reproduces the p1
records.
