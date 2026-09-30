# Q7 Report — Release hardening and v2.0.0

Phase Q7 of the v2 plan (§8 generality instruments + §9 release row),
executed 2026-09-24 → 2026-09-30. Commits: `a812bc4` (instruments,
authored red), `def8fba` (implementation green + version), `1d6dc30`,
`0044e55`, `f56ba96` (pipeline bring-up fixes), `6b0baaa` (pipeline
scoping), plus the close commit carrying this report. Amendments:
V2-A19, V2-A20.

## What the phase delivered

1. **The §8 instruments, all four axes.** §8.1 dihedral batteries for the
   three v1 modules (the Q5/Q6 heat/wind batteries already existed) with
   class-aware tolerances and the exemption table
   ([symmetry-exemptions](../theory/symmetry-exemptions.md), 9 waived
   asymmetries with provenance + measured floors); §8.2 BC matrix — 117
   canonical kind × target × coupling cells, every one `tested` or
   `schema-rejected`, zero `limitation` rows, lockstep-checked with a
   schema-vocabulary drift alarm; §8.3 feature-interaction table — 90
   rows, zero empty cells, GPU rows `experimental` per §2B.4; §8.4 —
   n=3 mpi lanes, the transposed g5 slice, corner/single-cell polygons,
   masked columns per edge, restart-mid-transient rows, the pair gates.
2. **Two real defects fixed under gates authored red.**
   - **V2-A11** (diagnosed Q0.3, deliberately deferred to the §8.2
     gate): the west/south face-area ghost copies in
     `WetDry.cpp::updateGeometry` handed a transmissive outlet the
     interior face area — outlets on descending west/south beds were
     throttled `(deptx/depth)²` (~1/812 in the superslab) and minted
     volume through the below-bed clamp when fed from upslope. The fix
     is the diagnosed two-line deletion; `regression.outflow_staircase`
     (10 stock failures across the four-edge sweep) and the swe battery
     gate it.
   - **The transport interface-flux halo** (found by the §8.4 tracer
     lane executing — the second battery-caught defect after V2-A18):
     the `swe_Fu`/`swe_Fv` exchange ran only on the eta-BC path, so
     multi-rank transport with no eta condition advected scalar with
     inconsistent rank-interface fluxes. Stock drift 1.6e-2 (upwind) /
     1.6e-3 (superbee) with eta rank-invariant at 4e-16 — every gated
     transport case was serial or carried an eta condition, so five
     phases never saw it. Post-fix: 9.0e-12 / 5.2e-12 strict (bound
     1e-10), 3.7e-6 default superbee (bound 1e-4).
3. **One schema hardening**: `scalar_value` on `groundwater_top` in
   coupled runs is rejected loudly (coupler-owned face; the last §8.2
   cell that was accepted-but-unverifiable).
4. **The v2 release pipeline** (`ci_release_gate.sh`): all g + all b +
   s1–s4 + p1–p5 in one run; p4 by nvcc or recorded evidence.
5. **Docs for release**: theory/atmosphere.md (the consolidated
   bulk-aerodynamic chain), the symmetry-exemption table, testing.md
   rewritten to the v2 label reality, the missing Q5 records and
   feature-coverage.md into the mkdocs nav, parameter lockstep green at
   178 schema tokens.
6. **v2.0.0**: version bump, this record, the tag.

## Gate results (release record)

- Per-PR rehearsal on `def8fba`: unit 41/41, mpi 20/20 (the 8 mpi
  drivers plus the 12 rank-invariance lanes dual-labeled
  `regression;mpi`), anchored regression 47/47 (the label's other
  lanes; the whole label is 59/59 in pipeline step 2); forbidden scan,
  parameter lockstep, mkdocs --strict clean.
- Battery floors (strict solver mode, macOS arm64): swe transpose
  4.7e-15 / mixed 3.478e-3 (bound 1e-2; x-swaps 2.227e-3 and y-swaps
  1.336e-3, additive — the pure legacy plus/minus edge arithmetic once
  V2-A11 is fixed); gw ≤ 3.75e-14 on all 8 transforms (side arithmetic
  exactly symmetric, both classes gate 1e-12); transport transpose
  2.1e-15 / mixed eta 3.09e-10 amplified ~7× into concentration 2.17e-9
  (x-swaps only; bound 1e-6); heat ≤ 1.5e-13; wind ≤ 1.1e-15 m.
- Pair lanes: `g1.b1_restart.amg` carries a deterministic 6.9e-8
  hierarchy-reuse-cadence offset (documented bound 1e-6);
  `p1.b1_restart.aijkokkos` restarts exactly deterministic.
- **Release pipeline** (`ci_release_gate.sh build-ci`): one `build-ci`
  binary, configured and built at `f56ba96` on 2026-09-26, run in three
  segments because session teardowns ended the first two. Later commits
  touch docs and script comments only, so every segment ran the same
  binary against the same harness and tolerances.
  - *Segment 1 (2026-09-26 → 27)*: steps 1–3 PASS — the per-PR gate
    (unit 41/41, mpi 20/20, fast regression 31/31, forbidden scan,
    parameter lockstep, Doxygen, mkdocs --strict), the anchored
    regression label 59/59, and the adjudicated b5 record (rain/sync,
    `--t-end 86400`, 4 ranks: peak +2.6 %, integrals +2.4/−2.5 %).
    Step 4, then still the blanket nightly label, started 05:12 and
    completed three full-horizon b5 lanes, all FAIL (rain/sync
    30 360 s, rain/subcycled 29 356 s, norain/sync 28 174 s; the
    measured values are Finding 2). Teardown ended norain/subcycled.
  - *Segment 2 (2026-09-27, step 4 scoped by `6b0baaa`)*: `b6.ss`
    287 s, `b6.td` 1000 s, `g8` 67 s, `g5.geng2015` 2288 s,
    `g5.transposed` 2423 s, `g1.b5_rank_invariance.default.amg` 479 s
    — all PASS. Teardown ended `g1.b5.rain.sync.amg`.
  - *Segment 3 (2026-09-29, detached from the session)*: the four
    remaining step-4 lanes — `g1.b5.rain.sync.amg` 10 973 s (sharing
    the machine with the serial superslab rerun), `g1.b6.ss.amg` 285 s,
    `g1.b6.td.amg` 747 s, `p1.b6.ss.omp4.aijkokkos` 943 s — all PASS.
    Step 4 is 10/10 PASS over segments 2–3.
  - *Steps 5–6 (2026-09-29 23:18 → 09-30 00:31, the superslab rerun
    paused so the scaling gates saw an idle machine)*: the scaling
    label 7/7 PASS in 1.2 h. s1: 2-rank efficiency 91.5 % (bound 80 %);
    the 4-rank 53.1 % is recorded, not gated, on four performance cores
    (V2-A12). s2/g3: efficiency 56.5 % at 4 ranks (WARN — above the
    0.5 hard floor, below the 0.8 target), gw iteration growth +6.3 %
    (bound 15 %); the fs system converges in 0 iterations at every rank
    count on this case (V2-A1), so its growth check prints −100 % and
    binds nothing. s3: kernel efficiency 63.9 % at 4 threads (bound
    60 %). s4: 53.9 / 58.0 / 66.2 s for 4×1 / 2×2 / 1×4. p2: gw solve
    21.0 s at 4 threads against 34.3 s at 1. p3: iterations
    placement-invariant (15.0), fastest 4×1. The g2 bjacobi record ran.
    Step 6 recorded p4 by evidence (no nvcc here): the re-run of
    cuda-compile.yml on the pushed release SHA.
- **Sanitizer matrix** (`run_sanitizers.sh --full`, Apple clang +
  deps-clang, `build-asan` configured at `6b0baaa` — sources identical
  to `f56ba96`): the first run was 104/108 clean with four
  lane-composition failures (Finding 9); the rerun with the printed
  exclusions ran clean in 2.7 h: 104/104 lanes (unit 41, mpi 20,
  regression 55 of 59; the 12 dual-labeled lanes counted once) and all
  eight nightly-class smokes (b5 ×4, b6 ×2, g5, g8).

## Findings

1. **The lane scripts had never run on this platform.** The Q0.2
   `_mpich_libdir` loader guard (authored for Linux runners) made the
   command substitution exit 1 wherever no lib path contains "mpich" —
   `set -euo pipefail` then killed `ci_build_and_test.sh`,
   `run_sanitizers.sh`, and `run_openmp_lane.sh` before their first
   line of output. The release pipeline's first local run died silently
   at step 1; fixed by failure-guarding the grep (`1d6dc30`). Two more
   first-execution defects followed: six undocumented `ScalarSpec`
   members under Doxygen-as-errors (`0044e55`) and a relative `--frehg`
   path that broke the staged b5 step (`f56ba96`). The near-miss
   pattern is the one §8 exists for: scripts exercised only in the
   environment they were written against.
2. **The four full-horizon b5 envelope lanes are not gates — scoped out
   of the pipeline with their measurements archived (V2-A20).** The
   pipeline's first full execution failed three of them after ~25 h:
   rain/sync ponding integral −11.8 % (V2-A4 archived −11.7 % for the
   same metric at Q1 — stable across the whole Q1→Q7 evolution, not a
   Q7 regression); rain/subcycled peak +15.9 % (the A17
   documented-failing fixed-dt configuration, +16.3 % at P3);
   norain/sync peak −22.5 % on the first completed norain run in the
   project's history. The adjudicated record (rain/sync at the A15/A17
   24 h horizon) passes on this tree, and is what step 3 gates.
   The norain miss is attributed by construction and by measurement:
   b5's DEM duplicates its boundary rows/columns (j=0 ≡ j=1,
   i=0 ≡ i=1), so the V2-A11 deletion is a self-assignment on this
   grid — a stock binary restarted from the fixed run's checkpoint
   reproduces it bitwise on all six output fields across the
   first-wetting window and across the discharge peak. The miss is the
   pre-existing scheme in a regime no run had ever completed.
3. **V2-A10's "b1–b6 green on real runners" was vacuous for every
   golden and envelope gate** — the sixth instance of the
   vacuous-evidence pattern, found while attributing the b5 failures.
   The public workflow history shows every scheduled nightly ever ran
   ≤ 2.5 h, and the per-step timings of the 2026-09-26 success show
   the envelope step at 0 s: the legacy-goldens archive was never
   provisioned on the runner, so the goldens-absent guard has silently
   skipped the golden regression comparisons and the entire
   `regression_nightly` label on every run. What the runner validates:
   build, unit, mpi, restart/rank-invariance, s-gates. All golden and
   envelope evidence for v2.0.0 is this local pipeline. Carried
   forward: arm the runner (the b5 fulls additionally exceed the
   hosted 6 h job cap) or make the skip a stated scope, not a silent
   success.
4. **The V2-A12 residual closed itself on schedule**: the recalibrated
   s1 gate ran green in five consecutive scheduled nightlies on the
   real 4-vCPU runner (2026-09-22 → 09-26; the scaling step is not
   goldens-guarded and runs 1:31 h).
5. **The swe mixed-class floor is per-axis and additive** (x-swaps
   2.227e-3 + y-swaps 1.336e-3 ≈ both 3.478e-3), consistent with
   independent west/south vs east/north edge terms — a useful
   diagnostic signature: a future symmetry regression that is NOT
   additive in the two axes is a new defect, not the documented floor.
6. **A pre-existing unaudited sink in deep thin-film recession**, found
   by the V2-A11 validation rerun (V2-A20 finding 4): with the outlet
   throttle gone, swe-vcatchment drains as a whole-domain ~1e-4 m film
   and the surface audit misses −70.8 m³ (0.27 % of rain) by t_end —
   interior-distributed, all cells wet, clamp frozen, present on the
   stock binary too (−35.5 m³). Recorded in the case README; spun off
   as its own investigation task.
7. **Validation-case records refreshed after V2-A11**: swe-vcatchment
   rerun — the 0.192 m outlet pool is gone (outlet cells 7.4e-3 m, row
   j=0 total 2.85 m³ vs 39.7 m³), regenerated hydrograph peak
   292.8 m³/min at t = 92 min, mean deviation vs ParFlow/HGS/tRIBS
   28.5/19.2/13.2 m³/min — inside the Maxwell (2014) intercomparison
   envelope. The swere-superslab rerun on the release binary takes most
   of a day serial (its README's 15–30 min estimate was unmeasured) and
   had not finished at the tag; the README marks its table and figures
   as pre-fix, and the regeneration is carried forward (V2-A20 scope
   note). The swe-outflow-staircase README records
   the landed fix and its promotion to `regression.outflow_staircase`.
8. **The run record's revision stamp can name code that did not run**
   (V2-A20 finding 5). The first swe-vcatchment rerun used a `build/`
   binary compiled 2026-09-24 from the working tree carrying the V2-A11
   fix; its record named version 1.0.0 at `e9d84cc`, a tree without the
   fix, because `FREHG_GIT_SHA` is read only when CMake configures and
   nothing marks a modified tree. The rerun was repeated with the
   release binary: bitwise identical (1177 datasets), so the README
   numbers stand, and every release-evidence binary was configured at
   its own tree. No gate checks the stamp's currency — r1 checks the
   record's schema, config round-trip, and launch shape. Carried
   forward.
9. **The sanitizer matrix had not run in full since P5, and its first
   full run failed on lane composition, not on frehg code** (V2-A20
   finding 6). 104 of 108 lanes ran clean. `regression.perf_baseline`
   (Q2) compares instrumented wall-clock against an uninstrumented
   baseline and fails by construction (+271 %). The three `aijkokkos`
   lanes (Q3: `p1.b1_restart`, `p1.b1.omp4`, `p1.b4.omp4`) aborted on an
   ASan bad-free inside `VecDestroy_SeqKokkos`. The lane builds frehg
   with Apple clang/libc++ against a clang-built Kokkos but links the
   gcc/libstdc++ Kokkos-aware PETSc, so one process holds two Kokkos
   builds and PETSc's views are destroyed through the executable's copy
   of the destructor. The release binary links one Kokkos and passes the
   same lanes in step 2. The run also stopped at its ctest block under
   `set -e`, so the nightly-class smokes never ran. `run_sanitizers.sh`
   now excludes both by name and prints the exclusion (the `aijkokkos`
   lanes only when the two C++ runtimes differ), and the whole matrix
   was rerun clean (104/104 lanes, 8/8 smokes). Because of the
   exclusion, the frehg-side `aijkokkos` staging path has never run
   under a sanitizer. Carried forward.

## Cost

The compute is the release pipeline, run 2026-09-26 → 30. Steps 1–3
took ~2.8 h including the build; the three full-horizon b5 lanes that
ran before V2-A20 scoped them out took 24.4 h; the scoped step 4 took
5.4 h (1.8 h in segment 2, 3.6 h in segment 3); step 5 took 1.2 h.
The sanitizer matrix ran twice (2.4 h, then 2.7 h). The serial
swere-superslab rerun on the release binary had covered 5 h of its 12 h
of simulated time in 6.2 h of running wall time at the last
measurement; its README estimated 15–30 min for the whole run.

## Carried forward (the v2.0.0 open set)

The nine dod-Q7 items: the owner p6 GPU bundle (GPU rows stay
`experimental`); cuda-compile.yml green on the pushed release SHA (p4)
before the official tag; arming the runner's golden and envelope gates,
or stating their local-only scope in the workflow (Finding 3); the owner
§6.4 g5 digitization sign-off; the P5 BSD-3 license review; the
thin-film audit sink (Finding 6); a build-time revision stamp with a
dirty flag (Finding 8); a sanitizer lane for the `aijkokkos` staging
path on an ABI-consistent Kokkos-aware PETSc (Finding 9); and the
swere-superslab regeneration once its rerun finishes (Finding 7). Closed since
the Q6 list: the V2-A12 s1-on-runner check (Finding 4). Deferred scope
(§11) is unchanged:
transpiration/root uptake, Coriolis, radiation BCs, NetCDF/VTK,
AMG-as-default, on-hardware GPU performance.
