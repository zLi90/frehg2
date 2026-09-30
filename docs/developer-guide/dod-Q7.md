# Definition of Done — Q7 (Release hardening, v2.0.0)

Per v2 development plan §8/§9 (Q7 row) and amendments V2-A19/V2-A20,
under the v1 plan's §11.4 DoD conventions. Every item names its
verification; commands were run green 2026-09-26 → 30 on macOS arm64
(gcc-16.2, MPICH 4.3, PETSc 3.25.1 + hypre 3.1.0, Kokkos 5.1.1 OpenMP
host backend; `FI_PROVIDER=tcp`, `UCX_TLS=tcp,self,sm`).

## Blocking (§9 Q7 row)

- [x] **One pipeline — all g + all b + s1–s4 + p1–p5**
      (`scripts/ci_release_gate.sh build-ci`, V2-A19 scope, step 4
      scoped per V2-A20): PASS, executed 2026-09-26 → 29 in three
      segments against one `build-ci` binary (configured and built at
      `f56ba96` on 2026-09-26; nothing the binary or the harness reads
      changed afterwards — later commits touch docs and script comments
      only). Session teardowns ended the first two segments: the first
      inside the four full-horizon b5 runs that V2-A20 then removed
      from the pipeline, the second after six of the ten scoped step-4
      lanes; each resumption ran only the lanes not yet executed
      (report-Q7.md lists every lane with its segment and time).
      Steps: the per-PR gate (strict build; unit 41/41 incl. the §8.2/
      §8.3 lockstep checkers; the 8 mpi drivers at n = 1/2/3/4; the
      fast regression subset 31/31; forbidden scan; parameter
      lockstep; Doxygen; mkdocs --strict),
      the full anchored regression label (59/59 — batteries, side
      sweeps, pair gates, g4/g6/g7/g9/g10, g1 amg/gamg lanes, p1
      lanes, g2, r2, perf baseline), the adjudicated b5 record
      (rain/sync, `--t-end 86400`, 4 ranks: PASS — peak +2.6 %,
      integrals +2.4/−2.5 %), the nightly regression label minus the
      never-adjudicated b5 full-horizon lanes (10/10 PASS — b6 ss/td,
      g8, g5 geng2015/transposed, the g1 AMG lanes for b5 rank
      invariance, b5 rain/sync, and b6 ss/td, and the p1 b6.ss
      aijkokkos lane), the
      nightly scaling label (7/7 PASS — s1 2-rank efficiency 91.5 %
      against 80 %, the 4-rank point recorded not gated on four
      performance cores per V2-A12; s2/g3 gw iteration growth +6.3 %
      with efficiency 56.5 %, a WARN above the 0.5 floor and below the
      0.8 target; s3 63.9 % at 4 threads; s4, p2, p3, and the g2
      bjacobi record), and p4 by evidence (no nvcc
      on this machine; `FREHG_P4_EVIDENCE` names the cuda-compile.yml
      requirement on the pushed release SHA — a release-checklist item
      below).
- [x] **§8.2 BC matrix, zero "unverified" rows** (`unit.bc_matrix`):
      117 canonical kind × target × coupling cells, 117 rows, every
      cell `tested` (ids cross-referenced against the gtest/ctest/
      validation universes) or `schema-rejected` (rejection + message
      pinned by a test), zero `limitation` rows. The checker also
      carries the schema-vocabulary drift alarm (a new BC kind cannot
      land without extending the matrix).
- [x] **§8.3 feature table, zero empty cells, GPU rows explicit**
      (`unit.feature_coverage`): 90 rows green; GPU rows are
      `experimental` exactly as §2B.4 defines (compile-verified,
      statically invariant-checked, CPU-physics-verified; on-device
      execution pending the owner p6 bundle) except the two
      CPU-verifiable rows (`p5` discipline, `aijkokkos` rehearsal),
      which are `supported`.

## Also required (§9 Q7 row)

- [x] **Sanitizer matrix over the new labels**
      (`scripts/run_sanitizers.sh --full`, Apple clang ASan+UBSan +
      deps-clang, `build-asan` at `6b0baaa`): clean, 104/104 lanes and 8/8 smokes — the
      unit/mpi/regression labels (now including the Q7 batteries, side
      sweeps, and pair gates) minus two printed lane-composition
      exclusions, plus the shortened nightly-class smokes (smoke-b5 ×4,
      smoke-b6 ×2, and the Q7-added smoke-g5/smoke-g8). The first run
      found the exclusions: 104/108 clean, `perf_baseline` failing its
      wall-clock comparison by construction and the three `aijkokkos`
      lanes aborting on the lane's two Kokkos builds (V2-A20 finding 6;
      not a frehg finding — the release binary links one Kokkos and
      passes those lanes in step 2).
- [x] **Docs**: theory chapters — atmosphere
      (`docs/theory/atmosphere.md`, new: the shared bulk-aerodynamic
      chain), heat (`docs/theory/temperature.md`, Q5), wind
      (`docs/theory/surface-water.md` §wind, Q6); the §8.1
      **symmetry-exemption table** (`docs/theory/symmetry-exemptions.md`,
      new: 9 waived asymmetries with provenance + the measured release
      floors); **parameter-table lockstep** green
      (`unit.parameter_docs` / `check_parameter_docs.py`); testing.md
      rewritten for the v2 label reality and the §8 instruments;
      feature-coverage.md and the Q5/Q7 phase records joined the mkdocs
      nav (strict build green).
- [x] **v2.0.0**: version bumped (CMakeLists `project VERSION`,
      CITATION.cff + date, README status), tagged `v2.0.0` on the close
      commit in the dev repo; synced dev → official per §1.2 (file
      copy; the owner commits and pushes).

## Q7 capability work (the §8 backfill fixes, V2-A19)

- [x] **V2-A11 fix landed** — the two west/south face-area ghost copies
      deleted (`WetDry.cpp::updateGeometry`); gated by
      `regression.outflow_staircase` (Outflow × {W,S,E,N}, both feed
      directions, absolute Manning/volume bands; 10 failures on stock →
      PASS) and the swe dihedral battery (stock mirror failures to
      5.2e-3 → post-fix mixed floor 3.478e-3 = the pure legacy edge
      arithmetic). b1–b6 unchanged: the full anchored tier and the b6
      nightly lanes green in the pipeline. The b5 exposure question
      from Q0.3 is closed by construction and by measurement rather
      than by the b5 nightly lanes (not gates, V2-A20): b5's DEM
      duplicates its boundary rows and columns, so the deletion is a
      self-assignment on that grid, and a stock binary restarted from
      the fixed run's checkpoints reproduces it bitwise on all six
      output fields over two 10 h windows (first wetting; the
      discharge peak). Of the two validation cases on the defective
      path, swe-vcatchment was rerun and regenerated; the
      swere-superslab rerun runs past the tag (carried-forward item 9).
- [x] **Unconditional `Fu`/`Fv` interface exchange** — the second
      battery-caught defect (V2-A19): stock multi-rank transport with
      no eta condition drifted 1.6e-2 (upwind) / 1.6e-3 (superbee)
      across decompositions, eta invariant at 4e-16 throughout;
      `regression.tracer_rank_invariance.{strict,default}` now PASS
      (post-fix invariance at machine precision — report-Q7.md records
      the achieved values); serial behavior and every eta-BC case
      bitwise unchanged.
- [x] **Schema hardening**: `scalar_value` × `groundwater_top` ×
      coupled rejected loudly (coupler-owned top;
      `ConfigTest.RejectsTemperaturePinOnCoupledGroundwaterTop`); the
      §8.2 backfill unit tests (~1300 lines: four-side kind coverage,
      corner/single-cell polygons, masked columns, decomposition,
      coupled budget closures).
- [x] **Gate-first (§6.1)**: every instrument authored red in
      `a812bc4` against stock source (staircase 10 failures, swe
      battery mirror failures, tracer lanes red, coverage checkers red
      by citing the then-unlanded backfill); green in the
      implementation commit with no gate weakened.

## Carried forward (open at v2.0.0, tracked)

1. **Owner p6 GPU bundle** (`scripts/gpu_acceptance.sh`): GPU rows stay
   `experimental` until a passing bundle returns; then the §8.3 rows
   flip to `supported` (a one-line table edit + the artifact under
   `docs/developer-guide/`).
2. **p4 on the release SHA**: cuda-compile.yml must be green on the
   official remote at the pushed v2.0.0 tree before tagging there (the
   workflows are unchanged since the Q0.2 bring-up at `b1ca0d0`; the
   device-facing sources changed in Q5–Q7, so the re-run is required,
   not a formality).
3. **Arm the runner's golden/envelope gates** (V2-A20 finding 2): the
   legacy-goldens archive has never been provisioned on the CI runner,
   so every scheduled nightly has skipped the golden regression
   comparisons and the envelope step (0 s) while reporting success —
   b1–b6 golden/envelope evidence is local-only at this release.
   Either provision the goldens on infrastructure that can hold them
   (the four b5 full-horizon lanes additionally exceed the hosted 6 h
   job cap and need a self-hosted runner or their own split jobs), or
   make the workflow state the local-only scope instead of skipping
   silently.
4. **Owner §6.4 digitization sign-off** (Q4 DoD item): the g5 overlay
   verification remains the one owner data touchpoint, unchanged by Q7.
5. **BSD-3 license review** flagged at P5 remains with the owner.
6. **The thin-film recession audit sink** (V2-A20 finding 4):
   pre-existing, measured and recorded in the swe-vcatchment README;
   investigation spun off as its own task.
7. **Build-time revision stamp** (V2-A20 finding 5): the run record's
   git SHA is captured when CMake configures, with no dirty flag, so an
   incrementally rebuilt binary can name a commit whose code it does
   not contain. Release evidence is unaffected (every evidence binary
   was configured at its tree); the fix is a build-time
   `git describe --always --dirty` header plus a test that a modified
   tree stamps `-dirty`. **Resolved after the tag** (merge `e6f82b9`;
   V2-A20 finding 5, Resolved note).
8. **Sanitize the aijkokkos staging path** (V2-A20 finding 6): this
   machine's sanitizer lane cannot run the `aijkokkos` lanes (its
   clang-built frehg and the gcc-built Kokkos-aware PETSc bring two
   Kokkos builds into one process), so the frehg-side
   `VecGetKokkosView` staging has not run under a sanitizer. It needs
   an ABI-consistent Kokkos-aware PETSc: clang-built in the clang
   prefix, or a Linux gcc ASan lane.
9. **Regenerate the swere-superslab record** (V2-A20 scope note): the
   V2-A11 rerun on the v2.0.0 binary started 2026-09-29 and takes most
   of a day serial, so it finishes after the tag. Then run
   `makeplot.py` (it rewrites both PDFs), replace the README's pre-fix
   table and secondary observations with the post-fix values, and
   commit. Until then the case's numbers and figures are pre-fix and
   must not feed the manuscript. **Resolved after the tag** (V2-A20
   scope note, Resolved note): the rerun finished in 13.1 h,
   `makeplot.py` regenerated both figures, and the README carries the
   post-fix record.

Closed since the Q5/Q6 lists: the **recalibrated-s1-on-runner check**
(V2-A12 residual) — five consecutive scheduled nightlies ran the
recalibrated s1 green on the 4-vCPU runner, 2026-09-22 → 09-26
(V2-A20 finding 3).

## Amendments this phase

V2-A19 (the §8 realization, the V2-A11 landing, the tracer
interface-flux defect, the release-pipeline scope) and V2-A20 (the
pipeline's first execution: b5 full-horizon scoping with the measured
numbers, the V2-A10 vacuous-runner-evidence correction, the V2-A12
residual closure, the thin-film sink, the configure-time revision
stamp, the sanitizer lane composition).
