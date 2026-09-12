# Definition of Done — Q2 (Run provenance and timing log)

Per v2 development plan §2A and §9, under the v1 plan's §11.4 DoD
conventions. Every item names the command that verifies it; commands were
run green on 2026-09-10/11 on macOS arm64 (gcc-15, Apple clang 15, MPICH
4.3, PETSc 3.25.1 + hypre 3.1.0, Kokkos 5.1.1; regression runs pin
`FI_PROVIDER=tcp`).

## Deliverables (v2 plan §2A.2) — all implemented and tested

- [x] **`run-record.yaml`** written by `io::RunRecord` into the HDF5
      output's directory: provenance (version/git SHA/build type/hostname,
      MPI ranks + [px, py] decomposition, OpenMP threads, Kokkos backend,
      ISO-8601 start/end, wall seconds, input path + SHA-256, restart
      parentage, `finished`), the **resolved** configuration, modules,
      per-BC summary (kind/target/value with series-file SHA-256, polygon
      vertex count + bounding box, global member cells), the merged timer
      tree, per-system solver telemetry, and the closing mass-audit
      budgets — every field asserted by `tools/check_run_record.py` on
      every gate run (r1).
- [x] **Rewritten at every output flush, finalized after run()** — six
      `flushRunRecord(false)` sites (one per output flush in the three
      loops and their exits) + `finalizeRunRecord()` in main after the
      timer report, write-then-rename so a crash cannot truncate the
      previous record. Verified: a completed run carries
      `finished: true` + the full `simulation` timer; the checker rejects
      unfinished records unless `--allow-partial`.
- [x] **Resolved-config serializer** (`resolvedConfigYaml`, the exact
      inverse of extraction; `frehg --resolve` prints it fenced by
      sentinels) — round-trip green on every gate case: the embedded
      configuration revalidates AND equals a re-resolve of the original
      input as parsed YAML (r1 check 4).
- [x] **New timer sections** completing the owner-requested breakdown:
      `init` (pre-loop restore/first outputs), `io/output`,
      `io/checkpoint`, `io/run_record`, `halo` (every exchange funnels
      through the one timed `exchangeSelected` path), `monitors`,
      `swe/begin_step` — plus the pre-existing swe/gw/transport/solve
      sections. Structured access via `Timer::merged()` shared with the
      stdout report (they cannot disagree).
- [x] **Solver + closure sections** feed from the same reductions as the
      stdout `solver summary` lines and the mass-audit monitors (one
      source of truth each).
- [x] **`run_scaling.py` reads the record** (stdout scrape kept as
      fallback only) — v2 plan §2A.4.

## Gates (v2 plan §2A.3, §9 Q2 row)

- [x] **r1 (per-PR, rides every lane):** `run_case` in the regression
      harness checks the record after every successful gate run —
      presence, schema, launch geometry (ranks/threads/decomposition),
      finalization, configuration round-trip. Green across the full
      per-PR sweep: **40/40** (`ctest -L '^unit$|^regression$' -j 1`,
      2026-09-11), which spans surface-only (b1/b4), gw-only (b2/b3),
      coupled (b5 lanes), coupled+transport (b6 lanes), restart
      (b1/b2/b5/b6 restart lanes; parentage verified in the restarted
      record), rank counts 1/2/4, and both amg/gamg solver lanes.
- [x] **r1 negative battery (§6.3):** `unit.run_record_checker` — a
      synthetic valid record passes; 15 corruptions (each missing
      key/field, launch mismatches, unfinished) each fail.
- [x] **r2 (per-PR):** `regression.r2_timer_coverage` — top-level
      sections tile the time loop at 1 rank: **100 %** measured
      (criterion ≥ 90 %); `halo` (0.139 s) and `io` nonzero at 4 ranks;
      run-record writing **0.01–0.03 %** of the loop (criterion ≤ 1 %).
- [x] **Physics untouched:** the bitwise restart lanes
      (b1/b2/b5/b6-restart) and both rank-invariance modes ran green in
      the same sweep — the new timers and the record writer perturb
      nothing observable.

## Standing criteria (v1 plan §11.2, unchanged)

- [x] Zero-warning build (gcc-15 `-Werror` full flag set) —
      `cmake --build build`.
- [x] Unit label green (12 entries incl. the new checker battery) —
      `ctest -L unit`.
- [x] Forbidden-pattern scan clean — `scripts/check_forbidden.sh`
      (the record's SHA-256 helper was rewritten iostream-only after the
      scan rejected `snprintf`).
- [x] Parameter-docs lockstep — `python3 scripts/check_parameter_docs.py`
      (no schema change in Q2; the run record is output, not input).
- [x] mkdocs strict build with the new output-reference chapter —
      `python3 -m mkdocs build --strict`.
- [x] Sanitizer lane (Apple clang ASan+UBSan) green over unit + mpi with
      the record writer exercised — `scripts/run_sanitizers.sh`
      (see report-Q2.md for the record).
