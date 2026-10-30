# DoD — Stage-3 phase T8 (default-flip evaluations, release pipeline, records)

Date: 2026-10-30. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §12; amendment S3-A11. Companion: `report-T8.md`.

| Check | State | Evidence |
|---|---|---|
| W8.1 reallocation DF: `build-s3-df` with `Reallocation::Off`, every test, unchanged tolerances — **14 fail** (b2/b3/g1 lanes, gw_orient, b5 rank/restart, b6 smoke timeout + subcycled + ss + td) ⇒ **flip rejected, `legacy` stays**, report filed | done | `stage3/df/realloc-report.md`, `w81-*.log` |
| W8.2 salt-store DF: same procedure with `SurfaceDryout::Store` — **all green** (60/60 regression; b6 ss 3.51 %/+5.47 %, td 4.82 %/+8.79 % inside their bands; td interface MAE 0.0389 m better than the 0.0430 m golden) ⇒ **flip KEPT: `surface_dryout` defaults to `store`**; combined run not owed (§2.7-5: only when both pass) | done | `stage3/df/w82-*.log`, `w82-b6{ss,td}-margins.log` |
| W8.3 pin update: `BITWISE_PINS` populated (`reallocation: legacy`, `surface_dryout: discard` — the §2.6 table was empty, a flip-survival defect); `et9-noet.yaml` pins discard in-file; 4 escaped unit tests pinned `Reallocation::Legacy`; `Config.cpp` absent-key paths now keep the Config.hpp initializer; `parameters.md` default column + theory + recipes updated; `ConfigTest.SurfaceDryoutDefaultIsStoreAfterTheT8Flip`; the ET+salinity rule fires only on explicit `discard` | done | S3-A11; 138/138 `records/T8/final-unit-mpi-s3.log` |
| W8.4 p1 lanes: `regression.s3.p1.{ob3,pd2coupled,et9}.omp4.{aij,aijkokkos}` (stage3_nightly; `--pd2-lane` added) — 6/6 green; et9's noET bitwise pin skipped under `--mat-type` by adjudication (recorded at 1 thread/aij; the tolerance criteria gate) | done | `w89-release-gate.log` + `w89-et9-p1-retest.log` |
| W8.5 forbidden scan + p5 negative compile | done | clean / green in the unit tier |
| W8.6 sanitizers with the stage3 labels: full matrix **184 tests, zero findings** (the one ctest failure was the et9 pin under instrumented codegen — `FREHG_SANITIZED_BUILD` guard added, retest green) | done | `w86-sanitizers.log`, `w86-et9-retest.log` |
| W8.7 coverage lockstep: +2 periodic-axis rejection cells (pd6 now probes BOTH axes), `unit.bc_matrix` green (122 rows = 122 cells), `unit.feature_coverage` green (127 rows), **zero empty cells** | done | checker outputs |
| W8.8 docs gate: doxygen 0 warnings, `mkdocs build --strict` green, parameter lockstep 235 tokens | done | release-gate step 1 + local |
| W8.9 `ci_release_gate.sh build-s3` (extended with the stage3/stage3_nightly tiers): steps 1–4 + b5 record green; resumed stage3_nightly + scaling_nightly green after the et9-pin adjudication | done | `w89-release-gate.log`, `w89-resume.log`, walls in `timings.md` |
| W8.10 study migration note (tide closure, periodic AX runs + generator field, ET recipe, DF outcomes, observability; no study file edited) | done | `stage3/records/T8/study-migration.md` |
| W8.11 completion report | done | `report-T8.md` |
| W8.12 snapshot `refs/stage3/T8`; full Stage-3 diff vs `refs/stage3/T0-start` | done | `records/T8/stage3-full.patch` |
| DF outcomes recorded (one kept, one reported) | done | S3-A11 + `stage3/df/` |
