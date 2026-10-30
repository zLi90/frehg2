# DoD — Stage-3 phase T4 (diagnostics and output)

Date: 2026-10-28. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §8; amendment S3-A7 (Field3 leg deferred to its first consumer, T6).

| Check | State | Evidence |
|---|---|---|
| dg1 green: rate-ledger identity measured **6.7e-16** over 12 intervals vs the authored 1e-12 (both sides sum the same per-step gain numbers); the V3 anti-vacuity reproduced (1 output frame with snapshot `seepage` = 0 on every wet cell while `seepage_rate` ≠ 0 — the subcycled dt = 1 / dtg ≤ 10 window) | done | `records/T4/green/dg1-green.log` |
| dg2 green: 61 decimated rows at `monitor_interval: 60` bitwise rows of the every-step tables (mass_audit, gw_mass_audit, point monitor p_corner); the 77 s + checkpoint-1800 lane exercises BOTH forced-row paths (49 rows = t0 + the 46 interval multiples + the checkpoint-forced 1800 + the t_end-forced 3600; schedule asserted exactly) | done | `records/T4/green/dg2-green.log` |
| dg-restart green: `seepage_rate`/`seepage`/eta/depth **bitwise** after restart at 1800 s (the `acc_<name>` fields + `acc_t_prev` scalar round trip; no S3-A6-class relaxation needed — the case carries no superbee/wind branch noise) | done | `records/T4/green/dg-restart-green.log` |
| `default_bitwise` green on the T4 binary (`monitor_interval` absent = a row every step; no rate variables = v2.0 checkpoint layout — bitwise vs the T0 hashes) | done | solo run (15.5 s) + the Tier P stage3 label |
| W4.1 `io/IntervalRates` registry (Field2 leg; S3-A7), per-cell accumulate → interval-mean write → reset; allocates only requested variables | done | `src/io/IntervalRates.{hpp,cpp}` |
| W4.2 `seepage_rate` = interval mean of the coupler's complete per-step gain (applied + bounce-back + vent), m/s positive upward; `seepage` kept and documented as a per-step sample | done | `Coupler::gain()`, `Simulation.cpp` registration/hook |
| W4.3 `output.monitor_interval`: accumulation split from row writing in all three loops; rows forced at checkpoint instants and t_end; restart recomputes the schedule from t0; resolved-config echo + r1 round trip | done | dg2 (incl. `check_run_record` on all lanes); `ConfigTest.MonitorIntervalAndSeepageRateParseRoundTripAndRejectMisuse` |
| Schema: `seepage_rate` in the surface whitelist with the requires-groundwater cross-field rule; `monitor_interval` whole seconds ≥ 0 | done | same ConfigTest; red logs `records/T4/red/` |
| Negative battery `unit.dg_gates.negative` green (identity on an exact synthetic bookkeeping pair; dropped-substep frame, vacuous V3 lane, 1e-15 column perturbation, missing t_end row, lost forced row all fail) | done | ctest |
| Red evidence (T3 binary rejects both target keys; all three gates "capability absent" FAIL) | done | `records/T4/red/` |
| W4.4 docs: `output.md` (snapshot-vs-rate semantics, decimation), `parameters.md` (both keys), `adding-an-output-variable.md` (the interval-rate pattern), `agents/postprocessing.md` (use `seepage_rate`, never snapshot `seepage`); coverage rows; doxygen + mkdocs --strict + all lockstep checkers clean | done | `unit.parameter_docs` (195 tokens), `unit.feature_coverage` (108 rows), `unit.bc_matrix` |
| Tier P green; S3-A1-scoped Tier N green (b5 24 h record + b6 ss/td) | done | `records/T4/timings.md` |
| Snapshot `refs/stage3/T4`; `phase.patch`/`phase.stat` | done | `records/T4/` |
