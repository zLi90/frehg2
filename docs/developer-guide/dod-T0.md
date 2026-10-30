# DoD — Stage-3 phase T0 (baseline, reference, harness)

Date: 2026-10-27. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §4.

| Check | State | Evidence |
|---|---|---|
| W0.4 diff `src/` vs `2026Frehg2/frehg2-2.0.0/src` empty (apart from `.DS_Store`) | done | `stage3/records/T0/baseline.md` |
| Four stagefix tolerance amendments reverted (b6 head floor 0.015, salt 0.10; b5 bjacobi 2e-3/3e-5) | done | `stage3/records/T0/baseline.md` |
| Tier P green, timings recorded | done (60/61 at `-j8`; `perf_baseline` failed under 8-way contention and **passed solo**, 123.6 s — a scheduling artifact, not a regression; see report §3) | `stage3/records/T0/tierP-regression.log`, `green/perf_baseline_solo.log` |
| Tier N green, timings recorded | done | `stage3/records/T0/tierN-nightly.log`, `tierN-b5record.log` |
| `stage3-T0.json` reference hashes recorded (13 cases, 1324 datasets) and `regression.s3.default_bitwise` green | done | `tests/regression/reference-hashes/stage3-T0.json`; report §4 |
| Override mechanism unit-tested (`--expect` fails on a missing/mismatched key) | done | `tools/test_run_record_checker.py` (4 new cases), `unit.run_record_checker` |
| W0.9 warning channel unit-tested and accepted by the run-record checker | done | `ConfigTest.ValidationWarningsDoNotInvalidate`; checker `warnings` case |
| `check_forbidden.sh` clean | done | report §5 |
| Docs gate (doxygen 1.17.0 `WARN_AS_ERROR`, `mkdocs build --strict`) clean | done | report §5 |
| Snapshot `refs/stage3/T0` exists | done | report §6 |
