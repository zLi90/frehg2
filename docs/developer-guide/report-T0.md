# Report — Stage-3 phase T0 (baseline, reference, harness)

Date: 2026-10-27. Companion: `dod-T0.md`. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §4.

## 1. Baseline reconciliation

- Snapshot `refs/stage3/T0-start` = `6b9f1562eaae` captured the pre-phase working tree
  (the owner's uncommitted 2026-10-05 work **plus** the 2026-10-26 stagefix prototype).
- The stagefix (`stage3/seeds/stagefix_v2.patch`, 14 files) was reverse-applied cleanly.
  This removed the hard-wired at-rest closure (Fix A), the stage-deposit removal (Fix B),
  the edge-diffusion audit booking, two unit tests, the doc hunks, **and the four
  tolerance amendments** (b6 `head_abs_floor` back to 0.015 m, `salt_mass_rel` back to
  0.10, b5 bjacobi default bounds back to 2.0e-3 / 3.0e-5). Those capabilities return in
  T2 in selectable form behind `inflow: {closure: ...}`, per owner decision O1.
- After the revert, `src/` is byte-identical to the study's production source
  `2026Frehg2/frehg2-2.0.0/src` (`.DS_Store` aside) — `stage3/records/T0/baseline.md`.
  The remaining 29 modified files against `main` (99fea72) are the owner's own
  2026-10-05 work, untouched.

## 2. What T0 added (no physics change)

| Piece | Files |
|---|---|
| Records root, snapshot tool | `stage3/` (outside the git tree), `stage3/tools/snapshot.sh` |
| Frozen reference binary | `stage3/ref-bin/frehg-T0` (+ `deps.txt`), stamp `99fea72cfe26-dirty` |
| Per-dataset bitwise hasher | `tools/h5_bitwise.py` (hash / compare, `--allow` globs for §2.6 declared exceptions, exit 77 off the recording toolchain) |
| Default-bitwise gate | `run_regression.py` subcommand `default-bitwise` (13 cases, `BITWISE_PINS` for flippable defaults), ctest `regression.s3.default_bitwise` (label `stage3_local`, `SKIP_RETURN_CODE 77`) |
| Override channel | `FREHG_STAGE3_OVERRIDES` applied at staging, module-gated; every applied override asserted in the run record via `check_run_record.py --expect` (anti-vacuity, plan §3.2) |
| Shared-harness battery | `scripts/test_s3_shared.py` → `unit.s3_shared` (hasher sensitivity incl. a 1e-15 perturbation failing, allow-glob semantics, toolchain skip, override gating) |
| W0.9 warning channel | `ValidationResult::warnings` + `Context::addWarning` (ConfigSchema), carried into `FrehgConfig::validationWarnings`, printed by `--validate` as `WARNING:` (still `VALID:`), logged at run start, persisted as the run record's `warnings:` list (`RunRecord::addWarning`); documented in `output.md` |
| Tolerance skeletons | `tests/regression/tolerances/stage3-{rl,ob,pd,dg,et}.yaml` (empty; filled at each phase's gate authoring) |
| Plan in-repo | `docs/developer-guide/FREHG2_STAGE3_DEVELOPMENT_PLAN.md` + mkdocs nav |

Unit-test additions: `ConfigTest.ValidationWarningsDoNotInvalidate`; run-record checker
battery gains 4 `--expect` cases and 1 `warnings`-list case.

## 3. Tier P (per-PR) measurements

`ctest -L '^regression$' -j8`, M3, `OMP_NUM_THREADS=1`:

- **60/61 passed, wall 1518 s ≈ 25.3 min** (`tierP-regression.log`).
- Unit + mpi tiers separately: 66/66, wall 356 s (`unit` 14.2 s·proc over 44 tests,
  `mpi` 2409 s·proc over 22).
- The one failure, `regression.perf_baseline`, showed **uniform** +53…+79 % across every
  timer section — the signature of core contention (the serial 404×220×25 perf case ran
  concurrently with 7 other tests on 4 P-cores), not of a code regression. **Solo rerun
  on the idle machine: PASS, 123.6 s** (`green/perf_baseline_solo.log`). Harness note
  for every later tier run: execute the `-j8` sweep with `-E perf_baseline`, then run
  `perf_baseline` serially — no tolerance was changed and no amendment is needed (the
  gate is unchanged; its documented measurement condition is an idle machine, the A23
  min-over-repeats protocol).

## 4. Reference hashes and the default-bitwise gate

- `tests/regression/reference-hashes/stage3-T0.json`: 13 cases, 1324 datasets, 131 kB,
  with the toolchain fingerprint (g++-16 16.2.0, arm64 Darwin).
- Compare mode (the actual gate): **PASS in 16.7 s wall** — all 13 cases rerun and
  bitwise-identical. Well inside the ≈ 40 s budget.
- Sensitivity proven in `unit.s3_shared`: a 1e-15 perturbation of one dataset fails the
  gate; an `--allow` glob exempts exactly the named datasets; a foreign toolchain skips
  (77) rather than fake-passing.

## 5. Gates and hygiene

- `scripts/check_forbidden.sh`: clean.
- Docs: `doxygen docs/Doxyfile` (1.17.0, warnings-as-errors): clean;
  `python3 -m mkdocs build --strict`: clean.
- Zero-warning build (`-Werror`, g++-16): clean (`ld` duplicate-library notes are
  pre-existing and come from the linker, not the compiler).

## 6. Tier N (phase close) measurements

Nightly batch (`ctest -L regression_nightly -E '^regression\.b5\.(rain|norain)\.'`,
`-j8`): **10/10 PASS, wall 11,510 s ≈ 3 h 12 min** (`tierN-nightly.log`). Per test
(contended — 61,855 s·proc over 11,510 s wall = 5.4× oversubscription of 4 P-cores):
g8 67 s; b6.ss/td 506/1475 s; g1.b6.ss/td.amg 517/1515 s; g1.b5_rank amg 394 s;
g5 + transposed 2846/2921 s; p1.b6 aijkokkos 2009 s; g1.b5.rain.sync.amg 10,599 s.

b5 rain/sync 24 h record, bjacobi, 4 ranks, idle machine (release step 3): **PASS,
8186 s wall** (`tierN-b5record.log`). Numbers equal the pre-stagefix record:
discharge 0/48 outside, peak +2.6 %, integral +2.4 %; ponding 0/54, −1.5 %, −2.5 %;
surface residual −10.0 m³ (pre-existing, bound 11.1), subsurface 3.9e-8 m³, clamp
creation 180 m³. The first attempt failed on a relative `--work` path (mpiexec ranks
resolve it from their own cwd) — harness note: always pass absolute `--work` when
invoking `run_regression.py` outside ctest.

Scheduling lesson for T1–T7 closes: run the nightly batch at `-j2` (long tests
first) and `perf_baseline` solo; a phase-scoped Tier N is proposed to the owner as
S3-A1 (undecided at T0 close).

## 7. Deviations

None. No amendment was required in T0.
