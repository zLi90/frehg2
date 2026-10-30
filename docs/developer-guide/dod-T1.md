# DoD — Stage-3 phase T1 (`groundwater.reallocation: legacy | disabled`)

Date: 2026-10-27. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §5; amendments S3-A2, S3-A3, S3-A4.

| Check | State | Evidence |
|---|---|---|
| rl1 unit gates green (exact-zero audit + ledger under forced saturation; bitwise retention-curve heads in a vadose column; schema round trip, `off`-token rejection, inert-key warnings) | done | `GwModule.ReallocationOffConservesWaterExactly`, `GwModule.ReallocationOffKeepsIsolatedHeadsOnTheRetentionCurve`, `ConfigTest.ReallocationModeParsesRoundTripsAndRejectsLegacyOnlyKeys`; `unit.all` green |
| rl2 green (disabled: ledger 3.3e-13, realloc exactly 0.0, θ bounds; both lanes: per-cell sums ≡ audit columns; legacy churn anti-vacuity 1.448 m³ ≥ 0.1·T; disabled vloss/substeps recorded per S3-A4) | done | `records/T1/green/rl-suite-verbose.log` |
| rl3 green (vadose: steadiness exactly 0, flux uniformity 9.3e-13, vloss exactly 0, substeps 1.000, realloc exactly 0) | done | same |
| rl4 green (anchor: disabled 1.366 vs legacy 38.70 over M = 5.80e4 — 28×; three-part bound) | done | same |
| rl5 green (vadose 8-orientation battery strict 1e-12 under `disabled`; coupled-salt strict rank invariance under `disabled`, r1 `--expect` proving the override) | done | same |
| rl-restart green (per-cell diagnostics checkpoint/restart bitwise, all five fields) | done | `records/T1/green/` (verbose ctest) |
| Negative battery `unit.rl_gates.negative` green (ledger two-way reference check; every criterion fails on its target corruption) | done | ctest |
| Red evidence recorded before the capability (rl2/rl3/rl4 "capability absent" on the T0 binary) | done | `records/T1/red/` |
| `default_bitwise` green (legacy default is bitwise v2.0) | done | tier log |
| W1.5 docs (theory post-allocation-modes section with the S3-A4 domain statement; parameters/output rows; coverage rows) + lockstep checkers + `check_forbidden` + doxygen/mkdocs strict | done | `unit.parameter_docs`, `unit.feature_coverage`, `unit.bc_matrix` green |
| W1.6 legacy pins listed | done | report §4 |
| W1.8 DF preview filed | done | `records/T1/df-preview.md` |
| W1.9 H2 evidence run filed | done | `stage3/evidence/T1-h2-2d/`, report §6 |
| Tier P green; S3-A1-scoped Tier N green; timings recorded | done | `records/T1/tierP-*.log`, `tierN-*.log`, `timings.md` |
| Snapshot `refs/stage3/T1`; `phase.patch`/`phase.stat` | done | `records/T1/` |
