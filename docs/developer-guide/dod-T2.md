# DoD — Stage-3 phase T2 (selectable inflow closure at prescribed-stage boundaries)

Date: 2026-10-28. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §6; amendment S3-A5 (gate adjudications at first green).

| Check | State | Evidence |
|---|---|---|
| ob1 green: discrete-chain entrance head, 3 lanes (external 1.8e-6 rel vs the chain, ratio 0.9996/0.9023 chain-exact; zero_gradient real-negative u_C = 2.2×u_A at R_C 0.002) | done | `records/T2/green/` (stage3 tier log) |
| ob2 unit green (Bernoulli entrance 0.90–1.15×v²/2g external; (v²−u_ext²)/2g with a moving reservoir; zero_gradient regime pinned in-test) | done | `SweModule.ExternalStageInflowPaysTheBernoulliEntranceHead`, `SweModule.ExternalReservoirVelocityEntersTheEntranceForm` |
| ob3 green (strict-solver ledgers 1.2e-8/2.5e-8 m³ over X = 3.0e4; per-cycle prism agreement 3.1e-5/7.6e-6 vs the 1e-3 floor, Fr² = 9.7e-5 measured in-gate) | done | same |
| ob4 unit green (external: stage holds exactly, deposit booked, closure 1e-10; zero_gradient: booked without moving, bitwise dynamics) | done | `CoupledModule.ExternalStageAbsorbsAndBooksCoupledDeposits`, `CoupledModule.ZeroGradientStageBooksDepositsWithoutMovingThem` |
| ob5 green (external KE monotone decay, influx 74 J = 6.2e-3·KE0; zero_gradient contrast 90× — the production-gyre mechanism reproduced and closed in 20×30×1200 s) | done | same |
| ob6 `default_bitwise` green with the declared §2.6 exceptions (`--allow monitor/{mass,transport,temperature}_audit`) | done | stage3 tier log |
| ob7 unit green (schema round trip; inflow on non-eta rejected; velocity under zero_gradient rejected) | done | `ConfigTest.InflowClosureParsesRoundTripsAndRejectsMisuse` |
| ob8 nightly record: ss and td PASS on the experiment criteria; golden residuals recorded (td: one beach cell 0.0196 m vs the 0.015 m zero_gradient floor; salt mass 9.88 %; interface MAE **improves**, 0.0381 m vs the golden's 0.0430 m) | done | `records/T2/green/ob8-rerun.log` |
| Negative battery `unit.ob_gates.negative` green (chain two-way reference check incl. the continuum Bernoulli+Manning closed form at 0.5 %; factor-two weight-1 entrance fails) | done | ctest |
| Red evidence (schema rejects `inflow` on the T1 binary) | done | `records/T2/red/` |
| Entrance-head diagnostics in the run record (`stage_ke_influx_J`, `stage_entrance_head_ratio`) + run-end warning | done | W2.6; ob5/ob1 records |
| W2.7 docs (theory section, parameters, output, case-setup recipe, coverage + bc-matrix rows incl. the extended canonical cells) + lockstep checkers + forbidden scan + mkdocs/doxygen | done | `unit.parameter_docs`, `unit.bc_matrix` (120 cells), `unit.feature_coverage` green |
| Tier P green; S3-A1-scoped Tier N green | done | `records/T2/tierP.log`, `tierN.log`, `timings.md` |
| Snapshot `refs/stage3/T2`; `phase.patch`/`phase.stat` | done | `records/T2/` |
