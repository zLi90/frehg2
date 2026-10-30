# dg — diagnostics and output gates (Stage-3 plan T4)

Per-PR gates for the interval-mean rate outputs (`seepage_rate`, the
`io::IntervalRates` registry) and monitor-row decimation
(`output.monitor_interval`). Harness subcommands `dg1`, `dg2`,
`dg-restart` in `tests/regression/run_regression.py`; pure checks in
`tests/regression/gate_diagnostics.py`; tolerances in
`tolerances/stage3-dg.yaml`; §6.3 battery in `scripts/test_dg_gates.py`.

| Gate | Case | Reference | Asserts |
|---|---|---|---|
| dg1 | `dg1-seepage.yaml` (4×4×8 coupled, SUBCYCLED dt = 1 s / dtg ≤ 10 s, y+ side head 0.3 m above the surface = exfiltration, ponded SW corner over a deep water table = infiltration, reallocation disabled, 1 h, outputs every 300 s) | the defining identity of an interval-mean rate: `Σ_cells rate·A·Δt_out = Δ(mass_audit.seepage)` per interval — both sides sum the SAME per-step gain numbers | identity ≤ 1e-12·max(\|Δ\|, 1e-12 m³) (measured 6.7e-16); **anti-vacuity** (the V3 reproduction): ≥ 1 output frame with snapshot `seepage` ≡ 0 on every wet cell while `seepage_rate` ≠ 0 |
| dg2 | dg1 in three lanes: every-step, `monitor_interval: 60`, and `monitor_interval: 77` + checkpoint 1800 (77 divides neither 1800 nor 3600 — the forced-row paths run by execution) | the every-step run's tables | every decimated row (all columns, every `/monitor` table) is BITWISE the every-step row at the same time; last row at t_end; the 77 s lane's row-time set equals the schedule exactly (t0 + multiples + checkpoint + t_end) |
| restart | dg1 + checkpoint at 1800 s | the continuous run | `seepage_rate` (and seepage/eta/depth) outputs bitwise after restart — the `acc_<name>` fields + `acc_t_prev` scalar round trip |

Every check is provably able to fail (`unit.dg_gates.negative`): a
dropped-substep rate frame, a vacuous V3 lane, a 1e-15 column
perturbation, a missing t_end row, a lost checkpoint-forced row.

**Red evidence (gate-first, plan §2.3):** the T3 binary rejects both
target keys (`output.variables.surface[3]: invalid value 'seepage_rate'`,
`output: unknown key 'monitor_interval'`) and all three gates fail
"capability absent" — `stage3/records/T4/red/`.
