# rl — reallocation-mode gates (Stage-3 plan T1)

Per-PR gates for `groundwater.reallocation: legacy | disabled` (Stage-3 plan §5;
owner decision O9). Harness subcommands `rl2`, `rl3`, `rl4` in
`tests/regression/run_regression.py`; tolerances in
`tolerances/stage3-rl.yaml`; §6.3 negative battery in
`scripts/test_rl_gates.py`. Each gate stages its case twice and runs a
`legacy` and a `disabled` lane by rewriting `groundwater.reallocation`.

| Gate | Case | Reference | Asserts |
|---|---|---|---|
| rl1 | unit tests (`tests/unit/test_gw.cpp`, `test_config.cpp`) | exact identities | `disabled`: realloc audit exactly 0.0, the step ledger closes to 1e-12, unsaturated heads sit bitwise on the van Genuchten retention curve h(θ); schema round trip + the inert-key warnings |
| rl2 | `rl2-tidal-box.yaml` (8×8×12 sand/silt checkerboard, side water tables 1.5/1.2 m, four 30-min infiltration pulses) | cumulative ledger identity ΔV = Δboundary_in − Δss_storage + Δrealloc − Δvloss (recomputed in-script from `/monitor/gw_mass_audit`) | `disabled`: realloc columns exactly 0, ledger ≤ `ledger_rel`·T (T = Σ\|boundary_in increments\|), \|Δvloss\| ≤ `vloss_rel`·T, θ within the soil-table bounds, substeps(disabled)/substeps(legacy) ≤ `substep_ratio`. Both lanes: Σ per-cell `realloc_created − realloc_dropped` and Σ `vloss` equal the audit columns to `percell_rel`. **Anti-vacuity**: the legacy lane's gross churn (created + dropped) ≥ max(`churn_floor_rel`·T, `churn_vloss_factor`·\|Δvloss_disabled\|) — a case too quiet to exercise reallocation fails its own gate. |
| rl3 | `rl3-steady-column.yaml` (1×1×40 sand, −2e-6 m/s infiltration, base pressure head 0) | the steady state is a uniform Darcy flux −2e-6 m/s | `disabled`: trailing-window θ span ≤ `steady_abs`, max_k \|qz_k − q_top\|/\|q_top\| ≤ `uniform_rel`, realloc exactly 0. `legacy`: churn and uniformity recorded, not gated. |
| rl4 | `x-orient/coupled-salt-base.yaml` + `transport.subsurface_update: consistent` | ledger bound | `disabled`: \|Δsubs_anchor\| ≤ max(`anchor_rel`·M_subs, `anchor_vloss_factor`·s_max·\|Δvloss\|) — with nothing created or dropped, the anchor reduces to the scalar on the clamp volume. `legacy` recorded. |
| rl5 | existing `gw-orient` dihedral battery and `rank-invariance-coupled-salt --mode strict`, run with `FREHG_STAGE3_OVERRIDES=groundwater.reallocation=disabled` | the run set itself | existing criteria; the r1 `--expect` hook proves the override reached every run |

Steadiness is asserted in-gate (rl3's trailing window); every check is
provably able to fail on corrupted synthetic input (`unit.rl_gates.negative`).

**Red evidence (gate-first, plan §2.3):** on the T0 baseline binary the
schema rejects `groundwater.reallocation` and the `realloc_*` output
variables, so rl2/rl3/rl4 fail "capability absent" —
`stage3/records/T1/red/`. Compiled unit gates land with the config change
(amendment S3-A2).
