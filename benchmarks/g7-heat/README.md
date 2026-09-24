# g7 — surface heat exchange, analytic (v2 plan §4.2; per-PR)

Surface-temperature gates against closed forms computed in-script
(`tests/regression/gate_heat.py`; §6.3 battery:
`scripts/test_g678_gates.py`).

| Case | Reference | Criterion |
|---|---|---|
| `g7a-edinger.yaml` | Edinger et al. (1968) equilibrium relaxation (still basin, T₀ = 30 °C, T_e = 20 °C, K_e = 30 W/m²/K, e-fold 1.614 d) | T(t) within 0.5% of the initial excess (measured 0.0079 K — the explicit-Euler bias exactly); per-step heat-ledger closure ≤ 1e-8 of the basin content (measured ~4e-16) |
| `g7a2-bulk.yaml` | offline root of the identical bulk chain (Tetens/Liu/L_v(T)/εσT⁴; hand-checked spot values Q(20) = +44.9, Q(22) = −122.2 W/m², root 20.55 °C) | the basin settles at the root within 0.05 °C (measured |d| < 1e-3) plus a trailing-steadiness bound |
| `g7b-channel.yaml` | van Genuchten (2013) river ADE + first-order exchange (u = 0.5 m/s, D_L = 5 m²/s, K_e = 25 W/m²/K, T_e = ambient) | measured u within 1% of nominal; steady profile L2 ≤ 1% of the inlet excess (0.01%); mid-channel transient ≤ 2% of the step vs the closed-form step response **convolved with the measured inlet history** (1.75%) — the inlet cell is a finite mixing volume, so its measured rise is the boundary signal (the measured-flux principle applied to the inlet) |

Bring-up findings (case comments):

- **The frictionless channel is a resonator** — a discharge inlet and an
  eta-Dirichlet outlet both reflect, so a from-rest start rings forever
  (±2% velocity standing wave after 7 transit times). The run
  initializes AT the uniform-flow steady state (`uu: 0.5`), which is
  exact; the V2-A16 friction trick is unavailable because the gated
  state has u ≠ 0.
- **Front resolution (§6.3 rule 3):** the breakthrough front width
  √(2Dt) ≈ 917 m at mid-channel; dx = 500 m (< 2 cells) measured 6.2%
  of the step, dx = 100 m (~9 cells) measures 1.75% — gate runs dx=100.

Harness: `run_regression.py g7`; tolerances in
`tests/regression/tolerances/g6-g8-heat.yaml`.
