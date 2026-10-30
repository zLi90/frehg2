# et — evapotranspiration gates (Stage-3 plan T5)

Per-PR gates for the `evapotranspiration` block (open-water + soil
channels, time-basis complementarity, the supply cap, the salt store and
held-seepage salt; plan §9, owner decisions O4/O6/O7). Harness
subcommands `et1`–`et9`, `et-orient`, `et-rank`, `et-restart` in
`tests/regression/run_regression.py`; pure references in
`tests/regression/gate_et.py`; tolerances in `tolerances/stage3-et.yaml`
(adjudications: amendment S3-A8); §6.3 battery in
`scripts/test_et_gates.py`; unit lane `AtmEt.*` in
`tests/unit/test_atm.cpp` (et0).

| Gate | Case | Reference / asserts |
|---|---|---|
| et0 | unit `AtmEt.*` | potentialBulk ≡ the v2.0 evaporationRate tree BITWISE at liu2006/r_s = 0; FAO-56 and log-profile hand values; the Lee–Pielke β closed form; the supply cap both signs; the θ_fc retention default |
| et1 | `et1-column.yaml` (uncoupled sand column, mild met) + `et1b-coupled.yaml` (the coupled D1 lane) | a: new block vs the v2.0 bulk path — same expression tree, bitwise through substep 15, then the fma-contraction class (bounds 1e-7; S3-A8.1). b: the D1 red — the v2.0 guard-skipped predictor leaves a 5872 m first-output head difference; the new path's ledger closes at 1.7e-13 |
| et2 | et1a column as a coupled run with a dry surface | coupled-dry ≡ uncoupled to 1e-10 (measured 8e-11 strict); anti-vacuity: depth ≡ 0 |
| et3 | `et3-pond-switch.yaml` (0.4 mm pond, et-met, dt 1 s / dtg 7 s), lanes bulk / prescribed / v2.0-negative | pond phase = E_w·t and the depth ledger to 1e-12 (measured 3e-15); wet+soil ≤ t; combined rate = E_w·t while α₁ = 1 (measured 0 / 3e-15); post-switch continuity 1.000000; the v2.0 lane's window-lag double count drains the pond early (removal ratio 0.26 vs the 1e-3 violation floor) |
| et4 | `et4-balance.yaml` (closed silt column, 12 h; S3-A8.3) | boundary_in = −soil_evaporation exactly; ΔS = boundary_in − ss + realloc − vloss at 4e-15; vloss ≡ 0 |
| et7 | `et7-cap.yaml` (θ = θ_r + 1e-4, demand 1e-5 m/s) | θ ≥ θ_r − 1e-14; unrealised = demand − applied exactly; vloss = 0 |
| et8 | a `et8a-soil-salt.yaml`; b/c `et8b-pond-store.yaml` (store / discard); d `et8d-held-seepage.yaml` | a: total mass exact, per-cell c·θ to the tiny-K class (6.7e-7 ≤ 1e-5; S3-A8.4), adjust/anchor ≤ 1e-12 M. b: s·h = 35·h₀ wet, store peak = 35·h₀·A at dry-out, release on the rain pulse, total (mass+store) 4e-15. c (real negative): discard loses exactly 35·h₀·A. d: the full audit identity 4e-15 with surf_adjust ≡ 0 EXACTLY (nothing deleted — D9), store > floor while held (S3-A8.5) |
| et9 | `et9-tidal-flat.yaml` + the harness noET lane | sub/supra/intertidal zone assertions; wet+soil ≤ t; both channels > 1e-6 m³; the salt audit identity (2.7e-9 ≤ 1e-7 default-solver class); substep ratio ≤ 1.2; the noET lane BITWISE vs its T5-start pin (1271 datasets — T5 default-neutrality by execution) |
| et-x | `orient-base.yaml` through the dihedral battery and n = 1/2/4 (trajectory-stabilized: 0.1 m bed steps, ≥ 0.05 m shoreline margins; S3-A8.9) | rank: machine precision (θ 1.3e-15, tracer 1.1e-13, store 5e-20). Battery: transpose exact (tracer 0.0, store 8e-20); mixed at the swe 1e-2 class, θ at 0.4 (front displacement under the documented legacy asymmetry) |
| restart | et3-prescribed (checkpoint mid-transition) + et8b (checkpoint with a full store) | all fields bitwise incl. the clocks and the store |

**Red evidence (gate-first, plan §2.3):** the T4 binary rejects every
case on the T5 target keys (`evapotranspiration`, `surface_dryout`,
`output.variables.atmosphere`, `salt_store`) — `stage3/records/T5/red/`.
