# ob — open-boundary inflow-closure gates (Stage-3 plan T2)

Per-PR gates for the selectable inflow closure at prescribed-stage
boundaries, `boundary_conditions[].inflow: {closure: zero_gradient |
external, velocity: {u, v}}` (Stage-3 plan §6/§13.A; owner decision O1/O2).
Harness subcommands `ob1`, `ob3`, `ob5` in
`tests/regression/run_regression.py`; pure references and checks in
`tests/regression/gate_openbc.py`; tolerances in
`tolerances/stage3-ob.yaml`; §6.3 negative battery in
`scripts/test_ob_gates.py`; nightly record `ob8` (b6 ss/td under
`external`).

| Gate | Case | Reference | Asserts |
|---|---|---|---|
| ob1 | `ob1-channel.yaml`, lanes A `external`, B `external` + u_ext = 0.1, C `zero_gradient` | the **discrete steady face chain** recomputed in-script from the code's own formulas (Momentum.cpp advection incl. the legacy double drag factor, WetDry.cpp max-of-two face geometry, the W2.3 conservative entrance), shot on Q to hit the east stage. Authoring self-check: per-face residuals 2.8e-16; lane A u_in 0.3052 m/s vs the independent continuum Bernoulli+Manning closed form 0.3067 (0.5 %); entrance ratio R_A = 0.9996 (the 2026-10-26 prototype measured 1.01). | A and B: face-velocity profile vs the chain; entrance-head ratio vs the chain's; trailing-window steadiness. C (the REAL negative lane): u_C ≥ 2·u_A (chain: 4.4×) and R_C ≤ 0.2 (chain: 0.0019). |
| ob2/ob4/ob7 | compiled unit gates (`test_swe.cpp`, `test_coupling.cpp`, `test_config.cpp`) | hand formulas / identities | stencil state per closure; deposit routing + `stage_deposit` booking; schema messages (land with the capability, amendment S3-A2) |
| ob3 | `ob3-tidal-basin.yaml`, both closures | ledger identity + the P14 momentum-order argument | per lane the surface ledger closes; between lanes the per-cycle tidal prism agrees to max(coef·Fr², floor) with Fr² measured in-gate (anti-vacuity: the bound's own scale is printed) |
| ob5 | `ob5-gyre.yaml`, both closures | energy accounting (the production-gyre mechanism in miniature, report_V1) | external: domain KE decays monotonically after the seed adjustment and the stage-face KE influx stays ≤ influx_rel·KE0. zero_gradient: imports ≥ 10× the external lane's KE (the real negative + anti-vacuity floor). Contingency C-OB5 (plan §16) holds the fallback. |
| ob6 | `regression.s3.default_bitwise` | T0 hashes | no `inflow` key ⇒ bitwise v2.0 (declared exceptions: the appended `mass_audit.stage_deposit` column and the `surf_boundary` audit fix, plan §2.6) |
| ob8 | b6 ss/td with `closure: external` on every eta BC (`stage3_nightly`) | the Kuan experiment | interface MAE within the existing experiment bounds; legacy-golden residuals recorded, not gated (the golden embeds the stage-deposit kick) |

Steadiness is asserted in-gate (ob1's trailing window; ob5's monotone
tail); every check is provably able to fail on corrupted input, including
the factor-two weight-1 entrance (`unit.ob_gates.negative`).

**Red evidence (gate-first, plan §2.3):** on the T1 baseline binary the
schema rejects the `inflow` block, so ob1/ob3/ob5 fail "capability
absent" — `stage3/records/T2/red/`.
