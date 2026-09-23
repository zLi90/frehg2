# g4 — analytic drawdown + evaporative concentration (v2 Q4, plan §3.3)

Per-PR gate, four sub-criteria on tiny closed basins (harness:
`tests/regression/run_regression.py g4`; criteria in
`tests/regression/gate_evap.py`, tolerances in `tolerances/g4-evap.yaml`,
§6.3 negative battery in `scripts/test_g45_gates.py`).

| Sub-case | Reference | State against pre-Q4 code |
|---|---|---|
| (a) `g4a-drawdown` | η(t) = η₀ − E·t (prescribed E) | **passes** (drawdown exact to ~1e-14 m) |
| (b) `g4b-concentration` | s(t) = s₀V₀/V(t), salt mass constant | **fails**: concentration exact, but prescribed evaporation induces a scalar leak through the domain-edge faces (measured 6.6e-3 of the salt mass over the run; zero when E = 0) — the stale edge ghosts make the solve produce pre-correction edge flow rates that the velocity walls cancel for water but not for the transport step's flux snapshot |
| (c) `g4c-bulk` | same, E from the bulk-aerodynamic closed form | **fails**: capability absent (the `atmosphere:`/`evaporation: bulk` keys are the Q4 target schema) |
| (d) `g4d-drain` | depth ≥ 0 through dry-out + audited shortfall | **fails**: the pre-Q4 ledger books the potential rate after the basin dries (closure error = 100 % of the initial volume) |

Authored failing per §6.1 (the Q2 "gates (red)" precedent); the Q4
capability PRs turn (b)-(d) green.

Two bring-up findings recorded here because no golden pins them:

1. **A "closed" basin needs explicit walls.** The default domain edge keeps
   the legacy transmissive behavior (b1's drain). Under evaporation the
   interior drops below the outside stage and water flows back *in* —
   unbounded volume growth (measured: η grows ×1.2/step at dt = 10 s).
   Every g4 case therefore closes its edges with a `kind: velocity`,
   value 0 condition.
2. **Prescribed evaporation was authored-unexercised in v1**: b1 sets
   `evaporation: 0`, and no other gate enables it, so neither finding was
   visible to the golden suite (the V2-A11 pattern on the evaporation axis).
