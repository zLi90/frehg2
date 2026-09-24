# g8 — Horton–Rogers–Lapwood convection onset (v2 plan §4.2; nightly, recommended)

The one deterministic pass/fail gate available for density-coupled flow:
a saturated porous slab heated from below brackets the analytic critical
Rayleigh number Ra_c = 4π² ≈ 39.48 — the gate that exercises the
thermal density term in the Darcy momentum, which g6/g7 never touch.
(The thermal Elder problem is explicitly rejected: ≥ 11 grid-dependent
steady states, plan §4.2.)

Ra in hydraulic-conductivity form: Ra = K β_T ΔT H_eff / α_e with
α_e = λ_b/(ρc)_w = 5e-7 m²/s. Pinning CELLS puts the isothermal planes
at the pinned-cell centres, so H_eff = H − dz = 0.95 m (V2-A17):

| Variant | K [m/s] | Ra_eff | Criterion | Measured |
|---|---|---|---|---|
| subcritical (harness rewrites ks*) | 7.5e-3 | 28.5 | seeded-mode amplitude decays below 10% of seed | 2e-6; trailing Nu = 1.0000 |
| supercritical (`g8-hrl.yaml`) | 1.375e-2 | 52.25 | trailing-half Nusselt > 1.05 persistently | Nu_min = 1.53; mode saturates at 2.8 K |

The 2H-wide box with insulated impermeable sides admits the critical
wavenumber a_c = π/H exactly (its m = 2 cosine mode); the committed
seed `input/ic_temperature.dat` (generator `input/gen_hrl_ic.py`) is
the conduction profile + 0.5 K of that mode. An earlier draft seeded
sin(πx/H), which is ORTHOGONAL to the admissible cosine mode (its
projection is exactly zero — the growth would have ridden the weaker
m = 3 component); caught by the §6.3 review, and the negative battery
now pins the orthogonality.

The Nusselt number is evaluated at the mid-plane from the temperature
and qz snapshots (advective + conductive vertical heat flux over the
conductive reference); the mode amplitude is the projection of the
horizontal-mean-removed temperature onto the seeded mode. Both
extractors carry §6.3 spot checks (conduction ⇒ Nu = 1 exactly;
synthetic roll ⇒ Nu > 1; seeded 0.37 K recovered to 1e-9).

Harness: `run_regression.py g8`; tolerances in
`tests/regression/tolerances/g6-g8-heat.yaml`.
