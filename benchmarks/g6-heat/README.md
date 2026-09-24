# g6 — subsurface heat, analytic (v2 plan §4.2; per-PR)

Temperature-transport gates against closed forms, all computed in-script
(`tests/regression/gate_heat.py`; §6.3 battery:
`scripts/test_g678_gates.py`). Realization decisions in plan amendment
V2-A17 — notably: every flow-driven sub-gate reads the ACHIEVED Darcy
flux from the output, asserts it within a window of the target, and
evaluates the closed form at the measured value.

| Case | Reference | Criterion |
|---|---|---|
| `g6a-bp.yaml` | Bredehoeft & Papadopulos (1965) steady profile | max-norm ≤ 1% of ΔT at measured Pe ∈ ~{−5,−1,0,+1,+5} (harness rewrites the bottom head per Pe around the hydrostatic baseline) |
| `g6b-ogata.yaml` | Ogata–Banks heat form, OGS parameter set ((ρc)_eff = 2e6, λ = 2.2, v_t = 1.5e-6 m/s) | superbee breakthrough ≤ 1% of the step at ~1/5/10/20/50 m for t ≥ 2e6 s; upwind recorded (1–5%), not gated |
| `g6c-stallman.yaml` | Stallman (1965) damped diel sinusoid, q_z ∈ {0, ±2e-6} m/s | fitted amplitude ratio ≤ 5% and phase lag ≤ 600 s at z = 0.08/0.18/0.28 m below the pinned top cell, relative to the top cell's own fit |
| `g6c-coupled.yaml` | same Stallman roots | the composed coupled variant (plan §4.2): an eta-pinned pond's wet-cell temperature is reset from the diel series and heat enters through the coupled seepage/conduction exchange — same criteria, zero new reference data |
| `orient-base.yaml` | the run's own identity | the §8.1 heat 8-orientation battery (strict-mode, tol 1e-12; measured ~1.4e-13) + the two-scalar restart and rank-invariance lanes |

Bring-up findings recorded in the case comments and V2-A17:

- **Steady horizon (g6a):** Pe ≠ 0 breaks the odd symmetry that kills
  the n = 1 mode at Pe = 0; τ₁ = L²/(α(π² + Pe²/4)) reaches 1.36e7 s at
  |Pe| = 1, so t_end = 8e7 s (a 2e7 s first cut left a 0.83 K transient
  bump — 8× the budget — symmetric in ±Pe).
- **Inlet-cell ambiguity (g6b):** the pinned inlet cell's half-thickness
  is a boundary-position ambiguity that decays with the front width;
  gate window t ≥ 2e6 s plus the dz₀ = 0.125 m mesh keep every point
  under the bound (dz₀ = 0.25 measured 2.08% at the first output; errors
  halved with dz₀ — representation, not scheme).
- **Orientation battery:** caught the legacy y+-only limiter side-ghost
  admission (the b6 sea-side rule) being applied to temperature — a side
  thermal Dirichlet was admitted on y+ but clipped on the other three
  sides (+0.29 K corner artifact via the corner-face spill of side
  polygons). Temperature now admits prescribed side ghosts on all four
  sides; salinity keeps the golden-pinned legacy rule
  (`SubsurfaceTransport.cpp`).

`input/gen_diel.py` and `input/gen_orient_ic.py` regenerate the
committed inputs (§6.2). The harness entry is
`run_regression.py {g6, heat-orient, heat-restart, rank-invariance-heat}`;
tolerances in `tests/regression/tolerances/g6-g8-heat.yaml`.
