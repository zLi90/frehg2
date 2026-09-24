# Q5 Completion Report — Temperature transport (surface + subsurface)

**Phase:** v2 Q5 (plan §4, gates g6/g7/g8 + the heat x-batteries)
**Dates:** 2026-09-24
**Commits:** `b67f9aa` (gates red) → implementation → close (this
report). Amendments V2-A17 (design realization) and V2-A18 (the
orientation-battery catch).

## Summary

Temperature is the model's second registered scalar. The plan §4.1
"two registered scalars, not a general N rewrite" is realized as a
`ScalarSpec` parameterization of the existing `ScalarSolver`: the
salinity spec reproduces the v1 solver byte-for-byte (field names,
checkpoint layout, arithmetic — the full b-tier and both restart lanes
prove it), and the temperature spec switches in the thermal physics:
the θ+κ retardation basis (SEAWAT form), conduction as
λ_eff/(ρc)_w in the dispersion tensor's molecular slot (bulk property,
not θ_s-scaled), donor-value top faces (heat travels with water),
no rain/evaporation dilution, cell-pinning top/bottom Dirichlet rows,
an all-sides limiter ghost admission, and the surface heat-exchange
source (Edinger equilibrium and full bulk modes) audited in a new
`surf_atmos` ledger column. The density coupling gained the thermal
term r_ρ = 1 + β_s·s − β_T·(T − T₀) with every coefficient moved from
compile-time constants to schema-validated configuration whose defaults
are the legacy values (the §4.1 compile-time-constant fix,
golden-neutral by construction).

All gates green at close:

| Gate | Criterion | Measured |
|---|---|---|
| g6a B&P Péclet sweep | ≤ 1% of ΔT at measured Pe ∈ {−5…+5} | 0.087 K worst (Pe ±5) vs 0.1 K; fluxes exact to 1e-19 |
| g6b Ogata–Banks (OGS params, superbee) | ≤ 1% of the step, 5 depths | 0.04–0.38%; upwind recorded 0.3–4.9% |
| g6c Stallman (q_z ∈ {0, ±2e-6}) | amplitude ≤ 5%, phase ≤ 600 s, 3 depths | ≤ 3.4% / ≤ 339 s |
| g6c-coupled (composed, plan §4.2) | same | ≤ 0.9% / ≤ 331 s — the coupled heat path matches the prescribed-BC physics |
| g7a Edinger + energy ledger | 0.5% of excess; ledger ≤ 1e-8/step | 0.0079 K (= the explicit-Euler bias exactly); ledger 4e-16 |
| g7a2 bulk equilibrium | settle at the offline root ± 0.05 °C | |d| < 1e-3 °C at T_e = 20.55 |
| g7b channel plume | u ± 1%; steady L2 ≤ 1%; transient ≤ 2% | 0.5000 exact; 0.01%; 1.75% |
| g8 HRL onset (nightly) | Ra_eff 28.5 conducts, 52.25 convects (Ra_c 39.48) | mode → 2e-6 of seed, Nu = 1.0000; Nu_min = 1.53 |
| heat 8-orientation battery | 1e-12 strict-mode | 1.5e-13 (T), 1.8e-15 (head) |
| two-scalar restart determinism | 1e-12 | exact-class pass |
| heat rank invariance 1/2/4 | strict 1e-12 / default 1e-4 | 4.3e-15 / 4.5e-8 |

## The findings that matter beyond this phase

### 1. The first battery-caught defect (V2-A18)

The heat 8-orientation battery failed its first execution at 1.4–2.9 K
on every y-involving transform. Root cause: the subsurface limiter
admits a prescribed side-ghost value on the **y+ side only** (the
legacy sea-side rule b6 pins), compounded by `BoundarySet`'s
corner-face spill (a side polygon spanning a corner hands the condition
every domain-edge face of the corner cell — load-bearing for b6's
one-cell-wide tank). A side thermal Dirichlet was limiter-admitted on
y+ but clipped on x−/x+/y−. Fix: the temperature instance admits
prescribed side ghosts on all four sides; salinity keeps the
golden-pinned legacy rule; the corner spill is documented, not changed.
After the fix the battery passes at 1.5e-13. This is the first x-gate
catch by *execution* (V2-A9…A12 were vacuous/unexercised-gate
findings) and seeds the §8.1 exemption table with its first two rows.

### 2. Gate-authoring physics caught at authoring time (§6.3 working as designed)

- **The drafted HRL seed was orthogonal to the critical mode** —
  sin(πx/H) has exactly zero projection on the admissible m = 2 cosine
  mode of the 2H box; the growth would have ridden the weaker m = 3
  component near its own (higher) threshold. Caught in review,
  regenerated, and the orthogonality is now itself a negative-battery
  expectation.
- **Gravity in the head sweep**: the drafted g6a head arithmetic
  centered the Pe sweep on dh = 0 — the hydrostatic no-flow baseline is
  h_p,bot = L, not 0. The V2-A17 measured-flux protocol (assert the
  achieved qz within a window, evaluate closed forms at the measured
  value) makes the gates structurally immune to this class; every
  measured flux landed within 1e-17–1e-19 of target.
- **Pinned-cell boundary geometry**: pinning CELLS puts the isothermal
  planes at the pinned centres — g6a gates on L_eff = L − dz, g8 on
  H_eff = H − dz (Ra_eff 52.25/28.5 for the configured 55/30), and g6b
  carries a ±dz₀/2 inlet-plane ambiguity that decays with the front
  width (gating window t ≥ 2e6 s + the dz₀ = 0.125 m mesh).

### 3. Physics/timescale bring-up traps (recorded in the case comments)

- **g6a steady horizon:** Pe ≠ 0 breaks the odd symmetry that
  annihilates the n = 1 transient at Pe = 0; τ₁ = L²/(α(π² + Pe²/4)).
  The first-cut t_end = 2e7 s left a 0.83 K bump symmetric in ±Pe
  (8× budget); 8e7 s ≥ 5.9 τ₁ everywhere.
- **The frictionless g7b channel is a resonator** (discharge inlet and
  eta outlet both reflect; nothing damps): a from-rest start rings
  forever at ±2% velocity. The run initializes at the exact
  uniform-flow steady state; the V2-A16 friction trick is unavailable
  because the gated state has u ≠ 0. The transient criterion convolves
  the closed-form step response with the **measured inlet history**
  (the inlet cell is a finite mixing volume — a 1000 s effective
  release delay that would otherwise eat 20% of the step at the
  steepest point).
- **dtg is the transport clock in gw-only runs**: halving `time.dt`
  without the `groundwater.timestep` clamps left the explicit thermal
  diffusion number at 1.4 → instant blow-up (caught by the gate within
  seconds; the case comments now pair the two).

### 4. Decisions of record

- Temperature is `modules.temperature` + a `temperature:` section (not
  a `transport.scalars` list) — V2-A17 rationale: byte-untouched
  salinity surface, thermal parameter names, and one `scalar:` selector
  on `scalar_value` instead of doubled §8.2 matrix rows.
- `scalar_value` on `groundwater_top/bottom` = cell-pinning Dirichlet,
  **temperature-only in v2.0** (the salinity rows stay schema-rejected
  pending a gate — §8.2 discipline).
- `atmosphere.surface_temperature` is now required exactly when a Q4
  evaporation consumer reads it; the Q5 heat exchange always evaluates
  at the local water temperature. The evaporation VOLUME pathways stay
  exactly Q4's (the latent term shares the formulas, not the plumbing).
- `atmosphere.wind_speed` and `surface_water.wind` stay separate (the
  Q6 open question): scalar bulk-transfer forcing vs vector momentum
  forcing; unifying would impose a vector→scalar convention across two
  independent schema surfaces.
- GLM still-air floor: `wind_speed_floor` default 0.5 m/s applied at
  `MetForcing::sample` for every bulk consumer (no gated case sits
  below it).
- μ(T) and the UNESCO ρ(s,T) polynomial are out of scope (recorded;
  linear ρ(T) is what g8 gates); polynomial deferred to §11.
- Rain enters at the cell's own temperature (measured into the anchor
  column); radiation *schemes* remain out of scope (prescribed
  shortwave/longwave series only), per the plan §10 scope fence.

## Cross-checks and neutrality

- **b6 unchanged with temperature off** and the whole b-tier: the
  refactor is spec-keyed; salinity keeps the "s_*" names, checkpoint
  layout, and arithmetic (guarded branches and exact +0.0 terms).
  Full per-PR tier green post-capability; b6 full-horizon nightly and
  the §4.3 g5 cross-check recorded in dod-Q5.md.
- **g5 cross-check (§4.3):** g5 rerun with the full temperature module
  active at uniform 20 °C (β_T = 2e-4 armed, (T−T₀) ≡ 0) — compared
  BITWISE against the Q4-mode run (stronger than "statistically
  identical"); `scripts/cross_check_g5_temperature.py`.
- **r1 lockstep:** every new schema key round-trips through
  `resolvedConfigYaml` (the emitter mirrors the module cross-rules);
  every gate run passes the r1 `check_run_record` hook.
- **p5:** `check_forbidden.sh` clean; new kernels follow the
  MemSpace-view idiom with no backend branches (one nvcc
  extended-lambda hazard — a lambda-capturing-lambda in the baroclinic
  rework — was flattened before it reached the CUDA lane).

## What Q5 leaves open

- The §8.1/§8.2 x-gate backfill (pre-existing): V2-A18 seeds the
  exemption table (salinity y+ admission; corner-face spill).
- V2-A11 (west/south `outflow`), owner p6 GPU bundle, owner §6.4
  digitization sign-off — unchanged from Q4/Q6.
- Remaining phase: Q7 (release).
