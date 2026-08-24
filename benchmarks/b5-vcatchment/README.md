# b5-vcatchment — Tilted-V catchment (Kollet et al. 2017)

Coupled surface–subsurface intercomparison case (plan §9 b5; the blocking
P3 gate). The model runs the symmetric half of the published tilted V —
one 100 m hillslope (5 % toward the channel, 2 % along it) plus half the
channel (rows j = 0..4) — exactly as the SERGHEI reference run did; gate
comparisons double the model discharge and storage
(`plot_discharge.py: 2.0 * data`).

Two scenarios (both gated, plan §9):

- `b5-vcatchment.yaml` — **rain**: 100 mm/h for 20 h, then 100 h of
  recession; discharge peaks near the full-domain rain rate ~1111 m³/h.
- `b5-vcatchment-norain.yaml` — **no-rain**: pure drainage/seepage from the
  initial water table (2 m below terrain); discharge is a few m³/h over
  120 h.

Both run in both coupling modes at the gate (`coupling.mode: sync`
committed; the harness patches `subcycled` — plan §9: the subcycled mode is
gated here because no subcycled golden exists anywhere).

## Reference data

`legacy/benchmarks/b5-vcatchment/{discharge,ponding}/*.csv` (the
development-side legacy archive; not distributed with this repository) —
digitized ParFlow/CATHY/HGS/Cast3M curves (time [h] vs discharge [m³/h] or
ponded storage [m³]; Cast3M has no rain-scenario discharge curve). The gate
(`tests/regression/run_regression.py b5`) applies the plan §9 envelope
rules; `tests/regression/tolerances/b5-vcatchment.yaml` documents them and
the measurement (mass-audit discharge, ±900 s smoothing).

## Authoring decisions adjudicated at the P3 gate

- **Raster orientation**: frehg2 reads the first data row as j = 0 while
  SERGHEI reads rasters north-up; the SERGHEI polygon/monitor coordinates
  flip as y → 55 − y (the rasters themselves are shared byte-for-byte).
  The channel is rows j = 0..4 with the outlet corner at (j = 0, i = 0).
- **Subsurface mesh**: `terrain_layers: uniform` (amendment A12) — the
  reference's terrain-parallel 5 m slab (SERGHEI `height: 5.0`). The legacy
  scaled-terrain wedge leaves the outlet columns ~0.5 m deep and the
  no-rain scenario without a water table at the channel.
- **Outlet**: a stage reservoir at the channel-end bed (kind eta, 0.0) on
  the five channel-end cells (i = 0, j = 0..4). The P0 provisional −10 m
  sink differs only in the (irrelevant) held value; the amendment-A2
  `outflow` kind deadlocks on the flat channel end (ghost stage = cell
  stage gives no gradient — measured: a 3 m lake filled the domain).
  Discharge is measured from the mass-audit boundary outflow (the P1
  finding: stored velocities under-read discharge).
- **The outlet faces must be a free outfall** (the strip roughness below):
  a face from a wet cell into a bed-level reservoir carries the full stored
  head as its local gradient, and the legacy point-implicit drag
  linearization (D = 1/(1 + ½ dt C_D |u^n|/h)) surge/stalls whenever the
  equilibrium drag number ½ dt C_D |u_eq|/h exceeds ~1 — for this waterfall
  geometry that means any face Manning n ≳ 0.35. Measured at n = 2.0: the
  outlet faces under-convey ~60× (0.002 m/s against a 0.92 m head), the
  channel impounds an artificial flat lake whose storage mimics the
  reference ponding, and rain-end releases it catastrophically (discharge
  spike +51 % of the envelope peak in ~20 min, ponding 671 → 169 m³).
  SERGHEI's raster lubricated its own outflow cells the same way
  (corner n = 0.001 for its edge-outflow scheme).
- **wetting_face_depth 1e-3**: sub-mm rain films must stay wet for
  infiltration (min_depth 1e-6) but must not advect momentum across
  hairline faces (SERGHEI's `dryDepth` 1e-6 is the min_depth analogue;
  b6's legacy wtfh was 5e-4).
- **Controller knobs** `courant_max: 10`, `dq: 0.05/0.1` (legacy defaults
  2 and 0.01/0.02): the legacy controller pins the common adaptive step
  near 0.4 s for the entire 20 h rain window (front cells passing near
  saturation carry dK/dθ ~ 1 m/s; capacity-limited ponds drain at ~Ks).
  The relaxed knobs hold it near 2 s; the envelope gate and the volume
  audits adjudicate the accuracy (the subsurface budget closes to 1e-10 of
  the initial storage either way).
- **Roughness authoring** (input/roughness.dat): the wall column keeps the
  raster's 16.0; the flow surfaces are reconstructed against the published
  envelope itself, anchored by measured runs:
  - **outlet strip (i = 0, j = 0..4): n = 0.1** — the free-outfall
    lubricant (see the outlet bullet; the raster's 16.0/0.001 belong to
    SERGHEI's own edge-outflow scheme, and its 0.001 corner is the same
    lubricant idea — 0.1 keeps the outfall out of the surge/stall regime
    while damping the near-frictionless runaway that 0.001 produces
    against a reservoir head).
  - **channel rows (j = 0..4, i ≥ 1): n = 6.0** — essentially the raster's
    own 6.26. The early rejection of 6.26 ("needs a 1.8 m channel") was
    measured *through the choked outlet*: with the free outfall the
    channel takes its kinematic depth (~0.3 m at the plateau) and the
    raster value is vindicated.
  - **hillslope: n = 1.45** (raster: 0.626). At 0.626 the measured ponding
    plateau is half the envelope (−44 %) with discharge dead-center — the
    storage deficit is distributed film depth, the one observable the
    hillslope roughness controls without touching the discharge plateau
    (which equals the rain rate in any steady state).
  - Method: two anchor runs (channel/hillslope 2.0/0.626 → plateau
    337 m³; 5.0/1.25 → 540 m³) calibrate the kinematic scaling
    (film depth ∝ n^0.6, predicted 543 vs measured 540); the third run
    (6.0/1.45, predicted ~600) landed at 596 vs the 605 m³ envelope mean
    (−1.5 %) with discharge unchanged (peak +2.6 %, 0/48 outside). This is
    reference reconstruction against the benchmark's defining curves, not
    model calibration — the discharge curve is insensitive to the storage
    knobs (its plateau is mass balance, its ramp/recession stayed
    dead-center across the whole reconstruction range), so the ponding
    reconstruction cannot be tuned to fake the discharge gate.
- The `reallocation_surplus` default (`drop`) stands: the coupled runs'
  reallocation activity is ~0.01 % of the rain volume (the A7 adjudication
  concern does not bind here).

The coupled exchange itself — including the supply-limited classification
this case forced (amendment A13) — is documented in
`docs/theory/exchange-flux.md`.
