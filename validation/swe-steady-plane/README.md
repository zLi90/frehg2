# swe-steady-plane — steady flow toward normal depth on a tilted plane

Two SWE-only cases in one directory, from SERGHEI
`ShallowWater/Analytical/steadyPlane/{slope2,slope8}`. Each is a 100×1 plane
(1000 m × 10 m, dx = dy = 10 m) with a constant unit discharge entering at the
high (west) end; the flow relaxes to a uniform steady depth.

- **`slope2.yaml`** — 2% slope, bed 19.9 → 0.1 m, outer `TRANSMISSIVE`,
  downstream **imposed depth** at the east. **Manning n = 0.033** (see the
  friction note below — SERGHEI's slope2 deck flags `friction: none`, which we
  override).
- **`slope8.yaml`** — **Manning n = 0.033**, bed 79.6 → 0.4 m (8% slope), outer
  `REFLECTIVE`, **free outflow** at the east. The imposed q and Manning n give a
  normal depth ≈ 0.01 m (q = (1/n)·h^(5/3)·√S, S = 0.08), confirming the SERGHEI
  bcvals.

## Conversion notes

- **bcvals convention.** SERGHEI `extbc.input` bcvals = [downstream h, inflow q,
  –]. slope2: `0.01 0.00196306 0.0`; slope8: `0.01 0.003973 0.0`.
  - inflow `bctype 7` → `kind: discharge`. **Unit conversion:** bcvals[1] is a
    *unit* discharge [m²/s] but frehg2's `kind: discharge` value is a *total*
    volumetric discharge [m³/s] over the inflow-polygon width. Width = ny·dy =
    1·10 = 10 m, so the frehg2 values are 0.00196306·10 = **0.0196306** (slope2)
    and 0.003973·10 = **0.03973** (slope8) m³/s.
  - slope2 outflow `bctype 6` (imposed depth 0.01 m) → `kind: eta`. **Judgment
    call:** SERGHEI prescribes a depth but frehg2's `kind: eta` is a
    free-surface *elevation*, so the downstream bed (0.1 m) is added →
    **eta = 0.11 m**. (Using the bare 0.01 would sit below the bed and dry the
    outlet.)
  - slope8 outflow `bctype 5` (free outflow) → `kind: outflow` (no value); the
    stage is extrapolated down the continued bed slope at the east faces.
- **Outer boundary.** slope2 `TRANSMISSIVE` and slope8 `REFLECTIVE` both reduce
  to the same explicit inflow+outflow polygons here (the single row has no
  transverse flux); the difference is only in how the un-prescribed east face is
  treated, captured by outflow `kind eta` (slope2) vs `kind outflow` (slope8).
- **Friction (slope2 deviation).** SERGHEI's slope2 deck sets `friction: none`
  but leaves `roughness: 0.033`, and its `extbc` comments the inflow as the
  *"normal specific discharge"* with the *"steady water depth"* downstream — i.e.
  q and h = 0.01 m are the Manning normal-flow values for n = 0.033 at S = 0.02.
  With friction genuinely off that pairing is degenerate (no normal depth exists;
  the plane drains dry — verified). We therefore run slope2 with **n = 0.033**,
  the roughness the discharge was derived from and the value slope8's deck
  applies. This is a deliberate deviation from the slope2 deck's `friction: none`
  flag, treated as a deck inconsistency.
- **Initial condition.** SERGHEI `initialMode h` = uniform depth 0.02 m over the
  sloping bed, so eta = bed + 0.02 was written as a field
  (`input/eta_slope2.dat`, `input/eta_slope8.dat`). Single-row strips → no DEM
  row-reversal needed (verbatim copies `input/dem_slope2.dat`,
  `input/dem_slope8.dat`).
- **Time step.** dt = 1.0 s for both; u = q/h is ~0.1 m/s (slope2) / ~0.4 m/s
  (slope8) over dx = 10 m, so the advective CFL is far from binding.

## Scheme caveat & observed accuracy

frehg2's surface solver is **semi-implicit** versus SERGHEI's explicit scheme.
Both cases are subcritical and steady. frehg2 reproduces the correct
qualitative solution — a **spatially uniform depth** across the plane (flat
water surface parallel to the bed, the signature of normal flow) — but the
converged depth/discharge sit ~20% off the SWASHES/Manning analytic
(h ≈ 0.012 vs 0.01 m at default settings). The offset traces to the ~1 cm flow
depth lying inside frehg2's `thin_layer_depth` friction-regularization band
(default 0.1 m ≫ h); the deeper MacDonald flume case, with h ≈ 0.75 m, matches
its target discharge to <1%. Treat these two as qualitative normal-depth
demonstrations rather than tight quantitative benchmarks.

## Cost (OMP_NUM_THREADS=1)

| case | cells | t_end / dt | cell-steps | wall | tier |
|---|---|---|---|---|---|
| slope2 | 100 | 20000 / 1.0 | 2.0e6 | ~1 s (measured) | LOCAL |
| slope8 | 100 | 3000 / 1.0  | 3.0e5 | <1 s (measured) | LOCAL |

Each variant writes its own output (`out/slope2.h5`, `out/slope8.h5`).

## Validate

```
cd validation/swe-steady-plane
OMP_NUM_THREADS=1 ../../build/src/frehg --validate slope2.yaml  # VALID:
OMP_NUM_THREADS=1 ../../build/src/frehg --validate slope8.yaml  # VALID:
```
