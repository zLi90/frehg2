# swe-vcatchment — surface-only V-catchment overland flow

Classic V-catchment (Di Giammarco / "tilted-V") overland-flow benchmark, from
SERGHEI `Hydrology/Analytical/vcatchment/input`. Two side planes drain into a
central channel that slopes to a single outlet at the south edge; uniform
rainfall produces a rising-then-falling outlet hydrograph compared against the
analytical (kinematic) solution (`ReferenceData/analytical.csv`).

> **Note:** this is the **surface-only** overland-flow case. It is distinct from
> the coupled surface–subsurface tilted-V benchmark **b5-vcatchment**; here
> there is no groundwater and no infiltration.

162×100 cells, dx = dy = 10 m (1620 m × 1000 m).

## Conversion notes

- **Rainfall.** SERGHEI 10.80 mm/h for 0–5400 s then dry, converted to m/s
  (÷3.6e6 = 3.0e-6 m/s) as a series `input/rain.dat`. frehg2 time series are
  piecewise-linear, so a 1-second step-down (5400 → 5400.01 s) reproduces the
  SERGHEI pulse; a trailing 0 beyond t_end closes the series.
- **Spatial Manning.** SERGHEI `roughness: file` (0.015 on the planes, 0.15 in
  the central channel) → `friction.coefficient.file` (`input/roughness.dat`).
- **Infiltration.** SERGHEI `InfModel: none` → no subsurface loss; this is a
  pure SWE run (frehg2 has no surface-only infiltration term).
- **Outlet.** SERGHEI outer `BCtype REFLECTIVE` (closed) with one `extbc`
  outlet `bctype 5` (free outflow) → `kind: outflow` polygon over the 3 channel
  cells (x 794–824 m) at the south edge. The outlet is the lowest bed (~0.096 m
  at x ≈ 805 m); after the DEM row-reversal it lands on j = 0.
- **Raster row order (important).** Genuine 2D relief → the SERGHEI north-up DEM
  **and** roughness rasters were row-reversed into `input/dem.dat` /
  `input/roughness.dat`. (The roughness happens to be y-invariant; reversed
  anyway for consistency.) The raster YLLCORNER (−0.2 m) is inert; grid origin
  is (0, 0), and the outlet polygon y-span was set to capture the j = 0 row.
- **Dry start.** SERGHEI `initialMode dry` → eta = −10 m (far below the bed).
- **min_depth.** SERGHEI `dryDepth 0` is below the schema minimum → 1e-6.
- **Time step.** dt = 2.0 s; overland velocities are small (~0.3 m/s) over
  dx = 10 m, so the advective CFL is far from binding.

## Scheme caveat

frehg2's surface solver is **semi-implicit** versus SERGHEI's explicit scheme.
This is subcritical overland flow, so the scheme difference is minor.

## Resolved limitation: the south outlet throttle (V2-A11, fixed in v2.0.0)

Through v1 this case's south-edge `kind: outflow` BC ran through the
defective west/south ghost rule (Q0.3 diagnosis, amendment V2-A11): the
boundary face was given the *interior* face area, gauged over the higher of
the two beds, and the transmissive volume goes as that area squared — so the
outlet was throttled `(0.0121/0.192)² ≈ 1/252` and the channel outlet cells
held a **0.192 m pool** against the 0.2 m upslope bed step (38.4 of the
39.7 m³ standing in row j = 0 at `t_end`), inflating the recession tail
(volume 7719 → 1151 → 690 → 451 → 318 → 253 m³ over the last outputs).

**v2.0.0 (the Q7 fix, commit `def8fba`) removes the throttle.** Rerun
2026-09-26: the outlet row drains freely — the channel outlet cells hold
7.4e-3 m (bed-step pools gone), row j = 0 carries 2.85 m³ total, and the
recession trace at the same output times is 462 → 342 → 285 → 245 → 224 m³.
The hydrograph against the Maxwell et al. (2014) intercomparison (regenerated
figure): frehg2 peak 292.8 m³/min at t = 92 min; mean |frehg2 − model| =
28.5 (ParFlow) / 19.2 (HGS) / 13.2 (tRIBS) m³/min. The behavior is gated by
[`../swe-outflow-staircase/`](../swe-outflow-staircase/README.md)
(`regression.outflow_staircase`, all four edges).

## Known limitation: a small unaudited sink in deep thin-film recession

With the throttle gone this case now exercises a regime no gate covers —
the whole domain draining as a ~1e-4 m film for hours — and the mass audit
does not fully close there: by `t_end` the identity
`volume = rain − evaporation − boundary_outflow + bc_inflow + clamped`
carries a **−70.8 m³ residual (0.27 % of rain)**, accumulated only after
t ≈ 15000 s (rain-phase and peak-flow windows close to a few m³). Measured
properties (2026-09-26): the `volume` column matches the field integral of
the depth snapshots exactly, so a real sink goes unaudited; all cells stay
wet (no dry-out truncation; the clamp column stays frozen); the loss is
interior-distributed and roughly constant per step (~0.011 m³/step). It is
**pre-existing, not the V2-A11 fix** — a stock (pre-fix) binary shows the
same signature (−35.5 m³; smaller only because the throttled outlet passes
about half the late thin-film flux). For calibration: b4's east-edge outflow
audit closes to 1.7e-4 m³ and the staircase's steady-state budget to
0.075 %. Affects late-recession storage numbers read from the audit, not
the hydrograph (which is d/dt of the audited `boundary_outflow` itself).

## Cost (OMP_NUM_THREADS=1, ~1e-7 s per cell-step)

16200 cells × (35000 / 2.0) = 2.8e8 cell-steps ≈ **~30 s — LOCAL**.

## Validate

```
cd validation/swe-vcatchment
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swe-vcatchment.yaml  # VALID:
```
