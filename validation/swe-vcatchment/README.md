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

## Known limitation: the outlet cannot drain below the upslope sill (V2-A11)

This case's `kind: outflow` BC is on the **south (-y)** edge, which runs
through the defective west/south ghost rule diagnosed under plan item Q0.3
(`src/swe/WetDry.cpp:158`, `Asy(0, i) = Asy(1, i)`): the boundary face is
given the *interior* face area, gauged over the higher of the two beds rather
than the outlet cell's own bed. The transmissive volume goes as that area
squared, so the outlet is throttled until it fills to the upslope sill.

Measured in `out/output.h5`: the channel cells `i = 80, 81` at `j = 0` hold
**0.192 m** of water against a 0.2 m bed step along y, at a throttle of
`(0.0121 / 0.192)^2 ~ 1/252`. That pool is 38.4 of the 39.7 m³ standing in
row `j = 0` at `t_end`, and it is what the long recession tail in the `volume`
column is draining (7719 -> 1151 -> 690 -> 451 -> 318 -> 253 m³).

The recession limb of the hydrograph and any late-time storage number from
this case are therefore affected. The rising limb and peak, which are set
upslope, are not. See plan amendment **V2-A11** for the derivation, the
verified two-line fix, and why it is deferred behind an x-gate;
[`../swere-superslab/README.md`](../swere-superslab/README.md) has the full
write-up, and [`../swe-outflow-staircase/`](../swe-outflow-staircase/README.md)
is the minimal reproducer with an analytic answer.

## Cost (OMP_NUM_THREADS=1, ~1e-7 s per cell-step)

16200 cells × (35000 / 2.0) = 2.8e8 cell-steps ≈ **~30 s — LOCAL**.

## Validate

```
cd validation/swe-vcatchment
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swe-vcatchment.yaml  # VALID:
```
