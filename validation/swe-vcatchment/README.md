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

## Cost (OMP_NUM_THREADS=1, ~1e-7 s per cell-step)

16200 cells × (35000 / 2.0) = 2.8e8 cell-steps ≈ **~30 s — LOCAL**.

## Validate

```
cd validation/swe-vcatchment
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swe-vcatchment.yaml  # VALID:
```
