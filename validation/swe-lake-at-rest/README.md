# swe-lake-at-rest — 2D well-balancing (C-property)

Still water over irregular 2D bathymetry, from SERGHEI
`ShallowWater/Analytical/LakeAtRest1` (200×200 cells, dx = dy = 1 m). The lake
surface sits at eta = 3.0 m; the tallest peaks (bed up to 4.0 m) emerge dry.
The test is the 2D C-property: velocities must stay at machine zero.

A sibling variant **LakeAtRest2** exists in the SERGHEI suite (same 200×200
grid, a different bathymetry field); only LakeAtRest1 is converted here — the
conversion is identical up to `input/dem.dat`.

## Conversion notes

- **Initial condition.** eta = bed + `hini` = 3.0 m wherever wet (verified
  constant) → `initial_conditions.surface.eta` constant 3.0.
- **Frictionless.** SERGHEI `friction: none` (its listed roughness 0.035 is
  inert with friction off) → Manning coefficient 0.0.
- **Boundaries.** `BCtype REFLECTIVE` → no `boundary_conditions` entries.
- **Raster row order (important).** This is genuine 2D relief, so the SERGHEI
  north-up DEM was **row-reversed** into `input/dem.dat` (frehg2 treats the
  first data row as j = 0). The raster XLLCORNER/YLLCORNER (150, 150) is inert;
  frehg2 places the grid at origin (0, 0).
- **Time step.** dt = 0.5 s (semi-implicit free surface; advective CFL
  non-binding at rest).

## Scheme caveat

frehg2's surface solver is **semi-implicit** versus SERGHEI's explicit
scheme; for still water the difference is minor.

## Cost (OMP_NUM_THREADS=1, ~1e-7 s per cell-step)

40000 cells × (100 / 0.5) = 8.0e6 cell-steps ≈ **0.8 s — LOCAL**.

## Validate

```
cd validation/swe-lake-at-rest
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swe-lake-at-rest.yaml  # VALID:
```
