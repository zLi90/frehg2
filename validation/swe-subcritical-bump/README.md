# swe-subcritical-bump — SWASHES steady subcritical flow over a bump

Steady subcritical flow over a parabolic bump (SWASHES §3.1 subcritical), from
SERGHEI `ShallowWater/SWASHES/steady_state/bumps/subcritical/input`. A 1000×4
strip (25 m × 0.1 m, dx = dy = 0.025 m). Constant unit discharge enters at the
west, a fixed downstream depth holds the east; the steady free-surface profile
is compared against the SWASHES analytical solution.

## Conversion notes

- **Inflow / outflow.** SERGHEI outer `BCtype TRANSMISSIVE`; the two open faces
  are the imposed BCs from `extbc.input` (bcvals `2.0 4.42 0.0` = [downstream h,
  inflow q, –]):
  - west `bctype 7` (imposed discharge) → `kind: discharge`. **Unit conversion:**
    SWASHES specifies a *unit* discharge q = 4.42 m²/s, but frehg2's
    `kind: discharge` value is a *total* volumetric discharge [m³/s] spread over
    the inflow-polygon width. Domain width = ny·dy = 4·0.025 = 0.1 m, so the
    frehg2 value is 4.42 × 0.1 = **0.442 m³/s**. (Confirmed by mass balance: the
    outlet recovers u = Q/(h·W) = q/h ≈ 2.2 m/s, and h·u along the channel
    reproduces q ≈ 4.42 m²/s.)
    Polygon at x ≈ 0.
  - east `bctype 6` (imposed depth 2.0 m) → `kind: eta`. The east bed is 0, so
    the stage equals the depth, **eta = 2.0 m**, polygon at x ≈ 25.
- **eta vs depth (judgment call).** SERGHEI `bctype 6` prescribes a *depth*;
  frehg2 `kind: eta` prescribes a free-surface *elevation* (verified in
  `SurfaceSolver.cpp`: the BC value is stage + a global datum offset, which is 0
  here since the bed ≥ 0). Here the two coincide because the downstream bed is 0.
- **Initial condition.** SERGHEI `initialMode file` ships a depth field that is
  exactly still water at the downstream level (eta = bed + hini = 2.0 constant),
  so `eta` is set constant 2.0; it relaxes to the steady profile under the BCs.
- **Frictionless.** `friction: none` → Manning coefficient 0.0.
- **min_depth.** SERGHEI `dryDepth 0` is below the schema minimum (must be > 0),
  so min_depth = wetting_face_depth = 1e-6.
- **Raster.** 1D strip (bump varies in x only) → DEM copied verbatim.
- **Time step.** dt = **0.001 s**. The semi-implicit free-surface solve diverges
  for this frictionless, started-from-rest transient at the SWASHES nominal
  0.01 s (CFL blows up around t = 0.5 s and the `fs_` linear system diverges);
  0.001 s is stable. At dt = 0.001 s the explicit advective CFL is ≈ 0.09
  (u = q/h ≈ 2.2 m/s, dx = 0.025 m).
- **End time.** t_end = **200 s**. This is a steady-state test; the profile is
  established by ~100 s (SERGHEI's nominal 1000 s is far past convergence). 200 s
  keeps the run local (~2–3 min) rather than ~11 min for the full 1000 s.

## Scheme caveat

frehg2's surface solver is **semi-implicit** versus SERGHEI's explicit
Roe/shock-capturing scheme. This case is subcritical and steady, so the scheme
difference is minor.

## Cost (OMP_NUM_THREADS=1, ~1.7e-7 s per cell-step measured here)

4000 cells × (200 / 0.001) = 8.0e8 cell-steps ≈ **~2.3 min — LOCAL** (measured).
The small dt is what makes this the most expensive SWE case in the collection;
the full SWASHES 1000 s horizon would be ~11 min (HPC territory) and is
unnecessary for a steady-state comparison.

## Validate

```
cd validation/swe-subcritical-bump
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swe-subcritical-bump.yaml  # VALID:
```
