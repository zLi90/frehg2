# swe-macdonald-subcritical — MacDonald/SWASHES subcritical flume

Steady subcritical flow in a varying-geometry flume with a MacDonald-type
analytic bed, from SERGHEI
`ShallowWater/SWASHES/steady_state/flumes/long-channel/subcritical/input`. A
1000×4 strip (1000 m × 4 m, dx = dy = 1 m). Constant unit discharge enters at
the west; a fixed downstream depth holds the east. Validated against the
MacDonald/SWASHES analytical steady profile.

## Conversion notes

- **Inflow / outflow.** SERGHEI outer `BCtype TRANSMISSIVE`; the imposed BCs
  from `extbc.input` (bcvals `0.748378 2.0 0.0` = [downstream h, inflow q, –]):
  - west `bctype 7` → `kind: discharge`. **Unit conversion:** SWASHES gives a
    *unit* discharge q = 2.0 m²/s, but frehg2's `kind: discharge` value is a
    *total* volumetric discharge [m³/s] over the inflow-polygon width. Width =
    ny·dy = 4·1 = 4 m, so the frehg2 value is 2.0 × 4 = **8.0 m³/s**. Polygon at
    x ≈ 0. (Verified: the outlet recovers q = h·u ≈ 1.99 m²/s, <1% from target.)
  - east `bctype 6` (imposed depth **0.748378 m**, confirmed) → `kind: eta`.
    **Judgment call:** frehg2 `kind: eta` is a free-surface *elevation*, so the
    east bed (0.00572192 m) is added → **eta = 0.75409992 m**.
- **Initial condition.** SERGHEI `initialMode file` ships the MacDonald depth
  field; eta = bed + hini is written to `input/eta.dat`.
- **Friction.** Manning n = 0.033.
- **Raster.** 1D strip (bed varies in x only) → DEM copied verbatim.
- **Time step.** dt = 0.2 s ≈ the explicit advective CFL limit (u = q/h ≈
  2.7 m/s over dx = 1 m → 0.5·dx/u ≈ 0.19 s); advection is explicit.

## Scheme caveat

frehg2's surface solver is **semi-implicit** versus SERGHEI's explicit
Roe/shock-capturing scheme. This case is subcritical and steady, and the flow
is ~0.75 m deep (well above the thin-layer friction band), so agreement is
good: the outlet discharge matches the imposed q to <1%.

## Cost (OMP_NUM_THREADS=1)

4000 cells × (2000 / 0.2) = 4.0e7 cell-steps ≈ **~7 s — LOCAL** (measured).

## Validate

```
cd validation/swe-macdonald-subcritical
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swe-macdonald-subcritical.yaml  # VALID:
```
