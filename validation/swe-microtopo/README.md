# swe-microtopo — rainfall runoff over microtopography

Rainfall-runoff over synthetic sinusoidal microtopography, from SERGHEI
`Hydrology/Analytical/microtopo/<variant>/input`. A short rain pulse generates
overland flow that drains to an east-edge outlet; the outlet hydrograph is
compared against the analytical/reference solution. Three variants (one YAML
each, sharing this README):

| variant | file | grid | dx = dy [m] | domain [m] | slope |
|---|---|---|---|---|---|
| S0.01 (WL1.7357, A0.01) | `microtopo-s0.01.yaml` | 92×35  | 0.086957 | 8.0 × 3.04 | 1% |
| S0.02 (WL0.4143, A0.03) | `microtopo-s0.02.yaml` | 209×79 | 0.038278 | 8.0 × 3.02 | 2% |
| S0.05 (WL2.0000, A0.03) | `microtopo-s0.05.yaml` | 88×33  | 0.090909 | 8.0 × 3.00 | 5% |

## Conversion notes

- **Rainfall.** SERGHEI 7.5 mm/h for 0–1800 s then dry, converted to m/s
  (÷3.6e6 = 2.0833e-6 m/s) as a series (`input/rain_s0.0*.dat`) with a
  1-second step-down and a trailing 0 beyond t_end. (The S0.01 source file lists
  only two points, relying on SERGHEI's stepwise hold; since frehg2 series are
  piecewise-linear, all three variants use the explicit constant-pulse idiom.)
- **Friction.** Uniform Manning n = 0.055.
- **Infiltration.** none → pure SWE run.
- **Outlet.** SERGHEI outer `BCtype REFLECTIVE` (closed) with one `extbc`
  outlet `bctype 5` (free outflow) → `kind: outflow` polygon on the east edge
  (spanning the full y-extent). The SERGHEI outlet polygons for S0.02/S0.05 are
  a few cells wide, but frehg2 applies `outflow` only at genuine domain-edge
  faces (`appendSideFaces`), so only the east column receives an outflow face —
  the extra width is harmless.
- **Raster row order (important).** Genuine 2D relief → each SERGHEI north-up
  DEM was row-reversed into `input/dem_s0.0*.dat` (frehg2 first data row = j=0).
- **Dry start.** SERGHEI `initialMode dry` → eta = −10 m (below all beds; the
  microtopo beds dip to about −0.4 m in S0.05).
- **Time step.** dt = 0.1 s for S0.01 and S0.05; **dt = 0.05 s for S0.02**,
  whose finer grid (dx ≈ 0.038 m) needs the smaller step for the explicit
  advection term. The free surface is semi-implicit.

## Scheme caveat

frehg2's surface solver is **semi-implicit** versus SERGHEI's explicit
Roe/shock-capturing scheme. This is subcritical rainfall-runoff, so the scheme
difference is minor.

## Cost (OMP_NUM_THREADS=1, ~1e-7 s per cell-step)

| variant | cells | t_end / dt | cell-steps | est. wall | tier |
|---|---|---|---|---|---|
| S0.01 | 3220  | 8000 / 0.1  | 2.6e8 | 26 s (measured)   | LOCAL |
| S0.02 | 16511 | 8000 / 0.05 | 2.6e9 | 6.8 min (measured) | LOCAL |
| S0.05 | 2904  | 8000 / 0.1  | 2.3e8 | 21 s (measured)   | LOCAL |

## Validate

```
cd validation/swe-microtopo
OMP_NUM_THREADS=1 ../../build/src/frehg --validate microtopo-s0.01.yaml  # VALID:
OMP_NUM_THREADS=1 ../../build/src/frehg --validate microtopo-s0.02.yaml  # VALID:
OMP_NUM_THREADS=1 ../../build/src/frehg --validate microtopo-s0.05.yaml  # VALID:
```
