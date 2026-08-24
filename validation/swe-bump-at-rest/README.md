# swe-bump-at-rest — well-balancing / C-property (still water over a bump)

Two SWE-only cases in one directory, both from the SERGHEI SWASHES
"bumps/rest" suite (`ShallowWater/SWASHES/steady_state/bumps/rest/{emerged,
immersed}`):

- **`emerged.yaml`** — lake surface at eta = 0.1 m; the parabolic crest
  (bed up to 0.2 m) pokes above the water and stays dry.
- **`immersed.yaml`** — lake surface at eta = 0.5 m; the bump is fully
  submerged.

Both are still water at rest over a variable bed on a 1000×4 strip
(25 m × 0.1 m, dx = dy = 0.025 m). The test is the C-property: a
well-balanced scheme keeps velocities at machine zero indefinitely (no
spurious currents from the bed-slope / pressure balance).

## Conversion notes

- **Initial condition.** SERGHEI ships the depth field `hini.input`; adding it
  to the bed gives a flat free surface (eta = 0.1 m emerged, 0.5 m immersed,
  both verified constant where wet). Encoded as `initial_conditions.surface.eta`
  constant. Where the bed rises above the lake level the depth is negative and
  frehg2 clamps the cell dry — the correct lake-at-rest behaviour.
- **Frictionless.** SERGHEI `friction: none`; frehg2 has no "none" law, so the
  Manning law is used with coefficient 0.0.
- **Boundaries.** SERGHEI `BCtype REFLECTIVE` (closed box) → no
  `boundary_conditions` entries.
- **Raster.** The bump varies only in x, so the DEM rows are identical and the
  SERGHEI north-up raster is copied verbatim (`input/dem.dat`, shared — the two
  cases have byte-identical DEMs).
- **Time step.** frehg2's free surface is solved **semi-implicitly**, so
  gravity waves are not step-limited and the advective CFL is non-binding at
  rest; dt = 0.1 s is used (well above the explicit acoustic limit ~0.013 s).

## Scheme caveat

frehg2's surface solver is **semi-implicit**, versus SERGHEI's explicit
Roe/shock-capturing scheme. For this still-water case the difference is minor;
both should preserve the rest state.

## Cost (OMP_NUM_THREADS=1, ~1e-7 s per cell-step)

| case | cells | t_end / dt | cell-steps | est. wall | tier |
|---|---|---|---|---|---|
| emerged  | 4000 | 10000 / 0.1 | 4.0e8 | 156 s (measured) | LOCAL |
| immersed | 4000 | 1000 / 0.1  | 4.0e7 | 29 s (measured)  | LOCAL |

## Validate

```
cd validation/swe-bump-at-rest
OMP_NUM_THREADS=1 ../../build/src/frehg --validate emerged.yaml   # VALID:
OMP_NUM_THREADS=1 ../../build/src/frehg --validate immersed.yaml  # VALID:
```
