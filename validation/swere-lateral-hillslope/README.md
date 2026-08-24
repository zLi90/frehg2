# swere-lateral-hillslope — Coupled lateral hillslope subsurface flow

Converted from the SERGHEI-SWE-RE reference case `case4/`
(`serghei-swe-re-testcases/case4`, SERGHEI supplement; not shipped): a Hydrus2D/analytic-referenced
coupled surface–subsurface **lateral hillslope** run over 5 years
(157 680 000 s). Lateral drainage happens in the saturated subsurface between
two fixed-head boundaries; the surface only handles light rainfall and any
exfiltration through the built-in coupler.

## Geometry and mesh

- One-cell-wide **x–z strip**: `nx = 40`, `ny = 1`, `nz = 60`; `dx = dy = 100 m`
  (`cellsize`), `dz = 0.25 m` (`dz_base`; `60 × 0.25 = 15.0 m = height`).
- The DEM (`dem.input`) slopes 3.9 → 0.0 m (a gentle 0.001 grade). The 15 m
  subsurface is a **terrain-parallel slab**, so `follow_terrain: true` +
  `terrain_layers: uniform` (the b5 / amendment-A12 idiom): each column carries
  its own 60 × 0.25 m slab below its bed. The inclined hillslope is thereby
  approximated on a rectilinear grid.
- `ny = 1`, so the rasters are x-only and need no north-up row reversal. The 3D
  fields are flat `(j*nx+i)*nz+k` → `i*nz+k`, with `k = 0` the top layer.

## Conversion decisions

- **gw scheme**: SERGHEI `gw_scheme 1` (PCA) → `pca` (the sole frehg2 scheme).
- **Surface forcing**: `sw.input` dry start, Manning `n = 0.03`. There is **no
  `extbc.input`**, so the surface has reflective walls (`BCtype REFLECTIVE`) and
  is driven only by `rainfall.input` (light, 0.017–0.108 mm/h) plus coupled
  seepage. Rain converted mm/h → m/s (`/3.6e6`), hours → seconds.
- **Lateral head BCs**: `gwbc.input` sets `dir 2` (−x, "left") and `dir 1` (+x,
  "right") to `bctype 2` = prescribed head, each with a *depth-varying* head
  file (`headleft/headright.input`): hydrostatic below a water table, clamped to
  −1.25 m in the unsaturated zone. frehg2's head BC has **no per-layer profile
  form**, so each is represented as a **hydrostatic head BC** (`head = eta − z`)
  with `eta` = the water-table elevation recovered from the file's saturated
  slope: `eta_left = −4.1 m`, `eta_right = −14.1 m`. These reproduce the files'
  bottom heads (6.88 / 0.78 m) to ~1 %; the −1.25 m unsaturated clamp is not
  reproduced (negligible lateral flux there).
- **Top coupling**: `gwbc.input` "surface" (`dir 6`, `bctype 8`) marks the
  coupled top boundary; frehg2's coupler owns that exchange → **no explicit BC**.
- **Initial conditions**: subsurface `initialMode 2` (file-h) → `groundwater.head`
  from a 3D flat field built from `head.input` (`k = 0` = top). Surface
  `initialMode dry` → `eta = −10 m` (below every bed).
- **specific_storage**: absent from the SERGHEI input; set `1e-5` (schema
  requires it; the b5 value).

## Caveats

- **Semi-implicit surface solver**: frehg2's SWE step is semi-implicit and
  differs from SERGHEI's explicit scheme.
- **Dynamic seepage coupling**: SERGHEI switches the top BC between
  ponding and seepage during the run; frehg2 approximates this with its
  built-in surface–subsurface coupler.
- **Lateral head profile**: the depth-varying fixed-head files are collapsed to
  a single hydrostatic water-table elevation per side (see above).

## Cost estimate

Sync coupling: both modules march the common adaptive step, so the 3D
groundwater cell count dominates. Cells = `40 × 1 × 60 = 2 400`. For a
quasi-static 5-year run the controller settles near `dt_max = 90 s`
(assumption), giving ≈ `1.5768e8 / 90 ≈ 1.75e6` steps and
`2 400 × 1.75e6 ≈ 4.2e9` gw cell-steps × `2e-6 s ≈ 2.3 h` (surface adds ~7 min).

**≈ 2–2.5 h wall time → HPC** (well over the 10-min local threshold; the long
horizon, not the small grid, is the cost). Run `--validate` only for authoring.

## Validation

```
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swere-lateral-hillslope.yaml
# -> VALID: swere-lateral-hillslope.yaml   (exit 0)
```
