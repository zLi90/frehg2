# swere-superslab — Coupled heterogeneous superslab infiltration/return flow

Converted from the SERGHEI-SWE-RE reference case `case5/`
(`serghei-swe-re-testcases/case5`, SERGHEI supplement; not shipped): a heterogeneous multi-soil
vertical-slab ("superslab") infiltration and return-flow test, coupled
surface–subsurface, over 12 h (43 200 s). Rainfall infiltrates a layered
profile whose low-conductivity lens/block force return flow to the surface
outlet.

## Geometry and mesh

- **2-D vertical (x–z) cross section**: `nx = 100` (`dem.input n_cols`),
  `ny = 1` (`n_rows`), `nz = 100` (`ndepth`); `dx = dy = 1.0 m` (`cellsize`).
  `soilID.input` confirms the shape — 100 rows (depth layers `k`) × 100 columns
  (x cells `i`).
- **dz = 0.05 m**. `subsurface.input` lists `height 5.0`, `ndepth 100`, but
  `dz_base 1.0` is inconsistent (`100 × 1.0 = 100 ≠ 5.0`). The slab is 5 m over
  100 layers, and `wt.input` puts the water table exactly at the domain bottom
  (`wt = bed − 5.0` for every column), confirming `dz = height/ndepth = 0.05 m`;
  `dz_base 1.0` is treated as a stale field.
- The DEM rises 0.0 → 9.9 m; the 5 m slab is terrain-parallel →
  `follow_terrain: true` + `terrain_layers: uniform` (b5 / A12 idiom). The
  sloped domain is approximated on a rectilinear grid.
- `ny = 1`, so no north-up row reversal. 3D fields are flat `(j*nx+i)*nz+k` →
  `i*nz+k`, with row `r` of `soilID.input` mapped to layer `k = r` (`k = 0` =
  top/surface, matching frehg2's vertical index).

## Conversion decisions

- **Soils**: `vg.input` holds 3 semicolon-separated soils → three `soil.types[]`.
  The 3D soil map is built from `soilID.input` into `input/soil_id.dat`:
  - id 0 — high-K **background**, `Ks 2.70e-3` (α 6.0, n 2.0, θr 0.02);
  - id 1 — a sloping low-K **lens** dipping right and down, `Ks 6.94e-6`
    (α 1.0, n 3.0, θr 0.03);
  - id 2 — a near-surface very-low-K **block** (cols ~40–60, top ~28 layers),
    `Ks 2.78e-7` (α 1.0, n 3.0, θr 0.03).
  All three share θs = 0.1.
- **gw scheme**: SERGHEI `gw_scheme 1` (PCA) → `pca`.
- **Surface forcing**: `sw.input` dry start, Manning `n = 0.0036`. Driven by
  `rainfall.input` (50 mm/h for 3 h, then dry to 12 h; mm/h → m/s `/3.6e6`,
  h → s) plus the coupled exchange.
- **Surface outlet BC**: `extbc.input` "outlet", `bctype 5`, direction (−1, 0)
  = −x edge (cell `i = 0`), `polygon.input` → `kind: outflow` (free/transmissive
  outfall), `target: surface`.
- **GW BCs**: `gwbc.input` has only "surface" (`dir 6`, `bctype 8`) = coupled
  top → coupler owns it, **no explicit BC**. All other faces are no-flow
  (`BCtype REFLECTIVE`) = frehg2 default, so no bottom/side BC is written.
- **Initial conditions**: subsurface `initialMode 4` (file-wt) →
  `groundwater.water_table` from `wt.input` (2D raster; table at the domain
  bottom). Surface `initialMode dry` → `eta = −10 m` (below every bed).
- **specific_storage**: absent from the SERGHEI input; set `1e-5` (schema
  requires it; the b5 value).

## Caveats

- **Semi-implicit surface solver**: frehg2's SWE step differs from SERGHEI's
  explicit scheme.
- **Dynamic seepage coupling**: SERGHEI switches the top BC between ponding and
  seepage during the run; frehg2 approximates this with its built-in coupler.
- **Strong heterogeneity**: the soil `Ks` ratio spans ~1e4, which stresses the
  coupled infiltration front and the adaptive step.

## Cost estimate

Sync coupling: the 3D groundwater cell count dominates. Cells =
`100 × 1 × 100 = 10 000`. This is a dynamic infiltration event with
`dt_max = 1.0 s`; the wetting front through the low-K block drives the common
step below 1 s during the 3 h rain (assume an effective `dt ≈ 0.5–1.0 s`),
giving ≈ 43 200–86 400 steps and `10 000 × ~6e4 ≈ 6e8` gw cell-steps ×
`2e-6 s ≈ 15–30 min` (surface is negligible by comparison).

**≈ 15–30 min wall time → HPC** (over the 10-min local threshold; the small
`dt_max` over the fine 10 000-cell grid is the cost). Run `--validate` only for
authoring.

## Validation

```
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swere-superslab.yaml
# -> VALID: swere-superslab.yaml   (exit 0)
```
