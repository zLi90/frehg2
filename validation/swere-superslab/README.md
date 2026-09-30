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

## Findings (Q0.3, 2026-09-21)

This case was flagged on 2026-08-30 as a broken run — "15 m³ of rain injected,
zero ponding/outflow/seepage" — and became plan item **Q0.3**. Both halves of
that flag are now resolved; the details and the derivation are in plan
amendment **V2-A11**.

### The "broken run" was a partially written file

`out/output.h5` was read by `makeplot.py` while the 12 h job was still
running. Early in a dry-start run ponding, outflow and seepage *are* all ~0,
which is exactly the degenerate-run signature the script tests for, so it
reported a mass-balance breakdown that did not exist. The completed file was
a 43 200 s run (86 401 `mass_audit` rows), its embedded config
byte-identical to `swere-superslab.yaml`, and its budget closed:

```
volume 0.1071 = rain 15.0247 - evap 0 - outflow 1.4311 + bc_inflow 0
              + seepage (-13.4870) + clamped 5.6e-11      (residual -5.4e-4 m3)
```

`makeplot.py` now checks `t_end` from the embedded config **before** applying
the degenerate test, so "not finished yet" can no longer be reported as
"broken". The residual formula was also wrong (it subtracted the signed
seepage and omitted `bc_inflow`/`clamped`), which is why it printed ~27 m³
unaccounted on a run that closes to 5e-4 m³.

### The rain source path — what Q0.3 gates for Q4 — is exact

`rain.dat` applies 1.388889e-05 m/s over 100 m² from t = 0 to 10 800 s,
ramping to 0 at 10 836 s: `1.3889e-5 x 100 x 10818 = 15.025 m3`, matching the
recorded `rain` of 15.02465 m³. Nothing about the rainfall or evaporation
source path is implicated.

### A real defect this case did expose: the west/south transmissive outlet

In the pre-fix run, the residual 0.107 m³ of surface storage at `t_end` was
a **0.1034 m pool standing in the outlet cell** (`i = 0`), level at
eta ~ 0.1035 — pinned just
below the upslope neighbour's bed at 0.1 m. The reference codes hold ~2e-4 m³
there, which is the Manning normal depth `h = (n q / sqrt(S))^(3/5) = 2.7e-4 m`
for this outlet.

The cause was the west-edge ghost rule at v1.x `src/swe/WetDry.cpp:144`
(`Asx(j, 0) = Asx(j, 1)`, ported from `shallowwater.c:1067-1099`). It hands
the boundary face the **interior** face area, which is gauged over the higher
of the two beds — the 0.1 m sill — instead of over the outlet cell's own bed.
The transmissive BC's released volume (`FreeSurface.cpp:311`) is proportional
to that area squared, so the face is throttled by
`(deptx(j,1) / depth(j,0))^2 ~ 1/812` here and the cell cannot discharge until
it fills to the sill.

Note the override *destroys* a correct value: the face kernel
(`WetDry.cpp:49-57`) already computes `deptx(j,0)` over the halo column, and
because the bed ghost is a zero-gradient copy (`SurfaceSolver.cpp:214`) and
the outflow ghost sets `eta(j,0) = eta(j,1) - drop`, that value is exactly the
outlet cell's own depth. Deleting the west and south copies is sufficient.

**Fixed in v2.0.0 (Q7).** Q0.3's deliverable was the diagnosis; the fix waited
behind the x-gate that the plan's §6.1 gate-first rule requires for BC kind ×
side (§1.3). That gate is `regression.outflow_staircase`, built from the 10×1
reproducer below: 10 failures on the stock code, PASS once the two ghost
assignments were deleted (`WetDry.cpp::updateGeometry`). V2-A11 has the
measurements; V2-A19 records the landing.

### Scope of the defect

Only the **west (-x)** and **south (-y)** edges were affected; east/north
transmissive faces use `Sxp`/`Syp`, the cell's own coefficient, and were
correct. Of the five cases in this repo that use `kind: outflow`, three are on
the east edge — including **b4-govindaraju, the only gate that exercised the
BC at all**, which is why b1-b6 never caught this. The two on the defective
path are both validation cases: this one (west) and `swe-vcatchment` (south),
where the same signature appeared independently — a 0.192 m pool against that
case's 0.2 m channel bed step, holding 38 of the 39.7 m³ left in the outlet
row (gone after the fix; see that case's README).

### The record (v2.0.0 release binary)

Rerun after the fix on the release binary (`2.0.0` at `f56ba96`,
2026-09-29 → 30); `makeplot.py` regenerated both figures from it. The
pre-fix column is the v1.x run the findings above come from. The
hydrograph ranges span ParFlow, CATHY and Cast3M; the digitized HGS
outflow series stays below 0.002 m³/h and is not in them.

| signal | frehg2 v2.0.0 | frehg2 pre-fix | reference range |
| --- | --- | --- | --- |
| outlet peak | 0.404 m³/h @ 8.52 h | 0.361 m³/h @ 8.93 h | 0.288-0.456 m³/h @ 7.50-8.44 h |
| return-flow onset | 6.54 h | 7.20 h | 4.78-6.51 h |
| rain-phase ponding | 1.10e-2 m³ over i = 40..60 | 1.10e-2 m³ | 5.36e-3 (HGS) - 6.66e-3 (Cast3M) |
| outlet storage, late | 1.40e-4 m³ | 1.07e-1 m³ | ~2e-4 m³ |

The outlet pool is gone: late outlet storage falls from 0.107 m³ to
1.40e-4 m³, the normal-depth scale the reference codes hold. The throttled
outlet had also delayed the hydrograph. Onset and peak now come 0.66 h and
0.41 h earlier, 0.03 h and 0.08 h after the latest reference code
(ParFlow, 6.51 h and 8.44 h), and the peak is inside the range; the curve
is closest to ParFlow's, with a peak about 10 % lower. Rain-phase ponding
over the low-K block is unchanged, since no water reaches the outlet
before 6.54 h. It stays 1.7-2.1x the reference range, in the right place.
The surface budget closes to +1.3e-3 m³ at `t_end` (8.5e-5 of the rain),
and 1.800 m³ leaves through the outlet against 1.431 m³ pre-fix.

A minimal 10x1 reproducer isolates the defect cleanly, and is in-tree at
[`../swe-outflow-staircase/`](../swe-outflow-staircase/README.md): fed at the
outlet cell, the stock code traps 0.1092 m³ and passes **zero** outflow until
t = 1001.5 s of a 4000 s run; with lines 144/158 removed the same case settles
at 2.039e-4 m depth — the analytic normal depth scale — and passes 0.3997 m³ of
the 0.400 m³ injected. The companion upslope-fed case shows the same wrong area
failing in the opposite direction (38.3 % over-drain, 39.3 % of the input
minted by the below-bed clamp).

### Secondary observations (not defects of this case)

- `gw_mass_audit` shows `ss_storage` and `vloss` each reaching ~202 m³
  (202.1 and −201.3 m³ at `t_end`; pre-fix ~202 each), ~13x the 15.0 m³
  of rain, very nearly cancelling. This is specific-storage churn in dry
  unsaturated soil — the P2 `Ss` + dry-clay runaway, already documented.
- `realloc_dropped` = 6.312 m³, ~48 % of the 13.2 m³ infiltrated (pre-fix
  6.176 m³ of 13.5 m³): the documented lateral-reallocation drop (the
  b2/b3 `reallocation_surplus` fork).

## Caveats

- **Semi-implicit surface solver**: frehg2's SWE step differs from SERGHEI's
  explicit scheme.
- **Dynamic seepage coupling**: SERGHEI switches the top BC between ponding and
  seepage during the run; frehg2 approximates this with its built-in coupler.
- **Strong heterogeneity**: the soil `Ks` ratio spans ~1e4, which stresses the
  coupled infiltration front and the adaptive step.

## Cost

Measured, not estimated. On the v2.0.0 release binary (Apple M3, serial,
sharing the machine with other pipeline lanes for most of the run), the
sync common step holds at `dt = 0.5 s` (86 400 steps), and the groundwater
iterations per step grow as the profile drains: 4 at t = 1 h, 14 at 3 h,
21 at 4.5 h, 24 from 7 h to 11 h (the `gw N it` field of the log's output
lines). The 12 h run took 13.1 h of wall time (2026-09-29 17:18 →
09-30 07:36, less a 1.2 h pause while the release pipeline's scaling gates
needed an idle machine). The pre-fix run on 2026-08-30 took about a day,
which is why its file was read mid-write (Findings). The 15–30 min this
section used to give was an unmeasured estimate.

**About 13 h serial on a laptop → HPC or a long unattended run.** Run
`--validate` only for authoring.

## Validation

```
OMP_NUM_THREADS=1 ../../build/src/frehg --validate swere-superslab.yaml
# -> VALID: swere-superslab.yaml   (exit 0)
```
