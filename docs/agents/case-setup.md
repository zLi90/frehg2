# Case-setup guide for AI assistants

**Purpose:** everything an LLM needs to create correct Frehg2 simulation
cases — the YAML configuration plus its input data files — for a user's
problem. Use this file **together with the complete key reference**
[`docs/user-guide/parameters.md`](../user-guide/parameters.md): that
table is CI-locked to the schema source and is authoritative for every
key, type, default, and rule. This guide carries what the table cannot:
the workflow, the conventions, and the pitfalls.

## 1. The workflow that produces reliable configs

1. **Copy the nearest example case** (§3) into a new directory — never
   start from blank YAML.
2. Edit it for the user's problem, consulting `parameters.md` for every
   key you touch.
3. **Iterate the validate loop until clean**:
   `frehg --validate case.yaml` → prints `VALID:` or `INVALID:` with one
   line per problem. Unknown keys get nearest-key suggestions; types,
   ranges, required/one-of groups, cross-field rules, and input-file
   existence are all checked. Fix and repeat.
4. **Smoke-run a shortened horizon** (temporarily set `time.t_end` to a
   few hundred steps): confirm it runs to `run complete`, then check the
   mass audit closes (see `docs/agents/postprocessing.md` §5) and the
   fields look physical.
5. Restore the real horizon and run.

`--validate` cannot catch dynamics (a too-large `dt` blowing up
mid-run); the smoke run is what catches that.

## 2. Universal conventions

- **SI units everywhere**: meters, seconds, m/s, m²/s, m³/s. Convert the
  user's data (mm/h rain → m/s, etc.) at file-creation time.
- **One absolute elevation datum** shared by `bottom_elevation`, initial
  `eta`, stage boundary values, and hydrostatic heads. Water depth is
  `eta − bottom`; a cell is dry where `eta ≤ bottom`. **A dry start is
  expressed as `eta` well below the bed** (e.g. `constant: -10` for a
  bed near 0) — this is intentional, not an error.
- Grid: regular Cartesian; cell `(i, j)` center at `x = (i+0.5)·dx`,
  `y = (j+0.5)·dy`; `i` runs west→east, `j` south→north; subsurface
  layer `k = 0` at the land surface, increasing **downward**, thickness
  `dz·dz_stretch^k`.
- Flat-file ordering: 2D rasters `j·nx + i` (i fastest); 3D fields
  `(j·nx + i)·nz + k`. Headered rasters: **first data row is `j = 0`**
  (file order, *not* ESRI north-up).
- `time.t_start`, `t_end`, `output_interval` are **whole seconds**
  (outputs are keyed by integer seconds in the HDF5 file).
- Every relative path in the YAML resolves against the **YAML file's
  directory**; output paths resolve against the run's working directory.

## 3. Which example to copy

| Problem type | Copy | Notes |
|---|---|---|
| overland flow / rainfall runoff | `benchmarks/b1-sw/` or the worked example in `docs/user-guide/worked-example.md` | Manning friction, rain series, outflow/eta outlet |
| Chézy-friction channel/plane | `benchmarks/b4-govindaraju/` | `friction.law: chezy` |
| still-water / wetting-drying tests | `validation/swe-lake-at-rest/`, `swe-bump-at-rest/` | well-balancing checks |
| 1D/2D infiltration, drainage, water table | `benchmarks/b2-gw/` (column), `validation/re-2d-open-lateral/` (2D + BCs) | groundwater-only |
| layered / heterogeneous soils | `benchmarks/b3-kirkland/` | `soil.types` list + 3D `soil.map.file` |
| coupled surface–subsurface catchment | `benchmarks/b5-vcatchment/` | the fully documented coupled reference |
| surface scalar transport | `validation/transport-surface-tracer-advection/` | upwind vs superbee pair |
| subsurface advection–dispersion | `validation/transport-ogata-banks-column/` | |
| variable-density / saltwater intrusion | `validation/transport-henry-saltwater-intrusion/` (simplest), `benchmarks/b6-kuan/` (full coupled tidal) | |

Every case directory has a README stating what it does and its expected
result — read it before adapting.

## 4. Run modes and required sections

At least one flow module must be on. `boundary_conditions` is always
optional (none = closed basin).

| Mode | `modules` | Additionally required sections |
|---|---|---|
| surface only | `surface_water: true` | `surface_water`, `initial_conditions.surface.eta` |
| groundwater only | `groundwater: true` | `groundwater`, `soil`, `initial_conditions.groundwater` (exactly one of `water_table`/`head`/`moisture`) |
| coupled | both true | all of the above; `coupling.mode` optional (`sync` default) |
| + transport | `transport: true` | `transport` section + `initial_conditions.transport.{surface,groundwater}` for each active grid |

## 5. Input data files

**Rasters** (bed, roughness, initial fields): either a flat
whitespace-separated list (`#` comments allowed) of exactly `nx·ny`
values in `j·nx + i` order, or a headered raster (`ncols`/`nrows`
required; `cellsize`, `nodata_value` optional; nodata → NaN). 3D fields
(`soil.map.file`, 3D ICs) are flat lists of `nx·ny·nz` values in
`(j·nx + i)·nz + k` order. A wrong count aborts with
observed-vs-expected — trust that message.

**Time series** (rain, evaporation, wind, time-varying BCs): two columns
`t value`, seconds + SI, strictly increasing time, `#` comments,
piecewise-linear, **clamped** to the end values outside the sampled
range. Step changes are expressed with close sample pairs
(`12000 5.5e-6` then `12001 0`).

Generate these files with small Python snippets in the case directory
(the worked example shows the pattern); keep them under `input/`.

## 6. Pitfalls (each of these has bitten someone)

1. `surface_water.min_depth` and `wetting_face_depth` are **required**
   (no defaults). Typical: `1.0e-6` both; use a larger
   `wetting_face_depth` (~1e-3) when sub-mm rain films should infiltrate
   but not advect momentum.
2. **`kind: discharge` takes TOTAL volumetric discharge [m³/s]** over
   the region, not unit discharge [m²/s]. Converting literature unit
   discharge: multiply by the inflow width. Getting this wrong over-feeds
   by the width factor and blows up.
3. **`kind: outflow` needs a continuing bed slope** at the boundary; on
   a *flat* channel end it deadlocks — prescribe a bed-level `kind: eta`
   stage there instead (the b5 outlet pattern).
4. Boundary polygons rasterize against **cell centers**; a point exactly
   on an edge counts as inside; a polygon selecting zero cells is fatal.
   Legacy-style index rectangles become 4-vertex polygons that
   comfortably overlap the intended centers.
5. **Coupled runs**: any configured `groundwater_top` condition must be
   `kind: flux` (the coupler owns the top head); `head`-kind top
   conditions are rejected. Rain falls on wet *and* dry cells.
6. **`coupling.mode: sync`** (default): `time.dt` is only the *initial*
   step — the run marches on a common adaptive step bounded by
   `groundwater.timestep.dt_min/dt_max`. Bound the step there, not via
   `dt`. **`subcycled`**: the surface marches at the *fixed* `time.dt` —
   the user must choose it to satisfy the surface CFL
   (`dt ≲ dx / √(g·h_max)` as a first estimate); too large degrades
   accuracy in flashy runoff.
7. Surface-only runs use fixed `time.dt`; watch the log's CFL warnings
   on the smoke run and shrink `dt` if they persist.
8. `groundwater.timestep`: `dt_init`, `dt_min`, `dt_max` all required;
   the adaptive controller lives within `[dt_min, dt_max]`. For
   fixed-step studies set all three equal.
9. **Transport constraints** (load-bearing; details in
   `validation/README.md`): prescribed subsurface boundary
   concentration enters on the **y+ (north) side only** — orient the
   geometry so scalar inflow comes from y+; the density law
   `r_ρ = 1 + 7.44e-4·s` is hardwired, so a target density contrast is
   set via the salinity surrogate `s = (Δρ/ρ)/7.44e-4`; configured
   dispersivities act on volumetric face fluxes, not Darcy velocities
   (not directly the physical dispersivity).
10. Surface `scalar_value` conditions: on a `discharge`-covered region
    they set the inflow concentration; elsewhere they hold the region's
    wet cells at the value (tide/stage salinity).
11. Shallow steady flows inside `friction.thin_layer_depth` (default
    0.1 m) sit in the legacy drag-regularization band — cm-scale normal
    depths can offset ~20 % from Manning theory; set `thin_layer_depth`
    below the expected depths when that matters.
12. Model scope limits — do not try to configure around them: no
    periodic boundaries, no Darcy–Weisbach friction, no
    surface-infiltration source term other than the groundwater coupling
    itself, one scalar, groundwater scheme is PCA only, no adaptive
    *surface* stepping outside sync-coupled mode.
13. `output.variables` lists must match enabled modules (`seepage` exists
    only in coupled runs; `concentration` needs groundwater+transport,
    `concentration_surface` needs surface+transport). Monitors record
    every time step at one cell — use them for hydrographs instead of
    frequent field output.
14. Restart: `output.checkpoint.interval` writes `/checkpoint/<t>`
    groups (one always at `t_end`); resuming needs the same physical
    configuration and reproduces the uninterrupted run bitwise.

## 7. Minimal skeleton (surface-water rain-runoff)

For orientation only — copy a real example instead (§3):

```yaml
simulation: {id: my-case}
domain:
  nx: 40
  ny: 5
  nz: 1
  dx: 1.0
  dy: 1.0
  dz: 1.0
  bottom_elevation: {file: input/dem.dat}
time: {dt: 0.2, t_end: 1200, output_interval: 60}
modules: {surface_water: true, groundwater: false, transport: false}
surface_water:
  friction: {law: manning, coefficient: {constant: 0.03}, thin_layer_depth: 0.1}
  min_depth: 1.0e-6
  wetting_face_depth: 1.0e-6
  rainfall: {series: {file: input/rain.dat}}
initial_conditions:
  surface: {eta: {constant: -1.0}}        # dry start: below the bed
boundary_conditions:
  - name: outlet
    region: {polygon: [[39.5, -1.0], [41.0, -1.0], [41.0, 6.0], [39.5, 6.0]]}
    target: surface
    kind: outflow
output:
  filename: out/output.h5
  variables: {surface: [eta, depth, uu, vv]}
  monitors: [{name: outlet, i: 39, j: 2, variables: [depth, uu]}]
```

Everything not shown has a documented default or is mode-dependent —
`parameters.md` is the complete truth.
