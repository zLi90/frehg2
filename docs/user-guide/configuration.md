# The configuration file

A run is driven entirely by one YAML file. It is validated against a **strict
v2 schema**: unknown keys are rejected (with spelling suggestions), types and
ranges are enforced, cross-field rules are checked, and referenced input
files must exist. Validate early and often:

```bash
build/src/frehg --validate my-case.yaml
```

Conventions: **units are SI** (meters, seconds, m/s); every relative path is
resolved against the directory of the configuration file; "one of" groups
require exactly one member. The complete, machine-checked key reference is
the [parameter reference](parameters.md) — the sections below walk through
each configuration block.

## Which sections a case needs

Every run needs `simulation`, `domain`, `time`, `modules`, and `output`.
The rest depends on which modules are on:

| Run mode | Additionally required sections |
|---|---|
| surface water only | `surface_water`, `initial_conditions.surface` |
| groundwater only | `groundwater`, `soil`, `initial_conditions.groundwater` |
| coupled | all of the above; `coupling` is optional (defaults to `sync`) |
| + transport | `transport`, `initial_conditions.transport.{surface,groundwater}` for each active grid |

`boundary_conditions` is always optional (a domain with none is a closed
basin — for groundwater, a sealed box).

## `simulation`, `domain`, `time`

```yaml
simulation:
  id: my-case                # short run identifier (required)
  title: My overland case    # free text (optional)

domain:
  nx: 200                    # global cells in i (x)   (required)
  ny: 10                     # global cells in j (y)   (required)
  nz: 1                      # subsurface layers; use 1 for SWE-only (required)
  dx: 0.5                    # cell size in i [m]      (required)
  dy: 0.5                    # cell size in j [m]      (required)
  dz: 1.0                    # top-layer thickness [m]; unused by SWE but required
  # dz_stretch: 1.0          # geometric layer growth factor (default 1.0)
  bottom_elevation:          # bed elevation — exactly one of:
    constant: 0.0            #   uniform bed [m]
    # file: input/dem.dat    #   or a raster (see the input-data page)
  follow_terrain: false      # terrain-following subsurface mesh (SWE: leave false)
  # terrain_layers: scaled   # terrain-mesh vertical rule: scaled (legacy wedge)
  #                          # or uniform (dz profile below each bed); needs
  #                          # follow_terrain: true
  decomposition:             # optional; default auto
    mpi_nx: auto
    mpi_ny: auto

time:
  dt: 0.1                    # surface time step [s] (required). Fixed — except in
                             # sync-coupled runs, where it is only the INITIAL value
                             # of the common adaptive step (see `coupling` below)
  t_start: 0                 # start time [s]               (default 0)
  t_end: 4800                # end time [s]                 (required, > t_start)
  output_interval: 120       # spatial-output cadence [s]   (required, whole seconds)
```

The grid is a regular Cartesian mesh; cell `(i, j)` has its center at
`x = (i + 0.5)·dx`, `y = (j + 0.5)·dy`. `i` runs west→east, `j` runs
south→north. `output_interval` and `t_end`/`t_start` are in whole seconds
because output times are keyed by integer seconds in the HDF5 file.

## `modules`

```yaml
modules:
  surface_water: true
  groundwater: true          # both on = coupled run
  transport: true            # scalar transport (needs a flow module)
```

At least one flow module must be on; enabling both runs the coupled
surface–subsurface model. Transport rides whichever flow modules are on
and requires a `transport:` section plus scalar initial conditions. See
`benchmarks/b2-gw/b2-gw.yaml` for a minimal groundwater-only case,
`benchmarks/b5-vcatchment/b5-vcatchment.yaml` for a complete coupled case,
and `benchmarks/b6-kuan/b6-kuan-td.yaml` for coupled flow with salt
transport and density feedback.

## `surface_water`

```yaml
surface_water:
  gravity: 9.81                    # [m/s²] (default 9.81)
  friction:
    law: manning                   # manning | chezy   (default manning)
    coefficient:                   # exactly one of:
      constant: 0.03               #   Manning n (or Chezy C if law: chezy)
      # file: input/roughness.dat  #   or a roughness raster
    thin_layer_depth: 0.1          # legacy hD exponent-switch depth [m]
  viscosity: {x: 1.0e-6, y: 1.0e-6}  # eddy viscosity [m²/s]
  min_depth: 1.0e-6                # wet/dry threshold [m]        (REQUIRED)
  wetting_face_depth: 1.0e-6       # face-wetting threshold [m]   (REQUIRED)

  rainfall:                        # optional; exactly one of constant/series
    constant: 5.0e-6               #   rate [m/s]
    # series: {file: input/rain.dat}
    # exclude: {polygon: [[x0,y0],[x1,y1],...]}  # region receiving no rain
  evaporation: {constant: 0.0}     # optional; [m/s]

  wind:                            # optional; off by default
    enabled: false
    cd: 0.0013                     # drag coefficient
    attenuation_depth: 5.0         # thin-layer attenuation depth [m]
    north_angle: 0                 # grid-to-north rotation [deg]
    speed: {constant: 0.0}         # or series: {file: ...}
    direction: {constant: 0.0}     # or series: {file: ...}
```

`min_depth` and `wetting_face_depth` are the only required surface-water
keys; everything else has a legacy-faithful default. `friction.law: chezy`
selects the drag law added by the upgrade plan (§5.8); with `chezy` the
coefficient is the Chézy `C`, otherwise it is Manning's `n`.

`rainfall.exclude.polygon` expresses the legacy hard-coded "no rain on the
last grid row" and is only needed to reproduce that specific legacy
convention (the b1 benchmark uses it). For a new case with uniform rainfall,
omit it. In coupled runs rain falls on every non-excluded cell, wet or dry
(dry cells start infiltrating immediately).

## `groundwater` and `soil`

Required whenever `modules.groundwater` is on (alone or coupled):

```yaml
groundwater:
  scheme: pca                      # sole allowed value
  use_full3d: true                 # false zeroes lateral K in unsaturated cells
  timestep:                        # the adaptive subsurface step dtg [s]
    dt_init: 0.1                   # initial dtg          (required)
    dt_min: 0.01                   # lower clamp          (required)
    dt_max: 10.0                   # upper clamp          (required)
    # dq_grow: 0.01                # flux-change growth threshold
    # dq_shrink: 0.02              # flux-change shrink threshold
    # courant_max: 2.0             # legacy Co_max
  specific_storage: 1.0e-5         # Ss [1/m]              (required)
  # reallocation_surplus: drop     # drop | redistribute (see parameters.md)

soil:
  types:                           # one or more van Genuchten soil types
    - name: loam
      ksx: 2.89e-6                 # saturated K [m/s] in i / j / k
      ksy: 2.89e-6
      ksz: 2.89e-6
      theta_s: 0.33                # saturated water content
      theta_r: 0.0                 # residual water content
      vg_alpha: 1.43               # van Genuchten alpha [1/m]
      vg_n: 1.56                   # van Genuchten n
      # aev: 0.0                   # saturation-cutoff head [m], <= 0
  map: {constant: loam}            # uniform type by name, or file: (3D id field)
```

The initial state comes from exactly one of three forms under
`initial_conditions.groundwater`: `water_table` (elevation of a hydrostatic
water table [m]), `head` (pressure head [m]), or `moisture` (water
content), each `constant` or `file`.

The subsurface mesh is the box `nz` layers deep below the domain: layer
thicknesses follow `dz`·`dz_stretch`^k; with `follow_terrain: true` the
mesh follows the bed instead of a flat top (see `domain.terrain_layers`
for the two vertical rules). Groundwater boundary conditions target
`groundwater_top`, `groundwater_bottom`, or `groundwater_side` regions
with kinds `head` or `flux` — including free drainage
(`value: {gravity: true}`) and hydrostatic side heads
(`value: {hydrostatic: {eta: ...}}`); see the `boundary_conditions`
section below.

## `coupling` (coupled runs)

```yaml
coupling:
  mode: sync                       # sync | subcycled   (default sync)
```

With both flow modules on, the coupler exchanges water each step exactly
where the legacy model did: wet surface cells drive the subsurface as a
ponded-head boundary (limited to the water actually available — thin rain
films cannot over-infiltrate), dry cells act as seepage faces, and
subsurface seepage returns to the surface as a depth change.

- **`sync`** (the legacy default): both modules march in lockstep on the
  *common adaptive step*. `time.dt` is only its initial value; the step
  then adapts within the `groundwater.timestep` bounds. Robust for every
  regime.
- **`subcycled`**: the surface marches at the *fixed* `time.dt` and the
  subsurface subcycles inside it with its adaptive dtg; exchange volumes
  are accumulated conservatively across the window. Choose `time.dt` to
  satisfy the surface CFL yourself — a too-large fixed step degrades
  accuracy in flashy overland-flow regimes (the b5 config documents a
  measured example).

Two coupled-run rules the validator enforces: configured `groundwater_top`
conditions must be `kind: flux` (the coupler owns the top head), and every
column must remain active in the surface system. The applied seepage rate
is available as the `seepage` output variable and the exchanged volume as
the `seepage` mass-audit column (see the [output reference](output.md)).

## `transport` (scalar transport)

```yaml
modules: {surface_water: true, groundwater: true, transport: true}
transport:
  scheme: {advection: superbee}    # upwind | superbee   (default upwind)
  surface_diffusivity: {x: 0.0, y: 1.0e-10}      # [m²/s]
  dispersion:                      # subsurface tensor
    longitudinal: 0.002            # [m]
    transverse: 4.0e-04            # [m]
    molecular: 1.0e-10             # [m²/s]
  # bounds: {min: 0, max: 35}      # optional clamp (default [0, unbounded))
initial_conditions:
  transport:
    surface: {constant: 35.0}      # with the SWE module
    groundwater: {constant: 0.0}   # with the GW module
```

One scalar (e.g. salinity) advances once per step after the flow update:
explicit finite-volume advection (first-order upwind or TVD superbee),
constant diffusion on the surface, the full anisotropic dispersion tensor
in the subsurface, and — in coupled runs — exchange through the seepage
(infiltrating water carries the surface concentration down, discharging
water carries the top-cell concentration up). Scalar boundary values use
`kind: scalar_value` conditions (below): on the surface they prescribe the
concentration of a region's wet cells (a stage salinity — pair the polygon
with your `eta` condition, or cover a whole tidal flat) or, when the
polygon lies on a `discharge` condition, the inflow concentration; on
`groundwater_side` they set the boundary concentration the inflow carries
and the density boundary sees.

Setting `groundwater.density_coupling.enabled: true` feeds the subsurface
scalar back into the Darcy fluxes as the legacy density/viscosity ratios
(r_rho = 1 + 0.000744 s, r_visc = 1/(1 + 0.0022 s)) — the baroclinic
pathway saltwater-intrusion problems need (`benchmarks/b6-kuan`). Outputs:
`concentration` (subsurface field), `concentration_surface` (surface
field), and the `/monitor/transport_audit` scalar-budget table
([output reference](output.md)).

## `initial_conditions`

```yaml
initial_conditions:
  surface:
    eta:  {constant: -10.0}   # initial free-surface elevation [m] (required for SWE)
    # eta:  {file: input/eta0.dat}
    uu:   {constant: 0.0}     # optional initial x-velocity
    vv:   {constant: 0.0}     # optional initial y-velocity
```

`eta` is the **absolute** free-surface elevation (same datum as
`bottom_elevation`). A cell is dry where `eta ≤ bottom`; a common way to
start a dry domain is to set `eta` well below the bed everywhere (e.g.
`-10.0` when the bed is near 0). Water depth is `eta − bottom`. Frehg2
applies the legacy internal elevation offset transparently — you always work
in absolute elevations, and outputs are in absolute elevations too.

When the groundwater module is on, `initial_conditions.groundwater` is
also required — exactly one of `water_table`, `head`, or `moisture`.

## `boundary_conditions`

An optional list. Each entry rasterizes a polygon region onto the grid and
applies a condition to a `target`: `surface` (SWE) or one of the
subsurface faces `groundwater_top`, `groundwater_bottom`,
`groundwater_side`:

```yaml
boundary_conditions:
  - name: outlet                    # unique name (required)
    region:
      polygon: [[21.8, -0.01], [22.0, -0.01], [22.0, 2.0], [21.8, 2.0]]
    target: surface
    kind: outflow                   # surface: eta | discharge | velocity | outflow
    # value: {constant: ...}        # required except for kind: outflow

  - name: drain                     # a groundwater example
    region: {polygon: [[-0.1, -0.1], [22.1, -0.1], [22.1, 2.1], [-0.1, 2.1]]}
    target: groundwater_bottom
    kind: flux
    value: {gravity: true}          # free drainage
```

Surface kinds:

| `kind` | Meaning | `value` |
|---|---|---|
| `eta` | prescribed free-surface elevation [m] | required (`constant` or `series`) |
| `discharge` | prescribed total discharge into the region [m³/s] | required |
| `velocity` | prescribed face-normal velocity [m/s] at domain-edge faces | required |
| `outflow` | free (transmissive) outflow — stage extrapolated down the continued bed slope at the domain-edge faces | **none** |
| `scalar_value` | concentration of the region's wet cells, or of a `discharge` inflow | required (needs transport) |

Groundwater kinds:

| `kind` | Meaning | `value` |
|---|---|---|
| `head` | prescribed pressure head [m] | `constant`, `series`, or `hydrostatic: {eta: ...}` (head = eta − z) |
| `flux` | prescribed Darcy flux [m/s] | `constant`, `series`, or `gravity: true` (free drainage) |
| `scalar_value` | boundary concentration at `groundwater_side` faces | required (needs transport) |

`outflow` is the right choice for a "let water leave here" open surface
boundary (it neither reflects nor holds back the flow) — but it needs a
continuing bed slope; on a *flat* channel end, prescribe a `kind: eta`
stage instead (the b5 outlet does this). A time-varying value uses
`value: {series: {file: bc.dat}}` instead of `constant`. The polygon must
have at least three vertices; a point exactly on an edge counts as inside.
In coupled runs `groundwater_top` conditions must be `kind: flux` (the
coupler owns the top head), and the hydrostatic side ghost follows the
*live* surface — the configured `hydrostatic: {eta: ...}` value only
matters in uncoupled runs.

## `output`

```yaml
output:
  filename: out/output.h5           # single HDF5 file (required)
  variables:
    surface: [eta, depth, uu, vv]   # subset of {eta, depth, uu, vv, seepage}
    groundwater: [hydraulic_head, water_content]
                                    # subset of {hydraulic_head, water_content, qx, qy, qz}
    transport: [concentration, concentration_surface]
  monitors:                         # optional point probes
    - name: outlet_q
      i: 199                        # global cell index in i (0 ≤ i < nx)
      j: 4                          # global cell index in j (0 ≤ j < ny)
      variables: [depth, uu]
  checkpoint:
    interval: 0                     # [s]; 0 = off. A final checkpoint is always written.
```

`variables.surface` / `variables.groundwater` / `variables.transport`
select which fields are written at each `output_interval` (each list
requires its module; `seepage` — the surface-applied exchange rate [m/s] —
exists in coupled runs only). Monitors record their listed variables
**every time step** at one cell, giving a high-resolution hydrograph.
Mass-balance tables are written automatically (see the
[output reference](output.md)).

## `restart` and `solver` (optional)

```yaml
restart:
  enabled: true
  file: out/output.h5   # a file containing a checkpoint at `time`
  time: 9000            # checkpoint time to resume from [s]

solver:
  # Per-system linear-solver selection (v2 plan §2.2). bjacobi-icc is the
  # v1 default; amg = hypre BoomerAMG (needs a hypre-enabled PETSc, the
  # scripts/ci_install_deps.sh build); gamg = PETSc smoothed aggregation.
  # AMG-class choices reuse their hierarchy across solves and rebuild on
  # the reuse_* policy (and, for the surface system, whenever the wet/dry
  # mask changes).
  surface:
    preconditioner: bjacobi-icc   # bjacobi-icc | amg | gamg
    rtol: 1.0e-8
    atol: 1.0e-14
    max_iterations: 500
    reuse_max_solves: 50          # amg/gamg: rebuild at least every N solves
    reuse_iteration_factor: 1.5   # amg/gamg: early rebuild trigger
  groundwater:
    preconditioner: bjacobi-icc
    max_iterations: 1000
  petsc_options_file: ""  # optional PETSc options file overriding the fs_/gw_ defaults
```

Resuming from a checkpoint reproduces the uninterrupted run bitwise — see
[Checkpoint and restart](restart.md).
