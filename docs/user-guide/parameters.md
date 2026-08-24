# Frehg2 YAML configuration reference (v2 schema)

Every recognized key, its type, whether it is required, and its default.
This table is checked against `src/core/ConfigSchema.cpp` by
`scripts/check_parameter_docs.py` (run in CI): a key present in the schema
but missing here — or vice versa — fails the build.

Conventions: units are SI (meters, seconds); paths are relative to the
configuration file; "one of" groups require exactly one member.

## simulation

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `simulation.id` | string | yes | — | short run identifier |
| `simulation.title` | string | no | "" | free-text description |

## domain

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `domain.nx` | int ≥ 1 | yes | — | global cells in i |
| `domain.ny` | int ≥ 1 | yes | — | global cells in j |
| `domain.nz` | int ≥ 1 | yes | — | subsurface layers |
| `domain.dx` | real > 0 | yes | — | cell size in i [m] |
| `domain.dy` | real > 0 | yes | — | cell size in j [m] |
| `domain.dz` | real > 0 | yes | — | top-layer thickness [m] |
| `domain.dz_stretch` | real > 0 | no | 1.0 | geometric layer growth (legacy `dz_incre`) |
| `domain.bottom_elevation.constant` | real | one of | — | uniform bed elevation [m] |
| `domain.bottom_elevation.file` | path | one of | — | bed elevation raster |
| `domain.follow_terrain` | bool | no | false | terrain-following subsurface mesh |
| `domain.terrain_layers` | `scaled` \| `uniform` | no | scaled | terrain-mesh vertical rule (amendment A12): `scaled` spans each column's nz layers from the local bed to the common box bottom (the legacy wedge — b6's tank); `uniform` stacks the configured `dz`·`dz_stretch`^k profile below each column's own bed (the terrain-parallel slab of the b5 reference). Requires `follow_terrain: true` |
| `domain.decomposition.mpi_nx` | int ≥ 1 or `auto` | no | auto | ranks along i |
| `domain.decomposition.mpi_ny` | int ≥ 1 or `auto` | no | auto | ranks along j |

## time

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `time.dt` | real > 0 | yes | — | surface step [s]; fixed, except in sync-coupled runs where it is the initial value of the common adaptive step (amendment A10) |
| `time.t_start` | real ≥ 0, whole seconds | no | 0 | start time [s] |
| `time.t_end` | real > 0 | yes | — | end time [s]; must exceed `t_start` |
| `time.output_interval` | real ≥ 1, whole seconds | yes | — | spatial output cadence [s] (§7 keys times by integer seconds) |

## modules

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `modules.surface_water` | bool | no | false | run the SWE module |
| `modules.groundwater` | bool | no | false | run the Richards module |
| `modules.transport` | bool | no | false | run scalar transport (needs a flow module) |

At least one of `surface_water`/`groundwater` must be true.

## surface_water (required section when the module is on)

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `surface_water.gravity` | real > 0 | no | 9.81 | [m/s²] |
| `surface_water.friction.law` | `manning` \| `chezy` | no | manning | drag law (`chezy` is the plan §5.8 addition) |
| `surface_water.friction.coefficient.constant` | real | one of | — | Manning n or Chezy C |
| `surface_water.friction.coefficient.file` | path | one of | — | roughness raster |
| `surface_water.friction.thin_layer_depth` | real > 0 | no | 0.1 | legacy `hD` exponent-switch depth [m] |
| `surface_water.viscosity.x` | real ≥ 0 | no | 1e-6 | eddy viscosity in i [m²/s] |
| `surface_water.viscosity.y` | real ≥ 0 | no | 1e-6 | eddy viscosity in j [m²/s] |
| `surface_water.min_depth` | real > 0 | yes | — | wet/dry threshold [m] (legacy `min_dept`) |
| `surface_water.wetting_face_depth` | real > 0 | yes | — | face-wetting threshold [m] (legacy `wtfh`) |
| `surface_water.wind.enabled` | bool | no | false | quadratic wind stress |
| `surface_water.wind.cd` | real ≥ 0 | no | 0.0013 | drag coefficient (legacy `Cw`) |
| `surface_water.wind.attenuation_depth` | real > 0 | no | 5.0 | thin-layer attenuation depth (legacy `CwT`) |
| `surface_water.wind.north_angle` | real | no | 0 | grid-to-north rotation [deg] |
| `surface_water.wind.speed.constant` | real | one of | — | wind speed [m/s] |
| `surface_water.wind.speed.series.file` | path | one of | — | wind speed series |
| `surface_water.wind.direction.constant` | real | one of | — | wind direction [deg] |
| `surface_water.wind.direction.series.file` | path | one of | — | wind direction series |
| `surface_water.rainfall.constant` | real | one of | 0 | rain rate [m/s] |
| `surface_water.rainfall.series.file` | path | one of | — | rain series |
| `surface_water.rainfall.exclude.polygon` | [x, y] list ≥ 3 | no | — | region receiving no rainfall. Expresses the legacy hardcoded skip of the last global row (`shallowwater.c:596`), which the b1 goldens embed for their outlet row; omit it for uniform rain |
| `surface_water.evaporation.constant` | real | one of | 0 | evaporation rate [m/s] |
| `surface_water.evaporation.series.file` | path | one of | — | evaporation series |

## groundwater (required section when the module is on)

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `groundwater.scheme` | `pca` | no | pca | sole allowed value (documents the design, plan §6) |
| `groundwater.use_full3d` | bool | no | true | false zeroes lateral K in unsaturated cells |
| `groundwater.timestep.dt_init` | real > 0 | yes | — | initial dtg [s] |
| `groundwater.timestep.dt_min` | real > 0 | yes | — | dtg lower clamp [s]; `dt_min ≤ dt_init ≤ dt_max` |
| `groundwater.timestep.dt_max` | real > 0 | yes | — | dtg upper clamp [s] |
| `groundwater.timestep.dq_grow` | real > 0 | no | 0.01 | flux-change growth threshold |
| `groundwater.timestep.dq_shrink` | real > 0 | no | 0.02 | flux-change shrink threshold |
| `groundwater.timestep.courant_max` | real > 0 | no | 2.0 | legacy `Co_max` |
| `groundwater.specific_storage` | real ≥ 0 | yes | — | Ss [1/m] |
| `groundwater.reallocation_surplus` | `drop` \| `redistribute` | no | `drop` | post-allocation surplus of saturation-adjacent cells: discard (legacy sweep, groundwater.c:1026 disabled the transfer; the b2 golden embeds it) or redistribute vertically into pore room (mass-conserving; b3 uses it) |
| `groundwater.density_coupling.enabled` | bool | no | false | baroclinic feedback (requires transport) |

## soil (required section when groundwater is on)

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `soil.types[].name` | string | yes | — | unique type name |
| `soil.types[].ksx` | real ≥ 0 | yes | — | saturated K in i [m/s] |
| `soil.types[].ksy` | real ≥ 0 | yes | — | saturated K in j [m/s] |
| `soil.types[].ksz` | real ≥ 0 | yes | — | saturated K in k [m/s] |
| `soil.types[].theta_s` | real ∈ (0, 1] | yes | — | saturated water content |
| `soil.types[].theta_r` | real ∈ [0, 1] | yes | — | residual water content |
| `soil.types[].vg_alpha` | real > 0 | yes | — | van Genuchten α [1/m] |
| `soil.types[].vg_n` | real > 1 | yes | — | van Genuchten n |
| `soil.types[].aev` | real ≤ 0 | no | 0 | saturation-cutoff head [m] (live outside `use_mvg` in legacy; plan §2.1) |
| `soil.map.constant` | string | one of | — | uniform type by name |
| `soil.map.file` | path | one of | — | 3D id field, flat `(j*nx+i)*nz+k`, ids index `soil.types` |

## coupling

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `coupling.mode` | `sync` \| `subcycled` | no | sync | surface–subsurface synchronization (plan §5.7). `sync`: both modules march in lockstep on the common adaptive step (initial `time.dt`, clamped to `groundwater.timestep` bounds — legacy solve.c:193, amendment A10); `subcycled`: fixed surface `time.dt` with the adaptive dtg subcycled inside it, exchange volumes accumulated conservatively (amendment A9) |

## initial_conditions

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `initial_conditions.surface.eta.{constant,file}` | one of | when SWE on | — | initial stage [m] |
| `initial_conditions.surface.uu.{constant,file}` | one of | no | — | initial u velocity |
| `initial_conditions.surface.vv.{constant,file}` | one of | no | — | initial v velocity |
| `initial_conditions.groundwater.water_table.{constant,file}` | one of three forms | when GW on | — | water-table elevation [m] |
| `initial_conditions.groundwater.head.{constant,file}` | one of three forms | when GW on | — | pressure head [m] |
| `initial_conditions.groundwater.moisture.{constant,file}` | one of three forms | when GW on | — | water content |
| `initial_conditions.transport.surface.{constant,file}` | one of | when transport+SWE | — | surface concentration |
| `initial_conditions.transport.groundwater.{constant,file}` | one of | when transport+GW | — | subsurface concentration |

## boundary_conditions (list; plan §5.6)

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `boundary_conditions[].name` | string | yes | — | unique name |
| `boundary_conditions[].region.polygon` | list of [x, y], ≥ 3 | yes | — | region vertices [m]; on-edge counts inside |
| `boundary_conditions[].target` | `surface` \| `groundwater_top` \| `groundwater_bottom` \| `groundwater_side` | yes | — | where the condition applies |
| `boundary_conditions[].kind` | `eta` \| `discharge` \| `velocity` \| `outflow` \| `head` \| `flux` \| `scalar_value` | yes | — | eta/discharge/velocity/outflow: surface only; head/flux: groundwater only; scalar_value needs transport. `outflow` is free (transmissive) outflow at domain-edge faces — the stage is extrapolated down the continued bed slope — and takes no `value`. `scalar_value` semantics (P4): target `surface` — members covered by a `discharge` condition inject that inflow's concentration (legacy `s_inflow`); every other member is a wet-cell Dirichlet (the legacy tide salinity, generalized to any region — dry cells stay at 0); target `groundwater_side` — the ghost concentration along the member edge (legacy `s_yp`/`s_ym`), feeding upwinded inflow, the limiter bound, and the baroclinic boundary face |
| `boundary_conditions[].value.constant` | real | one of | — | fixed value |
| `boundary_conditions[].value.series.file` | path | one of | — | time series |
| `boundary_conditions[].value.gravity` | `true` | one of | — | free drainage (kind flux; legacy `bctype_GW` code 3) |
| `boundary_conditions[].value.hydrostatic.eta` | real | one of | — | head = eta − z (kind head). In coupled runs the hydrostatic ghost follows the *live local surface* instead — the edge column's bed plus its current depth, density-scaled once transport+density coupling are active (the legacy enforce_head_bc rule; amendment A18) — and the configured eta applies only without a surface module |

## transport (required section when the module is on)

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `transport.scheme.advection` | `upwind` \| `superbee` | no | upwind | advection scheme |
| `transport.surface_diffusivity.x` | real ≥ 0 | no | 1e-10 | [m²/s] |
| `transport.surface_diffusivity.y` | real ≥ 0 | no | 1e-10 | [m²/s] |
| `transport.dispersion.longitudinal` | real ≥ 0 | no | 0 | [m] |
| `transport.dispersion.transverse` | real ≥ 0 | no | 0 | [m] |
| `transport.dispersion.molecular` | real ≥ 0 | no | 1e-10 | [m²/s] |
| `transport.bounds.min` | real | no | 0 | scalar lower bound (replaces the legacy [0, 200] clamp; plan §3.2) |
| `transport.bounds.max` | real | no | unbounded | scalar upper bound |

## output

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `output.filename` | string | yes | — | single HDF5 file (plan §7) |
| `output.variables.surface` | list ⊆ {eta, depth, uu, vv, seepage} | no | — | requires the SWE module |
| `output.variables.groundwater` | list ⊆ {hydraulic_head, water_content, qx, qy, qz} | no | — | requires the GW module |
| `output.variables.transport` | list ⊆ {concentration, concentration_surface} | no | — | requires transport; `concentration` is the subsurface scalar (with the GW module), `concentration_surface` the surface scalar (with the SWE module) |
| `output.monitors[].name` | string | yes | — | unique monitor name |
| `output.monitors[].type` | `point` | no | point | monitor kind |
| `output.monitors[].i` | int ≥ 0 | yes | — | global cell index in i (< nx) |
| `output.monitors[].j` | int ≥ 0 | yes | — | global cell index in j (< ny) |
| `output.monitors[].variables` | list, ≥ 1 | yes | — | variables of enabled modules |
| `output.checkpoint.interval` | real ≥ 0, whole seconds | no | 0 | 0 = off; final checkpoint always written |

## restart

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `restart.enabled` | bool | no | false | resume from a checkpoint |
| `restart.file` | string | when enabled | "" | checkpointed output file (must exist) |
| `restart.time` | real | no | 0 | checkpoint time to resume from [s] |

## solver

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `solver.petsc_options_file` | string | no | "" | PETSc options file (overrides the `fs_`/`gw_` defaults, plan §5.1) |

## runtime

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `runtime.gpu_aware_mpi` | `auto` \| `on` \| `off` | no | auto | pass device buffers to MPI or stage through host mirrors (plan §5.3) |
