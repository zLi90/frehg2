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
| `domain.periodic.x` / `domain.periodic.y` | bool | no | false | Stage-3 T3: periodic lateral topology per axis — a **domain** property (every module wraps: surface momentum/continuity, the groundwater matrix and fluxes, transport, baroclinic faces). Needs ≥ 3 cells on the axis; a BC whose polygon selects only faces on a periodic axis is a setup error, and non-periodic inputs across the seam warn (run-record `warnings`) |
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
| `modules.transport` | bool | no | false | run scalar (salinity) transport (needs a flow module) |
| `modules.temperature` | bool | no | false | run temperature transport (v2 Q5; needs a flow module) |

At least one of `surface_water`/`groundwater` must be true.


## evapotranspiration (Stage-3 T5)

Replaces BOTH legacy evaporation keys (`surface_water.evaporation`,
`groundwater.evaporation`) — the schema rejects mixing the forms. See
`docs/theory/evapotranspiration.md` for the channel semantics, the
time-basis complementarity, the supply cap and the salt rules.

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `evapotranspiration.potential` | `bulk` \| `penman_monteith` \| `prescribed` | yes | — | the shared potential mode. `bulk` needs the atmosphere block incl. `surface_temperature`; `penman_monteith` (Stage-3 T7) needs the atmosphere block incl. `shortwave_incoming` + `longwave_in` and is rejected with `modules.temperature`; `prescribed` needs `prescribed` |
| `evapotranspiration.prescribed` | series/constant [m/s] | with `prescribed` | — | the prescribed potential rate |
| `evapotranspiration.open_water` | bool | no | true | the surface channel (needs the SWE module): wet cells lose min(E_w·dt, depth) per step |
| `evapotranspiration.soil` | map | no | absent | presence enables the soil channel (needs the GW module) |
| `evapotranspiration.soil.region` | polygon | no | all columns | the soil-channel membership |
| `evapotranspiration.soil.limiter` | `alpha_barton` \| `beta_lee_pielke` \| `none` | no | alpha_barton | soil-rate limiter; `alpha_barton` acts inside the bulk humidity difference (bulk only — rejected with PM/prescribed) |
| `evapotranspiration.energy.albedo_water` | real [0, 1] | no | 0.06 | open-water albedo (penman_monteith, Stage-3 T7) |
| `evapotranspiration.energy.ground_heat_fraction.water` | real [0, 1] | no | 0.0 | G = fraction·R_n, open-water channel (penman_monteith) |
| `evapotranspiration.energy.ground_heat_fraction.soil` | real [0, 1] | no | 0.1 | G fraction, soil channel (penman_monteith) |
| `evapotranspiration.energy.ground_heat_fraction.canopy` | real [0, 1] | no | 0.1 | G fraction, canopy channel (penman_monteith) |

## vegetation (Stage-3 T6)

Static per-cell vegetation (needs the GW module and the ET block; with
salinity also `transport.subsurface_update: consistent`). See
`docs/theory/evapotranspiration.md` for the partition, Feddes stress and
the salt-exclusion rule.

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `vegetation.types[].name` | string | yes | — | unique type name |
| `vegetation.types[].lai` | real ≥ 0 | yes | — | leaf area index |
| `vegetation.types[].extinction` | real > 0 | yes | — | Beer–Lambert k; τ = exp(−k·LAI) scales the ground ET rates |
| `vegetation.types[].leaf_resistance` | real ≥ 0 | yes | — | r_l [s/m]; r_c = r_l/(0.5·LAI) |
| `vegetation.types[].root_depth` | real > 0 | yes | — | L [m] below the column's top surface (shallower columns renormalise, with a warning) |
| `vegetation.types[].root_profile` | `uniform` \| `linear` | no | uniform | root-density shape (linear = triangular, top-heavy) |
| `vegetation.types[].feddes.h1, h2` | real [m] | no | absent | anaerobic cutoff (together or not at all); OMIT for wetland species — α = 1 at any h ≥ h3, so transpiration continues while inundated (O7) |
| `vegetation.types[].feddes.h3, h4` | real [m] | yes | — | the dry ramp (α = 1 at h3, 0 at h4); h1 ≥ h2 ≥ h3 > h4 |
| `vegetation.types[].albedo` | real [0, 1] | with PM | — | canopy albedo (penman_monteith, Stage-3 T7: required on every type when the potential is PM) |
| `vegetation.map` | `{constant: <name>}` or `{file: <id raster>}` | yes | — | per-cell type; integer ids j·nx+i, **0 = bare**, k ≥ 1 = types[k−1] |
| `vegetation.lai.file` | path | no | absent | per-cell LAI override raster |

With salinity transport, ET requires `transport.surface_update: consistent`
(SWE) and `transport.subsurface_update: consistent` (GW); `transport.
surface_dryout` must stay `store` (the default since the T8 flip — an
explicit `discard` is rejected).

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
| `surface_water.wind.law` | `constant` \| `garratt` \| `smith-banke` \| `wu` \| `large-pond` | no | constant | v2 Q6 Cd(U₁₀) law: `constant` is the uncapped legacy `Cw`; the published laws — Garratt (1977) (0.75+0.067·U)e-3, Smith & Banke (1975) (0.63+0.066·U)e-3, Wu (1982) (0.8+0.065·U)e-3, Large & Pond (1981) piecewise — are capped at `cap` |
| `surface_water.wind.cd` | real ≥ 0 | no | 0.0013 | drag coefficient (legacy `Cw`; constant law only) |
| `surface_water.wind.cap` | real > 0 | no | 3.5e-3 | Cd cap for the U₁₀ laws (not the constant law) |
| `surface_water.wind.attenuation_depth` | real > 0 | no | 5.0 | thin-layer attenuation depth (legacy `CwT`) |
| `surface_water.wind.north_angle` | real | no | 0 | grid-to-north rotation [deg] (compass form only) |
| `surface_water.wind.speed.constant` | real | one of | — | wind speed [m/s] (compass form) |
| `surface_water.wind.speed.series.file` | path | one of | — | wind speed series |
| `surface_water.wind.direction.constant` | real | one of | — | wind direction [deg] |
| `surface_water.wind.direction.series.file` | path | one of | — | wind direction series; interpolated on the circle (unit-vector chord — 350°→10° passes through 0°, never 180°; sample finely enough that successive directions differ well below 180°) |
| `surface_water.wind.u10.constant` | real | one of | — | v2 Q6 grid-frame wind +x component [m/s]; the (u10, v10) pair replaces speed/direction (wrap-free) and ignores `north_angle` |
| `surface_water.wind.u10.series.file` | path | one of | — | u10 series |
| `surface_water.wind.v10.constant` | real | one of | — | wind +y component [m/s] |
| `surface_water.wind.v10.series.file` | path | one of | — | v10 series |
| `surface_water.rainfall.constant` | real | one of | 0 | rain rate [m/s] |
| `surface_water.rainfall.series.file` | path | one of | — | rain series |
| `surface_water.rainfall.exclude.polygon` | [x, y] list ≥ 3 | no | — | region receiving no rainfall. Expresses the legacy hardcoded skip of the last global row (`shallowwater.c:596`), which the b1 goldens embed for their outlet row; omit it for uniform rain |
| `surface_water.evaporation.constant` | real | one of | 0 | prescribed evaporation rate [m/s]; subtracted unconditionally over the non-excluded cells with the legacy dry clamp (b1 golden fidelity) |
| `surface_water.evaporation.series.file` | path | one of | — | prescribed evaporation series |
| `surface_water.evaporation.exclude.polygon` | [x, y] list ≥ 3 | no | — | region receiving no evaporation (prescribed mode; the rainfall-exclude symmetry, v2 Q4) |
| `surface_water.evaporation.mode` | `bulk` | no | — | v2 Q4 bulk-aerodynamic open-water evaporation: rate per step from the `atmosphere` block with water-surface q_g = q_sat(T_s), applied to wet cells only, at most the available depth per cell (no volume creation). Excludes `constant`/`series`/`exclude`; requires `atmosphere` |

## atmosphere (v2 Q4, plan §3.2)

Met forcing for the bulk-aerodynamic module (one vapor-physics module,
several consumers: open-water and soil evaporation now, the Q5 surface heat
exchange later). Every value is a `{constant: x}` or `{series: {file: p}}`
choice.

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `atmosphere.air_temperature` | series/constant | yes | — | [°C] |
| `atmosphere.surface_temperature` | series/constant | when a Q4 evap consumer is on | — | T_s [°C] for the evaporation consumers (surface bulk / soil bulk). The Q5 surface **heat** exchange never reads it — it evaluates the bulk chain at the local water temperature (V2-A17) |
| `atmosphere.pressure` | series/constant | yes | — | P₀ [kPa] |
| `atmosphere.specific_humidity` | series/constant | one of | — | q_a [−]; exactly one humidity form |
| `atmosphere.relative_humidity` | series/constant | one of | — | [−] of q_sat(air temperature) |
| `atmosphere.wind_speed` | series/constant | yes | — | U [m/s] for the bulk transfer functions. Configured independently of `surface_water.wind` (the vector momentum forcing) — Q5 decided against unifying them (a scalar-vs-vector convention would couple two independent schema surfaces); point both at the same series when physical consistency matters |
| `atmosphere.shortwave` | series/constant | no | 0 | absorbed shortwave radiation [W/m²] (Q5 bulk heat exchange; radiation *schemes* are out of scope, plan §10) |
| `atmosphere.longwave_in` | series/constant | no | 0 | incident longwave radiation [W/m²] (Q5 bulk heat exchange and the T7 Penman–Monteith potential; the emitted εσT⁴ term is computed) |
| `atmosphere.shortwave_incoming` | series/constant | no | 0 | **incident** shortwave R_s [W/m²] for the Penman–Monteith potential (Stage-3 T7) — each channel applies its own albedo; distinct from `shortwave`, which stays the Q5 *absorbed* SW |
| `atmosphere.wind_speed_floor` | real ≥ 0 | no | 0.5 | [m/s] still-air lower bound on the wind entering the bulk transfer functions (the GLM practice, plan §4.1) — applied to every bulk consumer at sampling |
| `atmosphere.aerodynamic_resistance.law` | `liu2006` \| `fao56` \| `log_profile` | no | liu2006 | r_a law of the ET consumers only (Stage-3 T5); the Q4/Q5 bulk consumers stay pinned to liu2006. `fao56` = 208/u2 with u2 = U·4.87/ln(67.8 z_ref − 5.42); `log_profile` = ln((z−d)/z0m)·ln((z−d)/z0h)/(κ²U) |
| `atmosphere.aerodynamic_resistance.reference_height` | real > 0 | no | 2.0 | z_ref [m] (fao56, log_profile) |
| `atmosphere.aerodynamic_resistance.z0m` | real > 0 | no | 0.01 | momentum roughness [m] (log_profile) |
| `atmosphere.aerodynamic_resistance.z0h` | real > 0 | no | 0.001 | scalar roughness [m] (log_profile) |
| `atmosphere.aerodynamic_resistance.displacement` | real ≥ 0 | no | 0.0 | displacement height d [m] (log_profile) |

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
| `groundwater.reallocation` | `legacy` \| `disabled` | no | `legacy` | Stage-3 T1: run the post-allocation step (`legacy`, the v2.0 sweep) or skip it (`disabled`: the corrector's flux-conservative θ stands, unsaturated heads re-derive from the retention curve, nothing is created or dropped — the `realloc`/`realloc_dropped` audit columns are exactly 0; the token is not `off`, which YAML 1.1 parsers read as the boolean false). With `disabled`, `reallocation_surplus` and `density_coupling.reallocation_gradient` are inert and warn |
| `groundwater.predictor_corrector` | `legacy` \| `consistent` | no | `consistent` (S3-A18 flip, 2026-11-04; the T11 opt-in was legacy) | Stage-3 T11 (O17): whether the PCA head solve (predictor) and the Darcy-flux θ update (corrector) share one face rule. `legacy` (golden-pinned) keeps the ported v1 split — the corrector's one-sided face Ks (x/y: the minus cell, z: the lower cell), the ponded top face averaged in the predictor only, dry seepage faces sealed in the predictor but live in the corrector, and the terrain-gravity term without density factors; at saturated Ks contrasts the reallocation restore turns the split into created/dropped water ([theory](../theory/groundwater.md)). `consistent` shares every face rule between the halves (each flank's own Ks; Ksz under ponding; the pressurized-seepage mode 3; the r²v¹ terrain-gravity form) — recommended for heterogeneous and terrain-following density-coupled runs. Bitwise `legacy` on homogeneous uncoupled soil |
| `groundwater.face_conductivity` | `arithmetic` \| `harmonic` | no | `arithmetic` | Stage-3 T11 (O17b; needs `predictor_corrector: consistent`): the interior face-conductivity mean. `arithmetic` is the v1 rule, 0.5 (K_a + K_b) of the flanks' own closures. `harmonic` takes the thickness-weighted harmonic mean of the saturated conductivities times the arithmetic mean of the relative permeabilities — exact for saturated series flow (the arithmetic interface face carries ~25× the exact conductance at a 100× contrast), seals a Ks = 0 flank exactly, and wetting fronts do not stall (the pc3–pc5 gates) |
| `groundwater.saturated_contact` | `legacy` \| `head` | no | `head` (S3-A19 flip, 2026-11-16; the T12 opt-in was legacy) | Stage-3 T12 (S3-A19): how the post-allocation sweep decides that a cell touches saturation, which selects the θ := θ(h) restore over the isolated h := h(θ) one. `legacy` (golden-pinned) is the v1 test θ ≥ θs, read in place while the same sweep restores the neighbours: the restored set follows the loop order, so results depend on the MPI decomposition and the thread schedule (b5: 2.4e-5 field-relative between 1 and 5 ranks after one step). `head` reads every neighbour from the pre-sweep state and counts it as saturated when θ ≥ θs or its head is above air entry (where θ(h) = θs) — rank-, thread- and order-invariant ([theory](../theory/groundwater.md)) |
| `groundwater.reallocation_surplus` | `drop` \| `redistribute` | no | `drop` | post-allocation surplus of saturation-adjacent cells: discard (legacy sweep, groundwater.c:1026 disabled the transfer; the b2 golden embeds it) or redistribute vertically into pore room (mass-conserving; b3 uses it) |
| `groundwater.density_coupling.enabled` | bool | no | false | baroclinic feedback (requires transport and/or temperature): r_rho = 1 + β_s·s − β_T·(T − T₀), r_visc = 1/(1 + β_sv·s) |
| `groundwater.density_coupling.beta_saline` | real ≥ 0 | no | 7.44e-4 | β_s [L/g] (v2 Q5: the legacy compile-time constant, now configurable; requires transport) |
| `groundwater.density_coupling.beta_saline_viscosity` | real ≥ 0 | no | 2.2e-3 | β_sv [L/g] (requires transport) |
| `groundwater.density_coupling.thermal_expansion` | real ≥ 0 | no | 0 | β_T [1/K] thermal expansion (v2 Q5; requires temperature; 0 = thermally passive) |
| `groundwater.density_coupling.reference_temperature` | real | no | 20 | T₀ [°C] for the thermal term (requires temperature) |
| `groundwater.density_coupling.reallocation_gradient` | `freshwater` \| `density` | no | freshwater | vertical head gradient of the post-allocation split. `freshwater` (legacy, golden-pinned) tests dh/dz − 1, which reads a hydrostatic saline column (dh/dz = r_ρ) as an upward gradient: when a saturated layer separates two cells in saturated contact, their restore surpluses are swapped through it every step and the scalar cannot follow (the Nueces coupled case lost 70 % of its infiltrated salt). `density` tests dh/dz − r_ρ, the gravity term of the Darcy flux — recommended with density coupling |
| `groundwater.evaporation.mode` | `bulk` | yes (in block) | — | v2 Q4 bulk-aerodynamic soil evaporation: potential rate from the `atmosphere` block, actual rate limited by the surface-layer soil relative humidity α₁ = min(1, 1.8·w_g/(w_g + 0.30)) (Geng & Boufadel Eq. 6; condensation passes through), applied as the top-face flux over the region. Requires `atmosphere`; uncoupled groundwater runs only; the region must not overlap a `groundwater_top` condition. Adds the cumulative `evaporation` column to `/monitor/gw_mass_audit` |
| `groundwater.evaporation.region.polygon` | [x, y] list ≥ 3 | yes (in block) | — | the evaporation zone |

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
| `soil.types[].theta_fc` | real in [0, 1] | no | θ(h = −3.3 m) from van Genuchten | field capacity for the `beta_lee_pielke` ET limiter (Stage-3 T5) |
| `soil.types[].albedo` | real in [0, 1] | with PM soil channel | — | soil surface albedo (penman_monteith, Stage-3 T7: required on every type when PM runs the soil channel) |
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
| `initial_conditions.temperature.surface.{constant,file}` | one of | when temperature+SWE | — | surface temperature [°C] |
| `initial_conditions.temperature.groundwater.{constant,file}` | one of | when temperature+GW | — | subsurface temperature [°C] |

## boundary_conditions (list; plan §5.6)

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `boundary_conditions[].name` | string | yes | — | unique name |
| `boundary_conditions[].region.polygon` | list of [x, y], ≥ 3 | yes | — | region vertices [m]; on-edge counts inside |
| `boundary_conditions[].target` | `surface` \| `groundwater_top` \| `groundwater_bottom` \| `groundwater_side` | yes | — | where the condition applies |
| `boundary_conditions[].kind` | `eta` \| `discharge` \| `velocity` \| `outflow` \| `head` \| `flux` \| `scalar_value` \| `scalar_cauchy` | yes | — | eta/discharge/velocity/outflow: surface only; head/flux: groundwater only; scalar_value and scalar_cauchy need transport. `outflow` is free (transmissive) outflow at domain-edge faces — the stage is extrapolated down the continued bed slope — and takes no `value`. `scalar_value` semantics (P4): target `surface` — members covered by a `discharge` condition inject that inflow's concentration (legacy `s_inflow`); every other member is a wet-cell Dirichlet (the legacy tide salinity, generalized to any region — dry cells stay at 0); target `groundwater_side` — the ghost concentration along the member edge (legacy `s_yp`/`s_ym`), feeding upwinded inflow, the limiter bound, and the baroclinic boundary face. `scalar_cauchy` (v2 Q4): zero-total-scalar-flux top condition (Geng & Boufadel 2015 Eq. 7) — water crosses the subsurface top face (whatever the flow-side top condition sends), scalar mass does not, and the top-cell limiter admits the exact evaporative concentration / infiltration dilution; homogeneous (takes no `value`), target `groundwater_top` in uncoupled groundwater runs only |
| `boundary_conditions[].scalar` | `salinity` \| `temperature` | no | `salinity` | which registered scalar a `scalar_value` condition prescribes (v2 Q5); valid on kind `scalar_value` only. With `scalar: temperature`, targets `groundwater_top`/`groundwater_bottom` become valid as **cell-pinning** Dirichlet rows: the member columns' top/bottom cells are re-imposed after every update (the b6 tide rule applied vertically — the g6/g8 gates). The salinity top/bottom rows stay schema-rejected in v2.0 (no gate exercises them; §8.2 discipline) |
| `boundary_conditions[].value.constant` | real | one of | — | fixed value |
| `boundary_conditions[].value.series.file` | path | one of | — | time series |
| `boundary_conditions[].inflow.closure` | `zero_gradient` \| `external` | no | `zero_gradient` | Stage-3 T2, kind `eta` only: the momentum state beyond a prescribed-stage boundary. `zero_gradient` = v2.0 (inflow imports the interior momentum — use where the boundary cuts a uniform flow); `external` = the stage cells are an external reservoir (entrance pays (v²−u_ext²)/2g; coupled deposits leave through the reservoir, booked as `stage_deposit`; see theory/surface-water.md) |
| `boundary_conditions[].inflow.velocity.u` / `.v` | constant \| series [m/s] | no | 0 (still) | the external reservoir's velocity components (`closure: external` only) |
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
| `transport.legacy_evap_allowance` | bool | no | false | v2 Q4: keep the legacy hardcoded +0.01/step limiter allowance on coupled dry evaporating columns instead of the exact in-step concentration factor (golden pinning) |
| `transport.subsurface_update` | `legacy` \| `consistent` | no | consistent when the resolved `surface_update` is conservative (its coupled books need it), else legacy | subsurface salinity update (requires groundwater). `legacy` (golden-pinned) advects with the last subsurface substep's fluxes over the surface dt but divides by that substep's flux volume, so subcycled runs with dtg ≠ dt lose or gain salt outside every audit column, surface steps without a substep reuse stale fluxes, the surface-interface diffusion acts only on steps with seepage, and post-allocation transfers are recorded on the receiver's face only. `consistent` moves the scalar with every volume the groundwater module moved: face volumes summed over all substeps (transfers on every face crossed), the previous transport step's water content as the mass basis, the interface diffusion on every wet step, and no one-sided interface gain. Identical to `legacy` up to rounding when dt = dtg and nothing is reallocated; recommended for coupled salinity runs |
| `transport.surface_dryout` | `discard` \| `store` | no | **store** (T8 default flip, S3-A11; v2.0 was discard) | a drying surface cell's remaining salt: `discard` (the v2.0 behavior, golden-pinned — lost into surf_adjust; token renamed from the plan's `delete` by S3-A8.10, the T0 forbidden-word rule) or `store` (Stage-3 T5/O4 — since T9 (O15) the PAIRED water–salt store — every thin-layer removal parks its water together with its salt, rain and evaporation are untouched, parked water is inert, and the store returns at rewetting or once it reaches `wetting_face_depth` (`theory/transport.md`); the `water_store`/`salt_store` outputs and the `store_volume`/`surf_store` audit columns expose it). The default-bitwise lanes and goldens pin `discard` explicitly |
| `transport.surface_update` | `legacy` \| `consistent` \| `conservative` | no | **conservative** (S3-A18 flip, 2026-11-04; v2.0–T10 default was legacy) | surface salinity update. `legacy` divides the transported mass by the previous step's flux volume and clamps moving cells only (golden-pinned, [theory quirk 2](../theory/transport.md)); under unsteady wetting/drying it concentrates closing thin films above every input value and books the mismatch into `surf_anchor`. `consistent` divides by the current step's flux volume (the flow rates that moved the mass) and applies the local min/max limiter in every wet cell — recommended for tidal and wind-driven wetting/drying. `conservative` further advects with the free-surface solve's own continuity face volumes — the water that actually moved, including across faces the velocity limiter zeroes after the solve — so salt crosses drying faces with its water, a uniform scalar stays uniform through wetting/drying, and `surf_anchor` falls to the solve-residual class ([theory](../theory/transport.md)); recommended with `surface_dryout: store` for coupled tidal wetting/drying |

## temperature (v2 Q5; required section when the module is on)

The second registered scalar (plan §4.1). The transport machinery
(advection schemes, limiter, ledgers, decomposition staging) is shared
with salinity; these are the per-scalar parameters and the thermal
physics: the subsurface mass basis is θ + κ with the retardation
κ = (1 − θ_s)·(ρc)_s/(ρc)_w, conduction enters the dispersion tensor's
molecular slot as λ_eff/(ρc)_w (a bulk property, **not** scaled by θ_s),
uncoupled top faces advect the donor value (heat travels with the
water), and rain/evaporation change volume without changing temperature.

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `temperature.scheme.advection` | `upwind` \| `superbee` | no | upwind | advection scheme |
| `temperature.surface_diffusivity.x` | real ≥ 0 | no | 1e-10 | [m²/s] (needs SWE) |
| `temperature.surface_diffusivity.y` | real ≥ 0 | no | 1e-10 | [m²/s] (needs SWE) |
| `temperature.thermal_conductivity` | real > 0 | when GW on | — | λ_eff [W/m/K] bulk effective thermal conductivity |
| `temperature.heat_capacity_water` | real > 0 | no | 4.184e6 | (ρc)_w [J/m³/K] |
| `temperature.heat_capacity_solid` | real > 0 | when GW on | — | (ρc)_s [J/m³/K] solid volumetric heat capacity (the retardation) |
| `temperature.dispersivity.longitudinal` | real ≥ 0 | no | 0 | thermal dispersivity [m] (needs GW) |
| `temperature.dispersivity.transverse` | real ≥ 0 | no | 0 | thermal dispersivity [m] (needs GW) |
| `temperature.bounds.min` | real | no | unbounded | optional lower bound [°C] (unlike salinity, both bounds default open) |
| `temperature.bounds.max` | real | no | unbounded | optional upper bound [°C] |
| `temperature.surface_exchange.mode` | `equilibrium` \| `bulk` | yes (in block) | — | surface heat flux (needs SWE). `equilibrium`: Q = −K_e (T − T_e) (Edinger). `bulk`: Q_net = Q_sw + ε(LW_in − σT_K⁴) − Q_lat − Q_sens with ε = 0.97 and the latent/sensible terms on the shared Q4 atm chain, evaluated at the local water temperature; requires `atmosphere`. Wet cells only; audited in the `surf_atmos` ledger column |
| `temperature.surface_exchange.equilibrium.temperature` | real | with mode equilibrium | — | T_e [°C] |
| `temperature.surface_exchange.equilibrium.coefficient` | real > 0 | with mode equilibrium | — | K_e [W/m²/K] |

Output variables: `temperature` (subsurface, °C) and
`temperature_surface` under `output.variables.temperature`;
`temperature_surface` is also a valid point-monitor variable. An active
temperature module writes the `/monitor/temperature_audit` heat ledger
(the transport-audit identity with the extra `surf_atmos` column).

## output

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `output.filename` | string | yes | — | single HDF5 file (plan §7) |
| `output.variables.surface` | list ⊆ {eta, depth, uu, vv, seepage, seepage_rate, film_removal_rate, water_store} | no | — | requires the SWE module; `seepage` and `seepage_rate` also require the GW module. `seepage` is a PER-STEP sample of the applied rate (m/s) — do not integrate it; `seepage_rate` is the interval MEAN of the complete coupled surface gain (applied seepage + bounce-back + vent deposits, m/s, positive upward) whose integral `sum(rate · area · Δt_out)` recovers the audit's `seepage` column exactly (Stage-3 T4). `film_removal_rate` is the interval-mean thin-layer removal per cell (m/s; Stage-3 T9, every mode); `water_store` is the paired store's parked water per cell (m³; Stage-3 T9, needs `transport.surface_dryout: store`) |
| `output.variables.groundwater` | list ⊆ {hydraulic_head, water_content, qx, qy, qz, realloc_created, realloc_dropped, vloss, root_uptake} | no | — | requires the GW module; the three Stage-3 cumulative diagnostics [m³] allocate their fields only when requested; `root_uptake` (T6) is the interval-MEAN root-uptake sink [1/s] on the T4 registry semantics (`realloc_*` are identically zero — and warn — with `reallocation: disabled`) |
| `output.variables.transport` | list ⊆ {concentration, concentration_surface, salt_store} | no | — | requires transport; `concentration` is the subsurface scalar (with the GW module), `concentration_surface` the surface scalar (with the SWE module) |
| `output.monitors[].name` | string | yes | — | unique monitor name |
| `output.monitors[].type` | `point` | no | point | monitor kind |
| `output.monitors[].i` | int ≥ 0 | yes | — | global cell index in i (< nx) |
| `output.monitors[].j` | int ≥ 0 | yes | — | global cell index in j (< ny) |
| `output.monitors[].variables` | list, ≥ 1 | yes | — | variables of enabled modules |
| `output.variables.atmosphere` | list ⊆ {evaporation_rate, soil_evaporation_rate, soil_evaporation_unrealised_rate, wet_time, soil_evaporation_time, transpiration_rate, transpiration_unrealised_rate} | no | — | the ET observables (Stage-3 T5): interval-mean rates [m/s] (the T4 mechanism) and the cumulative complementarity clocks [s] |
| `output.monitor_interval` | real ≥ 0, whole seconds | no | 0 | monitor/audit ROW-WRITE decimation (Stage-3 T4). 0 = a row every step (the v2.0 behavior, bitwise). > 0: audit columns still ACCUMULATE every step, but rows are written only at the first step reaching each interval multiple, at checkpoint instants, and at `t_end` — every written row is bitwise the row the every-step run writes at the same time. Point monitors sample at the same instants. On restart the schedule is recomputed from the restart time like the output schedule |
| `output.checkpoint.interval` | real ≥ 0, whole seconds | no | 0 | 0 = off (no checkpoints at all); > 0 also writes a final checkpoint at `t_end` |

## restart

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `restart.enabled` | bool | no | false | resume from a checkpoint |
| `restart.file` | string | when enabled | "" | checkpointed output file (must exist) |
| `restart.time` | real | no | 0 | checkpoint time to resume from [s] |

## solver

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `solver.surface.preconditioner` | `bjacobi-icc` \| `amg` \| `gamg` | no | bjacobi-icc | free-surface (`fs_`) preconditioner: the v1 block-Jacobi/ICC(0) default, hypre BoomerAMG (`amg`, needs a hypre-enabled PETSc), or PETSc GAMG (v2 plan §2.2) |
| `solver.surface.mat_type` | `aij` \| `aijkokkos` | no | aij | `fs_` PETSc matrix/vector backend: host `aij`, or `aijkokkos` running the solve through Kokkos Kernels on the build's execution space (the threaded-CPU/GPU lane; needs a Kokkos-enabled PETSc). Device builds force `aijkokkos` (v2 plan §2B.2) |
| `solver.surface.rtol` | real > 0 | no | 1e-8 | `fs_` relative tolerance |
| `solver.surface.atol` | real >= 0 | no | 1e-14 | `fs_` absolute tolerance |
| `solver.surface.max_iterations` | integer >= 1 | no | 500 | `fs_` iteration cap |
| `solver.surface.reuse_max_solves` | integer >= 0 | no | 50 | amg/gamg hierarchy reuse cadence: rebuild at least every this many solves (0 = every solve); ignored for bjacobi-icc |
| `solver.surface.reuse_iteration_factor` | real >= 1 | no | 1.5 | early rebuild when a solve exceeds this factor times the post-rebuild iteration count |
| `solver.groundwater.preconditioner` | `bjacobi-icc` \| `amg` \| `gamg` | no | bjacobi-icc | Richards (`gw_`) preconditioner (same choices; the amg defaults use the 3D anisotropy parameters, v2 plan §2.2.3) |
| `solver.groundwater.mat_type` | `aij` \| `aijkokkos` | no | aij | as `solver.surface.mat_type`, for `gw_` |
| `solver.groundwater.rtol` | real > 0 | no | 1e-8 | `gw_` relative tolerance |
| `solver.groundwater.atol` | real >= 0 | no | 1e-14 | `gw_` absolute tolerance |
| `solver.groundwater.max_iterations` | integer >= 1 | no | 1000 | `gw_` iteration cap |
| `solver.groundwater.reuse_max_solves` | integer >= 0 | no | 50 | as `solver.surface.reuse_max_solves`, for `gw_` |
| `solver.groundwater.reuse_iteration_factor` | real >= 1 | no | 1.5 | as `solver.surface.reuse_iteration_factor`, for `gw_` |
| `solver.petsc_options_file` | string | no | "" | PETSc options file (overrides every `fs_`/`gw_` default above, plan §5.1) |

## runtime

| Key | Type | Required | Default | Description |
|---|---|---|---|---|
| `runtime.gpu_aware_mpi` | `auto` \| `on` \| `off` | no | auto | pass device buffers to MPI or stage through host mirrors (plan §5.3) |
