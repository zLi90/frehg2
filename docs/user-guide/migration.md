# Migrating from legacy Frehg

Frehg2 replaces the legacy flat `key = value` input file with the strict v2
YAML schema. There is no automatic converter from legacy `input` files (the
legacy format has no reliable machine-readable contract); this page is the
by-hand mapping. For the *archival v1 YAML drafts* that circulated during
planning, `tools/migrate_yaml_v1_to_v2.py` converts automatically — see the
end of this page.

## Legacy parameter → v2 key map

This mirrors Appendix A of the
[upgrade plan](../developer-guide/FREHG2_UPGRADE_PLAN.md), which is the
binding version.

| Legacy key(s) | Frehg2 YAML | Notes |
|---|---|---|
| `NX NY dx dy dz dz_incre botZ` | `domain.{nx,ny,nz,dx,dy,dz,dz_stretch,bottom_elevation}` | `nz` was computed from the bathymetry range in legacy; explicit in Frehg2 |
| `mpi_nx mpi_ny use_mpi nthreads` | `domain.decomposition.{mpi_nx,mpi_ny}` (default auto) | `use_mpi`/`nthreads` dropped (always MPI + Kokkos) |
| `dt Tend dt_out` | `time.{dt,t_end,output_interval}` | `NT` was parsed-but-unused in legacy; no equivalent exists |
| `bath_file actv_file` | `domain.bottom_elevation.file` | the active mask folds into bathymetry nodata |
| `min_dept wtfh hD manning grav viscx viscy` | `surface_water.{min_depth,wetting_face_depth,friction.thin_layer_depth,friction.coefficient,gravity,viscosity.x/y}` | |
| `sim_wind Cw CwT windspd winddir north_angle` | `surface_water.wind.*` | |
| `q_rain rain_file / q_evap evap_file` | `surface_water.{rainfall,evaporation}` constant/series | the aerodynamic `evap_model` is dropped |
| `n_tide tide_locX/Y tide_file init_tide` | `boundary_conditions[]` kind `eta` with polygon regions | legacy rectangles become 4-vertex polygons |
| `n_inflow inflow_locX/Y inflow_file` | `boundary_conditions[]` kind `discharge` | |
| `bctype_GW[6] qtop qbot qyp qym htop hbot` | `boundary_conditions[]` targets `groundwater_top/bottom/side`, kinds `head`/`flux` | legacy face code 3 (free drainage) → `value: {gravity: true}` |
| `Ksx Ksy Ksz Ss wcs wcr soil_a soil_n aev` | `soil.types[].{ksx,ksy,ksz,theta_s,theta_r,vg_alpha,vg_n,aev}`, `groundwater.specific_storage` | heterogeneous soils via `soil.map.file` (new) |
| `use_full3d follow_terrain sync_coupling dt_adjust dt_max dt_min Co_max` | `groundwater.{use_full3d,timestep.*}`, `domain.follow_terrain`, `coupling.mode` | `iter_solve`, `use_corrector`, `post_allocate`, `n_substep` dropped (PCA always runs the corrector; post-allocation always on) |
| `n_scalar baroclinic superbee difux difuy disp_lon disp_lat init_s_* s_tide s_yp s_ym` | `modules.transport`, `groundwater.density_coupling`, `transport.*`, `scalar_value` BCs and `initial_conditions.transport` | single scalar in v1.0 |

Legacy features deliberately **not** carried over (diffusive wave, subgrid
topography, the Newton groundwater path, modified van Genuchten, the
aerodynamic evaporation model, and others) are listed with rationale in
[Removed features](../theory/removed-features.md). If your legacy case used
one of them, there is no v2 equivalent by design.

## Conventions that changed

- **Everything is SI at the interface.** Legacy inputs mixing units (mm/h
  rain, cm heads) must be converted; series and rasters are read as SI.
- **Elevations are absolute.** Legacy worked in offset elevations
  internally; Frehg2 applies that offset transparently. Initial `eta`,
  `bottom_elevation`, and all outputs use one datum.
- **Boundary regions are polygons over cell centers**, not index
  rectangles. A legacy `locX/locY` rectangle becomes a 4-vertex polygon in
  physical coordinates; a point exactly on the polygon edge counts as
  inside.
- **Output is one parallel HDF5 file** keyed by integer seconds, not
  per-variable ASCII files (`ascii_golden_to_h5.py` in `tools/` converts
  legacy ASCII outputs when you need to compare).
- **Solver failures are fatal.** Legacy silently continued on linear-solver
  divergence; Frehg2 aborts with the solver diagnostics.

## Migrating v1 YAML drafts

The v1 experimental YAML drafts are upgraded automatically:

```bash
python3 tools/migrate_yaml_v1_to_v2.py input.v1.yaml -o case.yaml
build/src/frehg --validate case.yaml
```

The rename table in the tool's own docstring is the authoritative contract;
any v1 key not covered there is a hard error — the tool refuses to guess,
so silent key loss is impossible.
