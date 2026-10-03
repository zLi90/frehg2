# Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `INVALID: ... (N errors)` at load | schema violation. Read each `- ...` line; unknown keys print a suggestion. Run `--validate` until clean. |
| Run is extremely slow on a small grid | thread launch overhead — set `OMP_NUM_THREADS=1` ([running](running.md#threads-and-performance)). |
| `transport: true` rejected | a flow module must be on, the `transport:` section and scalar initial conditions must exist, and `density_coupling` additionally requires transport and/or temperature. |
| `temperature: true` rejected | as for transport: a flow module, the `temperature:` section, and `initial_conditions.temperature` for each active grid; with groundwater on, `thermal_conductivity` and `heat_capacity_solid` are required ([configuration](configuration.md#temperature-heat-transport-v2)). |
| Missing `atmosphere` error | `surface_water.evaporation.mode: bulk`, `groundwater.evaporation`, and `temperature.surface_exchange.mode: bulk` all read the `atmosphere:` block ([configuration](configuration.md#atmosphere-meteorological-forcing-v2)). |
| `preconditioner: amg` or `mat_type: aijkokkos` fails at startup | `amg` needs a hypre-enabled PETSc and `aijkokkos` a Kokkos-enabled PETSc; build PETSc as the [installation guide](installation.md) describes, or use `gamg` / `aij`. |
| A salinity side boundary does not raise the concentration | only the y+ (north) side admits a prescribed subsurface salinity; on the other sides it is clipped to the interior range. Orient the case so salt enters from y+ ([configuration](configuration.md#transport-scalar-transport)). |
| Temperature pin on `groundwater_top` rejected | in coupled runs the coupler owns the top face; top/bottom temperature pins are for uncoupled groundwater runs (a bottom pin works in both). |
| `groundwater_top` condition rejected in a coupled run | the coupler owns the top head; configured top conditions must be `kind: flux` ([configuration](configuration.md#coupling-coupled-runs)). |
| Sync-coupled run ignores my `time.dt` | by design: `dt` is only the initial value of the common adaptive step; bound it via `groundwater.timestep.dt_min`/`dt_max`. |
| Subcycled coupled run is inaccurate / CFL warnings | the surface marches at the *fixed* `time.dt` in subcycled mode; reduce `time.dt` (or use `mode: sync`). |
| Value-count mismatch reading a raster | the file has the wrong number of values for `nx·ny`; check ordering (`j*nx + i`) and that header `ncols`/`nrows` match `nx`/`ny`. |
| MPI or HDF5 errors at output time | HDF5 was not built parallel, or Kokkos/PETSc/HDF5 were built against a different MPI. Rebuild the stack against one MPI. |
| CMake can't find Kokkos/PETSc/HDF5 | add their install prefix to `-DCMAKE_PREFIX_PATH`. |
| Build fails with compiler-wrapper errors | don't set `CMAKE_CXX_COMPILER=mpicxx`; set `CXX` to the real compiler ([installation](installation.md#a-note-on-compiler-wrappers)). |
| A test or run hangs at exit with ~100 % CPU | MPICH's libfabric *sockets* provider can wedge `MPI_Finalize` (seen on macOS). Run with `FI_PROVIDER=tcp` — the test harness pins this itself. |
| Domain looks dry everywhere | `eta` starts below `bottom` (a dry start is intentional); water appears once rain/inflow accumulates. Check `depth` output, not just `eta`. |
| Groundwater fields are NaN in some cells | those cells are outside the active subsurface (above a terrain-following column's bed); this is expected, mask them when plotting. |
