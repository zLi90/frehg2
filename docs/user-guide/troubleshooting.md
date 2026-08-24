# Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `INVALID: ... (N errors)` at load | schema violation. Read each `- ...` line; unknown keys print a suggestion. Run `--validate` until clean. |
| Run is extremely slow on a small grid | thread launch overhead — set `OMP_NUM_THREADS=1` ([running](running.md#threads-and-performance)). |
| `transport: true` rejected | a flow module must be on, the `transport:` section and scalar initial conditions must exist, and `density_coupling` additionally requires transport. |
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
