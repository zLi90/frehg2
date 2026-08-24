# Architecture

Frehg2 is a single executable (`frehg`) built from five static libraries
with strict layering — physics libraries never include each other, and
cross-module data flows through driver-owned coupled state:

```
frehg (main.cpp)
  └─ frehg::driver      time loops, module orchestration, output/monitor/
       │                checkpoint mediation, restart
       ├─ frehg::swe        surface-water module
       ├─ frehg::gw         groundwater module
       ├─ frehg::transport  scalar-transport module
       ├─ frehg::coupling   the surface–subsurface coupler and exchange
       └─ frehg::core       config, grid, halo, BC, linear solver, I/O,
            │               time series, timers, logging
            └─ Kokkos, PETSc, MPI, HDF5, yaml-cpp
```

| Directory | Library | Contents |
|---|---|---|
| `src/core` | `frehg::core` | `Types`, `Config`/`ConfigSchema`, `Grid`, `HaloExchanger`, `LinearSystem` (PETSc COO), `PetscSession`, `TimeSeries`, `Timer`, `Logger` |
| `src/bc` | (in core) | `Polygon`, `BoundarySet` (rasterization, per-BC sub-communicators), `BoundaryKinds` |
| `src/io` | (in core) | `Hdf5Output`, `Monitor`, `Checkpoint`, `GridDataReader` |
| `src/swe` | `frehg::swe` | `SurfaceSolver`, `Momentum`, `FreeSurface`, `WetDry`, `SurfaceSources` |
| `src/gw` | `frehg::gw` | `RichardsSolver`, `VanGenuchten`, `Predictor`, `Corrector`, `Reallocate`, `AdaptiveStep`, `TerrainMetric`, `Baroclinic` |
| `src/transport` | `frehg::transport` | `ScalarSolver`, `SurfaceTransport`, `SubsurfaceTransport`, `Dispersion`, `Limiters` |
| `src/coupling` | `frehg::coupling` | `Coupler`, `Exchange` |
| `src/driver` | `frehg::driver` | `Simulation` (the three time loops: surface-only, groundwater-only, coupled) |

## Rules the layering enforces

- **Physics libraries never include each other.** The SWE module knows
  nothing of the subsurface; the transport module reads the flow state
  through driver-wired view structs. All cross-module state (exchange
  fluxes, the scalar the baroclinic ratios read) is owned by the driver or
  the coupler and passed as views.
- **No `Kokkos::SharedSpace`/UVM anywhere.** Fields are device Views
  (`Field2`/`Field3` in `Types.hpp`, `(j, i[, k])` LayoutRight); host
  access happens exclusively via `create_mirror_view` + `deep_copy` at I/O
  boundaries.
- **Zero backend `#ifdef`s in physics code.** Backend differences are
  confined to `HaloExchanger` (the GPU-aware-MPI toggle — the only place it
  exists) and PETSc type strings. This is what makes the CUDA compile-only
  CI lane meaningful.
- **No silent failure.** Every PETSc/HDF5/MPI return code is checked
  (`FREHG_PETSC_CHECK`/`FREHG_H5_CHECK`); all failures route through
  `log::fatal`, which throws `frehg::FatalError` (never `exit()` in library
  code). Linear-solver divergence is fatal by design — the legacy model's
  silent continue was a bug, not behavior to preserve.
- **Memory discipline.** No raw `new`/`delete`/`malloc`/`free` — Kokkos
  Views and standard containers only, enforced by
  `scripts/check_forbidden.sh`.

## Control flow of a step

1. **Driver** advances time (fixed `dt`, or the common adaptive step in
   sync-coupled runs) and evaluates time series.
2. **Coupler** (coupled runs) classifies columns (wet-Dirichlet /
   supply-limited flux / seepage face) and hands the subsurface its top
   boundary; after the subsurface step it returns the accumulated exchange
   volume to the surface.
3. **SWE**: momentum sources (advection, viscosity, drag, wind) → implicit
   free-surface solve (`LinearSystem("fs_")`) → velocity update →
   wetting/drying → rain/evaporation.
4. **GW**: predictor (7-point head system, `LinearSystem("gw_")`) →
   corrector (face K, Darcy fluxes, θ update) → moisture reallocation →
   adaptive-step controller.
5. **Transport** (after the flow step): advective + dispersive/diffusive
   fluxes on both grids, monotonicity clipping, scalar exchange through the
   seepage, and (if enabled) the density/viscosity ratio update the next
   subsurface step reads.
6. **Driver** writes output/monitors/audits at output times and checkpoints
   at checkpoint times.

The equation-level description of each module, with legacy `file:line`
provenance for every preserved algorithm, is in the theory section
([surface water](../theory/surface-water.md),
[groundwater](../theory/groundwater.md),
[exchange](../theory/exchange-flux.md),
[transport](../theory/transport.md)).

## Error handling and logging

`log::info/warn/error` print rank-0-deduplicated messages; `log::fatal`
logs and throws `frehg::FatalError`, which `main` catches and converts to
`MPI_Abort` when running on more than one rank. Tests assert on
`EXPECT_THROW(..., frehg::FatalError)` — never on process exit.

## Timers

`Timer::start/stop` (or `Timer::Scoped`) auto-nest into paths
(`"simulation"`, `"swe/free_surface"`, `"groundwater/solve"` …);
`Timer::report(comm)` prints an MPI-merged min/mean/max table at the end of
every run. The [performance report](performance.md) is built from these.
