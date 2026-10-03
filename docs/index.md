# Frehg2

Frehg2 is a production-grade rewrite of Frehg 1.0 — a coupled semi-implicit
2D shallow-water / 3D mixed-form Richards / solute- and heat-transport
model — in C++20 with Kokkos (CPU/GPU-portable), MPI, PETSc, yaml-cpp
configuration, and single-file parallel HDF5 output. The core physics is
preserved exactly from the validated legacy code, and the six legacy
benchmarks are blocking validation gates. Version 2.0 adds new physics
and solver options, each with its own gates. The binding specifications
are the [upgrade plan](developer-guide/FREHG2_UPGRADE_PLAN.md) (v1) and
the [v2 development plan](developer-guide/FREHG2_V2_DEVELOPMENT_PLAN.md).

## What Frehg2 simulates

- **Surface water** — the θ-scheme 2D shallow-water solver: rainfall,
  prescribed or bulk-aerodynamic evaporation, wind stress (constant
  drag or Cd(U₁₀) laws), Manning or Chézy friction, wetting/drying, and
  an implicit free surface (PETSc CG).
- **Groundwater** — the mass-conservative PCA mixed-form Richards solver
  (Li et al. 2021): van Genuchten soils, adaptive subsurface stepping,
  head/flux/free-drainage/hydrostatic boundary conditions,
  terrain-following or partial-cell meshes, and bare-soil evaporation
  limited by soil moisture.
- **Coupled surface–subsurface flow** — both modules exchanging
  infiltration and seepage each step, in lockstep (`sync`) or subcycled
  mode, with a per-step conservation audit.
- **Salinity transport** — advected with upwind or TVD-superbee fluxes,
  diffused on the surface, dispersed anisotropically in the subsurface,
  exchanged through the seepage, and optionally fed back into the
  subsurface density/viscosity for baroclinic problems such as saltwater
  intrusion. Evaporation can concentrate it (salinization).
- **Heat transport** (v2) — temperature on both grids, with conduction,
  thermal dispersion, retardation by the solid matrix, surface heat
  exchange (equilibrium or bulk-aerodynamic with radiation), and optional
  thermal density feedback.
- **Atmospheric forcing** (v2) — one meteorological block feeding the
  open-water and bare-soil evaporation and the surface heat exchange.

Every mode has HDF5 field output, mass-audit tables, point monitors,
checkpoint/restart (restart reproduces the uninterrupted run, bitwise for
flow and salinity), and
a `run-record.yaml` provenance record. Each linear system can use
block-Jacobi/ICC, BoomerAMG or GAMG preconditioning, and with a
Kokkos-aware PETSc the solves thread too. The GPU path is
**experimental**: all kernels are Kokkos and the CUDA build compiles, but
device execution has not yet been validated.

## Where to start

| I want to… | Read |
|---|---|
| install and build | [Installation](user-guide/installation.md) |
| set up a run | [Configuration](user-guide/configuration.md), [Parameter reference](user-guide/parameters.md) |
| run and analyze | [Running a simulation](user-guide/running.md), [Output reference](user-guide/output.md) |
| reproduce the validation | [Benchmark walkthroughs](user-guide/benchmarks.md) |
| build a case from scratch | [Worked example](user-guide/worked-example.md) |
| port a legacy Frehg case | [Migrating from legacy Frehg](user-guide/migration.md) |
| understand the numerics | the [Theory](theory/surface-water.md) section |
| modify the code | the [Developer guide](developer-guide/architecture.md) |

## Provenance

The hydrodynamics follows Li & Hodges (2019a, b); the variably saturated
groundwater solver follows Li, Özgen-Xian & Maina (2021); the coupled
surface–subsurface salinity modeling that defines Frehg 1.0 — the
validated predecessor this rewrite preserves — is described in Li,
Hodges & Shen (2023, *J. Hydrol.* 618, 129268). The rewrite is
specified, gated, and documented by the
[upgrade plan](developer-guide/FREHG2_UPGRADE_PLAN.md), the
[v2 development plan](developer-guide/FREHG2_V2_DEVELOPMENT_PLAN.md), and
the per-phase records under the developer guide; full citations are in the
repository's `CITATION.cff`.
