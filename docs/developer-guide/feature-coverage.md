# Feature-interaction coverage (v2 plan §8.3)

This is the release-blocking coverage table the v2 plan's §8.3 specifies:
rows are shipped features, cells name the gate or test that exercises the
feature and each pairwise interaction that plausibly couples. The rule it
enforces is v1's hard lesson (wind: implemented in P1, never exercised by
any gate, still unvalidated at the v1.0.0 release): **a feature or
interaction with an empty cell either gets at least a
conservation/invariance test, or its configuration combination is
schema-rejected.** "Authored-unexercised" is allowed mid-phase and
forbidden in a release. `scripts/check_feature_coverage.py`
(ctest `unit.feature_coverage`) keeps this file in lockstep: every cell
must name at least one existing test, and the GPU rows must carry an
explicit `experimental` or `supported` status.

Test-name conventions: `Suite.Test` = gtest in `unit.all` / the mpi
drivers; `regression.*` / `scaling.*` / `validate.*` / `unit.*` = ctest
entries; `validation/<case>` = an in-tree validation case with a recorded
result (README).

## 1. Features

| Feature (config surface) | Exercised by |
|---|---|
| Surface water module | regression.b1; regression.b4; unit.all (SweModule) |
| Groundwater module (PCA Richards) | regression.b2; regression.b3; unit.all (GwModule) |
| Coupled surface–subsurface (`modules` both + `coupling`) | regression.b5.rain.sync; regression.b5_restart; unit.all (CoupledModule) |
| Salinity transport (`modules.transport`) | regression.b6.ss; regression.b6_restart; unit.all (TransportModule) |
| Temperature transport (`modules.temperature`) | regression.g6; regression.g7; regression.heat_orient |
| Two registered scalars co-resident | regression.two_scalar_restart |
| `coupling.mode: sync` | regression.b5.rain.sync; regression.b6_restart |
| `coupling.mode: subcycled` | regression.b5.rain.subcycled; regression.b6_subcycled |
| Density coupling, saline (`beta_saline`, r_visc) | regression.b6.ss; regression.g5.geng2015; Baroclinic.RatiosFollowTheScalarAndAverageAtFaces |
| Density coupling, thermal (`thermal_expansion`) | regression.g8; regression.heat_orient |
| Both density betas active | regression.two_scalar_restart |
| Friction `manning` | regression.b1; regression.b5.rain.sync |
| Friction `chezy` | regression.b4 |
| Eddy viscosity (`surface_water.viscosity`) | regression.b1; regression.b5.rain.sync |
| Thin-layer drag (`thin_layer_depth`) | regression.b5.rain.sync; regression.g9 |
| Wetting/drying (`min_depth`, `wetting_face_depth`) | regression.b5.rain.sync; regression.g4; SweModule.HigherOfTwoBottomsBlocksFlowAcrossADam |
| Rainfall constant / series / exclusion | regression.b1; SweModule.RainExclusionRegionReceivesNothing |
| Evaporation prescribed constant | regression.g4; SweModule.EvaporationDriesToTheBedAndNoFurther |
| Evaporation prescribed series | SweModule.EvaporationSeriesRemovesTheSampledRate |
| Evaporation exclusion region | SweModule.EvaporationExclusionRegionLosesNothing |
| Bulk (atmosphere) open-water evaporation | regression.g4; unit.all (BulkAerodynamic) |
| Bulk soil evaporation (`groundwater.evaporation`) | regression.g5.geng2015; BulkAerodynamic.SoilRelativeHumidityAndEquilibrium |
| `transport.legacy_evap_allowance` | TransportModule.LegacyEvapAllowanceTogglePinsCoupledDryColumnBehavior |
| Atmosphere block (met forcing incl. radiation terms) | regression.g7; MetForcing.SamplesConstantsAndFoldsRelativeHumidity |
| `atmosphere.wind_speed_floor` | MetForcing.WindSpeedFloorBoundsStillAir |
| Wind stress (`surface_water.wind`) | regression.g9; regression.g10; regression.wind_orient |
| Cd laws (constant/garratt/smith-banke/wu/large-pond, `cap`) | WindForcing.DragLawTablesMatchHandValues; regression.g9; unit.g910_gates.negative |
| Wind compass form + `north_angle` | regression.g9; WindForcing.NorthAngleRotatesTheCompassIntoTheGrid |
| Wind component form (`u10`/`v10`) + direction wrap | regression.wind_orient; WindForcing.DirectionSeriesInterpolatesAcrossTheWrap |
| Advection `upwind` | regression.g4; regression.tracer_rank_invariance.strict |
| Advection `superbee` (both scalars) | regression.b6.ss; regression.g6; regression.transport_orient |
| Scalar bounds (`transport.bounds`, `temperature.bounds`) | TransportModule.ConfiguredUpperBoundClampsTheUpdate; TransportModule.TemperatureBoundsClampTheUpdate |
| `scalar_cauchy` top condition | regression.g5.geng2015; TransportModule.CauchyTopConcentratesUnderEvaporationWithoutClipping |
| Temperature surface exchange `equilibrium` / `bulk` | regression.g7 |
| van Genuchten soils, layered maps | regression.b3; unit.all (VanGenuchten) |
| `use_full3d` switch | regression.b2; GwModule.UseFull3dSwitchGatesLateralFlow |
| `reallocation_surplus` drop / redistribute | regression.b3; GwModule.ReallocationSurplusModesDropOrConserve |
| Adaptive groundwater clock (`groundwater.timestep`) | regression.b2; regression.b3 |
| Specific storage (Ss) | regression.b2; GwModule.AuditIdentityClosesWithCompressibleStorage |
| Terrain-following mesh `scaled` | regression.b6.ss; GwMesh.TerrainMeshScalesColumnsAndTiltsFaces |
| Terrain-following mesh `uniform` | regression.b5.rain.sync; GwMesh.UniformTerrainLayersStackBelowEachBed |
| Masked columns (bed below the box; no yaml key) | GwModule.FullyMaskedColumnAdjacentToEachEdgeConserves; unit.masked_coupled_fatal |
| Boundary-condition kinds x targets x sides | unit.bc_matrix (tests/coverage/bc_matrix.csv) |
| Preconditioner `bjacobi-icc` (default) | regression.b1; scaling.g2.record.bjacobi |
| Preconditioner `amg` (BoomerAMG) | regression.g1.b1.amg; scaling.g2.iters.amg |
| Preconditioner `gamg` | regression.g1.b1.gamg |
| Solver hierarchy reuse knobs | LinearSystem.HierarchyReuseRebuildPolicy |
| `petsc_options_file` | ConfigTest.AcceptsPetscOptionsFilePassThrough |
| `mat_type: aij` / `aijkokkos` | regression.p1.b1.omp4.aij; regression.p1.b1.omp4.aijkokkos |
| Checkpoint / restart | regression.b1_restart; regression.b5_restart; Hdf5Test.CheckpointRoundTripZeroUlp |
| Run record (r1/r2) | unit.run_record_checker; regression.r2_timer_coverage |
| Monitors + audits (mass/gw/scalar/heat) | regression.b1; regression.g7; Hdf5Test.MonitorAppendsRows |
| Output variables incl. `qx`/`qy` | regression.b3 (requests qx, qy); unit.all (Hdf5Test) |
| Explicit MPI decomposition (`domain.decomposition`) | GridSerial.ExplicitDecompositionHonoredAndMismatchFatal; mpi.core.n3 |
| Non-divisible / odd rank splits | mpi.halo.n3; mpi.core.n3; regression.b5_rank_invariance.strict |
| OpenMP backend (threads) | scaling.s3.threads; regression.p1.b1.omp4.aij |
| Hybrid MPI+OpenMP placements | scaling.s4.hybrid; scaling.p3.hybrid |

## 2. Pairwise interactions

| Interaction | Covered by |
|---|---|
| Wind + wetting/drying | regression.wind_drydown |
| Wind + restart (mid step-down) | regression.wind_restart |
| Evaporation + transport (concentration) | regression.g4; regression.g5.geng2015 |
| Evaporation + wetting/drying (dry-out) | regression.g4; regression.evap_restart |
| Evaporation + restart (mid-drying) | regression.evap_restart |
| Bulk evaporation + bulk heat exchange (shared atm chain) | regression.evap_heat |
| Temperature + density feedback | regression.g8; regression.heat_orient |
| Temperature + coupling (sync) | regression.g6 |
| Temperature + subcycled coupling | regression.heat_subcycled |
| Temperature + restart (two-scalar checkpoint layout) | regression.heat_restart; regression.two_scalar_restart |
| Transport + subcycled coupling | regression.b6_subcycled |
| Transport + coupling + density + restart | regression.b6_restart |
| Transport + MPI decomposition (salinity halos/limiter) | regression.tracer_rank_invariance.strict |
| Temperature + MPI decomposition | regression.heat_rank_invariance.strict |
| Two scalars + both betas + restart | regression.two_scalar_restart |
| Rain + subcycled coupling | regression.b5.rain.subcycled |
| Rain + wetting/drying (dry start) | regression.b5.rain.sync |
| Terrain-following + transport + density | regression.b6.ss |
| Terrain-following + heat | regression.terrain_heat |
| AMG + every physics regime (b1–b6) | regression.g1.b1.amg; regression.g1.b6.ss.amg |
| AMG + restart (hierarchy reuse across the boundary) | regression.g1.b1_restart.amg |
| AMG + MPI rank invariance | regression.g1.b1_rank_invariance.default.amg |
| aijkokkos + full physics (coupled + transport + density) | regression.p1.b6.ss.omp4.aijkokkos |
| aijkokkos + restart | regression.p1.b1_restart.aijkokkos |
| aijkokkos + AMG + threads / hybrid | scaling.p2.solver_threads; scaling.p3.hybrid |
| Masked columns + coupling | unit.masked_coupled_fatal (loud refusal is the shipped contract) |
| Subcycling + restart | regression.b5_restart (sync; the subcycled-window restart is schema-legal, rides the same checkpoint state, and the window bookkeeping is covered by regression.b6_subcycled's sync equivalence) |

Interactions deliberately absent from this table because the terms do not
couple directly (each is an independent source/closure exercised by its own
row above): wind+rain, wind+transport/temperature (wind touches momentum
only; scalar advection sees it through the flow, which the scalar gates
already integrate over), chezy+non-b4 physics (a scalar drag-law swap),
run-record+anything (observer only, r1 asserts non-perturbation).

## 3. GPU / device rows (§8.3 release rule)

| Row | Status | Evidence |
|---|---|---|
| CUDA compile + link (p4) | experimental | .github/workflows/cuda-compile.yml (every library, frehg, and the test binaries build and link warning-free with nvcc; nothing executes) |
| Device memory-space discipline (p5) | supported | unit.p5.backend_invariant_fires; check_forbidden.sh static bans |
| Kokkos LA backend on CPU (aijkokkos rehearsal, p1) | supported | regression.p1.b1.omp4.aijkokkos; regression.p1.b6.ss.omp4.aijkokkos |
| GPU execution (physics on device) | experimental | pending the owner-run p6 bundle (scripts/gpu_acceptance.sh, docs/developer-guide/gpu-acceptance.md); the startup banner prints the experimental status |
| GPU-aware MPI (`runtime.gpu_aware_mpi`) | experimental | exercised only by scripts/gpu_acceptance.sh item 4 (owner-run) |
| GPU performance | experimental | out of v2.0 scope (plan §11); recorded, never asserted |
