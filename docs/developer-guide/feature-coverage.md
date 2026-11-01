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
| `transport.subsurface_update` (`legacy` / `consistent`) | TransportModule.ConsistentSubsurfaceUpdateClosesTheSubcycledLedger; TransportModule.ConsistentSubsurfaceUpdateReducesToLegacyWithoutSubcycling; ConfigTest.SubsurfaceUpdateAndReallocationGradientParseAndRoundTrip; ConfigTest.RejectsConsistentSubsurfaceUpdateWithoutGroundwater |
| `groundwater.density_coupling.reallocation_gradient` (`freshwater` / `density`) | TransportModule.DensityReallocationGradientKeepsInfiltratedSalt; ConfigTest.SubsurfaceUpdateAndReallocationGradientParseAndRoundTrip |
| `transport.surface_update` (`legacy` / `consistent`) | TransportModule.ConsistentSurfaceUpdateKeepsDrainingFilmsWithinTheInputRange; TransportModule.ConsistentSurfaceUpdateShrinksTheAnchorInSloshingFlow; ConfigTest.SurfaceUpdateParsesRoundTripsAndRejectsUnknownValues |
| `transport.surface_update: conservative` (Stage-3 T10) | SweModule.SolveFaceVolumesTelescope; TransportModule.SolveConsistentUniformStaysUniform; TransportModule.SolveConsistentDryingFaceCarriesSalt; regression.s3.sx1/sx2/sx4/sx5/sx-restart (per-PR), sx3/sx_rank/sx_seam (nightly) |
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
| `groundwater.reallocation` legacy / disabled (Stage-3 T1) | regression.s3.rl2 (both lanes, per-cell-sum identity, legacy churn anti-vacuity); regression.s3.rl3 (vadose steadiness/uniformity/cost bounds); regression.s3.rl4 (coupled-salt anchor reduction); regression.s3.rl5.orient_off (vadose 8-orientation battery, strict); regression.s3.rl5.coupled_salt_rank_off (strict rank invariance); GwModule.ReallocationOffConservesWaterExactly; GwModule.ReallocationOffKeepsIsolatedHeadsOnTheRetentionCurve; ConfigTest.ReallocationModeParsesRoundTripsAndRejectsLegacyOnlyKeys |
| Per-cell realloc/vloss diagnostics + checkpoint (Stage-3 W1.3) | regression.s3.rl2 (sums = audit columns); regression.s3.rl-restart (bitwise restart of the cumulative fields) |
| Inflow closure `zero_gradient` / `external` on kind eta (Stage-3 T2) | regression.s3.ob1 (discrete-chain entrance head, 3 lanes incl. the zero_gradient real-negative); regression.s3.ob3 (ledger + tidal prism, both closures); regression.s3.ob5 (no-free-energy vs zero_gradient contrast); regression.s3.ob8.b6.{ss,td} (nightly: b6 under external vs the experiment); SweModule.ExternalStageInflowPaysTheBernoulliEntranceHead; SweModule.ExternalReservoirVelocityEntersTheEntranceForm; ConfigTest.InflowClosureParsesRoundTripsAndRejectsMisuse |
| Stage-deposit routing + `stage_deposit` audit (Stage-3 W2.4, closure × coupled) | CoupledModule.ExternalStageAbsorbsAndBooksCoupledDeposits; CoupledModule.ZeroGradientStageBooksDepositsWithoutMovingThem; regression.b6_subcycled (the surface identity now closes with the booked term) |
| Periodic topology `domain.periodic` × surface/wind (Stage-3 T3) | regression.s3.pd1 (discrete wind–friction balance, strict + default-solver record + transposed periodic-y lane); regression.s3.pd4 (seam books no boundary flux); GridSerial.PeriodicWrapAndEdgeTopology |
| Periodic × transport (seam advection, superbee far neighbors) | regression.s3.pd2 (roll equivalence, upwind implicit via the slab + superbee scheme, seam-crossing anti-vacuity); regression.s3.pd-restart (the project restart standard: flow fields 1e-12 abs, tracer at its measured class — S3-A6) |
| Periodic × groundwater (matrix wrap, side codes, edge faces, baroclinic) | regression.s3.pd2 (gw lane); regression.s3.pd5 (confined Darcy, zero x flux, column-K invariance) |
| Periodic × coupled (seepage across the seam, external closure, consistent salinity both ways, terrain-following) | regression.s3.pd2 (pd2-coupled roll lane + `coupled_exchange_floor` anti-vacuity); regression.s3.pd3 (coupled 1/2/4 + 2 × 2 decompositions) |
| Periodic × MPI decompositions (self-wrap / same-rank-both-sides / wrapped ranks) | regression.s3.pd3 (n = 1/2/4, coupled also 2 × 2); mpi.halo.n4 (the periodic wrap case runs at every n) |
| Periodic axis-generality (x ↔ y under rotation, owner decision O3) | regression.s3.pd8 (gw + surface dihedral batteries with the W3.7 flag swap); pd1 transposed lane |
| Periodic schema/setup rules (§13.B) | regression.s3.pd6 (x-edge BC fatal, seam warning in log + record); ConfigTest (nx >= 3) |
| Interval-mean rate outputs `seepage_rate` (Stage-3 T4, IntervalRates registry) | regression.s3.dg1 (rate-ledger identity + the V3 snapshot-vs-rate anti-vacuity); regression.s3.dg-restart (acc_* checkpoint round trip, bitwise) |
| Monitor/audit row decimation `output.monitor_interval` (Stage-3 T4) | regression.s3.dg2 (bitwise row subset at 60 s; forced rows at checkpoint/t_end on the 77 s lane); regression.s3.default_bitwise (key absent = every step, bitwise) |
| ET open-water channel (Stage-3 T5 W5.3: wet-fraction clock, per-cell removal) | regression.s3.et3 (pond phase closed-form + the switch); regression.s3.et9 (tidal zones); AtmEt.BulkPotentialMatchesV20EvaporationRateBitwise |
| ET soil channel (W5.4/W5.5: dry-fraction split, supply cap, guard-free D1 fix) | regression.s3.et1 (v2.0-chain agreement + the coupled D1 red); regression.s3.et2 (coupled dry == uncoupled); regression.s3.et4 (water balance); regression.s3.et7 (cap books) |
| ET limiters and r_a laws (alpha_barton / beta_lee_pielke / none; liu2006 / fao56 / log_profile) | AtmEt.AeroResistanceLaws; AtmEt.LeePielkeBeta; AtmEt.ThetaFcDefaultIsRetentionAt330cm; unit.et_gates.negative (two-way chain check) |
| Time-basis complementarity (no double count across channels) | regression.s3.et3 (wet+soil <= t; the v2.0 window-lag real negative); regression.s3.et9 (3 tidal cycles, per cell) |
| Surface salt store `transport.surface_dryout: store` (O4, D6) | regression.s3.et8 (b: capture at dry-out + rewet release; c: the v2.0 delete loss as the real negative); salt_store output + surf_store audit column |
| Thin-layer removal books P1–P4 + `film_removed`/`film_pending`/`stage_volume` audit columns (Stage-3 T9 W9.1) | SweModule.ThinLayerRemovalsBookedPerPath; regression.s3.ws1 (P1 recorded, P4 thin vents floor; the discard twin's identity at every row + the film_removed-dropped negative); regression.s3.ws4 (P2, P3 per-lane books); regression.s3.default_bitwise (mass_audit T8 prefix bitwise) |
| Paired water–salt store × no ET, coupled seepage face (S3-A12 / O15) | regression.s3.ws1 (source bound 35 psu, I2, ws2 ledgers; red 494 psu on T8); regression.s3.ws1n (the owner's H2 repro; red 72 854 psu) |
| Paired store × legacy evaporation | regression.s3.ws3 (d: Σ(E − clamped) + V_s = V0, M_s = M0, V_s monotone while dry); regression.s3.ws4 (P3 lanes) |
| Paired store × ET open water | regression.s3.ws3 (a-ET bitwise vs T8; b: Σ E + V_s = V0); regression.s3.et8, et9, et9-veg, et-orient (W9.6 re-measured) |
| Paired store × rain | regression.s3.ws3 (a: store ≡ discard bitwise; c: the rewet closed form); SweModule.StoreReturnsAtWettingFaceDepth |
| Paired store × surface-only | regression.s3.ws3, ws4; TransportStore.* (unit) |
| Paired store × periodic | regression.s3.et9 (periodic-x tidal flat, explicit store); regression.s3.pd2 (pd2-coupled resolves store) |
| Paired store × temperature (R9) | regression.s3.ws4 (P2T/P3T: parked heat = T·V_s); TransportStore.TemperatureParksWithItsWater |
| Paired store × restart | regression.s3.ws-restart (checkpoint with ≥ 10 parked cells; water_store/salt_store/fields bitwise) |
| Paired store × rank | regression.s3.ws_rank (stable slope base: S3-A6 classes + Σ stores; ws1 at n = 2/4: source bound, I2, ledgers — S3-A13) |
| Subsurface ET salt (W5.8 exact top factor; the cauchy-style zero scalar flux) | regression.s3.et8 (a: per-cell c·theta identity, adjust/anchor at rounding) |
| Held-seepage salt (W5.9, D9) | regression.s3.et8 (d: full-identity closure, surf_adjust == 0 exactly, store > floor while held) |
| ET x coupled / external stage / periodic / salinity (the everything smoke) | regression.s3.et9 (+ the noET lane BITWISE vs its T5-start pin — T5 default-neutrality by execution) |
| Vegetation partition (Beer–Lambert tau on the ground channels; T_pot = (1 − tau)·E_x(r_c)) | regression.s3.et5 (dense canopy: E_soil < 1e-12·T, the bare contrast negative, the ALL-ZERO-map bitwise pin); AtmEt.BeerLambertLimits; AtmEt.CanopyResistanceAndRootIntegral |
| Feddes stress + the transpiration supply cap (O7 wetland form) | regression.s3.et7-t (below-h4 shutoff, theta floor, books exact); AtmEt.FeddesPiecewise |
| Root-uptake sink in the Richards step (predictor/corrector shared sign) | GwModule.SinkInSaturatedConfinedCellDropsHeadBySinkOverSs; regression.s3.et4-t (the extended closure identity) |
| Root fractions on stretched/masked/terrain columns + the root_uptake Field3 rate (S3-A7's leg, first consumer) | regression.s3.et6 (per-cell exactness, uniform + linear); GwEt.RootFractionsIntegrateToOneOnStretchedMaskedColumns |
| Transpiration salt exclusion (water out, salt stays; exact per-cell factor) | regression.s3.et8-t (c·theta identity + the rootFrac correlation) |
| Vegetation x coupled tidal / O7 inundated transpiration / restart | regression.s3.et9-veg (the marsh band; T during inundated instants); regression.s3.et-veg-restart (bitwise incl. root_uptake and the store) |
| Penman–Monteith potential (Stage-3 T7: q-form, analytic dq_sat/dT, per-channel albedo/G) | regression.s3.et10 (T_act = pm_q_form to 1e-10; FAO-56 Example-19 within delta_conv; the albedo negative); AtmEt.PmEquilibriumLimit; AtmEt.PmAerodynamicLimit; AtmEt.PmHighResistanceLimit; AtmEt.DqsatDtMatchesCentralDifference; unit.et_gates.negative (the FAO-56 transcription lanes) |
| PM x bulk commensurability (the energy-balance surface temperature) | AtmEt.PmCommensurateWithBulkAtEnergyBalanceTemperature (absolute-gap order 2.03, S3-A10) |
| PM x temperature module / PM x alpha_barton / PM without radiation keys | schema-rejected: ConfigTest.RejectsPenmanMonteithWithTemperatureModule; ConfigTest.RejectsPenmanMonteithWithoutRadiationPair; the alpha_barton rule (stage3/records/T7/red/pm-potential-validate.log) |
| PM x soil / vegetation albedo declarations | schema-rejected when absent: ConfigTest.RejectsPenmanMonteithWithoutTypeAlbedos |
| ET x orientation / rank / restart | regression.s3.et-orient (dihedral battery); regression.s3.et-rank (n = 1/2/4); regression.s3.et-restart (et3 mid-transition + et8b with a full store, bitwise) |
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
| `surface_update: consistent` + MPI decomposition | regression.tracer_rank_invariance.strict (superbee-consistent lane) |
| `subsurface_update: consistent` + subcycled coupling + density + MPI decomposition | regression.coupled_salt_rank_invariance.strict; regression.coupled_salt_rank_invariance.default |
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
run-record+anything (observer only, r1 asserts non-perturbation),
`surface_update`+restart (the consistent update reads only the current
step's flow rates and carries no state across the step boundary; the
checkpointed legacy flow-rate snapshots stay in the layout, unused),
`subsurface_update`+restart (the window volumes are consumed by every
transport step, so they are empty at each checkpoint, and the mass basis
is the restored θ; measured on x-orient/coupled-salt-base: restart
reproduces the straight run to rounding, 5e-16 in the flow, as the legacy
update does on this feature mix).

## 3. GPU / device rows (§8.3 release rule)

| Row | Status | Evidence |
|---|---|---|
| CUDA compile + link (p4) | experimental | .github/workflows/cuda-compile.yml (every library, frehg, and the test binaries build and link warning-free with nvcc; nothing executes) |
| Device memory-space discipline (p5) | supported | unit.p5.backend_invariant_fires; check_forbidden.sh static bans |
| Kokkos LA backend on CPU (aijkokkos rehearsal, p1) | supported | regression.p1.b1.omp4.aijkokkos; regression.p1.b6.ss.omp4.aijkokkos |
| GPU execution (physics on device) | experimental | pending the owner-run p6 bundle (scripts/gpu_acceptance.sh, docs/developer-guide/gpu-acceptance.md); the startup banner prints the experimental status |
| GPU-aware MPI (`runtime.gpu_aware_mpi`) | experimental | exercised only by scripts/gpu_acceptance.sh item 4 (owner-run) |
| GPU performance | experimental | out of v2.0 scope (plan §11); recorded, never asserted |
