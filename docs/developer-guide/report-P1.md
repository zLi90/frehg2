# Phase 1 development report — handoff for P2+

**Audience:** the developer (or agent session) implementing P2 and later
phases. Read together with `FREHG2_UPGRADE_PLAN.md` (binding spec; P1 added
amendments A1–A5 to its Amendment log), `report-P0.md` (foundation handoff),
and `dod-P1.md` (the verified Definition of Done). This report explains what
P1 built, what was decided and why, and what P2/P3 must know.

**Status:** P1 complete. All plan §10 P1 exit criteria green: b1 and b4
gates pass, b1 rank invariance holds (strict ≤ 1e-12), single-rank b1 runs
at 0.26× legacy serial, and every P0 criterion remains green (21/21 ctest,
zero-warning gcc + clang builds, ASan/UBSan clean, forbidden scan clean,
Doxygen clean).

---

## 1. What P1 added and where

- **`frehg::swe`** (`src/swe/`): the complete legacy θ-scheme SWE module.
  `SurfaceSolver` owns the state fields and orchestrates the two legacy
  phases (`solveFreeSurface` = `solve_shallowwater`, `updateVelocity` =
  `shallowwater_velocity`; the split is where P3's groundwater call slots
  in, mirroring `solve.c:51/:97`). `Momentum.cpp` (advection/viscosity/
  drag/wind/velocity update/uy-vx interpolation), `FreeSurface.cpp`
  (RHS, coefficients, prescribed-stage handling, COO fill + PETSc solve,
  boundary-flux audit), `WetDry.cpp` (depths, geometry, wetting limiter,
  velocity limiters, edge ghosts), `SurfaceSources.cpp` (rain/evaporation),
  `SweFormulas.hpp` (pointwise closures, unit-tested against hand values).
- **`frehg::driver`** (`src/driver/Simulation`): time loop, output/monitor/
  checkpoint mediation, restart, the `/monitor/mass_audit` volume-budget
  table (columns: volume, cumulative rain/evaporation/boundary-outflow/
  bc-inflow, reduced over ranks, one row per step). `frehg <config.yaml>`
  runs; extra CLI args go to PETSc.
- **Core additions:** `HaloExchanger::exchangeWithCorners` (amendment A1 —
  the uy/vx four-point stencil needs one diagonal neighbor; two sequential
  phases, wide j-rows carry settled i-halo columns); `Grid` auto-
  decomposition now respects extents (4 ranks on b1's 1×10 grid choose
  1×4, not MPI_Dims_create's 2×2 which left empty blocks); surface
  `velocity`/`outflow` kinds rasterize onto domain-edge faces in
  `BoundarySet`.
- **Schema (with docs + tests in lockstep):** `surface_water.rainfall.
  exclude.polygon` (A3), BC `kind: outflow` taking no value (A2), seepage-
  requires-groundwater cross-check, rainfall constant/series exclusivity.
- **Regression harness** (plan §8.3): `tests/regression/run_regression.py`
  (gates b1, b4, b1-restart, rank-invariance), `compare_h5.py` (element-wise
  `|Δ| ≤ max(abs_floor, rel·|ref|)`, always prints achieved-vs-allowed),
  `tools/ascii_golden_to_h5.py` (legacy ASCII → §7 HDF5 layout),
  `tests/regression/tolerances/*.yaml`. Goldens come from
  `../legacy/benchmarks` (CMake cache var `FREHG_LEGACY_BENCHMARKS`), never
  committed. ctest label `regression`; the mpi-labeled rank-invariance
  gates run b1 at 1/2/4 ranks in both PC modes.

## 2. Fidelity findings P2+ must not re-learn

1. **The elevation offset is part of the algorithm.** Legacy shifts all
   elevations by −min(bottom) (`initialize.c:142-147`). The dry-cell
   regularization row (`Sct = dx·dy`, RHS `η·dx·dy`) is not shift-invariant
   — wetting of dry cells depends on absolute η. Without the offset, b1's
   outlet row never wets and the fields are wrong by centimeters. Frehg2
   reproduces it (`SurfaceSolver::elevationOffset()`); outputs subtract it.
   P3's coupling (surface depth → subsurface head BC) must stay in the
   internal frame consistently.
2. **Stored velocities ≠ conserved fluxes.** `update_velocity` applies the
   drag factor twice (`shallowwater.c:745-746`, preserved). Discharge
   measured as `uu·depth` under-reads by ~20 % on b4; the mass-audit table
   is the mass-consistent discharge source. Any later gate that measures
   discharge should use the audit (b4's tolerance file documents this).
3. **Boundary faces see the pre-source ghost stage** (legacy call order):
   rain induces a small outward boundary-face velocity each step; a flat
   "closed" basin leaks slowly at east/north edges. Rank-invariant, tested
   (`SweModule.RainAccumulatesExactlyOnAClosedCell` uses a 1×1 domain,
   which is genuinely closed).
4. **Legacy b1 output labels are one step early** (`solve.c:121` trigger):
   golden `<t>` holds the state at `t − dt`. Absorbed by the b1 tolerances
   (worst achieved/allowed 0.23); if a later phase tightens tolerances,
   this shift is the floor.
5. **Legacy defects fixed on port** (§2.1 hazard class; full table in
   `docs/theory/surface-water.md`): uy/vx wrap/OOB indexing, per-rank
   discharge division, the north-edge `>` vs `>=` rank test, uninitialized
   ghost Fu/Fv. None affects b1/b4 goldens.
6. **The b4 outlet decision** (amendment A2): stage-sink at −10 m keeps the
   outlet column dry; plain open boundary retains half the rain. `kind:
   outflow` extrapolates the ghost stage down the continued bed slope.
   **P3 should revisit b5's provisional −10 m outlet against it.**

## 3. Decisions and conventions added in P1

- Surface BC application points: stage cells are enforced in
  `enforceSurfBc` and enter the matrix as identity rows with column
  elimination (SPD preserved; identical solution to legacy's unsymmetric
  form). Discharge/velocity/outflow modify the RHS in `assembleRhs`/
  `applyOutflowCorrections`; every RHS modification is mirrored into the
  audit so budget closure is exact by construction.
- Restart state for SWE = {eta, uu, vv, cflx, cfly} (+ time/step); all
  other fields are recomputed by `refreshDerivedState()`. b1
  run→checkpoint→restart matches uninterrupted to 1e-12
  (`regression.b1_restart`). P3 adds dtg and subsurface state.
- Output filename resolves against the working directory (regression runs
  copy the case dir to scratch); input paths resolve against the config
  file's directory.
- Strict rank-invariance mode = `-fs_pc_type jacobi -fs_ksp_rtol 1e-13
  -fs_ksp_atol 1e-16` on the command line (PETSc consumes argv).
  Measured: strict ≤ ~1e-13 at 1/2/4 ranks; default bjacobi/icc drifts to
  5.4e-5 over b1's 3600 steps → amendment A5 set the default-lane bound to
  1e-4 (1e-6 was pre-data).
- Wind direction series interpolate linearly (legacy held the previous
  sample); documented deviation, no benchmark exercises wind.

## 4. Performance (exit criterion + facts for P5)

b1 single-rank (3600 steps, 1×10 cells), macOS arm64:

| Configuration | Wall time |
|---|---|
| legacy serial (`-O3` numerics, parser at `-O0`*) | 1.11 s |
| frehg2, `OMP_NUM_THREADS=1` | **0.29 s** (0.26× legacy; bound was ≤ 2×) |
| frehg2, `OMP_NUM_THREADS=2` | 16 s |
| frehg2, default thread count | 34 s |

Timer report (OMP_NUM_THREADS=1): simulation 0.199 s = free_surface 0.099 s
(of which PETSc solve 0.019 s) + velocity 0.073 s.

*The legacy makefile is stale (`linsys.c` does not exist; `subroutines.c`
missing); the corrected file list builds, but `-O1` and above break the
`read_one_input` dangling-stack-buffer defect (`utility.c:41-59`) — `Tend`
parses as garbage and the run stops after one step. Parser files compiled
`-O0`, numerics `-O3`.

Benchmark-scale grids are kernel-launch-latency-bound (~35 Kokkos kernels ×
2 fences per step): thread counts > 1 only add OpenMP fork/join cost. The
regression harness pins `OMP_NUM_THREADS=1`. Real parallelism pays at b5
scale and on GPUs — measure there (P5), don't "fix" this now.

## 5. What P2 should do first

1. Re-run `scripts/ci_build_and_test.sh`; re-read plan §10 P2 and §9 b2/b3
   gates.
2. Add `src/gw/` as STATIC `frehg_gw` mirroring `src/CMakeLists.txt`'s
   frehg_swe block; physics libraries never link each other — the driver
   owns cross-module state (plan §4). `Simulation` currently fatals on
   `modules.groundwater: true` — replace that guard with the gw-only loop
   (legacy loops gw-only cases with the surface `dt` as the outer step;
   see report-P0 §3).
3. The 7-point subsurface system: same pattern as `buildCooPattern`/
   `fillAndSolve` in `FreeSurface.cpp` but with `gid3()` and prefix `"gw_"`
   (see plan §5.1 and report-P0 §6). `Grid::buildGlobalIds(ktop)` with the
   real bathymetry-derived ktop replaces the all-zeros mask the driver
   passes today.
4. Reuse the regression harness: add a `b2`/`b3` gate function to
   `run_regression.py` and a tolerances file; `ascii_golden_to_h5.py`
   needs a 3D variant (`/groundwater/<var>/<t>`, `(j·NX+i)·NZ+k`).
5. Remember the early b6 GW-only smoke comparison at P2 exit (plan §10 P2)
   — it sizes the P4 tolerance risk.

## 6. Open items and risks carried forward

- CI workflows remain unexercised (no remote); first push must watch them.
  The regression jobs need `FREHG_LEGACY_BENCHMARKS` pointed at a fetched
  goldens location on runners.
- MPICH on this machine intermittently spins in `MPI_Finalize` (OFI sockets
  provider; observed once, not reproducible, also implicated after an
  `MPI_Abort`). If a test hangs at exit with ~100 % CPU, sample the process
  before killing: the signature is `MPIDI_OFI_mpi_finalize_hook` polling.
- b5's outlet mapping and the subcycled-coupling gate arrive at P3; b5 runs
  will also be the first real (multi-rank, larger-grid) exercise of
  `exchangeWithCorners` on 2D decompositions beyond the halo unit test.
- The `velocity` BC kind has defined, unit-tested semantics (prescribed
  face-normal velocity at domain-edge faces) but no benchmark exercises it;
  if a later case needs distributed (non-edge) velocity forcing, that is a
  schema discussion, not a patch.
