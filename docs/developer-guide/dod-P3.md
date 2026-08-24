# Definition of Done — P3 (Surface–subsurface coupling)

Per upgrade plan §10 P3 and §11.4. Every item names the command that
verifies it; all commands were run green on 2026-08-11 on macOS arm64
(gcc-15.2, Apple clang 15, MPICH 4.3, PETSc 3.25.1, Kokkos 5.1.1, HDF5
1.14.5 parallel, deps at `/Users/zhili/Codes/local`). Plan amendments
A9–A17 (the exchange limit and subcycled accumulator, the sync common
adaptive step, coupled rain and the coupler-owned top, the uniform terrain
slab, the coupled-regime corrections adjudicated at the b5 gate, the
b5 rank-invariance lanes, the owner-directed shortened gate horizons, the
gate-adjudication follow-ups — bitwise coupled restart, output-label fix,
the audited below-bed clamp, the b5 free outfall — and the owner-directed
early conclusion of the b5 envelope battery) are logged at the end of
`FREHG2_UPGRADE_PLAN.md` in this commit.

## Deliverables (plan §10 P3) — all implemented and tested

- [x] `frehg::coupling` STATIC library (`src/coupling/`: Coupler,
      Exchange) wired into the driver's coupled loop — sync marches both
      modules on the common adaptive step (amendment A10), subcycled keeps
      the fixed surface dt with the legacy while-a-dtg-fits window
      (solve.c:72-75) — `cmake --build build` (zero warnings, `-Werror`).
- [x] `Exchange.cpp` / the coupled top branches: wet-cell Dirichlet head =
      surface depth with the saturated face conductivity (and the
      preserved predictor/corrector face-K asymmetry), infiltration
      limited to the available surface water per the §5.7 resolution
      (amendment A9 — no porosity factor; per-window budget), seepage → η
      per §5.7, the dry seepage face, the saturation bounce-back, the
      reallocation vent, the coupled evaporation correction, and the
      supply-limited classification (amendment A13; SERGHEI
      GwBC.h:277-288) —
      `build/tests/frehg_unit_tests --gtest_filter='CoupledModule.*'`
      (8 tests), provenance table in `docs/theory/exchange-flux.md`.
- [x] Coupled mass audit: the exchange accumulates through
      `/monitor/gw_mass_audit` + the mass_audit `seepage` column; the
      closure identity (every term measured: exchanged, bounce, vent,
      discarded, evaporated, in-transit) is asserted per step —
      `--gtest_filter='CoupledModule.PondedInfiltrationClosesMassBudget:CoupledModule.SubcycledWindowClosesMassBudget'`.
- [x] Checkpoint/restart on coupled runs: SWE {eta, etan, uu, vv, cflx,
      cfly} + GW {h, wc} + the coupler state {seep_accum, qss} + scalars
      {dtg, tgw_lag, dt_surface} — `ctest --test-dir build -R
      regression.b5_restart` (bitwise; etan and dt_surface exist so the
      refresh can reconstruct the stage-boundary correction's west/south
      edge-slot velocities, which live in halo slots the interior-only
      checkpoint cannot carry but feed the uy/vx interpolation).
- [x] `domain.terrain_layers: uniform` (amendment A12) with schema
      cross-check, docs, and mesh test —
      `--gtest_filter='GwMesh.UniformTerrainLayersStackBelowEachBed:ConfigTest.*TerrainLayers*'`.
- [x] Schema cross-check: `groundwater_top` conditions must be kind flux
      in coupled runs (amendment A11) —
      `--gtest_filter='ConfigTest.RejectsCoupledGroundwaterTopHead'`.
- [x] b5 configurations for both scenarios
      (`benchmarks/b5-vcatchment/b5-vcatchment{,-norain}.yaml`), the
      envelope gate, tolerances, and README adjudications —
      `ctest --test-dir build -R 'validate.b5|regression.b5'`.

## Exit criteria (plan §10 P3)

- [x] **b5 envelope gate — rain/sync PASS; battery concluded early by
      the owner (amendment A17)** on the A15 shortened horizon
      (`tests/regression/run_regression.py b5 --t-end 86400`, 4 ranks):
      discharge 0/48 reference times outside the envelope, peak +2.6 %
      (allowed 15 %), integral +2.4 % (allowed 10 %); ponding 0/54
      outside, peak −1.5 %, integral −2.5 %; surface budget residual
      −10 m³ of 1.11e4 (bound 1e-3) with the below-bed clamp measured at
      1.6 % of rain. rain/subcycled measured FAIL at the committed
      dt = 5 s (peak +16.3 %; diagnosis and the ~2 s recommendation in
      the config and A17); norain not gated. The full-horizon runs remain
      as `ctest --test-dir build -R 'regression.b5\.'`
      (regression_nightly label).
- [x] **Coupled mass closure < 1e-10/step** —
      `--gtest_filter='CoupledModule.*ClosesMassBudget'` (the §8.1
      three-column ponded-infiltration toy, sync and subcycled; every
      audit term measured, residual at rounding).
- [x] **b1/b2 regressions still pass with unchanged tolerance files** —
      `ctest --test-dir build -R 'regression.b1$|regression.b2$'`
      (b1 worst achieved/allowed 0.231 — identical to the P1 record —
      and mass closure 1.6e-11 of rain; b2 unchanged; b3/b4 also re-run
      green, b4 rel-L2 0.082, peak +0.3 %).
- [x] **Restart-at-half ≡ uninterrupted to 1e-12** —
      `ctest --test-dir build -R regression.b5_restart` (coupled, 4
      ranks, bitwise — every output label of the restarted leg agrees
      with the uninterrupted run to the last bit after the edge-slot
      reconstruction fix; `regression.b1_restart` and
      `regression.b2_restart` also bitwise).
- [x] **Rank invariance on b5-short** —
      `ctest --test-dir build -R regression.b5_rank_invariance`
      (amendment A14): strict one-step field agreement 5.7e-14 at
      1/2/4 ranks (bound 1e-12); default 600 s bulk volumes — subsurface
      1.6e-6/1.1e-5 (bound 3e-5), exchanged ≤ 4.9e-4, surface ≤ 1.7e-5
      of the exchange (bounds 2e-3) — with the per-field trajectory
      difference recorded (chaotic, up to ~1.2; the measured controller
      fork is documented in the amendment).
- [x] **P0/P1/P2 criteria still green** — `scripts/ci_build_and_test.sh`
      (zero-warning gcc + clang builds, all unit/mpi labels, b1–b4 gates,
      forbidden scan, parameter docs, Doxygen exit 0), ASan/UBSan lane
      clean (`scripts/run_sanitizers.sh`, clang lane per report-P0 §4 —
      no findings; the coupled 600 s rank-invariance runs execute under
      both sanitizers).

## Additional gates (plan §11)

- [x] Fidelity provenance table for the coupling —
      `docs/theory/exchange-flux.md` (legacy `file:line` → Frehg2
      function for every exchange piece, the §5.7 derivation, the
      supply-limited classification, and the coupled audit identity);
      `docs/theory/surface-water.md` gains the three P3 defect rows
      (dry-cell rows, dry-cell drag, four-edge audit) and
      `docs/theory/groundwater.md` the coupled-top cross-reference.
- [x] Newly-found dead code recorded — the `rss`/`r_async` monitor row in
      `docs/theory/removed-features.md`.
- [x] Parameter docs in lockstep — `python3 scripts/check_parameter_docs.py`
      (`domain.terrain_layers`, the coupling.mode semantics, time.dt).
- [x] Forbidden scan clean — `scripts/check_forbidden.sh`.

## Known facts recorded for later phases

- The b5 runs pin the coupled-regime corrections (amendment A13): the
  supply-limited exchange classification, the zero-depth continuity
  closure for dry cells, the min_depth-floored drag, and the four-edge
  audit. b6 (P4) is always capacity-limited, so its goldens see only the
  legacy paths — but its tidal waterline moves, so the A9 limit and the
  A13 closure both run there; the loose b6 tolerances absorb them.
- Coupled trajectories are chaotic at threshold density (amendment A14):
  any per-cell comparison across configurations must use bulk
  observables or expect O(1) scatter. P4's b6 comparisons are at
  quasi-steady state, which is exactly the robust regime.
- The `seepage` output/checkpoint field is the surface-applied rate
  [m/s] of the last coupled step (the legacy qss observable).
- In sync mode the marching dt is coupler-owned (`Coupler::nextDt`);
  time.dt is its initial value only (amendment A10). The subcycled
  subsurface clock lag is checkpointed (`tgw_lag`).
