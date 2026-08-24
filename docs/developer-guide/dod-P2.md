# Definition of Done — P2 (Groundwater module)

Per upgrade plan §10 P2 and §11.4. Every item names the command that
verifies it; all commands were run green on 2026-07-26 on macOS arm64
(gcc-15.2, Apple clang 15, MPICH 4.3, PETSc 3.25.1, Kokkos 5.1.1, HDF5
1.14.5 parallel, deps at `/Users/zhili/Codes/local`). Plan amendments
A6–A8 (per-cell zcell + adaptive-run checkpoint keys, the post-allocation
resolution with `groundwater.reallocation_surplus`, the b2 rank-invariance
variant) are logged at the end of `FREHG2_UPGRADE_PLAN.md` in this commit.

## Deliverables (plan §10 P2) — all implemented and tested

- [x] `frehg::gw` STATIC library (`src/gw/`: RichardsSolver, Predictor,
      Corrector, Reallocate, AdaptiveStep, TerrainMetric + VanGenuchten.hpp)
      wired into the driver's groundwater-only time loop (adaptive dtg as
      the outer step, legacy solve.c:43/groundwater.c:193 semantics) —
      `cmake --build build` (zero warnings, `-Werror` set).
- [x] Van Genuchten kernels with the live `aev` cutoff (θ(h) only; C and K
      cut off at h > 0) — `build/tests/frehg_unit_tests
      --gtest_filter='VanGenuchten.*'` (b2 Warrick soil and b6 sand at
      h ∈ {−10, −1, −0.02, −0.01, 0, +0.5} m, 1e-12 relative).
- [x] PCA predictor: 7-point `LinearSystem("gw_")` COO system with storage
      C(h) + Ss·θ/θs, legacy face-conductivity rules (incl. the saturated
      prescribed-head top face and the seal-pass ordering) —
      `ctest --test-dir build -R regression.b2` plus
      `--gtest_filter='GwModule.*'`.
- [x] Corrector: darcy_flux-faithful face fluxes (minus-cell Ks on x/y,
      lower-cell Ks on z), θ update by flux divergence with the Ss factor,
      final clamp with vloss accounting — same commands; the audit identity
      ΔV = boundary_in − ss_storage + realloc − vloss closes ≤ 1e-12
      (`GwModule.AuditIdentityClosesWithCompressibleStorage`).
- [x] `Reallocate` faithful to `groundwater.c:959-1424`, always on
      (plan Appendix A), deterministic per-column sweep, with the
      surplus-direction knob `groundwater.reallocation_surplus`
      (amendment A7) — `--gtest_filter='GwModule.Realloc*'` + the b2
      (drop) and b3 (redistribute) gates.
- [x] Adaptive dtg (`groundwater.c:1651-1727` semantics: dq 0.01/0.02
      thresholds, ×1.25/×0.75, dK/dθ Courant cap, [dt_min, dt_max] clamp,
      cross-rank min; top cells excluded) —
      `--gtest_filter='GwModule.AdaptiveStep*'`.
- [x] `use_full3d` (face-centric sealing, matrix and corrector) —
      `--gtest_filter='GwModule.UseFull3d*'`.
- [x] GW BC kinds: top/bottom `head`, `flux`; bottom `flux: gravity` (free
      drainage, with the legacy bctype-index fix); side `head`
      (hydrostatic and constant) and `flux` (positive = into the domain) —
      `--gtest_filter='GwModule.Prescribed*:GwModule.Gravity*:GwModule.Hydro*:GwModule.Side*'`
      + schema cross-checks (`ConfigTest.RejectsGravityValueOffTheBottom`,
      `RejectsHydrostaticValueOffTheSides`).
- [x] Subsurface HDF5 output per the §7 contract (hydraulic_head,
      water_content, per-area qx/qy/qz; per-cell `/groundwater/zcell/0` —
      amendment A6) + `/monitor/gw_mass_audit`; groundwater
      checkpoint/restart ({h, wc} + dtg) —
      `ctest --test-dir build -R 'regression.b2$|regression.b2_restart'`.
- [x] Subsurface mesh: bathymetry-derived ktop with the legacy partial-cell
      rules, terrain-following per-column dz and face metrics —
      `--gtest_filter='GwMesh.*'`; `Grid::buildGlobalIds` now receives the
      real mask in groundwater runs.
- [x] Regression harness extension: `b2`, `b3`, `b2-restart`,
      `rank-invariance-b2`, `b6-gw-smoke` gates; `compare_h5.py` dataset
      groups; `ascii_golden_to_h5.py` 3D converter; tolerance files
      `tolerances/b2-gw.yaml`, `tolerances/b3-kirkland.yaml` —
      `ctest --test-dir build -L regression`.

## Exit criteria (plan §10 P2)

- [x] **vG unit tests at 1e-12** —
      `build/tests/frehg_unit_tests --gtest_filter='VanGenuchten.*'`
      (reference values hand-evaluated in IEEE double with the code's
      operation order; a 40-digit cross-check agrees to each expression's
      double-precision conditioning).
- [x] **b2 gate passes** — `ctest --test-dir build -R regression.b2`:
      worst achieved/allowed ratios — hydraulic_head 0.40 against
      |Δ| ≤ max(1e-3 m, 1 %), water_content 1.8e-4 of the 0.005 absolute
      allowance; wetting front within 2.0 % of Warrick at 11700 s
      (allowed 5 %) at 1.00–1.01× the legacy golden's own error
      (allowed 1.2×).
- [x] **b3 gate passes** — `ctest --test-dir build -R regression.b3`:
      h = 0 contour max distance 0.160 m, RMS 0.073 m; h = −400 contour
      max 0.170 m, RMS 0.088 m (allowed 0.2 / 0.1 m); internal mass
      balance 0.3 % of the injected volume (allowed 0.5 %). Adjudications
      (contour units, Ss = 0, `reallocation_surplus: redistribute`) in
      `benchmarks/b3-kirkland/README.md`.
- [x] **b2 rank invariance** —
      `ctest --test-dir build -R regression.b2_rank_invariance`
      (the 4×4 replication, amendment A8): strict mode
      (`-gw_pc_type jacobi -gw_ksp_rtol 1e-13 -gw_ksp_atol 1e-16`, through
      t = 23400 s) max relative field difference at 1/2/4 ranks 3.3e-14
      (bound 1e-12 — the rank-invariance proof); default bjacobi/icc mode
      over the full horizon 1.4e-3 (bound 5e-3, data-derived like
      amendment A5: discrete saturation-front cell crossings amplify
      rank-layout rounding past t ≈ 23400 s even at machine-precision
      solves — strict mode itself measures 9.4e-4 on the full horizon).
- [x] **Closed-domain GW mass balance ≤ 1e-8 relative per step** —
      `--gtest_filter='GwModule.ClosedColumnConservesMass*'` (measured at
      rounding level; the audit identity is exact by construction).
- [x] **Debug-build θ ∈ [θr, θs] assertion, zero violations on b2/b3** —
      the `finalizeWaterContent` invariant check (active without NDEBUG)
      ran green through both full gates in a Debug build:
      `cmake -B build-debug -DCMAKE_BUILD_TYPE=Debug ...` then
      `run_regression.py b2|b3 --frehg build-debug/src/frehg ...`.
- [x] **Early b6 GW-only smoke comparison run and recorded** —
      `ctest --test-dir build -R regression.b6_gw_smoke` (also
      plan risk register): b6-kuan-ss reduced to groundwater-only
      (surface/transport off, monitors dropped, density coupling off) vs
      the syncV4 golden at t = 36000 s over 1224 cells:
      **|Δh| max 0.0118 m, mean 0.0057 m; |Δθ| max 0.0031, mean 2e-4;
      saturated zone exact in θ, |Δh| ≤ 0.0118 m.** Context: the golden is
      a coupled Newton run with salinity; the gap sits at the P4 head
      tolerance floor (max(0.01 m, 5 %)) before any coupling/transport —
      the PCA-vs-Newton risk to the b6 gate is small.
- [x] **P0/P1 criteria still green** — `scripts/ci_build_and_test.sh`
      (zero-warning gcc + clang builds, all unit/mpi labels, b1/b4 gates,
      forbidden scan, parameter docs, Doxygen exit 0), ASan/UBSan lane
      clean (`scripts/run_sanitizers.sh` on the clang lane per
      report-P0 §4).

## Additional gates (plan §11)

- [x] Fidelity provenance table for every preserved GW equation —
      `docs/theory/groundwater.md` (legacy `file:line` → Frehg2 function),
      including the preserved quirks (live aev, corrector-vs-predictor face
      K, top-cell controller exclusion, the consistency-restore economics,
      the seal-pass ordering for masked columns, the Ss feedback finding)
      and the eight §2.1-class correctness fixes (all in paths no benchmark
      configuration reaches).
- [x] Newly-found dead code recorded — five P2 rows in
      `docs/theory/removed-features.md` (allocate_recv machinery, repeat
      flag + stale check_room, the dead adaptive-step error estimator,
      terrain top-angle tables, Newton helpers).
- [x] Parameter docs in lockstep — `python3 scripts/check_parameter_docs.py`
      (`groundwater.reallocation_surplus` documented with both gates'
      adjudications).

## Known facts recorded for later phases

- The post-allocation audit columns (`realloc`, `realloc_dropped`) are the
  first diagnostic when a coupled (P3/P4) budget drifts; b5/b6 must choose
  their `reallocation_surplus` mode at their gates (amendment A7).
- Nonzero Ss + dry flat-retention soil + the consistency restore is a
  runaway feedback (b3 finding); b6's steep sand bounds it but P4 should
  watch the ss variant.
- A fully saturated, Ss = 0, all-Neumann groundwater system is singular
  (in legacy too): groundwater-only tests need a head anchor or Ss > 0.
- Masked columns (ktop > 0) keep a sealed predictor top face (legacy
  seal-pass order): first exercised by b5 at P3.
