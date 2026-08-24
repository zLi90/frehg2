# Definition of Done — P1 (Surface-water module)

Per upgrade plan §10 P1 and §11.4. Every item names the command that
verifies it; all commands were run green on 2026-07-25 on macOS arm64
(gcc-15.2, Apple clang 15, MPICH 4.1, PETSc 3.25.1, Kokkos 5.1.1, HDF5
1.14.5 parallel, deps at `/Users/zhili/Codes/local`). Plan amendments
A1–A5 (corner exchange, `kind: outflow`, `rainfall.exclude`,
waterfall-location removal, default-lane rank-invariance bound) are logged
at the end of `FREHG2_UPGRADE_PLAN.md` in this commit.

## Deliverables (plan §10 P1) — all implemented and tested

- [x] `frehg::swe` STATIC library (`src/swe/`: SurfaceSolver, Momentum,
      FreeSurface, WetDry, SurfaceSources + SweFormulas.hpp) and
      `frehg::driver` (`src/driver/Simulation`) with the SWE-only time
      loop; `frehg <config.yaml>` runs a simulation —
      `cmake --build build` (zero warnings, `-Werror` set).
- [x] Momentum source: upwind advection with the exact 0.5–0.7 CFL damping
      ramp (`shallowwater.c:166-169`), central eddy viscosity with the
      legacy face-area asymmetry, point-implicit Manning drag with the `hD`
      exponent switch, **Chezy law** (plan §5.8), wind stress with
      thin-layer attenuation —
      `build/tests/frehg_unit_tests --gtest_filter='SweFormulas.*'`
      (hand-computed values) plus the b1 field gate below.
- [x] Implicit free surface via `LinearSystem("fs_")`: one global SPD
      5-point system, COO fill, prescribed-stage rows with column
      elimination, dry-cell regularization, legacy elevation offset —
      `ctest --test-dir build -R regression.b1`.
- [x] Wetting/drying: min-depth, one-cell wetting limiter,
      higher-of-two-bottoms face depths, face-vanishing and dry-cell
      velocity limiters —
      `build/tests/frehg_unit_tests --gtest_filter='SweModule.*'`
      (dam-blocking, evaporation-clamp, offset tests) + b1 gate.
      `waterfall_location` was dropped as dead (amendment A4; row in
      `docs/theory/removed-features.md`).
- [x] Rain/evaporation sources with the configured exclusion region
      (amendment A3) —
      `build/tests/frehg_unit_tests --gtest_filter='SweModule.Rain*:SweModule.Evap*'`.
- [x] Surface BC kinds: `eta` (stage rows + mass-consistent edge-face
      velocity), `discharge` (global-count division), `velocity`
      (edge-face prescription), `outflow` (amendment A2);
      `scalar_value` is rasterized by the P0 BoundarySet and consumed at
      P4 — `build/tests/frehg_unit_tests
      --gtest_filter='SweModule.EtaCondition*:SweModule.Discharge*:SweModule.Velocity*:SweModule.Outflow*'`.
- [x] Monitors + `/monitor/mass_audit` volume-budget table; HDF5 surface
      output per the §7 contract — b1/b4 gates read them;
      `ctest --test-dir build -R 'unit.all'` (Monitor suite).
- [x] SWE checkpoint/restart (prognostic eta, uu, vv, cflx, cfly) —
      `ctest --test-dir build -R regression.b1_restart`
      (restart at t = 9000 s matches the uninterrupted run to 1e-12).
- [x] `HaloExchanger::exchangeWithCorners` (amendment A1) —
      `ctest --test-dir build -R mpi.halo` (corner values bitwise at
      1/2/4 ranks, both staging modes).
- [x] Schema additions in lockstep with docs and tests
      (`rainfall.exclude`, `kind: outflow`, seepage cross-check) —
      `ctest --test-dir build -R 'unit.all|unit.parameter_docs|validate.'`
      (17 + 6 crafted-invalid cases).
- [x] Regression harness (plan §8.3): `tests/regression/run_regression.py`,
      `compare_h5.py`, `tools/ascii_golden_to_h5.py`, per-case tolerance
      files; goldens fetched from `../legacy/benchmarks`, never committed —
      `ctest --test-dir build -L regression`.

## Exit criteria (plan §10 P1)

- [x] **b1 gate passes** — `ctest --test-dir build -R regression.b1`:
      worst achieved/allowed ratios — eta 0.063, depth 0.231, vv 0.131,
      uu 0.000 (identically zero, as golden) against |Δ| ≤ max(5e-4 m,
      1 %) for eta/depth and max(1e-3 m/s, 2 %) for velocities; cumulative
      mass-balance error 2.7e-11 of rainfall volume (allowed 1e-3). Golden
      note: legacy's output trigger fires one step early (solve.c:121),
      absorbed within tolerance (`tests/regression/tolerances/b1-sw.yaml`).
- [x] **b4 gate passes** — `ctest --test-dir build -R regression.b4`:
      hydrograph rel-L2 0.082 (≤ 0.15); peak +0.3 % (±10 %);
      time-to-99 %-peak −7.9 % (±10 %); outflow volume −0.10 % of rainfall
      (≤ 5 %) and −0.22 % of the reference volume (≤ 3 %). Metric
      definitions and the outflow-kind adjudication:
      `benchmarks/b4-govindaraju/README.md`.
- [x] **b1 rank invariance** —
      `ctest --test-dir build -R regression.b1_rank_invariance`:
      strict mode (`-fs_pc_type jacobi -fs_ksp_rtol 1e-13 -fs_ksp_atol
      1e-16`) max relative field difference 1/2/4 ranks ≤ 1e-12 (measured
      ~1e-13 — the rank-invariance proof); default mode (bjacobi/icc at
      production tolerances) measured 5.4e-5 against the amendment-A5
      bound of 1e-4 — plan §8.2's pre-data 1e-6 was amended from this
      measurement (preconditioner drift echoed through wet/dry threshold
      crossings over 3600 steps).
- [x] **Single-rank b1 wall time ≤ 2× legacy serial** — measured 0.29 s
      (`OMP_NUM_THREADS=1`) vs 1.11 s legacy serial (-O3 numerics; parser
      files at -O0 because of the legacy stack-buffer defect,
      utility.c:41-59) — 0.26× legacy. Timer report archived in
      `report-P1.md` §6.
- [x] **P0 criteria still green** —
      `ctest --test-dir build -L unit --output-on-failure` (16/16 suites,
      97 gtest cases), `ctest --test-dir build -L mpi` (6/6),
      `scripts/check_forbidden.sh`, `python3
      scripts/check_parameter_docs.py`, `doxygen docs/Doxyfile` (exit 0),
      zero-warning gcc-15 and Apple clang builds, ASan/UBSan clean
      (clang lane) — commands as in `dod-P0.md`.

## Additional gates (plan §11)

- [x] Fidelity provenance table for every preserved SWE equation —
      `docs/theory/surface-water.md` (legacy `file:line` → Frehg2
      function), including the preserved load-bearing quirks (elevation
      offset, double drag factor, pre-source ghost stage) and the four
      correctness fixes of the §2.1 hazard class.
- [x] Newly-found dead code recorded — three new rows in
      `docs/theory/removed-features.md` (waterfall_location flags,
      `Sct == 0` fallback, `rain_sum`).
- [x] Parameter docs in lockstep — `python3 scripts/check_parameter_docs.py`.

## Known facts recorded for later phases

- Benchmark-scale grids are kernel-launch-latency-bound: `OMP_NUM_THREADS=1`
  is the fastest CPU configuration (b1: 0.29 s vs 16 s at 2 threads); the
  regression harness pins it. Revisit at P5 scaling on b5.
- The b5 provisional stage-sink outlet (−10 m) should be re-examined
  against `kind: outflow` at the P3 gate (amendment A2 rationale).
- The mass-audit table is the mass-consistent discharge source; stored
  velocities carry the legacy double drag factor
  (docs/theory/surface-water.md) — P3/P4 gates that measure discharge
  should use the audit, as b4 does.
