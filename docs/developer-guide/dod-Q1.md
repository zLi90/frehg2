# Definition of Done — Q1 (Scalable preconditioner)

Per v2 development plan §2 (Q1) and §9, under the v1 plan's §11.4 DoD
conventions. Every item names the command that verifies it; commands were
run green on 2026-09-09/10 on macOS arm64 (gcc-15, Apple clang 15, MPICH
4.3, PETSc 3.25.1 **with hypre 3.1.0**, Kokkos 5.1.1, deps at
`/Users/zhili/Codes/local`; regression runs pin `FI_PROVIDER=tcp`).
Amendments V2-A1…V2-A4 are logged at the end of
`FREHG2_V2_DEVELOPMENT_PLAN.md` in this commit.

## Deliverables (v2 plan §2.2) — all implemented and tested

- [x] `solver.{surface,groundwater}` YAML block: per-system
      `preconditioner` (`bjacobi-icc` | `amg` | `gamg`), `rtol`, `atol`,
      `max_iterations`, `reuse_max_solves`, `reuse_iteration_factor`;
      strict schema validation (enum + unknown-key rejection) and docs in
      lockstep — `build/tests/frehg_unit_tests
      --gtest_filter='ConfigTest.Solver*'` and
      `python3 scripts/check_parameter_docs.py` (144 tokens documented).
- [x] Defaults reproduce v1 exactly (bjacobi-icc; fs 1e-8/1e-14/500, gw
      cap 1000) — `--gtest_filter='ConfigTest.SolverDefaultsReproduceV1'`;
      plain b1/b2 gates re-run green after the change
      (`ctest -R 'regression.b1$|regression.b2$'`).
- [x] hypre BoomerAMG path (`amg`): HMIS + ext+i + P_max 4; strength
      threshold 0.5 and one aggressive-coarsening level for the 3D `gw_`
      system, 0.25 / none for the 2D `fs_` system; plain-message fatal
      when PETSc lacks hypre; every default injected only-if-absent so
      the options file / command line always wins —
      `--gtest_filter='LinearSystem.PoissonSolvesUnderHypreAmg'`.
- [x] PETSc GAMG path (`gamg`): agg, nsmooths 1, threshold 0.02 —
      `--gtest_filter='LinearSystem.PoissonSolvesUnderGamg'`; unknown
      preconditioners are fatal
      (`--gtest_filter='LinearSystem.UnknownPreconditionerIsFatal'`).
- [x] Hierarchy-reuse policy (KSPSetReusePreconditioner): rebuild on the
      `reuse_max_solves` cadence, on iteration growth past
      `reuse_iteration_factor` × the post-rebuild count + 2, on
      `forceRebuild()`, and rebuild-and-retry-once on stale-hierarchy
      divergence (V2-A2); bjacobi-icc keeps the v1 rebuild-every-solve
      semantics — `--gtest_filter='LinearSystem.HierarchyReuse*:
      LinearSystem.DefaultPreconditionerRebuildsEverySolve'`.
- [x] Wet/dry-mask rebuild trigger: Allreduced mask checksum in
      `SurfaceSolver::fillAndSolve` forces a rebuild when the A13
      dry-closure pattern changes (collective by construction) —
      exercised by every coupled g1 lane (mask changes as b5 wets).
- [x] Dry-row audit (v2 plan §2.2.5), recorded in `report-Q1.md`:
      symmetric Dirichlet column elimination and the compressed gw
      pattern already satisfied the requirements; the one finding —
      literal-1.0 Dirichlet diagonals among cell-area-scale wet rows —
      fixed by scaling the prescribed-stage rows by the cell area
      (row solution unchanged; b1 golden gate re-verified).
- [x] Solver telemetry: per-system `SolverTelemetry` (solves, iteration
      mean/max, rebuilds, retries, timed KSPSetUp/KSPSolve split) and the
      rank-0 `solver summary fs|gw:` end-of-run lines — the g2/s-gate
      parsing contract (`scripts/run_scaling.py` SOLVER_LINE).
- [x] hypre in every build path: `scripts/ci_install_deps.sh` and
      `build_frehg2_hpc.sh` configure `--download-hypre`; CI dependency
      cache keys bumped (`…-petsc3.25.1-hypre-v2`); local PETSc rebuilt
      with hypre 3.1.0 (`grep PETSC_HAVE_HYPRE
      /Users/zhili/Codes/local/include/petscconf.h`).

## Gates (v2 plan §2.3, §9 Q1 row)

- [x] **g1 per-PR (10 lanes):** b1–b4 under `amg` and `gamg` pass their
      unchanged golden/metric criteria; default-mode rank-invariance
      b1/b2 under `amg` at 1/2/4 ranks —
      `ctest -R regression.g1 -LE regression_nightly` (100 % passed,
      2026-09-09).
- [x] **g2:** mean gw CG iterations 8.47 (1 rank) → 9.00 (8 ranks) on
      the A23 synthetic case under `amg`: ratio 1.06 ≤ 1.10; recorded
      count seeded in `tests/regression/tolerances/g2.yaml` (V2-A1: the
      case's fs system converges in 0 iterations, so the recorded gate
      binds gw; fs is exercised by g1) —
      `ctest -R scaling.g2.iters.amg`.
- [x] **g1 nightly:** b5 rain/sync under `amg` at 4 ranks on the parent's
      A17-approved 24 h record horizon (V2-A4), b6 ss and td under `amg`,
      and default-mode b5 rank invariance under `amg` with the
      solver-keyed data-derived bounds (V2-A4) —
      `ctest -R regression.g1 -L regression_nightly` (4/4 passed,
      2026-09-10: b5 discharge 0/48 + ponding 0/54 outside, integrals
      +2.5 %/−2.5 %; b6 ss MAE 0.0334 m / td 0.0390 m, salt mass
      +5.5 %/+8.7 %).
- [x] **s-gate harness live (§7.1/§7.2):** `scripts/run_scaling.py`
      `--gate strong|weak|threads|hybrid|iters|perf-baseline` with JSON
      artifacts; nightly ctest entries s1–s4 + the g2 bjacobi record
      armed under label `scaling_nightly` and run green locally —
      `ctest -L scaling_nightly -j 1` (results in `report-Q1.md` §4).
- [x] **Performance baseline seeded (§7.1.2):**
      `docs/developer-guide/perf-baseline.json` written by
      `run_scaling.py --gate perf-baseline --seed-baseline` (min over 3
      repeats) and the per-PR check passes —
      `ctest -R regression.perf_baseline`.

## Standing criteria (v1 plan §11.2, unchanged)

- [x] Zero-warning build, gcc-15 `-Werror` full flag set —
      `cmake --build build` (clean).
- [x] Unit + MPI labels green — `ctest -L unit` (11/11 entries incl. the
      schema↔docs lockstep and benchmark validation), `ctest -L mpi`
      (all green with the regression overlaps).
- [x] Per-PR regression label green including the new lanes (g1 per-PR,
      g2, perf-baseline) in one sweep — `ctest -L '^regression$'`
      (2026-09-10 record in `report-Q1.md`).
- [x] Forbidden-pattern scan clean — `scripts/check_forbidden.sh`.
- [x] Parameter-docs lockstep — `python3 scripts/check_parameter_docs.py`.
- [x] mkdocs strict build — `mkdocs build --strict`.
- [x] Sanitizer lane (Apple clang ASan+UBSan, deps-clang prefixes,
      halt_on_error) green over the full per-PR labels — 33/33 tests
      (unit 11, mpi 14, regression 8 incl. the restart and
      rank-invariance lanes) in 2840 s, `run_sanitizers.sh: clean`,
      2026-09-10. The amg/gamg PC paths are exercised by the
      LinearSystem unit tests inside the lane; hypre itself is
      uninstrumented (system dependency), so coverage is the
      frehg2/PETSc boundary.
