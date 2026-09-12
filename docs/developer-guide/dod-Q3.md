# Definition of Done — Q3 (Performance-portable parallelization)

Per v2 development plan §2B and §9, under the v1 plan's §11.4 DoD
conventions. Every item names the command that verifies it; commands were
run green on 2026-09-11 on macOS arm64 (gcc-15, MPICH 4.3, PETSc 3.25.1
+ hypre 3.1.0 + Kokkos Kernels, Kokkos 5.1.1 OpenMP host backend; all
runs pin `FI_PROVIDER=tcp`). GPU items are compile/statically verified
per §2B.4 — device execution is the owner-run p6 bundle.

## Deliverables (v2 plan §2B.2) — all implemented and tested

- [x] **B1 — device-consistent linear algebra.** `MATAIJ` is no longer
      hardcoded: `solver.<system>.mat_type: aij | aijkokkos` (schema +
      docs + resolved-config round-trip in the same PR, the Q2 lockstep);
      device builds (`MemSpace != HostSpace`) force `aijkokkos`/`kokkos`
      and reject host overrides with a fatal; the resolved type is read
      back from PETSc (`MatGetType`) and recorded. The solve stages
      RHS/solution through `VecGetKokkosView{Write,}` when the vector is
      Kokkos — zero-copy view alias on host, one device-side copy on GPU —
      and through the v1 host mirror otherwise. `petscvec_kokkos.hpp` is
      included first in the TU (PETSc hard requirement). Verified:
      64×64×25 synthetic, aij vs aijkokkos, max abs field diff 1.57e-9
      (inside rank-invariance-class rounding), identical gw iterations.
- [x] **B2 — backend-aware AMG defaults.** `coarsen_type` HMIS (host) /
      PMIS (device — the only GPU-capable coarsening); `ext+i` interp both
      (device-supported); on Kokkos lanes the smoother defaults to
      l1-scaled Jacobi (thread-invariant; hypre's hybrid host default
      drifts with OMP_NUM_THREADS — amendment V2-A7); aggressive-
      coarsening interp left to PETSc ≥ 3.25's device-aware selection.
      All via `setOptionDefault` — options file / command line still win.
      The **resolved** coarsening/smoother are read back from the options
      database post-`KSPSetFromOptions` and written into the run record's
      solver section (`amg_coarsen_type`, `amg_relax_type`).
- [x] **B3 — Kokkos-enabled PETSc in every build path.**
      `build_frehg2_local.sh`, `build_frehg2_hpc.sh`,
      `scripts/ci_install_deps.sh`: `--with-kokkos-dir
      --download-kokkos-kernels --with-openmp` (+ `--with-cuda=1` on the
      FREHG_CUDA lane); local script hard-fails if the installed PETSc
      lacks `PETSC_HAVE_KOKKOS_KERNELS`; CI cache keys bumped
      (`…-kokkoskernels-v3`).
- [x] **B4 — one GPU-aware-MPI switch.** `runtime.gpu_aware_mpi` remains
      the single toggle (halo exchanger); the p6 bundle exercises the
      solve under it on/off at 2+ ranks.
- [x] **B5 — build-system consistency.** `FREHG_BACKEND=auto|serial|
      openmp|cuda|hip` cache variable: requested backend verified against
      the found Kokkos, and a device Kokkos paired with a PETSc lacking
      `PETSC_HAVE_KOKKOS` fails **at configure time** with instructions.
- [x] **B6 — guardrails for absent hardware.** See p4/p5 below.

## Gates (v2 plan §2B.3)

- [x] **p1 — backend invariance (per-PR + nightly).**
      `ctest -R regression.p1` — b1 and b4 at OMP_NUM_THREADS=4 under
      both `aij` and `aijkokkos`: 4/4 PASS, unchanged gate criteria
      (harness: `run_regression.py --mat-type`). Nightly:
      `regression.p1.b6.ss.omp4.aijkokkos` (full coupled+transport+
      density physics) PASS. **p1 caught a live backend defect while
      closing:** Kokkos Kernels' IC(0) factorization produces a weaker
      factor on the ill-scaled Kuan matrix (dy = 0.05 / dz = 0.04 cells)
      — CG stalled at 2.1e-7 preconditioned residual, one order short of
      rtol, on both systems; under `aij` the same case converges. Fix:
      on Kokkos lanes the bjacobi-icc sub-factorization defaults to
      PETSc's host path (`sub_pc_factor_mat_solver_type petsc`), making
      the factors identical to the aij lane by construction while SpMV
      and vector work stay on Kokkos. The device-intended preconditioner
      is amg, where no sequential triangular solve exists.
- [x] **p2 — OpenMP solver thread scaling (nightly).**
      `ctest -R scaling.p2` PASS (788 s). On the 2.2M-cell case,
      aijkokkos+amg, 1 rank: gw solve 34.6 s → 21.8 s → 20.4 s at
      1/2/4 threads (**the solve threads**; hard assertion 4t < 1t),
      iterations 15/15/15 (±2 % bound; exactly 0 % drift), volumes
      identical to 1e-8. Negative battery: `scripts/test_p_gates.py`
      (ctest `unit.p_gates.negative`) — 8/8 defect injections behave.
- [x] **p3 — hybrid placements (nightly).** `ctest -R scaling.p3` PASS
      (520 s). 4 PEs on aijkokkos+amg: 4×1 = 55.5 s, 2×2 = 57.7 s,
      1×4 = 56.9 s — placement now costs ≤ 4 % (Q1's s4 measured 1×4 28 %
      slower; the Amdahl bottleneck is gone). Iterations 15 in all
      placements (≤ 10 % bound), volumes to 1e-8.
- [x] **p4 — device build integrity (CI).** `cuda-compile.yml` promoted:
      Kokkos+CUDA PETSc (`FREHG_CUDA=1` installer builds `--with-cuda`),
      `-DFREHG_BACKEND=cuda`, `FREHG_ENABLE_TESTS=ON` — compile **and
      link**, binary + test executables, zero warnings. Runs on push/PR
      (first execution needs Q0.2's real remote, tracked there).
- [x] **p5 — pointer discipline (per-PR).** (a) `check_forbidden.sh`
      rule: `Mat|Vec|KSP|PC` call receiving `.data()` outside
      `core/LinearSystem.cpp` fails the scan — clean; (b) compile-time
      invariant `PetscBackendConsistent<PETSC_HAVE_KOKKOS, MemSpace>`
      static-asserted in LinearSystem.cpp, with negative test
      `unit.p5.backend_invariant_fires` (WILL_FAIL build of the forbidden
      pair) — PASS; (c) debug-build first-solve `PetscMemType` assertion
      (device build must see device memory) + fatal on non-Kokkos vector
      types under a device build.
- [x] **p6 — GPU acceptance bundle (authored + dry-run; owner-executed on
      hardware).** `scripts/gpu_acceptance.sh` + protocol/criteria in
      `docs/developer-guide/gpu-acceptance.md`. CPU dry-run
      (`--allow-host`, both binaries = host build): backend SKIP as
      designed; invariance PASS (surface+subsurface volume 153541.295…
      identical, gw iters 21 = 21); multi-rank with `gpu_aware_mpi`
      on/off both PASS; bundle exit 0.

## Standing lanes (must not regress)

- [x] Unit tier: `ctest -L unit` 15/15 (includes the two new negative
      tests and the parameter-docs lockstep, which caught and forced the
      `mat_type` docs row).
- [x] Per-PR regression tier: `ctest -L '^regression$'` green
      post-change (b1–b4, g1 lanes, restart/rank-invariance, r1 hooks on
      every run, r2, perf-baseline).
- [x] `scripts/check_forbidden.sh` clean; `mkdocs build --strict` clean.
- [x] r1 lockstep: `mat_type` added to schema + `resolvedConfigYaml` +
      parameters.md in the same change; round-trip green on every gate
      case (r1 runs inside every regression lane).

## Documentation

- [x] `docs/user-guide/installation.md`: lane-selection table (§2B.2),
      `FREHG_BACKEND` row, GPU section rewritten to the `experimental`
      status contract.
- [x] `docs/user-guide/parameters.md`: `solver.*.mat_type` rows.
- [x] `docs/developer-guide/performance.md`: "Performance portability
      (v2 Q3)" chapter (what changed, why the solve never threaded
      before, gate map).
- [x] `docs/developer-guide/gpu-acceptance.md` (new, in mkdocs nav): the
      p6 protocol, build recipes for both binaries, criteria table,
      promotion procedure.
- [x] Device builds print the `experimental` status notice at startup
      (§2B.4).

## Known limits (stated, not finessed)

- GPU execution is **formally unverified** until the owner returns a
  passing p6 bundle from the GPU HPC; the lanes ship `experimental`
  (v2 plan §2B.4). This is a tracked open item with an owner action, not
  a silent gap.
- The HIP compile lane is not in CI (no ROCm toolchain available);
  `FREHG_BACKEND=hip` is supported by the build system and the same
  static guards apply.
- p2/p3 thresholds were measured on the fanless M3 dev machine (A23
  lottery); recalibration on real runners follows the amendment protocol.
- **The self-contained build scripts must build Kokkos shared**
  (`-DBUILD_SHARED_LIBS=ON`); a static Kokkos core is duplicated into the
  frehg executable and PETSc's `libkokkoskernels.dylib`, giving two runtime
  singletons and a first-`VecKokkos` segfault (amendment V2-A8). The
  p-gates ran against the shared `local/` prefix and never rebuilt from the
  scripts, so this packaging defect escaped; p4 is extended with a host
  self-contained-build smoke gate to close the gap.
- Amendments this phase: V2-A6 (phase insertion + s3 flat-solve
  amendment), V2-A7 (l1-scaled Jacobi on Kokkos lanes), V2-A8 (shared
  Kokkos in the build scripts + p4 build-provenance gate).
