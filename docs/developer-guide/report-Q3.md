# Q3 Report — Performance-portable parallelization

**Phase:** v2 Q3 (v2 plan §2B, inserted by amendment V2-A6)
**Completed:** 2026-09-11
**Machine:** macOS arm64 (M3, 8 cores), gcc-15, MPICH 4.3, PETSc 3.25.1
(+ hypre 3.1.0, Kokkos Kernels), Kokkos 5.1.1 (Serial + OpenMP host backends)

## What the phase delivered

One source tree now covers six deployment lanes — laptop serial/OpenMP/MPI,
HPC MPI and MPI+OpenMP hybrid, single- and multi-GPU — with the linear
solver backend-consistent in each. The headline defect fixed: **the PETSc
solve never threaded**. `MATAIJ` was hardcoded and PETSc's native matrix is
not OpenMP-parallel, so `OMP_NUM_THREADS` accelerated only the Kokkos
physics kernels while the solve stayed serial per rank (encoded in s3's
"solve asserted flat", now amended). Worse, the COO assembly handed a
`MemSpace` pointer to a host-typed matrix — correct on CPU builds where the
two spaces coincide, a latent wrong-memory defect on any device build, and
invisible to every existing lane because the CUDA lane never linked or ran.

## Measured results

**p2 — the solve now threads** (2.2M-cell synthetic, aijkokkos + amg, 1 rank):

| threads | gw solve [s] | gw iters | simulation [s] |
|---|---|---|---|
| 1 | 34.6 | 15 | 128.5 |
| 2 | 21.8 | 15 | 70.9 |
| 4 | 20.4 | 15 | 55.9 |

Iterations are *exactly* thread-invariant (the V2-A7 l1-scaled-Jacobi
smoother; hypre's hybrid host default drifts with thread count by
construction). Solve efficiency at 4 threads is 42 % — recorded, not gated:
AMG's V-cycle bandwidth-bound triangular work does not scale like the
physics kernels, and the honest hard assertions are "strictly faster than
1 thread" and "algebra unchanged".

**p3 — placement no longer matters much** (4 PEs, aijkokkos + amg):
4×1 = 55.5 s, 2×2 = 57.7 s, 1×4 = 56.9 s (≤ 4 % spread; Q1's s4 measured
1×4 28 % slower than 4×1). With the solve threaded, hybrid placement is now
governed by memory/NUMA, not by an unthreadable serial fraction — this is
the property that makes the one-rank-per-GPU model viable.

**p1 — backend choice does not change physics:** b1/b4 at 4 threads under
aij and aijkokkos pass the unchanged v1 tolerances; b6 (coupled + transport
+ density) passes nightly under aijkokkos. Direct A/B on the smoke case:
max abs field difference 1.57e-9.

**p1 earned its keep before the phase even closed.** The first b6-aijkokkos
run FAILED: CG stalled at 2.1e-7 preconditioned residual (rtol needs
1.9e-8) on the Kuan matrix — first the fs system, and with fs pinned to
host factors, then the gw system, while the identical case converges under
`aij`. Root cause: Kokkos Kernels' IC(0) factorization is measurably weaker
than PETSc's on this ill-scaled matrix (0.05 m × 0.04 m cells, density
coupling). Fix: on Kokkos lanes, bjacobi-icc keeps its sub-factorization on
PETSc's host path — factors identical to the aij lane by construction,
SpMV/vector work still on Kokkos. b1/b4 had passed under the Kokkos factors
(better-conditioned systems); only the full-physics nightly case exposed
it, which is exactly why the b6 lane is in the p1 set.

**p6 dry-run (CPU rehearsal):** all items PASS with volumes byte-identical
(153541.295… m³) and iterations equal; item 1 (device backend) SKIP as
designed under `--allow-host`.

## Decisions and their reasons

- **`aijkokkos` is opt-in on host, forced on device.** The host `aij` lane
  is the golden-pinned v1 behavior; nothing moves under existing configs.
  Device builds cannot correctly use host types, so they are not allowed to.
- **PMIS on device, HMIS on host** (BoomerAMG coarsening): PMIS is the only
  coarsening hypre implements on GPUs; `ext+i` interpolation (already the
  Q1 choice) is device-supported and stays everywhere.
- **l1-scaled Jacobi smoother on the Kokkos lanes (V2-A7):** costs ~1.7×
  iterations (8-9 → 15 on the gate case) but is thread-invariant and is the
  smoother family the device path uses — so the OpenMP lane rehearses the
  exact GPU algebra, which is the point of p1. Host-aij amg keeps hypre's
  default; the Q1 g2 record is untouched.
- **GPU ships `experimental` (§2B.4).** §8.3 forbids authored-unexercised
  features; CI has no GPU. The honest resolution is a named status —
  compile+link-verified (p4), statically invariant-checked (p5), CPU-physics
  -verified (p1–p3), printed at startup on device builds — cleared only by
  the owner-run p6 bundle. The promotion is one table edit plus committing
  the returned report.

## What the owner runs on the GPU HPC

1. Build the CPU reference and the device binary from the same commit
   (recipes in `docs/developer-guide/gpu-acceptance.md`; `build_frehg2_hpc.sh`
   carries the annotated configure line).
2. `scripts/gpu_acceptance.sh --frehg-gpu … --frehg-cpu … --mpiexec … --gpus N`
3. Send back `gpu-acceptance-report.txt` + the collected `record-*.yaml`.
   The script prints PASS/FAIL per item; a passing bundle closes the loop.

## Infrastructure notes for later phases

- Every new physics kernel (Q4–Q6) inherits p1/p5 automatically: the
  forbidden-pattern rule rejects raw `.data()` into PETSc outside
  `LinearSystem`, the compile-time invariant guards the backend pairing,
  and the p1 lanes run per-PR. This was the reason Q3 preceded the physics.
- `run_regression.py --mat-type` and `run_scaling.py --mat-type` compose
  with `--solver`; the run record now carries `mat_type`,
  `amg_coarsen_type`, `amg_relax_type` per system (resolved, not requested).
- `petscvec_kokkos.hpp` must be the first PETSc include in any TU using the
  Kokkos Vec API — PETSc enforces it with an #error.
- CI cache keys moved to `…-kokkoskernels-v3`; the first real-runner
  execution (Q0.2, still open) will rebuild all dependency caches.
- **Build Kokkos shared** (`-DBUILD_SHARED_LIBS=ON`, amendment V2-A8). A
  static Kokkos core is absorbed independently into the frehg executable and
  PETSc's `libkokkoskernels.dylib` → two runtime singletons → doubled
  `Kokkos::OpenMP::initialize` → segfault on the first `VecKokkos` access.
  This bit a fresh `build_frehg2_local.sh` build (frehg2-dev) even though
  every p-gate was green — the gates ran against the pre-existing *shared*
  `local/` prefix and never rebuilt the stack from the shipped script. p4 is
  extended with a host self-contained-build smoke gate so the build script is
  itself a gated artifact. Recovery from a static build: delete the Kokkos
  and PETSc installs from the deps prefix (PETSc embedded the static core and
  must rebuild too) and rerun; MPICH/HDF5/yaml-cpp are reused.

## Gate ledger

| Gate | Status | Evidence |
|---|---|---|
| p1 per-PR (b1/b4 × aij/aijkokkos @ 4 threads) | PASS 4/4 | ctest regression.p1.* |
| p1 nightly (b6 ss aijkokkos @ 4 threads) | PASS | ctest regression.p1.b6.ss.omp4.aijkokkos |
| p2 solver thread scaling | PASS | ctest scaling.p2 (788 s); build/scaling/p2.json |
| p3 hybrid placements | PASS | ctest scaling.p3 (520 s); build/scaling/p3.json |
| p4 device compile+link | armed in CI | cuda-compile.yml (first run = Q0.2) |
| p5 invariant + negative + scan | PASS | unit.p5.backend_invariant_fires; check_forbidden.sh |
| p2/p3 gate-logic negatives | PASS 8/8 | unit.p_gates.negative |
| p6 bundle | authored + CPU dry-run PASS | /tmp dry-run; owner action pending |
| Standing per-PR tier | PASS | ctest -L '^regression$' post-change |
| Unit tier | PASS 15/15 | ctest -L unit |
