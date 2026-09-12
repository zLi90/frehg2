# GPU acceptance protocol (gate p6, v2 plan §2B.3–§2B.4)

The CI machines have no GPU, so the GPU lanes ship as **experimental**:
compile+link-verified (p4), statically invariant-checked (p5), and
physics-verified on the equivalent CPU backends (p1–p3) — but never executed
on device hardware by CI. This document is the protocol that clears that
status. It is executed **by the owner on a GPU machine**; everything else in
Q3 was verified on the CPU lanes.

## Why CPU verification does not substitute

A Kokkos kernel that is correct and thread-invariant under the OpenMP backend
is very likely correct under CUDA — that is what p1–p3 establish. What OpenMP
*cannot* see is the host/device memory boundary: under OpenMP,
`MemSpace == HostSpace`, so every raw pointer and MPI buffer aliases host
memory and a memory-space mistake is invisible (the worked example: v1 handed
`MatSetValuesCOO` a `MemSpace` pointer while hardcoding host `MATAIJ` —
correct under OpenMP, wrong under CUDA, caught by nothing until Q3). p4/p5
cover that boundary statically; p6 is the one dynamic confirmation.

## Building the two binaries (same commit)

On the GPU machine, build the dependencies twice into separate prefixes:

1. **CPU reference:** `build_frehg2_hpc.sh` as shipped (Serial+OpenMP Kokkos,
   Kokkos-enabled PETSc).
2. **Device binary:** the same script with the step-3 Kokkos built with
   `-DKokkos_ENABLE_CUDA=ON -DKokkos_ARCH_<arch>=ON` and
   `CMAKE_CXX_COMPILER=<kokkos>/bin/nvcc_wrapper`, and the step-4 PETSc
   configure given `--with-cuda=1` (HIP: `-DKokkos_ENABLE_HIP=ON` /
   `--with-hip=1`). Frehg2 configures with `-DFREHG_BACKEND=cuda` — the
   CMake consistency check rejects a mismatched PETSc at configure time.

## Running the bundle

```bash
scripts/gpu_acceptance.sh \
  --frehg-gpu  <device-build>/src/frehg \
  --frehg-cpu  <cpu-build>/src/frehg \
  --mpiexec    $(command -v mpiexec) \
  --gpus       <number of GPUs to test> \
  --work       gpu-acceptance
```

One rank per GPU; if the MPI needs explicit device binding, export the
site's usual variable (e.g. `MPICH_GPU_SUPPORT_ENABLED=1`,
`CUDA_VISIBLE_DEVICES` per rank via the site wrapper) before launching.

## Items and criteria

The case is the s3/p2 synthetic grid (404×220×25 ≈ 2.2 M cells, coupled
SWE+Richards, amg on both systems, 30 fixed steps). Both binaries pin the
device-capable AMG configuration (PMIS + l1-scaled Jacobi) on the command
line so host-vs-device default divergence cannot masquerade as a defect.

| Item | Assertion | Bound |
|---|---|---|
| p6.backend | the device binary's run record reports a device `kokkos_backend` | must contain Cuda/HIP/SYCL |
| p6.invariance | final (surface+subsurface) volume, single GPU vs 1-rank CPU | ≤ 1e-8 relative |
| p6.invariance | gw mean iterations, single GPU vs 1-rank CPU | ≤ 10 % drift |
| p6.timing | single-GPU vs 1-rank-CPU simulation time | recorded, not gated |
| p6.multi-gpu | volumes at 2 (and `--gpus`) ranks, `gpu_aware_mpi` **on and off** | each ≤ 1e-8 vs single GPU |

The volume bound is the s3/s4 cross-configuration bound, reused unchanged.
The iteration bound is loose by design: PMIS coarsening is decomposition- and
platform-sensitive; a drift beyond 10 % indicates a real setup difference,
not rounding.

## Returning the result

Send back `gpu-acceptance-report.txt` plus the `record-*.yaml` files the
script collects (they carry the full provenance: backend, versions,
decomposition, timer tree, solver telemetry). On a passing bundle, the GPU
rows in the §8.3 feature-interaction table flip from `experimental` to
`supported` and the report is committed beside this file — that is the whole
promotion procedure. On a failing bundle, include the `*.log` files; the
first suspects are the site MPI's GPU-awareness (item 4 with
`gpu_aware_mpi: on`) and the rank→device binding.

## Dry-run (CI rehearsal)

The bundle itself is exercised on the CPU lanes before release (so the only
untested variable on the GPU machine is the hardware):

```bash
scripts/gpu_acceptance.sh --frehg-gpu build/src/frehg --frehg-cpu build/src/frehg \
  --mpiexec $(command -v mpiexec) --work /tmp/p6-dryrun --allow-host
```

`--allow-host` turns item 1 into SKIP; every other item executes the
identical code path.
