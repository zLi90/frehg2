# Running a simulation

Run by passing the configuration file to the executable:

```bash
build/src/frehg my-case.yaml
```

Output paths in the configuration are resolved against the **current working
directory**, so it is usually cleanest to run from the case directory:

```bash
cd benchmarks/b1-sw
OMP_NUM_THREADS=1 /path/to/build/src/frehg b1-sw.yaml
```

Frehg2 logs the effective configuration, then one line per output time
with the solver iteration counts (free-surface and/or subsurface,
depending on the enabled modules — sync-coupled runs also show the current
adaptive step), and finally a timer report.

## Parallel (MPI) runs

```bash
mpirun -np 4 build/src/frehg my-case.yaml
```

The domain is decomposed automatically over a 2D rank grid (respecting the
domain extents), or you can pin it with `domain.decomposition.mpi_nx/mpi_ny`.
Results are **rank-invariant**: the same case on 1, 2, or 4 ranks produces
the same fields (to round-off). The single HDF5 output file is written
collectively by all ranks.

## Extra PETSc arguments

Anything after the configuration path is handed to PETSc. For example, the
strict solver mode used for the bitwise rank-invariance proof:

```bash
build/src/frehg my-case.yaml -fs_pc_type jacobi -fs_ksp_rtol 1e-13 -fs_ksp_atol 1e-16
```

The free-surface system uses the `fs_` option prefix; the subsurface
system uses `gw_`.

## Threads and performance

Frehg2's host backend is Kokkos OpenMP. For **small (benchmark-scale)
grids**, the work per step is dominated by kernel-launch latency, so
`OMP_NUM_THREADS=1` is dramatically faster (b1 runs in ~0.3 s at one thread
versus tens of seconds at the default thread count). Set threads greater than
one only for large grids. For production:

```bash
OMP_NUM_THREADS=1 build/src/frehg my-case.yaml        # small grids
OMP_NUM_THREADS=8 mpirun -np 4 build/src/frehg big.yaml  # large grids
```

The measured strong-scaling behavior of the coupled solver (b5, 1→8 ranks)
and the per-module cost breakdown are in the developer guide's
[performance report](../developer-guide/performance.md).
