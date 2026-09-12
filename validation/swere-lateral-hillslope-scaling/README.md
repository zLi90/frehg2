# swere-lateral-hillslope-scaling — MPI strong-scaling case

A large, homogeneous 3D extension of
[`swere-lateral-hillslope`](../swere-lateral-hillslope/) (SERGHEI case4), built
to measure **MPI strong scaling**: a fixed-size problem run on an increasing
number of ranks, reporting speedup and parallel efficiency.

## What changed from the parent case

| | parent `swere-lateral-hillslope` | this case |
|---|---|---|
| x | nx 40, dx 100 m | **nx 400, dx 10 m** (same 4000 m width) |
| y | ny 1, dy 100 m | **ny 800, dy 10 m** (x-z plane duplicated in y) |
| z | nz 60, dz 0.25 m | nz 60, dz 0.25 m (unchanged) |
| cells | 2 400 | **19 200 000** |
| run length | 5 yr (157 680 000 s) | **1 hour (3 600 s)** |
| initial subsurface | 3D head file | **water table** (compact, see below) |

The **physics is unchanged**: the same 4000 m bed line (`elev = 3.95 − 0.001·x`),
the same van Genuchten soil, the same lateral head boundary conditions
(−4.1 m left, −14.1 m right), and the same rainfall + coupled-seepage forcing.
Because every y-slice is identical and the y-walls are no-flux, there is **no
y-gradient** — the widened domain reproduces the parent x-z hillslope on each
slice. That is exactly what you want for a scaling study: a big, uniform,
perfectly load-balanced workload whose per-rank cost is predictable.

## Files

```
swere-lateral-hillslope-scaling.yaml   the BASELINE run (decomposition 1 x 1)
gen_inputs.py                          regenerates the input rasters
input/dem.dat                          bed elevation   (ESRI raster, 400 x 800)
input/water_table.dat                  initial water table (ESRI raster, 400 x 800)
input/rain.dat                         constant rainfall time series
out/                                   run outputs (output.h5, monitors)
```

Regenerate the inputs (e.g. after editing the grid in `gen_inputs.py`):

```bash
python3 gen_inputs.py
```

## Design choices

- **`dy = 10 m`** (isotropic horizontal cells). Since the domain is homogeneous
  in y and the y-walls are no-flux, the physics is independent of `dy`; 10 m is
  chosen only so cells are square in plan. `ny = 800` factorises cleanly
  (2·4·5·5·... ) for a wide range of y-decompositions.

- **Water-table initial condition, not a 3D head file.** The parent case ships a
  full 3D pressure-head field; at this resolution that file would be ~150 MB of
  ASCII. Instead we hand frehg2 the water-table elevation
  (`initial_conditions.groundwater.water_table`, a ~3 MB 2D raster) and let
  `RichardsSolver` build the hydrostatic head `h = (WT + datum_offset) − z_cell`.
  The datum offset cancels against the mesh datum, so the **saturated-zone head
  matches the parent case exactly**; only the parent's −1.25 m unsaturated clamp
  is replaced by a clean hydrostatic profile above the table (a gentler start,
  which the solver relaxes anyway). The water table runs linearly from the two
  lateral head BCs: −4.1 m at the left edge to −14.1 m at the right.

- **Analytic inputs.** `gen_inputs.py` evaluates the parent bed line and water
  table at the refined cell centres, so the refined + duplicated domain is the
  *same* physical hillslope rather than an interpolation artefact. Both rasters
  are y-invariant (one nx-long row repeated ny times, in the reader's `(j·nx+i)`
  order — no ESRI north-up flip, matching the parent).

- **Minimal I/O.** `output_interval = t_end = 3600 s` writes fields once (at the
  end) so disk I/O stays out of the timing. Even so, one snapshot of the two 3D
  groundwater fields is ~0.47 GB — drop `hydraulic_head`/`water_content` from
  `output.variables.groundwater` if you want the run purely compute-bound.

## Running the strong-scaling study

The grid is **fixed** at 400 × 800 × 60 for every rank count. Strong scaling
changes **only** `domain.decomposition` (and the launch rank count); the input
rasters are shared across all runs. `mpi_nx · mpi_ny` must equal the number of
MPI ranks.

This file is the **baseline: `decomposition {mpi_nx: 1, mpi_ny: 1}`** — the
single-rank reference. Suggested ladder (y is the long axis, so decompose y
first; mix in x for large rank counts):

| ranks | decomposition | ranks | decomposition |
|---|---|---|---|
| 1  | `mpi_nx: 1, mpi_ny: 1`  | 16 | `mpi_nx: 1, mpi_ny: 16` |
| 2  | `mpi_nx: 1, mpi_ny: 2`  | 20 | `mpi_nx: 2, mpi_ny: 10` |
| 4  | `mpi_nx: 1, mpi_ny: 4`  | 32 | `mpi_nx: 2, mpi_ny: 16` |
| 8  | `mpi_nx: 1, mpi_ny: 8`  | 40 | `mpi_nx: 4, mpi_ny: 10` |

Copy this YAML per rank count (or edit `decomposition` in place) and launch with
the matching `-np`. On a workstation:

```bash
OMP_NUM_THREADS=1 mpirun -np 4 ../../build/src/frehg swere-lateral-hillslope-scaling.yaml
```

### Reading the result

At the end of every run frehg2 prints a **timer report**. The `simulation`
section is the whole time loop; its **`max[s]`** column (slowest rank) is the
strong-scaling metric:

```
timer report (4 ranks)
  section                        count      min[s]     mean[s]      max[s]
  simulation                         1      ...          ...       T_N   <- use this
  ...
```

With `T_1` from the baseline and `T_N` from the N-rank run:

- speedup    `S_N = T_1 / T_N`
- efficiency `E_N = T_1 / (N · T_N)`

Keep `OMP_NUM_THREADS=1` (one thread per rank) for every run so the comparison
is pure MPI, and run the rank counts back-to-back on an otherwise idle node.
The repo's `scripts/run_scaling.py` is hardwired to its own two cases and will
**not** drive this one — read the `simulation` timer directly.

## Cost and memory

- **Baseline (1 rank):** ~19.2 M cells in one PETSc/Kokkos address space. Expect
  several GB of RAM (PCA matrix + the Field3 working set) and a long wall time —
  it is meant to be the slow reference. If a single rank is impractical on your
  node, use the smallest rank count that fits as the reference and report
  efficiency relative to it (state which reference you used).
- **Output:** ~0.47 GB per field snapshot (see *Minimal I/O* above).

## HPC launch

Load the exact stack frehg2 was built against and launch with `mpirun`. On the
Intel cluster documented in `docs/agents/build.md` §11, that site's OpenMPI has
no Slurm PMI and a broken hwloc, so the working job is:

```bash
module purge
module load gcc cmake openmpi/4.1.6 hdf5 OpenBLAS/0.3.26
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export OMPI_MCA_rtc=^hwloc OMPI_MCA_hwloc_base_binding_policy=none
mpirun -np "$SLURM_NTASKS" /path/to/build/src/frehg \
    swere-lateral-hillslope-scaling.yaml
```

(On a cluster whose MPI has PMI/PMIx and a working hwloc, `srun` works and the
two `OMPI_MCA_*` lines are unnecessary.) Match `-n`/`SLURM_NTASKS` to
`mpi_nx · mpi_ny` in the YAML.
