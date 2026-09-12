# Performance report (P5)

All measurements: Apple M3 MacBook Air (4 performance + 4 efficiency
cores, fanless, 24 GiB unified memory), macOS, gcc-15 `RelWithDebInfo`
build, MPICH 4.3 (no rank binding exists on Darwin), PETSc 3.25.1, Kokkos
5.1.1 OpenMP host backend at `OMP_NUM_THREADS=1`, `FI_PROVIDER=tcp`.
Reproduce with:

```bash
python3 scripts/run_scaling.py --case synthetic --ranks 1 4 --repeats 3 \
    --frehg build/src/frehg --mpiexec mpiexec --repo . --work /tmp/scaling
```

## Methodology

Strong scaling is measured on **fixed-step runs**: the common adaptive
step is clamped (`dt_init = dt_min = dt_max`), so every rank count
marches the identical step sequence and does identical arithmetic work
(the adaptive controller forks across decompositions — amendment A14 —
which would otherwise contaminate wall-clock ratios). The gate value is
the `simulation` timer (the whole time loop, max across ranks).

On this machine, run-to-run variance is dominated by two platform
effects, both strictly additive: thermal DVFS (a fanless chip throttles
multi-core bursts within a minute) and rank placement (macOS may schedule
ranks onto E-cores; the slowest rank gates every collective; MPI offers
no binding on Darwin). `--repeats` therefore takes the **minimum** over
repeats — the standard estimator when noise can only slow you down. Every
attempt is published below (amendment A23).

## The scaling gate (amendment A23): 16× cells, 30-step burst

The b5 benchmark grid (139k cells) is too small to occupy four ranks on
this machine (see below), so the gate case is the same coupled physics —
rain infiltrating toward a water table; identical predictor / corrector /
reallocation / exchange kernels — on a 404×220×25 grid (2.2M cells),
30 fixed 2 s steps (long enough to amortize startup, short enough that a
4-core burst stays inside the thermal envelope).

Recorded gate artifact (`--ranks 1 4 --repeats 3`, min-selected):

| ranks | attempts [s] | min [s] | speedup | efficiency |
|---|---|---|---|---|
| 1 | 126.07, 127.05, 122.90 | 122.90 | 1.00 | 100 % |
| 4 | **38.89**, 48.69, 54.76 | 38.89 | 3.16 | **79.0 % — PASS (≥ 70 %)** |

Capability ensemble from the same case (back-to-back single runs on a
rested machine, order n1 → n4 → n4 → n2 → n8):

| ranks | simulation [s] | speedup | efficiency |
|---|---|---|---|
| 1 | 129.24 | 1.00 | 100 % |
| 4 | 38.10 (repeat 41.15) | 3.39 | **84.8 %** (repeat 78.5 %) |
| 2 | 72.22 | 1.79 | 89.5 % |
| 8 | 44.30 | 2.92 | 36.5 % |

The variance that motivates the min-over-repeats protocol, all on
identical work: 4-rank bursts spanned 38.1–62.0 s across the measurement
day, one sequence had 8 ranks outrun 4 (43.5 vs 62.0 s — placement, not
code), and *sustained* multi-minute 4-core loads (60-step horizon,
min-of-3 on a heat-soaked chip) capped at 51 % while 2-rank efficiency
held at 87–100 % throughout. Eight ranks always land on E-cores here and
are reported, not gated. On a bindable, uniform-core, actively cooled
node the script's exit-code gate applies unchanged and none of these
caveats exist.

## The b5 grid itself (informational)

`--case b5`: rain/sync, fixed dt = 2 s, 3600 steps (t = 7200 s), single
runs:

| ranks | simulation [s] | speedup | efficiency |
|---|---|---|---|
| 1 | 1070.60 | 1.00 | 100 % |
| 2 | 611.00 | 1.75 | 87.6 % |
| 4 | 553.04 | 1.94 | 48.4 % |
| 8 | 483.63 | 2.21 | 27.7 % |

Per-rank section spreads are tight (no load imbalance) and the plateau
hits the pure-Kokkos kernels and the PETSc solve alike — at 139k cells
there is simply too little work per rank per step (the launch-latency
regime report-P3/A15 recorded; the subsurface solve *rises* from 158 s at
2 ranks to 192 s at 8 on reduction latency).

## Per-module cost (b5, fixed 2 s step, serial)

| section | time [s] | share |
|---|---|---|
| simulation (3600 steps) | 1070.6 | 100 % |
| └ coupled_step | 1069.8 | 99.9 % |
| &nbsp;&nbsp;└ groundwater | 1066.4 | 99.6 % |
| &nbsp;&nbsp;&nbsp;&nbsp;└ groundwater/solve (PETSc CG) | 258.0 | 24.1 % |
| &nbsp;&nbsp;└ swe/free_surface | 2.4 | 0.2 % |
| &nbsp;&nbsp;&nbsp;&nbsp;└ free_surface/solve | 1.2 | 0.1 % |
| &nbsp;&nbsp;└ swe/velocity | 0.7 | 0.1 % |

At the common adaptive step the coupled b5 regime is entirely
subsurface-bound: ~76 % of runtime is the PCA machinery outside the
linear solve (predictor assembly, van Genuchten evaluations, corrector
fluxes, moisture reallocation, audits) and ~24 % the 7-point CG solve.
The surface module is negligible at this step size. Optimization effort,
if ever needed, should start in the subsurface kernels, and any change
there must preserve the term-by-term audits or re-gate (report-P4 §2).

## Practical guidance

- **Small (benchmark-scale) grids**: run `OMP_NUM_THREADS=1`, serial or
  few ranks — kernel-launch latency dominates (b1 runs in ~0.3 s at one
  thread vs tens of seconds at the default thread count; report-P1).
- **Large grids**: MPI ranks scale well while per-rank work stays large
  (≥ ~0.5M cells/rank here reached 79–85 % at 4 ranks); add OpenMP
  threads only once per-rank subdomains are large.
- Benchmark gate runtimes for planning: b1/b2 seconds, b3/b4 ~0.5 min,
  b6 ss ~6 min / td ~13 min serial, b5 shortened rain/sync gate ~1 h at
  4 ranks, b5 full horizon ~3.6 h per lane at 4 ranks (A15).

## Solver scalability (v2 Q1)

The v1 default preconditioner — CG with block-Jacobi/ICC(0) — is not
scalable: block count equals rank count, cross-block coupling is dropped,
and CG iterations grow with ranks. v2 adds YAML-selectable AMG paths
(`solver.<system>.preconditioner: amg` for hypre BoomerAMG, `gamg` for
PETSc smoothed aggregation) with an iteration-and-cadence-triggered
hierarchy-reuse policy (v2 plan §2.2.4) and per-system telemetry (the
end-of-run `solver summary` lines; parsed by `scripts/run_scaling.py`).

Standing gates (v2 plan §2.3, §7.2): g1 re-runs the b-gates under each
preconditioner (solver choice must not change physics); g2 asserts mean CG
iterations at 8 ranks stay within 1.10× the 1-rank count under AMG on the
A23 synthetic case; s1 keeps the 70 %-at-4-ranks strong-scaling criterion
permanent; s2 (weak), s3 (OpenMP kernel-time), and s4 (hybrid placement)
run nightly via `scripts/run_scaling.py --gate ...`. The per-module
perf-regression gate compares timers against
`docs/developer-guide/perf-baseline.json` (fail > +25 %).

The default preconditioner remains `bjacobi-icc` in v2.0: the b1–b6
goldens were gated under it, and flipping the default is a v2.1 decision
once g3 history exists on production hardware (v2 plan §2.2.1).

## Run provenance and timing record (v2 Q2)

The timer tree, solver telemetry, and mass-audit closure persist per run as
`run-record.yaml` beside the HDF5 output (see the user guide's output
reference for the format). New timer sections complete the breakdown:
`init` (pre-loop restore/first-output), `io/output`, `io/checkpoint`,
`io/run_record`, `halo` (all message passing — every exchange funnels
through one timed path), and `monitors`. Gates: r1 validates every
regression gate run's record including a byte-exact configuration
round-trip against `frehg --resolve`; r2 asserts the top-level sections
tile >= 90 % of the time loop (measured: 100 %), halo/io alive at 4 ranks,
and record writing <= 1 % of the loop (measured: 0.03 %).

## Performance portability (v2 Q3)

The execution backend is a build choice and the linear-algebra backend
follows from it (v2 plan §2B.2; lane table in the user guide's
installation chapter). What Q3 changed:

- `solver.<system>.mat_type: aijkokkos` runs the PETSc solve through Kokkos
  Kernels on the build's execution space — the path that threads the solve
  on CPU and is the *only* path on device builds (where it is forced). The
  RHS/solution handoff uses PETSc's Kokkos view API: zero copies on host,
  one device-side copy on GPU, never a device↔host round trip per solve.
- BoomerAMG defaults are backend-aware: HMIS coarsening on host, PMIS on
  device (the only GPU-capable coarsening); on the Kokkos lanes the smoother
  defaults to l1-scaled Jacobi because hypre's hybrid SOR/Jacobi host default
  is thread-count-dependent (measured: gw 8→15 iterations moving to l1, but
  identical at 1/2/4/8 threads, vs iteration drift under the default). All
  remain `setOptionDefault` — an options file or command line wins.
- Why the solve did not thread before Q3: PETSc's native `MATAIJ` is not
  OpenMP-threaded, so `OMP_NUM_THREADS` only ever accelerated the Kokkos
  physics kernels while the solve stayed serial per rank (the s3 gate's
  "solve asserted flat" encoded that design; amended by V2-A6).

Standing gates (v2 plan §2B.3): p1 re-runs b1/b4 (nightly: b6) at 4 threads
under both `aij` and `aijkokkos` — backend choice must not change physics;
p2 asserts the aijkokkos+amg gw solve actually threads (4-thread solve
< 1-thread) with thread-invariant algebra (iteration drift ≤ 2 %); p3 runs
4 PEs as 4×1/2×2/1×4 with volumes to 1e-8 and iterations within 10 % of
all-MPI; p4 compiles **and links** the full tree, tests included, under
Kokkos+CUDA in CI; p5 enforces pointer discipline (forbidden-pattern rule,
compile-time backend invariant with a WILL_FAIL negative test, debug
memtype assertions). GPU execution itself is gated by the owner-run p6
bundle (`scripts/gpu_acceptance.sh`, protocol in `gpu-acceptance.md`);
until it returns, GPU lanes are `experimental`.
