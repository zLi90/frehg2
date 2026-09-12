# Q1 phase report — Scalable preconditioner

**Scope:** v2 development plan §2 (Q1). **Closed:** 2026-09-10.
**Companion:** [dod-Q1.md](dod-Q1.md) (verification commands),
`FREHG2_V2_DEVELOPMENT_PLAN.md` amendment log (V2-A1…V2-A3).

## 1. What was delivered

The v1 solves — CG + block-Jacobi/ICC(0) hardcoded for both systems — are
now YAML-selectable per system (`solver.surface`, `solver.groundwater`):

- `bjacobi-icc` (default, golden-pinned v1 behavior, bit-for-bit
  constructor semantics preserved);
- `amg` — hypre BoomerAMG, the primary scalable path: HMIS coarsening,
  ext+i interpolation (P_max 4), strength threshold 0.5 + one
  aggressive-coarsening level on the anisotropic 3D `gw_` system
  (hypre's 3D guidance; thin layers dz ≪ dx make vertical coupling
  dominate the strength graph), 0.25 / none on the 2D `fs_` system;
- `gamg` — PETSc smoothed aggregation (agg, nsmooths 1, threshold 0.02):
  no external dependency, and the native device path if a Kokkos-enabled
  PETSc ever lands.

Every injected option uses the set-if-absent pattern, so
`solver.petsc_options_file` and command-line PETSc options override
everything, exactly as in v1.

**Hierarchy reuse (v2 plan §2.2.4).** AMG setup costs several solves and
the matrix is re-assembled every step, so amg/gamg freeze the hierarchy
(`KSPSetReusePreconditioner`) and rebuild on: the `reuse_max_solves`
cadence (default 50), iteration growth past `reuse_iteration_factor` ×
the post-rebuild count + 2, an external `forceRebuild()` — wired to an
Allreduced wet/dry-mask checksum in the surface solver, since an A13
dry-closure flip is a structural operator change — and a
rebuild-and-retry-once path when a *reused* hierarchy diverges (V2-A2;
fresh-hierarchy divergence stays fatal). All triggers are KSP-collective,
so every rank takes the same branch. `bjacobi-icc` keeps the v1
rebuild-every-solve semantics and ignores the reuse knobs.

**Telemetry (v2 plan §2.2.6).** `LinearSystem` accumulates per-system
solves, iteration mean/max, rebuilds, retries, and a timed
KSPSetUp/KSPSolve split; the driver prints rank-0
`solver summary fs|gw:` lines that `scripts/run_scaling.py` parses — the
g2/g3 measurement contract.

## 2. The dry-row audit (v2 plan §2.2.5)

Findings on the v1 assembly, recorded as required:

- **Symmetry: already satisfied.** Prescribed-stage (eta-BC) columns were
  already eliminated symmetrically into the right-hand side
  (`FreeSurface.cpp`, the P1 SPD design), and the A13 dry-cell closure
  keeps dry rows live and symmetric. The debug-build first-solve symmetry
  check stands guard.
- **No identity rows in the gw system at all:** masked cells have no rows
  (compressed gids), so the hazard class is absent there by construction.
- **One finding, fixed:** eta-BC rows carried a literal 1.0 diagonal among
  wet diagonals of cell-area scale — the poorly scaled decoupled row that
  perturbs AMG strength-of-connection coarsening. The prescribed-stage
  rows now scale the diagonal and right-hand side by the cell area; the
  row solution is unchanged, and the b1 golden gate (the only per-PR gate
  with eta-BC rows) re-verified green.

## 3. Gate record

| Gate | Result |
|---|---|
| g1 per-PR (b1–b4 × {amg, gamg}; rank-invariance default-mode b1/b2 × amg at 1/2/4 ranks) | **10/10 PASS**, criteria unchanged from plan §9 |
| g1 nightly (b5 rain/sync × amg at 4 ranks on the A17 24 h record horizon; b6 ss/td × amg; b5 rank-invariance default × amg) | **4/4 PASS** — b5: discharge 0/48 and ponding 0/54 points outside, peaks +2.7 %/−1.5 %, integrals +2.5 %/−2.5 % (the A17-record numbers reproduced); b6 ss interface MAE 0.0334 m (v1 record 0.0334), td 0.0390 m (v1 0.0389); salt mass +5.5 %/+8.7 % (allowed 10 %). Two adjudications were required and are logged as **V2-A4**: the b5 amg lane initially ran the never-adjudicated full horizon (only the ponding integral failed, −11.7 % vs 10 %, with all 165 pointwise checks green) and now matches the parent's approved `--t-end 86400` record; and the default-mode rank-invariance bounds are per-solver (bjacobi's data-derived A14 bounds don't transfer — solver choice alone at n=1 moves the observables past them; amg bounds set at ~2× the Q1 measurements: 4e-4 / 1.6e-2 / 3e-3). |
| g2 (iteration flatness, A23 synthetic, {1,2,4,8} ranks, amg) | **PASS**: gw mean 8.47 → 9.00, ratio 1.06 ≤ 1.10; recorded count seeded (`tolerances/g2.yaml`) |
| g2 record lane (bjacobi-icc, same case) | measured for the motivation record — §4 table |
| s1 strong (nightly) | PASS — §4 |
| s2 weak, amg (nightly) | PASS — §4 |
| s3 threads (nightly) | PASS — §4 |
| s4 hybrid (nightly) | PASS — §4 |
| perf-baseline | seeded (min over 3 repeats) and check green |

**Closing sweep (2026-09-10):** `ctest -L '^unit$|^regression$' -j 1` —
**38/38 passed** in one pass (unit 11; per-PR regression 27: b1–b4,
restarts, b6 smoke, strict+default rank invariance, the 10 g1 per-PR
lanes, g2, perf-baseline; 1891 s wall).

**Iteration counts per solve, A23 synthetic case (the headline numbers,
measured by the two g2 lanes):**

| ranks | gw mean, bjacobi-icc | gw mean, amg |
|---|---|---|
| 1 | 23.0 | 8.5 |
| 2 | 24.0 | 9.0 |
| 4 | 24.0 | 9.0 |
| 8 | 25.0 | 9.0 |

Block-Jacobi/ICC needs 2.6–2.8× the iterations of BoomerAMG on this case
and grows ~9 % over 1→8 ranks even while its blocks are still ~280k rows;
the growth compounds as blocks shrink at production rank counts (the
structural argument the literature quantifies — v2 plan §2.2 precedent),
while BoomerAMG is flat at 9 after the first-rank transient. This is the
documented motivation for Q1, now measured in CI rather than asserted.
(The case's fs system converges in 0 iterations — V2-A1 — so fs flatness
is exercised through g1.)

## 4. s-gate measurements (reference M3, min over 3 repeats)

All five `scaling_nightly` entries green locally (2026-09-09,
`ctest -L scaling_nightly -j 1`, 2319 s wall); JSON artifacts in
`build/scaling/` (uploaded by the nightly workflow on real runners).

- **s1 strong (bjacobi default, A23 synthetic):** 112.0 → 59.3 → 35.5 s
  over 1→2→4 ranks: efficiency 94.5 % at 2, **78.9 % at 4** (criterion
  ≥ 70 %; the P5-class measurement reproduced post-Q1).
- **s2 weak (amg, 202×110×25 per rank, 1→2→4 ranks):** efficiency 92.9 %
  at 2, **74.9 % at 4** — above the 0.5 hard floor, below the 0.8 target
  (WARN, the expected fanless-M3 regime); gw iterations 8.47 → 9.00
  (**+6.3 %**, gated ≤ 15 %); setup fraction ~1 % of solve time (the
  reuse policy amortizes as designed).
- **s3 threads (1 rank × {1,2,4} OpenMP, 2.2M cells):** simulation
  112.9 → 67.3 → 46.0 s; kernel time (simulation − solve) 92.9 → 25.7 s
  = **90 % kernel efficiency at 4 threads** (criterion ≥ 60 %); solve
  time flat at 20.0–20.3 s (the MPI-only solve, as asserted); mass-audit
  volume identical across thread counts within 1e-8 relative.
- **s4 hybrid (4 PEs as 4×1 / 2×2 / 1×4):** 41.7 / 44.5 / 53.4 s, all
  placements volume-identical within 1e-8 relative. 4×1 is fastest —
  consistent with the solve being MPI-only (more ranks = more parallel
  solve), which is the guidance `performance.md` gives users.

Thresholds are provisional until Q0.2 puts these lanes on real runners
(recalibration by amendment, the A23 discipline).

## 5. Decisions and deviations

- **Default stays `bjacobi-icc` in v2.0** (plan §2.2.1): the b-gate
  goldens were adjudicated under it; flipping the default is a v2.1
  decision once g3 history exists on production hardware.
- **V2-A1:** g2's recorded-count gate binds the gw system only (the
  synthetic case's fs solves are trivial); fs AMG behavior is covered by
  g1.
- **V2-A2:** stale-hierarchy divergence rebuilds and retries once before
  the fatal rule; fresh-hierarchy divergence unchanged from v1.
- **V2-A3:** 8-rank g2 measurements are valid on the 4P+4E machine —
  iteration counts are topology-independent; timing gates keep the A23
  4-rank discipline.
- **V2-A4:** the g1 b5 amg lane is held to the parent's A17-approved 24 h
  record horizon, and the default-mode b5 rank-invariance bounds are
  per-solver (data-derived for amg/gamg from the Q1 measurements, the
  A5/A8/A14 discipline). Full details and the archived full-horizon
  metrics in the plan's amendment log.
- The mkdocs strict build had rotted with mkdocs 1.6 (the v1 plan's
  `{#appendix-a/b}` custom anchors need `attr_list`, never enabled);
  fixed in `mkdocs.yml`.

## 6. Known limitations / notes for Q2+

- The **sanitizer lane scope** at Q1 close: ASan+UBSan (Apple clang lane)
  over the unit suite and the fast regression subset with amg/gamg paths
  exercised; the full A21 matrix re-runs at the next nightly. hypre is
  built un-instrumented (system dependency), so ASan interposition covers
  the frehg2/PETSc boundary, not hypre internals.
- The solve path remains **host MATAIJ, pure-MPI**: threads accelerate
  the Kokkos kernels only (s3 is honest about this — kernel-time gated,
  solve asserted flat). GPU AMG remains deferred (v2 plan §11) pending
  the GPU-validation debt: hypre-CUDA needs device matrices, GAMG +
  MATAIJKOKKOS is the likelier route.
- BoomerAMG operator-complexity gating (g3's ≤ 2.0 assertion) is armed in
  the weak-scaling lane only when `-gw_pc_hypre_boomeramg_print_statistics`
  is enabled; on the reference machine the setup fraction stayed low
  enough that complexity never bound. Real-runner calibration decides
  whether to tighten.
- For Q2+ authors: select `amg` in any new gate case whose grid exceeds
  the benchmark scale; the per-PR g1 lanes guarantee solver invariance,
  so gate tolerances must never be retuned per preconditioner.
