# Phase 2 development report — handoff for P3+

**Audience:** the developer (or agent session) implementing P3 and later
phases. Read together with `FREHG2_UPGRADE_PLAN.md` (binding spec; P2 added
amendments A6–A8 to its Amendment log), `report-P0.md` / `report-P1.md`
(earlier handoffs), and `dod-P2.md` (the verified Definition of Done). This
report explains what P2 built, what was decided and why, and what P3/P4 must
know.

**Status:** P2 complete. All plan §10 P2 exit criteria green: vG unit tests
at 1e-12, b2 and b3 gates pass, b2 rank invariance holds, closed-domain GW
mass balance ≤ 1e-8 per step, the debug-build θ ∈ [θr, θs] assertion is
violation-free on b2/b3, and the early b6 GW-only smoke comparison is
recorded (its numbers are in dod-P2.md §"smoke"). Every P0/P1 criterion
remains green.

---

## 1. What P2 added and where

- **`frehg::gw`** (`src/gw/`): the complete legacy PCA mixed-form Richards
  module. `RichardsSolver` owns the state fields and orchestrates the legacy
  `solve_groundwater` sequence; `VanGenuchten.hpp` (pointwise θ(h)/h(θ)/
  C(h)/K(h)/dK/dθ closures, unit-tested at 1e-12), `Predictor.cpp` (face
  conductivities, the 7-point `LinearSystem("gw_")` COO system, boundary
  ghost heads), `Corrector.cpp` (Darcy face fluxes, the θ update, moisture
  ghosts, the final clamp with the debug bounds assertion),
  `Reallocate.cpp` (the post-allocation consistency restore and
  redistribution — amendment A7), `AdaptiveStep.cpp` (the dtg controller),
  `TerrainMetric.cpp` (bathymetry, ktop masking with the legacy
  partial-cell rules, terrain-following metrics — feeds
  `Grid::buildGlobalIds` the real mask).
- **Driver:** groundwater-only runs march the outer loop on the adaptive
  dtg (legacy `solve.c:43` + `groundwater.c:193` semantics); outputs are
  labeled with the crossed `output_interval` boundary; checkpoints key to
  the crossed whole-second boundary with the exact state time and dtg in
  the header (amendment A6) — restart is bitwise (`regression.b2_restart`).
  Subsurface HDF5 output (`/groundwater/{hydraulic_head, water_content,
  qx, qy, qz}` — fluxes per unit area), geometry-aware
  `/groundwater/zcell/0` (amendment A6), and the `/monitor/gw_mass_audit`
  table (volume, boundary_in, ss_storage, realloc, realloc_dropped, vloss;
  the closure identity is exact by construction).
- **Schema (docs + tests in lockstep):** `groundwater.reallocation_surplus:
  drop | redistribute` (amendment A7); cross-checks pinning `value:
  {gravity}` to `groundwater_bottom` and `value: {hydrostatic}` to
  `groundwater_side`.
- **Regression harness:** gates `b2` (element-wise vs legacy goldens +
  the Warrick wetting-front metric), `b3` (digitized Kirkland contours,
  parsed from the legacy `makeplot.py`, + internal mass balance),
  `b2-restart`, `rank-invariance-b2` (4×4 replication — amendment A8), and
  the **recorded** `b6-gw-smoke` comparison; `compare_h5.py` gained a
  dataset-group parameter and `ascii_golden_to_h5.py` a 3D converter.

## 2. Fidelity findings P3+ must not re-learn

1. **The PCA post-allocation is where groundwater mass goes to die — watch
   it.** The legacy sweep (all benchmarks run `post_allocate = 0`) computes
   the redistribution volumes but has every transfer disabled: the
   over-saturation excess and the saturation-adjacent surplus are silently
   discarded, and the deficit direction *creates* water. On b2 this is a
   0.3 %-of-inflow effect the golden embeds (so the b2 gate pins the `drop`
   default); on b3 it ate ~25 % of the injected volume, starved Kirkland's
   perched saturation bulb, and made the h = 0 contour unreachable —
   b3 requires `reallocation_surplus: redistribute` (amendment A7 has the
   full story). **b5 and b6 must choose at their gates** and the
   `realloc`/`realloc_dropped` audit columns are the first thing to look at
   when a coupled budget drifts.
2. **Nonzero Ss + dry flat-retention soil + the consistency restore is a
   runaway feedback.** The predictor pulls the dry cell's head toward a wet
   neighbor, the corrector's Ss factor bleeds θ, h(θ) maps the loss to an
   even lower head, repeat — b3's Glendale flank walked to θr and
   h ~ −1e14 within 1200 s of simulation. Kirkland is storage-free
   physics, so b3 runs `specific_storage: 0`. b6's sand is steep-retention
   (bounded mismatch) but has Ss = 1e-5 — if P4's ss variant drifts dry
   anywhere, suspect this first (`docs/theory/groundwater.md`, "preserved
   quirks").
3. **A fully saturated, Ss = 0, all-Neumann groundwater system is
   singular** (C(h > 0) = 0 leaves no diagonal), in legacy exactly as here.
   The coupled P3 system always has wet-cell Dirichlet tops, but
   groundwater-only tests/cases need either an anchor (head BC) or Ss > 0.
4. **The corrector re-derives face conductivities from the predicted head
   with its own Ks selection** (x/y faces: the minus cell's Ks for both
   flanks; z faces: the lower cell's). At b3's sand/clay interfaces the
   predictor and corrector genuinely use different face K — preserved
   verbatim; don't "unify" it.
5. **Legacy is rank-variant in the subsurface; frehg2 is not.** One-sided
   interface conductivities/fluxes/geometry were replaced by two-sided
   evaluations from exchanged neighbor state (identical single-rank
   results). The post-allocation sweep is deterministic per column. If P3
   sees rank-dependent coupled results, the surface side or the coupling is
   the suspect, not `frehg::gw`.
6. **Masked columns (ktop > 0) keep a sealed top face in the predictor**
   (the legacy seal-pass ordering) while the corrector's boundary flux
   flows — b2/b3 have ktop = 0 everywhere, so **b5 is the first case that
   exercises this**; its gate adjudicates whether the quirk stands.
7. **Side-flux sign is now "positive into the domain"** (the legacy
   predictor and corrector disagreed at the y+ face; nothing exercised it).
   b6's inland freshwater inflow was converted to +2.77e-5 in both b6
   configs — P4 should not "fix" it back.
8. **The b3 contour units question (P0 report) is adjudicated:** the
   digitized h = 0 / h = −400 points match the simulated contours in the
   model's own head units — no rescaling (benchmarks/b3-kirkland/README.md).

## 3. Decisions and conventions added in P2

- Groundwater BC staging: polygon conditions rasterize to per-column /
  per-edge marker fields (`topCode/botCode/sideCode*` + values refreshed
  per step), decoding the legacy `bctype_GW` codes; defaults are no-flux.
  Prescribed-head tops use the saturated face conductivity; free drainage
  (`value: gravity`) is bottom-only; hydrostatic heads are side-only
  (schema cross-checks + tests).
- Face-field layout: `kx/qx(j, i, k)` is the x+ face of cell (j, i, k) with
  slot i = 0 the west boundary/interface face; `kzF/qzF` carry nz + 1
  planes (plane k = face above cell k; plane nz = the bottom face).
  Internal fluxes keep the legacy face-area factor [m³/s]; outputs divide
  by the face area.
- Restart state for GW = {h, wc} + the dtg scalar; everything else is
  re-derived at step start. Coupled P3 restart = SWE fields + these.
- The gw solver prefix is `gw_` (CG, bjacobi/icc, rtol 1e-8, atol 1e-14,
  maxit 1000 per plan §5.1); strict rank-invariance mode is
  `-gw_pc_type jacobi -gw_ksp_rtol 1e-13 -gw_ksp_atol 1e-16`.
- P4 hooks in place: the density/viscosity face-ratio fields (`r_rho*`,
  `r_visc*`) are multiplied in every formula and held at 1; `baroclinic_face`
  averaging lands with the transport module that produces non-unit values.
  `Vg/Vgn/Vgflux` cell-volume bookkeeping is *not* stored — all three are
  derivable from wc/wcn when transport needs them.

## 4. Performance

b2 single-rank (23,445 adaptive steps, 100 cells): 11 s wall
(`OMP_NUM_THREADS=1`; the P1 kernel-launch-bound findings apply unchanged —
launch count per gw step is ~25 kernels + 1 PETSc solve). b3 single-rank
(8,873 steps, 1,500 cells): 53 s. The b2 4×4 rank-invariance replication at
strict tolerances is the slowest P2 artifact (minutes per lane); the
regression harness keeps it in the `mpi`-labeled set.

## 5. What P3 should do first

1. Re-run `scripts/ci_build_and_test.sh`; re-read plan §10 P3 and §9 b5,
   plus amendments A2 (outflow kind; b5's provisional −10 m stage-sink
   outlet must be revisited) and A6–A8.
2. The coupling slots are prepared: `Simulation` fatals on
   surface+groundwater together — replace with the legacy `solve.c:25-166`
   sequence (surface free-surface solve → groundwater (lockstep or
   subcycled) → surface velocity update). `SurfaceSolver`'s two-phase split
   exists for exactly this.
3. The coupled top boundary replaces `topCode` for wet columns: Dirichlet
   ghost head = depth with saturated face Ksz, the infiltration limit
   (available surface water), and seepage → η per plan §5.7 (the
   `sim_shallowwater == 1` branches deliberately left out of
   `Predictor.cpp`/`Corrector.cpp`/`Reallocate.cpp` are marked "arrives
   with P3" at each site; `groundwater_flux:842-868` has the seepage
   bookkeeping to port).
4. Remember the elevation datum: `TerrainMetric` applies the same
   −min(bath) lift as the SWE module, so surface depth and subsurface
   heads already share a frame; the coupled ktop for the *surface* system
   under bathymetry masking is `ktop == nz` columns inactive (Grid handles
   both ids).
5. b5 will exercise: masked columns (finding 6), `exchangeWithCorners` on
   2D decompositions, the subcycled mode, and the reallocation-surplus
   choice. Measure before freezing tolerances.

## 6. Open items and risks carried forward

- CI workflows remain unexercised (no remote); the new regression jobs
  (b2/b3/b2-restart/rank-invariance-b2) join the fetch-goldens story from
  P1's report.
- The b6 GW-only smoke numbers (dod-P2.md) bound the PCA-vs-Newton +
  no-coupling gap at quasi-steady state. The comparison is intentionally
  loose context for P4's tolerance freeze, not a gate.
- MPICH on this machine intermittently wedges `mpiexec` at exit (OFI
  sockets; see report-P1). It bit twice during P2's rank-invariance runs —
  kill and rerun; the harness output is unaffected.
- The `groundwater.timestep.dq_grow/dq_shrink` defaults (0.01/0.02) and
  `courant_max` follow the legacy constants; no P2 benchmark stresses the
  Courant cap (b2's controller is dq-driven throughout). b5's coarser
  columns may.
