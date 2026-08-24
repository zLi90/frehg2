# Phase 3 development report — handoff for P4+

**Audience:** the developer (or agent session) implementing P4 and later
phases. Read together with `FREHG2_UPGRADE_PLAN.md` (binding spec; P3 added
amendments A9–A17 to its Amendment log), `report-P0.md` / `report-P1.md` /
`report-P2.md` (earlier handoffs), and `dod-P3.md` (the verified Definition
of Done). This report explains what P3 built, what was decided and why, and
what P4/P5 must know.

**Status:** P3 complete. The b5 envelope gate passed its rain/sync lane
on the owner-directed shortened horizon (amendment A15) and the owner then
concluded the battery early (amendment A17 — the reference envelope's own
uncertainty does not justify the remaining wall-hours; the full-horizon
runs remain in the suite as regression_nightly). Every other plan §10 P3
exit criterion is green at full strictness: coupled mass closure holds at
1e-10 per step, b1/b2 re-gate green with unchanged tolerance files,
coupled restart is bitwise at every output label, and the b5-short
rank-invariance lanes pass per amendment A14. Every P0/P1/P2 criterion
remains green.

---

## 1. What P3 added and where

- **`frehg::coupling`** (`src/coupling/`): `Coupler` (step orchestration,
  synchronization modes, restartable adaptive state, the coupled audit)
  and `Exchange.cpp` (window staging, the `subsurface_source` port, audit
  reductions). The coupler owns the exchange fields (seepage accumulator,
  window infiltration budget, deposit tally, the qss observable) and wires
  device views of them — plus the surface eta/depth — into the groundwater
  module (`gw::GwCoupling`), which mutates the surface state exactly where
  legacy did (the bounce-back and the reallocation vent write eta/depth
  mid-substep).
- **Coupled branches in `frehg::gw`**: the coupled top-face conductivity,
  matrix/RHS, corrector flux (with the A9 limit), the per-column
  capacity/supply/dry classification (A13), the top bookkeeping kernel
  (bounce-back, seepage accumulation, evaporation correction), the
  consistency-restore and vent branches in `Reallocate.cpp`. GwStepAudit
  gains the coupled terms (exchanged, bounce, vent, evaporated).
- **Driver:** `runCoupledLoop` — sync marches on the common adaptive step
  (A10; crossed-boundary output/checkpoint semantics of A6), subcycled
  keeps the fixed surface dt; the mass_audit table gains the `seepage`
  column; coupled checkpoints add {seep_accum, qss}, scalars
  {dtg, tgw_lag, dt_surface}, and the surface adds etan (both feed the
  restart refresh's edge-slot reconstruction, §2.7); the "seepage" output
  variable is the applied rate. The crossed-boundary output/checkpoint
  labels now advance past *every* boundary a step crossed — the old
  one-interval bump drifted the labels behind the state whenever
  interval < dt (latent since A6; no gated configuration hits it, a P3
  restart diagnostic with interval = 1 s did).
- **Mesh:** `domain.terrain_layers: uniform` (A12) — the terrain-parallel
  slab of the b5 reference; the legacy scaled-terrain wedge remains the
  default (b6).
- **Schema:** `terrain_layers` (+ cross-check requiring follow_terrain),
  the coupled `groundwater_top`-must-be-flux cross-check (A11).
- **Regression harness:** gates `b5` (envelope, --scenario × --coupling),
  `b5-restart`, `rank-invariance-b5` (A14 lanes);
  `tolerances/b5-vcatchment.yaml`; the four envelope runs carry the
  `regression_nightly` ctest label (hours each at 4 ranks — plan §11.2
  keeps b5/b6 out of the per-PR set).
- **Benchmarks:** `b5-vcatchment.yaml` reauthored (raster-orientation flip
  of the polygon/monitor coordinates, uniform terrain slab, stage-reservoir
  outlet, wetting_face_depth 1e-3, relaxed controller knobs — all
  adjudications in `benchmarks/b5-vcatchment/README.md`) and the
  `b5-vcatchment-norain.yaml` scenario added.

## 2. Fidelity findings P4+ must not re-learn

1. **The legacy coupled exchange is valid only in the capacity-limited
   regime.** b6 (deep ponds) lives there; b5's thin rain films do not.
   Applying the wet Dirichlet to a supply-starved column lets the
   consistency restore manufacture water at Ks-front speed (measured:
   16 m³ in one 7.8 s step). The A13 classification (from SERGHEI, the
   reference implementation and b5 envelope member) is now load-bearing:
   wet columns split capacity/supply on the saturated-Darcy demand
   estimate vs the available depth; supply columns feed the predictor and
   corrector the remaining depth as a flux, bounded by the top cell's
   pore room up to 0.9999 θs. **P4 must not "simplify" the classification
   away** — b6's tidal waterline will wander through both regimes.
2. **The §5.7 resolution extends to the infiltration limit** (A9): no
   porosity factor anywhere in the exchange. The b6 goldens were produced
   with the legacy wcs-factor limit; whatever deviation that causes at
   the moving waterline is part of the P4 tolerance budget (the plan
   already assigns §5.7 deviations to b6's loose gate).
3. **Legacy's rate-based seepage application over-draws when dtg < dt**;
   the volume accumulator (A9) is exact in both modes and reduces to the
   legacy product in sync mode. The dry-cell hold rule survives on the
   accumulator and is checkpointed (`seep_accum`).
4. **Sync coupling is the common adaptive step** (A10): legacy feeds the
   adapted dtg back into the surface dt. b6's committed dt = 0.05 s with
   dt_max = 0.5 is therefore the *initial* step there, not the horizon
   cost — plan the b6 runtimes accordingly.
5. **Coupled trajectories are chaotic at threshold density** (A14): the
   controller forks on single-cell θs-boundary rounding; even fixed-dtg
   runs drift to 2e-1 within 600 s. Bulk volumes are the rank-robust
   observables (subsurface volume 3e-7, exchange 3e-5 relative across
   1/2/4 ranks). Design any P4 cross-run comparison on bulk or
   quasi-steady observables — b6's gate (quasi-steady state + tidally
   averaged salt mass) already is.
6. **Three surface defects were exposed at b5's scale and fixed** (A13;
   `docs/theory/surface-water.md`): the frame-dependent dry-cell rows
   (replaced by the zero-depth continuity closure — this is also what
   makes wetting conservative and the matrix SPD), the never-written drag
   coefficient on dry cells (min_depth floor), and the two-edge-only
   volume audit (now four edges). b1/b4 re-gate green with unchanged
   tolerances (b1 worst achieved/allowed 0.231 — the P1 value exactly).
7. **The stage-boundary velocity correction's west/south writes land in
   halo slots** (`uu(j,0)`, `vv(0,i)`) that the §7 interior-only
   checkpoint cannot carry — yet they feed the uy/vx interpolation, so a
   zero-gradient refresh reconstruction seeds an O(1e-9)-class first-step
   difference that the outlet's point-implicit drag limit cycle amplifies
   into a persistent phase flip (measured: eta 1e-3, vv 3.7e-2, never
   decaying). The refresh therefore *reconstructs* exactly those writes
   from the restored eta^n (new checkpoint field `etan`), the completed
   step's dt (new coupled scalar `dt_surface`), and flow rates recomputed
   from the restored post-limiter velocities; the east/north branches
   write interior slots the checkpoint already restored and must stay
   skipped (re-running them over post-correction state perturbs them —
   b1's tide row pins this). Also: the refresh must not re-prescribe
   stages (it would wipe the rain film on stage cells) and must not
   recompute cflx/cfly (they belong to the completed step). With all
   three, coupled restart is bitwise at every output label; b1/b2
   restarts stay bitwise.
8. **The point-implicit drag linearization surge/stalls at drag numbers
   above ~1** (D = 1/(1 + ½ dt C_D |u^n|/h) is a lagged linearization; its
   fixed-point iteration turns into a 2-cycle when the equilibrium drag
   number is large). Consequence: a face from a wet cell into a bed-level
   stage reservoir — a waterfall face, whose local gradient is the whole
   stored head — chokes ~60× for any face n ≳ 0.35, impounding an
   artificial lake that releases catastrophically when the forcing stops
   (measured on b5: +51 % discharge spike at rain end). Any future case
   with a stage-reservoir outfall needs a small-n lubricant strip (b5 uses
   n = 0.1; SERGHEI's raster used 0.001 for its own scheme), or the
   linearization needs a self-consistent solve — deferred, since every
   benchmark regime is either low-drag-number (b1-b4) or lubricated (b5).
   Amendment A16.
9. **The below-bed clamp is a measurable legacy mass defect at b5's
   scale** (~1.6 % of rain; effectively zero in b1-b4): the mass_audit
   table now measures it (`clamped` column, amendment A16), so budget
   closure is exact bookkeeping with the defect reported as data.
10. **The masked-column (ktop > 0) predictor seal quirk remains
   unexercised**: b5's reference geometry required the uniform slab
   (ktop = 0 everywhere). It will stay unexercised through v1.0 unless a
   case with a regular mesh over non-flat bathymetry couples.

## 3. Decisions and conventions added in P3

- The coupled top boundary is coupler-owned (A11): wet ⇒ Dirichlet depth,
  dry ⇒ seepage face (legacy bctype_GW[5] = 2 with qtop = 0 — b6's
  committed configuration); configured `groundwater_top` flux conditions
  feed qtop (guards + evaporation correction preserved); head kinds are
  rejected in coupled runs.
- Coupled rain falls on every non-excluded cell (A11); the legacy
  wet-cell-only branch discarded rain over dry land.
- The b5 outlet is a bed-level stage reservoir (kind eta, 0.0). The A2
  `outflow` kind deadlocks on flat channel ends (ghost stage = cell
  stage, ghost velocity = quiescent pool → a 3 m lake filled the domain);
  it remains correct for sloped outlets (b4). Discharge is measured from
  the mass-audit table; its eta-BC accounting is exact by construction.
- `CouplingAudit` + the extended `GwStepAudit` carry every exchange term;
  the closure identity is written out in `docs/theory/exchange-flux.md`
  and asserted by `CoupledModule.*ClosesMassBudget` at 1e-10 per step in
  both modes.
- The "seepage" output/checkpoint variable is the surface-applied rate of
  the last coupled step [m/s] (the legacy qss observable).

## 4. The b5 gate (plan §9)

**Verdict (amendments A15/A17):** rain/sync PASS on the 24 h shortened
horizon; battery concluded early by the owner after that lane. Numbers of
the gated run (4 ranks, the committed configuration exactly):

| metric | model vs envelope | tolerance |
|---|---|---|
| discharge, band | 0/48 reference times outside | all inside ±0.1·peak |
| discharge, peak | 1.13e3 vs 1.10e3 m³/h (+2.6 %) | ±15 % |
| discharge, integral (2.7–23.8 h) | +2.4 % | ±10 % |
| ponding, band | 0/54 reference times outside | all inside ±0.1·peak |
| ponding, peak | 596 vs 605 m³ (−1.5 %) | ±15 % |
| ponding, integral (2.3–23.8 h) | −2.5 % | ±10 % |
| surface budget residual | −10 m³ of 1.11e4 rain | 1e-3 |
| subsurface budget residual | 3.9e-8 m³ | 1e-3 |
| below-bed clamp creation (reported) | 180 m³ = 1.6 % of rain | measured, not gated |

rain/subcycled measured FAIL at the committed time.dt = 5 s (discharge
peak +16.3 %, integral +10.4 %; ponding and budgets pass) — the fixed
surface step is the whole story (CFL ~10 excursions, clamp creation 8 %
of rain); the config records the diagnosis and the ~2 s recommendation.
norain was stopped mid-flight by the A17 decision.

**What the adjudication established** (full detail in the benchmark
README and amendment A16): the outlet must be a free outfall — with the
strip at channel roughness the point-implicit drag limit cycle chokes the
outfall ~60×, impounds an artificial 0.92 m lake, and releases it at rain
end as a +51 % discharge spike; the flow-surface roughness was
reconstructed against the envelope's own storage–discharge arithmetic by
measured kinematic scaling (anchors 337/540 m³ at two configs, h ∝ n^0.6
predicted 543/600, landed 540/596 vs the 605 m³ envelope mean), landing on
channel n = 6.0 — essentially the raster's own 6.26, vindicated once the
outlet artifact was removed — hillslope 1.45, outfall strip 0.1.

Runtime facts (macOS arm64, 4 ranks, OMP_NUM_THREADS=1): the rain-window
common step settles near 1.8 s under the relaxed controller knobs
(courant_max 10, dq 0.05/0.1 — the legacy defaults pin it near 0.4 s, see
the config README); the four envelope runs are nightly-class (hours each).
The `regression_nightly` label separates them from the per-PR set.

## 5. What P4 should do first

1. Re-run `scripts/ci_build_and_test.sh`; re-read plan §10 P4 and §9 b6,
   plus amendments A9–A17 (A9's limit and A13's classification both act
   at b6's moving waterline; A16's drag-limit-cycle finding matters for
   any stage-reservoir outfall).
2. The transport module slots into the driver loop after the velocity
   update (legacy solve.c:112-116); the coupler already aggregates the
   per-window gw audits P4's scalar exchange will extend (the legacy
   sseepage uses qss, scalar.c:132-144 — the `Coupler::seepageRate` field
   is exactly that observable).
3. The r_rho/r_visc face-ratio hooks are still held at 1 and multiplied
   everywhere (P2); `baroclinic_face` lands with transport. The b6
   hydrostatic side BC uses a *constant* reference stage (P0 decision) —
   decide at the b6 gate whether the tidal variant needs it to follow the
   tide series (legacy adds the live surface depth, enforce_head_bc:777).
4. b6 runs sync (golden syncV4); with A10 that means the whole model
   marches on the adaptive dtg from 0.05 s within [0.05, 0.5] — measure
   the ss-variant runtime early.
5. The b6 gW-only smoke numbers (dod-P2.md) plus the A9/A13 deviations
   size the golden-comparison risk; the experiment remains the primary
   gate (plan §9).

## 6. Open items and risks carried forward

- CI workflows remain unexercised (no remote); the b5 nightly label needs
  a scheduled runner with ~a working day of budget, or scenario sharding.
- The four-edge audit slightly changes the *measured* boundary_outflow of
  cases with open west/south edges; b1–b4 gates are unaffected (verified),
  but any later tolerance tightening should re-baseline from current
  numbers.
- MPICH's intermittent wedge (reports P1/P2) reappeared once during P3:
  an overnight envelope battery stalled mid-run for 8 h with all ranks
  idle. The gate's `--reuse-output` flag exists exactly for this — the
  partial HDF5 output stays readable and the envelope checks clip to the
  covered window, so a wedged run's completed hours are not lost.
- The default-lane b5 rank-invariance bounds (A14) are data-derived with
  ~3-30× headroom; the harness prints achieved values on every run so
  they can be revisited.
- The b5 rain/subcycled lane fails at the committed time.dt = 5 s and is
  diagnosed but not re-gated (A17): re-running it at dt ~2 s (and the two
  norain lanes) is the first step if the b5 envelope is ever
  re-adjudicated — e.g. with fresh reference data replacing the digitized
  curves.
