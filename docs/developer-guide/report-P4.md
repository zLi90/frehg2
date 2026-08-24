# Phase 4 development report — handoff for P5

**Audience:** the developer (or agent session) implementing P5. Read
together with `FREHG2_UPGRADE_PLAN.md` (binding spec; P4 added amendments
A18–A20 to its Amendment log), the earlier `report-P0..P3.md` handoffs,
and `dod-P4.md` (the verified Definition of Done). This report explains
what P4 built, what was decided and why, and what P5 must know.

**Status:** P4 complete. Both b6 variants pass the plan §9 gate with the
primary criteria at full strictness — the tidal variant's salt-interface
agreement with the Kuan experiment (MAE 0.0389 m) is *better than the
legacy golden's own* (0.0430 m) — with the golden-comparison tolerances
data-adjudicated per amendment A20. Every P0–P3 criterion remains green,
including bitwise coupled restart with the transport state.

---

## 1. What P4 added and where

- **`frehg::transport`** (`src/transport/`): `ScalarSolver` (setup, scalar
  BCs, step orchestration, end-of-step snapshots, restart refresh),
  `SurfaceTransport.cpp` (scalar_shallowwater), `SubsurfaceTransport.cpp`
  (scalar_groundwater + advective/dispersive face fluxes + ghosts),
  `Dispersion.cpp` (the tensor), `Limiters.hpp` (tvd_superbee). The module
  owns the scalar fields and reads the flow state through driver-wired
  view structs (plan §4 layering: physics libraries never include each
  other).
- **Baroclinic activation in `frehg::gw`**: `Baroclinic.cpp`
  (update_rhovisc + baroclinic_face; named constants
  `kBaroclinicBetaRho/Visc`), `attachScalar`, and the coupled live
  hydrostatic side ghosts (amendment A18) — recomputed per substep from
  live state so restarts stay bitwise.
- **Core:** `HaloExchanger::exchangeWithCorners` now serves 3D fields (the
  A1 protocol unchanged; its "no 3D corners" guard repealed — the
  dispersion cross terms read diagonal subsurface scalars).
- **Driver:** transport construction/wiring (including
  `gw_->attachScalar`), the transport step in all three loops, the
  `/monitor/transport_audit` table, transport output variables
  (`concentration`, `concentration_surface`), checkpoint fields {s_surf,
  s_subs, s_fu_old, s_fv_old, s_dzz_top}, and the restart ordering
  (surface refresh → transport refresh → gw refresh — the baroclinic
  ratios and coupled ghosts read the restored scalar's ghosts).
- **Schema:** `concentration_surface` output variable; surface
  `scalar_value` conditions generalized (wet-cell Dirichlet regions +
  discharge-paired inflow concentrations — amendment A19/A20).
- **Regression harness:** gates `b6` (`--variant ss|td`, golden fields +
  the Kuan-experiment interface MAE + tidally averaged salt mass + bounds),
  `b6-restart`; `tolerances/b6-kuan.yaml`; ctest entries `regression.b6.*`
  (regression_nightly) and `regression.b6_restart` (per-PR, added to
  `ci_build_and_test.sh`); every harness run pins `FI_PROVIDER=tcp`.
- **Docs:** `docs/theory/transport.md` (provenance + preserved quirks +
  A19 resolutions + the scalar budget identity), the groundwater.md
  baroclinic/A18 section, four removed-features rows, parameters.md and
  README updates, the b6 README provenance/gate record.

## 2. Fidelity findings P5 must not re-learn

1. **The td golden ran the legacy wet-cell salinity override.** Every wet
   surface cell of every td golden output is exactly 35.00 psu (the ss
   golden shows a smooth ramp) — the commented-out "Kuan 2019" block
   (scalar.c:258-262) was active in the td golden build. Without it the
   upper saline plume never forms and the td primary gate is unreachable
   (measured: interface MAE 0.114 m vs the 0.107 m cap). The td config's
   `sea-surface-salinity` condition expresses it explicitly; do not remove
   it, and do not add one to ss (its golden and physics are ramp-formed).
2. **The legacy surface exchange-diffusion coefficient is unwritten
   memory** (`smap->dz`, map.c:49; reads at scalar.c:137,145). A
   fresh-heap zero makes the term ±∞ and the limiter pins the film to the
   local extremum — the ss golden's near-shore ramp is that fixed point.
   Frehg2 uses 2·Dzz/dz3d(top) (amendment A19); the residual near-shore
   and wedge-fringe deviation is what the A20 tolerance adjudication
   absorbs (salinity outsiders 3.51 %/4.82 % vs the 6 % allowance, all in
   the 1-2-cell band of the ~30 psu front).
3. **The transport ledger is legacy-lagged and the audit measures it.**
   The surface flux volume uses the *previous* step's flow rates
   (volume_by_flux runs before update_velocity) and the subsurface flux
   volume predates the reallocation/clamp adjustments; under accelerating
   flow Σ s·V drifts by exactly the audited anchor terms. Conservation at
   1e-8/step holds in the quasi-steady/diffusive regimes (where the plan
   gates it); the closure identity holds to rounding everywhere. Any P5
   performance rework of the transport step must preserve the term-by-term
   audit or re-gate b6.
4. **Transport restart state beyond the scalars:** the pre-correction
   flow-rate snapshots (`s_fu_old`/`s_fv_old` — the end-of-step Fu/Fv are
   *pre*-stage-correction values a restart cannot rebuild from the
   checkpointed post-correction velocities) and the carried top-cell
   dispersion coefficient (`s_dzz_top` — its source fluxes are not restart
   state). Both are checkpointed; b6 restart is bitwise at every output.
5. **The superbee stencil needs data one cell beyond the halo.** Legacy
   silently degraded interface faces to upwind (rank-local guards);
   Frehg2 stages shifted far-neighbor fields and corner-exchanges the
   subsurface scalar so decomposed runs reproduce the serial golden
   stencil. The b6 gate itself runs serial (1×68; the goldens are serial).
6. **The dispersion tensor is calibrated on volumetric face fluxes**
   [m³/s] — legacy never divides the face area out (scalar.c:977-994), and
   b6's dispersivities embed that scale. Converting to Darcy velocities
   would silently rescale D by ~Ax/Ay/Az and break the gate.
7. **The coupled hydrostatic side ghost follows the live surface** (A18):
   `(bed − z_c)·r_face + depth·r_face` at plus sides (minus sides leave
   the depth term unscaled — the preserved legacy asymmetry), recomputed
   per substep. The configured `hydrostatic: {eta: ...}` value only
   matters without a surface module (the b6_gw_smoke path).

## 3. Decisions and conventions added in P4

- Surface `scalar_value` semantics (parameters.md): discharge-covered
  members inject the inflow concentration; all other members are wet-cell
  Dirichlet cells (the tide rule generalized to regions). Groundwater-side
  members set the boundary concentration (ghost value) that the inflow
  advection, the limiter bound, and the baroclinic boundary face read.
- The scalar budget identity (docs/theory/transport.md):
  Δ(surface mass) = exchange + source + boundary + adjust + anchor and
  Δ(subsurface mass) = −exchange + boundary + adjust + anchor, every term
  in `/monitor/transport_audit` (cumulative), the P3 "defects are data"
  convention.
- The b6 interface metric: nearest 0.5-crossing per experimental point on
  the legacy plotting script's fine grid (the td columns are non-monotone:
  upper plume + deep wedge; the point set traces both interfaces).
  Identical extraction for model and golden.
- b6 runs serial in the gate (like its goldens); the transport module is
  rank-invariant by construction — measured 1.0e-14/1.3e-14 max relative
  difference at 1/2/4 ranks over 200 strict fixed-dtg coupled+transport
  steps of b6-ss (concentration, both grids, head, eta) — but no b6 rank
  lane was added (the plan requires none; the coupled A14 lanes already
  exercise decomposition).

## 4. The b6 gate (plan §9) — record

| metric | ss | td | allowed |
|---|---|---|---|
| golden head, max abs diff | 0.0121 m | 0.0129 m | max(0.015 m, 5 %) — A20 |
| golden salinity outsiders (max(1 psu, 10 %)) | 3.51 % | 4.82 % | 6 % — A20 |
| interface MAE vs experiment (model) | **0.0334 m** | **0.0389 m** | 0.107 m and 1.5× golden |
| interface MAE vs experiment (golden) | 0.0245 m | 0.0430 m | — |
| tidally averaged salt mass vs golden | +5.5 % | +8.8 % | 10 % |
| salinity range, all outputs, both grids | [0, 35] | [0, 35] | [0, 35] |

Runtimes (serial, OMP_NUM_THREADS=1, FI_PROVIDER=tcp): ss ~6 min, td
~13 min — the A10 sync common step rides the 0.5 s dtg cap instead of the
legacy fixed 0.05 s, so b6 is no longer nightly-class in cost (the label
is kept per plan §11.2; P5 may promote it to the per-PR set).

## 5. What P5 should do first

1. Re-run `scripts/ci_build_and_test.sh` and re-read plan §10 P5; all six
   gates are now green locally (`ctest -L 'regression|regression_nightly'`
   — the b5 envelope runs remain the only hours-long pieces).
2. The sanitizer matrix (P5 deliverable) extends the clang lane to *all*
   test labels; the P4 lane covered unit+mpi+b1/b2 as before.
3. The performance report can reuse `/monitor` and the timer tree; note
   OMP_NUM_THREADS=1 remains fastest at benchmark scales and the b5
   strong-scaling study is specified at 1→8 ranks.
4. The CUDA compile-only lane: the transport module is plain Kokkos
   (no new PETSc surface); the one platform-sensitive piece is the
   `parallel_reduce` with multiple reducers in the transport kernels.
5. Carried risks: CI workflows still unexercised (no remote); the
   libfabric sockets provider wedges MPI_Finalize on this machine —
   FI_PROVIDER=tcp everywhere (A19); the b5 subcycled rain lane remains
   diagnosed-not-regated at the committed dt (A17).

## 6. Open items and risks carried forward

- The b6 element-wise salinity comparison is fringe-band limited by
  construction (a ≤ 1-2-cell front displacement flips ~30 psu cells);
  if the gate is ever re-adjudicated, an interface-aware metric (e.g.
  band-shifted comparison) would be more informative than a tighter
  exceedance fraction.
- The masked-column (ktop > 0) coupled paths (report-P3 item 10) remain
  unexercised through P4; the transport kernels handle kTop generically
  but no gate covers them.
- Transport with subcycled coupling is implemented with the legacy
  structure (last-substep fluxes over the surface dt) but no benchmark
  gates it (b6 is sync; the unit closure tests run it implicitly through
  the coupled toy in sync mode only).
- `initial_conditions.transport` file inputs and multi-scalar (the schema
  reserves the list form) remain unexercised beyond unit coverage of the
  constants.
