# Scalar transport (frehg::transport)

The Frehg2 transport module preserves legacy Frehg's explicit FV scalar
transport (upgrade plan §3.1 item 4). This page is the P4 provenance table
required by plan §11.1 rule 7 and §12: every preserved piece names the
legacy `file:line` it reproduces and the Frehg2 function that implements
it, followed by the legacy quirks the port carries deliberately, the
defined resolutions of the legacy undefined/rank-dependent behaviors
(amendment A19), and the scalar budget identity the transport audit
asserts.

## Transport step

The scalar advances once per surface step, after the flow modules (legacy
`solve.c:112-116`; `ScalarSolver::step` called by the driver after the
coupled/surface/groundwater step): the surface scalar first
(`scalar_shallowwater`, `scalar.c:25-298`; `ScalarSolver::stepSurface`),
then the subsurface scalar (`scalar_groundwater`, `scalar.c:303-498`;
`ScalarSolver::stepSubsurface`). Both use the surface step dt and the flow
state of the completed step; in subcycled coupling the subsurface advection
uses the window's *last* substep's Darcy fluxes over the whole surface dt —
the legacy structure exactly (only sync mode is golden-pinned; b6 runs
sync).

## Provenance

| Piece | Legacy provenance | Frehg2 implementation |
|---|---|---|
| Surface scalar mass ledger s·Vsn + explicit increments | `scalar.c:42, 111-119` | `stepSurface` increment kernel |
| Vsn: the previous velocity phase's volume | `update_subgrid_variable`, `shallowwater.c:1050` | `vsn_` snapshot at the end of each transport step |
| Surface advection: upwind/superbee face values, current flow rates | `scalar.c:45-112` | increment kernel; `tvdSuperbee` (`Limiters.hpp` = `tvd_superbee`, `subroutines.c:541-561`) |
| Flux volume from the *previous* step's flow rates | `volume_by_flux`, `shallowwater.c:996-1032` (called before `update_velocity`, `:112-113`) | `fuOld_`/`fvOld_` snapshots; `vflux_` in the increment kernel |
| Surface diffusion (constant difux/difuy) | `scalar.c:114-119` | increment kernel |
| Seepage scalar exchange (upwind concentration + two-point diffusion) | `scalar.c:121-152` | increment kernel (`sseepage_`); the subsurface subtracts the same mass (`scalar.c:371-375`) |
| Surface limiter: wet-stencil extrema, velocity gate, sentinel guards | `scalar.c:154-179, 206-232` | increment + update kernels |
| Rain/evaporation dilution | `scalar.c:243-254` | update kernel |
| Ghost/tide/dry rules in order | `scalar.c:263-288` | ghost, Dirichlet, dry passes |
| Inflow concentration mass source, limiter release | `scalar.c:180-195` | `surfaceInflow_` lists (paired with discharge conditions) |
| Subsurface scalar mass ledger s·Vgn | `scalar.c:322` | `stepSubsurface` increment kernel (`wcn`·V) |
| Subsurface advection (upwind/superbee on volumetric face fluxes) | `advective_flux`, `scalar.c:501-716` | increment kernel |
| Coupled top face advects through the seepage term only | `scalar.c:700-702` | top-face rules in the increment kernel |
| Dispersion tensor (molecular·θs + α_L/α_T on volumetric fluxes) | `dispersion_tensor`, `scalar.c:958-1003` | `updateDispersionTensor` (`Dispersion.cpp`) |
| Dispersive face fluxes with cross terms | `dispersive_flux`, `scalar.c:719-847` | `dispFluxX/Y/ZLower` device helpers |
| Live y-edge dispersive doublings | `scalar.c:348-350` | increment kernel (`jjp`/`jjm` × 2 at global y edges) |
| Coupled top-interface dispersive flux (wet surface only) | `scalar.c:352-363, 834-845` | top-face `jkm` branch |
| Subsurface flux volume Vgflux | `volume_by_flux_subs`, `groundwater.c:1634-1644` | recomputed in the update kernel |
| Subsurface limiter: conductive-stencil extrema, qtop allowance | `scalar.c:378-434, 452-458` | increment + update kernels |
| Scalar bounds (legacy hard [0, 200] → configurable) | `scalar.c:460-483` | `transport.bounds` (plan §3.2; default [0, ∞), negative clamp at the lower bound) |
| Side/top/bottom scalar ghosts, s_surfkP export | `enforce_scalar_bc`, `scalar.c:851-917` | `enforceSubsurfaceBc` |
| Density/viscosity ratios r_rho/r_visc | `update_rhovisc`, `scalar.c:921-955` | `gw::RichardsSolver::updateBaroclinicFaces` (`Baroclinic.cpp`; constants `kBaroclinicBetaRho/Visc`) |
| Face ratio means and boundary rules | `baroclinic_face`, `groundwater.c:338-415` | `Baroclinic.cpp` (see docs/theory/groundwater.md) |

## Preserved legacy quirks (golden-pinned; do not "fix")

1. **The dispersion tensor runs on volumetric face fluxes** [m³/s], not
   Darcy velocities: legacy never divides the face-area factor out of
   qx/qy/qz (`scalar.c:977-994`), so the configured dispersivities are
   calibrated against face-scaled fluxes. b6's `disp_lon = 0.002` was tuned
   for this scale; the goldens embed it.
2. **The surface flux volume lags the velocity update by one step**:
   `volume_by_flux` runs before `update_velocity` (`shallowwater.c:112-113`),
   so the concentration denominator uses the previous step's flow rates
   while the advective increments use the current ones. Under accelerating
   flow the ledger and the actual volume disagree by the audited
   `surf_anchor` term; in quasi-steady flow (every gated regime) it
   vanishes.
3. **The subsurface ledger re-anchors against the post-reallocation θ**:
   Vgflux predates the reallocation/clamp adjustments; the difference is
   the audited `subs_anchor` term.
4. **The coupled top-interface dispersive flux is one-sided**: the
   subsurface gains `dispersive_flux` at the kM ghost (`scalar.c:362`)
   *and* loses the seepage term's diffusive part, while the surface only
   gains the seepage term — the interface gain has no surface counterpart.
   Measured into `subs_boundary` (negligible at b6's dispersivities).
5. **Surface scalar advects across closed boundary faces**: the closed-edge
   fold keeps the water in, but the boundary-face velocities are live
   computed values and the scalar upwinds across them against the
   zero-gradient ghost. Measured into `surf_boundary`.
6. **The limiter pins single-conductive-neighbor cells**: a cell whose only
   conductive neighbor is one face (e.g. the box-bottom cell) has
   s_min = s_max = that neighbor's value and is dragged to it whenever its
   own update leaves the interval (`scalar.c:460-463`). Measured into
   `subs_adjust`.
7. **The rain dilution ignores the rain's own timing nuances**: legacy
   dilutes by `rain > 0 ? rain : -evap` with the cell's wet area
   (`scalar.c:243-254`); Frehg2 keeps the structure but uses the per-cell
   masked rain (the legacy global rate ignored its own exclusion rows — no
   gated case rains on transport).

## Defined resolutions of legacy undefined/rank-dependent behavior (A19)

1. **The surface exchange-diffusion denominator** `smap->dz` is allocated
   and never written in legacy (`map.c:49`; read only at
   `scalar.c:137,145`) — an uninitialized read of the plan §2.1 hazard
   class. Its garbage value scales the film↔top-cell diffusive exchange
   arbitrarily (a fresh-heap zero makes it ±∞, which the limiter then pins
   to the local extremum — the b6 ss golden's smooth near-shore salinity
   ramp is this pinning's fixed point). Frehg2 uses the top cell's
   thickness (`2 Dzz / dz3d(top)`), the physically defined two-point
   coefficient. The residual b6 deviation sits in the near-shore film and
   the wedge toe and is adjudicated in the gate tolerances (amendment A20).
2. **Superbee far-neighbor guards** were rank-local (`smap->ii != nx-1`
   compares against the *rank's* extent), silently degrading interface
   faces to first-order upwind — a rank-count-dependent stencil. Frehg2
   guards against the global edges and stages shifted far-value fields so
   decomposed runs reproduce the serial golden stencil.
3. **Dispersive cross terms read diagonal neighbors** that legacy never
   exchanged (stale interface corners; unset domain-edge diagonal ghosts
   read fresh-heap zeros). Frehg2 corner-exchanges the subsurface scalar
   (the A1 two-phase protocol, extended to 3D fields) and keeps
   domain-edge diagonal ghosts at their deterministic ghost values.
4. **Dead in-branch boundary doublings** of `dispersive_flux`
   (`scalar.c:786, 831`) are unreachable behind the ghost `actv == 0`
   guards and are dropped (the live `jjp`/`jjm` doublings at
   `scalar.c:348-350` are kept).
5. **The scalar-mass ghost writes** (`sm_surf`/`sm_subs` edge copies and
   exchanges) are dropped: no live code reads them (the top-ghost write at
   `scalar.c:333` even indexes a 3D array with a 2D index).
6. **The uncoupled top-face advection asymmetry** (superbee runs leak the
   donor value through head-condition tops while upwind runs zero the
   inflow side, `scalar.c:653-697`) reduces to the derived branch outcomes,
   implemented directly.

## The scalar budget identity

Per step, with every term measured and reduced into
`/monitor/transport_audit` (cumulative columns; the P3 audit convention —
defects are data, not residuals):

```
Δ(Σ s·dept·A) = exchange + surf_source + surf_boundary + surf_adjust
                + surf_anchor
Δ(Σ s·θ·V)    = -exchange + subs_boundary + subs_adjust + subs_anchor
```

`exchange` is the seepage scalar mass entering the surface (the same value
leaves the subsurface); `surf_source` combines inflow mass and the
prescribed-salinity resets; the `adjust` terms are the limiter/bounds
clips; the `anchor` terms are quirks 2 and 3; the `boundary` terms are
quirks 4 and 5 plus the true open-boundary advective/dispersive transfers.
The closure holds to rounding on every configuration
(`TransportModule.*ClosesTheScalarBudget`); on quasi-steady/diffusive
closed domains the non-boundary terms vanish and Σ s·V is conserved to
better than 1e-8 per step (`TransportModule.*ConservesScalarMass` — the
plan §10 P4 exit criterion).

## Surface scalar boundary semantics (plan §5.6 extension, A19)

A surface `scalar_value` condition's member cells split by carrier: cells
covered by a `discharge` condition inject the inflow concentration
(`current_inflow · dt · s_in / n_cells`, limiter released —
`scalar.c:180-195`); all other members are wet-cell Dirichlet cells (the
legacy tide rule, `scalar.c:275-282`, generalized to any region). The b6
tidal variant prescribes the flooding seawater over the whole tank: its
golden held every wet surface cell at exactly 35 psu at every output —
the legacy wet-cell salinity override (`scalar.c:258-262`, the commented
"Kuan 2019" block) was active in the golden run — expressed in Frehg2 as
an explicit whole-tank condition (amendment A20; b6 README).
