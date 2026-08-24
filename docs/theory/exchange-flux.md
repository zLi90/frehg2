# Surface–subsurface exchange (frehg::coupling)

The Frehg2 coupler preserves legacy Frehg's flux-based coupling (upgrade
plan §3.1 item 3) and carries the plan §5.7 resolution of the legacy
seepage-unit defect. This page is the P3 provenance table required by plan
§11.1 rule 7 and §12: every preserved piece names the legacy `file:line` it
reproduces and the Frehg2 function that implements it, followed by the §5.7
conservation derivation, the coupled-regime corrections adjudicated at the
b5 gate (amendment A13), and the audit identity the coupled mass test
asserts.

## Coupled step

One coupled step (legacy `solve`, `solve.c:25-166`;
`coupling::Coupler::step`):

1. surface free-surface phase (`solve_shallowwater`, `solve.c:51-52`);
2. depth refresh so the subsurface sees the post-solve depth
   (`update_depth` at the end of `solve_shallowwater`, `shallowwater.c:87`);
3. groundwater phase — **sync**: one lockstep subsurface step of the common
   dt (`solve.c:68-69`); **subcycled**: adaptive-dtg substeps while a whole
   dtg fits into the surface step (`solve.c:72-75`), the subsurface clock
   carrying its lag across steps;
4. seepage application to η (`subsurface_source`, `shallowwater.c:639-683`;
   `Coupler::applySeepage`);
5. surface velocity phase (`shallowwater_velocity`, `solve.c:97-98`).

Sync mode marches **both** modules on the common adaptive step: legacy
feeds the adapted dtg back into the surface dt (`solve_groundwater:193`),
starting from `time.dt` (`solve.c:37`) and clamped to
`groundwater.timestep` bounds. Frehg2 adapts between whole steps so one
consistent dt covers both phases; legacy swapped dt mid-step between the
eta solve and the velocity update, a §2.1-class inconsistency (amendment
A10).

## Exchange provenance

| Piece | Legacy provenance | Frehg2 implementation |
|---|---|---|
| Wet-column Dirichlet top: ghost head = surface depth | `enforce_head_bc`, `groundwater.c:795-799` | the depth enters the top-face formulas directly (`Predictor.cpp` coupled branches; no k-ghost plane exists) |
| Wet-column predictor face conductivity 0.5 (Ksz + K(h)) | `compute_K_face`, `groundwater.c:276-284` (legacy's scalar `param->Ksz` generalizes to the top cell's own Ksz) | `RichardsSolver::computeFaceConductivity` top rule |
| Wet-column corrector face conductivity = saturated Ksz (the preserved predictor/corrector asymmetry) | `darcy_flux`, `subroutines.c:102-103` | `RichardsSolver::computeFluxes` mode-1 branch |
| Top-face Darcy flux with the depth as ghost head | `darcy_flux`, `subroutines.c:71-75, 87-108` | `computeFluxes` coupled top branch |
| Infiltration limited to the available surface water | `darcy_flux`, `subroutines.c:164-169` (**without** the spurious porosity factor — §5.7, amendment A9) | `computeFluxes` mode-1 limit against the window budget |
| Dry-column seepage face: infiltration prohibited, exfiltration free, configured qtop with moisture guards | `darcy_flux`, `subroutines.c:171-180`; matrix/rhs `groundwater.c:641-663` | `computeFluxes` mode-0 branch; `assembleSystem` coupled dry branch |
| Coupled default dry-top = the bctype_GW[5] = 2 semantics | b6's committed legacy input (the only coupled legacy configuration); amendment A11 | coupler-owned top; configured `groundwater_top` conditions must be kind flux (schema cross-check) |
| Saturation bounce-back (downward flux into a roomless top cell under a dry surface returns to η at once) | `groundwater_flux`, `groundwater.c:844-851` | `RichardsSolver::applyCoupledTopBookkeeping` |
| Seepage accumulation across substeps with the dry-cell hold | `groundwater_flux`, `groundwater.c:852-857`; `reset_seepage` in `subsurface_source` `shallowwater.c:657-675` | `applyCoupledTopBookkeeping` + `Coupler::applySeepage` (volume form — amendment A9) |
| Evaporation correction (positive qtop leaving a moist, surface-dry column does not pond) | `groundwater_flux`, `groundwater.c:858-863` | `applyCoupledTopBookkeeping` |
| Seepage → η with the below-bed clamp | `subsurface_source`, `shallowwater.c:652-658` | `Coupler::applySeepage` (clamp remainder audited) |
| Wet surface counts as saturated contact in the consistency restore | `check_adj_sat`, `groundwater.c:1061` | `Reallocate.cpp` phase-1 coupled branch |
| Reallocation vent onto the surface (any amount when wet; only above min_depth when dry, smaller vents discarded) | `allocate_send`, `groundwater.c:1185-1207` | `Reallocate.cpp` phase-2 coupled branch (drops audited) |
| Coupled rain on every non-excluded cell | amendment A11 (legacy `evaprain` coupled branch `shallowwater.c:602-612` rained only on wet cells, discarding rain over dry land; no golden pins it) | `SurfaceSolver::evapRain`, unchanged from surface-only runs |

## The §5.7 seepage-unit resolution (binding)

The Darcy flux q at the coupled top face is a volumetric flux per unit plan
area [m/s]. The corrector updates subsurface storage as Δθ = q dtg / dz per
cell, i.e. Δ(Σ θ dz) = q dtg per unit area; volume conservation across the
interface therefore requires the surface to change by exactly Δη = ∓q dtg —
**no porosity factor** (the `· wcs` variant in the legacy comments at
`shallowwater.c:643-646` double-counts porosity and breaks the coupled
budget by a factor θs).

The same argument governs the infiltration limit (amendment A9): "flux
cannot exceed surface water available" means |q| dtg ≤ available depth. The
legacy bound `vseep = |q| dtg wcs` (`subroutines.c:166-168`) let the
subsurface withdraw up to depth/θs while the surface lost at most the
depth (the below-bed clamp discarded the difference) — creating water
whenever the limit bound. With the factor dropped, the limit-bound exchange
hands the subsurface exactly the pond and lands the surface on the bed
(`CoupledModule.InfiltrationLimitedToAvailableWater`).

**Subcycled conservation (amendment A9).** Legacy accumulated the exchange
as a *rate* (Σ q_i / Az over substeps) and applied it as rate × surface dt —
consistent only when dtg ≥ dt (the stale rate re-applies each surface step)
and over-drawing by the subcycle count when dtg < dt. The Frehg2 accumulator
carries the exchanged *volume per unit area* (Σ q_i dtg_i / Az); the surface
applies the accumulator itself. Sync mode (dtg = dt) reproduces the legacy
product exactly. The infiltration limit debits a per-window available-depth
budget seeded from the post-solve depth and credited by every coupled
deposit — the live legacy `dept` the limit read, made consistent across
substeps.

## Supply-limited exchange (amendment A13)

Legacy applies the wet-column Dirichlet unconditionally. When the
saturated-Darcy demand exceeds the available surface water (a thin rain
film over conductive soil — b5's entire hillslope), the predictor still
pulls the top-cell head to the surface depth, the consistency restore then
fills θ to θ(h) ≈ θs, and the reallocation vents the corrector's overfill
back to the surface: the infiltration front races at Ks speed funded by
*manufactured* water (measured on b5: 16 m³ created in one 7.8 s step). No
legacy benchmark reaches this regime — b6, the only coupled legacy case, is
always deep-ponded.

Frehg2 adopts the reference implementation's classification (SERGHEI
`GwBC.h:277-288`, the same author's evolution of this exchange; SERGHEI's
results are the b5 envelope member):

- **mode 1 (capacity-limited)**: the saturated-conductivity demand estimate
  `q_est = 2 Ks (h_top − depth)/dz − Ks` can be funded
  (−q_est dtg ≤ available): the legacy Dirichlet path, with the A9 limit as
  the backstop;
- **mode 2 (supply-limited)**: the predictor and corrector both take the
  remaining surface water as a top flux (SERGHEI `GwBC.h:605-609, :408`),
  bounded by the top cell's pore room up to the legacy effective-saturation
  mark 0.9999 θs — one large adaptive substep can neither over-pressurize a
  filling column nor land it exactly on the singular C(h) = 0 state
  (SERGHEI's CFL-sized surface steps keep the same formula bounded
  implicitly; the plan §5.2 adaptive dtg does not);
- **mode 0 (dry)**: the legacy seepage face.

## Coupled-regime surface corrections (amendment A13)

Two legacy surface defects sit below the b1/b4 tolerance floor but are
fatal at b5's scale; both are §2.1-class fixes recorded in
`docs/theory/surface-water.md`:

- **dry-cell rows**: legacy held dry cells at η with live off-diagonal legs
  whenever a face velocity was nonzero (`shallowwater.c:353-360`) — a
  frame-dependent pull on the neighbors' absolute stage (its accidental
  wetting mechanism). Frehg2 closes dry cells as zero-depth continuity rows
  (free-surface area = cell area at the bed); the higher-of-two-bottoms
  face depths guarantee live legs only carry flux *into* a dry cell, so
  wetting is implicit, conservative, symmetric, and translation-invariant.
- **drag on dry cells**: legacy wrote the drag coefficient only where
  Vs > 0 (`shallowwater.c:984`), leaving faces against persistently dry
  cells undragged (b5's outlet cells reached km/s). Dry cells evaluate the
  same law at the min_depth floor; wet-cell values are bit-identical.

The volume audit also counts the explicit fluxes of all four domain edges
(P1 audited east/north only; b5's west-edge outflow outlet is the first
case where the omission is material).

## The coupled audit identity

Per surface step, with every term measured
(`CouplingAudit`/`GwStepAudit`/`SurfaceStepAudit`; asserted at 1e-10 of the
total volume by `CoupledModule.*ClosesMassBudget`, plan §5.7/§8.1):

```
Δ(V_surf + V_gw + V_in_transit) =
    rain − evaporation − boundary_outflow + bc_inflow + clamped  (surface)
  + (gw_boundary_in + exchanged + vent)                        (external gw)
  + bounce − discarded − evaporated_from_seepage
  − ss_storage + realloc − vloss                               (gw internal)
```

`V_in_transit` is the held seepage (the dry-cell hold rule);
`exchanged`/`vent`/`bounce` are the three surface-deposit paths;
`discarded` is the seepage-application clamp remainder (infiltration that
would draw eta below the bed); `clamped` is the free-surface phase's
measured below-bed clamp creation (amendment A16). The identity is exact by
construction — physical non-conservation (the consistency restore, the
reallocation drops, the clamps) appears in the measured terms, never in the
residual.
