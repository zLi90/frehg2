# Groundwater module (frehg::gw)

The Frehg2 groundwater module preserves legacy Frehg's mass-conservative
predictor–corrector ("PCA") mixed-form Richards solver (upgrade plan §3.1
item 2; the Newton and modified-Picard schemes are dropped, §3.2). This page
is the P2 provenance table required by plan §11.1 rule 7 and §12: every
preserved equation names the legacy `file:line` it reproduces and the Frehg2
function that implements it, followed by the preserved quirks, the
deliberate deviations, and the correctness fixes.

## Scheme summary

Per subsurface step of size dtg (legacy `solve_groundwater`,
`groundwater.c:57-199`; `gw::RichardsSolver::step`):

1. save h<sup>n</sup> and θ<sup>n</sup>, evaluate C(h<sup>n</sup>), refresh
   halos and boundary ghosts;
2. **predictor** — arithmetic-mean face conductivities from h<sup>n</sup>;
   one global 7-point linear system for h with storage
   (C(h<sup>n</sup>) + Ss θ<sup>n</sup>/θs) (SPD; PETSc CG, options prefix
   `gw_`, rtol 1e-8 as legacy `SetRTCAccuracy`);
3. **corrector** — face conductivities from the predicted h, Darcy face
   fluxes, θ from the flux divergence with the Ss compressibility factor;
4. **post-allocation** — the θ/h consistency restore and gradient-directed
   redistribution of over-saturation excess (always on, plan Appendix A);
5. final θ clamp to [θr, θs] with volume-loss accounting;
6. the adaptive dtg controller (grow 1.25 / shrink 0.75 on the max flux
   change, ∂K/∂θ Courant cap, [dt_min, dt_max] clamp, cross-rank min).

Groundwater-only runs march the outer loop on the adaptive dtg itself
(legacy feeds dtg back into `param->dt`, `groundwater.c:193` with
`solve.c:43`); outputs are labeled with the crossed multiple of
`output_interval`, holding the state within one dtg of the label — the
legacy trigger semantics (`solve.c:121-126`).

## Equation provenance

| Piece | Legacy provenance | Frehg2 implementation |
|---|---|---|
| θ(h) van Genuchten retention with the **live aev cutoff** | `compute_wch`, `subroutines.c:280-312` (`aev` outside the `use_mvg` guard at `:292`) | `VanGenuchten.hpp` (`waterContentFromHead`) |
| h(θ) inversion with the θr + 1e-7 floor | `compute_hwc`, `subroutines.c:314-346` | `VanGenuchten.hpp` (`headFromWaterContent`) |
| C(h) specific capacity (cutoff at h > 0, **not** aev) | `compute_ch`, `subroutines.c:349-376` | `VanGenuchten.hpp` (`capacityFromHead`) |
| K(h) Mualem conductivity (cutoff at h > 0, capped at Ks) | `compute_K`, `subroutines.c:378-405` | `VanGenuchten.hpp` (`conductivityFromHead`) |
| dK/dθ for the Courant cap (0.9999 θs limiter) | `compute_dKdwc`, `subroutines.c:451-483` | `VanGenuchten.hpp` (`conductivityDerivative`) |
| Arithmetic-mean face conductivities, each flank with its own state and Ks; head-BC edges take the interior K; no-flux edges sealed | `compute_K_face`, `groundwater.c:202-292` | `Predictor.cpp` (`computeFaceConductivity`) |
| Saturated top face under a prescribed head (face K = Ksz) | `groundwater.c:287`, `subroutines.c:102-105` | `Predictor.cpp` / `Corrector.cpp` top-face rules |
| Coupled top boundary (wet Dirichlet / seepage face / supply-limited flux) | `groundwater.c:276-284, 552-557, 635-663, 795-799, 835-868`, `subroutines.c:71-75, 102-105, 160-180` | the P3 coupled branches; the full provenance table is in [exchange-flux.md](exchange-flux.md) |
| `use_full3d = 0`: lateral K sealed at unsaturated cells | `compute_K_face:295-307` | `Predictor.cpp` (face-centric form; see deviations) |
| Seal of inactive-cell faces, ordered last | `compute_K_face:310-316` | `Predictor.cpp` (`gw_face_k_seal_*` kernels, same ordering) |
| Matrix coefficients (storage + six legs; head edges doubled at half spacing; flux edges removed from the diagonal) | `groundwater_mat_coeff`, `groundwater.c:503-569` | `Predictor.cpp` (`assembleSystem`) |
| Right-hand side (storage, gravity, terrain terms, boundary legs) | `groundwater_rhs`, `groundwater.c:572-682` | `Predictor.cpp` (`assembleSystem`) |
| Linear solve (legacy LASPack CG+SSOR at rtol 1e-8 → PETSc CG at the same rtol, plan §5.1) | `solve_groundwater_system`, `groundwater.c:729-748` | `core::LinearSystem("gw_")`, `Predictor.cpp` (`fillAndSolve`) |
| Boundary ghost heads (zero-gradient, prescribed, hydrostatic sides) | `enforce_head_bc`, `groundwater.c:751-813` | `Predictor.cpp` (`enforceHeadBc`) |
| Darcy face fluxes (x/y from the minus cell's Ks, z from the lower cell's; boundary overrides; the saturated-top and free-drainage rules) | `darcy_flux`, `subroutines.c:27-195`; `groundwater_flux`, `groundwater.c:816-841` | `Corrector.cpp` (`computeFluxes`) |
| θ from flux divergence with the Ss compressibility factor | `update_water_content`, `groundwater.c:903-930` | `Corrector.cpp` (`updateWaterContent`) |
| Moisture ghosts (saturated at head-BC sides) | `enforce_moisture_bc`, `groundwater.c:933-956` | `Corrector.cpp` (`enforceMoistureBc`) |
| Post-allocation: consistency restore, saturated-adjacency test, gradient split, vertical send walks | `reallocate_water_content:959-1041`, `check_adj_sat:1044-1071`, `check_head_gradient:1074-1150`, `allocate_send:1153-1291` | `Reallocate.cpp` (`reallocateWaterContent`) |
| Final θ clamp with volume-loss accounting | `solve_groundwater:163-189` | `Corrector.cpp` (`finalizeWaterContent`) |
| Adaptive dtg controller | `adaptive_time_step`, `groundwater.c:1651-1727` | `AdaptiveStep.cpp` (`adaptTimeStep`) |
| Terrain-following mesh and metrics; regular-mesh partial cells; ktop masking | `build_subsurf_map`, `map.c:196-616` | `TerrainMetric.cpp` |
| Initial conditions (moisture / head / water table) | `ic_subsurface`, `initialize.c:998-1162` | `RichardsSolver.cpp` (`applyInitialConditions`) |
| Density/viscosity face ratios (hooks held at 1 until P4) | `baroclinic_face`, `groundwater.c:338-414`; `update_rhovisc`, `scalar.c:921` | `r_rho*/r_visc*` fields in every formula |

## Preserved load-bearing quirks

- **The aev cutoff applies to θ(h) only.** C(h) and K(h) cut off at h > 0.
  Between aev and 0 a cell is saturated in θ but still on the Mualem K
  curve — live in b1/b2/b6 and pinned by `test_vangenuchten`.
- **The corrector's face conductivities are not the predictor's.**
  `darcy_flux` re-evaluates K from the predicted head with its own Ks
  selection: x/y faces use the *minus* cell's saturated conductivity for
  both flanks; z faces use the *lower* cell's. Heterogeneous-soil faces
  (b3's sand/clay interfaces) genuinely differ between the two phases.
- **Top cells are excluded from the adaptive controller** (`istop` guard,
  `groundwater.c:1672`), and the per-cell flux sums divide by the cell's own
  face areas, not the per-face areas (`:1674-1675`).
- **The consistency restore discards or creates water by design.** Cells
  adjacent to saturation take θ := θ(h_predicted) whether that gains or
  loses volume relative to the conserved corrector value
  (`groundwater.c:1034`); the audit's `realloc` column makes the moved
  volume observable, and `groundwater.reallocation_surplus: redistribute`
  (amendment A7) re-enables the legacy send path for the loss direction.
- **Masked columns (ktop > 0) keep a sealed top face in the predictor.**
  The legacy seal loop runs after the top-face rule and zeroes it again;
  preserved verbatim. (b5 was expected to adjudicate this at P3, but its
  reference geometry required the uniform terrain slab — amendment A12 —
  under which every column keeps ktop = 0; the quirk remains unexercised
  by any benchmark.)
- **PCA compressibility/consistency feedback on dry, flat-retention soils.**
  With Ss > 0, a persistent mismatch between the predicted head and the
  restored h(θ) pumps θ downward (the Ss term of the corrector) — on b3's
  Glendale clay this runs away (P2 finding). Kirkland is an unsaturated
  problem and its SERGHEI-style reference carries no storage term, so b3
  runs with Ss = 0. Watch the b6 sand at P4 (steep retention, so the
  mismatch stays bounded, but the mechanism exists whenever Ss > 0).

## Deviations from legacy (all rank-invariance or correctness driven)

| Deviation | Legacy behavior | Rationale |
|---|---|---|
| Interface faces evaluate two-sided from exchanged neighbor state | One-sided values that disagreed across ranks (`compute_K_face:238-262`, `darcy_flux` ghost collapse, interface geometry copied from own cells `map.c:509-527`) | The two ranks sharing a face must see the same conductivity/flux (plan §8.2 rank invariance); single-rank results identical |
| `use_full3d` sealing is face-centric and also applied to the corrector fluxes | Cell-ordered writes (races under parallel sweep); the corrector re-derived lateral K and leaked | Face-centric = the single-rank legacy semantics; a leaking corrector contradicts the documented 1D-column mode (§2.1). Unobservable in any benchmark (b2 is a single column) |
| Post-allocation phases: classify (parallel) then per-column sequential walks; fresh pore room | Flat-order interleaved sweep with pre-corrector room (`check_room` before `update_water_content`), decomposition-dependent and able to overfill past θs | Deterministic, rank-invariant, and the θ ≤ θs invariant the P2 exit criterion asserts (amendment A7) |
| Lateral post-allocation fractions dropped and audited | Order-dependent single-neighbor delivery, remainder dropped silently (retry disabled at `:983`) | Column-local walks per plan §5.2 ("the Reallocate sweep is sequential per column"); the drop is observable in `/monitor/gw_mass_audit` |
| Side-flux sign: positive = into the domain, both edges, both phases | Predictor treated positive y+ flux as outflow while the corrector treated it as inflow (`groundwater_rhs:612` vs `darcy_flux:155`); never exercised (b6's qym is the only user and its magnitude enters both variants) | One convention, documented in parameters.md; b6's configs carry the converted value |
| Output qx/qy/qz are per-unit-area fluxes [m/s] | Legacy wrote the face-area-multiplied fluxes [m³/s] | User-facing units; the b2 gate compares head/θ only |
| zcell output is the full per-cell `(j·NX+i)·NZ+k` dataset | (P0 wrote one uniform level table) | b3's `makeplot.py` pins the §7 layout |

## Legacy defects fixed on port (§2.1 hazard class; none exercised by any benchmark configuration)

| Defect | Legacy location | Fix |
|---|---|---|
| Prescribed-flux bottom RHS term missing dtg/dz on qbot | `groundwater_rhs`, `groundwater.c:631` | Dimensionally consistent term mirroring the top-flux form |
| Prescribed-head bottom RHS enters with the wrong sign | `groundwater.c:633` (`+= Gzp·h` instead of `-=`) | Dirichlet sign, matching the y-edge treatment |
| Bottom free drainage tests the *top* boundary code | `darcy_flux`, `subroutines.c:190` (`bctype_GW[5]` in the bottom branch) | Free drainage keyed to the bottom condition (`value: gravity`, target `groundwater_bottom`) |
| Ghost soil parameters read uninitialized | `Ksx/vga/vgn/wcs/wcr` ghosts never written (`ic_subsurface` fills `n3ci` of `n3ct`) | Domain-edge ghosts mirror the interior; interface ghosts halo-exchanged |
| `check_head_gradient` adds a cell index to a thickness | `groundwater.c:1110` (`dz3d[ii] + icjckM[ii]`) | 0.5 (dz(k) + dz(k−1)) |
| Prescribed-head bottom face spacing uses uninitialized ghost dz | `groundwater_mat_coeff:525` at k = nz−1 | Half-cell spacing (the Dirichlet-at-face convention the top uses) |
| Masked-column (ktop > 0) top-face flux stored the below-face value | `groundwater_flux:839` + `darcy_flux:71` (`kk[ic] == 0` test) | The top face of *any* column takes the back-z rule with its boundary overrides |
| Bottom ghost head copied from the column *top* | `enforce_head_bc:793` (`kPou = h[kMin]`) | Ghost used only by the prescribed-head bottom (set to that value); otherwise unused |

## Baroclinic activation (P4)

With `modules.transport` and `groundwater.density_coupling.enabled` on, the
driver wires the transport scalar into the module
(`RichardsSolver::attachScalar`) and every step opens by evaluating the
legacy density/viscosity ratios (`update_rhovisc`, `scalar.c:928-930`):
r_rho = 1 + 0.000744 s and r_visc = 1/(1 + 0.0022 s), the named constants
`kBaroclinicBetaRho`/`kBaroclinicBetaVisc` (plan §3.1 item 5), over every
cell including ghosts — the ghost scalar carries the transport module's
side Dirichlet values, so a saline head boundary densifies its boundary
face. The face means and the boundary rules reproduce `baroclinic_face`
(`groundwater.c:338-415`; `Baroclinic.cpp`):

- interior/interface faces: arithmetic mean of the two flanks;
- minus-side domain edges: mean of the interior cell and the ghost, with a
  head condition taking the ghost's own ratio one-sided (`:407-410`);
- plus-side domain edges: 1 unless a head condition takes the ghost's ratio
  one-sided (`:361-364`); legacy left non-head plus edges at 1 because the
  `actv` pair guard fails there — preserved;
- the face above each column's top cell: mean of the top cell and the
  surface scalar in coupled runs (`:387-389` with the `enforce_scalar_bc`
  ghost), one-sided otherwise; the bottom boundary face stays at 1;
- legacy defined the head-BC rules only for y faces (its x sides had no
  head conditions); the x rules mirror the y rules per the P2 x-side
  generalization.

The predictor matrix/RHS, the corrector fluxes, and the θ update already
multiplied every hook (P2 "present-but-inactive"); activation changes no
formula. The hydrostatic side ghosts also become live (amendment A18): in
coupled runs the ghost head follows the legacy live form
(`enforce_head_bc:769-788`) — `(bed − z_c) r_face + depth r_face` at plus
sides, `(bed − z_c) r_face + depth` at minus sides (the legacy asymmetry,
`:777` vs `:787`) — recomputed from live state at every substep start so
restarts reproduce it bitwise (legacy refreshed it post-solve, carrying
mid-step depth one step; the half-step freshness difference is a §2.1-class
consistency choice recorded in A18). Uncoupled runs keep the configured
reference form `(value − z_c) r_face`.

## Mass audit

`/monitor/gw_mass_audit` records, per step (cumulative, reduced over ranks):
`volume` (Σ θ dx dy dz), `boundary_in` (net inflow through domain-boundary
faces), `ss_storage` (water moved into compressible storage by the Ss factor
of the θ update), `realloc` (net θ change by the post-allocation step),
`realloc_dropped` (the part of `realloc` that is pure loss), and `vloss`
(the final clamp). The identity

    ΔV = boundary_in − ss_storage + realloc − vloss

closes to rounding every step by construction (`GwModule.AuditIdentity…`
asserts ≤ 1e-12), and a fully closed domain conserves volume to better than
1e-8 relative per step (plan §10 P2; `GwModule.ClosedColumnConservesMass`).
