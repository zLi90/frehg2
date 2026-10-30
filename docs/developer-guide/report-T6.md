# Report — Stage-3 phase T6 (vegetation and transpiration)

Date: 2026-10-29. Companion: `dod-T6.md`. Plan §10/§13.C; owner decision O7 (transpiration continues while inundated). Amendment S3-A9 (seven entries) records the adjudications and a major pre-existing defect this phase's restart lane uncovered.

## 1. What landed

| Piece | Where |
|---|---|
| `vegetation` block: type table (LAI, extinction, r_l, root depth/profile, Feddes h1–h4 with the wetland no-h₁h₂ form, T7 albedo), id map (0 = bare), optional LAI raster; §13.C rules (needs GW + the ET block; consistent subsurface update with salinity; Feddes ordering; name/map referential checks) | `Config.*`, `ConfigSchema.cpp` (W6.1) |
| `beerLambert` (exactly 1 at LAI 0), `canopyResistance`, `feddes` (piecewise, O7 wetland form), `rootFractionIntegral` (uniform/linear closed forms) | `atm/EtFormulas.hpp` (W6.2) |
| `setupVegetation` (host, static): per-cell τ/r_c/Feddes staging, closed-form root fractions on the terrain-following/masked column measured from the ktop face, shallow-column renormalisation with a warning count, all-bare-map inert mode; `updateTranspiration` per substep: T_pot = (1 − τ)·E_x(r_c), the Feddes- and supply-capped volumetric sink S_k (NO inundation test — O7), demand/applied books, the sink dtg window; τ scales the soil channel in-kernel and the open-water channel via an attached factor field | `gw/Evapotranspiration.cpp`, `swe/SurfaceSources.cpp` (W6.3) |
| The sink in the Richards step via **separate post-kernels** (S3-A9.7): the predictor adds −dtg·rRho·S·V to the assembled RHS, the corrector subtracts (dtg·rRho·S·V)/coeff from θ and books the ss-audit delta — the pinned kernels stay byte-identical | `Predictor.cpp`, `Corrector.cpp` (W6.4) |
| `sinkWin_` window; transport `vgflux −= sinkWin` with the exact per-cell factor (combined vap + sink on ET top cells) — water out, salt stays | `RichardsSolver.cpp`, `SubsurfaceTransport.cpp` (W6.5) |
| `gw_mass_audit` + `transpiration`/`transpiration_unrealised` (closure ΔV = b_in − ss + re − vloss − T), `/atmosphere` transpiration rates, the **Field3 IntervalRates leg** (S3-A7's deferral closed) with `groundwater: root_uptake` [1/s], checkpoint `acc_root_uptake`, run-record `transpiration_m3` | `Simulation.*`, `io/IntervalRates.*` (W6.6) |
| Docs: theory sections, parameters, input-data map ids, the §6d marsh recipe, 6 coverage rows | W6.7 |

## 2. The S3-A9 find: a coupled-transport restart defect older than Stage 3

The et-veg-restart lane blew up (c → 1e30 in one step), and the controls
peeled every T6 feature away: the SAME failure reproduces on plain et9
with no vegetation, no ET block, under `legacy` AND `consistent` updates,
periodic or walled — **the base coupled transport restart has been broken
since it existed**. Mechanism: a subcycled window can take zero
groundwater substeps right after a restart (the restored gw clock lags),
and the transport then reads the last substep's per-substep state — the
volumetric fluxes, the face conductivities (the limiter's conductive-
stencil gates) and `dtgLast` — which the checkpoint never carried. As
allocation zeros, every limiter gate stayed closed and the first
transport step clipped the entire subsurface scalar to its empty-window
seed (subs_adjust = −M). b6-restart never caught it because its dtg ≤ dt
guarantees a substep every window. **Fix**: checkpoints with transport
now carry `qx/qy/qzf_last`, `kx/ky/kzf_last` and the scalar `dtg_last`,
reproducing the documented stale-flux legacy quirk bitwise; the minimal
control and all 14 restart lanes (10 ctest + et/dg/pd harness) are green
on the fixed binary, and runs without transport keep the prior checkpoint
layout byte-identical.

## 3. The fma-contraction rule, now a doctrine

Touching the pinned predictor-assembly / θ-update kernels — even with a
guarded, never-taken branch — changed the compiler's contraction and
broke `default_bitwise` on b2 (the S3-A8.1 class, second measurement).
The W6.4 sink therefore applies in separate post-kernels that only launch
when vegetation is active; the default path is byte-identical **by
construction**, not by hope, and `default_bitwise` plus the et9 noET pin
(1271 datasets) re-verified it on the final binary.

## 4. Physics results

The partition behaves: under a τ = 9e-14 canopy the soil channel shuts
off (E_soil 14 orders below T) and the depletion spreads over the root
zone (top-cell share 0.051 vs the bare contrast's 0.53) — et5's three
lanes. Root fractions are exact on stretched, masked, zone-cutting
columns (et6, 3e-10 against the `root_uptake` output — the Field3
registry leg's first exercise). Feddes stress and the supply cap keep
the books closed to the last bit (et7-t), transpiration excludes salt
exactly (et8-t: c·θ conserved per cell, salinity rising with rootFrac at
correlation 0.9999931), and the marsh band transpires through inundation
(O7 observed on the tidal flat) while bare cells carry exactly zero.

## 5. Tier P / Tier N (S3-A1 scope)

T6 touched groundwater + transport + the coupling/output layers: Tier N =
the b5 24 h record + b6 ss/td + g5 ×2 + g8 alongside the anchored Tier P
sweep, g2, unit/mpi/stage3 and perf solo. Results and walls:
`stage3/records/T6/timings.md`.
