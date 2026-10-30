# Report — Stage-3 phase T5 (Feature C core: coupled evaporation with salt conservation)

Date: 2026-10-29. Companion: `dod-T5.md`. Plan §9/§13.C; owner decisions O4 (the salt store), O6 (PM deferred to T7), O7 (recorded for T6). Amendment S3-A8 (nine entries) records the first-green adjudications and case redesigns.

## 1. What landed

| Piece | Where |
|---|---|
| The `evapotranspiration` block (potential bulk/prescribed; penman_monteith validates-fails until T7), `atmosphere.aerodynamic_resistance` (liu2006/fao56/log_profile), `soil.types[].theta_fc`, `transport.surface_dryout` (discard/store — the T8 DF constant; token per S3-A8.10), the `/atmosphere` output section, all §13.C cross-field rules, resolved-config echo, r1 round trip | `Config.*`, `ConfigSchema.cpp` (W5.1) |
| `atm/EtFormulas.hpp`: `aeroResistance` (liu2006 forwards to the v2.0 formula), `potentialBulk` (≡ the v2.0 tree at liu2006/r_s = 0, asserted bitwise), `limiterBeta` (Lee & Pielke 1992), `supplyCap` (both signs; the 0.9999 θ_s condensation room), `thetaFcFromRetention` | W5.2 |
| The open-water channel: per-cell removal min(E_w·dt, depth), the wet-fraction clock (`wetFrac`, `wetClock`, the step-start snapshot the coupler reads) | `swe/SurfaceSources.cpp`, `SurfaceSolver.*` (W5.3) |
| The soil channel (`gw/Evapotranspiration.cpp`): per-substep potential × limiter × **dry fraction** (the wet-clock split — time-basis complementarity by construction), supply-capped; a coupled wet column stands down; demand/applied/unrealised books | W5.4 |
| Guard-free `qVap` on both sides of the substep (predictor RHS + corrector flux — the D1 fix) and the exact accumulator strip + channel-volume audit (D2) | `gw/Predictor.cpp`, `gw/Corrector.cpp` (W5.5) |
| The window vapour volume (`windowVapourVolume`) feeding the W5.8 exact top-cell concentration factor; ET top faces pass water but no scalar (the cauchy relation; replaces the D5 flag on ET columns) | W5.6, `transport/SubsurfaceTransport.cpp` (W5.8) |
| The surface salt store (O4/D6): dry-cell capture (positive mass only — a robustness find), rewet release after limiter/bounds/dilution, the `surf_store` audit column, `salt_store` output, `s_store` checkpoint; per-cell D8 dilution `(rainAdded − evapRemoved)·A` | `transport/SurfaceTransport.cpp`, `ScalarSolver.*` (W5.7) |
| Held-seepage salt (D9): `exfNew` (this step's new exfiltration) debits the top cell as the water leaves; held shares land in the store; the applied backlog carries no second copy | `coupling/Exchange.cpp`, `Coupler.*` (W5.9) |
| Audits (`soil_evaporation`, `soil_evaporation_unrealised` as a boundary_in sub-account), `/atmosphere` interval rates + clocks (the T4 registry), checkpoint fields, run-record closure lines | `Simulation.*`, (W5.10) |
| Docs: `theory/evapotranspiration.md` + parameters/configuration/output/transport/atmosphere/case-setup §6c + 9 coverage rows | W5.11 |

## 2. The defect inventory, closed by measurement

Every D-item the plan targeted is now a measured red-vs-green pair:
**D1** — et1b: the v2.0 guard-skipped predictor leaves a **5872 m**
first-output head error on a saturated coupled column; the guard-free
path's ledger closes at 1.7e-13. **D2** — the audit books the channel's
own volume (et4: boundary_in = −soil_evaporation exactly). **D4** — et7:
θ ≥ θ_r with unrealised = demand − applied to the last bit. **D5/D8** —
et8a/et3. **D6** — et8b/c: `discard` loses exactly 35·h₀·A; `store`
conserves to 4e-15 through dry-out and rewet. **D9** — et8d:
`surf_adjust ≡ 0` through a hold-apply seepage cycle. **The window-lag
double count** — et3's v2.0 lane removes volume at 0.26 of E_w·t (the
double-counted supply drains the pond early); the ET path holds the
combined-rate identity at 0 (bulk) / 2.9e-15 (prescribed) with
`wet_time + soil_time ≤ t` exact.

## 3. Default-neutrality, twice over

`default_bitwise` is green, and et9's noET lane — the full coupled
tidal-flat case WITHOUT the ET block — is **bitwise against its T5-start
pin (1271 datasets)** recorded on the T4 binary before any T5 source
change: every T5 edit (predictor/corrector branches, transport update
passes, coupler fields) is proven inert when the block is absent, by
execution on a case that exercises all the touched code paths.

## 4. The trajectory-stability doctrine, applied again

The S3-A6 twin-control tool earned its keep twice: (1) et9's rank lane —
the 1e-12 bed twin reproduces the rank deltas exactly (tracer 137.9 vs
137.8), so the tide front is chaotic and the strict identities moved to a
stabilized base; (2) the orient/rank base itself took three authorings
(drying rim → front-chaotic; 0.04 m bed steps → shoreline flips under
exactly the 1e-2 legacy mixed class; 0.1 m steps + low-K soil → stable).
On the stable base, rank invariance holds at **machine precision**
(θ 1.3e-15, store 0.0) and the battery's transpose lanes are **exact**
(tracer 0.0, store 0.0, θ 0.0) — T5 adds no asymmetric arithmetic. The
mixed lanes are recorded, not gated (S3-A8.9): the measured 5.9e-2
depth / 20-psu shoreline deltas are the DOCUMENTED legacy
live-east/north-edge-face asymmetry (the id lane is provably static at
1.5e-4 while the rotated pond drains a shoreline row), a pre-existing
surface property gated by its own swe-orient battery.

## 5. Remaining scope notes

The vG dry singularity (et4's 3-day horizon) and saturated-regime
Ss-released salt (et8a's Ss > 0 anchor) are pre-existing domain limits of
the Richards/θ-basis formulations, documented in the case headers — the
S3-A4 applicability class, unchanged by T5. PM (O6) and vegetation (O7)
arrive in T7/T6.

## 6. Tier P / Tier N (S3-A1 scope)

T5 touched surface + groundwater + transport + coupling: Tier N = the b5
24 h record + b6 ss/td + g5 ×2 + g8, alongside the anchored Tier P sweep,
g2, unit/mpi/stage3 (adds ~14 et entries) and perf solo. Results and
walls: `stage3/records/T5/timings.md`.
