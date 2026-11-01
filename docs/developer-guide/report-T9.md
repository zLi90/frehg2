# Report — Stage-3 phase T9 (paired water–salt store for thin-layer removals)

Date: 2026-10-31. Companion: `dod-T9.md`. Plan §12A (added by S3-A12,
owner decision O15); amendment S3-A13. Owner commits (O10).

## 1. Why T9 existed

The T8 default flip made `transport.surface_dryout: store` (O4's
salt-only store) the shipped default. On the owner's v30 study, a coupled
2-D tidal beach without ET, it blew up: 72 850 psu surface salinity, 9 039
psu in the soil and 47 CFL warnings within 0.64 of a tide (S3-A12). The
`discard` twin stayed at 35 psu. The mechanism: reallocation vents put
sub-`min_depth` films on a dry seepage face, the wetting limiter deleted
their water unbooked, and the store kept their salt. Each rewetting then
released the accumulated salt into a 1e-5 m film. O15 replaced the store
with a **paired** one: thin-layer removals park their water together with
its salt, rain and evaporation are untouched, and parked water is inert.

## 2. What landed

| Piece | Where |
|---|---|
| Per-path booking of every thin-layer removal (P1 isolated, P2 small, P3 clamp film, P4 thin vents, the sweep) in every mode, as read-only passes beside the pinned kernels; `mass_audit` appends `film_removed`, `film_pending`, `stage_volume` [, `store_volume`, `store_release`]; run-record per-path lines; the `film_removal_rate` output | `WetDry.cpp`, `SurfaceSources.cpp`, `Simulation.cpp` (W9.1) |
| The paired store: P1–P3 park instead of deleting; the end-of-step sweep after the last η change (after `applySeepage` and stage routing) parks every 0 < η − z_b ≤ `min_depth` film; the R5/R6 release; `water_store` output; `swe_water_store` checkpoint | `SurfaceSolver.*`, `Coupler.cpp` (W9.2) |
| Transport: the R8 store-in predicate (no water after the surface step); R3 parks c·ΔV at the window-capped concentration; R4 keeps only the evaporative crust (≤ c·V_evap) and the held-seepage credit (S3-A13 bounds); R5 releases on the combined basis; R9 parks heat with its water (`t_store`, `temperature_audit.surf_store`) | `SurfaceTransport.cpp`, `ScalarSolver.*` (W9.3) |
| Under the store every dry-cell vent deposits (the sweep then parks it with its salt); `vent_thin_m3`. `discard` and legacy unchanged | `Reallocate.cpp` (W9.4) |
| `water_store` schema rule + ConfigTest pin; docs: transport/surface-water theory, output columns and identities, parameters, the agent guide; 10 coverage rows; study-migration addendum | W9.5 |

Without the store (and in `discard`) every field is bitwise T8, and
`default_bitwise` plus the et9 noET lane prove it. The `mass_audit` table
only appends columns, and its T8 prefix is pinned bitwise per case (ws5).

## 3. The gates

| Gate | Result |
|---|---|
| ws0 (8 unit tests + negative battery) | green; the R5 mix exact, a vBasis-only release fails |
| ws1 (per-PR, compressed S3-A12) | T8 red 494 psu → T9 **35.000 psu** everywhere; I2, ledgers, anti-vacuity, CFL ≤ 1.5 |
| ws1n (the owner's repro, nightly) | T8 red 72 854 psu → T9 **35.000 psu** everywhere; water 8.9e-13, salt 5.6e-10 |
| ws2 (ledgers, in the lanes) | store and discard identities close at every row (≤ max(1e-10 X, 5e-11 m³)) |
| ws3 (rain/evaporation invariance) | store ≡ discard on the never-drying basin; ET basin bitwise vs T8; pond dry-out and rewet closed forms exact (1e-15 class) |
| ws4 (path coverage) | P2/P3 books exact, pairing, R9 heat |
| ws-x restart / rank | restart bitwise with 24 parked cells; rank per S3-A13 |
| ws5 (pins) | `default_bitwise` + et9 noET + the mass_audit T8 prefix |

## 4. The W9.1 attribution

The S3-A12 vanished volume, 1.12e-3 m³ exchanged onto dry cells in zoom A,
is **P2's**. The vents land after the solve as 1.0–1.7e-5 m films.
`updateDepth` reports them dry, and the small-depth rule deletes them at
the next solve, 88 % within the same or the following second. With
`film_removed`, `film_pending` and `stage_volume` booked, the coupled
identity closes at every audit row to 1e-11 of the exchange scale. Without
the film terms it misses by 5.5e-3 (`stage3/records/T9/attribution.md`).

## 5. What W9.6 found

- **R4 had to be bounded.** Keeping the whole unpaired remainder whenever
  evaporation took any water turned a dry-front advection remainder into a
  21 psu store on a 10 psu case. The crust is now at most c·V_evap and the
  held keep at most the step's held credit (S3-A13). True evaporative
  dry-out is unchanged.
- **et8d's T5 lane was vacuous on its release half.** On the 1 × 2 strip
  the pond fell below its 5 mm `min_depth` across both cells, and T8
  deleted it every step: the surface held salt in 0 of 241 frames while the
  store banked salt-only. The lane is now one column, where D9 is exercised
  (220 wet frames) and exact on both binaries.
- **Owner item: a transport/solve mismatch at drying faces (pre-existing).**
  When a cell ends a step dry, the velocity limiter zeroes its outflow face
  after the solve has already moved water across it. The transport
  therefore cannot move that water's salt: the water arrives salt-free, and
  the stranded salt beyond R3's cap is booked to `surf_adjust`. At the
  study's `min_depth` (1e-5 m) this is film-scale: ws1n's `surf_adjust` is
  −0.077 psu·m³ versus −0.46 under discard. At a pond-scale `min_depth`
  (et8d's 1 × 2 strip, 5 mm) it removed 24 % of the exfiltrated salt. Every
  earlier mode had the same mismatch but deleted or banked whatever it
  touched. The fix (driving the transport's face volumes from the solve's
  continuity) is a scheme change outside T9.
- Every other store-resolving T5–T8 gate (et9, et9-veg, et-orient,
  et-restart, et8a/b, the restart lanes, b6 ss/td) is green with
  unchanged tolerances.

## 6. W9.7 — the default under the paired semantics: `store` stays

The §2.7 re-evaluation ran every tier with tolerances unchanged on the
build whose default is the paired store: unit/mpi/stage3 153/153, the
regression sweep 60/60, g2, perf, the Tier N nightly set 5/5, the
`stage3_nightly` lanes 10/10 (ws1n, ws_rank, the s3.p1 lanes, ob8) and
the b5 24 h record. **All green: `transport.surface_dryout: store` stays
the shipped default**, now meaning the paired store, and no
`stage3/df/store-T9-report.md` is owed.

The decisive cases:

- **ws1n** (the S3-A12 repro) is at 35.000 psu, against the T8 store's
  72 854 psu.
- **b6 ss/td** stay inside their bands under both modes, with the paired
  store at least as good as `discard`. On ss the tidal salt mass is +5.35 %
  versus +5.47 % and the salinity band 3.59 % versus 3.51 %; td is
  identical to `discard` at the printed precision. The store parked
  2.0–2.9e-3 m³ in these runs, so it was not vacuous.

The v30 study matrix can drop its interim `discard` pin
(`stage3/records/T9/study-migration-addendum.md`).

## 7. Records

Timings and walls: `stage3/records/T9/timings.md`. Red evidence:
`records/T9/red/`; green logs: `records/T9/green/`; the attribution:
`records/T9/attribution.md`. Snapshot `refs/stage3/T9`, phase patch
`records/T9/phase.patch` (vs `refs/stage3/T9-start`).

## 8. Residual limitations

- The drying-face transport/solve mismatch (§5) is pre-existing and
  recorded for the owner. It is film-scale at the study's `min_depth`.
- ws1/ws1n log bounded CFL warnings (≤ 1.25) at R6 return events
  (S3-A13).
- Coupled tidal beaches are trajectory-chaotic in both modes, so
  cell-level rank and perturbation comparisons are not meaningful there.
  The ws gates assert the invariant guarantees instead (source bound, I2,
  ledgers).
