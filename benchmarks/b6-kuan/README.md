# b6-kuan — Kuan et al. tidal saltwater intrusion (lab sandbox)

Two configurations share this directory's `input/` data:

- `b6-kuan-ss.yaml` — no-tide variant (seaward stage constant at 0.15 m);
  golden: `legacy/benchmarks/b6-kuan/out-ss-syncV4/`.
- `b6-kuan-td.yaml` — tidal variant (stage from `input/tide.dat`);
  golden: `legacy/benchmarks/b6-kuan/out-td-syncV4/`.

## Recorded golden/config discrepancy (plan §2.3, §9, risk register)

The archived legacy input file (`legacy/benchmarks/b6-kuan/input`, a
development-side artifact not distributed with this repository) sets
`sync_coupling = 0` (asynchronous/subcycled coupling) and points `foutput` at
an `out-ss-async10V4/` directory — but the golden datasets actually archived
are the **sync** runs (`out-ss-syncV4/`, `out-td-syncV4/`). Both Frehg2
configurations therefore pin `coupling.mode: sync` to match the goldens, not
the committed input.

Additionally, the goldens were generated with the legacy Newton iteration
(`iter_solve = 1`), which the plan drops (§3.2): Frehg2 revalidates b6 under
PCA with the relaxed golden tolerance, and the Kuan experimental
salt-interface data is the primary gate (§9).

## Input provenance

- `input/bathymetry.dat` — verbatim copy of `ex3_kuan_input/bath`
  (68 values, one per j row).
- `input/tide.dat` — `ex3_kuan_input/tide1` converted from alternating
  time/value lines to two-column `t value` (129600 samples, 2 s spacing).
- Initial conditions are constants from the legacy input (`init_eta`,
  `init_wt_abs`, `init_s_surf`, `init_s_subs`); the `*_ic` files in
  `ex3_kuan_input/` correspond to `eta_file/h_file/wc_file = 0` in the
  committed input and are not used.

## The td sea-surface salinity condition (amendment A20)

Every wet surface cell of every td golden output is exactly 35.00 psu,
while the ss golden shows a smooth discharge-diluted near-shore ramp: the
td golden was generated with the legacy wet-cell salinity override active
(the commented-out "Kuan 2019" block, `scalar.c:258-262` — `s_surf := 35`
on wet cells). Physically it idealizes the flooding water as the ocean
reservoir at seawater salinity; without it the fresh beach discharge
dilutes the advancing flood film and the upper saline plume — the
experiment's defining tidal feature — never forms. `b6-kuan-td.yaml`
expresses it as the explicit `sea-surface-salinity` condition (a
whole-tank surface `scalar_value`, i.e. a wet-cell Dirichlet);
`b6-kuan-ss.yaml` keeps the tide-cell-only value matching its own golden.

## P4 gate record (2026-08-21, serial, FI_PROVIDER=tcp)

Both variants PASS (tolerances/b6-kuan.yaml; adjudications in amendment
A20):

| metric | ss | td | allowed |
|---|---|---|---|
| golden head, max abs diff | 0.0121 m | 0.0129 m | max(0.015 m, 5 %) |
| golden salinity, cells outside max(1 psu, 10 %) | 3.51 % | 4.82 % | 6 % |
| interface MAE vs Kuan experiment (model) | 0.0334 m | 0.0389 m | 0.107 m and 1.5x golden |
| interface MAE vs Kuan experiment (golden) | 0.0245 m | 0.0430 m | — |
| tidally averaged salt mass vs golden | +5.5 % | +8.8 % | 10 % |
| salinity range over all outputs | [0, 35] | [0, 35] | [0, 35] |

The salinity outsiders sit in the 1-2-cell fringe band of the ~30 psu
front (wedge toe and intertidal strip); the resolved legacy defects
(amendment A19 — most visibly the uninitialized surface exchange-diffusion
coefficient whose limiter-pinned fixed point the ss golden's near-shore
ramp is) and the PCA-vs-Newton scheme swap displace the steep front by
<= 1-2 cells. Runtime: ~6 min (ss) / ~13 min (td) serial on the reference Apple M3 CPU
(scale to your hardware) — the sync common step rides the 0.5 s dtg cap, not the 0.05 s
legacy fixed step.
