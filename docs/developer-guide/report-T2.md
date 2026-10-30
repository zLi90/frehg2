# Report — Stage-3 phase T2 (selectable inflow closure at prescribed-stage boundaries)

Date: 2026-10-28. Companion: `dod-T2.md`. Plan §6/§13.A; owner decisions O1 (zero_gradient stays the default, bitwise) and O2 (deposit routing tied to `external`). Amendment S3-A5 records the first-green adjudications.

## 1. What landed

| Piece | Where |
|---|---|
| `boundary_conditions[].inflow: {closure: zero_gradient \| external, velocity: {u, v}}` (kind eta only; u/v constant or series, default still) | `Config.*`, `ConfigSchema.cpp` (+cross-field rules), `BoundarySet.*` (per-step component evaluation), resolved-config echo |
| External-reservoir momentum: static reservoir-face masks, per-step reservoir velocity state, `ExtStencil<Ext>` compile-time dispatch — the `Ext = false` instantiation folds to the v2.0 expression trees (verified bitwise by ob6) | `SurfaceSolver.*`, `Momentum.cpp` |
| Conservative entrance: self-advection across a reservoir face is d(u²/2), so a steady entrance pays (v² − u_ext²)/2g exactly — the prototype's 0.5-weight finding generalized to a moving reservoir; cross-advection and uy/vx interpolation read the reservoir state | `Momentum.cpp` (W2.3) |
| Stage-deposit routing (O2): booked in `mass_audit.stage_deposit` under both closures; removed through the reservoir before the velocity update under `external` only — zero_gradient dynamics stay bitwise (the v2.0 kick preserved, now measured) | `WetDry.cpp` (`routeStageDeposits`), `Coupler.cpp`, `Simulation.*` |
| Edge-diffusion audit booking (the prototype's third finding, audit-only) | `SurfaceTransport.cpp` (W2.5) |
| Entrance-head diagnostics: run-record `stage_ke_influx_J` + `stage_entrance_head_ratio`, run-end warning when zero_gradient inflow imports energy without paying the head | `FreeSurface.cpp`, `Simulation.cpp` (W2.6) |
| Gates ob1–ob8, `gate_openbc.py` (the discrete steady chain), `stage3-ob.yaml`, `unit.ob_gates.negative`; closure-mismatch-across-a-face made fatal at setup; "no reservoir face" warning | `benchmarks/ob-openbc/`, `tests/regression/`, `scripts/test_ob_gates.py` |

## 2. The reference design (what makes ob1 self-contained)

The ob1 reference is a **discrete steady face chain** built from the code's own
update formulas (max-of-two face geometry, point-implicit drag with the legacy
double factor, the conservative entrance), marched face-by-face at a trial Q and
shot on Q to hit the downstream stage. Authoring self-checks (§6.3): per-face
residuals 2.8e-16 at the chain's own solution, and the external lane agrees with
an independent continuum Bernoulli+Manning closed form to 0.5 %. First green
(strict solver): the model matches the chain to **1.8e-6 relative** on every
face of both external lanes, with entrance ratios chain-exact (0.9996 / 0.9023).

## 3. The ob5 result (the production-gyre mechanism closed)

The 2026-10 campaign gyre in miniature: a seeded along-shore jet by a constant
stage. Under `external` the domain KE decays monotonically and the stage faces
import 74 J (0.6 % of KE0) over 1200 s; under `zero_gradient` the same faces
import **6671 J — 90×** (the gate's real negative lane and its anti-vacuity
proof). With the walls of the production layout this is the free-energy budget
of report_V1 reproduced in a 3.4 s gate.

## 4. ob8 (b6 under `external`) — the owner-facing record

Both variants **pass the experiment criteria**; the interface MAE *improves*
(ss 0.0325 vs 0.0334; td 0.0381 vs 0.0389 zero_gradient; golden's own 0.0430).
Golden-referenced residuals are recorded, not gated (S3-A5): the legacy golden
embeds the stage-deposit kick (0.90 m³ of the td seepage lands on the stage
row), so `external` moves one beach cell to 0.0196 m against the golden's
0.015 m floor and the tidally-averaged salt mass to 9.88 % (zero_gradient:
8.79 %). These are the data for a future owner decision on flipping the b6
reference; **no tolerance of the zero_gradient b6 gate was changed** (unlike
the 2026-10-26 prototype, whose four loosened bounds T0 reverted).

## 5. S3-A5 adjudications (first green; measured values in `stage3-ob.yaml`)

1. ob1/ob3 run the **strict solver**: the gated identities are momentum-chain
   and ledger identities, and the default rtol 1e-8 leaves solver-class
   residuals (measured: ledger −6.6e-4 m³ default vs 2.5e-8 strict).
2. ob1 harness face indexing corrected (flat index k = face k+1).
3. ob3 `ledger_rel` 1e-12 → 1e-11 (1440-step accumulation at X = 3.0e4 m³;
   12× headroom over the measured 2.5e-8).
4. ob5 `influx_rel` 1e-3 → 2e-2: the seeded slosh exchanges
   pressure-accelerated water through the entrance, which carries real KE
   (measured 6.2e-3); the discriminator is the 10× contrast (measured 90×).
5. ob8 gates the experiment criteria and records the golden residuals (§4).
6. `check_bc_matrix.py` canonical cells extended with the eta-inflow axis
   (120 rows).

## 6. Tier P / Tier N (S3-A1 scope)

Tier P: regression 60/60 at `-j8` (2546 s) + `perf_baseline` solo PASS
(127 s) + unit/mpi/stage3 84/84 (351 s; T2 adds ob1/ob3/ob5 and four unit
gates). Tier N (S3-A1 scope — T2 touched the surface/eta regime): b6.ss +
b6.td + g8 3/3 (697 s at `-j2`); ob8 both variants PASS on the experiment
criteria (664 s); b5 24 h record in `records/T2/tierN-b5record.log`.
Pipeline note: the regression step's tests 1–38 ran on the T2 binary and
the rest on the T3-converted rebuild (the 58-site EdgeTopology
conversion) — green across the boundary plus the step's own
`default_bitwise` PASS is the conversion's default-neutrality evidence.
Details: `stage3/records/T2/timings.md`.

## 7. Deviations

S3-A5 (above). The coupled surface identity in `tests/unit/test_coupling.cpp`
now subtracts `stageDeposit` (exact under both closures: zero under
zero_gradient-with-walls cases, the booked removal under external).
`regression.b6_subcycled`'s identity and tolerance are untouched — under the
default closure the dynamics are bitwise and the new column is informational.
