# Report — Stage-3 phase T7 (Penman–Monteith potential)

Date: 2026-10-29. Companion: `dod-T7.md`. Plan §11; amendment S3-A10 (six entries).

## 1. What landed

| Piece | Where |
|---|---|
| `atmosphere.shortwave_incoming` (incident R_s [W/m²], distinct from the Q5 **absorbed** `shortwave`, whose meaning is unchanged); `MetSample::shortwaveIncoming`; `evapotranspiration.energy` (albedo_water 0.06, ground_heat_fraction water 0 / soil 0.1 / canopy 0.1); `soil.types[].albedo` and `vegetation.types[].albedo` with has-flags; the W5.1 "not available before T7" rule removed; the §13.C PM rules (radiation pair required, PM + temperature rejected, alpha_barton-with-PM rejected, per-type albedo required under PM — S3-A10.5) | `Config.*`, `ConfigSchema.cpp`, `MetForcing.hpp` (W7.1) |
| `penmanMonteithRate(T_air, P, q_a, R_s, LW_in, albedo, gFrac, r_a, r_s)` — the q-form sharing the bulk chain's Tetens/λ/ρ_a/c_pa/ε — and the analytic `dqsatdT(T, P)` | new `atm/PenmanMonteith.hpp` (W7.2) |
| Dispatch per channel: open water (r_s = 0, `albedo_water`, water G), soil (r_s = 0 then the beta/none limiter, the **top soil type's** albedo per column, soil G), canopy (r_s = r_c, the **vegetation type's** albedo per column, canopy G); staged per-column albedo fields beside θ_fc/r_c | `SurfaceSolver.cpp`, `Evapotranspiration.cpp` (W7.3) |
| Docs: the PM theory section (derivation, commensurability, the FAO convention difference), parameters, 4 coverage rows | W7.4 |

## 2. The FAO-56 anchor

The owner's FAO-56 PDF was transcribed — all 65 printed intermediates of
the Chapter-4 worked examples with page references — into
`benchmarks/et-evapotranspiration/reference/fao56-examples.yaml`
(Example 18, Uccle daily, ET₀ 3.88 ≈ 3.9 mm/d; Example 19, N'Diaye
hourly, 0.63 and 0.00 mm/h; the plan's recalled values confirmed from
source). `gate_et.fao56_daily()/fao56_hourly()` (FAO's own e-form chain)
reproduce every intermediate to the published precision, with each value
perturbation-tested (S3-A10.2 records the precision-floor subtlety on 11
coarse-printed keys and the full-precision-u2 requirement). Against that
chain the model's q-form measures its convention gap:
**δ_conv = 0.016166 mm/h (2.58 %)** — inside the plan's expected 1–3 %
class — which sets et10's tolerance (× 1.2 by S3-A10.1).

## 3. The gates

- **et0 part 3** (unit): dq_sat/dT analytic vs central difference
  (1e-8); the equilibrium limit **exact** (and its perturbed-s
  negative); the aerodynamic limit exact and G-independent; the
  high-resistance limit (S3-A10.4: the suppression ratio is
  (s+1)/(s+1+r_s/r_a), so the 1e-12 bound is a math-limit check at
  r_a = 0.1); commensurability with the bulk mode at the energy-balance
  surface temperature — **absolute-gap order 2.03** (gated ≥ 1.8),
  relative order 0.98 recorded (S3-A10.3: the secant-vs-tangent error is
  first order relative; both rates vanish linearly, so their gap is
  second order absolute).
- **et10** (model level): the FAO reference grass column (LAI 2.88,
  extinction 11 ⇒ τ ≈ 1.7e-14, r_c = 69.4 s/m, albedo 0.23,
  G = 0.1·R_n, `fao56` r_a at z_ref 2, radiation closing FAO's budget
  exactly) under et6's exactness conditions. T_act over the hour equals
  `gate_et.pm_q_form()` at **1.7e-16** (bound 1e-10), the rate
  0.64310 mm/h sits 0.0131 from FAO's published 0.63 (bound 0.0244),
  and the albedo 0.23 → 0.20 negative misses the identity by 3.1e-2.
- **Schema pins** (unit): PM accepted with the radiation pair; rejected
  without it, with `modules.temperature`, and with undeclared soil /
  vegetation albedos.

## 4. Bitwise discipline

The PM dispatch lives in the host-side potential evaluations and inside
the T5/T6 feature-gated kernels (`gw_et_update`, the transpiration
sink) — never in the pinned Richards kernels. `default_bitwise` is green
on the final binary, and the T5/T6 spot gates (et1's bitwise bulk-chain
identity, et4-t, et5, et6, et8-t) confirm the new in-kernel branches did
not perturb the bulk expression trees (the S3-A8.1/A9.7 fma class did
not fire — these kernels launch only when the ET block is active, so the
default path never sees them).

## 5. Tier P / Tier N (S3-A1 scope)

T7 touched the atmosphere/ET layers (groundwater + surface consumers):
Tier N = the b5 24 h record + b6 ss/td + g5 ×2 + g8 alongside the
anchored Tier P sweep, g2, unit/mpi/stage3 and perf solo. Results and
walls: `stage3/records/T7/timings.md`.
