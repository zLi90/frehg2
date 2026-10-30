# Report — Stage-3 phase T4 (diagnostics and output)

Date: 2026-10-28. Companion: `dod-T4.md`. Plan §8; amendment S3-A7 (the registry's Field3 leg deferred to its first consumer).

## 1. What landed

| Piece | Where |
|---|---|
| `io::IntervalRates` — a registry of named per-cell depth accumulators: contributors add each step's deposit, each output write emits `acc/(t − t_prev)` (an interval-MEAN rate, m/s) and resets; only requested variables allocate; checkpoint `acc_<name>` + scalar `acc_t_prev` | `src/io/IntervalRates.{hpp,cpp}` (W4.1) |
| `seepage_rate` (surface, coupled only): interval mean of the coupler's complete per-step `gain` — applied seepage + bounce-back + vent deposits — positive upward; the write normalizes over the TRUE elapsed state time so adaptive loops that cross an output label stay exact | `Coupler::gain()`, `Simulation.cpp` (W4.2) |
| `output.monitor_interval` (whole seconds, 0/absent = every step): audit ACCUMULATION every step unchanged, row WRITES decimated to the first step reaching each multiple, forced at checkpoint instants and t_end; point monitors sample at the same instants; restart recomputes the schedule from t0 | `Config.*`, `ConfigSchema.cpp`, `Simulation.*` all three loops (W4.3) |
| Docs: snapshot-vs-rate semantics, decimation, the interval-rate pattern for future variables, agent guidance ("integrate `seepage_rate`, never snapshot `seepage`") | `output.md`, `parameters.md`, `adding-an-output-variable.md`, `agents/postprocessing.md` (W4.4) |
| Gates dg1/dg2/dg-restart, `gate_diagnostics.py`, `stage3-dg.yaml`, `unit.dg_gates.negative`, `ConfigTest.MonitorIntervalAndSeepageRate…` | `benchmarks/dg-diagnostics/`, `tests/regression/`, `scripts/test_dg_gates.py` |

## 2. The defining identity (dg1)

An interval-mean rate is only trustworthy if it integrates back to the
audited volume. dg1 holds `Σ_cells rate·A·Δt_out = Δ(mass_audit.seepage)`
per output interval on a coupled box with BOTH signs of exchange active
(a side head 0.3 m above the surface exfiltrating, a ponded corner
infiltrating over a deep water table). Both sides sum the very same
per-step `gain` numbers — the rate field per cell over steps, the audit
per step over cells — so only summation-order rounding separates them:
**measured 6.7e-16 against the authored 1e-12** (first green, no
adjudication). The anti-vacuity lane reproduces the V3 finding that
motivated the feature: with dt = 1 s and dtg ≤ 10 s most surface steps see
no subsurface substep, so at ≥ 1 output instant the per-step snapshot
`seepage` is identically 0 on every wet cell while `seepage_rate` ≠ 0.

## 3. Decimation is a row-subset property (dg2)

Accumulation and row-writing are now separate concerns: a decimated run's
audit row at time t is BITWISE the every-step run's row at t (61/61 rows
at `monitor_interval: 60`, three tables). A third lane beyond the plan's
two (noted in S3-A7) sets `monitor_interval: 77` with a checkpoint at
1800 s — 77 divides neither 1800 nor 3600, so the checkpoint-forced and
t_end-forced row paths are proven by execution (49 rows, the schedule set
asserted exactly), not left to coincide with the 60 s multiples.
`default_bitwise` is the other half of the guarantee: the key absent means
every step, bitwise vs the T0 hashes.

## 4. Restart (dg-restart)

The accumulators are prognostic: `acc_<name>` fields plus the `acc_t_prev`
scalar checkpoint and restore, and a missing scalar under requested rate
outputs is fatal with a named message. The gate is **bitwise** as authored
— restart at 1800 s reproduces every `seepage_rate`/`seepage`/eta/depth
frame exactly. (The S3-A6 surface-restart relaxation was not needed: this
case has no superbee/wind branch-amplification noise.) Runs without rate
outputs keep the v2.0 checkpoint layout byte-identical.

## 5. S3-A7 (scope note, not a tolerance change)

The registry lands Field2-only; the plan's Field3 (1/s) leg has no
consumer before T6's `root_uptake` (W6.6) and authored-unexercised code is
this project's documented failure class (V2-A9..A12). It is added in T6
with the gate that exercises it. No tolerance was adjudicated anywhere in
T4; no zero_gradient or legacy behaviour was touched (dynamics identical
by `default_bitwise`).

## 6. Tier P / Tier N (S3-A1 scope)

Tier P: anchored `^regression$` sweep at `-j8` minus `perf_baseline` +
`scaling.g2.iters.amg` + `perf_baseline` solo + unit/mpi/stage3 at `-j4`
(adds dg1/dg2/dg-restart, `unit.dg_gates.negative`, the ConfigTest and
`validate.dg1-seepage`). Tier N (T4 touched the io/driver/coupling output
layers; dynamics untouched by `default_bitwise`): the b5 24 h record +
b6 ss/td. Results and walls: `stage3/records/T4/timings.md`.
