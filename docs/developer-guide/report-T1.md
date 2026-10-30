# Report — Stage-3 phase T1 (`groundwater.reallocation: legacy | disabled`)

Date: 2026-10-27. Companion: `dod-T1.md`. Plan §5; amendments S3-A2 (compiled gates land with the config change), S3-A3 (`disabled`, not `off` — the YAML 1.1 boolean trap), S3-A4 (semantics: four designs measured, design 3 adopted with a vadose domain of applicability; gate rescope).

## 1. What landed

| Piece | Where |
|---|---|
| `groundwater.reallocation: legacy \| disabled` (enum `Reallocation`, the T8 DF constant) | `Config.hpp/.cpp`, `ConfigSchema.cpp`, resolved-config echo, run-record closure line |
| `disabled` semantics (S3-A4 design 3): the legacy sweep minus every non-conservative action — no sends, no drops, no θ := θ(h) creation; the conservative isolated-cell h := h(θ) restore kept verbatim; the over-saturation excess goes to the audited finalize clamp | `gw/Reallocate.cpp` (mode branches inside the legacy kernel; file header carries the domain-of-applicability statement and the S3-A4 history) |
| Inert-key warnings under `disabled` (surplus / reallocation_gradient / realloc_* outputs) through the W0.9 channel | `ConfigSchema.cpp` |
| Per-cell cumulative diagnostics `realloc_created` / `realloc_dropped` / `vloss` [m³]: allocated only when requested, booked at every creation/discard site of both sweep phases, output + checkpoint | `Reallocate.cpp`, `RichardsSolver.*`, `Corrector.cpp` (vloss field), `Simulation.cpp` |
| Gates rl1–rl5 + rl-restart, `gate_realloc.py` pure checks, `stage3-rl.yaml` tolerances (every bound carries its authoring measurement), `unit.rl_gates.negative` | `benchmarks/rl-realloc/`, `tests/regression/`, `scripts/test_rl_gates.py`, `tests/CMakeLists.txt` |
| Docs: theory "Post-allocation modes" section, parameters/output rows, coverage rows | `docs/` |

## 2. The S3-A4 finding (the phase's main technical result)

The plan's W1.2 semantics ("every unsaturated cell re-derives h from the retention
curve") is **not viable**: four designs were implemented and measured (table in
the amendment). The adopted `disabled` is exactly conservative and deterministic
everywhere, clean in vadose regimes, and — by measurement — the PCA scheme has no
h↔θ consistency mechanism near saturation other than the legacy, water-creating
θ := θ(h) restore. In saturated regimes `disabled` sheds/mints clamp volume
(audited per cell as `vloss`) and amplifies rounding chaotically. A conserving
redesign of the saturated-contact restore is plan contingency **C-RL — an owner
decision**, with these measurements as its data:

- saturated gw-orient block: per-cell |vloss| up to 0.33 m³ in 0.0625 m³ cells; O(1) orientation divergence (deterministic per run);
- rl2 saturated box: vloss 5.743 m³ over T = 2.380 m³; substeps 562 vs 503 (1.12×);
- vadose rl3: vloss exactly 0; substeps ratio 1.000; steadiness exactly 0; flux uniformity 9.3e-13;
- vadose rl5 battery: strict 1e-12 PASS under `disabled`; coupled-salt strict rank lane PASS;
- rl4 coupled-salt anchor: 38.70 → 1.366 over M = 5.80e4 (28× reduction).

## 3. Gate-first record

rl2/rl3/rl4 red ("capability absent") on the T0 binary: `records/T1/red/`.
First green-attempt failures that reshaped the gates (all in S3-A4):
the `off` YAML-boolean round-trip corruption (S3-A3), the over-forced rl2 pulse
(2e-5 m/s > silt Ksat — both modes shed ~90 % of it; now 8e-7), the rl3 ledger
absolute floor, the rl4 three-part anchor bound, and the churn anti-vacuity floor
re-anchored to the legacy lane's own throughput.

## 4. W1.6 legacy pins

Tests asserting legacy sweep semantics, pinned to `Reallocation::Legacy` with the
plan-§2.7 comment: `GwModule.ReallocationSurplusModesDropOrConserve`,
`GwModule.ReallocateRestoresIsolatedHeadConsistency`. (All other reallocation
mentions in the unit tier assert mode-independent identities or coupled audit
plumbing and run on the default.)

## 5. Tier P / Tier N (S3-A1 scope)

Tier P: regression 60/60 at `-j8` (2566 s) + `perf_baseline` solo PASS
(128 s) + unit/mpi/stage3 77/77 (376 s). Tier N (S3-A1 scope, `-j2`):
b6.ss 316 s, b6.td 740 s, g5 1858 s, g5-transposed 1810 s, g8 52 s — 5/5;
b5 24 h record PASS with the T0 figures (discharge 0/48 outside, ponding
0/54). Details: `stage3/records/T1/timings.md`.

## 6. W1.8 DF preview and W1.9 H2 evidence

Full analysis: `stage3/records/T1/df-preview.md`. Headlines: under forced
`disabled`, 13/60 regression gates and b6 ss/td fail — all in the S3-A4
saturated-regime classes (b2/b3 goldens embed the legacy restore; coupled
water-table cases slow and drift; saturated batteries amplify rounding).
The W1.9 evidence run on the study's own `c_H2_2D_U00` (2 tidal cycles):
`legacy` completes and reproduces the report_V2 creation signature
(realloc net +2.80 m³ vs physical net boundary exchange −2.66 m³, gross
drop 17.47 m³, per-cell sums exact); **`disabled` diverges at
t = 1189 s** (fs solve indefinite after a fringe-driven CFL blow-up —
a loud abort, no corrupt output). Consequence for O9: unless C-RL lands a
conserving saturated-contact restore before T8, the default-flip
evaluation will keep `legacy`, and `disabled` ships as the documented,
gated vadose-domain option.

## 7. Snapshot note

`refs/stage3/T1` also carries T2's pre-authored gate layer (authored during
the T1 pipeline's wall time, per the gate-first rule): `benchmarks/
ob-openbc/` (ob1/ob3/ob5 cases + inputs + README), `gate_openbc.py` (the
discrete-chain reference, battery-validated: lane A 0.3052 m/s vs the
continuum closed form 0.3067, entrance ratio 0.9996; zero-gradient contrast
4.39x), `scripts/test_ob_gates.py`, and the filled `stage3-ob.yaml`. No
Stage-2 source or harness file is touched by them; they are red by
construction until T2's W2.1 lands.

## 8. Deviations

S3-A2, S3-A3, S3-A4 (see the plan's amendment log). No tolerance of any
pre-existing gate was changed.
