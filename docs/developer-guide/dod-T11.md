# DoD — Stage-3 phase T11 (predictor–corrector consistency)

Date: 2026-11-04. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §12C (added by S3-A16, owner decision O17); amendment S3-A17. Companion: `report-T11.md`.

| Check | State | Evidence |
|---|---|---|
| W11.0 corrections BEFORE any gate: the corrector verified ρ-weighted (`updateWaterContent`), the D-M/D-V dichotomy withdrawn, M4 narrowed to the terrain-gravity term = the S1 force (one P5 change closes both); pc4 1-cell wall redesign; pc7 identity-form gate (the raw Σs·θV assert would re-assert booked limiter terms) | done | S3-A16/A17 |
| W11.1 gates authored first: `benchmarks/pc-predictor-corrector/` (13 cases + `input/`), `tests/regression/gate_pc.py` pure checks, `scripts/test_pc_gates.py` negative battery, `tolerances/stage3-pc.yaml`, harness gates pc1–pc8/pc-restart/pc-rank, 8 per-PR + 3 nightly + 8 validate ctest entries | done | `records/T11/` |
| **Red** on the frozen frehg-T10 (tokens stripped): rl2 churn/T 2.91; pc2-H2 creation 6.97 m³/m/cycle; pc2-H0 residual 0.021; pc3 arithmetic interface 3.4× the exact series flux; pc4 arithmetic wall leak 5.36e-6 m³/s; pc6a restore = 36 % of the ponded infiltration (M2); pc6b predictor seals the artesian face, restore mints 0.33 m³/h (M3); pc8 salty slab creeps at 1.01e-2·Ks·A (S1), fresh control 1.1e-13; pc7a flat lock θV drift exactly 0 | done | `records/T11/red/`, `measured-red-T10.md` |
| W11.2 shared face rule `src/gw/FaceConductivity.hpp`; corrector `computeFluxesImpl<PC,HARM>` takes each flank's own Ks on all three axes (M1); predictor/corrector templated compile-time per the fma rule (S3-A8.1), PC=false trees literal legacy | done | `Predictor.cpp`, `Corrector.cpp` |
| W11.3 M2: ponded top face Ksz in BOTH halves (O17c). M3: pressurized dry column (h_top > rZp·dz/2) gets the Dirichlet seepage leg in both halves, mode 3, clamp booked to `pc_clamp_m3` (run record, MPI-reduced) | done | `Predictor.cpp`, `Corrector.cpp`, `Simulation.cpp` |
| W11.4 P5 (M4+S1): the terrain-slope gravity term in the z-gravity form — K·vXp·rXp²·sinα predictor RHS, K·vXp·rXp·sinα corrector flux — a resting salty slab on a 10 % slope stays at rest | done | pc8 green |
| W11.5 tokens `groundwater.predictor_corrector: legacy\|consistent` (default legacy), `groundwater.face_conductivity: arithmetic\|harmonic` (+ schema rule harmonic⇒consistent); parse/resolve/describe; docs: `theory/groundwater.md` (T11 section + quirk rewrite), `theory/exchange-flux.md` (mode 3), `parameters.md` (240 tokens), `agents/case-setup.md` §9b, coverage (+2 rows, 140) | done | checkers green by exit code |
| **pc0** (4 unit tests: hetero column restore books nothing under both face rules, legacy churns — red control in-tree, homogeneous consistent ≡ legacy bitwise K3 `EXPECT_EQ`, ponded top face) + `unit.pc_gates.negative` | done | 267/267 |
| **pc1** green: rl2 lane under `consistent` — churn 3.061e-3 = 5.86e-3 of T (cap 0.01; red 2.91); gw ledger 1.1e-13; anti-vacuity 320 saturated Ks-contrast faces | done | `work/third-pc1.log` |
| **pc2** green (nightly, v30 H2 beach at T/4): gross created 0.0003 m³/m/cycle (cap 0.01; red 6.97; the report's Ks-only patch 0.045); ledger 5.4e-3 vs the 1e-2 coupled-subcycled class | done | `green/pc2.log` |
| **pc2-h0** green: ledger 4.0e-3 (stock 0.021 created → 0.0005, recorded) — the H0 residual was M3's | done | `green/pc2-h0.log` |
| **pc3** green: every interior face at its exact series resistance (harmonic, z/stretch/x); steady |dh| = 0.0; arithmetic lane records a 25.5× interface conductance (anti-vacuity) | done | `work/second-pc3.log` |
| **pc4** green: max \|q\| through the Ks = 0 wall = 0.0 exactly (harmonic); arithmetic leaks 5.36e-6 m³/s (anti-vacuity) | done | `work/first-pc4.log` |
| **pc5** green: wetting front crosses the 100× interface, θ below it 0.469 (> 0.2) | done | `work/third-pc5.log` |
| **pc6** green: (a) M2 — K1 at saturation: last-quarter churn exactly 0.0; ratio consistent/legacy 0.0433 (≤ 0.05); (b) M3 — artesian seepage 5.0000e-6 vs analytic 5e-6 m³/s (rel 4.1e-14), restore 4.3e-14 of the exchanged | done | `work/third-pc6.log` |
| **pc7** green (O17d): lock + heat — gw ledger ≤ 2.7e-15 m³, scalar identity ≤ 4.0e-14·M, K1 front churn ≤ 1.9e-12 of V | done | `work/second-pc7.log` |
| **pc8** green (S1): salty slab max \|q\| = 1.42e-13 of Ks·A (red 1.01e-2); fresh control 1.15e-13 | done | `work/first-pc8.log` |
| **pc-restart** bitwise (checkpoint incl. `realloc_dropped`); **pc-rank** n = 1/2/4 strict: head ≤ 5.3e-15 | done | `work/first-pc-restart.log`, `first-pc-rank.log` |
| "M5" edge-leg suspicion adjudicated NOT a defect: v2 corner rule explains the probe's 3×; predictor edge leg exactly 2c·dtg (matrix dump); rl2 head-row churn is the chord class, inside pc1's cap | done | `work/m5-adjudication.md`, S3-A17 |
| Bitwise: `default_bitwise` + et9 noET pins green on the T11 binary (default `legacy` untouched) | done | stage3 tier |
| W11.8 override re-runs (15 lanes): rl3/pd2/pd5/dg1/et8/et9-veg/ws4/sx2/b6-gw-smoke PASS; rl2 (the M1-encoding floor), rl4 (anchor-ratio premise), et9 (noET pin), ws1/ws1n (anti-vacuity floors + 1.07e-9 identity slip; `pc_clamp_m3` 3000× larger — booked consistently), sx1 (2.12e-6 vs 1e-6 psu bound) fail BY DESIGN/adjudicated | done | `work/w118-*.log`, S3-A17 |
| W11.9 flip evaluation (comparison build, default → `consistent`): default-bitwise 11/13 lanes bitwise un-pinned (only gw-orient + coupled-salt change); 13/13 WITH the new `BITWISE_PINS` entry; legacy-golden gates b3/b5/b6 + et9 on the flipped binary adjudicated face by face ⇒ report, owner decides (item 8c) | done | `work/w119*-*.log`, report §6 |
| S3-A17 amendment (plan + mirror); dod/report; mkdocs nav | done | plan §17 |
| Tier P (unit\|mpi 129/129, stage3 56/56, regression 60/60, g2) / Tier N (nightly 5/5, stage3_nightly 16/16 incl. pc2/pc2-h0/pc-rank, b5 24 h record PASS); perf_baseline deferred — both attempts contended (W11.9 b5, then an external SERGHEI run); rerun on a quiet machine | done | `records/T11/timings.md` |
| Snapshot `refs/stage3/T11`; `phase.patch` vs `refs/stage3/T11-start` (961b9884) | done | `records/T11/` |
