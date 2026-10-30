# DoD — Stage-3 phase T5 (Feature C core: coupled evaporation with salt conservation)

Date: 2026-10-29. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §9; amendment S3-A8 (nine first-green adjudications and case redesigns, each with measured data).

| Check | State | Evidence |
|---|---|---|
| et0 green: `AtmEt.*` 5/5 — potentialBulk ≡ the v2.0 evaporationRate tree BITWISE (liu2006, r_s = 0, both over u and α); FAO-56/log-profile hand values to 1e-12; Lee–Pielke β closed form; the supply cap both signs; θ_fc = θ(h = −3.3 m) | done | ctest `frehg_unit_tests` |
| et1 green: lane a holds the new/v2.0 chain agreement at the fma-contraction class (bitwise through substep 15; measured peak θ 1.3e-9 ≤ 1e-7 — S3-A8.1); lane b (the coupled D1 red, S3-A8.2) measures a 5872 m guard-mismatch head difference with the new path's ledger at 1.7e-13 | done | `records/T5/green/et1-green.log` |
| et2 green: coupled-dry ≡ uncoupled at 8.0e-11 (≤ 1e-10 strict); depth ≡ 0 anti-vacuity | done | `records/T5/green/et2-et4-green.log` |
| et3 green (bulk + prescribed + the v2.0 real negative): pond phase 2.5e-15/3.1e-15; wet+soil ≤ t EXACT; combined rate 0 / 2.9e-15 (through the α₁ = 1 horizon, S3-A8.6); post-switch continuity 1.000000; the v2.0 window-lag lane drains the pond early (removal ratio 0.26 — a 74 % violation vs the 1e-3 floor) | done | `records/T5/green/et3-green.log` |
| et4 green (12 h lane, S3-A8.3): boundary_in = −soil_evaporation at 0; ΔS identity 4.2e-15; vloss ≡ 0; the 3-day vG dry singularity documented (a pre-existing Richards domain limit, the S3-A4 class) | done | same log |
| et7 green: θ ≥ θ_r − 4e-6 above the floor; unrealised = demand − applied at 0; vloss = 0 | done | same log |
| et8 green (a/b/c/d): a total mass EXACT, per-cell c·θ 6.7e-7 (the tiny-K class, S3-A8.4), adjust/anchor ≤ 7e-17/6e-16; b s·h = 35 h₀ wet, store peak = 35 h₀ A, total 4.3e-15; c (real negative) discard loses exactly 35 h₀ A; d full identity 4.2e-15 with **surf_adjust ≡ 0 exactly** (D9), store 0.449 while held (S3-A8.5 driver redesign) | done | `records/T5/green/et8-green.log` |
| et9 green: zone assertions exact (subtidal soil ≡ 0, supratidal wet ≡ 0); complementarity 0.0; both channels active (4.2/2.5 m³); salt identity 2.7e-9 (default-solver class, S3-A8.7); **the noET lane BITWISE vs its T5-start pin (1271 datasets)** — T5 default-neutrality by execution | done | `records/T5/green/et9-et-restart-green.log` |
| et-x green: rank invariance at MACHINE precision (θ 1.3e-15, tracer 1.1e-13, store 0.0 over n = 1/2/4); battery transpose EXACT (tracer 0.0, store 0.0, θ 0.0); mixed lanes recorded (the legacy live-edge-face asymmetry, S3-A8.9, measured on a provably static base) | done | `records/T5/green/et-rank-green.log`, `et-orient-green.log` |
| et restart green: et3 mid-transition + et8b with a full store — eta/depth/clocks/tracer/store all BITWISE | done | `records/T5/green/et9-et-restart-green.log` |
| D1–D9 removed on the ET path (D3 subsumed, D7 out of scope); legacy paths untouched: g4 + g5 + `default_bitwise` green on the T5 binary | done | `records/T5/legacy-check.log` |
| Negative battery `unit.et_gates.negative` green (the bulk chain two independent ways; t_d; the salt closed forms; every checker fails its corruption) | done | ctest |
| Red evidence (the T4 binary rejects every case on the T5 target keys; harness red runs) | done | `records/T5/red/` |
| W5.11 docs: `theory/evapotranspiration.md` (new, in nav), parameters (216 tokens), configuration, output, transport, atmosphere, case-setup §6c, coverage rows (117); doxygen + mkdocs --strict + all lockstep checkers clean | done | checker outputs |
| Tier P green; S3-A1-scoped Tier N green (b5 24 h record + b6 ss/td + g5 ×2 + g8 — T5 touched surface + groundwater + transport + coupling) | done | `records/T5/timings.md` |
| Snapshot `refs/stage3/T5`; `phase.patch`/`phase.stat` | done | `records/T5/` |
