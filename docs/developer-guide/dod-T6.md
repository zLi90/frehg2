# DoD — Stage-3 phase T6 (vegetation and transpiration)

Date: 2026-10-29. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §10; amendment S3-A9 (seven entries: six first-green adjudications and the pre-existing coupled-transport restart defect found and fixed).

| Check | State | Evidence |
|---|---|---|
| et0 (part 2) green: `AtmEt.BeerLambertLimits` / `FeddesPiecewise` (every segment + the wetland no-h₁h₂ form) / `CanopyResistanceAndRootIntegral`; `GwEt.RootFractionsIntegrateToOneOnStretchedMaskedColumns` (1×2, stretch 1.2, masked column, L cuts a cell; uniform + linear); `GwModule.SinkInSaturatedConfinedCellDropsHeadBySinkOverSs` (dh = −S·dtg/Ss to 1e-12, θ = θ_s to 1e-14 — the predictor/corrector sign fixed by test) | done | ctest 10/10 |
| et4-t green: ΔS = boundary_in − ss + realloc − vloss − **transpiration** at 2.7e-14; both channels active (T 5.7e-3 / E 1.5e-2 m³); the drop-T negative misses by 0.27 | done | `records/T6/green/` |
| et5 green: (i) the ALL-ZERO vegetation map **bitwise** vs vegetation-absent (3 groups); (ii) dense canopy (τ = 9e-14): E_soil 1.1e-14 vs T 9.8e-2, top-cell share 0.051 ≤ 0.10; (iii, real negative) the bare 2 h contrast concentrates at the top (share 0.53 ≥ 0.5) — S3-A9.3 | done | same |
| et6 green: per-cell Δθ_k = uptake_k·t (the `root_uptake` Field3 rate — S3-A7's deferred leg, first consumer) at 3.1e-10 ≤ 1e-9 (the K = 1e-15 class), uniform + linear, on the stretched masked geometry; the uniform lane's shared T_pot/L sink value exact — S3-A9.1 | done | same |
| et7-t green: below-h₄ shutoff (T ≡ 0, unrealised = demand, E_soil > 0); θ ≥ θ_r − 1e-14 under an unfundable demand with T + unrealised = demand at 0; vloss ≡ 0; T_act ≤ T_pot per row | done | same |
| et8-t green: total salt exact (2e-16); per-cell c·θ 3.0e-8 (tiny-K class); adjust/anchor ≤ 1.3e-16; salinity-rise ~ rootFrac correlation 0.9999931 ≥ 1 − 1e-4 (S3-A9.4) | done | same |
| et9-veg green: complementarity exact; T > 0 on every marsh cell (min 2.0e-3) **including inundated instants (O7 observed)**; T ≡ 0 outside the band; cumulative T 0.156 m³ | done | same |
| et-veg-restart green: eta/θ/root_uptake/salt_store/clocks/T-rate all **bitwise** across a mid-cycle checkpoint — after the S3-A9 restart fix | done | same |
| **S3-A9 restart defect fixed**: a zero-substep post-restart window read stale (allocation-zero) per-substep fluxes/conductivities/dtg and wiped the whole subsurface scalar (subs_adjust = −M in one step; reproduced with no T6 feature, legacy AND consistent, since the first coupled-transport restart; b6-restart's dtg ≤ dt masked it). Checkpoints now carry `q*_last`, `k*_last`, `dtg_last`; the minimal control restarts **bitwise**; ALL restart gates re-ran green (10 ctest lanes + et/dg/pd harness lanes) | done | `records/T6/green/all-restart-gates.log` |
| W6.4 sink in separate post-kernels (S3-A9.7): the pinned assembly/update kernels stay byte-identical — `default_bitwise` green on the final binary (an in-kernel guarded branch had broken it on b2, the twice-measured fma class) | done | solo run + stage3 label |
| T5 gates still green; vegetation-absent bitwise (et5(i) + the et9 noET pin, 1271 datasets, on the final binary); g4 green | done | `all-t6-gates-final.log` |
| W6.7 docs: theory partition/Feddes/roots/salt/O7 sections, parameters (229 tokens), input-data (map ids), case-setup §6d marsh recipe, coverage (123 rows); doxygen + mkdocs --strict + forbidden + bc_matrix clean | done | checker outputs |
| Tier P green; S3-A1-scoped Tier N green | done | `records/T6/timings.md` |
| Snapshot `refs/stage3/T6`; `phase.patch`/`phase.stat` | done | `records/T6/` |
