# DoD — Stage-3 phase T3 (periodic lateral topology)

Date: 2026-10-28. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §7; amendment S3-A6 (per-field-class bounds and the coupled-lane redesign, each with its measured datum).

| Check | State | Evidence |
|---|---|---|
| pd1 green: discrete wind–friction balance (strict lane spreads ≤ 1e-13, u vs the bisection balance 0.07446 m/s; default-solver translation asymmetry recorded; transposed periodic-y lane green — O3 axis symmetry by execution) | done | `records/T3/green/` (stage3 tier log) |
| pd2 green, 3 lanes: surface roll identity at the S3-A6 field classes (uu 7.3e-8 ≤ 1e-6, tracer 7.3e-6 ≤ 1e-4, seam-crossing anti-vacuity peak 7.55); gw lane strict (1e-13 ≤ 1e-11); coupled lane (saturated tidal basin — uu 1.5e-7, head 2.1e-11, c_gw 2.0e-10; exchange anti-vacuity 0.284 psu ≥ 0.05) | done | same |
| pd3 green: wrapped-rank invariance n = 1/2/4 (+ coupled 2 × 2) — gw strict ≤ 5e-13 vs the 1e-12 bound; surface/coupled at the pd2 field classes (uu 7.0e-8/1.4e-7), with the walled control measuring ~10× larger deviations than the periodic lane on every field | done | same |
| pd4 green: the seam books no boundary flux (ledger closure; unit test `boundaryOutflow == 0.0` exactly on a periodic channel) | done | `SweModule.PeriodicSeamCarriesNoBoundaryFlux`, stage3 tier log |
| pd5 green: confined Darcy — head x-spread/curvature 3e-15, layer-shape parallelism (hydrostatic per-layer profiles), qx 4.5e-19 vs qy 2.0e-7 (denormal class), column-K lane invariant | done | stage3 tier log |
| pd6 green: nx < 3 rejected; x-edge-only BC fatal with the W3.4 message; the W3.5 seam warning fires on the non-periodic bump lane; `disabled`/periodic round trips pinned | done | same; `ConfigTest.Periodic*` |
| pd7 `default_bitwise` green (no `periodic` key ⇒ bitwise vs the T0 hashes, §2.6 declared exceptions only) | done | stage3 tier log |
| pd8 green, 2 batteries with the W3.7 flag swap: gw (rl5 vadose base, 1e-14 class vs the 1e-12 strict bound) and surface (`pd8-orient-base.yaml`, scalar fields — transpose at the pd2 classes, mixed at the x-orient swe class; measured ≤ 1.2e-4) | done | same |
| pd-restart green at the project restart standard (S3-A6): flow fields ABSOLUTE 1.9e-15 ≤ 1e-12; tracer 4.0e-6 rel ≤ 1e-4 — the walled control measures strictly larger on every field | done | same |
| The 58-site EdgeTopology conversion table fully converted; W3.8 forbidden rule (raw global-edge tests) clean on the converted tree | done | `scripts/check_forbidden.sh` clean; site table in `report-T3.md` |
| Conversion default-neutrality: 10/10 rl/ob/bitwise gates green across the T2→T3 binary boundary; `default_bitwise` PASS on the converted tree | done | `records/T2/timings.md` pipeline note; stage3 tier log |
| mpi.halo.n{1,2,4} gain `PeriodicWrapFillsHalosFromTheFarEdge`; pd3 runs under `FREHG_EXPECT_RANKS` | done | mpi tier log |
| Negative battery `unit.pd_gates.negative` green (discrete balance vs an independent fixed point; roll/ledger/linear-head/zero-flux checkers; S3-A6 lanes: measured classes pass, defect classes fail; coupled-exchange hollow-lane negative) | done | ctest |
| Red evidence (schema rejects `domain.periodic` on the T2 binary, all five original cases; the S3-A6-authored lanes fail on the same key) | done | `records/T3/red/` |
| W3.9 docs (grid-indexing, halo-exchange, configuration, parameters, theory seam sections, symmetry-exemptions, case-setup torus recipe, coverage rows incl. periodic × coupled) + lockstep checkers + forbidden scan + mkdocs/doxygen | done | `unit.parameter_docs`, `unit.bc_matrix` (120 cells), `unit.feature_coverage` (106 rows) green |
| Tier P green; S3-A1-scoped Tier N green (T3 touched surface + groundwater + transport) | done | `records/T3/timings.md` |
| Snapshot `refs/stage3/T3`; `phase.patch`/`phase.stat` | done | `records/T3/` |
