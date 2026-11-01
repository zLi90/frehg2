# DoD — Stage-3 phase T10 (solve-consistent surface transport)

Date: 2026-10-31. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §12B (added by S3-A14, owner decision O16); amendment S3-A15. Companion: `report-T10.md`.

| Check | State | Evidence |
|---|---|---|
| W10.0 gates authored first: `benchmarks/sx-solve-consistent/` (sx1 uniform beach, sx2 drying strip, sx4 front slope, sx5 uniform heat + inputs), `tests/regression/gate_sx.py` pure checks, `scripts/test_sx_gates.py` negative battery, `tolerances/stage3-sx.yaml`, harness gates sx1/sx2/sx3/sx4/sx5/sx-restart/sx-rank, 13 ctest entries | done | `records/T10/` |
| **Red** on the frozen T9 binary: all 7 harness gates exit 1 (token schema-rejected, capability absent); measured `consistent` baselines — sx1 uniform breach 31.6 psu, anchor 1.4e-3·M; sx2 stranding −0.109 = 24 % of the exfiltrated salt; sx4 drift 1.26e-1; ws1n anchor +28.4 = 2.4e-2·M | done | `records/T10/red/` |
| W10.1 C1: `recordSolveFaceVolumes()` — Q per face post-solve (the boundary-audit expression; velocity-BC and outflow faces included), allocated/run only under the token | done | `FreeSurface.cpp`; unit `SweModule.SolveFaceVolumesTelescope` |
| W10.2 C2/C3/C4/C6 in `stepSurfaceImpl<bool QVOL>` (compile-time; the fma rule's third measurement — runtime ternaries broke the et9 pin, template fixed it); coupler `gain`/`bounce` books; stage-cell volume + limiter skip; pass-through far-upwind substitution + inflow-mix infiltration; paired overdraft park (donor cap measured and rejected, C-SX) | done | S3-A15; pins green |
| W10.3 token: Config parse/resolve/describe, schema enum + ET rule extension + the conservative→subsurface-consistent rule; docs: `theory/transport.md` (diagnosis, Q, C-books, J1–J4), `theory/surface-water.md` (step order + recording note), `parameters.md`, `agents/case-setup.md`, coverage row | done | checkers green |
| **sx0** (3 unit tests: telescope, uniform-stays-uniform, drying-face-carries-salt) + `unit.sx_gates.negative` + 4 validates | done | 263/263; battery green |
| **sx1** green: J1 1.9e-7 (≤ 1e-6), J2 anchor 9.2e-13·M (≤ 1e-10), C4 adjust 3.2e-5·M (≤ 1e-4, booked), film_removed > 0 | done | `records/T10/green/sx1.log` |
| **sx2** green (redesigned — flat strip vacuous, tilted pond, dt 2 s; S3-A15): ledger identity 6.4e-14·M, `surf_adjust` 0 exactly, pit store 0.049 psu·m³ | done | `green/sx2.log` |
| **sx3** green (ws1n, production solver): anchor 4.7e-9·M (≤ 1e-7, the production-residual class; red 2.4e-2·M), I2, both ws2 identities | done | `green/sx3.log` |
| **sx4** green: identity 3.6e-9·M, extrema in [5, 35] ± 1e-9, front moves | done | `green/sx4.log` |
| **sx5** green (needed the tide-temperature BC + the C6 exchange fix): heat anchor 1.9e-14·M, |T−20| ≤ 1.9e-7 | done | `green/sx5.log` |
| **sx-x** green: restart bitwise (ckpt 1800); rank n = 1/2/4 strict ≤ 5.8e-12 with the 1e-12-bed twin control (scalar bound 1e-7, measured 1.28e-9, S3-A15); seam (pd2 coupled + override) | done | `green/sx-restart.log`, `green/sx-rank.log`, trio log |
| Bitwise: `default_bitwise` + et9 noET + et9-veg green on the T10 binary (`legacy`/`consistent` untouched) | done | stage3 tier |
| W10.4 override re-runs: ws1/ws1n/ws4/et8/et8-t/et9-veg PASS; ws3 (a-ET) and et9 noET bitwise sub-checks fail BY DESIGN (consistent-pinned lanes; adjudicated) | done | `work/w104-*.log`, S3-A15 |
| W10.5 before→after table: ws1n anchor 2.4e-2·M → 4.7e-9·M; et8d stranding 24 % → 0; b6 ss/td PASS with tolerances unchanged (td salt mass +9.89 %, cap 10 %); v30 H2 smoke 35.000 psu everywhere, anchor 5.2e-9·M | done | report §5, `work/w105-*.log` |
| W10.6 flip evaluation: comparison build — 12/13 default-bitwise lanes bitwise, only the consistent-using coupled-salt lane and the et9 noET pin change ⇒ report, owner decides (item 7b) | done | `work/w106-*.log` |
| S3-A15 amendment (plan + mirror); dod/report; mkdocs nav | done | plan §17 |
| Tier P (unit\|mpi 119/119, stage3 61/61, regression 60/60, g2, perf) / Tier N (nightly 5/5, stage3_nightly 13/13 incl. the sx trio, b5 24 h record PASS) | done | `records/T10/timings.md` |
| Snapshot `refs/stage3/T10`; `phase.patch` vs `refs/stage3/T10-start` | done | `records/T10/` |
