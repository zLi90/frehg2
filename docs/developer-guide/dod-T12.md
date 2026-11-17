# DoD — Stage-3 phase T12 (saturated-contact classification)

Date: 2026-11-16. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §12D (owner decision O19); amendment S3-A19. Companion: `report-T12.md`.

| Check | State | Evidence |
|---|---|---|
| Trigger reproduced: the frozen `frehg-T11` b5 one-step strict lane is rank-dependent on the dev machine at n = 3/5/6 (2.76e-6 / 2.37e-5 / 2.76e-6). Every diverging cell sits on a subdomain seam column (i = 34; i = 81–83) | done | `records/T12/red/sc1-legacy-measured.log` |
| W12.0 gates authored first: `tolerances/stage3-sc.yaml`, `tests/regression/gate_sc.py` pure checks, `scripts/test_sc_gates.py` battery, harness gates sc1/sc2, ctest entries (`regression.s3.sc1` labelled `stage3;mpi`, `.sc2` `stage3`, `unit.sc_gates.negative`) | done | `tests/CMakeLists.txt`, `run_regression.py` |
| **Red** on `frehg-T11`: sc1/sc2 capability absent (the key is rejected); sc0's head assertions on the legacy rule fail 36 times (both GwSc.Head* tests); the battery passes on its synthetic inputs | done | `records/T12/red/` |
| W12.0 static audit of all 162 kernels: no other in-place neighbour read. Two scatter races found (`swe_velocity_bc_rhs`, `swe_outflow_rhs`) and reported, not fixed | done | `records/T12/kernel-audit.md` |
| W12.1 key `groundwater.saturated_contact: legacy \| head`: parse, schema, resolved-config round trip (r1), `parameters.md`, coverage row; checkers green by exit code | done | `unit.parameter_docs`, `unit.feature_coverage`, `unit.run_record_checker` |
| W12.2 snapshot `wcSnap_`/`hSnap_` (head only) plus the read aliases `wcR`/`hR` in `gw_reallocate_classify`. The `sat` predicate is θ ≥ θs ∨ (head ∧ h > h_ae). No template and no `if constexpr` (nvcc rule moot). Header invariance claim corrected | done | `Reallocate.cpp`, `RichardsSolver.*` |
| **sc0** green (3 GwSc unit tests: mirror symmetry bitwise on x/y/z; whole fringe restored; legacy one-way chain recorded as anti-vacuity) | done | `unit.all` |
| **sc1** green: head ≤ 1.04e-14 at n = 2, 3, 4, 5, 6, 8; restore 1.911e-2 m³; legacy anti-vacuity 2.374e-5 at n = 5 | done | `records/T12/green/sc1-sc2-verbose.log` |
| **sc2** green: 1 vs 4 OpenMP threads, 5 fields bitwise | done | same |
| `default_bitwise` green under `legacy` before the flip | done | `records/T12/green/stage3_local.log` |
| W12.3 `head` override on 55 Stage-3 lanes: 54 pass; et9's noET pin breaks by design. The et9 case files pin `legacy` | done | `records/T12/w12.3/`, S3-A19 |
| W12.4 flip (`Config.hpp` initializer → `Head`; `BITWISE_PINS` += `saturated_contact: legacy`): unit\|mpi 131/131, stage3 58/58, regression 60/60, Tier N 5/5, stage3_nightly 16/16, b5 24 h PASS, perf PASS, tolerances unchanged → flip kept | done | `records/T12/flip/` |
| Two S3-A18 recalibrations reversed with data: gw-orient `tol_abs_mixed` 5e-7 → 1e-12 (head 1.8e-15; legacy 7.6e-8); b5 default-lane bounds 1.5e-4/4e-3 → the A14 3e-5/2e-3 (head n=4 3.2e-6 / 1.7e-4) | done | `records/T12/w12.3/gw-orient-*.log`, `records/T12/flip/b5-rank-restored-bounds.log` |
| W12.5 docs: `theory/groundwater.md` (new section, deviations row), `symmetry-exemptions.md` gw row, `testing.md`, `agents/case-setup.md` 9c, `parameters.md` (default `head`), V2-A9 correction note, plan §12D/§17-9/S3-A19, mkdocs nav | done | doxygen and mkdocs `--strict` exit 0 |
| Snapshot `refs/stage3/T12`; `phase.patch` vs `refs/stage3/T11` | done | `records/T12/` |
| Sync to official by file copy (the owner commits and pushes); x86 CI is the final check | done | `records/T12/sync.md` |
