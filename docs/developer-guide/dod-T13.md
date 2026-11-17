# DoD — Stage-3 phase T13 (side-face right-hand-side scatters)

Date: 2026-11-17. Plan: `FREHG2_STAGE3_DEVELOPMENT_PLAN.md` §12E; amendment S3-A20. Companion: `report-T13.md`.

| Check | State | Evidence |
|---|---|---|
| Start condition: `refs/stage3/T12` exists (d329fe4); `records/T12/flip/pipeline-status.txt` reads `flip pipeline done`; the working tree equals the T12 snapshot tree (a5601d8) | done | `git for-each-ref refs/stage3` |
| Frozen reference binary `stage3/bin/frehg-T12` (+ `.deps.txt`) | done | `stage3/bin/` |
| W13.0 seam: the two kernels moved, unchanged, into `scatterVelocityBcRhs` / `scatterOutflowRhs`. The refactor alone is bitwise equal to `frehg-T12` on b4 (0/33 datasets) and keeps default_bitwise and et9 green | done | `records/T13/baseline/bitwise-et9-seam.log`, `evidence/fp-contract-b4.log` |
| W13.0 gates authored first: `SweRs.*` and `unit.rs0` (`unit;stage3`, 8 threads, fails on skip); `gate_rs.py`, `tolerances/stage3-rs.yaml`, the harness gate `rs1`, `regression.s3.rs1`; `scripts/test_rs_gates.py` and `unit.rs_gates.negative` | done | `tests/CMakeLists.txt`, `run_regression.py` |
| **Red** rs0 on the T12 arithmetic: 16/16 launches lossy at 8 threads; 800/800 over 50 repeats at 4 and at 8 threads; the cases skip at 1 and 2 threads | done | `records/T13/red/rs0-*.log`, `rs0-repeat-summary.txt` |
| **Red** rs1 on `frehg-T12`: 3/3 gate runs fail (audit 1.2e-8 to 6.0e-8; up to 22 datasets differ). Evidence lanes: g4a departs from serial in 10 of 12 runs at 5/6 threads, 0 of 3 at 3 threads, 0 at 4/8 threads | done | `records/T13/red/rs1-T12-rep*.log`, `evidence/straddle-g4a-T12.log`, `evidence/g4a-T12-audit-columns.txt` |
| W13.1 `Kokkos::atomic_add` on both writes | done | `FreeSurface.cpp` |
| **Green** rs0: 0/800 lossy launches at 4, 8 and 16 threads | done | `records/T13/green/rs0-repeat-summary.txt` |
| **Green** rs1: 22 datasets bitwise, audit 4.4e-16 at 5 and 6 threads; 6/6 runs; g4a 20/20 runs bitwise at 3/5/6/7 threads | done | `records/T13/green/rs-gates-verbose.log`, `rs1-repeat5.log`, `evidence/straddle-g4a-atomic.log` |
| Battery: 26 checks pass | done | `unit.rs_gates.negative` |
| default_bitwise and the et9 pins green on the atomic binary | done | `records/T13/green/rs0-bitwise-et9.log` |
| C-BIT (§2.6): serial b4 differs from T12 at contraction level (η ≤ 3.4e-15 m); under `-ffp-contract=off` the atomic and in-place forms are bitwise. Recorded in S3-A20 | done | `records/T13/evidence/fp-contract-b4.log`, `b4-serial-seam-vs-atomic.txt` |
| `scripts/check_forbidden.sh` clean; `-Werror` build | done | — |
| Tier P on `build-race` (`-L '^(unit\|mpi\|regression\|stage3\|stage3_local)$' -j8`): 234/235 in 37 min. The one failure, `regression.perf_baseline`, is a timing failure under `-j8` load (+96 % in the untouched groundwater module). Run alone, as the T12 flip records do, it passes: groundwater +19 % WARN at load 6, `swe/free_surface` +0.2 % | done | `records/T13/tierP.log`, `tierP-perf-alone.log`, `timings.md` |
| Docs: plan §12E, §1/§14/§15/§17-10, S3-A20, `report-T13.md`, `testing.md`, coverage row, mkdocs nav | done | — |
| Tier N on `build-race` (final source; `FreeSurface.cpp.o` byte-identical to the Tier P build): b6 ss/td, g8, g5 ×2 5/5 (2 757 s); `stage3_nightly` 16/16 (5 121 s); b5 24 h record PASS (11 069 s, concurrent; equal to the T12 record at report precision) | done | `records/T13/tierN/` |
| Snapshot `refs/stage3/T13`; `phase.patch` vs `refs/stage3/T12`; sync to official by file copy (the owner commits and pushes) | done | `records/T13/phase.patch`, `records/T13/sync.md` |
