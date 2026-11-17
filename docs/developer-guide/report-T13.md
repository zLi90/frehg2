# Report — Stage-3 phase T13 (side-face right-hand-side scatters)

Date: 2026-11-17. Plan §12E, amendment S3-A20. DoD: `dod-T13.md`.

## 1. Why T13 existed

T12's static audit of every Kokkos kernel (`records/T12/kernel-audit.md`)
found two scatter races in the free-surface right-hand side:
- `swe_velocity_bc_rhs`, which does `Srhs(j, i) += add`;
- `swe_outflow_rhs`, which does `Srhs(j, i) -= volume`.

Both loop over a side-face BC list. `appendSideFaces` pushes one entry per
domain-edge face, so a corner cell holds two consecutive entries and an
nx = 1 or ny = 1 axis cell two or three. Two entries of one launch then
update one slot unprotected. The kernel's reduction books both
contributions, so a lost update shows up in the mass audit as volume that
was booked but never moved.

## 2. When the race fires

For lists shorter than 40 × threads, Kokkos' OpenMP `RangePolicy` gives
each thread ⌈n/T⌉ consecutive entries. A duplicated pair can race only when
a chunk starts at its second entry:

| Case (list) | Entries | Thread counts where a pair straddles |
|---|---|---|
| g4a (walls on all four sides) | 14 | 3, 5, 6 (not 4, 7, 8) |
| b4 (east outlet strip) | 12 | ≥ 12 |
| x-orient swe-orient (west wall, south part) | 3 | ≥ 3, but both entries there carry 0 |

Every existing threaded lane runs at 4 or 8 threads, and the corner cases
are gated at one thread. That is why nothing saw the race. On the frozen
`frehg-T12` binary, g4a departs from serial:

| Threads | Runs departing | Size |
|---|---|---|
| 5 | 4 of 6 | depth/η up to 3.8e-9 m |
| 6 | 6 of 6 | surface closure residual 1.8e-8 to 6.0e-8 m³ (serial: 7.6e-13) |
| 3 | 0 of 3 | — |

`bc_inflow` matches serial to 7e-12: the volume is booked, but it never
reaches `Srhs`. On CUDA, adjacent entries run in one warp, so a loss is
close to certain. CUDA is not measured here.

## 3. What landed

- **The fix.** Both writes are now `Kokkos::atomic_add(&Srhs(j, i), ±v)`.
- **The refactor.** The two kernels moved into free functions,
  `scatterVelocityBcRhs` and `scatterOutflowRhs` (`SurfaceSolver.hpp`,
  `FreeSurface.cpp`), which the two `SurfaceSolver` methods call. This lets
  a unit test drive them with a crafted list. The refactor alone is
  bitwise equal to `frehg-T12` on b4.
- **No configuration change.** No key, schema or default changed.

## 4. Gates (`stage3-rs.yaml`)

| Gate | Red on the T12 arithmetic | Green (atomic) |
|---|---|---|
| rs0 (`unit.rs0`, 8 threads): both free functions over a 65 536-entry list on one cell, with dyadic data so every partial sum is exact | 16/16 launches lose updates (velocity kept 4 907–7 287 of 20 480); 800/800 over 50 repeats at 4 and 8 threads | 0/800 lossy at 4, 8 and 16 threads |
| rs1 (`regression.s3.rs1`): g4a at 1, 5 and 6 threads; surface fields and monitors bitwise, audit ≤ 1e-12 m³ | 3/3 gate runs fail; audit 1.2e-8 to 6.0e-8 | 6/6 runs pass: 22 datasets bitwise, audit 4.4e-16 |
| battery (`unit.rs_gates.negative`) | — | 26 checks |

Anti-vacuity:
- `SweRs.SideFaceListsRepeatCornerAndAxisCells` pins that the production
  builder repeats corner and axis cells, as consecutive runs.
- `unit.rs0` fails if its threaded cases skip. They skip below 4 threads,
  which is why `unit.all` at one thread does not run them.
- rs1 rebuilds the list from the case YAML and asserts that a duplicated
  cell straddles a chunk boundary at every lane's thread count.

## 5. Serial bitwise, and the one place it is not (C-BIT)

`default_bitwise` (`stage3_local`) and the et9 pins stay green. Serial b4
does not stay bitwise on the dev toolchain:
- **What moved.** At t = 600 s, 22 of 33 datasets differ from `frehg-T12`:
  η by up to 3.4e-15 m, uu by up to 4.7e-14 m/s.
- **The cause.** GCC on aarch64 contracts `Srhs -= s * drop` into one
  fused multiply-subtract. The atomic's argument rounds `s·drop` first.
- **The check.** Built with `-ffp-contract=off`, the atomic and in-place
  forms agree bitwise (0/33).
- **Scope.** Only outflow conditions with a nonzero bed drop, on
  FMA-contracting builds. The velocity kernel has no contractible product,
  and g4a, swe-orient and default_bitwise are unchanged. No gate pins an
  outflow case bitwise, and Tier P is green.

## 6. For the owner

- **Accept or redirect the C-BIT change.** Accept the contraction-level
  movement of serial outflow runs on aarch64, or ask for the deterministic
  alternative recorded in S3-A20: one thread applies each run of a cell's
  consecutive entries in list order. That alternative keeps the T12
  rounding on FMA toolchains, and it also makes threaded runs bitwise to
  serial when two entries of one cell both carry volume. The atomic leaves
  those runs to the thread order, at rounding level. No shipped case has
  such a pair.
- **Tier P.** 234/235 at `-j8` in 37 min. `perf_baseline` failed under
  that load and passed alone: `swe/free_surface` +0.2 %, so the atomics
  cost nothing measurable.
- **Tier N (phase close, 2026-11-17).** b6 ss/td, g8 and g5 ×2 5/5;
  `stage3_nightly` 16/16; b5 24 h record PASS, its discharge and ponding
  numbers equal to T12's record at report precision.

## 7. State

- **Dev tree.** Uncommitted (O10). The change is
  `src/swe/{SurfaceSolver.hpp,FreeSurface.cpp}`, plus the gates,
  tolerance file, battery, ctest entries and docs.
- **Builds.** Built and tested in `build-race`. `build-race-nofma` exists
  only for the C-BIT evidence.
- **Records.** `stage3/records/T13/` holds `red/`, `green/`, `evidence/`
  (the thread lanes, the fp-contract test, the audit attribution),
  `tierP.log` and `phase.patch`.
- **Snapshot.** `refs/stage3/T13`.
- **Official repo.** Synced by file copy 2026-11-17 (`records/T13/sync.md`),
  uncommitted; the owner commits and pushes.
