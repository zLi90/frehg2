# Report — Stage-3 phase T12 (saturated-contact classification)

Date: 2026-11-16. Plan §12D (O19), amendment S3-A19. DoD: `dod-T12.md`.

## 1. Why T12 existed

The owner pushed the T11 tree (official ee265cb). The `build` lane went red
on the x86-64 runner. Its only failure was
`regression.b5_rank_invariance.strict`: one coupled b5 step, strict solve,
n=1 vs n=2 differing by 3.704e-7 against the V2-A9 bound 1e-9. On the dev
machine the same lane passed at 1.0e-14 / 5.7e-14, so the failure looked
platform-specific.

It is not. The post-allocation classification (`gw_reallocate_classify`)
decides "saturated neighbour" by reading the neighbour's θ (and, for the
gradient split, h) **in place**, while the same kernel writes θ and h:
- **Restore chains in loop order.** A fringe cell (head above air entry,
  θ just below θs after the corrector) is restored to θs and then looks
  saturated to the cells visited after it. The restore chains along the
  ascending loop index (Gauss–Seidel).
- **Chains break at seams.** The halo copy is pre-sweep, so the chain
  breaks at every subdomain's west/south seam.
- **Threads race.** Under threads the reads are a race.
- **Clamp threshold.** θ ≥ θs thresholds a clamped value at its clamp,
  so deep-saturated cells sit on it to within the solve residual.

On the dev machine the T11 binary passes at n = 2/4/8 only by luck. At
n = 3/5/6 it measures 2.76e-6 / 2.37e-5 / 2.76e-6, with every diverging
cell on a seam column. The defect dates from the v1 port. T11's flux
changes multiplied the fringe cells until the x86 runner tripped at n=2.

## 2. What landed

- **The key.** `groundwater.saturated_contact: legacy | head`, default
  **`head`** after the W12.4 flip.
- **The `head` classification.** It reads every neighbour from a
  pre-sweep snapshot of θ and h (halos included). A neighbour counts as
  saturated when θ ≥ θs or h > h_ae, which is the state θ := θ(h)
  assigns. No new constant.
- **`legacy` is the v1 test, bit for bit.** It is pinned by
  `BITWISE_PINS` and the two et9 case files.

Candidates measured and rejected (b5, one step, §12D table):
- a snapshot with the exact θ test: still rank-dependent, through the
  θs threshold;
- a snapshot with θ ≥ θs(1 − 1e-10): invariant only under the strict
  solver.

The adopted rule changes 5 473 cells by at most 2.4e-4 per step on b5,
and lowers the step's net restore volume by 18 %.

## 3. Gates (`stage3-sc.yaml`)

| Gate | Result |
|---|---|
| sc0 (GwSc ×3) | mirror symmetry bitwise on x/y/z; whole fringe restored; the legacy one-way chain recorded. Red: 36 failures on the legacy rule |
| sc1 | head ≤ 1.04e-14 at n = 2, 3, 4, 5, 6, 8; restore 1.91e-2 m³; legacy anti-vacuity 2.37e-5 at n = 5 |
| sc2 | 1 vs 4 threads bitwise |
| battery | 13 checks; every corruption caught |

sc1 is labelled `stage3;mpi`, so the official CI's `ctest -L mpi` step
runs it on the x86 runner.

## 4. W12.3 / W12.4

- **W12.3.** Under the `head` override, 54 of 55 Stage-3 per-PR lanes
  pass. et9's noET pin breaks by design and is pinned `legacy`.
- **W12.4.** With `head` as the default and tolerances unchanged, every
  tier passes (S3-A19 table), so the flip stays.

Legacy-golden movement:

| Gate | Quantity | Before (legacy rule) | `head` |
|---|---|---|---|
| b3 | distances | identical to report precision | same |
| b5 24 h | discharge integral | +2.4 % | +2.8 % |
| b5 24 h | ponding integral | −2.5 % | −2.3 % |
| b6-ss | salt mass | 6.66 % | 6.79 % |
| b6-td | salt mass | 9.89 % | 9.20 %, away from its 10 % bound |

## 5. Two earlier adjudications were measuring this defect

- **gw-orient mixed class (S3-A18 finding 4).** The rot180 /
  antitranspose signal (wc 1.056e-3 at t = 300, head 7.6e-8 at t_end) was
  diagnosed as a scheme-independent transient whose decay the consistent
  pair slows, and the bound was relaxed to 5e-7. The signal is the
  classification: a mirror reverses the sweep direction, a transpose does
  not. Under `head` it measures 4.4e-16 at t = 300 and ≤ 1.8e-15 at
  t_end. The bound is back at 1e-12.
- **b5 default-lane rank bounds (S3-A18 finding 3).** These were widened
  on seam-affected numbers. Under `head`, n=2 measures 5.1e-9 subsurface
  (it was 3.06e-5). The A14 bounds are restored.
- **V2-A9's "platform-rounding flip".** It was this defect's seam
  footprint. The correction is appended to the v2 plan. The b5 strict
  bound stays at 1e-9 until the x86 runner reports sc1 under `head`.

## 6. Found on the way, not fixed here

The static kernel audit (`records/T12/kernel-audit.md`) found no other
in-place neighbour read. It did find two scatter races:
`swe_velocity_bc_rhs` and `swe_outflow_rhs` read-modify-write `Srhs` for
duplicated corner entries of side-face BC lists. They are exact in serial
and lose an update under threads. They are proposed as a separate task.

## 7. State at close

- **Dev tree (uncommitted, per O10).** T12 is in the dev tree, snapshot
  `refs/stage3/T12`.
- **Official repo.** Synced by file copy (`records/T12/sync.md`); the
  owner commits and pushes.
- **Expected CI result.** The `build` lane's b5 strict test should pass
  on x86: it runs under the `head` default, and sc1 gates the same
  scenario at 1e-12.
- **`cuda-compile`.** It stays red for the separate T11 nvcc
  first-capture issue, which is out of T12's scope. T12 adds no
  `if constexpr`.
