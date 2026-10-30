# Report — Stage-3 phase T8 (default-flip evaluations, release pipeline, records)

Date: 2026-10-30. Companion: `dod-T8.md`. Plan §12; amendment S3-A11.
T8 is the closing phase: the two §2.7 default-flip evaluations, the p1
lanes for the Stage-3 kernels, the sanitizer matrix and release pipeline
extended to the Stage-3 labels, the coverage lockstep, the study
migration note, and the final records. Owner commits (O10).

## 1. The default-flip outcomes (§2.7, O9/O4)

| Flip | Verdict | Where |
|---|---|---|
| `groundwater.reallocation` → `disabled` (W8.1, O9) | **FAILS — `legacy` stays.** 14 tests fail with unchanged tolerances. Decisive physics: b6's interface MAE vs the Kuan experiment is 13× worse under `disabled` (ss 0.329 m vs the 0.025 m golden), tidal salt mass +140 %/+213 %, b3's mass balance opens to 12.47 % (bound 0.5 %), the dtg controller stalls (b6_gw_smoke timeout). `disabled` ships as the documented vadose-domain option. | `stage3/df/realloc-report.md` |
| `transport.surface_dryout` → `store` (W8.2, O4) | **PASSES — KEPT.** 60/60 regression, 138/138 unit/mpi/stage3 (after the pin fix below), perf green. The decisive ebb-drying b6 cases stayed inside their golden bands with real margin: ss salinity 3.51 % of cells outside (≤ 6 %), salt mass +5.47 % (≤ 10 %); td 4.82 % and +8.79 %. td's interface MAE vs the experiment (0.0389 m) is BETTER than the discard-mode golden's own 0.0430 m — the store is a fidelity improvement on the tidal lane, not a tolerated drift. | `stage3/df/w82-*.log` |

Per §2.7-5 the combined run is owed only when both flips pass
individually; it is not. W8.3: the default-bitwise pins keep explicit
legacy values (below), `parameters.md`'s default column, the theory
section, and `ConfigTest.SurfaceDryoutDefaultIsStoreAfterTheT8Flip`
(shipped-default + explicit-discard round trip) landed; the ET+salinity
store rule now fires only on an EXPLICIT `discard`.

## 2. What the evaluations caught: three flip-survival defects (fixed)

1. **`BITWISE_PINS` was empty.** The §2.6 table — "v2.0-equivalent pins
   injected so the default-bitwise rule survives the T8 flips" — was
   never populated by T1/T5. Under the flipped binaries,
   `regression.s3.default_bitwise` (b2, g6c-coupled) and et9's noET pin
   lane went non-bitwise: those gates were testing *the shipped
   default*, not *the v2.0 configuration*. Now populated
   (`reallocation: legacy`, `surface_dryout: discard`, module-gated,
   printed per lane), and `et9-noet.yaml` pins `discard` in-file.
2. **Four unit tests escaped the W1.6 explicit-pin rule** and asserted
   legacy reallocation semantics through the default (two CoupledModule
   seepage/vent tests, two TransportModule ledger tests). Pinned
   `Reallocation::Legacy` at their config sites.
3. **The `Config.cpp` parse sites had hardcoded fallbacks**
   (`"legacy"`/`"discard"`), so §2.7's "change the single constant in
   Config.hpp" would silently never reach a run. The absent-key path now
   keeps the struct initializer — verified with `--resolve` on both DF
   binaries; the r1 resolved-config hook is what proved each flip
   reached every run (§2.7-6).

Two further comparison-definition adjudications of the same family
(S3-A11): et9's noET bitwise pin is recorded at 1 thread on the default
aij path, so it is skipped under `--mat-type` (the W8.4 p1 lanes, OMP 4
+ backend override) and under `FREHG_SANITIZED_BUILD=1` (instrumented
codegen under the same compiler — the h5_bitwise fingerprint is compiler
identity only); every tolerance-based et9 criterion still gates in both
contexts and passes.

## 3. Release-pipeline scope extensions (W8.4–W8.9)

- **W8.4**: `regression.s3.p1.{ob3,pd2coupled,et9}.omp4.{aij,aijkokkos}`
  (label `stage3_nightly`; the harness gained `--pd2-lane coupled`). All
  six green with each gate's own criteria.
- **W8.5**: `check_forbidden.sh` clean over the full Stage-3 tree (the
  MemSpace/backend rules cover `atm/PenmanMonteith.hpp` and the T3–T7
  kernels); the p5 negative compile stays green in the unit tier.
- **W8.6**: `run_sanitizers.sh` now runs the `stage3`/`stage3_local`
  labels in all three matrix tiers. Full matrix: **184 tests, zero
  sanitizer findings** (the single ctest failure was the et9 pin
  definition above, green once guarded).
- **W8.7**: BC matrix +2 canonical cells (`edge-bc × periodic-{x,y}-axis`,
  schema-rejected, both probed by regression.s3.pd6 — the y-axis probe
  added here); `check_bc_matrix` OK (122 rows/cells);
  `check_feature_coverage` OK (127 rows, all T1–T7 features and
  interactions present, zero empty cells).
- **W8.8**: doxygen 0 warnings; `mkdocs build --strict` green;
  `check_parameter_docs` OK (235 schema tokens).
- **W8.9**: `ci_release_gate.sh` gained the Stage-3 tiers (step 2b:
  `stage3|stage3_local`; step 4 rider: `stage3_nightly`) — a release
  pipeline that skipped them would ship Stage 3 release-unexercised (the
  §8.3 v1-wind lesson applied to the pipeline itself). Run on the final
  build: steps 1–4 green (unit 77, mpi 11, fast regressions 31, per-PR
  regression 61 incl. the 6 p1 lanes, stage3 39, the b5 24 h record, the
  nightly tier 10); the first pass aborted at the et9 p1 pin (adjudicated
  above), and the resumed run re-executed stage3_nightly and
  scaling_nightly to completion (`w89-release-gate.log`,
  `w89-resume.log`). p4: no nvcc on this machine — Stage 3 CHANGED
  device-facing sources, so cuda-compile.yml must be re-verified green on
  the pushed Stage-3 SHA before any official tag (owner action at the
  sync).

## 4. The b5 record, once more

The release gate's b5 rain/sync 24 h record on the final T8 binary
(store default, legacy reallocation): discharge peak +2.6 %, ponding
−1.5 %, subsurface residual 3.9e-8 m³ — numerically identical to every
Stage-3 record since T2. The kept flip does not touch it (b5 carries no
transport), and the rejected flip never reached it (reverted).

## 5. Study migration (W8.10)

`stage3/records/T8/study-migration.md`: the tide BC gains
`inflow: {closure: external}` (kills the report_V9 open-boundary
momentum import); AX runs switch to `domain.periodic.x` with the x-wall
BCs and `*_xe` corner soil types removed and an x-periodic random field
from the generator; the ET/vegetation/PM recipe; `store` needs no key
now; the new observability outputs. No study file was edited.

## 6. Residual limitations

- The reallocation default remains `legacy` with its known net-creation
  term (report_V2; +2.80 m³ per 2 tidal cycles on c_H2_2D_U00) —
  measurable per run via `realloc_created`/`realloc_dropped`; C-RL
  redesign directions recorded for the owner in the DF report.
- GPU execution stays experimental pending the owner p6 bundle; p4 needs
  the CI re-verification above.
- The et9 noET bitwise pin (and the default-bitwise hashes) are defined
  on this machine's g++-16 uninstrumented serial build only; other
  contexts gate on tolerances.
- s1's 70 %@4 headroom clause remains honestly unenforced on the
  hardware in hand (v2 plan note, unchanged by Stage 3).

## 7. Amendments this phase

S3-A11 (plan §18): both DF verdicts with measured values, the three
flip-survival defects, the two pin-definition adjudications, the
pipeline scope extensions, the pd6 y-axis probe note.
