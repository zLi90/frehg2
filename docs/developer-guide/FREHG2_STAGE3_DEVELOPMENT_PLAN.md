# FREHG2 Stage-3 Development Plan (implementation-ready)

Date: 2026-10-27. Status: **awaiting owner approval; no source changed.**

This plan supersedes, for execution, `FREHG2_STAGE3_PLAN_2026-10-27.md`. That file stays as the assessment record: premise corrections P1–P15 and the reasoning behind each design choice. This document is the binding specification. It covers:
- what to build, in which files and in which order;
- which gates are authored first, with their exact configuration, reference, tolerance, budget and negative check;
- what "done" means for each phase.

Line numbers refer to the **T0 baseline tree**. That is the dev repo `~/Codes/frehg2-upgrade/frehg2` after step W0.3 reverse-applies the uncommitted stagefix. Its `src/` is then expected to be byte-identical to `2026Frehg2/frehg2-2.0.0/src`: on 2026-10-27 the two trees differed only in the seven stagefix source files, and the reverse patch applies cleanly. W0.4 verifies the identity. Line numbers drift as work proceeds, so always grep for the quoted code.

Notation:
- **W** = work item.
- **G** = gate.
- Phases are **T0–T8**.
- Gate families:
  - `rl*`: reallocation
  - `ob*`: open boundary (Feature A)
  - `pd*`: periodic (Feature B)
  - `dg*`: diagnostics/output
  - `et*`: evapotranspiration (Feature C)
- ctest names are `regression.s3.<gate>` unless stated.
- ctest labels:
  - `stage3`: per-PR tier;
  - `stage3_nightly`: phase-close tier;
  - `stage3_local`: bitwise, dev-machine only.

---

## 0. Owner decisions adopted

| # | Decision (owner, 2026-10-27) | How this plan implements it |
|---|---|---|
| O1 | Keep `zero_gradient` as the default inflow closure | Default unchanged. The four loosened tolerances in the dev tree are reverted at T0 (W0.3). Every existing case is bitwise unchanged unless it opts in. |
| O2 | Stage-deposit routing tied to `closure: external` (accepted with O1) | `external`: deposits on stage cells are removed before the velocity update and booked. `zero_gradient`: booked only (dynamics bitwise). §6 W2.4. |
| O3 | Periodicity axis-generic (accepted) | `domain.periodic: {x, y}`. Both axes are implemented and gated (§7). |
| O4 | Salt concentrates in the remaining water; **when the water is completely removed, the salt is stored and released when the cell becomes wet again** | Subsurface: E and T remove water only, and salt stays in the pore water (the supply cap keeps θ ≥ θ_r > 0, so the water is never completely removed). Surface: a new per-cell **salt store**. A cell that dries deposits its whole remaining salt mass into the store; the store is released in full into the cell's water at the first step the cell is wet again. Key `transport.surface_dryout: delete \| store` (§9 W5.7). Default `delete` until the T8 default-flip evaluation (§2.7, owner item §17-1). The ET block with salinity requires `store`. |
| O5 | Vapour-flux family (accepted) | `bulk \| penman_monteith \| prescribed`, one shared aerodynamic resistance, Beer–Lambert parallel two-source partition, Feddes uptake (§9–§11). |
| O6 | Defer PM reading the simulated temperature | The schema rejects `penman_monteith` with `modules.temperature: true`. Bulk with temperature keeps v2.0 behaviour (documented). |
| O7 | Transpiration continues while inundated | The T sink is independent of the wet/dry classification. Feddes h₁/h₂ (optional) give a too-wet cutoff for non-wetland types. |
| O8 | Dry-out residue | Superseded by O4 (salt store). |
| O9 | User can **disable reallocation**; complete the gates with it disabled; **if all tests pass, disabled becomes the default** | `groundwater.reallocation: legacy \| off` (T1, §5). Every Stage-3 gate with a groundwater module runs `off`, except the declared legacy lanes. The default flips at T8 under the default-flip protocol (§2.7), only if every test passes with unchanged tolerances. |
| O10 | Work in `frehg2-upgrade`; the owner commits at the end | No commits by the agent. Phase boundaries are recorded as snapshot refs outside any branch (§2.2). The owner's index, branch and stash are never touched. |
| O11 | Monitor decimation in scope | `output.monitor_interval` (T4). |
| O12 | Keep per-PR / nightly tiering | §2.5. |
| O13 | Feddes h₃/h₄ required per vegetation type, no physiology defaults | Schema-required. h₁/h₂ are optional, and omitting both disables the anaerobic reduction. |
| O14 | Study reruns | Out of this plan. T8 delivers a case-migration note for the study. |
| O15 (2026-10-31, revises O4) | Thin-layer removals park their **water together with its salt** in a per-cell store, released together at rewetting. Rain and evaporation (their additions, removals, dilution and concentration) are not affected by the store. Parked water is **inert**: it does not flow, infiltrate or evaporate, ET included. A new phase with its own gates implements it. | Phase T9 (§12A), amendment S3-A12. The salt-only store of O4 becomes the special case V_s = 0 (salt left by complete evaporative dry-out). The interim default of `transport.surface_dryout` is an open owner item (§17-6). |
| O16 (2026-10-31) | Surface salt (and heat) moves **with the water the free-surface solve moved**: the transport's face volumes are the solve's own continuity face volumes, not velocities rebuilt afterwards from the new depths and then limited. | Phase T10 (§12B), amendment S3-A14. |

---

## 1. Scope and phase map

| Phase | Content | Depends on | Effort (agent-days, predicted) |
|---|---|---|---|
| **T0** | Baseline reconciliation, reference binary and hashes, harness skeleton, override mechanism, records | — | 0.5–1 |
| **T1** | `groundwater.reallocation: off`; per-cell reallocation/vloss diagnostics; conservation gates; first default-flip preview | T0 | 1.5–2 |
| **T2** | Feature A: selectable inflow closure (`zero_gradient` default, `external`), stage-deposit routing, entrance-head diagnostics | T0 | 1.5–2 |
| **T3** | Feature B: periodic topology (x and y) | T2 (shares `WetDry`/`FreeSurface` edge code) | 2–3 |
| **T4** | Interval-mean flux outputs (`seepage_rate` …), monitor decimation | T0 | 1 |
| **T5** | Feature C core: ET block (`bulk`, `prescribed`), two evaporation channels with time-basis complementarity, single vapour flux, supply cap, surface salt store, transport windows | T1, T3 (et9 uses periodic x), T4 | 2.5–3 |
| **T6** | Vegetation and transpiration (Feddes sink, root fractions, Beer–Lambert split, salt exclusion) | T5 | 2–3 |
| **T7** | Penman–Monteith closure (the r_a laws ship in T5), incoming-shortwave input, FAO-56 verification | T6 | 1.5–2 |
| **T8** | Default-flip evaluation (reallocation, salt store), full pipeline, coverage tables, sanitizer, docs, records | all | 1.5–2 |
| **T9** (added by S3-A12) | Paired water–salt store for thin-layer removals (O15): booking of every surface thin-layer removal, the paired store and its release, coupled sub-`min_depth` vents, re-adjudication of the store-resolving gates, DF re-evaluation of the `store` default | T8 | 2.5–3.5 |
| **T10** (added by S3-A14) | Solve-consistent surface transport (O16): the transport advects with the free-surface solve's continuity face volumes, so salt and heat move with the water; post-solve manipulations (limiters, removals, sources, exchange) keep their own explicit books | T9 | 1.5–2.5 |

Execution order: T0 → T1 → T2 → T3 → T4 → T5 → T6 → T7 → T8. The total is ≈ 14–19 agent-days. T9 follows the closed T8 (S3-A12, 2026-10-31); T10 follows the closed T9 (S3-A14, 2026-10-31).

---

## 2. Working rules (binding)

### 2.1 Environment

| Item | Value |
|---|---|
| Source tree | `~/Codes/frehg2-upgrade/frehg2` (git, branch `main`; the owner's 2026-10-05 uncommitted work is part of the baseline) |
| Build dir | `frehg2/build-s3` (gitignored via `build-*/`); T8 adds `build-s3-df` |
| Configure | `cmake -S . -B build-s3 -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16 -DCMAKE_PREFIX_PATH=/Users/zhili/Codes/local -DKokkos_DIR=/Users/zhili/Codes/local/lib/cmake/Kokkos -DHDF5_DIR=/Users/zhili/Codes/local/cmake -Dyaml-cpp_DIR=/Users/zhili/Codes/local/lib/cmake/yaml-cpp -DMPIEXEC_EXECUTABLE=/Users/zhili/Codes/local/bin/mpiexec -DFREHG_ENABLE_TESTS=ON -DFREHG_WERROR=ON` (copied from `build-stagefix/CMakeCache.txt`) |
| MPI launcher | **always** `/Users/zhili/Codes/local/bin/mpiexec` (MPICH 4.3.0). The PATH `mpiexec` is Homebrew Open MPI and silently degrades the world to size 1. |
| Run environment | `OMP_NUM_THREADS=1 OMP_PROC_BIND=false FI_PROVIDER=tcp UCX_TLS=tcp,self,sm` (the ctest `_regress_env`) |
| Legacy goldens | `~/Codes/frehg2-upgrade/legacy/benchmarks` (`FREHG_LEGACY_BENCHMARKS` default) |
| Records root | `~/Codes/frehg2-upgrade/stage3/` (outside the git tree; layout in §2.2) |
| Read-only | `2026Frehg2/` (whole study tree, including `frehg2-2.0.0/`, `cases_dev/`, `outputs/`). Copy anything needed into `stage3/`. |

### 2.2 Repository discipline (no commits)

1. Never `git commit`, `git stash`, `git checkout -- <file>`, `git reset`, or modify the index. Edits go to the working tree only.
2. **Phase snapshot**: at the start of T0 and at the end of every phase, run `stage3/tools/snapshot.sh <label>`. The script (authored in W0.1) does this:
   ```bash
   set -euo pipefail
   cd ~/Codes/frehg2-upgrade/frehg2
   tmp=$(mktemp)
   export GIT_INDEX_FILE="$tmp"
   git read-tree HEAD
   git add -A                      # respects .gitignore; private index only
   tree=$(git write-tree)
   c=$(git commit-tree "$tree" -p HEAD -m "stage3 snapshot $1 (not on any branch)")
   git update-ref "refs/stage3/$1" "$c"
   unset GIT_INDEX_FILE; rm -f "$tmp"
   echo "$1 -> $c"
   ```
   The phase diff for review is `git diff refs/stage3/T<n-1> refs/stage3/T<n>`. It is written to `stage3/records/T<n>/phase.patch` together with `--stat`. The owner can delete all refs with `git for-each-ref --format='%(refname)' refs/stage3 | xargs -n1 git update-ref -d`.
3. **Records layout:**
   ```
   stage3/
     tools/snapshot.sh
     seeds/stagefix_v2.patch            # copy of the 2026-10-26 prototype (T2 seed)
     ref-bin/frehg-T0, ref-bin/deps.txt # frozen T0 binary + otool -L
     records/T<n>/red/*.log             # gate-first failure evidence
     records/T<n>/green/*.log           # final gate output
     records/T<n>/tierP.log, tierN.log, timings.md, phase.patch, phase.stat
     evidence/<name>/                   # investigation runs (not gates)
     df/                                # T8 default-flip evaluation
   ```

### 2.3 Gate-first protocol

For each phase:
1. Register the phase's gates (cases, tolerance file, harness subcommand, ctest entry, negative battery) **before** writing capability code.
2. Run them on the pre-capability binary. Save the failure (or "capability absent" from `--validate`) to `records/T<n>/red/`.
3. Implement until green, then save green logs.

Rules for registration:
- Gates of **future** phases are not registered until their phase starts, so the per-PR tier stays green between phases.
- Every gate passes the v2 §6.3 self-checks at authoring:
  - (1) the reference is checked against ≥ 2 independent spot values;
  - (2) a negative check proves the gate can fail;
  - (3) a refinement study exists where the metric is discretisation-driven.
- **Anti-vacuity rule (binding, from V2-A9…A12).** Every gate must contain an assertion that its scenario actually exercises the feature. For example: rl2's legacy lane must show reallocation activity above a floor; ob1's `zero_gradient` lane must differ from `external` by ≥ 2×; et3 must contain both channels. A gate whose scenario can pass with the feature disconnected is defective.

### 2.4 Process rules

- **Amendment log** `S3-A<n>`: appended to the in-repo copy of this plan (`docs/developer-guide/FREHG2_STAGE3_DEVELOPMENT_PLAN.md`, created in W0.6) for every deviation. That includes any tolerance change, gate redesign, scope change or contingency use. Each entry carries the measured data. No silent tolerance edits.
- **One `dod-T<n>.md` and one `report-T<n>.md`** per phase in `docs/developer-guide/`, plus mkdocs nav entries.
- **Forbidden patterns:** `scripts/check_forbidden.sh` stays clean. Note that it bans the words TODO, FIXME, XXX, HACK, WIP, placeholder, stub, "not implemented", "for now" and temporary/temporarily **anywhere in `src/` including comments**. T3 adds a new rule (§7 W3.9).
- Zero-warning builds (`-Werror`).
- **p1/p5 green for every new kernel**: the backend-invariance lanes are in T8.
- **Lockstep checkers stay green on every work item:**
  - `unit.parameter_docs`: every new key is documented in `docs/user-guide/parameters.md`;
  - `unit.bc_matrix`;
  - `unit.feature_coverage`: new rows in `docs/developer-guide/feature-coverage.md`;
  - `unit.run_record_checker`;
  - r1: every schema change extends `resolvedConfigYaml` (`src/io/RunRecord.cpp`) in the same work item, or the round-trip fails.
- New benchmark YAMLs are added to the `validate.<case>` foreach list in `tests/CMakeLists.txt`.
- Code style: match the surrounding code (doxygen `///` headers, provenance comments, Kokkos `MDRangePolicy` idiom, `Field2/Field3` views captured by value).

### 2.5 Re-verification tiers

| Tier | Command | When | Measured wall (M3, 8 cores) |
|---|---|---|---|
| **P** (per work item / "per PR") | `ctest --test-dir build-s3 -L '^(unit\|mpi\|regression\|stage3\|stage3_local)$' -j8 --output-on-failure` | After every work item that touches `src/`, `tests/` or benchmarks | Existing part ≈ 23 min (1360 s, `suite93` log) + g1/p1 lanes 3.6 min if run; Stage-3 adds ≤ 5 min |
| **N** (phase close) | Tier P, plus release steps 3–4: `run_regression.py b5 --scenario rain --coupling sync --t-end 86400 --ranks 4 …` and `ctest -L regression_nightly -E '^regression\.b5\.(rain\|norain)\.'`, plus `ctest -L stage3_nightly` | End of every phase, before the snapshot | ≈ 56 min for the 7-test nightly batch (measured) + b6 ss/td bjacobi + b5 24 h record (timed at T0) |
| **R** (release) | `scripts/ci_release_gate.sh build-s3` (adds `scaling_nightly`) + `scripts/run_sanitizers.sh` over new labels + doxygen + `mkdocs build --strict` | T8 | timed at T8 |

A phase is closed only with Tier N green and its timings recorded in `records/T<n>/timings.md`.

### 2.6 Default-bitwise rule

With **no Stage-3 key**, every case in the default-bitwise list (§3.3) produces HDF5 field and monitor datasets bitwise equal to the T0 reference hashes.

Declared, allowed exceptions:
1. Newly appended audit columns: `mass_audit.stage_deposit` (T2) and later appended columns.
2. `transport_audit.surf_boundary` values: the T2 edge-diffusion audit fix (audit-only).
3. New run-record lines.

Any other difference is a defect. From T1 on, the harness injects the "v2.0-equivalent pins" (`groundwater.reallocation: legacy`, and from T5 `transport.surface_dryout: discard` (token per S3-A8.10)), so the rule survives the T8 default flips.

**Comparison caveat (FP contraction).** If a bitwise lane differs by ≤ 1e-14 relative only in a kernel edited in the current phase, rebuild both trees with `-ffp-contract=off` to test whether FMA contraction caused it. Record the outcome by amendment, never by editing hashes silently.

### 2.7 Default-flip (DF) protocol (implements O9; proposed also for O4)

1. **Preview (T1, non-binding).** Run the YAML-driven tiers with `FREHG_STAGE3_OVERRIDES="groundwater.reallocation=off"` (§3.2). Record which existing gates fail in `records/T1/df-preview.md`. This informs the owner early and changes nothing.
2. **Evaluation (T8, binding).**
   - Change the single default constant in `src/core/Config.hpp`, e.g. `Reallocation reallocation = Reallocation::Off;`.
   - Build `build-s3-df` and run **every** test: unit, mpi, Tier N and `stage3`.
   - Tolerances stay unchanged.
3. **Rule.** If **all** pass, the flip stays: amendment, docs (`parameters.md` default column, theory page), and the resolved-config default. If **any** fails, revert the one-line change, keep `legacy`, and write `stage3/df/report.md`. That report gives, per failing gate, the measured value against its bound and a diagnosis: is the `off` result closer to the analytic/experimental reference than `legacy`, or is the gate's golden itself reallocation-dependent? The owner then decides.
4. **Unit tests asserting legacy reallocation semantics** are pinned to `Reallocation::Legacy` explicitly in T1 (W1.6). The flip therefore only affects tests that do not test reallocation itself. This is legitimate because `legacy` remains available.
5. **Salt store (O4).** The same protocol applies with `transport.surface_dryout: store`, evaluated after the reallocation evaluation. If both flips pass individually, a combined run must also pass all tests before both are kept. If the combined run fails, keep the reallocation flip only, per O9's priority, and report.
6. The DF evaluation must prove the change reached every run. The r1 hook asserts that each run record's resolved config carries the flipped value (W0.5).

### 2.8 Gate configuration rule for O9

Every Stage-3 gate case with `modules.groundwater: true` sets `groundwater.reallocation: off` explicitly. The only exceptions are:
- the declared legacy lanes (rl2/rl3 legacy lanes);
- the default-bitwise lanes, which pin `legacy`.

---

## 3. Shared infrastructure (built in T0)

### 3.1 Harness layout

| File | Content |
|---|---|
| `tests/regression/gate_realloc.py`, `gate_openbc.py`, `gate_periodic.py`, `gate_diag.py`, `gate_et.py` | **Pure** reference and check functions (no HDF5 I/O), as in `gate_wind.py`. Constants mirror the source (`AIR_DENSITY = 1.225`, `WATER_DENSITY = 998.0`, `GRAVITY = 9.81`, Garratt law). |
| `tests/regression/run_regression.py` | New subcommands `rl`, `ob`, `pd`, `dg`, `et`, `default-bitwise`, `s3-restart`, each with a `--lane` option. HDF5 plumbing and staging reuse `stage_case`, `run_case`, `rewrite_config`, `run_dihedral_battery`. |
| `tests/regression/tolerances/stage3-{rl,ob,pd,dg,et}.yaml` | Every numeric bound, with a comment citing the gate and the authoring measurement |
| `scripts/test_{rl,ob,pd,dg,et}_gates.py` | Negative batteries: every check function fails on corrupted synthetic input; reference spot values. Registered as `unit.<x>_gates.negative`. |
| `benchmarks/rl-realloc/`, `ob-openbc/`, `pd-periodic/`, `dg-diagnostics/`, `et-evapotranspiration/` | Cases + `README.md` in the `benchmarks/g9-wind/README.md` template (reference, tolerance citation, steadiness rule, refinement study, negative battery) |
| `tests/regression/reference-hashes/stage3-T0.json` | §3.3 |
| `tools/h5_bitwise.py` | `hash <h5> --out json` and `compare <h5> <json> --allow <dataset-glob>...` |

### 3.2 Override mechanism (`FREHG_STAGE3_OVERRIDES`)

- `run_regression.py` reads the environment variable `FREHG_STAGE3_OVERRIDES` (`;`-separated `dotted.key=value`).
- `apply_solver_override()` (already called by `stage_case`) additionally applies each override to every staged YAML whose modules make the key meaningful:
  - `groundwater.*` needs `modules.groundwater`;
  - `transport.*` needs `modules.transport`.
- `groundwater.reallocation_surplus` and `groundwater.density_coupling.reallocation_gradient` are left in place: with `off` they only raise the W1.1 warning. This is deliberate, so that the T8 default flip runs every existing YAML (b3, g5 and the study cases carry these keys) without edits.
- `check_run_record()` passes `--expect key=value` to `tools/check_run_record.py` for each active override. The checker asserts that the resolved config carries it, so an override that missed a run fails that run (anti-vacuity).

### 3.3 Default-bitwise gate (`regression.s3.default_bitwise`, label `stage3_local`)

- **Case list** (fast per-PR cases; each ≤ ~10 s): `b1-sw/b1-sw`, `b2-gw/b2-gw`, `g4-evap/g4a-drawdown`, `g4-evap/g4b-concentration`, `g4-evap/g4d-drain`, `g6-heat/g6c-coupled`, `g7-heat/g7b-channel`, `g9-wind/g10-seiche`, `x-orient/{swe,gw,transport}-orient-base` (identity orientation only), `x-orient/conveyor-base`, `x-orient/coupled-salt-base`.
- **Procedure:** stage, inject the v2.0-equivalent pins (§2.6), run, then `tools/h5_bitwise.py compare` against `stage3-T0.json` with the declared allow-list.
- **Coverage of the list:** SWE with eta, coupled with eta and salt, groundwater-only reallocation, transport, heat.
- **Hash scope:** SHA-256 per dataset of `/surface/*`, `/groundwater/*`, `/transport/*`, `/temperature/*`, `/monitor/*` (point monitors and audit tables).
- **Toolchain fingerprint** in the JSON: compiler path and version, `CMAKE_BUILD_TYPE`, PETSc/Kokkos versions, `uname -m`. A mismatch exits 77 (`SKIP_RETURN_CODE 77`): the gate is meaningful only on the dev toolchain.
- **Predicted cost:** ≈ 40 s.

### 3.4 Restart gate (`regression.s3.restart`, grows per phase)

- Each lane runs a case continuously and with a restart at mid-time (`output.checkpoint.interval` + `restart:`), then requires bitwise equality of all output datasets at t_end.
- Lanes are added by T3 (periodic surface), T4 (interval accumulators), T5 (dry-time, store, ET) and T6 (vegetated ET).
- Cost: ≤ 3 s per lane.

---

## 4. T0 — Baseline, reference, harness

**Objective.** A clean, measured starting point that equals the study's production source. The work also creates every shared tool the later phases use.

| W | Action | Detail |
|---|---|---|
| W0.1 | Records root | `mkdir -p ~/Codes/frehg2-upgrade/stage3/{tools,seeds,ref-bin,records,evidence,df}`; author `tools/snapshot.sh` (§2.2) |
| W0.2 | Start snapshot | `stage3/tools/snapshot.sh T0-start` (captures the tree including the stagefix and the 10-05 work) |
| W0.3 | Remove the prototype | `cp 2026Frehg2/outputs/debug_v1/impl_v1/stagefix_v2.patch stage3/seeds/`; `git apply -R --check stage3/seeds/stagefix_v2.patch` (verified clean 2026-10-27); then `git apply -R stage3/seeds/stagefix_v2.patch`. This reverts Fix A/B, the edge-diffusion audit fix, the two stagefix unit tests, the docs hunks and **all four tolerance amendments**. |
| W0.4 | Prove the baseline | `diff -r src "/Users/zhili/Documents/Work/个人/论文/2026Frehg2/frehg2-2.0.0/src"` must report only `.DS_Store`. `grep -n head_abs_floor tests/regression/tolerances/b6-kuan.yaml` must read 0.015 and `salt_mass_rel` 0.10. The b5 bjacobi bounds in `run_regression.py` must read 2.0e-3 / 3.0e-5 / 2.0e-3. Record the outputs in `records/T0/baseline.md`. |
| W0.5 | Harness | §3.1–§3.3: override mechanism + `check_run_record.py --expect`; `tools/h5_bitwise.py`; empty `stage3-*.yaml` tolerance files with headers; `unit.<x>_gates.negative` battery skeletons. The skeletons for rl/ob/pd/dg/et contain only the shared-helper self-tests until their phase. |
| W0.6 | Plan in repo | Copy this file to `docs/developer-guide/FREHG2_STAGE3_DEVELOPMENT_PLAN.md`; add mkdocs nav entries ("Stage-3 development plan (binding spec)" and a "Stage-3 phase records" group). Copies of this plan in `frehg2-upgrade/` are not edited further. |
| W0.7 | Build + Tier N | Configure `build-s3` (§2.1), build, run Tier P and Tier N. Expected green (it equals the study's production source). A failure here is pre-existing: record it, diagnose it, and report it to the owner before T1. Time each tier, plus b6.ss, b6.td and the b5 24 h record individually, into `records/T0/timings.md`. |
| W0.8 | Reference binary + hashes | `cp build-s3/src/frehg stage3/ref-bin/frehg-T0`; `otool -L … > deps.txt`; run the §3.3 list, hash into `tests/regression/reference-hashes/stage3-T0.json`, register `regression.s3.default_bitwise` (green by construction) |
| W0.9 | Schema warning channel | `ValidationResult` (`core/Config.hpp:533`) gains `std::vector<std::string> warnings`; the schema context gains `addWarning(where, msg)`. `frehg --validate` prints each as `WARNING: <where>: <msg>` and still prints `VALID:`. At run start each is logged with `log::warn` and written to a new run-record list `warnings:` (`io/RunRecord.cpp`; `tools/check_run_record.py` accepts the key; `tools/test_run_record_checker.py` gains a case). Unit test `ConfigTest.WarningsDoNotInvalidate`. Every "warning" in this plan (§13) uses this channel; runtime-only warnings (seam, root truncation, entrance head) call `RunRecord::addWarning` as well as `log::warn`. |
| W0.10 | Close | `snapshot.sh T0`; `dod-T0.md`, `report-T0.md` (baseline proof, timings, hash list) |

**DoD T0:**
- [ ] W0.4 diff is empty;
- [ ] Tier P and Tier N are green, with timings recorded;
- [ ] `stage3-T0.json` is committed in the tree (as a file; the owner commits);
- [ ] the override mechanism is unit-tested: a synthetic run record missing the override fails `check_run_record.py --expect`;
- [ ] the W0.9 warning channel is unit-tested and accepted by the run-record checker;
- [ ] `check_forbidden.sh` is clean;
- [ ] snapshot `refs/stage3/T0` exists.

---

## 5. T1 — Reallocation option (`groundwater.reallocation: legacy | off`)

**What "reallocation" is here.** It is the PCA post-allocation step `RichardsSolver::reallocateWaterContent()` (`src/gw/Reallocate.cpp`), called once per Richards step after the corrector (`RichardsSolver.cpp:713`, `RichardsSolver::step`). It does four things:
1. It sends over-saturation excess (θ ≥ θs) up or down the column along the head gradient and **drops the lateral share** (`Reallocate.cpp:212-226`).
2. In cells in saturated contact, it **restores θ := θ(h_pred)**. That **creates** water when θ(h) > θ (`:271-298`) and drops or redistributes the surplus otherwise.
3. In isolated unsaturated cells, it re-derives h := h(θ) (`:265-270`).
4. It drops remainders at column ends and sub-`min_depth` vents (`:380-403`, `:422`).

The final clamp `finalizeWaterContent()` (`Corrector.cpp:554`) is **not** reallocation. It stays in both modes, and its volume is audited as `vloss`.

**`off` semantics:**
- The corrector's flux-conservative θ stands.
- Every active cell with θ < 0.9999 θs takes h := `headFromWaterContent(soil, θ)` (`VanGenuchten.hpp:65`).
- Cells at or above 0.9999 θs keep the predicted h.
- `finalizeWaterContent()` then clamps θ into [θr, θs], booking `vloss` as today.
- Nothing is sent, restored, created or dropped. `realloc`, `realloc_dropped` and `cplVent` are exactly 0.

### Work items

| W | File(s) | Exact change |
|---|---|---|
| W1.1 | `core/Config.hpp`, `Config.cpp`, `ConfigSchema.cpp`, `io/RunRecord.cpp` | In `GroundwaterConfig`, add `enum class Reallocation { Legacy, Off };` and `Reallocation reallocation = Reallocation::Legacy;`. This one line is the T8 default-flip constant; give it a comment saying so. Schema: `{"reallocation", enumeration(false, {"legacy", "off"})}` in the `groundwater` map. Cross-field **warnings** (W0.9 channel; not errors, so that existing YAMLs still validate when T8 flips the default): `groundwater.reallocation_surplus: has no effect with reallocation: off (nothing is transferred)`, and `groundwater.density_coupling.reallocation_gradient: has no effect with reallocation: off`. Echo it in `resolvedConfigYaml`. |
| W1.2 | `gw/RichardsSolver.hpp/.cpp`, `gw/Reallocate.cpp` | Add member `bool reallocOff_`. In `step()`, replace `reallocateWaterContent();` with `if (reallocOff_) { restoreHeadFromWaterContent(); } else { reallocateWaterContent(); }`. New `RichardsSolver::restoreHeadFromWaterContent()` in `Reallocate.cpp`: a cell-parallel kernel over active cells, `if (wc < 0.9999 * thetaS) h = headFromWaterContent(soil, wc);`. It sets `audit_.reallocAdjust = audit_.reallocDropped = audit_.cplVent = 0.0`. `qzPathCorr_` is already zeroed at `RichardsSolver.cpp:711`; `ventWin_` gets no deposit. Before relying on it, grep `room_\|sendUp_\|sendDown_` to confirm that nothing outside `Reallocate.cpp` reads them. Update the file header to describe both modes. |
| W1.3 | `Reallocate.cpp`, `RichardsSolver.hpp`, `driver/Simulation.cpp`, `ConfigSchema.cpp` | Per-cell cumulative diagnostics, allocated **only when requested** as output (so default runs are untouched): `reallocCreated_` and `reallocDroppedCell_` (Field3, m³). Phase 1 books `max(0, after − before)·V` as created and the dropped lateral/surplus share as dropped on the cell. Phase 2 books the column-end and discarded-vent remainders on the **sending** cell k. New groundwater output variables: `realloc_created`, `realloc_dropped`, and `vloss` (the existing cumulative `vloss_` field). Add them to `groundwaterVariables()` and `groundwaterVarMeta` (units `m3`, long names). Checkpoint them as `realloc_created`, `realloc_dropped`, `vloss_cum` when requested. |
| W1.4 | `Simulation.cpp` | Run-record closure line `groundwater: reallocation: <legacy\|off>` |
| W1.5 | docs | `theory/groundwater.md`: new section "Post-allocation modes" with the four legacy actions above, the `off` definition, the anchor consequence (W1.7) and the default-flip rule. Also `parameters.md`, `configuration.md`, and `output.md` (new variables). |
| W1.6 | `tests/unit/*.cpp` | Pin every unit test that asserts reallocation semantics to `cfg.groundwater.reallocation = Reallocation::Legacy`, with the comment `// asserts the legacy post-allocation (kept selectable; Stage-3 §2.7)`. Find them with `grep -nE 'reallocAdjust\|reallocDropped\|cplVent\|ReallocationSurplus\|reallocat\|vent' tests/unit/*.cpp`. List the pinned tests in `report-T1.md`. |
| W1.7 | — | Documentation only: under `off` with `transport.subsurface_update: consistent`, `subs_anchor` reduces to the scalar on the `vloss` volume (gate rl4) |
| W1.8 | `stage3/evidence/T1-df-preview/` | DF preview (§2.7-1): `FREHG_STAGE3_OVERRIDES="groundwater.reallocation=off" ctest -L '^regression$'` plus b6.ss/b6.td in the background. Table of pass/fail and measured value against bound in `records/T1/df-preview.md`. Non-blocking. |
| W1.9 | `stage3/evidence/T1-h2-2d/` | Evidence run (not a gate): copy `2026Frehg2/cases_dev/matrix/c_H2_2D_U00/` (read-only source), set `t_end` to 2 tidal cycles (89424 s), request `realloc_created, realloc_dropped, vloss`, and run legacy and off. Report net creation, exfiltration, `subs_anchor` and seepage for both. Expected ≤ 1 h each; measure. |

### Gates (authored first)

| Gate | Configuration | Asserts | Tolerance (`stage3-rl.yaml`) | Budget | Negative / anti-vacuity |
|---|---|---|---|---|---|
| **rl1** unit (`test_gw.cpp`): `GwModule.ReallocationOffConservesWaterExactly`, `GwModule.ReallocationOffKeepsUnsaturatedHeadsOnTheRetentionCurve`; `ConfigTest.ReallocationModeParsesRoundTripsAndRejectsLegacyOnlyKeys` | 1×1×20 column: sand over silt (sand Ks 1e-4, θs .41, θr .05, α 4, n 2.5; silt Ks 1e-6, θs .45, θr .10, α 1, n 1.5), Ss 1e-5, top flux −5e-6 m/s, bottom head with the water table at mid-column; 200 steps, dtg 10 s | Every step under `off`: `reallocAdjust == 0.0` and `reallocDropped == 0.0` exactly. Δvolume = boundary_in − ss_storage − vloss to 1e-14 relative. Unsaturated cells: h == headFromWaterContent(θ) bitwise. Saturated-class cells have finite h ≥ h(0.9999 θs). Schema messages as in W1.1. | exact / 1e-14 | < 0.1 s | Legacy lane of the same column shows reallocAdjust ≠ 0 in ≥ 1 step (anti-vacuity) |
| **rl2** `rl-realloc/rl2-tidal-box.yaml` | 8×8×12, dx = dy = 1 m, dz = 0.25 (no stretch), flat base, column 3 m. Sand/silt (rl1 parameters) in a 2×2×3-block checkerboard (generated `input/soil_map.dat`, `(j*nx+i)*nz+k`). y−: `head` `hydrostatic: {eta: 1.5}`. y+: `head` `hydrostatic: {eta: 1.2}`. Top: `groundwater_top` `flux` series `input/top_flux.dat` (negative = infiltration, the b3 convention): −2e-5 m/s for the first 1800 s of every 3600 s, 0 otherwise, written as a piecewise-linear `time value` series with 1 s ramps. Bottom no-flux, IC water table 1.3 m, Ss 1e-5, dtg ∈ [0.1, 30] s, t_end 14400 s (four wetting fronts crossing the heterogeneity into the water table). Lanes `legacy` and `off`. | **off:** cumulative `realloc` = 0 and `realloc_dropped` = 0 (exact). Ledger \|Δvol − (bin − ss + realloc − vloss)\| ≤ 1e-12·T (T = Σ\|boundary_in\| increments). \|vloss_cum\| ≤ ε_v·T. θ ∈ [θr, θs] at outputs. substeps(off)/substeps(legacy) ≤ 1.5. **Both:** Σ per-cell `realloc_created − realloc_dropped` and Σ `vloss` equal the audit columns to 1e-12 (the dg3 check moved here). | ε_v = 1e-4 initially; adjudicated at authoring to ≥ 10× headroom over the measured value, else contingency C-RL (§16) | ≈ 6 s (2 lanes) | **Anti-vacuity:** legacy gross churn (created + dropped) ≥ max(1e-6·T, 10·\|vloss_off\|), else redesign the case |
| **rl3** `rl-realloc/rl3-steady-column.yaml` | 1×1×40, dz 0.05 (2 m), sand; top `flux` −2e-6 m/s (infiltration); `groundwater_bottom` `head` constant 0.0 (pressure head 0: water table at the base); IC hydrostatic; t_end 518400 s (6 d), dtg ≤ 600 s. Lanes `legacy` and `off`. | **off:** steadiness max\|Δθ\| over the last 10 % ≤ 1e-10; flux uniformity max_k \|qz_k − q_top\|/\|q_top\| ≤ 1e-8; realloc = 0. **legacy:** uniformity and churn recorded in `report-T1.md` (they show whether the restore distorts a steady front). | 1e-10 / 1e-8 | ≈ 3 s | Synthetic trace with a 1e-6 flux kink fails the uniformity check |
| **rl4** coupled salt anchor | `x-orient/coupled-salt-base.yaml` with `transport.subsurface_update: consistent` (set by the gate if absent); lanes `legacy` and `off` | **off:** \|subs_anchor_cum\| ≤ max(1e-12·M_subs, 1.0001·s_max·\|vloss_cum\|). **legacy:** recorded. | as stated | ≈ 2 s | Synthetic anchor 1e-6·M fails |
| **rl5** batteries under `off` | `regression.s3.rl5.gw_orient_off` (existing `gw-orient` battery) and `regression.s3.rl5.coupled_salt_rank_off` (existing strict coupled-salt rank invariance), both with `ENVIRONMENT FREHG_STAGE3_OVERRIDES=groundwater.reallocation=off` | Existing criteria (`tolerances/x-orient.yaml`; strict rank bound) | existing | ≈ 6 s | The r1 `--expect` assertion proves the override reached each run |

`unit.rl_gates.negative`: the closure, exact-zero, uniformity, anchor and anti-vacuity checks each fail on corrupted synthetic input.

**DoD T1:**
- [ ] rl1–rl5 green;
- [ ] `default_bitwise` green with the legacy pin;
- [ ] W1.6 pins listed;
- [ ] DF preview table filed;
- [ ] evidence run filed (or its measured runtime recorded if aborted at > 2 h);
- [ ] docs and lockstep checkers green;
- [ ] Tier N green;
- [ ] snapshot `T1`.

---

## 6. T2 — Feature A: selectable inflow closure at prescribed-stage boundaries

**Seed.** `stage3/seeds/stagefix_v2.patch`: its Fix A masks and stencils, Fix B removal and audit, edge-diffusion audit fix, and unit tests. T2 re-applies that logic in selectable form. Do not `git apply` the seed wholesale: it hard-wires the at-rest closure and loosens tolerances.

**Semantics.**
- `zero_gradient` (default): v2.0 arithmetic, bitwise. The upstream state of an inflow face is the mass-consistent reconstruction (`WetDry.cpp:364-403`).
- `external`: the stage cell is an **external reservoir** with velocity (u_ext, v_ext), default (0, 0).
  - Every **reservoir face** enters the advection stencils at the reservoir's state, and self-advection across it uses the conservative difference of u²/2.
  - A reservoir face is one whose two cells are both external-stage cells, or one external-stage cell and one cell outside the global domain.
  - Coupled deposits on the stage cell leave through the reservoir before the velocity update.
  - Flow rates, the continuity back-solve and the free-surface system are unchanged.

### Work items

| W | File(s) | Exact change |
|---|---|---|
| W2.1 | `Config.hpp/.cpp`, `ConfigSchema.cpp`, `RunRecord.cpp` | In `BoundaryConditionConfig`, add `enum class InflowClosure { ZeroGradient, External }; InflowClosure inflowClosure = InflowClosure::ZeroGradient; SeriesOrConstant inflowU, inflowV;` (constant 0 default). Schema `bcEntry`: `{"inflow", map({{"closure", enumeration(true, {"zero_gradient", "external"})}, {"velocity", map({{"u", seriesOrConstant(false)}, {"v", seriesOrConstant(false)}})}})}`. Cross-field errors are in §13.A. Echo in the resolved config. |
| W2.2 | `bc/BoundarySet.*`, `swe/SurfaceSolver.hpp/.cpp` | Carry closure and the u/v series on `BoundaryCondition` (evaluated with the same `value(t)` machinery as the stage). In `buildBoundaryLists`: (a) `stageClosure_` (Field2<int>: 0 none, 1 zero_gradient, 2 external), halo-exchanged. (b) Setup fatal if two stage cells sharing a face carry different closures. (c) Static masks `resFaceX_`/`resFaceY_` (Field2, 1.0 on external-reservoir x+/y+ faces; the seed's `reservoirFace` lambda restricted to closure 2), `exchangeWithCorners`. (d) Warning `eta boundary '<name>': closure external has no reservoir face (a one-cell-thick interior stage line); it behaves as zero_gradient`. Per step, where the stage `current` values are refreshed: fill `uResX_(j,i)` and `vResY_(j,i)` on reservoir faces with u_ext and v_ext of the adjacent external-stage cell (the x-face carries u, the y-face carries v), then `exchangeWithCorners`. Allocate and run none of this when no BC is external (`hasExternal_ = false`). |
| W2.3 | `swe/Momentum.cpp` | In `momentumSource()` and `interpolateVelocity()`, capture `const bool ext = hasExternal_;`. When `!ext`, the existing expressions must stay textually identical. When `ext`: `uA(jc,ic) = resX(jc,ic) ? uResX(jc,ic) : uu(jc,ic)` (likewise `vA`). **x self-advection:** if `u > 0 && resX(j,i-1)`, use `(u*u - uE*uE)/(2*dx)` with `uE = uResX(j,i-1)`; else if `u < 0 && resX(j,i+1)`, use `(uE*uE - u*u)/(2*dx)` with `uE = uResX(j,i+1)`; else the original upwind form. **Cross term** (advection of u by `vxc`): the original form with `uA` at `(j±1,i)`. **advY:** symmetric (self term with `resY`/`vResY` at `(j∓1,i)`; cross term with `vA` at `(j,i±1)`). **`interpolateVelocity`:** `uy`/`vx` from `uA`/`vA`. The CFL damping, viscous terms, drag and wind are unchanged. Document in the file header that u_ext = 0 reproduces the seed's 0.5 weight (entrance head v²/2g) and u_ext = u reproduces no entrance loss. |
| W2.4 | `SurfaceSolver.*`, `coupling/Coupler.cpp`, `driver/Simulation.*` | `void SurfaceSolver::routeStageDeposits(const Field2<real_t>& gain)`. For external-stage cells with `gain != 0`, set `eta -= gain` and book `gain·A` into `audit_.stageDeposit`. For zero_gradient stage cells, book `gain·A` only, without changing η. Call it in `Coupler::step` between `applySeepage(dt)` and `surface_.updateVelocity()`. Add a `mass_audit` column `stage_deposit` (coupled runs only, appended last), widen the `recordMassAudit` arrays 7 → 8, and add run-record closure line `stage_deposit_m3`. Update the coupled closure helper in `tests/unit/test_coupling.cpp` and the `b6_subcycled` surface identity in `run_regression.py` to include the new term (seed hunks). The b6_subcycled tolerance stays unchanged; record the new residual. |
| W2.5 | `transport/SurfaceTransport.cpp` | Re-apply the seed's audit-only booking of diffusive mass through global-edge faces into `surf_boundary` (declared allowed difference, §2.6) |
| W2.6 | `swe/FreeSurface.cpp` (`accumulateBoundaryFluxes`), `Simulation.cpp` | Entrance-head diagnostics, per eta BC. For each face into a stage cell whose step volume is **inflow** (with face area A_f > 0): v = \|vol\|/(dt·A_f), Q = \|vol\|/dt, ΔKE = ½·998·v²·Q·dt; accumulate `ke += ΔKE` and `rNum += ΔKE·(η_stage − η_interior)/(v²/(2g))`. At output flush, reduce over the BC sub-communicator and write run-record lines `stage_ke_influx_J` and `stage_entrance_head_ratio` (= rNum/ke, or `n/a`) per BC. At run end, for a zero_gradient BC with ratio < 0.5 and max v > 0.1 m/s, log the warning in §13.A. |
| W2.7 | docs | `theory/surface-water.md`: new section "Inflow closure at prescribed-stage boundaries" (definitions, formulas, reservoir faces, deposit routing, diagnostics, provenance = 2026-10-26 diagnosis `report_V1.md`, prototype measurements). Also `parameters.md`, `configuration.md`, `output.md`, `agents/case-setup.md` (recipe: still sea or lake → `external`; moving sea → `external` with velocity; channel cut → `zero_gradient`), and `feature-coverage.md` rows (inflow closure × surface-only / coupled / transport / wind). |

### Gates (authored first)

| Gate | Configuration | Reference | Asserts | Tolerance (`stage3-ob.yaml`) | Budget | Negative / anti-vacuity |
|---|---|---|---|---|---|---|
| **ob1** `ob-openbc/ob1-channel.yaml` | nx = 24, ny = 1, dx = dy = 1 m, flat bed, IC η = 1.0; west stage cell η = 1.005, east stage cell η = 1.000; Manning n = 0.01, viscosity 0; dt = 0.25 s (CFL ≤ 0.4), t_end 2400 s, strict solver. Lanes **A** `external`; **B** `external` with `velocity: {u: 0.1}`; **C** `zero_gradient`. | **Discrete steady chain** solved in `gate_openbc.py`. Unknown Q, marched face by face with the code's exact discrete relations: face depth = max(η_l, η_r) − max(b_l, b_r); effh = A_sx/V_sx with V_sx = ½(V_l + V_r) (`WetDry.cpp`); steady momentum u = (u − dt·adv)·D² − g·dt·effh·Δη·D with D = 1/(1 + ½dt·C_D·\|u\|·A_szx/V_sx) and C_D = g n²/h^{1/3}; conservative entrance (W2.3); root-find Q to hit the east stage. Self-check: a one-step update of a given state in Python equals the model's to 1e-14. Continuum Torricelli/Bernoulli is recorded. | A and B: max over faces \|u − u_chain\|/u_chain; steadiness: max\|Δu\| over the last 10 % of t; run-record ratio R against R_chain. C: u_C ≥ 2·u_A and R_C ≤ 0.2. | 1e-6; 1e-9 m/s; R within 0.05 | ≈ 2 s | C is itself the real negative (predicted u_C ≈ 1.5 vs u_A ≈ 0.31 m/s). A synthetic weight-1 trace (u = √(gΔη)) fails the A check. |
| **ob2** unit (`test_swe.cpp`) | `SweModule.ExternalClosureSelfAdvectionIsConservative`, `SweModule.ExternalClosureTangentialStateEntersCrossAdvection`, `SweModule.ZeroGradientStencilIsUnchanged`: 6×4 domains, one `momentumSource()` call on prescribed fields | Hand formulas of W2.3 | advX/advY at first interior faces; zero_gradient = v2.0 formula | 1e-15 relative | < 0.1 s | Masks swapped: fails |
| **ob3** `ob-openbc/ob3-tidal-basin.yaml` | 10×20, dx = dy = 10 m, flat bed 2 m deep, north stage row η = 0.2 sin(2πt/3600), n = 0.02, dt = 5 s, 2 cycles; lanes zero_gradient and external | Ledger identity | Per lane: surface ledger closure ≤ 1e-12 of cumulative \|stage exchange\|. Between lanes: cumulative stage inflow per cycle differs by ≤ 1e-3 relative. | as stated | ≈ 2 s | Anti-vacuity: max Fr² recorded and < 1e-3 (the regime of the bound) |
| **ob4** unit (`test_coupling.cpp`) | `CoupledModule.ExternalStageRemovesAndBooksCoupledDeposits` (seed test adapted), `CoupledModule.ZeroGradientStageBooksDepositsWithoutRemovingThem`: 4×4×6 coupled, north stage row, exfiltrating head beneath | Identities | external: after routing, η(stage) == η_bc exactly; stage_deposit == Σ gain·A; coupled closure ≤ 1e-10·V. zero_gradient: η(stage) == η_bc + gain; stage_deposit == Σ gain·A. | exact / 1e-10 | < 0.5 s | Routing call removed: η kick ≠ 0 (fails) |
| **ob5** `ob-openbc/ob5-gyre.yaml` | 20×30, dx = dy = 1 m, bed from 0.5 m to 3 m depth toward the north stage row (η = 0 const), n = 0.025, initial along-shore jet uu = 0.3 m/s within 5 m of the stage row, 1200 s; lanes external and zero_gradient | Energy diagnostics | external: domain KE monotone decreasing after 60 s; stage KE influx ≤ 1e-3 of the initial domain KE | as stated | ≤ 30 s | zero_gradient KE influx ≥ 10× external (anti-vacuity). If authoring cannot produce this contrast within 30 s, amendment: replace with the nightly `regression.s3.ob5.d1_record` (shortened D1 case from `impl_v1/d1_case_patched.yaml`, 4 ranks, 1500 s; external max\|u\| at `p_low_xc` ≤ 1e-4 m/s, zero_gradient ≥ 0.1 m/s) |
| **ob6** | `regression.s3.default_bitwise` | T0 hashes | bitwise with the declared exceptions | — | (shared) | Any external lane would differ |
| **ob7** unit (`test_config.cpp`) | `ConfigTest.InflowClosure*` | — | §13.A messages; resolved-config round trip | — | < 0.1 s | — |
| **ob8** nightly record (`stage3_nightly`) | b6 ss/td with `inflow: {closure: external}` on every eta BC | Kuan experiment | Interface MAE against the experiment ≤ the existing b6 bound. Legacy-golden residuals recorded, not gated (the golden contains the deposit kick). | existing `b6-kuan.yaml` experiment bounds | ≈ 15–20 min each | — |

`unit.ob_gates.negative`: chain solver spot values (two independent evaluations: shooting and a direct Newton on the stacked system); each check fails on corrupted traces.

**DoD T2:**
- [ ] ob1–ob7 green, ob8 filed;
- [ ] `default_bitwise` green (exceptions: `stage_deposit` column, `surf_boundary`);
- [ ] b6_subcycled still green, new identity residual recorded;
- [ ] docs and coverage rows;
- [ ] Tier N green;
- [ ] snapshot `T2`.

---

## 7. T3 — Feature B: periodic lateral topology

**Semantics.** `domain.periodic.x: true` makes the faces between global columns nx−1 and 0 interior faces for **every** module: surface momentum and continuity, groundwater matrix and fluxes, transport and baroclinic terms. There is no x boundary any more. The same applies to `y`.

### Work items

| W | File(s) | Exact change |
|---|---|---|
| W3.1 | `Config.hpp/.cpp`, `ConfigSchema.cpp`, `RunRecord.cpp` | `DomainConfig`: `bool periodicX = false, periodicY = false;`. Schema: `{"periodic", map({{"x", boolean(false)}, {"y", boolean(false)}})}` in `domain`. Cross-field: `domain.periodic.x: needs nx >= 3 (got N)` (same for y). Echo in the resolved config; run-record line `topology: periodic_x/periodic_y`. |
| W3.2 | `core/Grid.hpp/.cpp` | Store the flags and wrap the neighbours at `Grid.cpp:67-70`: `rankWest_ = (pi_ > 0) ? rank_ - 1 : (perX_ ? rank_ + (px_ - 1) : MPI_PROC_NULL)`; `rankEast_ = (pi_ + 1 < px_) ? rank_ + 1 : (perX_ ? rank_ - (px_ - 1) : MPI_PROC_NULL)`; `rankSouth_ = (pj_ > 0) ? rank_ - px_ : (perY_ ? rank_ + px_ * (py_ - 1) : MPI_PROC_NULL)`; `rankNorth_ = (pj_ + 1 < py_) ? rank_ + px_ : (perY_ ? rank_ - px_ * (py_ - 1) : MPI_PROC_NULL)`. Add `struct EdgeTopology { int nx, ny; bool perX, perY; }` with `KOKKOS_INLINE_FUNCTION` predicates `west(gi)` (`!perX && gi == 0`), `east(gi)` (`!perX && gi == nx - 1`), `south(gj)`, `north(gj)`, `outsideX(gi)` (`!perX && (gi < 0 \|\| gi >= nx)`), `outsideY(gj)`, plus `Grid::edgeTopology()`, `periodicX()` and `periodicY()`. **No change** to `HaloExchanger` or the gid exchange. Verified by reading: receives use tag `200 + dir` from `neighbor_[dir]`, and sends use `200 + (dir ^ 1)`. With `px = 1` the west send to self lands in the east halo; with `px = 2`, the two messages to the same rank are told apart by tag. `Grid.cpp:100-103` Sendrecv tags 100/101 behave the same. Document this in `developer-guide/halo-exchange.md`. |
| W3.3 | the site table below | Convert every global-index edge test to an `EdgeTopology` predicate. Capture `const EdgeTopology topo = grid_.edgeTopology();` in each kernel. Rank-based sites need no edit, but the gates must prove each one is halo-covered. |
| W3.4 | `bc/BoundarySet.cpp` (`appendSideFaces`, ctor) | Never attach a face normal to a periodic axis. Count per BC the faces skipped for that reason. If a BC ends with zero members and skipped > 0, fatal `boundary condition '<name>': its polygon selects only faces on the periodic x edges, which domain.periodic.x makes interior; remove the condition`. Add the per-side face counts and the skipped count to the run-record BC summary. This is the generic rule "attach edge faces only on axes that have a boundary"; the v2.1 degenerate-axis fix (R1) later becomes one more condition here. |
| W3.5 | `Simulation.cpp` (after inputs are read) | Seam diagnostic for `bottom_elevation`, the soil map, the Manning map and 2-D/3-D IC rasters on each periodic axis. Compute s = max\|Δ\| across the seam and m = the largest interior neighbour jump. If s > m, warn `domain.periodic.x: bottom_elevation jumps by <s> m across the x seam (largest interior jump <m> m); the seam is a discontinuity, not a periodic continuation` and add it to the run-record `topology.warnings`. |
| W3.6 | `gw/TerrainMetric.cpp:342,360` | The boundary-face slope skip becomes `topo.outsideX(gi) \|\| topo.outsideX(gi + 1)`. On a periodic axis the seam face slope is computed from the wrapped (halo-exchanged, `:152`, `:306`) neighbour. |
| W3.7 | `tests/regression/run_regression.py` (`run_dihedral_battery`) | For the transforms `transpose`, `antitranspose`, `rot90` and `rot270`, swap `domain.periodic.x` and `.y` in the variant YAML |
| W3.8 | `scripts/check_forbidden.sh` | New rule "raw global-edge test (use EdgeTopology)", via `scan_src_excluding` with regex `\b(gi\|gj\|giMinus\|gjMinus\|iGlob\|jGlob)[[:space:]]*(==\|<\|>=)[[:space:]]*(0\|nxGlobal - 1\|nyGlobal - 1\|nxGlobal\|nyGlobal)\b` and exemption `core/Grid\.hpp`. Run it once on the converted tree; it must report zero hits. |
| W3.9 | docs | `developer-guide/grid-indexing.md` (EdgeTopology, wrap formulas); `halo-exchange.md`; `configuration.md`, `parameters.md`; theory pages (a seam is an interior face; no closed-edge leg fold; no stage reconstruction on a periodic axis); `symmetry-exemptions.md` (a periodic axis has no minus/plus edge class); `agents/case-setup.md` (uniform-coast recipe: drop the x walls and the `*_xe` soil types, generate random fields periodically, e.g. FFT synthesis on the torus); `feature-coverage.md` rows (periodic × surface / groundwater / coupled / transport / heat / wind / eta-external). |

**Site table** (T0 tree; action codes: **AUTO** = rank-based, follows the wrap, must be halo-covered; **TOPO** = convert to an `EdgeTopology` predicate; **NONE** = unaffected):

| File:line | Code | Action |
|---|---|---|
| `core/Grid.cpp:67-70` | neighbour assignment | W3.2 wrap |
| `core/Grid.cpp:270-273, 308-311`; `core/HaloExchanger.cpp:36, 265` | gid / halo unpack guards | AUTO |
| `bc/BoundarySet.cpp:39-50` | `iGlob == 0` / `== nx-1`, `jGlob …` (side faces) | TOPO + W3.4 |
| `swe/SurfaceSolver.cpp:210-213` | bed edge ghosts | AUTO (exchanged at `:209`) |
| `swe/FreeSurface.cpp:64-67` | η edge ghosts | AUTO |
| `swe/FreeSurface.cpp:263-272` | closed-edge leg fold `sct -= sxm/sxp/sym/syp` | TOPO (on a periodic axis the leg stays live; its COO column becomes valid through the wrapped gid) |
| `swe/FreeSurface.cpp:469-478` | boundary-outflow audit at global edges | TOPO |
| `swe/WetDry.cpp:74-75, 90-91` | `oneColumn` / `oneRow` | NONE (nx, ny ≥ 3 enforced on periodic axes) |
| `swe/WetDry.cpp:135-158` | west/south geometry ghosts | AUTO |
| `swe/WetDry.cpp:263-266, 413-414` | velocity edge ghosts, ghost corners | AUTO |
| `swe/WetDry.cpp:385-400` | stage-cell boundary-face reconstruction | TOPO (no reconstruction on a periodic axis) |
| `transport/SurfaceTransport.cpp:133-136` | `westCell…northCell` (superbee guards `:155-187`, boundary booking, W2.5 audit) | TOPO |
| `transport/SurfaceTransport.cpp:463-466` | scalar edge ghosts | AUTO |
| `transport/SubsurfaceTransport.cpp:324-327` | edge-cell flags | TOPO |
| `transport/SubsurfaceTransport.cpp:766-769` | ghosts / side BC | AUTO |
| `gw/RichardsSolver.cpp:28-31` | soil ghost mirror | AUTO (soil halos exchanged at `:263`) |
| `gw/RichardsSolver.cpp:887-896` | boundary-inflow audit | TOPO |
| `gw/Predictor.cpp:190-199` | edge K faces by side code | TOPO |
| `gw/Predictor.cpp:280, 293` | masked-column seal `gi >= 0 && gi < nxGlobal && gid < 0` | TOPO (`!topo.outsideX(gi) && gid < 0`) |
| `gw/Predictor.cpp:390-399, 415-424, 464-489` | head-edge leg doubling, side-flux RHS, edge RHS fold | TOPO |
| `gw/Predictor.cpp:592-595` | head ghosts | AUTO |
| `gw/Corrector.cpp:126-127, 155, 203-204, 227` | edge face-flux branches | TOPO |
| `gw/Corrector.cpp:509-512` | moisture ghosts | AUTO |
| `gw/Baroclinic.cpp:99, 110, 134, 144` | edge face density ratios | TOPO |
| `gw/TerrainMetric.cpp:22-60` | host ghost fills | AUTO |
| `gw/TerrainMetric.cpp:342, 360` | boundary-face slope skip | W3.6 |
| `gw/Reallocate.cpp` | no global-edge test (`i0`/`nxGlobal` are unused) | NONE |

### Gates (authored first)

| Gate | Configuration | Reference / assertion | Tolerance (`stage3-pd.yaml`) | Budget | Negative / anti-vacuity |
|---|---|---|---|---|---|
| **pd1** `pd-periodic/pd1-alongshore.yaml` (+ transposed `pd1-alongshore-T.yaml`) | 8×4, dx = dy = 10 m, flat bed, η = 0.5 m, Manning n = 0.05, `thin_layer_depth` 0.1, wind Garratt cap 3.5e-3 with u10 = 10, v10 = 0; periodic x; y edges closed with `kind: velocity value 0` (g9 practice); dt = 2 s, t_end = 4400 s (≈ 40 t_r); strict solver; monitor `uu` at (4, 2). Transposed: 4×8, periodic y, v10 = 10, x walls. | **Discrete closed form** (from `Momentum.cpp`, steady and uniform, with Δη = 0): g n² u² (1 + a/2)/h^{4/3} = ρ_a C_d (W − u)²/(ρ_w h), a = ½·dt·g n² u/h^{4/3}, with C_d recomputed by `gate_wind.cd_garratt(10)`; solve by bisection. Predicted u = 0.07446 m/s (bisection, 2026-10-27); t_r = h^{4/3}/(2 g n² u) ≈ 109 s. The continuum √(τh^{1/3}/(ρ_w g n²)) is recorded. Asserts: η spread over the grid; `uu` spread; \|u − u_disc\|/u_disc; steadiness (monitor span over the last 540 s). Default-solver lane: η spread bound = measured × 10 (bjacobi-ICC is not translation-symmetric; the T10c diagnosis measured 1e-6 class). | η spread ≤ 1e-12 m; u spread ≤ 1e-13 relative; u ≤ 1e-10 relative; span ≤ 1e-13 m/s | ≈ 2 s | Real negative lane: periodic off with x walls gives u ≈ 0 plus setup (must fail the u check) |
| **pd2** roll / seam invariance | **surface** `pd2-surface.yaml`: 12×6, dx = dy = 5 m, bed bump b(x) = 0.2·exp(−((x−20)/8)²) (raster), η = 1.0, n = 0.03, wind u10 = 8, v10 = 2, tracer s = 10 in columns i = 10–11, y walls, periodic x, 600 s at dt = 1 s, superbee and upwind sub-lanes. **gw** `pd2-gw.yaml`: 12×6×8, dx = dy = 1, dz = 0.25, K(x,y) = 1e-5·(1 + 0.5 sin(2πx/12)·cos(πy/6)) as harness-generated soil types, θs .41, θr .05, α 2, n 1.5, y− head η = 1.5, y+ head η = 1.0, top flux −1e-7, IC water table 1.25, Ss 1e-5, reallocation off, 1 day. **coupled** `pd2-coupled.yaml`: 12×8×6, dx = dy = 2, bed sloping in y plus an x bump, north tide stage (closure external, 0.2 m, 3600 s), salinity with consistent updates, reallocation off, 2 h at dt = 2 s. | Each lane runs as is and with **every x-dependent input rolled by k ∈ {3, 6}** columns. Assert roll(output_k, −k) == output_0 for all fields (η, u, v, s, h, θ, c). | 1e-11 strict initially; adjudicated to ≥ 10× headroom over the measured value | ≈ 20 s | The gate-first red run (before W3.3) must fail; a synthetic seam glitch fails the checker; anti-vacuity: the tracer crosses the seam (max s at the seam column > 1 during the run) |
| **pd3** rank invariance | pd2 bases at `mpi_nx` = 1 (self-wrap), 2 (same rank on both sides), 4; plus 2×2 for the coupled lane | Strict rank-invariance metric as in `gate_rank_invariance_tracer` | 1e-12 local (1e-9 on CI per V2-A9) | ≈ 15 s | r1 `--ranks` check proves the world size (`FREHG_EXPECT_RANKS`) |
| **pd4** seam is not a boundary | `pd4-tidal-periodic.yaml` (ob3 basin, periodic x, no x walls, closure external); unit `SweModule.PeriodicSeamCarriesNoBoundaryFlux` (periodic channel, uniform IC uu = 0.3, y walls, 10 steps) | Ledger closure ≤ 1e-12. Unit: `audit.boundaryOutflow == 0.0` exactly, and volume constant to 1e-14. | as stated | ≈ 2 s | The unit test on the pre-W3.3 build shows a nonzero outflow (red evidence) |
| **pd5** groundwater analytic | `pd5-confined.yaml`: 6×10×5, dx = dy = 1, dz = 0.5, surface at 2.5 m, fully saturated (y− head η = 3.0, y+ head η = 2.8, both above the surface), Ss 1e-5, 3600 s, periodic x. Lanes: uniform K, and K varying by column (3 periodic values). | Hydraulic head linear in y; x-face fluxes zero; q_y per column = K·0.2/L | \|H − H_lin\| ≤ 1e-10 m; \|q_x\| ≤ 1e-14·\|q_y\| | ≈ 2 s | Synthetic |
| **pd6** setup and schema negatives | `regression.s3.pd6` runs `pd6-x-wall.yaml` (periodic x plus an x-edge `velocity` BC) and asserts nonzero exit and the W3.4 message. `ConfigTest.Periodic*` cover nx < 3 and the round trip. The pd2 surface log must contain the W3.5 seam warning (its bump has a 0.002 m seam jump). | message substrings | — | ≈ 1 s | — |
| **pd7** | `default_bitwise` | — | bitwise | (shared) | — |
| **pd8** orientation battery | `pd-periodic/pd8-orient-base.yaml`: 8×8, periodic x, y walls, wind u10 = 6, v10 = 3, bed bump, tracer; 8 transforms (W3.7 flag swap), strict | Existing battery machinery | existing `x-orient.yaml` classes | ≈ 6 s | Flag swap disabled: the rot90 variant fails |
| **restart** | `s3-restart` lane `pd2-surface` (periodic) | bitwise | exact | ≈ 2 s | — |

`unit.pd_gates.negative`: discrete closed form spot values (bisection against a fixed-point iteration), roll checker, ledger checker. MPI unit tests: `mpi.halo.n{1..4}` gain `HaloTest.PeriodicWrapFillsHalosFromTheFarEdge`, which asserts on a periodic grid that each halo column equals the far-edge interior column, for px = 1, 2, 4.

**DoD T3:**
- [ ] pd1–pd8 green;
- [ ] W3.8 rule clean;
- [ ] the site table is fully converted, with every TOPO site cited in `report-T3.md`;
- [ ] `default_bitwise` green;
- [ ] docs;
- [ ] Tier N green;
- [ ] snapshot `T3`.

---

## 8. T4 — Diagnostics and output

| W | File(s) | Exact change |
|---|---|---|
| W4.1 | new `io/IntervalRates.hpp/.cpp`, `Simulation.*` | A registry of named per-cell accumulators (Field2/Field3 volume, m³, or depth, m) with `add(name, view)` contributor hooks. On each output write: write `acc/(t − t_prev)` (per cell area for 2-D rates, giving m/s; per cell volume for 3-D, giving 1/s), then reset. Allocate only the requested variables. Checkpoint `acc_<name>` (Field) and scalar `acc_t_prev`. |
| W4.2 | `coupling/Exchange.cpp` | Output `seepage_rate` (surface, coupled only): interval mean of the coupler's per-step `gain` (applied seepage + bounce-back + vent deposits) divided by the cell area, in m/s, positive upward. By construction, Σ rate·A·Δt equals the Δ of `mass_audit.seepage`. This replaces snapshot `seepage` for analysis; `seepage` stays, documented as a per-step sample (`report_V3.md`). |
| W4.3 | `Config.*`, `ConfigSchema.cpp`, `Simulation.cpp`, `RunRecord.cpp` | `output.monitor_interval` (integer seconds ≥ 0; 0 or absent means every step). Separate audit **accumulation** (every step, unchanged) from row **writing** (decimated). Write rows at the first step with t ≥ the next multiple of the interval, at checkpoint times and at t_end. Point monitors sample at the same instants. On restart, the next row time is recomputed from t0 like the output schedule. |
| W4.4 | docs | `output.md` (rate semantics, decimation), `parameters.md`, `adding-an-output-variable.md` (the interval-rate pattern), `agents/postprocessing.md` (use `seepage_rate`, never snapshot `seepage`, for exchange metrics) |

| Gate | Configuration | Asserts | Tolerance (`stage3-dg.yaml`) | Budget | Negative / anti-vacuity |
|---|---|---|---|---|---|
| **dg1** `dg-diagnostics/dg1-seepage.yaml` | Coupled 4×4×8, subcycled dt = 1 s, dtg ≤ 10 s; one side head above the surface (exfiltration) and a ponded corner (infiltration); reallocation off; 1 h, outputs every 300 s, `seepage` and `seepage_rate` | Per interval: \|Σ rate·A·Δt − Δ(mass_audit.seepage)\| ≤ 1e-12·max(\|Δ\|, 1e-12 m³) | 1e-12 | ≈ 2 s | **Anti-vacuity:** in ≥ 1 output frame the snapshot `seepage` is identically 0 on every wet cell while `seepage_rate` ≠ 0 (reproduces V3) |
| **dg2** decimation | dg1 run twice: every step vs `monitor_interval: 60` | Every decimated row (time and all columns of every point monitor and audit table) equals the every-step row at the same time, bitwise; the last row is at t_end | exact | ≈ 3 s | Synthetic row shift fails |
| **restart** | `s3-restart` lane dg1 (checkpoint at 1800 s) | `seepage_rate` outputs bitwise | exact | ≈ 2 s | — |

`unit.dg_gates.negative`: checker sensitivity.

**DoD T4:**
- [ ] dg1, dg2 and the restart lane green;
- [ ] `default_bitwise` green (`monitor_interval` absent means every step);
- [ ] docs;
- [ ] Tier N green;
- [ ] snapshot `T4`.

---

## 9. T5 — Feature C core: coupled evaporation with salt conservation

**Semantics (no vegetation yet).**
- An `evapotranspiration` (ET) block replaces both legacy evaporation keys.
- **Open-water channel:** every cell with depth > 0 at the surface source step loses min(E_w·dt, depth), where E_w is the open-water rate of the selected potential mode.
- **Soil channel:** every exposed column loses the soil rate E_s through the top face, weighted by the **dry-time fraction** of each groundwater substep and capped by the top cell's supply.
- **Complementarity holds on a time basis by construction.** Per cell, the time covered by the surface channel plus the time covered by the soil channel never exceeds elapsed time.
- **Salt:** evaporation removes water only.
  - Surface: a drying cell's remaining salt moves into a per-cell **salt store**, released into the water at the first wet step (O4).
  - Subsurface: the top cell keeps its salt (exact concentration factor).
  - Seepage held over a dry cell takes its salt into the store with it (D9).

**Inventory of defects this phase removes** on the ET path (numbering from the 2026-10-16 plan; the legacy paths are untouched):
- D1: predictor/corrector guard mismatch.
- D2: audit sums the whole top flux.
- D3: condensation and strip threshold.
- D4: no supply cap.
- D5: post-step-depth transport flag.
- D6: dry-out salt deletion.
- D8: rain/evap either/or in transport.
- D9: held-seepage salt.
- Window-lag double counting (assessment §1.C caveat 4).

D7 (solubility) remains out of scope.

### Work items

| W | File(s) | Exact change |
|---|---|---|
| W5.1 | `Config.hpp/.cpp`, `ConfigSchema.cpp`, `RunRecord.cpp` | New `struct EvapotranspirationConfig { bool enabled; enum class Potential { Bulk, PenmanMonteith, Prescribed } potential; SeriesOrConstant prescribed; bool openWater = true; bool soil = false; std::vector<std::array<real_t,2>> soilPolygon; enum class Limiter { AlphaBarton, BetaLeePielke, None } limiter; /* T7 */ EnergyConfig energy; };` in `FrehgConfig`. `AtmosphereConfig` gains `enum class AeroResistance { Liu2006, Fao56, LogProfile }` with `refHeight` (2.0), `z0m` (0.01), `z0h` (0.001), `displacement` (0.0). `SoilType` gains optional `thetaFc` (absent → θ(h = −3.3 m) from van Genuchten). `TransportConfig` gains `enum class SurfaceDryout { Delete, Store } surfaceDryout = Delete;` (the T8 DF constant). Schema and cross-field rules: §13.C. Resolved-config echo. PM keys are parsed in T7; until then `potential: penman_monteith` fails validation with `evapotranspiration.potential penman_monteith: not available before phase T7` (removed in T7). |
| W5.2 | new `src/atm/EtFormulas.hpp` (host/device, pure, doxygen with provenance) | `aeroResistance(law, U, params)`: liu2006 = `94.909·U^-0.9036` (the existing formula, moved; `BulkAerodynamic.hpp` keeps a forwarding call so its consumers stay bitwise); fao56 = `208/U2` with U2 = U·4.87/ln(67.8 z_ref − 5.42) (FAO-56 eq. 47); log_profile = `ln((z−d)/z0m)·ln((z−d)/z0h)/(0.41²·U)`. `potentialBulk(T_s, P, U, q_a, r_a, r_s, α)` = ρ_a(T_s,P)·(α·q_sat(T_s) − q_a)/(ρ_w·(r_a + r_s)); at r_a = Liu, r_s = 0 it is **the same expression tree** as `atm::evaporationRate` (asserted bitwise in et0). `limiterBeta(θ, θ_fc)` = θ ≥ θ_fc ? 1 : 0.25(1 − cos(πθ/θ_fc))² (Lee & Pielke 1992). `supplyCap(E, θ, θr, θs, dz, dtg)`: evaporation min(E, max(0, θ − θr)·dz/dtg); condensation max(E, −max(0, 0.9999 θs − θ)·dz/dtg). |
| W5.3 | `swe/SurfaceSources.cpp`, `SurfaceSolver.hpp/.cpp` | ET branch of `evapRain()` (active iff `et.enabled && et.openWater`). Per step E_w comes from the potential mode: bulk with r_s = 0, α = 1; prescribed = series. Per cell, after rain is added: `removed = (depth > 0) ? (E_w > 0 ? fmin(E_w·dt, depth) : E_w·dt) : 0` (negative = condensation onto water), `eta -= removed`; `evapRemoved_(j,i) = removed` [m]; `rainAdded_(j,i) = rain·dt·mask`; `wetFrac_(j,i) = (depth > 0) ? ((E_w > 0) ? fmin(1, depth/(E_w·dt)) : 1) : 0`; `wetClock_(j,i) += wetFrac·dt`; also `wetClockStepStart_ = wetClock` before the increment. The existing clamp follows unchanged; its volume stays in `clamped`. Book Σ removed·A into `audit_.evapVolume` (existing `evaporation` column). `currentEvaporation()` keeps returning the uniform E_w (legacy transport path only). |
| W5.4 | `gw/RichardsSolver.hpp`, new `gw/Evapotranspiration.cpp` (add to `src/CMakeLists.txt` `frehg_gw` list explicitly) | `setupEvapotranspiration(config)`: soil mask (region or all), `etKtop_`, θ_fc per cell, `qVap_` (Field2, m/s), `wetClockLast_` (Field2, s), `soilTime_` (Field2, s), audit fields. Setup fatal if the region overlaps a `groundwater_top` BC (the `RichardsSolver.cpp:539-544` message). `updateEvapotranspiration(t, dtg)` is called in `step()` **after** `classifyCoupledTop(dtg)` (uncoupled: after `updateBoundaryValues`). Per substep, per soil column: (1) E_s = potential(r_s = 0) × limiter: `alpha_barton` uses α₁(w_g) inside the humidity difference with w_g the existing surface extrapolation of `updateBulkEvaporation` (copy the lines verbatim; et1 asserts bitwise); `beta_lee_pielke` multiplies by β(w_g, θ_fc); `none` uses factor 1. (2) **Dry fraction**: uncoupled f_dry = 1; coupled, with a = tSub − dtg and b = tSub, `W_b = wetClockStepStart + wetFrac·(b − tStepStart)`, `f_wet = clamp((W_b − wetClockLast)/dtg, 0, 1)`, `wetClockLast = W_b`, f_dry = 1 − f_wet. (3) Applied flux `qVap = (mode == 0 \|\| !coupled) ? supplyCap(f_dry·E_s, θⁿ_top, …) : 0`. (4) `soilTime += (qVap != 0) ? f_dry·dtg : 0`. (5) Audit: demand D = f_dry·E_s·A·dtg, applied = qVap·A·dtg, unrealised = D − applied. `tStepStart = t − dt` is passed through a new `GwCoupling::stepStart` set by `Coupler::step` before the substep loop; the views `wetClockStepStart`, `wetFrac` are attached through `GwCoupling`. |
| W5.5 | `gw/Predictor.cpp`, `gw/Corrector.cpp` | Add `Field2<int> etCol = etColumn_` (1 on ET soil columns; empty view when ET is off, guarded by `const bool et = etActive_`). **Predictor** top branch: if `et && etCol(j,i)`, use `qtop = qVap(j,i)` with **no guard**, in place of the topValue guard pair (`Predictor.cpp:532-537`), for coupled mode 0 and uncoupled. **Corrector** `computeFluxes` mode-0 branch (`:330-336`) and the uncoupled flux branch: if `et && etCol`, `q += qVap·az` (no guard). **Corrector** `applyCoupledTopBookkeeping` (`:446-450`): if `et && etCol`, replace the strip with `seepAccum -= qVap·dtg; sumV += qVap·dtg·az;` for any sign and with no depth test (qVap ≠ 0 only in mode 0). Legacy paths are byte-identical when `!et`. |
| W5.6 | `RichardsSolver.cpp` (`accumulateTransportWindow`, `enableTransportWindow`, `resetTransportWindow`) | `vapWin_` (Field2, m³) `+= qVap·az·dtg` per substep. Exposed as `windowVapourVolume()` |
| W5.7 | `transport/SurfaceTransport.cpp`, `ScalarSolver.hpp/.cpp`, `Simulation.cpp` | **Salt store** (salinity only, `surfaceDryout == Store`). Field2 `sStore_` (scalar·m³), halo-free (cell-local). In the update pass: (a) non-wet branch (`!wet`, i.e. `vflux ≤ 0 \|\| dept ≤ 0`): replace `sumAdjust -= sm` with `store += sm; sumStoreIn += sm`. (b) Wet branch, **after** the limiter, bounds and dilution: if `store > 0`: `value += store/vBasis; sumStoreOut += store; store = 0`. Then the anchor as today. (c) Per-cell volumes under ET: `vre = (rainAdded − evapRemoved)·A` replaces the either/or rule (D8) and honours masks (`rainAdded` already includes the rain mask). New transport audit **instantaneous** column `surf_store` (Σ store, appended last). Identity: Δ(surf_mass + surf_store) = exchange + surf_source + surf_boundary + surf_adjust + surf_anchor. Output variable `salt_store` (transport, 2-D). Checkpoint `s_store`. Without `store`: arithmetic unchanged (bitwise). |
| W5.8 | `SubsurfaceTransport.cpp` (update pass, `:661-690` block) | ET top cell (`et && etCol && k == kTop`, needs `consistentSub`): the exact factor uses the window vapour volume, `f = (vgflux + vapWin)/vgflux` for vapWin > 0 (evaporation; for condensation f = vgflux/(vgflux − vapWin) applied to `lo`). The window anchors on `sOld` as the cauchy path does. This replaces the `coupledEvapTop` post-step-depth flag on ET columns (D5). The legacy `coupledEvapTop` branch is kept for non-ET columns. |
| W5.9 | `coupling/Exchange.cpp`, `Coupler.*`, `SurfaceTransport.cpp` | Held-seepage salt (ET with `store` only). In `beginWindow` record `heldPrev_ = seepAccum` (positive part). In `applySeepage`, set per cell `exfNew_ = max(0, a − heldPrev)` and `released_ = (applied && heldPrev > 0) ? heldPrev : 0`. Surface transport exchange (consistent exchange path), exfiltration part: mass = `exfNew·A·sKp` debited from the subsurface top cell every step (as the water left it). If applied, it is credited to the surface cell's mass. If held, it is credited to `sStore_` (booked as exchange plus store-in). The released held volume carries **no** exchange salt (its salt is already in the store, which W5.7(b) releases because the cell is now wet). The diffusive interface term and infiltration are unchanged. |
| W5.10 | `Simulation.cpp`, `RunRecord.cpp` | Audits: `mass_audit` keeps `evaporation` (= open-water removal). `gw_mass_audit` appends `soil_evaporation` and `soil_evaporation_unrealised` (ET only). Closure doc: Δvolume = boundary_in − ss_storage + realloc − vloss, with `soil_evaporation` a sub-account of `boundary_in`. New output section `output.variables.atmosphere` (2-D, group `/atmosphere`): `evaporation_rate` (interval-mean open-water removal, m/s; T4 mechanism), `soil_evaporation_rate`, `soil_evaporation_unrealised_rate`, `wet_time`, `soil_evaporation_time` (cumulative s). Checkpoint: `et_wet_clock`, `et_wet_clock_last`, `et_soil_time`, `s_store`. Run-record closure lines `open_water_evaporation_m3`, `soil_evaporation_m3`, `soil_evaporation_unrealised_m3`, `salt_store_final`. |
| W5.11 | docs | New `theory/evapotranspiration.md` (mkdocs nav under Theory): channels, time-basis complementarity with the W5.4 formula, supply cap, limiters, r_a laws, salt store, held-seepage salt, audit identities, deviations from v2.0 (D1–D9 list). Also `parameters.md`, `configuration.md`, `output.md`, `atmosphere.md` (r_a selector), `transport.md` (store), `agents/case-setup.md` (evaporating tidal-flat recipe), and `feature-coverage.md` rows (ET × coupled / uncoupled / surface-only / transport / periodic / wind / external stage / restart / rank). |

### Gates (authored first) — `benchmarks/et-evapotranspiration/`, `tolerances/stage3-et.yaml`

Shared met (`et-met`): T_air 25 °C, T_s 30 °C, P 101.3 kPa, RH 0.3, U 8 m/s, `wind_speed_floor` 0.5. The gate recomputes E_w in Python from `gate_et.potential_bulk` with Liu r_a. Predicted value ≈ 1.3e-6 m/s; the README records the computed value. Mild met (`et-met-mild`, for et1/et2 so that the supply cap never binds): T_air 20 °C, T_s 20 °C, RH 0.6, U 2 m/s.

| Gate | Configuration | Reference / asserts | Tolerance | Budget | Negative / anti-vacuity |
|---|---|---|---|---|---|
| **et0** (part 1) unit `test_atm.cpp` | `AtmEt.BulkPotentialMatchesV20EvaporationRateBitwise` (r_s = 0, α, Liu); `AtmEt.AeroResistanceLaws` (Liu at U = 1, 2, 8; FAO 208/u2 with the z_ref conversion; log profile against a hand value); `AtmEt.LeePielkeBeta` (0, ½θ_fc, θ_fc, θ_fc+); `AtmEt.SupplyCap` (both signs); `AtmEt.ThetaFcDefaultIsRetentionAt330cm` | Hand / Python values recorded in the test comments | 1e-15 relative (bitwise where stated) | < 0.1 s | Perturbed coefficients fail (`unit.et_gates.negative` mirrors each formula in Python) |
| **et1** new block ≡ old path | 1×1×30 column, sand (rl1 sand), dz 0.02, uncoupled, top bulk E with `et-met-mild`, salinity off, 1 day, dtg 60 s fixed, reallocation off. Lane **a**: unsaturated start (h = −0.5 m uniform); `evapotranspiration{bulk, soil, alpha_barton}` vs `groundwater.evaporation{mode: bulk}`. Lane **b**: saturated start. | a: all outputs bitwise. b: the old path shows the D1 predictor skip (first-step θ differs); the new path's ledger closes exactly, with soil_evaporation = Σ qVap·A·dtg. | bitwise / 1e-14 | ≈ 3 s | Anti-vacuity: lane a has unrealised ≡ 0 (cap never binds; else bitwise is not expected and the case must be milder); lane b differs from the old path at step 1 (else the D1 case is not reached; redesign) |
| **et2** coupled dry ≡ uncoupled | et1a column as the subsurface of a coupled 1×1 run with bed 0.5 m above the water table, no rain, `coupling: sync` with dt = 60 s | E(t), θ(z,t), h(z,t) against et1a new-block | 1e-10 relative | ≈ 3 s | Anti-vacuity: surface depth ≡ 0 at every output and f_dry ≡ 1 |
| **et3** switch | 1×1 coupled, bed 0 m, 1×1×20 sand column saturated with the water table at the surface, impermeable sides/bottom, Ss 1e-5, strict solver, pond h₀ = 0.4 mm, `et-met`, subcycled dt = 1 s, dtg fixed 7 s, 910 s, outputs every 70 s (= 10·dtg, so each output coincides with a substep end) with `atmosphere: [evaporation_rate, soil_evaporation_rate, wet_time, soil_evaporation_time]`. Lane **prescribed** with E = 1.3e-6 constant. | (1) Pond phase (t < t_d = h₀/E_w ≈ 308 s): cumulative open-water evaporation = E_w·t to 1e-12 relative, and depth(t) + cumulative exchange/A = h₀ − E_w·t to 1e-12·h₀ (the exchange term is the audited seepage; it is ≈ 0 but not asserted zero). At dry-out: h₀ = open-water removed + \|clamped\| + exchange/A to 1e-12·h₀. (2) Per output, per cell: wet_time + soil_evaporation_time ≤ t + 1e-9 (**no double count**). (3) Cumulative open-water + soil evaporation = E_w·t to 1e-10 relative at every output while α₁ = 1 (θ_s = 0.41 ≥ 0.375; asserted from the θ output). (4) Continuity: the soil rate in the first full frame after t_d over E_w is within 1e-3 of 1. (5) Ledger closes. | as stated | ≈ 3 s | **Real negative:** the v2.0 lane (surface `evaporation: {mode: bulk}` + `groundwater_top` `flux` = E_w) must violate (3) by > 1e-3 (window-lag double count). Anti-vacuity: ≥ 1 frame with 0 < wet fraction < 1. |
| **et4** water balance (E lane; T lane in T6) | Closed 1×1×30 silt column, uncoupled, ET soil, 3 days, reallocation off | ΔS = −∫E; vloss ≡ 0; unrealised booked = demand − applied | 1e-12 relative | ≈ 3 s | Drop the soil column from the identity: fails |
| **et7** supply cap (E part) | 1×1×10 column at θ = θr + 1e-4 in the top cell, very high E (prescribed 1e-5 m/s), 1 h | min θ ≥ θr − 1e-14; unrealised = ∫(f_dry·E − qVap) exactly; vloss = 0 | exact / 1e-14 | ≈ 1 s | Synthetic trace with θ < θr, or with unrealised ≠ demand − applied, fails the checker (`unit.et_gates.negative`) |
| **et8** salt (E lanes) | **a** soil: closed column, salinity 10, density off, ET soil, K tiny (1e-12) so redistribution is negligible. **b** pond: 1×1 surface-only + transport, 5 mm at 35 psu, open water to dryness, then rain 2 mm at t = 7200 s, `store`. **c** = b with `delete`. **d** held seepage: 2×1 coupled, cell 1 bed below the side head (seepage face, dry surface), `store`. | a: total salt constant 1e-12; per cell c·θ·V = c⁰θ⁰V; subs_adjust and subs_anchor ≤ 1e-12·M. b: s(t)·h(t) = 35·h₀ to 1e-12 until dry; at dry-out surf_mass = 0 and surf_store = 35·h₀·A; first wet step s = 35·h₀/h_rain; Σ salt constant 1e-12. d: total salt constant 1e-12; top-cell subs_adjust ≤ 1e-12·M; store > 0 while held. | as stated | ≈ 5 s | c: loss = 35·h₀·A shows in surf_adjust (the gate asserts that c **fails** the conservation check — the real negative). d with `delete`: subs_adjust > 0. |
| **et9** coupled tidal smoke | `et9-tidal-flat.yaml`: 4×24, dx = dy = 2 m, bed z = 0.40 − 0.02·y (to −0.56), north stage η = 0.3 sin(2πt/3600) with `inflow: {closure: external}`, periodic x, Manning 0.03. Subsurface 12 layers (dz 0.1, stretch 1.1), rl1 sand. Inland y− `head` hydrostatic η = 0.2, seaward y+ `head` hydrostatic η = 0. Salinity: sea 35, inland 0, IC 35 (consistent surface/subsurface updates, `store`). `et-met`, reallocation off, dt = 2 s, dtg ≤ 10 s, 3 cycles, outputs every 60 s. Lane **noET** (same without the block). | Subtidal cells (bed < −0.36): Σ soil_evaporation ≡ 0. Supratidal (bed > 0.36): Σ open-water evaporation ≡ 0. Intertidal (−0.2 ≤ bed ≤ 0.2): ≥ 1 wet→dry transition per cycle per cell (from `wet_time` increments). wet_time + soil_evaporation_time ≤ t + 1e-9 everywhere. Water ledger (surface + gw) ≤ 1e-10 of cumulative \|stage exchange\|. Salt ledger (surface + store + subsurface) ≤ 1e-10·M. Substeps ET/noET ≤ 1.2. noET lane bitwise against its T5-start hash. | as stated | ≈ 20 s | Anti-vacuity: open-water and soil totals both > 1e-6 m³ |
| **et-x** batteries | `et-evapotranspiration/orient-base.yaml` (8×8 coupled, corner pond plus sloping bed, ET, salinity, `store`, 600 s) through the 8-transform battery, strict; `regression.s3.et_rank` (et9 base at n = 1/2/4, strict) | Existing battery classes; strict rank bound | 1e-12 local | ≈ 15 s | r1 `--expect` |
| **restart** | `s3-restart` lanes et3 (checkpoint during the transition) and et8b (checkpoint while the store is full) | bitwise | exact | ≈ 4 s | — |

`unit.et_gates.negative`: every Python reference (bulk chain against two independent evaluations; the t_d closed form; the salt closed forms) and every checker on corrupted traces.

**DoD T5:**
- [ ] et0 (part 1), et1–et4, et7 (E), et8 (a–d), et9, et-x and restart green;
- [ ] g4/g5 and `default_bitwise` green (legacy keys untouched);
- [ ] docs and coverage rows;
- [ ] Tier N green;
- [ ] snapshot `T5`.

---

## 10. T6 — Vegetation and transpiration

**Semantics.**
- **Static per-cell vegetation** with Beer–Lambert transmission τ = exp(−k·LAI).
- **Ground rates** (open water under canopy, or exposed soil) are scaled by τ.
- **Potential transpiration:** T_pot = (1 − τ)·E_x(r_s = r_c), with canopy resistance r_c = r_l/(0.5·LAI). Here E_x is the active potential function: bulk uses T_s for the canopy too; prescribed uses E_pot and ignores r_s.
- **Transpiration is a non-compensated volumetric sink** over root-zone cells: S_k = α_F(h_k)·b_k·T_pot/Δz_k, capped by supply. It runs regardless of inundation (O7).
- **Salt:** water leaves, salt stays (exact factor).

### Work items

| W | File(s) | Exact change |
|---|---|---|
| W6.1 | `Config.*`, `ConfigSchema.cpp`, `RunRecord.cpp` | `struct VegetationType { name; lai; extinction; leafResistance; rootDepth; enum RootProfile {Uniform, Linear} profile; bool hasH12; real_t h1, h2, h3, h4; real_t albedo /*T7*/; }`, `VegetationConfig { types; mapFile; laiFile (optional) }`. Map ids: integer, **0 = bare**, k ≥ 1 refers to `types[k−1]`. Schema: §13.C. |
| W6.2 | `atm/EtFormulas.hpp` | `beerLambert(k, lai)` = exp(−k·lai) (exactly 1.0 at lai = 0); `canopyResistance(r_l, lai)`; `feddes(h, h1, h2, h3, h4, hasH12)` piecewise linear (0 above h1; ramps h1→h2 and h3→h4; 1 on [h3, h2]; 0 at or below h4; without h1/h2, 1 for h ≥ h3); `rootFractionIntegral(profile, L, d_a, d_b)` closed forms (uniform (min(d_b,L) − d_a)/L; linear [(2Ld − d²)/L²] evaluated from d_a to min(d_b, L); 0 for d_a ≥ L) |
| W6.3 | `gw/Evapotranspiration.cpp` | Setup (host, static): per column, z_top = top face of the `ktop` cell. Per active cell, depth interval [d_a, d_b] from `bot3d`/`dz3d`; `rootFrac_(j,i,k)` = integral / Σ integral (renormalised when the column is shallower than L, with a warning count). τ per cell from type LAI or the override raster. Per substep (`updateEvapotranspiration`): `tPot_(j,i)` = (1 − τ)·E_x(r_c); for each root cell, `sink_(j,i,k) = min(feddes(hⁿ)·rootFrac·tPot/dz3d, (θⁿ − θr)/dtg)` [1/s]. Audit T_act = Σ sink·V·dtg; unrealised = Σ (rootFrac·tPot·A·dtg − sink·V·dtg). The soil channel and open-water channel rates are multiplied by τ (open water via a per-cell factor field read by W5.3; exactly 1.0 when vegetation is absent, so W5 outputs stay bitwise). |
| W6.4 | `Predictor.cpp`, `Corrector.cpp` | One sink array read by both. Predictor RHS for active cells: add the volume term **−dtg·rRho·sink·V** with the same sign convention as the corrector numerator. Corrector `updateWaterContent` numerator: `− dtg·rRho·sink·V`. Guard: `const bool tr = transpirationActive_;` with an empty view when off (bitwise). Unit test `GwModule.SinkInSaturatedConfinedCellDropsHeadBySinkOverSs`: one saturated cell (θ = θs, h > 0), Ss > 0, impermeable faces, density off, constant sink S, one step: h − hn = −S·dtg/Ss to 1e-12 relative and θ = θs to 1e-14. The corrector identity θ_new = (θs − dtg·S)/(1 + Ss·Δh/θs) equals θs exactly only when the predictor carries the same sink with the same sign, so this test fixes the predictor sign. |
| W6.5 | `RichardsSolver.cpp` (window), `SubsurfaceTransport.cpp` | `sinkWin_(j,i,k) += sink·V·dtg`. Transport update: `vgflux −= sinkWin` (water out, salt stays); exact factor for any cell with sinkWin > 0: `f = (vgflux + sinkWin)/vgflux`, widen `hi *= f`, anchor on `sOld`. Top cell: `f = (vgflux + vapWin + sinkWin)/vgflux`. Needs `consistentSub` (schema). |
| W6.6 | `Simulation.cpp` | `gw_mass_audit` appends `transpiration`, `transpiration_unrealised`. Closure identity becomes Δvolume = boundary_in − ss_storage + realloc − vloss − transpiration (documented and asserted by et4). Outputs: `atmosphere: transpiration_rate`, `transpiration_unrealised_rate` (2-D); `groundwater: root_uptake` (3-D, interval-mean sink in 1/s). Run-record line `transpiration_m3`. |
| W6.7 | docs | `theory/evapotranspiration.md` (partition, Feddes, root mapping on the terrain-following mesh, the salt-exclusion statement, O7), `parameters.md`, `input-data.md` (vegetation map ids), `case-setup.md` (marsh recipe: h₁/h₂ omitted for wetland species), coverage rows |

### Gates (authored first)

| Gate | Configuration | Asserts | Tolerance | Budget | Negative / anti-vacuity |
|---|---|---|---|---|---|
| **et0** (part 2) unit | `AtmEt.BeerLambertLimits`, `AtmEt.FeddesPiecewise` (each segment plus the no-h₁h₂ form), `GwEt.RootFractionsIntegrateToOneOnStretchedMaskedColumns` (1×2 grid, dz 0.05 stretch 1.2, bed 0 / −0.12 with `follow_terrain: false` so column 2 is masked; L = 0.33 cuts a cell; uniform and linear), `GwModule.SinkInSaturatedConfinedCellDropsHeadBySinkOverSs` | closed forms | 1e-14 | < 0.1 s | perturbation |
| **et4** (T lane) | et4 column with vegetation (LAI 2, k 0.5, r_l 100, L 0.4, linear; Feddes h₃ −10, h₄ −150) | ΔS = −∫(E + T); audit T + E = realised ET; unrealised exact | 1e-12 | ≈ 3 s | Drop T from the identity: fails |
| **et5** partition and location | (i) et2 column with a vegetation map of all zeros. (ii) 1×1×50 silt column (rl1 silt, dz 0.02), uniform roots L = 0.4 m, LAI 30, k = 1 (τ = 9e-14), Feddes α ≡ 1 range, 1 day. (iii) Contrast lane: same column bare, `potential: prescribed` equal to (ii)'s T_pot, limiter `none`. | (i) bitwise against et2. (ii) E_soil < 1e-12·T; the top cell's share of the total depletion ≤ 2·b_top (b_top = 0.02/0.4 = 0.05). (iii) top cell's share ≥ 0.5. | bitwise / as stated | ≈ 5 s | Lane (iii) is the real negative: a top-face removal at the same rate fails (ii)'s criterion |
| **et6** root distribution | 1×2 columns of et0's geometry, K = 1e-15 (no redistribution), Feddes h₃ = −1000, h₄ = −2000 (α_F ≡ 1), 10 substeps | Per cell Δθ_k·Δz_k = rootFrac_k·T_pot·Δt (uniform and linear lanes); Σ = T_act = T_pot | 1e-12 | < 1 s | Fractions by cell centre (synthetic): fails |
| **et7** (T part) stress | (a) column at h = −5 m (below h₄ = −2): T_act ≡ 0, unrealised = T_pot, E_soil > 0. (b) Per substep T_act ≤ T_pot. (c) Root zone at θr + 1e-4 with high T_pot: θ ≥ θr − 1e-14, unrealised exact. | exact | as stated | ≈ 3 s | — |
| **et8** (T lane) salt | et8a with vegetation and E off (`soil` omitted; T only) | Per cell c·θ·V = const; total salt constant; subs_adjust, subs_anchor ≤ 1e-12·M; salinity rises in root cells ∝ rootFrac (correlation over cells = 1 to 1e-10) | 1e-12 | ≈ 2 s | T removing salt (synthetic): conservation fails |
| **et9** (vegetated lane) | et9 with a marsh band 10 ≤ y < 20 m (LAI 2.5, k 0.5, r_l 100, L 0.3, Feddes without h₁/h₂, h₃ −4, h₄ −80) | et9 assertions plus: T > 0 in the band including inundated steps (O7); T ≡ 0 outside | as et9 | ≈ 20 s | — |
| **restart** | `s3-restart` lane et9-vegetated (checkpoint mid-cycle) | bitwise | exact | ≈ 4 s | — |

**DoD T6:**
- [ ] all T6 gates green;
- [ ] T5 gates still green, with the T5 bitwise outputs unchanged when vegetation is absent;
- [ ] docs;
- [ ] Tier N green;
- [ ] snapshot `T6`.

---

## 11. T7 — Penman–Monteith potential

**Formulation** (specific-humidity form, sharing the bulk chain's functions so the two modes are commensurable):
- Δ_q = dq_sat/dT evaluated at T_air: the analytic derivative of `saturatedSpecificHumidity` (Tetens).
- λ = `latentHeatOfVaporization(T_air)`.
- s = λΔ_q/c_pa, with c_pa = 1005.
- ρ_a = `airDensity(T_air, P)`.
- A = (1 − albedo)·R_s + ε(LW_in − σT_K⁴) − G, with ε = 0.97 (`kWaterEmissivity`), T_K = T_air + 273.15, and G = fraction·R_n.

```
λ ρ_w E_x(r_s) = [ s·A + λ ρ_a (q_sat(T_air) − q_a)/r_a ] / [ s + 1 + r_s/r_a ]
```

The surface albedo and G fraction are chosen per channel:
- water: `energy.albedo_water` and fraction `water`;
- soil: the soil type's `albedo` and fraction `soil`;
- canopy: the vegetation type's `albedo` and fraction `canopy`.

Soil limiter: `beta_lee_pielke` or `none` (alpha_barton is rejected with PM).

| W | File(s) | Exact change |
|---|---|---|
| W7.1 | `Config.*`, `ConfigSchema.cpp`, `RunRecord.cpp`, `atm/MetForcing.hpp` | `atmosphere.shortwave_incoming` (R_s, W/m²; distinct from the Q5 `shortwave` = absorbed, whose meaning is unchanged); `MetSample::shortwaveIncoming`. `EnergyConfig { albedoWater = 0.06; gFractionWater = 0, gFractionSoil = 0.1, gFractionCanopy = 0.1; }`. Remove the W5.1 "not available before T7" rule; add the §13.C PM rules. |
| W7.2 | new `atm/PenmanMonteith.hpp` | `penmanMonteithRate(T_air, P, q_a, R_s, LW_in, albedo, gFrac, r_a, r_s)` as above, plus `dqsatdT(T, P)` (analytic) |
| W7.3 | `Evapotranspiration.cpp`, `SurfaceSources.cpp` | Dispatch `E_x` by potential mode: open water (r_s = 0, water albedo), soil (r_s = 0 × limiter, soil albedo), canopy (r_c, canopy albedo) |
| W7.4 | docs | PM section in `theory/evapotranspiration.md` (q-form derivation, commensurability statement, the FAO-56 convention difference), `parameters.md` |

| Gate | Configuration | Asserts | Tolerance | Budget | Negative |
|---|---|---|---|---|---|
| **et0** (part 3) | unit `test_atm.cpp` | `AtmEt.PmEquilibriumLimit`: q_sat(T_air) = q_a gives λρ_wE = s·A/(s+1) exactly. `AtmEt.PmAerodynamicLimit`: A = 0. `AtmEt.PmHighResistanceLimit`: r_s = 1e12 gives E < 1e-12·E(0). `AtmEt.PmCommensurateWithBulkAtEnergyBalanceTemperature`: for A ∈ {50, 100, 200, 400} W/m², solve T_s* from A = λρ_a(q_sat(T_s)−q_a)/r_a + ρ_a c_pa(T_s − T_a)/r_a by bisection, with ρ_a at T_air; \|E_PM − E_bulk(T_s*)\|/E against (T_s* − T_a) has a fitted order ≥ 1.8; the shipped bulk's air-density-at-T_s convention difference is recorded, not gated. `AtmEt.DqsatDtMatchesCentralDifference` (1e-8). | as stated | < 0.1 s | Perturbed s fails |
| **et0** (part 4) | Python `scripts/test_et_gates.py::fao56` | `gate_et.fao56_hourly()` and `gate_et.fao56_daily()` reproduce the FAO-56 Chapter 4 worked examples for **hourly** data (N'Diaye, recalled ≈ 0.63 mm/h for 14–15 h) and **daily** data (Brussels 6 July, recalled ≈ 3.9 mm/d). Every intermediate (Δ, γ, e°, e_a, R_ns, R_nl, R_n, G, ET₀) is checked to the published precision. **The agent fetches FAO Irrigation and Drainage Paper 56 (fao.org, document X0490E), transcribes inputs, intermediates and page numbers into `benchmarks/et-evapotranspiration/reference/fao56-examples.yaml`, and treats the recalled values above as untrusted.** | published precision (last printed digit ± 1) | < 1 s | Each transcribed value perturbed by 1 % fails |
| **et10** model-level PM | 1×1×20 grass column (vegetation LAI = 2.88, the FAO 24·h for h = 0.12 m grass; extinction k = 11 so τ = e^{−31.7} ≈ 1.7e-14; r_l = 100 → r_c = 69.4 s/m, the FAO reference value), Feddes α ≡ 1, wet soil, `aerodynamic_resistance: fao56` (z_ref 2), constant met = FAO hourly example inputs. LW_in is set so that ε(LW_in − σT_K⁴) = −R_nl,FAO and R_s = R_s,FAO, albedo 0.23, G fraction 0.1. Run 3600 s. | (1) Model T_act/h equals `gate_et.pm_q_form()` (the W7.2 formula in Python with the model's constants) to 1e-10. (2) Model against FAO's published ET₀: \|Δ\| ≤ δ_conv + rounding, where δ_conv = \|PM_q(model constants) − ET₀,FAO chain\| computed in et0 part 4 (expected 1–3 %; the measured value goes to the tolerance file with 1.2× headroom by amendment). | as stated | ≈ 2 s | Albedo 0.23 → 0.20 fails (1) |

**DoD T7:**
- [ ] et0 parts 3–4 and et10 green;
- [ ] FAO transcription file committed with page references;
- [ ] PM + temperature rejected (unit);
- [ ] docs;
- [ ] Tier N green;
- [ ] snapshot `T7`.

---

## 12. T8 — Default-flip evaluation, release pipeline, records

| W | Action | Pass condition |
|---|---|---|
| W8.1 | **Reallocation DF** (§2.7). Build `build-s3-df` with `Reallocation reallocation = Reallocation::Off;`, then run Tier N, `stage3`, `stage3_nightly`, unit and mpi. | All green, with unchanged tolerances: keep the flip and write amendment S3-A(n) with the measured values. Otherwise revert the constant and file `stage3/df/realloc-report.md` (owner decides). |
| W8.2 | **Salt-store DF** (§2.7-5), same procedure with `SurfaceDryout surfaceDryout = SurfaceDryout::Store;`. If W8.1 kept its flip, also run the combined build. | Same rule. b6 ss/td are the decisive cases (their surface cells dry at ebb); report surface and pore salinity against their bounds. |
| W8.3 | Pin update | If either flip is kept, the default-bitwise pins (§2.6) keep their explicit legacy values, so `default_bitwise` stays meaningful. `docs/user-guide/parameters.md` default columns and `resolvedConfigYaml` defaults are updated. |
| W8.4 | **p1 lanes for new kernels**: `regression.s3.p1.{ob3,pd2coupled,et9}.omp4.{aij,aijkokkos}` (label `stage3_nightly`), compared against the serial run with the existing default-mode per-solver bounds (V2-A4) | green |
| W8.5 | **p5**: `scripts/check_forbidden.sh` (MemSpace/backend rules already cover the new files); `unit.p5` negative compile still green | green |
| W8.6 | **Sanitizers**: `scripts/run_sanitizers.sh` with the `stage3` label added to its label list | clean |
| W8.7 | **Coverage lockstep**: `tests/coverage/bc_matrix.csv` rows for eta × {zero_gradient, external} × {uncoupled, coupled} and for edge BCs on periodic axes (schema-rejected); `feature-coverage.md` rows for every T1–T7 feature and interaction listed in the work items | `unit.bc_matrix`, `unit.feature_coverage` green, zero empty cells |
| W8.8 | **Docs gate**: `doxygen` (zero warnings) + `mkdocs build --strict` | clean |
| W8.9 | **Release pipeline**: `scripts/ci_release_gate.sh build-s3` (the final default build) | PASS; wall time recorded |
| W8.10 | **Study migration note** `stage3/records/T8/study-migration.md` | For the `2026Frehg2` matrix: which keys to add (`inflow: {closure: external}` on the tide BC; `domain.periodic.x` for AX runs, with the x walls and `*_xe` soil types removed; the ET block; `reallocation`/`surface_dryout` if not default), the generator changes (periodic random field), and the expected behaviour changes. No study file is edited. |
| W8.11 | **Completion report** `docs/developer-guide/report-T8.md` | Gate inventory with final measured values, timings per tier, amendments list, DF outcomes, residual limitations |
| W8.12 | Final snapshot `T8`; full diff `git diff refs/stage3/T0-start` written to `records/T8/stage3-full.patch` | Owner commits |

**DoD T8:**
- [ ] W8.1–W8.12 done;
- [ ] `ci_release_gate.sh` PASS;
- [ ] zero unverified BC-matrix rows and zero empty coverage cells;
- [ ] DF outcomes recorded (kept or reported).

---

## 12A. T9 — Paired water–salt store for thin-layer removals (added by S3-A12, O15)

Status: awaiting owner approval (2026-10-31). §2 applies unchanged: gate-first, records under `stage3/records/T9/`, snapshot `T9`, no commits. The section is numbered 12A so that the existing §13–§18 references (some quoted in source comments) stay valid.

**Removal-path inventory** (grep the quoted code; line numbers drift):

| Path | Where | What it removes | Audit at T8 |
|---|---|---|---|
| P1 | `swe/WetDry.cpp` `cflLimiter`, isolated-cell rule `eta(j, i) = bottom(j, i); cflActive(j, i) = 1.0;` | water the solve put on a dry cell that has no wet neighbour | none |
| P2 | `swe/WetDry.cpp` `cflLimiter`, `if (diff > 0.0 && diff < minDepth)` | films thinner than `min_depth` after the solve | none |
| P3 | `swe/SurfaceSources.cpp` `swe_evaprain_clamp`, `if (diff < minDepth)` | the sub-`min_depth` remainder after rain and evaporation (positive part) | `clamped` (signed, mixed with the below-bed lift) |
| P4 | `gw/Reallocate.cpp` vent branch, the else of `if (cplDept(j, i) > 0.0 \|\| d > minDepth)` | sub-`min_depth` vents onto a dry cell (the water never reaches η) | gw ledger (`realloc_dropped`) |

Not a removal: `updateDepth`'s `if (d <= minDepth) d = 0.0;` keeps the film in η but reports `dept` = 0, so the transport sees the cell dry while P2 removes the film only at the next solve (the R8 mismatch). Not thin-layer removals, unchanged by T9: the below-bed lift `swe_eta_clamp` (`FreeSurface.cpp`, creates volume, `clamped`), the infiltration remainder cut in `Coupler::applySeepage` (`discarded`), and the IC clamp.

**Design rules (binding).** Per-cell store state: water V_s [m³] and salt M_s [psu·m³] (the existing `sStore_`). The paired store is active iff the salinity store is active (surface water + transport + `surface_dryout: store`). Otherwise every T9 code path is inert, so runs without transport (the b5 class) are untouched.

- **R1 — rain and evaporation act first and are untouched.** The following are applied and booked exactly as at T8:
  - rain, evaporation (the ET open-water channel or the legacy `surface_water.evaporation` keys) and condensation;
  - their dilution or concentration of the cell's salt.

  The store receives only the positive leftover that a removal path P1–P4 takes after those sources have acted. A cell that never reaches a removal path is bitwise unchanged from T8.
- **R2 — the store never funds evaporation.** Evaporative demand beyond the water present stays in its existing audited path (`clamped`, the ET supply cap). Evaporation never debits V_s.
- **R3 — water and salt park together.** Every positive removal ΔV (P1–P4) moves into V_s in the same step as its salt c_park·ΔV moves into M_s.
  - c_park is the cell's limited post-source concentration (W9.3 fixes the transport basis), never above the cell's limiter-window maximum.
  - Dry-branch mass above c_park·ΔV is booked in `surf_adjust`, as `discard` books all of it at T8.
- **R4 — salt without water.** Salt that reaches a dry cell with no paired water in the step is kept in M_s only when:
  - (a) evaporation removed water from that cell in that step (the O4 crust of complete evaporative dry-out); or
  - (b) it is W5.9 held-seepage salt, paired with the coupler's held volume (`seepAccum`, positive part).

  Any other unpaired salt follows the `discard` path (`surf_adjust`, audited).
- **R5 — release is mixing.** At the first step the cell holds surface water, V_s returns to η and M_s returns to the water in the same step. The transport divides by vBasis + V_rel, so c = (M + M_s)/(V + V_s), and the anchor term stays the existing basis-lag correction rather than becoming a salt source. Without evaporation, c is a convex combination of source concentrations.
- **R6 — return without rewetting.** A dry cell whose V_s/A ≥ `surface_water.wetting_face_depth` returns its store at the next surface step, where its faces can carry the water. Below that depth, parked water stays parked.
- **R7 — parked water is inert (O15).** It does not flow, infiltrate, exchange heat or evaporate. Under ET it is outside the open-water channel and `wet_time`; the soil channel sees the cell as exposed, so complementarity is unchanged.
- **R8 — one predicate.** The transport's store-in predicate is "no surface water on the cell after the surface step" (η = z_b after the W9.2 sweep). It replaces T8's `!wet` (`vflux ≤ 0 || dept ≤ 0`). A cell holding water but no flux volume (wetted from dry by a vent, seepage or rain) keeps its salt in that water.
- **R9 — temperature.** With `modules.temperature`, the parked water's heat parks with it under R3/R5: the shared `ScalarSolver` store is enabled for the temperature instance, without R4. Returned water therefore does not dilute heat.

**Invariants (asserted by the gates).**
- **I1 (pairing):** R8 and R3 hold at every transport step.
- **I2 (source bound):** without evaporation, M_s/(V_s + held) never exceeds the maximum source concentration the cell has seen, so no released film exceeds it either.
- **I3 (conservation):** the identities below close to rounding at every audit row.
- **I4 (invariance):** lanes driven only by rain/evaporation, and lanes that never dry, are bitwise equal to T8.

**Identities (W9.1/W9.2).**
- **`discard`** (and every default-bitwise pin):
  Δvolume = rain − evaporation − boundary_outflow + bc_inflow + seepage + clamped − stage_deposit − film_removed,
  with film_removed = P1 + P2. P3 stays inside `clamped`, so that existing column remains bitwise.
- **`store`:**
  Δ(volume + store_volume) = rain − evaporation − boundary_outflow + bc_inflow + seepage + clamped − stage_deposit,
  with Δstore_volume = film_removed − store_release. Here film_removed counts P1, P2, the positive part of P3, and the W9.2 sweep (which parks the P4 vents); `clamped` keeps only the created below-bed part.
- **Salt (both modes):** Δ(surf_mass + surf_store + subs_mass) closes on the T5 ledger right-hand side as today.

### Work items

| W | File(s) | Exact change |
|---|---|---|
| W9.0 | records | Copy the S3-A12 repro (inputs, YAMLs, `analyze.py`, `vent_check.py`) from `/Users/zhili/Documents/Work/个人/论文/2026Frehg2/outputs/bugreport_dryout_store_v1/repro` into `stage3/evidence/dryout-store/`. Register every T9 gate (cases, `tolerances/stage3-ws.yaml`, harness subcommand `ws`, ctest entries, `scripts/test_ws_gates.py`), run them on the T8 binary, and save red logs to `records/T9/red/`. |
| W9.1 | `swe/WetDry.cpp`, `swe/SurfaceSources.cpp`, `swe/SurfaceSolver.hpp/.cpp`, `driver/Simulation.cpp`, `io/RunRecord.cpp` | Book every positive thin-layer removal per cell and per path, in all modes. A per-step Field2 `filmRemoved_` [m³] is filled by P1 and P2 (and, in store mode, by the positive part of P3 and the W9.2 sweep). New appended `mass_audit` column `film_removed` (cumulative). Run-record closure lines `film_removed_isolated_m3`, `film_removed_small_m3`, `film_removed_clamp_m3` (positive part of P3; informational under `discard`) and `film_removed_sweep_m3`. η arithmetic is unchanged outside store mode: fields stay bitwise, and `default_bitwise` differs only by the appended column. Attribute the S3-A12 vanished volume (1.1e-3 m³ in zoom A) to its path and record it in `report-T9.md`. |
| W9.2 | `swe/SurfaceSolver.*`, `swe/WetDry.cpp`, `swe/SurfaceSources.cpp`, `coupling/Coupler.*` | Store mode only. New Field2 `waterStore_` (V_s, m³) and per-step `parkedVol_`. P1, P2 and the positive part of P3 move their volume into `waterStore_` instead of deleting it. **End-of-step sweep**, after the last η change of the surface step (after `applySeepage`): every cell with 0 < η − z_b ≤ `min_depth` is set to η = z_b and its volume parked, which makes R8's predicate exact. Release (R5/R6) happens at the same point: a cell with η − z_b > `min_depth` and V_s > 0, or a dry cell with V_s/A ≥ `wetting_face_depth`, gets η += V_s/A, booked as `store_release`. `mass_audit` appends `store_volume` (instantaneous, Σ V_s) and `store_release` (cumulative). New output `surface: water_store` [m³]; checkpoint field `swe_water_store`. |
| W9.3 | `transport/SurfaceTransport.cpp`, `ScalarSolver.hpp/.cpp` | Store mode only. Store-in predicate R8 (reads `parkedVol_` and η − z_b). Dry-branch salt per R3 (c_park on the limited post-dilution basis, excess to `surf_adjust`) and R4 (evaporation crust and W5.9 held salt kept; other unpaired salt to `surf_adjust`). Release per R5, dividing by vBasis + V_rel. Consistent-exchange vent salt is credited only together with its water: the vent lands in η and is either wet or parked by W9.2 in the same step. Temperature instance (R9): the same paired store without R4. The W5.9 held path is unchanged. |
| W9.4 | `gw/Reallocate.cpp` | Store mode only. The vent branch deposits every dry-cell vent (the `d > minDepth` condition is dropped) and records it in the vent window, so its salt leaves the top soil cell through the existing vent path. The W9.2 sweep then parks the vent with its salt in the same surface step. Under `discard`/legacy nothing changes: `realloc_dropped` keeps the discard and the gw ledger stays bitwise. |
| W9.5 | `Config*`, `ConfigSchema.cpp`, docs | No new key: `surface_dryout: store` now means the paired store. New schema rule for the `water_store` output (§13.C). Docs: `theory/transport.md` (paired store, R1–R9, I1–I4); `theory/surface-water.md` (removal paths P1–P4, `film_removed`); `parameters.md`; `output.md` (columns, identities); `agents/case-setup.md` (replaces the no-ET store caveat); `feature-coverage.md` rows for store × {no ET, legacy evaporation, ET, rain, coupled, surface-only, periodic, temperature, restart, rank}; `stage3/records/T9/study-migration-addendum.md` (the v30 matrix can drop its explicit `discard` once W9.7 keeps `store`). |
| W9.6 | T5–T8 gates that resolve `store` | Re-run all of them and adjudicate every change by amendment with data: et8b's closed form gains V_r (S3-A12); et8d, et9, et9-veg and et-x are re-measured; b6 ss/td are run under `store`. No tolerance edit without an amendment. |
| W9.7 | `build-s3`, plus `build-s3-df` for the comparison | §2.7 re-evaluation of `store` as the shipped default under the paired semantics: unit, mpi, Tier N, `stage3` and `stage3_nightly`, with tolerances unchanged; b6 ss/td and ws1n are decisive. All green → `store` stays the default. Any failure → the default constant becomes `discard`, `stage3/df/store-T9-report.md` is written per §2.7-3, and the owner decides. |
| W9.8 | records | `dod-T9.md`; `report-T9.md` (gate inventory with measured values, the W9.1 path attribution, timings); mkdocs nav; snapshot `T9`; `records/T9/phase.patch`. |

### Gates (authored first) — `benchmarks/ws-water-salt-store/`, `tolerances/stage3-ws.yaml`

| Gate | Configuration | Reference / asserts | Tolerance | Budget | Negative / anti-vacuity |
|---|---|---|---|---|---|
| **ws0** unit (`test_transport.cpp`, `test_swe.cpp`) | `TransportStore.PairedReleaseMixesExactly`; `TransportStore.ParkedConcentrationWithinLimiterWindow`; `TransportStore.UnpairedSaltToAdjustWithoutEvaporation`; `TransportStore.UnpairedSaltKeptAfterEvaporativeDryout`; `TransportStore.NegativeResidualNeverStored`; `TransportStore.TemperatureParksWithItsWater`; `SweModule.ThinLayerRemovalsBookedPerPath` (P1, P2, P3 on constructed fields); `SweModule.StoreReturnsAtWettingFaceDepth` | Closed forms: c = (M + M_s)/(V + V_s); booked volume = removed volume | 1e-14 | < 0.1 s | `unit.ws_gates.negative` mirrors each closed form. A release that divides by vBasis alone fails. |
| **ws1** source bound (per-PR) | Small coupled section, authored at W9.0, with the S3-A12 ingredients: a dry seepage face fed by reallocation vents, a heterogeneous window, tidal stage with `external` closure, consistent updates, density coupling, no ET | (1) Max surface and subsurface salinity ≤ c_max·(1 + 1e-12) at every frame and monitor row, with c_max = 35. (2) Zero CFL warnings and zero clip warnings. (3) I2 on the store: M_s ≤ c_max·V_s per cell. | as stated (measured at authoring; any relaxation by amendment) | ≤ 30 s | **Red required:** the T8 binary must exceed c_max; otherwise the case is redesigned, and if no ≤ 30 s case reproduces it, C-WS applies. Anti-vacuity: Σ parked volume and the number of release events stay above floors set at authoring (1/10 of the measured value, by amendment). |
| **ws1n** source bound (nightly, `stage3_nightly`) | The S3-A12 repro `store.yaml`: 2-D H2 beach, 0 → 28 800 s | As ws1 | As ws1 | ≈ 10 min | Red: 72 850 psu and CFL 81 on the T8 binary. Anti-vacuity: parked volume ≥ 1e-4 m³ and ≥ 100 release events (re-set at authoring). |
| **ws2** ledgers | Asserted on the ws1, ws1n, ws3 and ws4 lanes, plus a `discard` lane of ws1 | Water: the mode's identity (above) at **every** audit row, not only at MSL-stage instants. Salt: the Δ(surf_mass + surf_store + subs_mass) identity. | 1e-10 of (cumulative \|stage exchange\| + rain + evaporation) | (inside the lanes) | Red: `film_removed` is absent on the T8 binary; the S3-A12 gap (1.2–3.2 % of the exchange) is the documented size. Negative: dropping `film_removed` from the identity fails on the `discard` lane. |
| **ws3** rain/evaporation invariance (O15) | **(a)** 4×4 closed basin, 0.2 m of water at 35 psu, prescribed rain series plus legacy prescribed evaporation, never dries; lanes `store` and `discard`. **(a-ET)** The same basin with the ET open-water channel plus rain, under `store`. **(b)** 1×1 pond, 5 mm at 35 psu, ET open water to dryness (et8b geometry). **(c)** Rain h_r after (b). **(d)** As (b) with legacy prescribed evaporation (no ET). | (a): `store` ≡ `discard` bitwise. (a-ET): bitwise against its T8-binary hash (I4). (b) and (d): Σ evaporation + V_s,final = V₀ (1e-12·V₀); V_s,final ≤ `min_depth`·A; M_s = M₀ (1e-12); V_s never decreases while the cell is dry (R2, R7; read from `water_store`). (c): first wet step c = M₀/(V_s + h_r·A) (1e-12); total salt constant. | as stated | ≈ 5 s | Anti-vacuity: in (a), salinity moves by > 1 psu through dilution and concentration; in (b), the cell dries. Negative: a store debited by evaporation fails (b). |
| **ws4** path coverage | One small lane per path: P1 (isolated solve-wetted cell), P2 (ebbing film on a slope), P3 (evaporation leaving a sub-`min_depth` film), P4 (sub-`min_depth` vent on a dry coupled seepage face); P2 repeated with `modules.temperature` (R9) | Per lane: the path's parked volume = its removed volume (1e-12); salt is paired at the cell's concentration (1e-12); the release mixes (R5); heat parks with the water (R9). | as stated | ≈ 10 s | Anti-vacuity: each lane's own path fires (above a floor). A lane whose path parks nothing fails. |
| **ws-x** rank and restart | `regression.s3.ws_rank`: ws1 at `mpi_ny` = 1/2/4, strict (`stage3_nightly`). `s3-restart` lane ws1: checkpoint taken while ≥ 10 cells hold V_s > 0. | Rank: the S3-A6 field classes, with Σ V_s and Σ M_s to 1e-12 relative. Restart: the S3-A6 restart standard, with `water_store` bitwise. | as stated | ≈ 90 s nightly / ≈ 4 s | r1 `--expect` |
| **ws5** pins | `regression.s3.default_bitwise` and the et9 noET pin lane (both pin `discard`) | Bitwise, except the appended `film_removed` column (§2.6 exception 1) | exact | (existing) | — |

**DoD T9:**
- [ ] red evidence for ws1n, ws1 (or C-WS) and ws2 in `records/T9/red/`;
- [ ] ws0–ws5 and ws-x green;
- [ ] W9.1 path attribution of the S3-A12 volume recorded;
- [ ] W9.6 adjudications recorded by amendment;
- [ ] W9.7 DF outcome recorded (kept or reported);
- [ ] docs and coverage rows;
- [ ] Tier N green;
- [ ] snapshot `T9`.

---

## 12B. T10 — Solve-consistent surface transport (added by S3-A14, O16)

Status: awaiting owner approval (2026-10-31). §2 applies unchanged: gate-first, records under `stage3/records/T10/`, snapshot `T10`, no commits.

**Diagnosis (code order, measured at T9).** In one surface step:

1. `solveFreeSurface` moves water across every face. The face volume is
   Q_x(j, i) = dt·Asx·Ex − Sxp·(η_{i+1} − η_i)^{n+1}. It is built from the
   **start-of-step** face areas and the new η, and is the same expression
   `accumulateBoundaryFluxes` already uses for the stage and edge faces.
   That is the water that actually moved.
2. `cflLimiter` and `evapRain` (P1–P3), then (coupled) the groundwater
   substeps, `applySeepage`, the stage routing and the T9 sweep/release.
3. `updateVelocity` recomputes depths and face areas from the **new** η,
   rebuilds u = (Ex − g·dt·effh_new·Δη)·Dx, and `applyVelocityLimiters`
   zeroes faces (a sender below `wetting_face_depth`, vanished faces,
   wetting-limited cells). Then Fu = u·Asx_new.
4. The transport advects salt with that Fu (`surface_update: consistent`
   also builds its flux volume from it).

So the transport's face volumes are not the ones that moved the water.
Reordering steps 3 and 4 would not fix this, because Fu is rebuilt on the
new geometry whether or not the limiter runs. The mismatch shows up in two
books:

- `surf_anchor`, the re-anchoring of the concentration from the transport's
  volume onto the real depth, creates or destroys salt. On ws1n (the
  owner's repro) it reaches **+28.4 psu·m³ = 2.4 % of the surface salt
  after 0.64 T**, the same in both modes (store +28.39, discard +28.23).
- At drying faces, water the solve moved arrives salt-free, and its salt is
  stranded in the drying cell, where it is discarded or capped into
  `surf_adjust`: −0.077 psu·m³ on ws1n, and 0.109 psu·m³ = 24 % of the
  exfiltrated salt on the T5 1 × 2 et8d strip (S3-A13).

**Design rules (binding).** Active under the new token
`transport.surface_update: conservative`. `legacy` and `consistent` stay
bitwise, so every existing gate and the default-bitwise lanes are untouched.
W10.6 then evaluates whether `consistent` should resolve to it.

- **C1 — record the solve's face volumes.** Right after `fillAndSolve`,
  before `cflLimiter`, store Q_x and Q_y [m³ per step] for every face from
  the assembled continuity: interior faces, open edges, stage faces and
  periodic seams. The face slots match Fu/Fv, and the west/south interface
  slots are halo-exchanged like Fu. One face, one number: the volume cell i
  sends east is exactly what cell i + 1 receives from the west.
- **C2 — salt moves with that water.** The advective mass increment is
  Σ Q·c_face, with the upwind or superbee face value chosen by the sign of
  Q and the far-neighbour guards unchanged. The flux volume is
  vflux = vsn + Σ Q, which is the post-solve, pre-manipulation water
  volume. Diffusion keeps its face areas.
- **C3 — every post-solve change keeps its own explicit book.** P1/P2/P3
  removals (T9: parked with their salt under store, `surf_adjust` under
  discard), rain and evaporation (the dilution volumes), seepage, vents and
  the held exchange (the exchange term), stage routing (`stage_deposit`, at
  the cell's concentration), and the T9 sweep and release. The velocity
  limiters then act only on the **next** step's momentum and never remove
  salt from this step's transport.
- **C4 — overdrawn cells.** Where the solve drains a cell below its bed
  (the below-bed lift `swe_eta_clamp`, booked in `clamped`), the lifted
  water is salt-free. The donor's outflow salt is capped at its available
  mass, and any negative remainder follows the existing negative-residual
  path (`surf_adjust`). Its size is measured and recorded.
- **C5 — the anchor becomes a basis check.** With C1–C3, vflux plus the
  booked post-solve volumes equals dept·A after the step, so `surf_anchor`
  must fall to rounding. A residual anchor is a defect, not a correction.
- **C6 — temperature.** The temperature instance reads the same Q (its own
  dilution rule unchanged).

**Invariants.**
- **J1:** a uniform scalar stays uniform to rounding through wetting,
  drying and limiter events.
- **J2:** `surf_anchor` ≤ 1e-12·M, and `surf_adjust` contains only the
  booked P1–P3/R4/C4 terms.
- **J3:** the boundary salt flux equals c_face times the very volumes in
  `mass_audit.boundary_outflow`, face by face.
- **J4:** salt that leaves a drying cell with its water arrives in the
  neighbour.

### Work items

| W | File(s) | Exact change |
|---|---|---|
| W10.0 | records, cases | Gate-first: author the sx cases and tolerances (`benchmarks/sx-solve-consistent/`, `tolerances/stage3-sx.yaml`), the harness subcommands `sx*`, the ctest entries and `scripts/test_sx_gates.py`; red logs on the T9 binary to `records/T10/red/` |
| W10.1 | `swe/FreeSurface.cpp`, `SurfaceSolver.*` | C1: a post-solve pass that fills `Qx_`/`Qy_` from the assembled continuity (read-only beside the pinned solve; exchanged; inert unless the token is `conservative`). Unit: telescoping, and Δη·A = −ΣQ on a source-free step |
| W10.2 | `transport/SurfaceTransport.cpp`, `ScalarSolver.*`, `driver/Simulation.cpp` (wiring) | C2/C3/C4/C6 behind the new token; the legacy and consistent expression trees unchanged (separate branch / kernel per the S3-A8.1/A9.7 fma rule) |
| W10.3 | `Config*`, `ConfigSchema.cpp`, docs | `transport.surface_update: conservative` (§13.C); resolved echo; `theory/transport.md` (solve-consistent update, J1–J4), `theory/surface-water.md` (the step order), `parameters.md`, `agents/case-setup.md` |
| W10.4 | gates | Re-run the T5–T9 store/ET/ws lanes under `conservative` (they keep `consistent` in their files): adjudicate every change by amendment with data |
| W10.5 | records | The anchor/adjust before–after table on ws1n, the et8d 1 × 2 strip, b6 ss/td and the v30 H0/H2 smokes in `report-T10.md` |
| W10.6 | DF | §2.7-style evaluation of making `consistent` resolve to `conservative` (owner decides: the v30 matrix uses `consistent`) — every tier, tolerances unchanged; b6 ss/td and ws1n decisive |
| W10.7 | records | `dod-T10.md`, `report-T10.md`, mkdocs nav, snapshot `T10`, `records/T10/phase.patch` |

### Gates (authored first) — `benchmarks/sx-solve-consistent/`, `tolerances/stage3-sx.yaml`

| Gate | Configuration | Asserts | Tolerance | Negative / anti-vacuity |
|---|---|---|---|---|
| **sx0** unit | `SweModule.SolveFaceVolumesTelescope` (interior faces cancel; Σ Q = −Δη·A per cell on a source-free step); `TransportModule.SolveConsistentUniformStaysUniform`; `TransportModule.SolveConsistentDryingFaceCarriesSalt` | closed forms | 1e-14 | a transport using Fu instead of Q fails the drying-face test |
| **sx1** uniform scalar through drying | ws1's coupled beach at uniform 35 psu (surface, soil, sea BC) under `conservative` | J1: c ≡ 35 at every frame and monitor row; J2: anchor and unbooked adjust ≤ 1e-12·M | as stated | red on T9 (`consistent`): anchor/adjust measured nonzero; anti-vacuity: drying-face and limiter events counted above a floor |
| **sx2** drying-face salt | the T5 et8d 1 × 2 strip (`min_depth` 5 mm) under `conservative` | J4: `surf_adjust` == 0 (no stranded salt), total salt exact, the neighbour receives the salt the solve's water carried | 1e-12 | red on T9: −0.109 psu·m³ |
| **sx3** ws1n anchor | the owner's repro under `conservative` | J2 on a production run: anchor ≤ ε·M (set at authoring by amendment); the ws1 source bound, I2 and ws2 ledgers still hold | as stated | red on T9: anchor +28.4 psu·m³ |
| **sx4** non-uniform advection across a front | a tracer step crossing a wetting/drying front on a slope | mass exact, no new extrema, J3 boundary salt = c × boundary volumes | 1e-12 | — |
| **sx5** temperature | sx1 with `modules.temperature` | uniform T stays uniform; heat ledger closes without an anchor | 1e-12 | — |
| **sx-x** | rank (mpi_ny 1/2/4, on a trajectory-stable base per S3-A13), restart (bitwise), periodic seam (pd2-coupled under `conservative`) | the S3-A6 classes; restart bitwise | as stated | — |

**DoD T10:**
- [ ] red evidence (sx1–sx3) in `records/T10/red/`;
- [ ] sx0–sx5 and sx-x green;
- [ ] W10.4 adjudications recorded by amendment;
- [ ] W10.6 outcome recorded (kept or reported);
- [ ] docs and coverage rows;
- [ ] Tier N green;
- [ ] snapshot `T10`.

---

## 13. Schema specification (all phases)

### 13.A Feature A (T2)

```yaml
boundary_conditions:
- name: tide-eta
  region: {polygon: [[-1.0, 79.625], [81.0, 79.625], [81.0, 81.0], [-1.0, 81.0]]}
  target: surface
  kind: eta
  value: {series: {file: input/tide.dat}}
  inflow:                        # optional; kind eta only
    closure: external            # zero_gradient (default) | external
    velocity:                    # optional, external only; default u = v = 0
      u: {constant: 0.0}         # [m/s] along +x (constant | series)
      v: {constant: 0.0}         # [m/s] along +y
```

| Rule | Error message |
|---|---|
| `inflow` on kind ≠ eta | `boundary_conditions[i].inflow: applies to kind eta only (it sets the momentum of the external reservoir behind a prescribed stage)` |
| `velocity` with `zero_gradient` | `boundary_conditions[i].inflow.velocity: only valid with closure external (zero_gradient takes the inflow momentum from the interior)` |
| Adjacent stage cells, different closures (setup fatal) | `boundary conditions '<a>' and '<b>': stage cells share a face but declare different inflow closures` |
| External BC with no reservoir face (setup warning) | `eta boundary '<name>': closure external has no reservoir face (a one-cell-thick interior stage line); it behaves as zero_gradient` |
| Run-end warning (zero_gradient, R < 0.5, max v > 0.1 m/s) | `eta boundary '<name>': inflow entered at up to <v> m/s while paying <R>·v²/2g of entrance head (stage_ke_influx <E> J); if the water beyond the boundary is still, use inflow: {closure: external}` |

### 13.B Feature B (T3)

```yaml
domain:
  periodic: {x: true, y: false}  # optional; default false/false
```

| Rule | Message |
|---|---|
| nx < 3 (or ny) | `domain.periodic.x: needs nx >= 3 (got <n>)` |
| BC selecting only periodic-axis faces (setup fatal) | `boundary condition '<name>': its polygon selects only faces on the periodic x edges, which domain.periodic.x makes interior; remove the condition` |
| Seam discontinuity (warning) | `domain.periodic.x: <field> jumps by <s> across the x seam (largest interior jump <m>); the seam is a discontinuity, not a periodic continuation` |

### 13.C Reallocation, ET, vegetation, salt store (T1, T5–T7)

```yaml
groundwater:
  reallocation: off              # legacy (v2.0) | off (T8 may flip the default)

atmosphere:                      # the single met block
  air_temperature: {constant: 25.0}
  surface_temperature: {constant: 30.0}      # potential bulk only
  pressure: {constant: 101.3}
  relative_humidity: {constant: 0.3}         # or specific_humidity
  wind_speed: {constant: 8.0}
  wind_speed_floor: 0.5
  shortwave_incoming: {series: {file: input/rs.dat}}   # NEW (T7); penman_monteith
  longwave_in: {series: {file: input/lw.dat}}          # penman_monteith
  aerodynamic_resistance:        # NEW (T5); ET consumers only
    law: liu2006                 # liu2006 (default) | fao56 | log_profile
    reference_height: 2.0        # fao56, log_profile [m]
    z0m: 0.01                    # log_profile [m]
    z0h: 0.001                   # log_profile [m]
    displacement: 0.0            # log_profile [m]

evapotranspiration:              # NEW (T5)
  potential: bulk                # bulk | penman_monteith (T7) | prescribed
  prescribed: {series: {file: input/epot.dat}}   # potential prescribed only [m/s]
  open_water: true               # default true
  soil:                          # omit = no soil channel
    region: {polygon: [[0, 0], [80, 0], [80, 80], [0, 80]]}   # optional, default all
    limiter: alpha_barton        # alpha_barton (bulk only) | beta_lee_pielke | none
  energy:                        # penman_monteith only (T7)
    albedo_water: 0.06
    ground_heat_fraction: {water: 0.0, soil: 0.1, canopy: 0.1}

vegetation:                      # NEW (T6), optional
  types:
  - name: marsh
    lai: 2.5
    extinction: 0.5
    leaf_resistance: 100.0       # r_l [s/m]
    root_depth: 0.30             # [m] below the column's top surface
    root_profile: linear         # uniform | linear
    feddes: {h3: -4.0, h4: -80.0}   # [m]; optional h1 >= h2 >= h3 (anaerobic cutoff)
    albedo: 0.20                 # penman_monteith only
  map: {file: input/veg.dat}     # integer ids, j*nx+i; 0 = bare, k = types[k-1]
  lai: {file: input/lai.dat}     # optional per-cell override

soil:
  types:
  - {name: mud, ..., albedo: 0.15, theta_fc: 0.38}   # NEW optional keys

transport:
  surface_dryout: store          # NEW (T5): delete (v2.0) | store (O4)
  surface_update: conservative   # NEW token (T10, O16): legacy | consistent | conservative
                                 # T9 (S3-A12, O15): store = the paired water–salt store

output:
  monitor_interval: 60           # NEW (T4) [s]; absent/0 = every step
  variables:
    surface: [seepage_rate, water_store]           # NEW (T4; water_store T9)
    groundwater: [realloc_created, realloc_dropped, vloss, root_uptake]   # NEW (T1, T6)
    transport: [salt_store]                        # NEW (T5)
    atmosphere: [evaporation_rate, soil_evaporation_rate, soil_evaporation_unrealised_rate,
                 transpiration_rate, transpiration_unrealised_rate,
                 wet_time, soil_evaporation_time]  # NEW section (T5/T6)
```

| Rule | Message |
|---|---|
| `reallocation: off` with `reallocation_surplus` / `reallocation_gradient` (**warning**) | `groundwater.reallocation_surplus: has no effect with reallocation: off (nothing is transferred)` / `groundwater.density_coupling.reallocation_gradient: has no effect with reallocation: off` |
| ET block with `surface_water.evaporation` or `groundwater.evaporation` | `evapotranspiration: replaces surface_water.evaporation and groundwater.evaporation; give one form` |
| `groundwater.evaporation` coupled (kept) | existing message + ` — coupled soil evaporation is configured through the evapotranspiration block` |
| ET without the atmosphere block (bulk/PM) | `evapotranspiration.potential <mode>: requires the atmosphere block` |
| bulk without `surface_temperature` | `evapotranspiration.potential bulk: requires atmosphere.surface_temperature` |
| PM without `shortwave_incoming` / `longwave_in` | `evapotranspiration.potential penman_monteith: requires atmosphere.<key>` |
| prescribed without `prescribed` | `evapotranspiration.potential prescribed: requires evapotranspiration.prescribed` |
| `alpha_barton` with PM or prescribed | `evapotranspiration.soil.limiter alpha_barton acts inside the bulk humidity difference; use beta_lee_pielke or none with potential <mode>` |
| PM + `modules.temperature` | `evapotranspiration: penman_monteith with modules.temperature is not available (the heat module's latent flux uses the bulk rate at the simulated temperature; the combination is ungated)` |
| bulk + `modules.temperature` + r_a law ≠ liu2006 | `atmosphere.aerodynamic_resistance: with modules.temperature the evaporated volume and the heat module's latent flux must share the liu2006 resistance` |
| `soil` without `modules.groundwater` | `evapotranspiration.soil: needs modules.groundwater` |
| `open_water` without `modules.surface_water` | `evapotranspiration.open_water: needs modules.surface_water` |
| ET + salinity without the consistent updates of the modules present | `evapotranspiration with salinity transport requires transport.surface_update: consistent (with modules.surface_water) and transport.subsurface_update: consistent (with modules.groundwater)` |
| ET + salinity with `surface_dryout: delete` | `evapotranspiration with salinity transport requires transport.surface_dryout: store (salt left by evaporating water must be kept)` |
| `surface_dryout: store` without a surface module or transport | `transport.surface_dryout: store needs modules.surface_water and modules.transport` |
| `surface_dryout: store` with `transport.bounds.max` (warning) | `transport.bounds.max clips concentrations above <max>; salt released from the store into thin films may be clipped (booked in surf_adjust)` |
| `output.variables.surface: water_store` without the store (T9) | `output.variables.surface: water_store needs transport.surface_dryout: store (with modules.surface_water and modules.transport)` |
| `vegetation` without the ET block or without `modules.groundwater` | `vegetation: has no effect without evapotranspiration` / `vegetation: transpiration needs modules.groundwater` |
| Feddes ordering | `vegetation.types[i].feddes: needs h1 >= h2 >= h3 > h4 (h1, h2 optional, given together)` |
| Map id without a type (setup fatal) | `vegetation.map: id <k> at cell (<j>, <i>) has no type (ids are 1..<n>, 0 = bare)` |
| Root depth beyond the column (setup warning) | `vegetation '<name>': root depth exceeds the column at <n> cells; truncated and renormalised` |
| θ_s < 0.375 with `alpha_barton` and `open_water` (warning) | `soil type '<name>': alpha_1(theta_s) = <a> < 1; the evaporation rate jumps by <(1-a)%> when a cell switches from inundated to exposed` |
| `output.variables.atmosphere` entries without their feature | `output.variables.atmosphere: <var> needs <evapotranspiration.soil / vegetation>` |
| `realloc_*` outputs with `reallocation: off` (warning) | `output.variables.groundwater: realloc_created/realloc_dropped are identically zero with reallocation: off` |

---

## 14. Gate inventory and budget

Predicted serial costs (1 rank unless the gate is a rank test, `OMP_NUM_THREADS=1`, M3). Every case is < 60 s; the largest are ob5 (≤ 30 s), et9 and et9-veg (≈ 20 s each).

| Family | ctest tests (label `stage3` unless noted) | Predicted wall |
|---|---|---|
| shared | `regression.s3.default_bitwise` (`stage3_local`), `regression.s3.restart` (lanes from T3–T6, T9) | 40 s + 18 s |
| rl | `unit` (rl1), `regression.s3.rl2`, `.rl3`, `.rl4`, `.rl5.gw_orient_off`, `.rl5.coupled_salt_rank_off`, `unit.rl_gates.negative` | 17 s |
| ob | `unit` (ob2, ob4, ob7), `regression.s3.ob1`, `.ob3`, `.ob5`, `unit.ob_gates.negative`; `regression.s3.ob8.b6_{ss,td}` (`stage3_nightly`) | 35 s (nightly ≈ 35 min) |
| pd | `regression.s3.pd1`, `.pd2`, `.pd3`, `.pd4`, `.pd5`, `.pd6`, `.pd8`, `mpi.halo` periodic case, `unit.pd_gates.negative` | 50 s |
| dg | `regression.s3.dg1`, `.dg2`, `unit.dg_gates.negative` | 5 s |
| et | `unit` (et0 parts 1–3, et6), `regression.s3.et1`…`et10`, `.et_x`, `.et_rank`, `unit.et_gates.negative` (incl. FAO transcription) | 85 s |
| ws (T9) | `unit` (ws0), `regression.s3.ws1`, `.ws3`, `.ws4`, `unit.ws_gates.negative` (ws2 and ws5 assert inside these lanes and the existing pins); `regression.s3.ws1n`, `.ws_rank` (`stage3_nightly`) | 45 s (nightly ≈ 12 min) |
| sx (T10) | `unit` (sx0), `regression.s3.sx1`, `.sx2`, `.sx4`, `.sx5`, `unit.sx_gates.negative`; `regression.s3.sx3`, `.sx_rank` (`stage3_nightly`) | ≈ 40 s (nightly ≈ 12 min) |
| **Stage-3 per-PR total** | | **≈ 4.9 min serial** with T9 (< 5 min; ≈ 4.1 min through T8); ≈ 1.7 min at ctest `-j8` |
| Existing per-PR tier (measured, T0 re-measures) | 93 tests | ≈ 23 min wall |
| Existing nightly tier (measured, partial) | 7-test batch | ≈ 56 min + b6 bjacobi ss/td + b5 24 h (T0 measures) |

---

## 15. Documentation deliverables by phase

| Phase | Pages |
|---|---|
| T0 | in-repo plan, mkdocs nav, `testing.md` (stage3 labels, override mechanism, default-bitwise gate) |
| T1 | `theory/groundwater.md` (post-allocation modes), `parameters.md`, `configuration.md`, `output.md` |
| T2 | `theory/surface-water.md` (inflow closure, deposit routing, diagnostics), `parameters.md`, `configuration.md`, `output.md`, `agents/case-setup.md`, `feature-coverage.md`, `bc_matrix.csv` |
| T3 | `grid-indexing.md`, `halo-exchange.md`, `configuration.md`, `parameters.md`, all theory pages (seam statement), `symmetry-exemptions.md`, `agents/case-setup.md`, coverage |
| T4 | `output.md`, `adding-an-output-variable.md`, `agents/postprocessing.md`, `parameters.md` |
| T5–T7 | new `theory/evapotranspiration.md`, `theory/atmosphere.md`, `theory/transport.md`, `parameters.md`, `configuration.md`, `input-data.md`, `output.md`, `agents/case-setup.md`, `agents/postprocessing.md`, coverage |
| T8 | `report-T8.md`, `parameters.md` defaults (if flipped), `AGENTS.md` (v2.0 → Stage-3 feature list), `README.md` feature list, `CITATION.cff` untouched (the owner versions) |
| T9 | `theory/transport.md` (paired store), `theory/surface-water.md` (removal paths, `film_removed`), `parameters.md`, `output.md`, `agents/case-setup.md`, coverage, `report-T9.md`, study-migration addendum |
| T10 | `theory/transport.md` (solve-consistent update, J1–J4), `theory/surface-water.md` (step order, face volumes), `parameters.md`, `agents/case-setup.md`, coverage, `report-T10.md` |

Every phase also writes `dod-T<n>.md` and `report-T<n>.md`.

---

## 16. Contingencies (each use requires an amendment with data)

| ID | Trigger | Response |
|---|---|---|
| C-RL | rl2 under `off` shows \|vloss\| > ε_v·T (the corrector overshoots θs or θr without reallocation) | Measure where (output `vloss`). If it is the over-saturation excess of saturated cells, report to the owner before changing anything: options are a smaller `dt_max`, or an `off`-mode excess-to-neighbour rule (a new design, owner decision). The DF evaluation then reports it. |
| C-OB5 | ob5 cannot show a ≥ 10× contrast in ≤ 30 s | Use the nightly D1 record (§6 ob5 row) |
| C-PD | pd2 roll invariance fails only under the default solver | Gates run strict (as designed). Record the default-solver measurement; do not loosen strict bounds. |
| C-ST | et8/et9 show solver failure or non-physical density feedback when the store releases into thin films | Report with data. Candidate (owner decision, O4 refinement): release proportional to depth up to a release depth. **Triggered 2026-10-30 outside the gates (S3-A12); the candidate is insufficient (406 psu, S3-A12) and is superseded by T9 (O15).** |
| C-WS | ws1 cannot reproduce the T8 red (salinity above c_max) in ≤ 30 s | ws1 becomes a nightly-only lane beside ws1n (amendment); the per-PR tier keeps ws0 and ws2–ws4 |
| C-INERT | R7 (inert parked water under ET) moves an ET gate beyond its W9.6 closed-form adjustment | Report with data. The owner's alternative (an O15 refinement) lets parked water evaporate, booked as evaporation |
| C-SX | A post-solve change cannot be given an explicit salt book under C3 (an unbooked η change found by J2), or C4's capped outflow grows beyond film scale on a gate | Report with data; the change is booked (new audit term) or the cap is redesigned before T10 closes |
| C-FAO | FAO-56 PDF not retrievable | Owner supplies the PDF; et0 part 4 and et10 (2) wait; et10 (1) proceeds |
| C-BIT | `default_bitwise` differs by ≤ 1e-14 only in an edited kernel | §2.6 caveat procedure |
| C-TIME | Any Stage-3 case > 60 s, or the suite > 5 min | Coarsen or shorten the case with a refinement note; else reclassify it to `stage3_nightly` by amendment (the g5 precedent) |

---

## 17. Owner review items (the only owner touchpoints)

1. **Plan approval** (this document), including the proposal to apply the O9 default-flip protocol to the salt store (§2.7-5).
2. **Each phase close:** review `records/T<n>/phase.patch`, the DoD and the report. The agent proceeds to the next phase without waiting unless the owner asks for phase-by-phase approval (please state your preference).
3. **Decisions on demand:** any contingency in §16; any DF failure report (T8); any tolerance amendment.
4. **FAO-56 PDF** if C-FAO triggers.
5. **Final commit** of the tree after T8 (`records/T8/stage3-full.patch` and the snapshot refs provided). Optionally delete the `refs/stage3/*` refs afterwards (§2.2).
6. **T9 (S3-A12), open:** (a) approve phase T9 (§12A); (b) the interim default of `transport.surface_dryout` until W9.7 passes — the agent recommends reverting it to `discard` through §2.7, because the shipped `store` default blows up coupled tidal runs without ET. Decided 2026-10-31 (O15): the paired water–salt store, rain and evaporation untouched, parked water inert under ET.
7. **T10 (S3-A14), open:** (a) approve phase T10 (§12B); (b) at W10.6, whether `transport.surface_update: consistent` should resolve to the solve-consistent update (the v30 matrix uses `consistent`; until then it opts in with `conservative`).

---

## 18. Amendment log (Stage-3)

### S3-A1 (T0 close, 2026-10-27) — §2.5: phase-scoped Tier N at intermediate phase closes

**Measured driver.** The T0 full Tier N cost 3 h 12 min wall for the nightly batch
(10/10 PASS, 61,855 s·proc at `-j8` on 4 P-cores = 5.4× oversubscription) plus
8,186 s for the b5 24 h bjacobi record — about 5½ h per phase close, nine times.

**Change (owner-approved 2026-10-27).** At the T1–T7 phase closes, Tier N runs:
the b5 24 h bjacobi record (release step 3), b6 ss + td, and the nightly lanes
whose module the phase touched (g5 + transposed and g8 for groundwater-touching
phases; g5 again at T5–T7). The AMG and backend duplicate lanes
(`g1.b5.rain.sync.amg`, `g1.b6.*.amg`, `g1.b5_rank_invariance.default.amg`,
`p1.b6.ss.omp4.aijkokkos`) run only at T0 (done, green) and T8. No gate's
physics or tolerance changes — only when untouched lanes re-run. Scheduling:
`-j2` with the long tests first; `perf_baseline` excluded from `-j8` sweeps and
run solo (its T0 `-j8` failure was uniform +53…+79 % contention; solo PASS
123.6 s).

### S3-A3 (T1, 2026-10-27) — §5/§13.C: the disabled-reallocation token is `disabled`, not `off`

`off` is a YAML 1.1 boolean: PyYAML — the gate harness, the run-record checker, and
any user postprocessing — reads it as `false`, which corrupted the r1 round trip on
the first rl execution (`reallocation: false` rejected by the schema). The YAML
token is `disabled`; the C++ enum stays `Reallocation::Off`. The schema rejects
`off` explicitly, and `ConfigTest.ReallocationModeParsesRoundTripsAndRejectsLegacyOnlyKeys`
pins the rejection. Later Stage-3 keys must avoid every YAML 1.1 boolean literal
(y/yes/no/on/off/true/false) as an enum value.

### S3-A4 (T1, 2026-10-27) — §5 W1.2: the disabled-mode semantics, four designs measured

The plan's W1.2 ("every unsaturated cell re-derives h from the retention curve")
is not viable. Four designs were implemented and measured on the rl gates, the
gw-orient battery (a block saturated to the surface) and the coupled-salt strict
rank lane:

| # | Design | Measured outcome |
|---|---|---|
| 1 | h := h(θ) for every θ < 0.9999 θs cell | Orientation branch flips at the cutoff on water-table cells: 7e-4 m head deltas (bound 1e-12), ~2.1× substeps |
| 2 | As 1, blended into the predicted head over the last 0.1 % of saturation | Pins the capillary fringe to the unsaturated branch; the water table cannot pass: 17 m head deltas, θ deltas 2.5e-2, 2.6× substeps, vloss 13 m³ against 0.37 m³ of forcing |
| 3 | Legacy sweep minus every non-conservative action (no sends, no drops, no θ := θ(h)); the conservative isolated-cell h := h(θ) restore kept verbatim | Exactly conservative and deterministic everywhere; **in saturated regimes** the clamp mints water at draining cells and deletes at filling ones (±0.33 m³ per 0.0625 m³ cell) and rounding amplifies chaotically (O(1) orientation divergence); **in vadose regimes** clean (rl3 steadiness exactly 0, flux uniformity 9.3e-13, vloss exactly 0, substeps 1.000) |
| 4 | θ keeps compressed storage above θs; h from a continuous extended curve with the Ss line h = (θ−θs)/Ss | The 1/Ss noise amplification blows up (1e19 deltas) — pressure must come from the solve, never from θ |

**Adopted: design 3**, with an explicit domain of applicability (vadose-dominated
problems), documented in `theory/groundwater.md` and the `Reallocate.cpp` header.
The saturated-regime behaviour is a property of the PCA scheme — the legacy
θ := θ(h) restore is its consistency mechanism there, and a conserving
replacement is a scheme redesign that plan contingency C-RL reserves for the
owner. Gate rescope that follows from this:

- rl2 (saturated heterogeneous box): the disabled lane asserts the conservation
  identities (ledger, exact-zero realloc, θ bounds, per-cell sums) and **records**
  vloss and substeps as the owner's C-RL data (measured: vloss 5.743 m³ over
  T = 2.380 m³; substeps 562 vs 503). The vloss and substep **bounds** move to
  rl3 (vadose), where the mode is applicable (measured: 0 and 1.000).
- rl2's legacy-churn anti-vacuity floor is churn ≥ 0.1·T(legacy) (measured
  churn/T = 2.9; the earlier vloss-linked floor tied it to an out-of-domain
  number).
- rl5's orientation battery runs on a new vadose base
  (`benchmarks/rl-realloc/rl5-orient-base.yaml`), passing strict 1e-12; the
  saturated x-orient gw block is outside the mode's domain and no strict battery
  can pass there under any consistency-free scheme. The coupled-salt strict rank
  lane stays as specified and passes.
- rl3's ledger bound gains a 1e-12 m³ absolute floor (rounding over 864 audit
  rows at T = 0.19 m³; measured residual 4e-13).
- rl4's bound becomes three-part (≤ 1e-4·M_subs, ≤ legacy/5, legacy ≥ 6e-5·M
  anti-vacuity); measured: disabled 1.366 vs legacy 38.70 over M = 5.80e4 — a 28×
  reduction, not the clamp-only level the plan's W1.7 note predicted (the coupled
  exchange basis also contributes).

**Implication for the T8 default flip (O9):** in saturated regimes `disabled` is
not a drop-in replacement; the DF evaluation's b2/b3/b6 lanes will measure exactly
that. The owner decides at T8 with the DF data and the C-RL options.

### S3-A5 (T2, 2026-10-28) — §6/§13.A: first-green adjudications of the ob gates

All physics landed as specified (the ob2 unit gates passed on their first
execution); six adjudications at first green, each with its measured datum in
`tolerances/stage3-ob.yaml` or the gate code:

1. ob1 and ob3 run the strict solver (jacobi, rtol 1e-13): their gated
   identities are the momentum chain and the ledger, and the default
   bjacobi-icc rtol 1e-8 leaves solver-class residuals (ob3 ledger measured
   −6.6e-4 m³ default vs 2.5e-8 m³ strict over X = 3.0e4 m³).
2. ob1 harness face indexing: flat output index k holds face k+1.
3. ob3 `ledger_rel` 1e-12 → 1e-11 (1440-step rounding accumulation; 12×
   headroom; the 1e-6 m³ synthetic hole still fails).
4. ob5 `influx_rel` 1e-3 → 2e-2: the seeded slosh exchanges
   pressure-accelerated water through the entrance faces, which carries real
   kinetic energy (measured 6.2e-3 of KE0); the discriminating assertion is
   the zero_gradient contrast (gate 10×, measured 90×).
5. ob8 gates the EXPERIMENT criteria and records the golden-referenced
   residuals: the legacy golden embeds the stage-deposit kick, so `external`
   legitimately departs from it (td: one beach cell 0.0196 m vs the 0.015 m
   floor; salt mass 9.88 %) while improving against the experiment
   (MAE 0.0381 m vs the golden's own 0.0430 m). No zero_gradient tolerance
   changed.
6. `check_bc_matrix.py` canonical cells extended with the eta-inflow axis.

### S3-A15 (T10, 2026-10-31) — first-green adjudications of the conservative update

**W10.2 — the fma rule's third measurement.** The first implementation
used runtime ternaries (`qvol ? Q : Fu`) with value-identical false
branches. `regression.s3.et9` (noET pin, consistent + superbee +
subcycled) FAILED bitwise while `default_bitwise` passed — the ternary
perturbed the compiler's contraction choices on the trajectory-chaotic
lane, exactly the S3-A8.1/A9.7 class. Cure (same as T2 Momentum `Ext`):
`stepSurfaceImpl<bool QVOL>` with `if constexpr`, so QVOL=false folds to
the original expression trees. Both pins 2/2 green after the fix. Rule
confirmed a third time: runtime flags do not enter golden-pinned kernels,
even value-identical ones.

**C3 adjudications (measured on sx1, each change moving a measured
residual class):**
1. *Coupled deposits* enter the flux volume as the coupler's applied
   `gain` field (post-clamp/post-hold), not the qss re-derivation, whose
   clamp/hold approximation the anchor used to absorb (no-tide sx1
   variant: uniform to 7.4e-9 with gain, 2.1e-1 with qss). A new
   `bounce` share field pairs the Corrector's direct deposits with
   top-cell salt.
2. *Prescribed-stage cells*: flux volume = dept·A (identity row, not a
   face balance; the reservoir's routing is otherwise unbooked —
   anchor 2.6 psu·m³ → 1.3e-9 on sx1) and the stencil limiter is
   skipped there (it minted salt onto the reservoir's fresh volume:
   adjust 2.5 → 3.2e-2).
3. *Pass-through fronts*: a start-of-step-dry cell the solve fills AND
   drains in one step hands its stale-zero upwind faces the same-axis
   far-upwind value, and infiltrates at its inflow-mix concentration
   (the J1-breaking seed: raw 60.3 psu at the swash front → 35.00005;
   subsurface uniformity 2.2e-3 → 1.9e-4).

**C4 — the donor cap is REJECTED; the overdraft parks paired.** The plan's
C4 (donor outflow capped at available mass) was implemented and measured:
it starves the real receiving cells of their water's salt — sx1's J1
breach grew from 6.5e-2 to 8.2 psu. This is the C-SX contingency firing:
redesigned. An overdrawn cell (solve drained it below bed; the below-bed
lift restored dept from minted, salt-free volume; its face volumes net
negative at park time) parks its lifted water PAIRED at its own limited
concentration; the ledger difference stays on the negative-residual path
into `surf_adjust` — now a BOOKED C4 column, not an unbooked class.
Measured on sx1: 4.3e-2 psu·m³ = 3.2e-5 of M over 2 h; the sx1/sx4 gates
bound it at 1e-4·M.

**Gate reworks at first green (all recorded in `stage3-sx.yaml` with
measurements):**
- sx2's flat split-pond was VACUOUS — the exfiltration never cleared
  min_depth, the surface stayed dry for the whole horizon and no face
  ever flowed. Redesigned: 4 mm shelf/pit tilt, an 8 mm 10-psu pond
  draining through the connecting face over the et8d exfiltration
  (dt 2 s: the pond's gravity-wave CFL at dt 60 was 17 and the run blew
  up — the et8d original was surface-dry and never felt it). Anti-vacuity
  moved to the pit's salt_store (its water parks below min_depth before
  it can read as concentration; measured 0.049 psu·m³).
- Total-salt checks are LEDGER-IDENTITY checks (total − Σ booked
  cumulative columns), not raw constancy: sx2's bottom head feeds 10-psu
  groundwater in through `subs_boundary` (3.6e-2 psu·m³, booked;
  identity residual 6.4e-14 of M).
- sx5 needed a `tide-temperature` stage BC (salinity had `tide-salt`;
  without it the stage cells sink the heat toward 0 — max |T−20| = 20).
- J1 tolerance: 1e-6 abs (measured 1.9e-7 — overdraft pairing +
  production-rtol residuals on swash films; strict-solve no-tide variant
  7.4e-9).
- sx3's anchor bound is 1e-7·M (measured 4.7e-9·M: the production-solver
  residual accumulated over 14 400 steps — C5's solve-residual class;
  the strict-solver sx1 measures 9.2e-13·M). The sx-rank 1e-12-bed twin's
  scalar control binds at 1e-7 (measured 1.28e-9: superbee amplifies the
  bed perturbation ~1e3 through the drying front; a chaotic lane is O(1),
  so the control still discriminates).

**W10.4/W10.5 adjudications.** Under the conservative override the T5–T9
lanes pass except the two BITWISE sub-checks pinned to `consistent`
trajectories (ws3's a-ET hash, et9's noET 1 271-dataset pin) — they fail
by design under a scheme-changing override and stay green without it.
Decisive before→after: ws1n anchor 2.4e-2·M → 4.7e-9·M; et8d-strip
stranding 24 % → 0 exactly; b6 ss/td PASS with tolerances unchanged
(td tidal salt mass +9.89 % vs the 10 % cap — moved 1.1 points toward
the cap, recorded); the v30 H2 one-cycle smoke holds 35.000 psu
everywhere (surface and subsurface) with anchor 5.2e-9·M.

**W10.6 (measured; owner item 7b).** With `consistent` resolving to
`conservative` in a comparison build: 12/13 `default_bitwise` lanes stay
bitwise — only `coupled-salt` (the one lane whose config uses
`consistent`) changes — and the et9 noET pin breaks, as it must. The
flip is clean for default-token cases but re-pins every
`consistent`-using golden lane, and the v30 matrix runs `consistent`:
reported, owner decides.

### S3-A14 (post-T9, 2026-10-31) — new phase T10: the surface transport moves salt with the water the solve moved (O16)

**Finding (T9 W9.6, S3-A13 owner item, extended).** The surface transport
does not advect with the face volumes the free-surface solve used. The
solve moves Q = dt·Asx·Ex − Sxp·Δη^{n+1} on the start-of-step face areas.
The transport instead uses Fu = u·Asx rebuilt by `updateVelocity` from the
new depths and face areas, after `applyVelocityLimiters` has zeroed some
faces (§12B diagnosis). This is a sequencing defect, as the owner put it:
the water is moved first, then the faces are limited, then the salt is
moved. Reordering alone would not cure it, because the rebuilt flow rates
also differ in geometry.

Measured on T9: `surf_anchor` +28.4 psu·m³ on ws1n, which is 2.4 % of the
surface salt by 0.64 T and identical in store and discard. At drying faces
the salt is stranded: −0.077 psu·m³ on ws1n, and 24 % of the exfiltrated
salt on the T5 1 × 2 et8d strip. Every Stage-3 mode carries it, as did
v2.0's `consistent` update.

**Decision (owner, 2026-10-31) → O16.** Salt and heat move with the water
the solve moved. This is implemented as phase T10 rather than an amendment
to T9, because it changes the surface transport scheme. It needs its own
red evidence, gates and a §2.7-style decision on whether `consistent`
should resolve to it. It ships behind a new token,
`transport.surface_update: conservative`, so every existing gate and the
default-bitwise lanes stay bitwise until that decision.

### S3-A13 (T9, 2026-10-31) — §12A: first-green adjudications of the ws gates, two identity columns, and the et8d redesign

**W9.1 — the per-row identity needs two more instantaneous columns.**
The plan's identities hold only at mean-stage instants as written (S3-A12's
"misses 1.2–3.2 % at MSL-stage instants"). Closing them at EVERY row (ws2)
needs:
1. `stage_volume` — the water on prescribed-stage cells. `boundary_outflow`
   books the flux INTO those cells (external reservoirs), so their storage is
   outside the identity's domain; on et9-noet it rides the tide by 3.8–13.4
   m³ and, unbooked, misses the per-row identity by 3.2× the exchange scale.
2. `film_pending` — films 0 < η − z_b ≤ `min_depth` still in η (reported
   dry by `updateDepth`, removed by P2 at the next solve: the R8 mismatch).
   Identically 0 under the paired store after the sweep.
With both columns (and `film_removed`), the et9-noet discard identity
closes at all 5 401 rows to 7.1e-12 of the exchange scale (strict solve);
without the film terms it misses by 5.5e-3. Prescribed-stage cells are
never booked as film removals (the reservoir).

**W9.1 attribution.** The S3-A12 zoom-A volume (1.118e-3 m³ exchanged onto
dry-to-dry cells) is P2's: the vents land as 1.0–1.7e-5 m films after the
solve, `updateDepth` reports them dry, and P2 deletes them at the next
solve. 88 % of that volume is removed within the same 1-s interval or the
next, and the rest leaves laterally through faces. P1 is 1.2e-11 m³
(solver-noise films) and P3 is 0 (`records/T9/attribution.md`).

**ws1 authored** (`benchmarks/ws-water-salt-store/ws1-source-bound.yaml`, 64
× 16, 7 200 s, generator `input/gen_inputs.py`). The S3-A12 ingredients are
compressed: a 12.5 % beach, the H2 log-K window, the external-closure tide,
density coupling, no ET. Red on the T8 binary: surface 494.4 psu,
subsurface 35.006 psu, 36 of 121 frames over 35 (the discard control holds
35.000). C-WS is not needed. ws1n is the owner's repro verbatim (inputs in
`ws1n/input`); red 72 854 psu surface, 9 039 psu subsurface, 47 CFL warnings
(max 80.95). Green on T9: both lanes at 35.000 psu everywhere (surface,
subsurface, monitors).

**ws1/ws1n criterion (2): "zero CFL warnings" becomes a cap of 1.5.** R6
returns a store of ≥ `wetting_face_depth` onto a dry slope cell. Its first
step runs on the lagged point-implicit drag of a dry cell, the same one-step
wetting burst any dry-cell wetting produces: the b5 rain record logs 23 884
such warnings (max 83.8) from rain alone. Measured on T9: max 1.247 (ws1,
36 warnings) and 1.248 (ws1n, 55 warnings). The first event was traced at
1-s resolution (t = 592 s: V_s 2.48e-4 m³ returned, face velocity 1.24 m/s
for one step). The T8 red's 47 warnings reached 81. "Zero clip warnings" is
kept, but the model has no clip-warning class (no `log::warn` on clipping);
the bound itself is asserted by criterion (1).

**ws1 anti-vacuity floors** are 1/10 of the authoring measurement: max ΣV_s
8.80e-4 m³, 44 release rows, 1.83e-3 m³ of P4 thin vents (the vent path
fires). P1 is recorded (8.6e-11 m³). ws1n keeps the plan's floors (1e-4 m³
parked, measured 1.09e-3; 100 release rows, measured 158).

**ws2 tolerances.** Water: `max(1e-10 X, 5e-11 m³)`. The absolute floor is
the coupled lanes' solve-residual class: the ws1 store and discard twins
jump by the same ~2e-11 m³ at the same rows (worst 3.5e-11 m³), which
exceeds 1e-10 of ws1's small exchange scale (0.25 m³). Measured ws1n
8.9e-13 relative. Salt: 1e-7 of M, the et9 reduced-sum class (ws1 2.0e-9,
ws1n 5.6e-10). The film_removed-dropped negative fails on the discard twin,
as required.

**ws3(a) bitwise** is read as "every field dataset bitwise and every audit
table's shared columns bitwise". The store lane's audits append
`store_volume`/`store_release`/`surf_store`. (a-ET) is bitwise against its
T8-binary hash (`reference-hashes/stage3-ws3aet-T8.json`), with the
mass_audit table excepted per ws5.

**ws3(c) closed form gains the concurrent evaporation step.** The ET
channel keeps evaporating during the rain step, so the first wet step has
c = M0 / (V_s + (h_r − E dt) A). The plan's form assumed none. Measured
exact to 1e-15.

**ws4 lanes.** P1 rides ws1's books (a solve wetting an isolated cell needs
flow) and P4 rides ws1's dry seepage face. The P2 lanes (a draining 10 %
slope) pair at the consistent update's uniformity class: M_s/(35 V_s) is
within 1.3e-9 (bound 1e-8). The P3 lanes are 1 × 1: on a 3 × 1 flat basin,
evaporation drives the legacy ghost-stage edge velocity, booked in
`surf_boundary`, which blurs the 1e-12 forms (measured 6.7e-10, identical
on the T8 binary). R9 on P3T is exact (parked heat = 20 V_s). On the P2T
slope it is within 0.87 %, because the temperature instance's draining film
carries the legacy lagged-update concentration (20.11 °C in a wet cell);
bound 1e-2.

**ws-rank adjudicated (S3-A6 doctrine).** ws1 is trajectory-chaotic in BOTH
modes. A 1e-12 bed perturbation moves η by 6.8e-4 and flips a cell's
surface salinity by 35 psu, and the discard twin is as rank-sensitive (n=1
vs n=2: 6.5e-5, 35 psu). So:
- (A) the identity classes gate on the trajectory-stable ws4 P2 slope at
  mpi_ny = 1/2/4, with its 1e-12 twin as the stability control (1e-12
  class). Results: η 7.3e-15, concentration 1.1e-13, per-cell water_store
  1.4e-15 m³, zero parking-timing flips. The store totals agree to 5.4e-11
  relative. The plan's 1e-12 is below the strict solve's rounding at the
  ~1e-5 m³ parked scale, so `store_sum_rel` is 1e-9.
- (B) ws1 at n = 2 and 4 gates the rank-invariant guarantees (source bound,
  I2, both ws2 ledgers). All green.

**ws5.** `default_bitwise` allows `monitor/mass_audit`, and its T8 prefix is
pinned bitwise per case. `reference-hashes/stage3-T8-massaudit.json` (12
cases including et9-noet) was recorded with `--record-massaudit` on the
frozen T8 binary. The et9 noET lane compares the same prefix.

**R4 bounded (found by W9.6 on et8d).** The first implementation kept the
whole unpaired remainder whenever evaporation had taken any water that
step. On et8d's 1 × 2 strip, a superbee/limiter dry-front remainder then
became a 21 psu store (source 10 psu) and re-released at 11.6 psu, which the
limiters clipped. R4 now keeps (a) at most the salt the evaporated water
carried, c V_evap at the limited concentration, and (b) at most this step's
held-exchange credit while the coupler holds the water. Everything else
follows the discard path, as R3 states. In a true evaporative dry-out
c V_evap equals the whole remainder, so the ws3 closed forms are unchanged.

**et8d redesigned to 1 × 1 (W9.6).** On the T5 1 × 2 strip the applied
pond split across both cells below `min_depth` (5 mm). At T8 the wetting
limiter deleted it every step, unbooked: surface salt appeared in 0 of 241
frames, the store grew to 0.449 psu·m³ salt-only (the S3-A12 class), and
the lane passed `surf_adjust == 0` vacuously. On one column the D9 lane is
exercised: the pond persists in 220 of 241 frames, the store releases into
it, `surf_adjust` == 0 exactly and `subs_adjust` 5.9e-15, identically on the
T8 and T9 binaries. et8b's s·h = 35 h0 identity is unchanged: its first wet
step mixes the parked V_r = 1.79e-6 m³ (S3-A12's predicted closed-form
change shows in the concentration, not in the asserted product).

**Owner item — a transport/solve mismatch at drying faces (pre-existing,
exposed by T9).** When the sending cell of a face ends the step dry, the
velocity limiter zeroes that face's flow rate after the free-surface solve
has already moved water across it (`dept < wetting_face_depth`). The
transport's advective flux is then zero, so the moved water arrives
salt-free and its salt stays in the drying cell, where R3's cap books the
excess as `surf_adjust` (the discard path). On the 1 × 2 et8d strip with
`min_depth` 5 mm that is 0.109 psu·m³ (24 % of the exfiltrated salt), plus
1.4e-2 of downstream `subs_adjust` from the salt-free water infiltrating.
At the study's 1e-5 m `min_depth` it is film-scale: ws1n's final
`surf_adjust` is −0.077 psu·m³ against −0.46 under discard. T8 had the
identical mismatch but deleted (discard) or banked salt-only (store) every
film it touched. A fix routes the transport's face volumes from the solve's
continuity, which is a scheme change outside T9's scope.

### S3-A12 (post-T8, 2026-10-31) — §9 W5.7 / O4 / §16 C-ST: the salt-only store breaks coupled tidal runs without ET; paired water–salt store, new phase T9

**Finding (outside the gates; the 2026Frehg2 v30 study, 2026-10-30).** A coupled
2-D tidal beach without ET (1 × 160 × 40, 12.5 % beach plus 1 % terrace, stage BC
`closure: external`, M2 tide 1 m, 35 psu, consistent updates, density coupling
with the `density` gradient, `reallocation: legacy` + `redistribute`, subcycled
0.5 s / ≤ 5 s, `min_depth` 1e-5 m, 25-type heterogeneous K) runs clean under
`discard` and blows up under the shipped default `store`. The two lanes differ
in that one key only; the packaged lanes reproduce the original study runs
bitwise. Repro and analysis scripts:
`/Users/zhili/Documents/Work/个人/论文/2026Frehg2/outputs/bugreport_dryout_store_v1/` (read-only here; W9.0 copies the
inputs into `stage3/evidence/dryout-store/`).

| 0 → 28 800 s (0.64 T) | `store` (T8 default) | `discard` |
|---|---|---|
| CFL warnings (first / max) | 47 (26 300 s / 81) | 0 |
| max surface / subsurface salinity | **72 850 / 9 039 psu** | 35.00 / 35.00 psu |
| first surface value > 35 psu | 24 060 s (film 1e-5 m) | never |
| `surf_anchor` / `surf_adjust` at the end | −1.57e7 / +1.58e7 psu·m³ | +28.2 / −0.46 psu·m³ |

The homogeneous-K variant does not blow up, but its films reach 70.6 psu (no
wind) and 92.7 psu (onshore wind). With Manning 0.060 instead of 0.025 the H2
case still reaches 4 390 psu.

**Mechanism (measured at 1 s, cell y = 31.75 m, 23 400–24 600 s).**
1. The cell is dry in 91 % of frames (176 wet↔dry switches). Its salt store
   grows 57× (3.1e-4 → 1.76e-2 psu·m³; a 35-psu film at `min_depth` on the
   0.5 m² cell holds 1.75e-4) in 69 increments, all while η − z_b = 0 exactly
   at both ends of the 1-s interval, with no decrement.
2. The growth equals the top-cell salinity (19–21 psu) times the coupled
   exchange booked on the cell (`seepage_rate`, which includes vent deposits):
   8.777e-3 vs 8.771e-3 psu·m³, per-interval correlation 1.000. The exchange
   arrives as reallocation vents onto the dry seepage-face cell,
   1.0–1.7e-5 m each, just above `min_depth`.
3. The vented water leaves η through a removal **no audit column books**:
   1.1e-3 m³ booked onto cells that read dry before and after vanished over
   the 1 200 s, while `mass_audit.clamped` moved 1.1e-11 m³. The two
   `cflLimiter` rules (`WetDry.cpp`: the isolated-cell rule
   `eta(j, i) = bottom(j, i); cflActive(j, i) = 1.0;` and the small-depth rule
   `if (diff > 0.0 && diff < minDepth)`) are the surface removals with no
   audit term. Code reading gives two ways the salt is kept without its water,
   both left for W9.1 to attribute:
   - the consistent exchange credits vent salt with no depth guard
     (`rate += vent(j, i) / (cellArea * dt) * sKp(j, i);`);
   - the store-in predicate `!wet` (`vflux ≤ 0 || dept ≤ 0`) can fire on a
     cell whose water is still in η.
4. Each rewetting releases the whole store into the first film:
   - zoom A: films of 1.0–1.9e-5 m at 3 200–3 430 psu;
   - zoom B (26 100–26 700 s): 5 900–13 700 psu, and 727 000 psu at
     y = 33.25 m, t = 26 301 s. `surf_anchor` jumps +150 psu·m³ in one step at
     26 301.5 s, and the first CFL warning is at 26 300 s.
5. Density coupling then carries the salt into the soil.

**Why the §16 C-ST candidate is not enough.** Releasing in proportion to depth
up to a release depth spreads the release in time, but the store still holds
salt without water. The zoom-B store (0.101 psu·m³ on 0.5 m²) released into
h = `wetting_face_depth` (5e-4 m) gives 406 psu, still 11.6 × the source. The
defect is the missing water, not the release rate.

**Independent ledger defect (all modes).** The coupled surface identity of the
v30 one-cycle smokes misses 1.2 % (H0) to 3.2 % (H2) of the exchanged volume at
MSL-stage instants. This is consistent with item 3 and affects `discard` as
well: T9 W9.1 books it.

**Why T8 did not see it.** The W8.2 lanes (b6 ss/td, the et gates) have no
reallocation vents onto a dry seepage face, and no lane bounds surface salinity
by the source maximum without ET.

**Owner decisions (2026-10-31) → O15.**
- Thin-layer removals park water and salt together, cell by cell.
- Rain and evaporation must not be affected by the store.
- Parked water is inert, under ET too.
- Implemented as a new phase with its own gates: T9, §12A.

O4's salt-only store survives as the V_s = 0 case: the salt left by complete
evaporative dry-out.

**Consequences for existing gates (adjudicated in T9, W9.6, by amendment with data):**
- et8b's first-wet-step closed form gains the parked sub-`min_depth` remainder:
  s = 35·h₀·A/(V_r + h_rain·A), with V_r ≤ `min_depth`·A.
- et8d, et9 and et-x are re-measured.
- `default_bitwise` and the et9 noET pin lane pin `discard` and change only by
  the appended `mass_audit` column (§2.6 exception 1).

**Open (owner, §17-6).** Whether the shipped default reverts to `discard`
through §2.7 until W9.7 passes. The study's v30 matrix pins `discard` meanwhile.

### S3-A11 (T8, 2026-10-30) — §2.7/§12: the default-flip evaluations, and the flip-survival defects they caught

**W8.1 (reallocation, O9): the flip FAILS — the default stays `legacy`.**
`build-s3-df` with `Reallocation::Off`: 14 tests fail with unchanged
tolerances (b2, b3, g1.b{2,3}.{amg,gamg}, gw_orient,
b5_rank_invariance.default, b5_restart, b6_gw_smoke timeout,
b6_subcycled, b6.ss, b6.td) — the T1 preview's prediction
(`records/T1/df-preview.md`) confirmed under the binding protocol.
Decisive physics: b6's interface MAE vs the EXPERIMENT is 13× worse
under `disabled` (ss 0.329 m vs the 0.025 m golden; td 0.307 m vs
0.043 m), tidally averaged salt mass +140 %/+213 %, b3's mass balance
opens to 12.47 % of inflow (bound 0.5 %), and the dtg controller stalls
(b6_gw_smoke timeout — the T1-preview ~4× slowdown). Per-gate numbers
and diagnoses: `stage3/df/realloc-report.md`. `disabled` stays the
documented vadose-domain option (exactly conservative there, rl3/rl5).

**W8.2 (salt store, O4): the flip PASSES and is KEPT —
`transport.surface_dryout` now defaults to `store`.** The same suite on
`SurfaceDryout::Store` (reallocation reverted): 60/60 regression, b6
ss/td (the decisive ebb-drying cases) inside their golden bands, perf
green, 137/138 unit/mpi/stage3 — the single failure was et9's noET
bitwise lane missing its §2.6 pin (below), green again once pinned. Per
§2.7-5 no combined run is owed (it is conditioned on BOTH flips passing
individually). W8.3 updates: the parameters.md default column, the
theory section, `ConfigTest.SurfaceDryoutDefaultIsStoreAfterTheT8Flip`
(the shipped-default pin + the explicit-discard round trip), and the
ET + salinity store rule now fires only on an EXPLICIT `discard` (an
absent key resolves to store).

**The flip-survival defects (all fixed in this phase):**

1. `BITWISE_PINS` — the §2.6 "v2.0-equivalent pins injected so the rule
   survives the T8 default flips" — was never populated by T1/T5. Under
   the flipped binaries, `regression.s3.default_bitwise` (b2,
   g6c-coupled) and et9's noET pin lane went non-bitwise: the gates were
   testing "the shipped default", not "the v2.0 configuration". Fixed:
   the table now carries `groundwater.reallocation: legacy` and
   `transport.surface_dryout: discard` (module-gated, printed per lane),
   and `et9-noet.yaml` pins `surface_dryout: discard` in-file (it is
   hashed against T5-start references).
2. Four unit tests escaped the W1.6 explicit-pin rule and asserted
   legacy reallocation semantics through the shipped default
   (CoupledModule.SeepageHeldOnDrySurfaceUntilMinDepth,
   CoupledModule.ReallocationVentDepositsOnWetSurface,
   TransportModule.DensityReallocationGradientKeepsInfiltratedSalt,
   TransportModule.ConsistentSubsurfaceUpdateClosesTheSubcycledLedger).
   Now pinned `Reallocation::Legacy` at their config sites.
3. The `Config.cpp` parse sites carried their own hardcoded
   `"legacy"`/`"discard"` fallbacks, so §2.7's "change the single
   constant in Config.hpp" would silently never reach a run. The
   absent-key path now keeps the struct initializer (verified with
   `--resolve` on both DF binaries); the r1 resolved-config echo is what
   proved the flip reached every run (§2.7-6).

**Release-pipeline scope (W8.6/W8.9):** `run_sanitizers.sh` and
`ci_release_gate.sh` predated Stage 3 and ran no stage3 label; both now
include `stage3`/`stage3_local` (and the release gate `stage3_nightly`)
— the §8.3 lesson (v1 wind: implemented, never exercised, shipped
unvalidated) applied to the pipeline itself.

**W8.4 adjudication:** the et9 p1 lanes' first run failed on et9's noET
BITWISE pin — the T5-start hashes are recorded at 1 thread on the
default aij path, so the comparison is undefined under OMP 4 + a
backend override (threaded reductions + different solver algebra). The
pin is skipped under `--mat-type` (printed), like perf_baseline under
sanitizers; et9's tolerance-based criteria (complementarity, zone
anti-vacuity, the salt-audit identity, the substep ratio) still gate
the p1 lanes, and both backends pass them (retest
`records/T8/w89-et9-p1-retest.log`).

**W8.6 adjudication (same pin class):** the full sanitizer matrix
(`run_sanitizers.sh --full` with the new stage3 labels) came back
183/184 — the one failure again et9's noET bitwise pin: sanitizer
instrumentation changes codegen under the SAME compiler, and the
h5_bitwise toolchain fingerprint is compiler identity only. The harness
now skips bitwise-pin comparisons under `FREHG_SANITIZED_BUILD=1`
(exported by run_sanitizers.sh); et9's tolerance-based criteria pass on
the sanitized binary (identical to 4 significant figures), retest
`records/T8/w86-et9-retest.log`. No sanitizer FINDING anywhere in the
matrix — both failures of the day were comparison-definition classes,
not defects.

**W8.7:** the BC matrix gains the periodic-axis rejection cells
(`edge-bc × periodic-{x,y}-axis`, schema-rejected, probed by
regression.s3.pd6 — the y-axis probe added here; its polygon must select
interior columns only, since a full-span wall legitimately attaches
corner faces of the non-periodic axis and defeats the probe).

### S3-A10 (T7, 2026-10-29) — §11: PM first-green adjudications and the measured convention gap

1. **δ_conv measured (et0 part 4 → et10's tolerance).** On the FAO-56
   Example-19 daytime lane (14.00–15.00 h), |PM_q(model constants) −
   ET₀,FAO chain| = **0.016166 mm/h (2.58 %)** — inside the plan's
   expected 1–3 % class. `tolerances/stage3-et.yaml` records
   `delta_conv_mm_h: 0.0194` (1.2× headroom) plus
   `fao_rounding_mm_h: 0.005` (FAO prints ET₀ to 0.01 mm/h). The model's
   et10 column lands at 0.64310 mm/h, |Δ| = 0.0131.

2. **The FAO transcription's sensitivity floor.** The plan's "each
   transcribed value perturbed by 1 % fails" cannot hold where the
   published precision is coarser than 1 % of the value (e.g. the
   hourly wind term printed 0.17, the night-lane entries printed to
   0.01): a 1 % nudge sits inside the last-digit ± 1 acceptance by
   construction. The battery perturbs those keys by 2.5× the precision
   floor instead and prints the list (11 of 65 compared keys); the
   remaining 54 use the literal 1 % rule. Also: the daily chain must be
   fed u2 at FAO's full precision (Eq. 47 on 10 km/h at 10 m →
   2.0776…), not the printed 2.078 — the printed wind_factor 6.450 is
   otherwise missed by one unit in the last digit; and `et0_rounded`
   (3.9, a display rounding of 3.88) is excluded from the chain
   comparison.

3. **Commensurability is second order in the ABSOLUTE gap, first order
   relative.** For the energy-balance surface temperature T_s* the
   secant-vs-tangent algebra gives L_PM/L_exact = 1 − δ/(s+1) with
   δ = λ q″·ΔT/(2 c_pa) — a first-order RELATIVE error (measured fitted
   slope 0.98 over A ∈ {50, 100, 200, 400} W/m², ΔT 0.64–4.7 K,
   relative gaps 0.4–3.1 %). The plan's "fitted order ≥ 1.8" holds for
   the absolute gap |E_PM − E_bulk(T_s*)| (measured slope **2.03**),
   since both rates also vanish linearly in ΔT; the unit gate fits the
   absolute gap and this entry records both slopes. The test drives the
   balance from q_a = q_sat(T_air) (zero aerodynamic baseline): with
   q_a < q_sat the baseline flux can exceed A and put T_s* below T_air
   — the balance still closes but the fit's abscissa would cross zero.

4. **PmHighResistanceLimit needs r_a ≤ ~0.25.** The suppression ratio is
   (s + 1)/(s + 1 + r_s/r_a); at 25 °C (s ≈ 2.9) the plan's
   E < 1e-12·E(0) bound at r_s = 1e12 requires r_a below ~0.25 s/m. The
   unit gate runs the pure math limit at r_a = 0.1 (physical r_a values
   land at ~1e-10 — still a vanishing rate, just not 1e-12).

5. **The soil/vegetation albedo is schema-required under PM** (every
   soil type exactly when the soil channel is on; every vegetation type
   when vegetation is present) rather than silently defaulting to 0 —
   a 0-albedo PM rate overstates A by (1−α)⁻¹ ≈ 30 % at FAO's grass
   albedo and would pass no gate that could catch a forgotten key. The
   type structs carry has-flags (`hasAlbedo`), the staged per-column
   fields default to 0 only where the schema has already guaranteed PM
   is not the consumer.

6. **et10's exactness lane measured at 1.7e-16** (T_act vs
   `gate_et.pm_q_form` over the hour, K = 1e-15/Ss = 0 exactness
   conditions per et6) — bound kept at the plan's 1e-10; the albedo
   0.23 → 0.20 negative misses by 3.1e-2.

### S3-A9 (T6, 2026-10-29) — §10: first-green adjudications, and a pre-existing coupled-transport restart defect found and fixed

**The restart defect (found by the et-veg-restart lane, present since the
first coupled transport restart):** a SUBCYCLED window can take ZERO
groundwater substeps in the first step after a restart (the restored gw
clock lags the surface clock), and the transport then reads the LAST
SUBSTEP'S per-substep state — the volumetric fluxes `qx_/qy_/qzF_`, the
face conductivities `kx_/ky_/kzF_` (the limiter's conductive-stencil
gates) and `dtgLast` — which the checkpoint never carried. On a restart
these were allocation zeros: every extrema gate stayed closed, the
limiter windows kept their seeds (lo = +bound, hi = 0 for salinity) and
the first transport step clipped the ENTIRE subsurface scalar to zero
(measured: subs_adjust = −M in one step; the et9 control reproduced it
with vegetation absent, with the ET block absent, under `legacy` AND
`consistent` updates, periodic or walled — a defect of the base coupled
transport restart). b6-restart never caught it because its dtg ≤ dt
guarantees ≥ 1 substep per window. **Fix:** coupled/gw transport
checkpoints now carry `qx_last`, `qy_last`, `qzf_last`, `kx_last`,
`ky_last`, `kzf_last` and the scalar `dtg_last` (restored into the gw
solver), reproducing the documented stale-flux legacy quirk bitwise;
runs without transport keep the prior layout byte-identical. The control
restarts BITWISE after the fix; all restart gates re-ran green on the
fixed binary. (Known remaining nuance, recorded not fixed: a TIME-VARYING
groundwater BC's `topValue_/sideValue_` would also be stale-from-
construction across a zero-substep restart window; every current gate
case uses constant gw BCs, and the surface BCs re-evaluate per step.)

Gate adjudications at first green:

1. **et6 asserts dz-free** (per-cell Δθ_k = uptake_k·t against the
   `root_uptake` output; the audit total; the uniform lane's shared
   T_pot/L sink value with the zone-cutting cell excluded at a strict
   98 % filter — the regular mesh trims partial top cells so no
   raster-side dz ladder is reliable) with **Ss = 0** (with Ss > 0 part
   of the sink is funded from compressible storage and the θ identity is
   not exact) and the K = 1e-15 class bound 1e-9 (measured 3.1e-10; the
   T5 et8a precedent — K = 0 is outside the vG solver).
2. **An all-bare vegetation map is legal** (warn + inert module): the
   et5(i) pin REQUIRES running the block with no vegetated cell and
   asserting bitwise against vegetation-absent (measured: 3 groups
   bitwise).
3. **et5's bare-contrast lane runs 2 h, not the day**: on the full
   horizon the silt redistribution smears the top-face removal over the
   whole column (measured: a uniform θ_r end state, vloss 0.11 — the dry
   singularity again), erasing the locality the lane contrasts. At 2 h
   the share is 0.53 (≥ the 0.5 criterion) vs the canopy lane's 0.051
   (≤ 2·b_top = 0.10) on its own full-day horizon.
4. **et8-t's proportionality bound is 1 − 1e-4** (measured correlation
   0.9999931): the multiplicative concentration factor is exp-nonlinear
   in the uptake (c1/c0 = θ0/θ1), so the linear correlation carries a
   quadratic linearization deviation; the exact statement is the
   per-cell c·θ identity (measured 3.0e-8 at the tiny-K class).
5. The audit column indices of the T lanes depend on which channels are
   configured (the columns are appended in soil-then-transpiration
   order); the gates read them per configuration.
6. The coupled audit aggregation (`Coupler::step`) now carries
   trDemand/trApplied like the T5 et terms — caught by et9-veg's
   anti-vacuity reading T = 0 from the audit while the rate outputs
   showed live transpiration.
7. **The W6.4 sink lives in SEPARATE post-kernels**, not inside the
   predictor-assembly / θ-update kernels: inserting even a guarded
   branch into those pinned kernels perturbed the compiler's fma
   contraction and broke `default_bitwise` on b2 (the S3-A8.1 class,
   now twice measured). The predictor adds the sink to the assembled
   RHS vector afterwards; the corrector subtracts (dtg·rRho·S·V)/coeff
   from θ afterwards (one extra rounding vs the in-kernel form — no
   golden pins the vegetated path, and the W6.4 unit gate's 1e-12/1e-14
   hold) and accumulates the ss-audit delta. The default path is
   byte-identical by construction and `default_bitwise` re-verified.

### S3-A8 (T5, 2026-10-28/29) — §9: first-green adjudications and case redesigns of the et gates

The capability landed as specified (D1–D9 all removed on the ET path;
`default_bitwise` green and the et9 noET lane BITWISE against its T5-start
pin — 1271 datasets — so every T5 edit is default-neutral by execution).
The adjudications, each with its measured datum in `stage3-et.yaml` or the
case header:

1. **et1 lane a is fma-contraction class, not bitwise.** The new and v2.0
   chains are the same expression tree, and the first 15 substeps ARE
   bitwise; the build's `-ffp-contract=fast` then contracts the two
   kernels differently (±1 ulp in the rate, first measured at substep 16:
   θ 5.6e-17). Through the vG front it peaks transiently (θ 1.3e-9 at 4 h)
   and decays (1.8e-11 at day end). Bounds 1e-7 (≥25× the peak).
2. **et1 lane b (the D1 red) is a COUPLED lane** (`et1b-coupled.yaml`):
   the v2.0 moisture-guard pair exists only in the coupled mode-0
   branches, so the plan's uncoupled saturated lane cannot reach D1.
   Measured: the guard-skipped predictor leaves a 5872 m first-output head
   difference against the guard-free path (the Ss response to the
   unseen sink), and the new path's ledger closes to 1.7e-13.
3. **et4 horizon 3 d → 12 h.** The plan's horizon drove the silt column
   into the vG dry singularity (θ → θ_r everywhere; measured: Ss phantom
   release 0.148 m³, vloss 0.038 m³ — the books stay exact but the state
   is outside the Richards domain, the S3-A4 applicability class; the
   legacy guarded path reaches the same state). The 12 h lane holds
   ΔS = boundary_in − ss + realloc − vloss at 4.2e-15 with vloss ≡ 0;
   the cap anti-vacuity lives in et7 (books exact at 1e-14).
4. **et8a pins the tiny-K limit, not K = 0** (K = 0 is outside the vG
   solver: 1e29 blowup). K = 1e-12 measures a 6.7e-7 per-cell
   redistribution class over the 2 h horizon — per-cell bound 1e-5,
   while the TOTAL mass and adjust/anchor hold 1e-12 (measured 0 / 7e-17
   / 6e-16). Ss and horizon trimmed to keep the top cell off the dry
   singularity.
5. **et8d redesigned twice**: both side-head authorings produced no
   measurable exfiltration; the working driver is a bottom pressure head
   0.2 m above hydrostatic (K = 1e-6, ~1e-4 m/step, holds ~40 steps under
   min_depth 5e-3). The conservation assertion is the FULL audit identity
   (measured 4.2e-15) plus `surf_adjust == 0` EXACTLY (v2.0's deletion
   channel never fires — the D9 observable); the bottom influx
   (subs_boundary 0.047) and the θ-basis re-anchor (subs_anchor 0.040)
   are booked data. A robustness fix fell out: the store captures only
   POSITIVE dry-cell mass (a negative advective residual against a dry
   front stays in surf_adjust as in legacy; first measured as a −19
   psu·m³ store on the violent first authoring).
6. **et3's criteria (1)/(3) read the audit columns**, not reconstructed
   state (the saturated column funds the soil flux from Ss storage,
   invisible in θ), and criterion (3) is asserted through the LAST output
   where the surface-extrapolated top moisture keeps α₁ = 1 (the plan's
   own premise; the bulk lane's top cell crosses the 0.375 mark near the
   end — measured rate decay 1.6514e-6 → 1.6413e-6). Measured at the
   asserted points: pond phase 2.5e-15/3.1e-15, combined rate 0 (bulk,
   t ≤ 700) and 2.9e-15 (prescribed, full horizon), post-switch
   continuity ratio 1.000000, the v2.0 window-lag negative 0.26 (a 74 %
   deficit — the double-counted supply empties the pond early) vs the
   1e-3 floor.
7. **et9's salt assertion is the full audit identity at 1e-7** (measured
   2.7e-9 under the default solver over 5400 steps — the reduced-sum
   rounding class; the strict et8 lanes hold the same identity at
   4e-15). The substep-ratio criterion reads the run records' solver
   `solves` counts.
8. `evapotranspiration.potential: penman_monteith` fails validation with
   the plan's "not available before phase T7" message (removed in T7).
9. **The et-x identity lanes run on a trajectory-stable base** (the S3-A6
   twin-control precedent). The first orient base (eta −0.20, a drying
   thin rim) was front-chaotic — the 1e-12 twin reproduced its rank
   deltas exactly (tracer 1.4e2 of a 35-psu field, wet_time 10 s) — so
   the strict per-field identities moved to a stabilized base (eta
   −0.05: every shoreline cell ≥ 0.01 m from the threshold while the
   600 s ET removal is 6e-5 m). Results: rank invariance at MACHINE
   precision (θ 1.3e-15, tracer 1.1e-13, store 5e-20 over n = 1/2/4);
   the battery's transpose lanes exact (tracer 0.0, store 8e-20, θ
   9e-15). The battery's six MIXED lanes are RECORDED, not gated: the
   documented legacy live-east/north-edge-face asymmetry
   (symmetry-exemptions.md) lets the rotated pond drain its shoreline
   row — measured on a provably static base (the id lane moves 1.5e-4 =
   the ET removal; the rot90 lane 5.9e-2 = one shoreline depth; the ICs
   mirror to 0.0) — so shoreline wet/dry status (tracer 20, θ 0.18) is a
   pre-existing surface property whose own gate is the swe-orient
   battery; T5's symmetry assertion is the exact transpose class. et9
   keeps the wet/dry-front physics under its trajectory-robust
   assertions (ledgers, zones, complementarity, the noET bitwise pin).
10. **The YAML token is `discard`, not `delete`** (the S3-A3 class): the
    T0 forbidden-word rule bans the literal `delete` anywhere in `src/`
    (raw-memory discipline), which a string enum token would trip. The
    semantics, default, and the T8 DF flip are unchanged; `delete` is
    rejected by the schema like any unknown token.

### S3-A7 (T4, 2026-10-28) — §8 W4.1: the registry's Field3 leg is deferred to its first consumer

`io::IntervalRates` lands with the Field2 (m/s) leg only. W4.1 specifies a
Field2/Field3 registry, but no T4 or T5 work item reads a 3-D rate — the
first consumer is T6's `root_uptake` (W6.6) — and authored-unexercised code
is this project's documented failure class (V2-A9..A12, the CUDA lane):
four of Stage-2's five process failures were code or gates that existed
without ever being executed. The Field3 (1/s, per cell volume) leg is
added in T6 alongside the gate that exercises it. Recorded in the
IntervalRates.hpp header and `adding-an-output-variable.md`. Everything
else in T4 landed as specified with no tolerance adjudications: dg1's
identity measured 6.7e-16 against the authored 1e-12, the dg2 subset and
schedule checks are exact, and dg-restart is bitwise as authored (the dg1
case has no superbee/wind branch noise, so the S3-A6 surface-restart
relaxation was not needed here). dg2 gained a third lane beyond the plan's
two (`monitor_interval: 77` + checkpoint 1800: 77 divides neither 1800 nor
3600, so the checkpoint- and t_end-forced row paths are exercised by
execution, not by coincidence with the 60 s multiples).

### S3-A6 (T3, 2026-10-28) — §7: pd2-surface, pd3-surface and pd-restart hold the scheme's own class, not bitwise

The periodic conversion itself is clean: the branch-free lanes hold the
authored strict bounds (pd2-gw roll identity measured 1e-13 vs the 1e-11
bound; pd3-gw rank invariance 5e-13 at n = 2/4 vs 1e-12; pd-restart flow
fields measured 1.9e-15 **absolute** — ulp class). The surface lane's
tracer/velocity deviations are the pre-existing scheme's exponential branch
amplification of solver rounding, proven by a **walled control**: the same
pd2-surface case with `periodic` removed and x walls added, run on the same
binary through the same harness lanes, measures deviations STRICTLY LARGER
than the periodic lane on every field —

- restart pair: flow abs 1.3e-13 vs periodic 1.9e-15; tracer rel 6.6e-6
  vs 4.0e-6;
- rank study n = 1/2/4 (strict jacobi+Krylov still reduces dot products
  across ranks, so n changes summation order): uu 6.6e-7 vs 7.0e-8,
  tracer 1.4e-5 vs 3.7e-6, eta 1.7e-9 vs 1.6e-10 — ~10× worse walled.

Two further diagnostics localize the mechanism: the pd2 roll deltas are
spread over ALL 12 columns (uu per-column max 2e-8–7e-8, maximum NOT at the
seam columns), and they grow exponentially (uu 3.2e-10 @ 60 s → 2.2e-9 @
120 → 2.1e-8 @ 300 → 7.3e-8 @ 600 s) — the signature of superbee/wet-dry/
relative-wind branch flips seeded by the 1-ulp solver summation-order
change that rolling the unknown ordering (or restarting the Krylov state,
or changing the rank count) produces. Adjudication:

- **pd2** splits the roll bound per field class: the gw lane keeps the
  authored `roll_abs` 1e-11; the surface flow fields gate at
  `roll_surface_abs` 1e-6 (measured max 7.3e-8, 14× headroom, 4 orders
  below a defect-class O(1e-2) seam error) and the tracer at
  `roll_tracer_abs` 1e-4 (measured 7.3e-6 on scale 7.5). The run length
  stays 600 s — it is what carries the seam-crossing anti-vacuity floor.
- **pd3** keeps `rank_abs` 1e-12 on the gw lane and adopts the same field
  classes on the surface lane (`rank_surface_abs` 1e-6, measured 7.0e-8;
  `rank_tracer_abs` 1e-4, measured 3.7e-6); the threshold logic is the
  pd2 check, battery-covered there.
- **pd-restart** adopts the project's restart standard (the wind_restart
  class; no existing gate holds a limited tracer bitwise across restart):
  flow fields `flow_abs` 1e-12 (measured 1.9e-15; even the walled
  control's 1.3e-13 passes, a restart-order defect fails by ≥ 8 orders)
  and the tracer `tracer_rel` 1e-4 (measured 4.0e-6).
- The battery (`scripts/test_pd_gates.py`) gains four S3-A6 lanes: the
  measured classes pass, defect-class corruptions (1e-4 seam column on
  flow, 1e-2 on tracer, 1e-9 restart-order flow error, 1e-2 rel tracer
  restart error) fail.
- **pd2-coupled** (the plan's third pd2 lane) is authored on a SATURATED
  tidal basin after two instructive failures, both diagnosed by the
  1e-12-twin control (perturb one bed cell by 1e-12 and rerun — if the
  twin reproduces the roll deltas, the trajectory, not the seam, is the
  source): (1) a marginally-drying row over a vg_n = 1.23 soil was
  chaotic with Lyapunov time ~40 s (twin and roll both plateau at eta
  3e-4, uu 1e-2, head 2.0 m); (2) a vadose gap under a permanently wet
  surface carried a singular wetting-front event (dh/dθ at the front)
  that retained O(1e-2) subsurface-tracer memory of any seed in both
  pairs. No trajectory identity can gate either regime; the vadose seam
  is pd2-gw's (clean, flat-terrain) lane, and a dedicated control proved
  terrain-following + periodic itself exact (roll deltas 6e-15 head /
  1e-20 flux on a terrain-following pd2-gw variant). The adopted lane —
  always-wet bed, saturated column, Ss = 1e-3 tidal loading, external
  closure, consistent salinity both ways — measures roll deltas within
  the field classes (uu 1.5e-7, head 2.1e-11, c_gw 2.0e-10, c_surface
  9.0e-8) with a real exchange signal (max |Δc_gw| = 0.284 psu), gated
  by the `coupled_exchange_floor` 0.05 anti-vacuity check.
- **pd8** adds the plan's surface battery (`pd-periodic/pd8-orient-base.yaml`:
  periodic-x torus bump, oblique wind rotated with the domain by the
  harness — the dihedral coord_op's linear part applied to (u10, v10) —
  tracer slab, y walls) beside the gw lane that ran first; it compares the
  SCALAR fields (eta, depth, tracer — the swe-orient convention: uu/vv are
  staggered vector components that transform with sign flips and a
  component swap, which the scalar push cannot express; a velocity
  asymmetry feeds eta within steps). `transpose` gates at the pd2 field
  classes (measured eta 6.5e-11, tracer 3.6e-7), the six mixed transforms
  at the x-orient swe class 1e-2 (the y walls' legacy minus/plus edge
  arithmetic; measured ≤ 1.2e-4).

### S3-A2 (T1, 2026-10-27) — §2.3: compiled unit gates cannot be merged red

The gate-first rule is executable verbatim for YAML-driven gates (they fail
`--validate` as "capability absent" on the pre-phase binary — the recorded red
evidence). Compiled unit tests referencing new config members cannot exist
before the members do without breaking the `-Werror` build, so each phase's
unit gates land in the same work item as the config change, with the schema
rejection of the new keys captured in the YAML gates' red logs standing as the
phase's red evidence. This is the Q4–Q6 practice (`test_atm.cpp` did not exist
before Q4), recorded here so it is a rule rather than a habit.

