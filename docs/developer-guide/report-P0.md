# Phase 0 development report — handoff for P1+

**Audience:** the developer (or agent session) implementing P1 and later
phases, starting with no memory of the P0 session. Read this together with
`FREHG2_UPGRADE_PLAN.md` (the binding spec, one directory above the repo)
and `dod-P0.md` (the verified Definition of Done). This report explains
*what exists, how it fits together, which conventions are binding, and what
was decided along the way* — so later phases don't rediscover or
accidentally contradict it.

**Status:** P0 complete. Initial commit `6c51c3f` (2026-07-18), branch
`main`, local repository at `/Users/zhili/Codes/frehg2-upgrade/frehg2`
(no remote configured). All P0 exit criteria verified green; every claim
below is enforced by a test or script named in `dod-P0.md`.

---

## 1. What exists and where

One static library `frehg::core` (sources in `src/core`, `src/bc`,
`src/io`) and the `frehg` executable (`src/main.cpp`, currently
`--validate` only — the simulation loop arrives with the physics phases).
The physics libraries `frehg::{swe,gw,transport,driver}` do **not** exist
yet; `src/CMakeLists.txt` is where P1/P2 add them.

### core/

| Component | What it does | Facts P1+ must not re-derive |
|---|---|---|
| `Types.hpp` | `real_t = double`; `Field2<T>`/`Field3<T>` = `Kokkos::View` `(j,i[,k])` LayoutRight in `MemSpace`; `HostField2/3`; `frehg::FatalError` | No `SharedSpace`/UVM anywhere; host access only via mirrors at I/O boundaries |
| `Logger` | `log::info/warn/error`; `log::fatal(msg)` logs then **throws `FatalError`** (never exits) | `log::msg() << ...` builds messages; info prints rank 0 only; `main` catches FatalError → `MPI_Abort` when size > 1. Tests assert on `EXPECT_THROW(..., frehg::FatalError)` |
| `Timer` | static registry; `start/stop(name)` auto-nests into paths (`"solve/assembly"`); `Timer::Scoped`; `report(comm)` = MPI-merged min/mean/max table | `stop` with a non-innermost name is fatal (enforces nesting) |
| `PetscSession` | RAII: MPI_Init_thread → Kokkos::initialize → PetscInitialize; reverse on destruction; optional PETSc options file | Construct exactly once, first thing in `main`/test mains. `session.comm()` is the world comm |
| `TimeSeries` | strictly-increasing `(t, v)` samples; linear interp; clamps at both ends; **stateless** (binary search — no cursor, so checkpoints need no series state) | `fromFile`: two-column `t value` text, `#` comments. Non-monotonic time, extra columns, empty file = fatal |
| `Config.hpp/.cpp` | typed config structs mirroring the v2 schema; `validateConfigFile(path)` → `ValidationResult` (collects *all* errors); `loadConfig(path)` (fatal on invalid) materializes every default; `describeConfig` for the log header | `FrehgConfig::resolvePath()` anchors file refs at the YAML's directory. `cfg.rawText` carries the config for HDF5 embedding |
| `ConfigSchema.cpp` | the declarative schema tree + walker + cross-field checks. **This file is the schema's single source of truth** | New keys ⇒ update schema here **and** `docs/user-guide/parameters.md` **and** a test, same PR (`unit.parameter_docs` fails otherwise). yaml-cpp gotcha: `operator[]` on an undefined node **throws** — use the internal `sub(node, key…)` helper for chained lookups |
| `Grid` | 2D block decomposition over (i,j) (full z-columns per rank, plan §5.2); non-divisible blocks `N/P + (p < N%P)`; auto layout via `MPI_Dims_create` (larger rank count goes to the larger extent); `buildGlobalIds(ktop)` builds compressed `gid2/gid3` (PetscInt, −1 = inactive) with halo gids exchanged; z-spacing `dz·stretch^k` | Rank map: `pi = rank % px`, `pj = rank / px`. gid order: rank blocks, then (j, i, k), k innermost — matches PETSc row ownership *and* the §7 output flattening. `ktop == nz` ⇒ column fully inactive in both 2D and 3D. `xCenter(i) = (i+0.5)dx` — polygons rasterize against these |
| `HaloExchanger` | register fields once (`add(name, field)`); `exchangeAll()` = one coalesced message per neighbor; `exchange({names})` targeted; Irecv×4 → pack kernels → fence → Isend×4 → Waitall → unpack | Constructor takes `gpuAware` bool — the **only** place that toggle exists. Host-staged path deep-copies through persistent host mirrors. Domain-edge halos are never written (BC system owns them). Tags pair `dir ^ 1` |
| `LinearSystem` | PETSc Mat/Vec/KSP with COO assembly: `setPattern(rows, cols)` once (`MatSetPreallocationCOO`; −1 indices ignored — fold boundary legs into diagonal/RHS), `setValues(view)` per solve, `solve(rhs, x)` → `SolveStats` | Defaults: KSP CG, PC bjacobi + icc(0) via options DB (respects user overrides), rtol 1e-8 / atol 1e-14. **Divergence is fatal** (legacy's silent continue was a bug, §5.1). Duplicate COO indices are summed (tested). `isSymmetric(tol)` via ‖A−Aᵀ‖∞ (MatIsSymmetric rejects MPIAIJ). Debug builds check symmetry on first solve |

### bc/

- `Polygon`: ray-cast point-in-polygon; **on-edge and on-vertex = inside**
  (fixed convention, unit-tested with convex/concave/edge/vertex cases).
- `BoundarySet`: rasterizes each configured BC at construction —
  footprint targets (surface / groundwater_top / groundwater_bottom) select
  every owned cell whose center is inside; `groundwater_side` selects only
  domain-edge cells and tags the edge face(s). Per-BC sub-communicator
  (`MPI_Comm_split`) over ranks with members. A polygon selecting zero
  cells globally is **fatal**. `BcFace` ordering documents the decoded
  legacy `bctype_GW` faces: [x+, x−, y+, y−, bottom, top], legacy codes
  0 = no-flux, 1 = head, 2 = flux, 3 = free drainage.

### io/

- `GridDataReader`: headered rasters (`ncols/nrows/xllcorner/yllcorner/
  cellsize/nodata_value`, any case, `n_cols` accepted) **or** flat lists;
  count must match exactly; nodata → NaN. **Binding row convention: first
  data row is j = 0 (file order), NOT ESRI north-up** — chosen because the
  legacy references were consumed that way. 3D fields are flat lists in
  `(j·nx+i)·nz+k` order.
- `Hdf5Output`: single parallel file, §7 layout exactly (`/surface/<var>/<t>`
  etc., `<t>` = integer seconds; NaN above ktop; units/time/long_name attrs;
  `/frehg2` root attrs incl. full config text and git SHA). Each rank writes
  its (j,i) block as one regular 1D hyperslab. `writeDistributed2/3` +
  attribute helpers are the shared raw layer Monitor/Checkpoint build on.
- `Monitor`: extendable `(time, values…)` table per probe; rows buffered on
  the owning rank; `flush()` collective.
- `Checkpoint`: `/checkpoint/<t>/…` bit-exact (0-ULP round-trip tested at
  1/2/4 ranks); header carries t, step, named scalars (e.g. `dtg`).

### tools/, benchmarks/, scripts/, .github/

- `tools/migrate_yaml_v1_to_v2.py`: v1→v2 migration; the rename table in its
  docstring is the migration contract. Unknown v1 keys are a hard error.
  Notable: `time.max_step` is **dropped** (legacy `NT` was parsed but unused).
  Tested by `tools/test_migrate_yaml.py` (ctest `unit.migrate_yaml`), which
  also asserts the committed b1/b2 configs embed the migration output
  (documented per-case exceptions where the legacy `input` file beat the
  v1 draft).
- `benchmarks/`: seven validated configs — `b1-sw`, `b2-gw`, `b3-kirkland`,
  `b4-govindaraju`, `b5-vcatchment`, `b6-kuan-{ss,td}` — with converted
  input data (provenance in YAML headers + READMEs). Goldens are *not*
  committed; the P1+ regression harness fetches them from
  `../legacy/benchmarks/`.
- `scripts/`: `check_forbidden.sh` (§11.3 gate), `check_parameter_docs.py`
  (schema↔docs lockstep; also ctest `unit.parameter_docs`),
  `ci_build_and_test.sh` (the §11.2 per-PR gate), `run_sanitizers.sh`,
  `ci_install_deps.sh` (CI-side Kokkos/PETSc builds).
- `.github/workflows/`: `build.yml` (gcc+clang matrix), `sanitize.yml`
  (ASan+UBSan), `regression-nightly.yml` (cron; regression label grows in
  P1+). **Never run in anger yet — no remote exists.** First push must
  watch them.

## 2. Binding conventions (do not change without a plan amendment)

1. **Field layout**: `(j, i[, k])`, halo width 1 in j and i, none in k,
   k unit-stride, k = 0 at the land surface increasing downward. Interior
   local indices 1..nyLocal / 1..nxLocal.
2. **Flattening** (files *and* HDF5): 2D `j·NX+i`, 3D `(j·NX+i)·NZ+k` —
   pinned by b3's `makeplot.py`; asserted element-wise in tests.
3. **Active masking**: single source `ktop(j,i)`; surface cell active iff
   `ktop < nz`; dry-but-active surface cells keep matrix rows.
4. **Error handling**: no exit/abort in library code; every PETSc/HDF5 call
   checked (`FREHG_PETSC_CHECK` / `FREHG_H5_CHECK` macros in the .cpp
   files); all failures route through `log::fatal`.
5. **Units**: SI (m, s); series and rasters in SI already (benchmark inputs
   were converted at data prep, e.g. mm/h → m/s).
6. **Forbidden patterns** (§11.3) are live: no TODO/FIXME/"for now"-style
   markers, no `\bnew\b`/`\bdelete\b` even in comments (write "fresh",
   "removes", …), `printf`/`std::cout` only in Logger.cpp/main.cpp, and the
   dropped-feature tripwire (never name identifiers after `difuwave`,
   `subgrid`, `use_mvg`, …) — run `scripts/check_forbidden.sh` early and
   often.

## 3. Decisions made in P0 that later phases inherit

- **Seven configs, not six**: b6 ships as two variants (ss/td) because §9
  gates both; they share `benchmarks/b6-kuan/input/`.
- **b6 runs `coupling.mode: sync`** although the committed legacy input says
  async — the goldens are the syncV4 runs (risk-register item; recorded in
  `benchmarks/b6-kuan/README.md`).
- **`soil.map.file` is a 3D id field** (`(j·nx+i)·nz+k`): b3's soil map is
  genuinely z-layered (file rows = k levels), b5's legacy map was a 2D slab
  replicated over nz (verified identical across all 25 slabs at data prep).
- **b2 fidelity choices**: where the archival v1 draft disagreed with the
  committed legacy `input` (dtg limits, initial moisture 0.033), the legacy
  input won; exceptions are listed in `tools/test_migrate_yaml.py`.
- **Provisional BC mappings to revisit at the gates**: b4/b5 "free outflow"
  outlets are expressed as `kind: eta` held far below the bed (−10 m) —
  SERGHEI used transmissive-outflow codes that have no direct v2 kind. The
  P1 (b4) and P3 (b5) gates decide whether this stands or the schema needs
  an outflow kind. Similarly b6's seaward groundwater head uses
  `value: {hydrostatic: {eta: 0.15}}` with a *constant* reference stage;
  whether the tidal variant needs the reference to follow the tide series
  is a P4 question.
- **b3 initial head −500 m** matches the uniform legacy `head.input`
  verbatim; whether the digitized h = 0 / h = −400 contours are meters is
  adjudicated by the P2 gate, not by the config.
- **Surface `dt` for gw-only runs** (b2: 1e-4 s legacy-faithful; b3: 10 s
  chosen) is the outer-loop step; P2 defines how the driver loops
  gw-only cases.
- **`time.max_steps` does not exist** in v2; don't reintroduce it.

## 4. Toolchain and environment (macOS arm64 dev machine)

- **Configure gcc lane** (primary):
  `cmake -B build -DCMAKE_PREFIX_PATH=/Users/zhili/Codes/local
  -DCMAKE_CXX_COMPILER=g++-15
  -DMPI_CXX_COMPILER=/Users/zhili/Codes/local/bin/mpicxx`.
  **Never set `CMAKE_CXX_COMPILER=mpicxx`**: the wrapper injects a plain
  `-I/Users/zhili/Codes/local/include`, which defeats CMake's SYSTEM-header
  treatment and `-Wconversion -Werror` then fails *inside Kokkos headers*.
- **Local PETSc 3.25.1 has NO Kokkos backend** (plan §5.1 assumed it did).
  The COO API is backend-neutral, but `-mat_type aijkokkos` cannot be
  smoke-tested here; `aij` is the tested path. A Kokkos-enabled PETSc is a
  P5/CI concern.
- **Clang + sanitizer lane**: gcc on darwin/arm64 ships no ASan/UBSan
  runtimes, so the local sanitizer lane is Apple clang:
  Kokkos (Serial) + yaml-cpp rebuilt with clang live in
  `/Users/zhili/Codes/frehg2-upgrade/deps-clang` (sibling of the repo, not
  committed). Configure with `MPICH_CXX=/usr/bin/clang++` in the
  environment, `-DCMAKE_CXX_COMPILER=/usr/bin/clang++`, and
  `CMAKE_PREFIX_PATH="…/deps-clang;/Users/zhili/Codes/local"`, plus
  `-DFREHG_SANITIZE=ON`. Two PUBLIC defines on `frehg::core` make the
  gcc-built PETSc headers clang-safe: `PETSC_SKIP_REAL___FLOAT128`,
  `PETSC_SKIP_CXX_COMPLEX_FIX` — keep them.
- **`-ffp-contract=off` on test targets only**: gcc fuses multiply-adds on
  arm64 by default, so "expected" values computed at different inlined call
  sites can differ by 1 ULP from stored values — this broke the 0-ULP HDF5
  tests until the flag was added. Any *new* test target doing bitwise
  double comparisons needs the same flag (see `tests/CMakeLists.txt`).
- Kokkos backends here: OpenMP + Serial (gcc lane runs OpenMP; clang lane
  Serial). PetscInt is 32-bit; `Grid` guards against overflow.
- Doxygen 1.17 via Homebrew; `doxygen docs/Doxyfile` must stay exit-0
  (`WARN_AS_ERROR = FAIL_ON_WARNINGS`, undocumented public API is an
  error — document every member you add, including `///<` on struct
  fields).

## 5. Test inventory (all green at commit `6c51c3f`)

`ctest` in `build/`: 16 tests. Labels: `unit` (10) — GoogleTest suite
(`Logger/Timer`, `TimeSeries`, `Config` incl. 17-case invalid battery,
`GridBlocks/GridSerial`, `LinearSystem` Poisson/symmetry/COO-fold,
`Polygon/BoundarySet`, `GridDataReader`, `Hdf5/Monitor/Checkpoint`),
migration test, parameter-docs lockstep, 7 × `frehg --validate`;
`mpi` (6) — `frehg_mpi_halo` and `frehg_mpi_core` at n = 1, 2, 4
(bitwise halos both staging modes, gid coverage/permutation, NX=101 blocks,
rank-spanning polygons + sub-comms, parallel HDF5 layout, parallel
checkpoint). The `regression` label is empty until P1.

## 6. What P1 (and P2) should do first

1. Re-read plan §10 P1/P2 deliverables and §9 gates; re-run
   `scripts/ci_build_and_test.sh` to confirm the foundation is still green.
2. Add `src/swe/` (or `src/gw/`) as a STATIC library in `src/CMakeLists.txt`
   with explicit source lists, linking `frehg::core` and
   `frehg::compiler_flags`; physics libraries never include each other
   (plan §4) — cross-module data flows through driver-owned state.
3. The intended assembly pattern for the free-surface system: build COO
   row/col arrays once from `grid.gid2()` (5-point; out-of-domain legs get
   column −1 and their coefficient folded), then per step fill the values
   view in a Kokkos kernel and call `setValues` + `solve` on a
   `LinearSystem(comm, "fs_", grid.activeCount2Local(),
   grid.activeCount2Global())`. The 7-point subsurface system is the same
   with `gid3()` and prefix `"gw_"`.
4. The regression harness (`tests/regression/compare_h5.py`,
   `tools/ascii_golden_to_h5.py`, tolerance files) is **deliberately not in
   P0** — it lands with the first gated benchmark (plan §8.3), including
   the golden-fetch script (goldens never committed).
5. Every physics PR must cite the legacy `file:line` it reproduces (§11.1
   rule 7 — Appendix B has the authoritative list) and update the theory
   guide's equation table in the same PR (`docs/theory/` currently holds
   only `removed-features.md`; the equation tables start with P1).
6. If you drop anything newly found dead, add a row to
   `docs/theory/removed-features.md` in the same PR (binding, §3.2).

## 7. Open items and risks carried forward

- CI workflows are untested against real GitHub runners (no remote yet);
  expect one shakedown iteration on first push (apt package names, PETSc
  download URL, cache keys).
- The b6 early GW-only smoke comparison (plan §10 P2 exit) sizes the
  PCA-vs-Newton golden gap **before** P4 tolerances freeze — don't skip it.
- Restart determinism (1e-12 vs uninterrupted) is a P3 gate; P0 only
  guarantees the storage layer is bit-exact.
- `Hdf5Output::writeGridMeta` writes global axes replicated from every rank
  (fine at benchmark scales); revisit only if profiling ever flags it.
