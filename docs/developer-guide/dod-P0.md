# Definition of Done — P0 (Foundation)

Per upgrade plan §10 and §11.4. Every item names the command that verifies
it; all commands were run green on 2026-07-18 on macOS arm64 (gcc-15.2,
Apple clang 15, MPICH 4.1, PETSc 3.25.1, Kokkos 5.1.1, HDF5 1.14.5 parallel,
deps at `/Users/zhili/Codes/local`).

## Deliverables (plan §10 P0) — all implemented and tested

- [x] Repository, CMake targets and flags — `cmake -B build
      -DCMAKE_PREFIX_PATH=/Users/zhili/Codes/local
      -DCMAKE_CXX_COMPILER=g++-15
      -DMPI_CXX_COMPILER=/Users/zhili/Codes/local/bin/mpicxx &&
      cmake --build build -j` (zero warnings under
      `-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Werror`).
      `frehg::core` + `frehg` exist now; `frehg::{swe,gw,transport,driver}`
      are added by their phases per the plan's phase gating.
- [x] `PetscSession` (RAII MPI+Kokkos+PETSc), `Types`, `Logger`, `Timer` —
      `ctest --test-dir build -R unit.all` (Logger/Timer suites).
- [x] `Config` + `ConfigSchema` + `--validate` —
      `ctest --test-dir build -R 'unit.all|validate.'`
      (ConfigTest: acceptance + 17-case invalid battery with message checks).
- [x] `migrate_yaml_v1_to_v2.py` + migrated benchmark YAMLs —
      `ctest --test-dir build -R unit.migrate_yaml`; the seven configs
      (b1..b5, b6-ss, b6-td) live under `benchmarks/`.
- [x] `Grid` (compressed gids, non-divisible blocks, ktop masking) —
      `ctest --test-dir build -R 'unit.all|mpi.core'`
      (incl. NX=101 / 4-rank block math and masked-column compression).
- [x] `HaloExchanger` — `ctest --test-dir build -R mpi.halo`
      (bitwise ghosts at 1/2/4 ranks; both GPU-aware and host-staged modes;
      targeted exchange; domain-edge halos untouched).
- [x] `LinearSystem` (COO, `fs_`/`gw_`-style prefixes) —
      `ctest --test-dir build -R unit.all` (LinearSystem suite).
- [x] `TimeSeries` — `ctest --test-dir build -R unit.all`.
- [x] `Polygon`/`BoundarySet` (rasterization, sub-communicators) —
      `ctest --test-dir build -R 'unit.all|mpi.core'`
      (on-edge = inside; polygon spanning rank boundaries).
- [x] `Hdf5Output`, `Monitor`, `Checkpoint` —
      `ctest --test-dir build -R 'unit.all|mpi.core'`.
- [x] CI workflows — `.github/workflows/{build,sanitize,regression-nightly}.yml`
      (gcc+clang matrix, ASan/UBSan lane, nightly regression; runnable via
      `scripts/ci_build_and_test.sh` and `scripts/run_sanitizers.sh`).
- [x] `scripts/check_forbidden.sh` — `scripts/check_forbidden.sh` (clean).

## Exit criteria (plan §10 P0)

- [x] All §8.1 non-physics unit tests pass —
      `ctest --test-dir build -L unit --output-on-failure` (10/10).
- [x] Halo bitwise-exact at 1/2/4 ranks —
      `ctest --test-dir build -L mpi` (6/6: halo + parallel core at n=1,2,4).
- [x] Poisson ‖x−x*‖∞ < 1e-9 —
      `build/tests/frehg_unit_tests --gtest_filter='LinearSystem.*'`
      (1D with Dirichlet boundary fold, COO duplicate summation, 2D
      manufactured solution; symmetry checks included).
- [x] HDF5 round-trip 0 ULP —
      `build/tests/frehg_unit_tests --gtest_filter='Hdf5Test.*'`
      (bit-exact field and checkpoint round-trips; §7 flattening
      j*NX+i and (j*NX+i)*NZ+k asserted element-wise; NaN masking).
- [x] `frehg --validate` passes all benchmark YAMLs —
      `ctest --test-dir build -R validate.` (7/7 configs), and fails
      correctly on the invalid battery —
      `build/tests/frehg_unit_tests --gtest_filter='ConfigTest.Rejects*'`
      (17 crafted-invalid cases, message content asserted).
- [x] Zero-warning build, gcc **and** clang —
      gcc-15: `cmake --build build` (0 diagnostics);
      Apple clang 15: `MPICH_CXX=/usr/bin/clang++ cmake -B build-clang
      -DCMAKE_PREFIX_PATH="…/deps-clang;/Users/zhili/Codes/local"
      -DCMAKE_CXX_COMPILER=/usr/bin/clang++ … && cmake --build build-clang`
      (0 diagnostics; Kokkos/yaml-cpp rebuilt with clang in `deps-clang`
      because libstdc++/libc++ cannot mix).
- [x] ASan/UBSan clean —
      `ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1
      ctest --test-dir build-clang` (16/16 under both sanitizers).
      Note: gcc on macOS arm64 ships no sanitizer runtimes, so the local
      lane uses the clang toolchain; the Linux/gcc lane runs in CI
      (`.github/workflows/sanitize.yml`, `detect_leaks=1`).
- [x] Forbidden-pattern scan clean — `scripts/check_forbidden.sh`.
- [x] Doxygen `WARN_AS_ERROR` clean — `doxygen docs/Doxyfile`
      (WARN_AS_ERROR=FAIL_ON_WARNINGS, WARN_IF_UNDOCUMENTED=YES, exit 0).

## Additional P0 gates introduced by plan §11

- [x] Parameter docs in lockstep with the schema —
      `python3 scripts/check_parameter_docs.py`
      (also a ctest: `unit.parameter_docs`).
- [x] `docs/theory/removed-features.md` seeded from plan §3.2.
- [x] b6 golden/config sync-vs-async discrepancy recorded —
      `benchmarks/b6-kuan/README.md` (risk-register item).

## Known environment facts recorded for later phases

- The local PETSc 3.25.1 was built **without** Kokkos support: the COO
  assembly API used by `LinearSystem` is backend-neutral per plan §5.1, but
  `-mat_type aijkokkos` cannot be smoke-tested on this machine. The
  CPU-default `aij` path is what all local tests exercise.
- Kokkos host backends here: OpenMP + Serial (gcc build runs OpenMP;
  the clang sanitizer build runs Serial).
