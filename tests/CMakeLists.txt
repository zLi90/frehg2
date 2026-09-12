# Frehg2 test suite: GoogleTest unit tests (label "unit") and plain MPI
# drivers run at 1/2/4 ranks (label "mpi"). Regression harness cases join in
# P1+ (label "regression").

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
  URL_HASH SHA256=1f357c27ca988c3f7c6b4bf68a9395005ac6761f034046e9dde0896e3aba00e4
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# ---------------------------------------------------------------------------
# Unit tests (single rank)
# ---------------------------------------------------------------------------
add_executable(frehg_unit_tests
  unit/main.cpp
  unit/test_logger_timer.cpp
  unit/test_timeseries.cpp
  unit/test_config.cpp
  unit/test_grid.cpp
  unit/test_linearsystem.cpp
  unit/test_polygon.cpp
  unit/test_griddatareader.cpp
  unit/test_hdf5.cpp
  unit/test_swe.cpp
  unit/test_vangenuchten.cpp
  unit/test_gw.cpp
  unit/test_coupling.cpp
  unit/test_transport.cpp)
target_link_libraries(frehg_unit_tests PRIVATE frehg::core frehg::swe frehg::gw
                                               frehg::transport frehg::coupling GTest::gtest
                                               frehg::compiler_flags)
target_compile_definitions(frehg_unit_tests PRIVATE
  FREHG_REPO_DIR="${CMAKE_SOURCE_DIR}")

add_test(NAME unit.all COMMAND frehg_unit_tests)
set_tests_properties(unit.all PROPERTIES LABELS "unit" PROCESSORS 1)

# p5 negative test (v2 plan §2B.3): the backend-consistency invariant in
# LinearSystem.hpp must REJECT a device MemSpace paired with a PETSc that
# lacks Kokkos. The target is excluded from the normal build; the test
# builds it on demand and passes only when compilation FAILS.
add_executable(frehg_compile_fail_backend EXCLUDE_FROM_ALL unit/compile_fail_backend.cpp)
target_link_libraries(frehg_compile_fail_backend PRIVATE frehg::core frehg::compiler_flags)
add_test(NAME unit.p5.backend_invariant_fires
         COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}"
                 --target frehg_compile_fail_backend)
set_tests_properties(unit.p5.backend_invariant_fires PROPERTIES
  LABELS "unit" PROCESSORS 1 WILL_FAIL TRUE)

# The tests compare doubles bitwise (0-ULP round-trip contracts); disable FP
# contraction so "expected" values computed at different inlined call sites
# cannot differ by a fused multiply-add.
target_compile_options(frehg_unit_tests PRIVATE -ffp-contract=off)

# ---------------------------------------------------------------------------
# MPI tests (1, 2, and 4 ranks)
# ---------------------------------------------------------------------------
add_executable(frehg_mpi_halo mpi/test_halo.cpp)
target_link_libraries(frehg_mpi_halo PRIVATE frehg::core frehg::compiler_flags)
target_compile_options(frehg_mpi_halo PRIVATE -ffp-contract=off)

add_executable(frehg_mpi_core mpi/test_parallel_core.cpp)
target_link_libraries(frehg_mpi_core PRIVATE frehg::core frehg::compiler_flags)
target_compile_options(frehg_mpi_core PRIVATE -ffp-contract=off)

foreach(nranks 1 2 4)
  add_test(NAME mpi.halo.n${nranks}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nranks}
            $<TARGET_FILE:frehg_mpi_halo>)
  set_tests_properties(mpi.halo.n${nranks} PROPERTIES LABELS "mpi" PROCESSORS ${nranks})

  add_test(NAME mpi.core.n${nranks}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nranks}
            $<TARGET_FILE:frehg_mpi_core>)
  set_tests_properties(mpi.core.n${nranks} PROPERTIES LABELS "mpi" PROCESSORS ${nranks})
endforeach()

# ---------------------------------------------------------------------------
# Tooling tests: the v1 -> v2 migration script against the archival drafts
# ---------------------------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME unit.migrate_yaml
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/test_migrate_yaml.py
          --repo ${CMAKE_SOURCE_DIR} --frehg $<TARGET_FILE:frehg>)
set_tests_properties(unit.migrate_yaml PROPERTIES LABELS "unit")

# Schema <-> parameter-docs lockstep check (plan §5.4, §11.1 rule 3).
add_test(NAME unit.parameter_docs
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/check_parameter_docs.py)
set_tests_properties(unit.parameter_docs PROPERTIES LABELS "unit")

# Run-record checker sensitivity (v2 plan §2A gate r1, §6.3 negative tests):
# a synthetic valid record passes; every structural corruption fails.
add_test(NAME unit.run_record_checker
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/test_run_record_checker.py
          --checker ${CMAKE_SOURCE_DIR}/tools/check_run_record.py)
set_tests_properties(unit.run_record_checker PROPERTIES LABELS "unit")

# p2/p3 gate-logic negative battery (v2 plan §6.3): each p-gate assertion is
# provably able to fail on synthetic result sets.
add_test(NAME unit.p_gates.negative
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/test_p_gates.py)
set_tests_properties(unit.p_gates.negative PROPERTIES LABELS "unit")

# ---------------------------------------------------------------------------
# Benchmark configuration validation (part of the P0 exit criteria)
# ---------------------------------------------------------------------------
foreach(case
    b1-sw/b1-sw
    b2-gw/b2-gw
    b3-kirkland/b3-kirkland
    b4-govindaraju/b4-govindaraju
    b5-vcatchment/b5-vcatchment
    b5-vcatchment/b5-vcatchment-norain
    b6-kuan/b6-kuan-ss
    b6-kuan/b6-kuan-td)
  get_filename_component(case_name ${case} NAME)
  add_test(NAME validate.${case_name}
    COMMAND $<TARGET_FILE:frehg> --validate
            ${CMAKE_SOURCE_DIR}/benchmarks/${case}.yaml)
  set_tests_properties(validate.${case_name} PROPERTIES LABELS "unit")
endforeach()

# ---------------------------------------------------------------------------
# Regression gates (plan §8.3, §9): P1 gates b1 and b4. Goldens are fetched
# from the sibling legacy checkout and never committed. OpenMP thread counts
# are pinned: these are launch-latency-bound miniature grids.
# ---------------------------------------------------------------------------
set(FREHG_LEGACY_BENCHMARKS "${CMAKE_SOURCE_DIR}/../legacy/benchmarks" CACHE PATH
    "Directory holding the legacy benchmark goldens (plan §8.3)")
set(_regress ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/regression/run_regression.py)
set(_regress_args
    --frehg $<TARGET_FILE:frehg>
    --mpiexec ${MPIEXEC_EXECUTABLE}
    --repo ${CMAKE_SOURCE_DIR}
    --legacy ${FREHG_LEGACY_BENCHMARKS})
# Benchmark-scale grids are launch-latency-bound: one OpenMP thread is the
# fastest configuration (the P1 performance note in report-P1.md has the
# measurements), and it keeps the gates deterministic. FI_PROVIDER pins
# libfabric to its tcp provider: the default sockets provider intermittently
# wedges MPICH's OFI finalize hook on this platform (the P1/P3 reports'
# known issue; at P4 it became persistent — post-completion, results
# unaffected — and the tcp provider exits cleanly).
set(_regress_env "OMP_NUM_THREADS=1;OMP_PROC_BIND=false;FI_PROVIDER=tcp")

add_test(NAME regression.b1
  COMMAND ${_regress} b1 ${_regress_args} --work ${CMAKE_BINARY_DIR}/regression/b1)
set_tests_properties(regression.b1 PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 900)

add_test(NAME regression.b2
  COMMAND ${_regress} b2 ${_regress_args} --work ${CMAKE_BINARY_DIR}/regression/b2)
set_tests_properties(regression.b2 PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 900)

add_test(NAME regression.b3
  COMMAND ${_regress} b3 ${_regress_args} --work ${CMAKE_BINARY_DIR}/regression/b3)
set_tests_properties(regression.b3 PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 1800)

add_test(NAME regression.b4
  COMMAND ${_regress} b4 ${_regress_args} --work ${CMAKE_BINARY_DIR}/regression/b4)
set_tests_properties(regression.b4 PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 1800)

add_test(NAME regression.b1_restart
  COMMAND ${_regress} b1-restart ${_regress_args}
          --work ${CMAKE_BINARY_DIR}/regression/b1-restart)
set_tests_properties(regression.b1_restart PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 1200)

add_test(NAME regression.b2_restart
  COMMAND ${_regress} b2-restart ${_regress_args}
          --work ${CMAKE_BINARY_DIR}/regression/b2-restart)
set_tests_properties(regression.b2_restart PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 1200)

# Plan §10 P2 exit: the early b6 groundwater-only smoke comparison against
# the syncV4 steady-state golden. Recorded (always passes), sizing the P4
# tolerance risk; the measured deviations are archived in dod-P2.md.
add_test(NAME regression.b6_gw_smoke
  COMMAND ${_regress} b6-gw-smoke ${_regress_args}
          --work ${CMAKE_BINARY_DIR}/regression/b6-gw-smoke)
set_tests_properties(regression.b6_gw_smoke PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 1800)

foreach(mode strict default)
  add_test(NAME regression.b1_rank_invariance.${mode}
    COMMAND ${_regress} rank-invariance --mode ${mode} ${_regress_args}
            --work ${CMAKE_BINARY_DIR}/regression/rank-${mode})
  set_tests_properties(regression.b1_rank_invariance.${mode} PROPERTIES
    LABELS "regression;mpi" ENVIRONMENT "${_regress_env}" PROCESSORS 4 TIMEOUT 1800)

  add_test(NAME regression.b2_rank_invariance.${mode}
    COMMAND ${_regress} rank-invariance-b2 --mode ${mode} ${_regress_args}
            --work ${CMAKE_BINARY_DIR}/regression/rank-b2-${mode})
  set_tests_properties(regression.b2_rank_invariance.${mode} PROPERTIES
    LABELS "regression;mpi" ENVIRONMENT "${_regress_env}" PROCESSORS 4 TIMEOUT 1800)

  add_test(NAME regression.b5_rank_invariance.${mode}
    COMMAND ${_regress} rank-invariance-b5 --mode ${mode} ${_regress_args}
            --work ${CMAKE_BINARY_DIR}/regression/rank-b5-${mode})
  set_tests_properties(regression.b5_rank_invariance.${mode} PROPERTIES
    LABELS "regression;mpi" ENVIRONMENT "${_regress_env}" PROCESSORS 4 TIMEOUT 3600)
endforeach()

# ---------------------------------------------------------------------------
# P3 coupled gates (plan §9 b5, §10 P3). The four full-horizon envelope runs
# are nightly-class (plan §11.2 keeps b5/b6 out of the per-PR set): hours of
# wall time each at 4 ranks. The restart gate is short and per-PR-sized.
# ---------------------------------------------------------------------------
add_test(NAME regression.b5_restart
  COMMAND ${_regress} b5-restart ${_regress_args} --ranks 4
          --work ${CMAKE_BINARY_DIR}/regression/b5-restart)
set_tests_properties(regression.b5_restart PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" PROCESSORS 4 TIMEOUT 3600)

foreach(scenario rain norain)
  foreach(coupling sync subcycled)
    add_test(NAME regression.b5.${scenario}.${coupling}
      COMMAND ${_regress} b5 --scenario ${scenario} --coupling ${coupling}
              ${_regress_args} --ranks 4
              --work ${CMAKE_BINARY_DIR}/regression/b5-${scenario}-${coupling})
    set_tests_properties(regression.b5.${scenario}.${coupling} PROPERTIES
      LABELS "regression_nightly" ENVIRONMENT "${_regress_env}" PROCESSORS 4
      TIMEOUT 43200)
  endforeach()
endforeach()

# ---------------------------------------------------------------------------
# P4 transport gates (plan §9 b6, §10 P4). The two full-horizon b6 runs stay
# out of the per-PR set (plan §11.2 runs b6 nightly); the restart gate is
# short and per-PR-sized. b6 is serial like the legacy golden runs (a 1x68
# tank; decomposition adds nothing at this size).
# ---------------------------------------------------------------------------
add_test(NAME regression.b6_restart
  COMMAND ${_regress} b6-restart ${_regress_args}
          --work ${CMAKE_BINARY_DIR}/regression/b6-restart)
set_tests_properties(regression.b6_restart PROPERTIES LABELS "regression"
  ENVIRONMENT "${_regress_env}" TIMEOUT 1800)

foreach(variant ss td)
  add_test(NAME regression.b6.${variant}
    COMMAND ${_regress} b6 --variant ${variant} ${_regress_args}
            --work ${CMAKE_BINARY_DIR}/regression/b6-${variant})
  set_tests_properties(regression.b6.${variant} PROPERTIES
    LABELS "regression_nightly" ENVIRONMENT "${_regress_env}" TIMEOUT 43200)
endforeach()

# ---------------------------------------------------------------------------
# Sanitized builds run 3-10x slower than the timeouts above assume; the
# sanitizer lane configures with a multiplier instead of maintaining a
# second set of numbers (plan §10 P5, scripts/run_sanitizers.sh).
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# v2 Q1 gates (v2 plan §2.3, §7): g1 solver invariance, g2 iteration
# flatness, the s-gate scalability matrix, and the per-module perf
# baseline. Gate-first (v2 plan §6.1): these entries land red and are
# turned green by the Q1 solver work.
# ---------------------------------------------------------------------------
set(_scaling ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/run_scaling.py
    --frehg $<TARGET_FILE:frehg> --mpiexec ${MPIEXEC_EXECUTABLE}
    --repo ${CMAKE_SOURCE_DIR})
set(_scaling_env "FI_PROVIDER=tcp;OMP_PROC_BIND=false")

# g1 per-PR lanes: the fast gates (b1-b4) under each new preconditioner, and
# the default-mode rank-invariance lanes under amg. b5/b6 amg lanes are
# nightly-class like their parents.
foreach(pc amg gamg)
  foreach(gate b1 b2 b3 b4)
    add_test(NAME regression.g1.${gate}.${pc}
      COMMAND ${_regress} ${gate} --solver ${pc} ${_regress_args}
              --work ${CMAKE_BINARY_DIR}/regression/g1-${gate}-${pc})
    set_tests_properties(regression.g1.${gate}.${pc} PROPERTIES
      LABELS "regression" ENVIRONMENT "${_regress_env}" TIMEOUT 1800)
  endforeach()
endforeach()

foreach(variant "rank-invariance;b1" "rank-invariance-b2;b2")
  list(GET variant 0 subcommand)
  list(GET variant 1 short)
  add_test(NAME regression.g1.${short}_rank_invariance.default.amg
    COMMAND ${_regress} ${subcommand} --mode default --solver amg ${_regress_args}
            --work ${CMAKE_BINARY_DIR}/regression/g1-rank-${short}-amg)
  set_tests_properties(regression.g1.${short}_rank_invariance.default.amg PROPERTIES
    LABELS "regression;mpi" ENVIRONMENT "${_regress_env}" PROCESSORS 4 TIMEOUT 1800)
endforeach()

add_test(NAME regression.g1.b5_rank_invariance.default.amg
  COMMAND ${_regress} rank-invariance-b5 --mode default --solver amg ${_regress_args}
          --work ${CMAKE_BINARY_DIR}/regression/g1-rank-b5-amg)
set_tests_properties(regression.g1.b5_rank_invariance.default.amg PROPERTIES
  LABELS "regression_nightly" ENVIRONMENT "${_regress_env}" PROCESSORS 4 TIMEOUT 3600)

foreach(scenario rain)
  # The parent gate's approved record is the A15/A17 shortened 24 h horizon
  # (the release pipeline gates b5 with --t-end 86400); g1 holds the amg
  # lane to the parent's record, not to the never-adjudicated full horizon
  # (amendment V2-A4).
  add_test(NAME regression.g1.b5.${scenario}.sync.amg
    COMMAND ${_regress} b5 --scenario ${scenario} --coupling sync --solver amg
            --t-end 86400 ${_regress_args} --ranks 4
            --work ${CMAKE_BINARY_DIR}/regression/g1-b5-${scenario}-amg)
  set_tests_properties(regression.g1.b5.${scenario}.sync.amg PROPERTIES
    LABELS "regression_nightly" ENVIRONMENT "${_regress_env}" PROCESSORS 4
    TIMEOUT 43200)
endforeach()

foreach(variant ss td)
  add_test(NAME regression.g1.b6.${variant}.amg
    COMMAND ${_regress} b6 --variant ${variant} --solver amg ${_regress_args}
            --work ${CMAKE_BINARY_DIR}/regression/g1-b6-${variant}-amg)
  set_tests_properties(regression.g1.b6.${variant}.amg PROPERTIES
    LABELS "regression_nightly" ENVIRONMENT "${_regress_env}" TIMEOUT 43200)
endforeach()

# g2 (per-PR, v2 plan §2.3): iteration flatness at {1,2,4,8} ranks on the
# A23 synthetic case under amg; block-Jacobi's growth is recorded nightly
# alongside (the documented motivation, not a gate).
add_test(NAME scaling.g2.iters.amg
  COMMAND ${_scaling} --gate iters --solver amg --ranks 1 2 4 8
          --tolerances ${CMAKE_SOURCE_DIR}/tests/regression/tolerances/g2.yaml
          --json ${CMAKE_BINARY_DIR}/scaling/g2-amg.json
          --work ${CMAKE_BINARY_DIR}/scaling/g2-amg)
set_tests_properties(scaling.g2.iters.amg PROPERTIES
  LABELS "regression" ENVIRONMENT "${_scaling_env};OMP_NUM_THREADS=1"
  PROCESSORS 8 TIMEOUT 3600)

add_test(NAME scaling.g2.record.bjacobi
  COMMAND ${_scaling} --gate none --solver bjacobi-icc --ranks 1 2 4 8
          --json ${CMAKE_BINARY_DIR}/scaling/g2-bjacobi.json
          --work ${CMAKE_BINARY_DIR}/scaling/g2-bjacobi)
set_tests_properties(scaling.g2.record.bjacobi PROPERTIES
  LABELS "scaling_nightly" ENVIRONMENT "${_scaling_env};OMP_NUM_THREADS=1"
  PROCESSORS 8 TIMEOUT 3600)

# s-gates (nightly, v2 plan §7.2): strong, weak, threads, hybrid.
add_test(NAME scaling.s1.strong
  COMMAND ${_scaling} --gate strong --ranks 1 2 4 --repeats 3
          --json ${CMAKE_BINARY_DIR}/scaling/s1.json
          --work ${CMAKE_BINARY_DIR}/scaling/s1)
set_tests_properties(scaling.s1.strong PROPERTIES
  LABELS "scaling_nightly" ENVIRONMENT "${_scaling_env};OMP_NUM_THREADS=1"
  PROCESSORS 4 TIMEOUT 7200)

add_test(NAME scaling.s2.weak.amg
  COMMAND ${_scaling} --gate weak --solver amg --ranks 1 2 4 --repeats 3
          --json ${CMAKE_BINARY_DIR}/scaling/s2.json
          --work ${CMAKE_BINARY_DIR}/scaling/s2)
set_tests_properties(scaling.s2.weak.amg PROPERTIES
  LABELS "scaling_nightly" ENVIRONMENT "${_scaling_env};OMP_NUM_THREADS=1"
  PROCESSORS 4 TIMEOUT 7200)

add_test(NAME scaling.s3.threads
  COMMAND ${_scaling} --gate threads --threads 1 2 4 --repeats 3
          --json ${CMAKE_BINARY_DIR}/scaling/s3.json
          --work ${CMAKE_BINARY_DIR}/scaling/s3)
set_tests_properties(scaling.s3.threads PROPERTIES
  LABELS "scaling_nightly" ENVIRONMENT "${_scaling_env}" PROCESSORS 4 TIMEOUT 7200)

add_test(NAME scaling.s4.hybrid
  COMMAND ${_scaling} --gate hybrid --repeats 3
          --json ${CMAKE_BINARY_DIR}/scaling/s4.json
          --work ${CMAKE_BINARY_DIR}/scaling/s4)
set_tests_properties(scaling.s4.hybrid PROPERTIES
  LABELS "scaling_nightly" ENVIRONMENT "${_scaling_env}" PROCESSORS 4 TIMEOUT 7200)

# Per-module perf regression vs the tracked baseline (v2 plan §7.1.2).
add_test(NAME regression.perf_baseline
  COMMAND ${_scaling} --gate perf-baseline --repeats 3
          --baseline ${CMAKE_SOURCE_DIR}/docs/developer-guide/perf-baseline.json
          --work ${CMAKE_BINARY_DIR}/scaling/perf-baseline)
set_tests_properties(regression.perf_baseline PROPERTIES
  LABELS "regression" ENVIRONMENT "${_scaling_env};OMP_NUM_THREADS=1" TIMEOUT 1800)

# r2 — timer coverage and record overhead (v2 plan §2A.3): direct children of
# the time loop cover >= 90 % at 1 rank; halo/io sections nonzero at 4 ranks;
# writing the run record costs <= 1 % of the loop. Reads the run record, not
# stdout. (r1 has no entry of its own: every regression gate's run_case now
# checks the record in-line.)
add_test(NAME regression.r2_timer_coverage
  COMMAND ${_scaling} --gate timer-coverage --repeats 1
          --work ${CMAKE_BINARY_DIR}/scaling/timer-coverage)
set_tests_properties(regression.r2_timer_coverage PROPERTIES
  LABELS "regression" ENVIRONMENT "${_scaling_env};OMP_NUM_THREADS=1"
  PROCESSORS 4 TIMEOUT 1800)

# ---------------------------------------------------------------------------
# v2 Q3 p-gates (v2 plan §2B.3): performance portability. Gate-first (§6.1).
# p4 (device compile+link) lives in .github/workflows/cuda-compile.yml; the
# p5 negative test is unit.p5.backend_invariant_fires above; p6 is the
# owner-executed bundle (scripts/gpu_acceptance.sh).
# ---------------------------------------------------------------------------

# p1 — backend invariance (per-PR): b1 and b4 under {OpenMP-4 + aij} and
# {OpenMP-4 + aijkokkos} must pass the unchanged gate criteria (the existing
# regression.b1/b4 entries are the {1-thread + aij} corner). This is the CPU
# rehearsal of the exact linear-algebra path a GPU build takes.
foreach(mt aij aijkokkos)
  foreach(gate b1 b4)
    add_test(NAME regression.p1.${gate}.omp4.${mt}
      COMMAND ${_regress} ${gate} --mat-type ${mt} ${_regress_args}
              --work ${CMAKE_BINARY_DIR}/regression/p1-${gate}-${mt})
    set_tests_properties(regression.p1.${gate}.omp4.${mt} PROPERTIES
      LABELS "regression" ENVIRONMENT "FI_PROVIDER=tcp;OMP_NUM_THREADS=4;OMP_PROC_BIND=false"
      PROCESSORS 4 TIMEOUT 3600)
  endforeach()
endforeach()

# p1 nightly extension: b6 (coupled + transport + density, the full physics)
# under aijkokkos at 4 threads.
add_test(NAME regression.p1.b6.ss.omp4.aijkokkos
  COMMAND ${_regress} b6 --variant ss --mat-type aijkokkos ${_regress_args}
          --work ${CMAKE_BINARY_DIR}/regression/p1-b6-ss-aijkokkos)
set_tests_properties(regression.p1.b6.ss.omp4.aijkokkos PROPERTIES
  LABELS "regression_nightly" ENVIRONMENT "FI_PROVIDER=tcp;OMP_NUM_THREADS=4;OMP_PROC_BIND=false"
  PROCESSORS 4 TIMEOUT 43200)

# p2 — OpenMP solver thread scaling (nightly): on the aijkokkos + amg path
# the gw solve must actually thread (4-thread solve < 1-thread solve) with
# thread-invariant algebra (iteration drift <= 2 %, volumes to 1e-8).
add_test(NAME scaling.p2.solver_threads
  COMMAND ${_scaling} --gate p2-solver-threads --threads 1 2 4 --repeats 3
          --json ${CMAKE_BINARY_DIR}/scaling/p2.json
          --work ${CMAKE_BINARY_DIR}/scaling/p2)
set_tests_properties(scaling.p2.solver_threads PROPERTIES
  LABELS "scaling_nightly" ENVIRONMENT "${_scaling_env}" PROCESSORS 4 TIMEOUT 7200)

# p3 — hybrid MPI+OpenMP placements (nightly): 4 PEs as 4x1 / 2x2 / 1x4 on
# aijkokkos + amg; volumes to 1e-8, iteration drift <= 10 % vs all-MPI,
# fastest placement recorded.
add_test(NAME scaling.p3.hybrid
  COMMAND ${_scaling} --gate p3-hybrid --ranks 4 --repeats 3
          --json ${CMAKE_BINARY_DIR}/scaling/p3.json
          --work ${CMAKE_BINARY_DIR}/scaling/p3)
set_tests_properties(scaling.p3.hybrid PROPERTIES
  LABELS "scaling_nightly" ENVIRONMENT "${_scaling_env}" PROCESSORS 4 TIMEOUT 7200)

set(FREHG_TEST_TIMEOUT_SCALE "1" CACHE STRING
    "Multiplier applied to every test TIMEOUT (sanitizer builds set 5)")
if(NOT FREHG_TEST_TIMEOUT_SCALE EQUAL 1)
  get_property(_frehg_tests DIRECTORY PROPERTY TESTS)
  foreach(_frehg_test ${_frehg_tests})
    get_test_property(${_frehg_test} TIMEOUT _frehg_timeout)
    if(_frehg_timeout)
      math(EXPR _frehg_timeout "${_frehg_timeout} * ${FREHG_TEST_TIMEOUT_SCALE}")
      set_tests_properties(${_frehg_test} PROPERTIES TIMEOUT ${_frehg_timeout})
    endif()
  endforeach()
endif()
