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
