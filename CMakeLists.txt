cmake_minimum_required(VERSION 3.23)

project(frehg2
  VERSION 1.0.0
  DESCRIPTION "Frehg2: coupled surface-water / groundwater / transport model"
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
endif()

option(FREHG_WERROR "Treat compiler warnings as errors" ON)
option(FREHG_SANITIZE "Enable Address and UndefinedBehavior sanitizers" OFF)
option(FREHG_ENABLE_TESTS "Build the test suite" ON)
option(FREHG_CLANG_TIDY "Run clang-tidy during compilation" OFF)
option(FREHG_GPU_AWARE_MPI "Compile-time default for GPU-aware MPI staging (runtime-overridable)" OFF)

# ---------------------------------------------------------------------------
# Dependencies (explicit find_package; no submodules, no bundled sources)
# ---------------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(Kokkos 5.1 REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
if(NOT HDF5_IS_PARALLEL)
  message(FATAL_ERROR "Frehg2 requires a parallel (MPI) HDF5 build; the HDF5 found at "
                      "'${HDF5_INCLUDE_DIRS}' was built without MPI support.")
endif()
find_package(yaml-cpp REQUIRED)
find_package(PETSc 3.20 REQUIRED)

include(FrehgCompilerFlags)
include(FrehgSanitizers)

if(FREHG_CLANG_TIDY)
  find_program(FREHG_CLANG_TIDY_EXE NAMES clang-tidy REQUIRED)
  set(CMAKE_CXX_CLANG_TIDY "${FREHG_CLANG_TIDY_EXE};--warnings-as-errors=*")
endif()

# Embed the git revision into the binary for output-file provenance.
find_package(Git QUIET)
set(FREHG_GIT_SHA "unknown")
if(GIT_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" -C "${CMAKE_CURRENT_SOURCE_DIR}" rev-parse --short=12 HEAD
    OUTPUT_VARIABLE FREHG_GIT_SHA_RAW
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE FREHG_GIT_RC)
  if(FREHG_GIT_RC EQUAL 0)
    set(FREHG_GIT_SHA "${FREHG_GIT_SHA_RAW}")
  endif()
endif()

add_subdirectory(src)

if(FREHG_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
