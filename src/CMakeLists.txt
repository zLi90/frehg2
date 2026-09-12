# Frehg2 libraries and executable. Source lists are explicit by policy
# (Section 4 of the upgrade plan); source globbing is forbidden and CI-checked.
#
# frehg::core carries configuration, grid/indexing, halo exchange, boundary
# conditions, linear solvers, I/O, timers, and logging. frehg::swe (P1) is
# the surface-water module, frehg::gw (P2) the groundwater module,
# frehg::transport (P4) the scalar-transport module, frehg::coupling (P3)
# the surface-subsurface coupler, and frehg::driver the time loop. Physics
# libraries never link each other — cross-module data flows through
# driver-owned state mediated by the coupler and the transport wiring
# struct (plan §4).

add_library(frehg_core STATIC
  core/Logger.cpp
  core/Timer.cpp
  core/PetscSession.cpp
  core/TimeSeries.cpp
  core/Config.cpp
  core/ConfigSchema.cpp
  core/Grid.cpp
  core/HaloExchanger.cpp
  core/LinearSystem.cpp
  bc/Polygon.cpp
  bc/BoundarySet.cpp
  io/GridDataReader.cpp
  io/Hdf5Output.cpp
  io/Monitor.cpp
  io/Checkpoint.cpp
  io/RunRecord.cpp)
add_library(frehg::core ALIAS frehg_core)

target_include_directories(frehg_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(frehg_core
  PUBLIC
    Kokkos::kokkos
    MPI::MPI_CXX
    PETSc::PETSc
    hdf5::hdf5
    yaml-cpp::yaml-cpp
  PRIVATE
    frehg::compiler_flags)

target_compile_definitions(frehg_core
  PUBLIC
    FREHG_GIT_SHA="${FREHG_GIT_SHA}"
    FREHG_VERSION="${PROJECT_VERSION}"
    FREHG_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
    # Frehg2 is double-precision only (plan §4); skipping PETSc's optional
    # __float128 surface and C++ complex operator shims keeps a gcc-built
    # (libquadmath-aware) PETSc usable from compilers without __float128,
    # e.g. clang on arm64. Frehg2 uses neither quad precision nor
    # PetscComplex.
    PETSC_SKIP_REAL___FLOAT128
    PETSC_SKIP_CXX_COMPLEX_FIX
    $<$<BOOL:${FREHG_GPU_AWARE_MPI}>:FREHG_GPU_AWARE_MPI_DEFAULT=1>)

add_library(frehg_swe STATIC
  swe/SurfaceSolver.cpp
  swe/Momentum.cpp
  swe/FreeSurface.cpp
  swe/WetDry.cpp
  swe/SurfaceSources.cpp)
add_library(frehg::swe ALIAS frehg_swe)
target_link_libraries(frehg_swe
  PUBLIC frehg::core
  PRIVATE frehg::compiler_flags)

add_library(frehg_gw STATIC
  gw/TerrainMetric.cpp
  gw/RichardsSolver.cpp
  gw/Predictor.cpp
  gw/Corrector.cpp
  gw/Reallocate.cpp
  gw/AdaptiveStep.cpp
  gw/Baroclinic.cpp)
add_library(frehg::gw ALIAS frehg_gw)
target_link_libraries(frehg_gw
  PUBLIC frehg::core
  PRIVATE frehg::compiler_flags)

add_library(frehg_transport STATIC
  transport/ScalarSolver.cpp
  transport/SurfaceTransport.cpp
  transport/SubsurfaceTransport.cpp
  transport/Dispersion.cpp)
add_library(frehg::transport ALIAS frehg_transport)
target_link_libraries(frehg_transport
  PUBLIC frehg::core
  PRIVATE frehg::compiler_flags)

add_library(frehg_coupling STATIC
  coupling/Coupler.cpp
  coupling/Exchange.cpp)
add_library(frehg::coupling ALIAS frehg_coupling)
target_link_libraries(frehg_coupling
  PUBLIC frehg::core frehg::swe frehg::gw
  PRIVATE frehg::compiler_flags)

add_library(frehg_driver STATIC
  driver/Simulation.cpp)
add_library(frehg::driver ALIAS frehg_driver)
target_link_libraries(frehg_driver
  PUBLIC frehg::core frehg::swe frehg::gw frehg::transport frehg::coupling
  PRIVATE frehg::compiler_flags)

add_executable(frehg main.cpp)
target_link_libraries(frehg PRIVATE frehg::driver frehg::compiler_flags)
