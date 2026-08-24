# FindPETSc.cmake — locate PETSc through pkg-config and expose it as the
# imported target PETSc::PETSc.
#
# PETSc installs a pkg-config file (PETSc.pc) but no CMake package config, so
# this wrapper bridges the two worlds. Point CMAKE_PREFIX_PATH (or
# PKG_CONFIG_PATH) at the PETSc installation prefix.
#
# Result variables:
#   PETSc_FOUND, PETSc_VERSION
# Imported target:
#   PETSc::PETSc

find_package(PkgConfig REQUIRED)

set(_frehg_saved_pkg_config_path "$ENV{PKG_CONFIG_PATH}")
foreach(_prefix IN LISTS CMAKE_PREFIX_PATH)
  if(EXISTS "${_prefix}/lib/pkgconfig")
    set(ENV{PKG_CONFIG_PATH} "${_prefix}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  endif()
endforeach()

pkg_check_modules(_frehg_petsc QUIET IMPORTED_TARGET PETSc)
if(NOT _frehg_petsc_FOUND)
  pkg_check_modules(_frehg_petsc QUIET IMPORTED_TARGET petsc)
endif()

set(ENV{PKG_CONFIG_PATH} "${_frehg_saved_pkg_config_path}")

if(_frehg_petsc_FOUND)
  set(PETSc_VERSION "${_frehg_petsc_VERSION}")
  if(NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc INTERFACE IMPORTED)
    target_link_libraries(PETSc::PETSc INTERFACE PkgConfig::_frehg_petsc)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
  REQUIRED_VARS _frehg_petsc_FOUND
  VERSION_VAR PETSc_VERSION)
