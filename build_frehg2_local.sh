#!/usr/bin/env bash
# build_frehg2_local.sh — build EVERY Frehg2 prerequisite from source into a
# private ./frehg-deps prefix, then build Frehg2 itself. No system package
# manager, no modules, no root: only a C/C++ compiler, cmake, make, curl,
# python3, and zlib headers (present on stock macOS and most Linux distros).
#
# RUN THIS ONE COMMAND from the Frehg2 top-level directory (the one holding
# CMakeLists.txt), e.g. inside frehg2-dev/:
#   bash build_frehg2_local.sh
#
# What it builds, in dependency order, all pinned and reused when already
# installed (delete ./frehg-deps to force a full rebuild):
#   1. MPICH 4.3.0            (MPI; C/C++ only, no Fortran)
#   2. HDF5 1.14.6            (parallel, compiled with the private mpicc)
#   3. yaml-cpp 0.8.0
#   4. Kokkos 5.1.1           (Serial + OpenMP host backends)
#   5. PETSc 3.25.1           (+ downloaded f2cblaslapack and hypre
#                              BoomerAMG — the v2 solver 'amg' path)
#   6. Frehg2                 (Release; executable at ./build/src/frehg)
# and finishes by validating the b1-sw benchmark and the Kuan
# saltwater-intrusion validation configs.
#
# Options (environment variables):
#   BUILD_JOBS=<n>       parallel build jobs        (default: CPU count)
#   FREHG_WERROR=ON|OFF  -Werror for Frehg2         (default: OFF — the
#                        zero-warning guarantee is pinned to gcc-15/clang-15;
#                        other compiler versions may warn harmlessly)
#   FREHG_ENABLE_TESTS=ON|OFF                        (default: OFF)
#   CC=... CXX=...       compilers (default: newest real GCC found, else
#                        clang; with clang, Kokkos OpenMP needs libomp —
#                        the script falls back to a Serial-only Kokkos and
#                        says so)
#
# After a successful build, run the Kuan case with:
#   cd validation/transport-kuan-saltwater-intrusion
#   OMP_NUM_THREADS=1 FI_PROVIDER=tcp \
#     ../../frehg-deps/bin/mpiexec -n 1 ../../build/src/frehg kuan-ss.yaml
#   OMP_NUM_THREADS=1 FI_PROVIDER=tcp \
#     ../../frehg-deps/bin/mpiexec -n 1 ../../build/src/frehg kuan-td.yaml
# (FI_PROVIDER=tcp pins MPICH's libfabric provider — the sockets default
# intermittently wedges MPI finalize on macOS; amendment A19. OMP_NUM_THREADS=1
# is fastest at benchmark grid sizes — report-P1.)

set -Eeuo pipefail
umask 022
ROOT_DIR="$PWD"
PREFIX="$ROOT_DIR/frehg-deps"
CACHE="$PREFIX/src-cache"
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="$ROOT_DIR/frehg2-build-$STAMP.log"
trap 'status=$?; echo; echo "BUILD FAILED at line $LINENO (exit $status). Full log: $LOG"; exit "$status"' ERR
exec > >(tee -a "$LOG") 2>&1

say() { printf '\n==================== %s ====================\n' "$*"; }
die() { echo "ERROR: $*"; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Required command '$1' is not installed."; }

MPICH_VERSION=4.3.0
HDF5_VERSION=1.14.6
YAMLCPP_VERSION=0.8.0
KOKKOS_VERSION=5.1.1
PETSC_VERSION=3.25.1

[[ -f "$ROOT_DIR/CMakeLists.txt" ]] || die "Run this script from the Frehg2 top-level directory containing CMakeLists.txt."
[[ -d "$ROOT_DIR/benchmarks/b1-sw" ]] || die "benchmarks/b1-sw is missing; this is not a complete Frehg2 checkout."
need cmake; need make; need curl; need python3; need tar
mkdir -p "$PREFIX" "$CACHE"

if command -v nproc >/dev/null 2>&1; then JOBS="${BUILD_JOBS:-$(nproc)}";
else JOBS="${BUILD_JOBS:-$(getconf _NPROCESSORS_ONLN)}"; fi

say "Frehg2 fully-private build (all prerequisites from source)"
echo "Repository:  $ROOT_DIR"
echo "Prefix:      $PREFIX"
echo "Build log:   $LOG"
echo "Build jobs:  $JOBS"

# ---------------------------------------------------------------------------
# 0/6: pick compilers. Prefer a real GCC (newest first) because the local
# zero-warning record and the Kokkos-OpenMP path are proven there; fall back
# to the system compiler. Do NOT set CC/CXX to MPI wrappers — CMake locates
# MPI itself, and wrapper -I flags defeat SYSTEM-header treatment.
# ---------------------------------------------------------------------------
say "0/6: select compilers"
if [[ -z "${CC:-}" || -z "${CXX:-}" ]]; then
  for v in 15 14 13 12; do
    if command -v "g++-$v" >/dev/null 2>&1; then CC="gcc-$v"; CXX="g++-$v"; break; fi
  done
fi
CC="${CC:-cc}"; CXX="${CXX:-c++}"
command -v "$CC" >/dev/null 2>&1 || die "C compiler '$CC' not found."
command -v "$CXX" >/dev/null 2>&1 || die "C++ compiler '$CXX' not found."
echo "CC:  $($CC --version | head -1)"
echo "CXX: $($CXX --version | head -1)"
KOKKOS_OPENMP=ON
if "$CXX" --version | head -1 | grep -qi clang; then
  # Apple/LLVM clang has no bundled libomp; probe for it.
  if ! echo 'int main(){return 0;}' | "$CXX" -fopenmp -x c++ - -o /tmp/frehg-omp-probe 2>/dev/null; then
    KOKKOS_OPENMP=OFF
    echo "NOTE: $CXX cannot link -fopenmp (no libomp); Kokkos will be built"
    echo "      Serial-only. Install a real GCC (e.g. brew install gcc) and"
    echo "      rebuild for the threaded backend."
  fi
  rm -f /tmp/frehg-omp-probe
fi
export CC CXX

fetch() {  # fetch URL FILE
  [[ -f "$2" ]] || curl --fail --location --retry 3 --output "$2" "$1"
}
extract_to() {  # extract_to TARBALL DEST-DIR (skips if DEST-DIR exists)
  [[ -d "$2" ]] && { echo "Reusing source tree: $2"; return 0; }
  # Extract into a scratch directory and move whatever single top-level
  # directory appears — robust to archives whose member names carry a
  # leading "./" (the HDF5 GitHub release tarballs do).
  local tmp="$CACHE/.extract-$$"
  rm -rf "$tmp"; mkdir -p "$tmp"
  tar -xf "$1" -C "$tmp"
  local entries=()
  while IFS= read -r line; do entries+=("$line"); done \
    < <(find "$tmp" -mindepth 1 -maxdepth 1 ! -name '.DS_Store')
  if [[ ${#entries[@]} -eq 1 && -d "${entries[0]}" ]]; then
    mv "${entries[0]}" "$2"
    rm -rf "$tmp"
  else
    mv "$tmp" "$2"   # flat archive: the scratch dir itself is the tree
  fi
  [[ -d "$2" ]] || die "Could not extract $1."
}
cmake_install() {  # cmake_install SRC BUILD [args...]
  local source="$1" build="$2"; shift 2
  cmake -S "$source" -B "$build" -G "Unix Makefiles" "$@"
  cmake --build "$build" --parallel "$JOBS"
  cmake --install "$build"
}

export PATH="$PREFIX/bin:$PATH"
export CMAKE_PREFIX_PATH="$PREFIX${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PREFIX/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
if [[ "$(uname -s)" == "Darwin" ]]; then
  export DYLD_LIBRARY_PATH="$PREFIX/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
else
  export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

# ---------------------------------------------------------------------------
say "1/6: MPICH $MPICH_VERSION (C/C++ only)"
if [[ -x "$PREFIX/bin/mpicc" && -x "$PREFIX/bin/mpiexec" ]]; then
  echo "MPICH already installed; reusing it."
else
  fetch "https://www.mpich.org/static/downloads/$MPICH_VERSION/mpich-$MPICH_VERSION.tar.gz" \
        "$CACHE/mpich-$MPICH_VERSION.tar.gz"
  extract_to "$CACHE/mpich-$MPICH_VERSION.tar.gz" "$CACHE/mpich-$MPICH_VERSION"
  pushd "$CACHE/mpich-$MPICH_VERSION" >/dev/null
  ./configure --prefix="$PREFIX" --disable-fortran --enable-shared \
              CC="$CC" CXX="$CXX"
  make -j "$JOBS"
  make install
  popd >/dev/null
fi
echo "mpicc: $("$PREFIX/bin/mpicc" -show | head -1)"

# ---------------------------------------------------------------------------
say "2/6: HDF5 $HDF5_VERSION (parallel, against the private MPICH)"
if [[ -x "$PREFIX/bin/h5pcc" ]] || "$PREFIX/bin/h5cc" -showconfig 2>/dev/null | grep -qi "parallel hdf5: yes"; then
  echo "Parallel HDF5 already installed; reusing it."
else
  fetch "https://github.com/HDFGroup/hdf5/releases/download/hdf5_$HDF5_VERSION/hdf5-$HDF5_VERSION.tar.gz" \
        "$CACHE/hdf5-$HDF5_VERSION.tar.gz"
  extract_to "$CACHE/hdf5-$HDF5_VERSION.tar.gz" "$CACHE/hdf5-$HDF5_VERSION"
  pushd "$CACHE/hdf5-$HDF5_VERSION" >/dev/null
  ./configure --prefix="$PREFIX" --enable-parallel --enable-shared \
              --disable-fortran --disable-cxx \
              CC="$PREFIX/bin/mpicc"
  make -j "$JOBS"
  make install
  popd >/dev/null
fi

# ---------------------------------------------------------------------------
say "3/6: yaml-cpp $YAMLCPP_VERSION and Kokkos $KOKKOS_VERSION (OpenMP=$KOKKOS_OPENMP)"
if [[ -f "$PREFIX/lib/cmake/yaml-cpp/yaml-cpp-config.cmake" || -f "$PREFIX/lib64/cmake/yaml-cpp/yaml-cpp-config.cmake" ]]; then
  echo "yaml-cpp already installed; reusing it."
else
  fetch "https://github.com/jbeder/yaml-cpp/archive/refs/tags/$YAMLCPP_VERSION.tar.gz" \
        "$CACHE/yaml-cpp-$YAMLCPP_VERSION.tar.gz"
  extract_to "$CACHE/yaml-cpp-$YAMLCPP_VERSION.tar.gz" "$CACHE/yaml-cpp-$YAMLCPP_VERSION"
  cmake_install "$CACHE/yaml-cpp-$YAMLCPP_VERSION" "$CACHE/build-yaml-cpp" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DYAML_BUILD_SHARED_LIBS=ON -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
fi
if [[ -f "$PREFIX/lib/cmake/Kokkos/KokkosConfig.cmake" || -f "$PREFIX/lib64/cmake/Kokkos/KokkosConfig.cmake" ]]; then
  echo "Kokkos already installed; reusing it."
else
  fetch "https://github.com/kokkos/kokkos/archive/refs/tags/$KOKKOS_VERSION.tar.gz" \
        "$CACHE/kokkos-$KOKKOS_VERSION.tar.gz"
  extract_to "$CACHE/kokkos-$KOKKOS_VERSION.tar.gz" "$CACHE/kokkos-$KOKKOS_VERSION"
  # BUILD_SHARED_LIBS=ON is load-bearing: frehg links Kokkos and so does
  # PETSc's downloaded kokkos-kernels. A static Kokkos core gets absorbed
  # into BOTH the frehg executable and libkokkoskernels.dylib, giving two
  # copies of Kokkos's global runtime state in one process -- the singleton
  # initializes twice ("Kokkos::OpenMP::initialize" prints twice) and the
  # first VecKokkos access segfaults. One shared libkokkoscore.dylib keeps a
  # single runtime shared across the executable and every dylib.
  cmake_install "$CACHE/kokkos-$KOKKOS_VERSION" "$CACHE/build-kokkos" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_CXX_STANDARD=20 \
    -DKokkos_ENABLE_SERIAL=ON "-DKokkos_ENABLE_OPENMP=$KOKKOS_OPENMP" \
    -DKokkos_ENABLE_TESTS=OFF
fi

# ---------------------------------------------------------------------------
say "4/6: PETSc $PETSC_VERSION (+ f2cblaslapack + hypre BoomerAMG + Kokkos backend)"
# Inherited PETSc environment variables poison configure (it aborts if
# PETSC_DIR points anywhere but the tree being built) and make; this build
# is fully self-contained, so drop them.
unset PETSC_DIR PETSC_ARCH
if [[ -f "$PREFIX/lib/petsc/conf/petscvariables" || -f "$PREFIX/lib64/petsc/conf/petscvariables" ]]; then
  echo "PETSc already installed; reusing it."
else
  fetch "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-$PETSC_VERSION.tar.gz" \
        "$CACHE/petsc-$PETSC_VERSION.tar.gz"
  extract_to "$CACHE/petsc-$PETSC_VERSION.tar.gz" "$CACHE/petsc-$PETSC_VERSION"
  pushd "$CACHE/petsc-$PETSC_VERSION" >/dev/null
  # --with-kokkos-dir + --download-kokkos-kernels (v2 plan §2B.2 B3): PETSc
  # reuses the Kokkos built in step 3, enabling -mat_type aijkokkos /
  # solver.*.mat_type=aijkokkos (the threaded-CPU and GPU solver lane).
  # --with-openmp also threads the downloaded hypre.
  ./configure --prefix="$PREFIX" \
    --with-cc="$PREFIX/bin/mpicc" --with-cxx="$PREFIX/bin/mpicxx" --with-fc=0 \
    --with-debugging=0 \
    --download-f2cblaslapack --download-hypre \
    --with-kokkos-dir="$PREFIX" --download-kokkos-kernels \
    --with-openmp=1 \
    COPTFLAGS=-O2 CXXOPTFLAGS=-O2
  make -j "$JOBS" all
  make install
  popd >/dev/null
fi
grep -q "PETSC_HAVE_HYPRE" "$PREFIX/include/petscconf.h" \
  || die "PETSc was installed without hypre; delete $PREFIX and rerun."
grep -q "PETSC_HAVE_KOKKOS_KERNELS" "$PREFIX/include/petscconf.h" \
  || die "PETSc was installed without Kokkos Kernels (mat_type aijkokkos needs it); delete $PREFIX and rerun."
# A reused PETSc must still have its Kokkos Kernels runtime present in the
# prefix -- petscconf.h advertises KOKKOS_KERNELS even after the library is
# gone. This catches "Kokkos rebuilt but PETSc reused" (e.g. the static->shared
# migration of V2-A8): the old libkokkoskernels.dylib is removed with the old
# Kokkos, libpetsc dangles, and the failure is otherwise a cryptic dyld abort
# at first run instead of an actionable message here.
[[ -f "$PREFIX/lib/libkokkoskernels.dylib" || -f "$PREFIX/lib/libkokkoskernels.so" ]] \
  || die "PETSc's Kokkos Kernels runtime is missing from $PREFIX/lib (Kokkos was likely rebuilt without PETSc). Delete '$PREFIX/lib/petsc' and '$PREFIX/src-cache/petsc-$PETSC_VERSION', then rerun to rebuild PETSc against the current Kokkos."

# ---------------------------------------------------------------------------
say "5/6: configure and compile Frehg2"
rm -rf "$ROOT_DIR/build"
# HDF5_ROOT + NO_FIND_PACKAGE_CONFIG_FILE: a Homebrew HDF5 ships CMake
# config files, which find_package prefers over module mode — bypassing
# HDF5_PREFER_PARALLEL and selecting the serial build. Forcing module mode
# rooted at the private prefix parses our parallel h5pcc instead.
cmake -S "$ROOT_DIR" -B "$ROOT_DIR/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" \
  -DMPI_HOME="$PREFIX" \
  -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" \
  -DHDF5_ROOT="$PREFIX" \
  -DHDF5_NO_FIND_PACKAGE_CONFIG_FILE=TRUE \
  -DHDF5_PREFER_PARALLEL=TRUE \
  -DFREHG_ENABLE_TESTS="${FREHG_ENABLE_TESTS:-OFF}" \
  -DFREHG_WERROR="${FREHG_WERROR:-OFF}"
cmake --build "$ROOT_DIR/build" --parallel "$JOBS"
EXE="$ROOT_DIR/build/src/frehg"
[[ -x "$EXE" ]] || die "The expected executable was not created: $EXE"
echo "SUCCESS: executable created: $EXE"

# ---------------------------------------------------------------------------
say "6/6: validate benchmark and Kuan validation inputs"
export OMP_NUM_THREADS=1
export FI_PROVIDER=tcp
"$EXE" --validate "$ROOT_DIR/benchmarks/b1-sw/b1-sw.yaml"
if [[ -d "$ROOT_DIR/validation/transport-kuan-saltwater-intrusion" ]]; then
  ( cd "$ROOT_DIR/validation/transport-kuan-saltwater-intrusion" && \
    "$EXE" --validate kuan-ss.yaml && "$EXE" --validate kuan-td.yaml )
fi

echo
echo "================================================================"
echo "FREHG2 BUILD COMPLETE"
echo "Executable:        $EXE"
echo "Dependency prefix: $PREFIX   (mpiexec: $PREFIX/bin/mpiexec)"
echo "Build log:         $LOG"
echo
echo "Run the Kuan saltwater-intrusion case:"
echo "  cd validation/transport-kuan-saltwater-intrusion"
echo "  OMP_NUM_THREADS=1 FI_PROVIDER=tcp \\"
echo "    $PREFIX/bin/mpiexec -n 1 $EXE kuan-ss.yaml"
echo "  (then kuan-td.yaml; plot with: python3 makeplot.py)"
echo "================================================================"
