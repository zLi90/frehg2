#!/usr/bin/env bash
#SBATCH --job-name=frehg2-build
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=frehg2-build-%j.out
#SBATCH --error=frehg2-build-%j.err
#
# Frehg2 automated HPC build script (SLURM)
#
# WHAT THIS SCRIPT DOES
#   1. Finds and loads usable module files for CMake, a C++20 compiler,
#      MPI, and (when available) parallel HDF5.
#   2. Builds missing user-space dependencies in $HOME/frehg-deps:
#      yaml-cpp 0.8.0, Kokkos 5.1.1 (shared, Serial+OpenMP), parallel
#      HDF5 1.14.6 if necessary, and PETSc 3.25.1 built Kokkos-aware
#      (--with-kokkos-dir --download-kokkos-kernels --with-openmp) plus
#      hypre BoomerAMG. No sudo or administrator privileges are needed.
#   3. Configures and compiles Frehg2 in build/ and confirms that
#      build/src/frehg exists and is executable.
#   4. Runs a schema validation and a short SLURM benchmark verification.
#
# BEFORE FIRST USE
#   A. Copy this file to the TOP-LEVEL Frehg2 repository directory, i.e.
#      the directory that contains CMakeLists.txt, src/, benchmarks/, and
#      scripts/.  Do NOT put it in build/.
#   B. Log in to the cluster and, from that directory, submit:
#         sbatch build_frehg2_slurm.sh
#   C. Watch progress:
#         squeue -u "$USER"
#         tail -f frehg2-build-<job-id>.out
#      On success, the executable is: build/src/frehg
#
# IMPORTANT SITE-SPECIFIC SETTINGS
#   The automatic module search is intentionally conservative. If it cannot
#   select correct modules, first inspect the names supplied by your centre:
#         module spider cmake gcc openmpi mpich hdf5
#      Then submit with the exact full module names, for example:
#         sbatch --export=ALL,FREHG_MODULE_CMAKE=cmake/3.28.3,\
#         FREHG_MODULE_COMPILER=gcc/13.2.0,FREHG_MODULE_MPI=openmpi/4.1.6,\
#         FREHG_MODULE_HDF5=hdf5/1.14.3-openmpi build_frehg2_slurm.sh
#   Or edit the four FREHG_MODULE_* variables immediately below. Do not use
#   an MPI wrapper (mpicxx) as Frehg2's CMAKE_CXX_COMPILER: this script uses
#   the underlying real gcc/g++ (or clang/clang++) instead.
#
#   A build needs outbound HTTPS access to download source dependencies.
#   If compute nodes have no Internet access, use a centre-provided HDF5,
#   Kokkos, PETSc and yaml-cpp module, or ask the support team how to obtain
#   approved source tarballs / a build allocation with Internet access.
#
# OPTIONAL OVERRIDES
#   FREHG_PREFIX=/path/to/user/prefix       dependency installation location
#   FREHG_MODULE_CMAKE=cmake/3.27.9         exact full module name
#   FREHG_MODULE_COMPILER=gcc/13.2.0        exact full module name
#   FREHG_MODULE_MPI=openmpi/4.1.6          exact full module name
#   FREHG_MODULE_HDF5=hdf5/1.14.3-openmpi   exact full module name
#   FREHG_JOBS=4                            number of simultaneous compilers
#   FREHG_ENABLE_TESTS=1                    also configure network-fetched tests
#   FREHG_WERROR=OFF                        only if a new compiler gives warnings
#
# This is a CPU build with Kokkos Serial+OpenMP host backends and a
# Kokkos-aware PETSc, so both the physics kernels and the linear solve
# thread under OMP_NUM_THREADS (mat_type aijkokkos). Frehg2's CUDA/HIP
# path ships experimental and GPU-unverified; use a separate, site- and
# GPU-specific device workflow for that.

set -Eeuo pipefail
shopt -s nullglob
trap 'rc=$?; echo "ERROR: failed at line $LINENO (exit $rc). Read the log above; no successful build was claimed." >&2; exit "$rc"' ERR

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX="${FREHG_PREFIX:-$HOME/frehg-deps}"
SRC_CACHE="$PREFIX/src-cache"
BUILD_ROOT="$ROOT_DIR/build"
JOBS="${FREHG_JOBS:-${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-4}}}"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || JOBS=4

FREHG_MODULE_CMAKE="${FREHG_MODULE_CMAKE:-}"
FREHG_MODULE_COMPILER="${FREHG_MODULE_COMPILER:-}"
FREHG_MODULE_MPI="${FREHG_MODULE_MPI:-}"
FREHG_MODULE_HDF5="${FREHG_MODULE_HDF5:-}"

say() { printf '\n========== %s ==========%s\n' "$*" ""; }
die() { echo "ERROR: $*" >&2; exit 1; }
need_command() { command -v "$1" >/dev/null 2>&1 || die "Required command '$1' is unavailable."; }
version_ge() { # true if $1 >= $2; works for numeric dotted versions
  [[ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" == "$2" ]]
}

[[ -f "$ROOT_DIR/CMakeLists.txt" ]] || die "Run this script from the top-level Frehg2 repository directory. CMakeLists.txt was not found beside this script."
mkdir -p "$PREFIX" "$SRC_CACHE"

say "1 of 8: Record environment and locate module command"
printf 'Host: %s\nDate: %s\nRepository: %s\nPrefix: %s\nSLURM job: %s\n' "$(hostname)" "$(date -Is)" "$ROOT_DIR" "$PREFIX" "${SLURM_JOB_ID:-not submitted through sbatch}"
if ! type module >/dev/null 2>&1; then
  # A few clusters require sourcing this only in non-interactive jobs.
  [[ -r /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  [[ -r /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
fi
type module >/dev/null 2>&1 || die "The environment-modules command is unavailable. Ask your HPC support team which shell initialization file enables modules in SLURM jobs."
module purge

# Return the first complete module name matching a regular expression. The
# module command commonly writes its list to stderr, hence 2>&1.
find_module() {
  local regex="$1" x
  while IFS= read -r x; do
    x="${x%%(*}"; x="${x// /}"
    [[ "$x" =~ ^$regex$ ]] && { printf '%s\n' "$x"; return 0; }
  done < <(module -t avail 2>&1 | sed '/^$/d' | sort -Vr)
  return 1
}
load_exact_or_find() {
  local variable="$1" requested="$2" regex="$3" selected=""
  if [[ -n "$requested" ]]; then
    selected="$requested"
  else
    selected="$(find_module "$regex" || true)"
  fi
  [[ -n "$selected" ]] || die "Could not find a suitable $variable module automatically. Run: module spider cmake gcc openmpi mpich hdf5 ; then resubmit with FREHG_MODULE_${variable}=<exact-module-name>."
  echo "Loading $variable module: $selected"
  module load "$selected" || die "Could not load $variable module '$selected'. Supply its exact full module name through FREHG_MODULE_${variable}."
}

# Module names vary by site. These patterns cover usual Module/Lmod names,
# but explicit FREHG_MODULE_* values above are always preferred.
load_exact_or_find CMAKE "$FREHG_MODULE_CMAKE" '(cmake|CMake)/[0-9].*'
load_exact_or_find COMPILER "$FREHG_MODULE_COMPILER" '(gcc|gnu|GCC|clang|llvm)/[0-9].*'
load_exact_or_find MPI "$FREHG_MODULE_MPI" '(openmpi|OpenMPI|mpi/openmpi|mpich|MPICH|intel-mpi|impi)/[0-9].*'

need_command cmake
need_command mpicc
need_command mpicxx
need_command curl
CMAKE_VERSION="$(cmake --version | awk 'NR==1{print $3}')"
version_ge "$CMAKE_VERSION" 3.23 || die "CMake $CMAKE_VERSION is too old; Frehg2 needs CMake >= 3.23. Load a newer CMake module."

# Deliberately derive real compiler paths before PETSc builds. Frehg2 must
# not use mpicc/mpicxx as its CMake compiler.
if command -v g++ >/dev/null 2>&1 && [[ "$(g++ --version 2>&1 | head -1)" =~ (GCC|g\+\+) ]]; then
  REAL_CXX="$(command -v g++)"; REAL_CC="$(command -v gcc)"
elif command -v clang++ >/dev/null 2>&1; then
  REAL_CXX="$(command -v clang++)"; REAL_CC="$(command -v clang)"
else
  die "No real gcc/g++ or clang/clang++ was available after loading the compiler module."
fi
COMPILER_VERSION="$($REAL_CXX --version | head -1)"
echo "CMake: $CMAKE_VERSION"
echo "Real C compiler: $REAL_CC"
echo "Real C++ compiler: $REAL_CXX ($COMPILER_VERSION)"
echo "MPI compiler wrapper: $(command -v mpicxx)"
mpiexec --version 2>&1 | head -5 || true

say "2 of 8: Load and verify a parallel HDF5 module if one exists"
if [[ -n "$FREHG_MODULE_HDF5" ]]; then
  echo "Loading requested HDF5 module: $FREHG_MODULE_HDF5"
  module load "$FREHG_MODULE_HDF5"
else
  HDF_CANDIDATE="$(find_module '(hdf5|HDF5|parallel-hdf5|hdf5-parallel)/[0-9].*' || true)"
  if [[ -n "$HDF_CANDIDATE" ]]; then
    echo "Trying HDF5 module: $HDF_CANDIDATE"
    module load "$HDF_CANDIDATE" || true
  else
    echo "No HDF5 module found; parallel HDF5 will be built below."
  fi
fi

USE_MODULE_HDF5=0
if command -v h5pcc >/dev/null 2>&1 && h5pcc -showconfig 2>/dev/null | grep -Eqi 'Parallel HDF5:[[:space:]]+yes'; then
  USE_MODULE_HDF5=1
  echo "Usable parallel HDF5 detected: $(command -v h5pcc)"
else
  echo "A usable parallel HDF5 was not detected; building HDF5 1.14.6 into $PREFIX."
fi

say "3 of 8: Define safe source-build helpers"
fetch_tarball() {
  local url="$1" archive="$2" directory="$3"
  if [[ ! -d "$directory" ]]; then
    [[ -f "$archive" ]] || curl --fail --location --retry 3 --output "$archive" "$url"
    tar -xf "$archive" -C "$SRC_CACHE"
  fi
}
cmake_build_install() {
  local source="$1" build="$2"; shift 2
  cmake -S "$source" -B "$build" -G "Unix Makefiles" "$@"
  cmake --build "$build" --parallel "$JOBS"
  cmake --install "$build"
}

say "4 of 8: Build yaml-cpp and Kokkos in your user prefix"
# Building these locally guarantees the same compiler family as Frehg2.
YAML_SRC="$SRC_CACHE/yaml-cpp-0.8.0"
fetch_tarball "https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz" "$SRC_CACHE/yaml-cpp-0.8.0.tar.gz" "$YAML_SRC"
if [[ ! -f "$PREFIX/lib/cmake/yaml-cpp/yaml-cpp-config.cmake" && ! -f "$PREFIX/lib64/cmake/yaml-cpp/yaml-cpp-config.cmake" ]]; then
  cmake_build_install "$YAML_SRC" "$SRC_CACHE/build-yaml-cpp-0.8.0" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_C_COMPILER="$REAL_CC" -DCMAKE_CXX_COMPILER="$REAL_CXX" \
    -DYAML_BUILD_SHARED_LIBS=ON -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
else
  echo "yaml-cpp already installed in $PREFIX; reusing it."
fi

KOKKOS_SRC="$SRC_CACHE/kokkos-5.1.1"
fetch_tarball "https://github.com/kokkos/kokkos/archive/refs/tags/5.1.1.tar.gz" "$SRC_CACHE/kokkos-5.1.1.tar.gz" "$KOKKOS_SRC"
if [[ ! -f "$PREFIX/lib/cmake/Kokkos/KokkosConfig.cmake" && ! -f "$PREFIX/lib64/cmake/Kokkos/KokkosConfig.cmake" ]]; then
  # BUILD_SHARED_LIBS=ON is load-bearing: frehg links Kokkos and so does
  # PETSc's downloaded kokkos-kernels. A static Kokkos core is absorbed
  # into BOTH the frehg executable and libkokkoskernels, giving two copies
  # of Kokkos's runtime singleton in one process -- it initializes twice
  # ("Kokkos::OpenMP::initialize" prints twice) and the first VecKokkos
  # access segfaults. One shared libkokkoscore keeps a single runtime.
  cmake_build_install "$KOKKOS_SRC" "$SRC_CACHE/build-kokkos-5.1.1" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_C_COMPILER="$REAL_CC" -DCMAKE_CXX_COMPILER="$REAL_CXX" \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_CXX_STANDARD=20 -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON \
    -DKokkos_ENABLE_TESTS=OFF
else
  echo "Kokkos already installed in $PREFIX; reusing it."
fi

say "5 of 8: Obtain a matching parallel HDF5"
if [[ "$USE_MODULE_HDF5" -eq 0 ]]; then
  HDF5_SRC="$SRC_CACHE/hdf5-hdf5_1.14.6"
  fetch_tarball "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.6.tar.gz" "$SRC_CACHE/hdf5-1.14.6.tar.gz" "$HDF5_SRC"
  if [[ ! -x "$PREFIX/bin/h5pcc" ]]; then
    cmake_build_install "$HDF5_SRC" "$SRC_CACHE/build-hdf5-1.14.6" \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_C_COMPILER="$(command -v mpicc)" -DCMAKE_CXX_COMPILER="$(command -v mpicxx)" \
      -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_CPP_LIB=ON -DHDF5_BUILD_TOOLS=ON \
      -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TESTING=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF
  else
    echo "User-built HDF5 already installed in $PREFIX; reusing it."
  fi
  export PATH="$PREFIX/bin:$PATH"
  export HDF5_ROOT="$PREFIX"
fi
h5pcc -showconfig 2>/dev/null | grep -Eqi 'Parallel HDF5:[[:space:]]+yes' || die "HDF5 is not parallel. Frehg2 requires an MPI-enabled HDF5 built with the currently loaded MPI."

say "6 of 8: Build PETSc in your user prefix"
PETSC_SRC="$SRC_CACHE/petsc-3.25.1"
fetch_tarball "https://gitlab.com/petsc/petsc/-/archive/v3.25.1/petsc-v3.25.1.tar.gz" "$SRC_CACHE/petsc-3.25.1.tar.gz" "$PETSC_SRC"
if [[ ! -f "$PREFIX/lib/petsc/conf/petscvariables" && ! -f "$PREFIX/lib64/petsc/conf/petscvariables" ]]; then
  pushd "$PETSC_SRC" >/dev/null
  # PETSc correctly uses MPI wrappers for ITS dependency build. This is not
  # the Frehg2 CMake compiler setting, which remains the real compiler.
  # --with-kokkos-dir + --download-kokkos-kernels reuse the shared Kokkos
  # built above so the free-surface/groundwater solves run on the Kokkos
  # backend and thread under OMP_NUM_THREADS (mat_type aijkokkos, v2 Q3);
  # --with-openmp also threads the downloaded hypre BoomerAMG (solver amg).
  ./configure --prefix="$PREFIX" --with-cc="$(command -v mpicc)" --with-cxx="$(command -v mpicxx)" \
    --with-fc=0 --with-debugging=0 --download-f2cblaslapack --download-hypre \
    --with-kokkos-dir="$PREFIX" --download-kokkos-kernels --with-openmp=1
  make -j "$JOBS" all
  make install
  popd >/dev/null
else
  echo "PETSc already installed in $PREFIX; reusing it."
fi
grep -q "PETSC_HAVE_HYPRE" "$PREFIX/include/petscconf.h" 2>/dev/null \
  || grep -q "PETSC_HAVE_HYPRE" "$PREFIX"/lib*/petsc/conf/petscvariables 2>/dev/null \
  || die "PETSc was installed without hypre (solver 'amg' needs it). Delete '$PREFIX/lib/petsc' and '$PREFIX/src-cache/petsc-3.25.1', then rerun."
grep -q "PETSC_HAVE_KOKKOS_KERNELS" "$PREFIX/include/petscconf.h" 2>/dev/null \
  || die "PETSc was installed without Kokkos Kernels (mat_type aijkokkos needs it). Delete '$PREFIX/lib/petsc' and '$PREFIX/src-cache/petsc-3.25.1', then rerun to rebuild PETSc against the shared Kokkos."
# A reused PETSc must still have its Kokkos Kernels runtime present:
# petscconf.h advertises KOKKOS_KERNELS even after the library is gone.
# This catches "Kokkos rebuilt but PETSc reused" -- the old libkokkoskernels
# is removed with the old Kokkos, libpetsc dangles, and the failure is
# otherwise a cryptic loader abort at first run instead of an actionable
# message here.
[[ -f "$PREFIX/lib/libkokkoskernels.so" || -f "$PREFIX/lib64/libkokkoskernels.so" \
   || -f "$PREFIX/lib/libkokkoskernels.dylib" || -f "$PREFIX/lib64/libkokkoskernels.dylib" ]] \
  || die "PETSc's Kokkos Kernels runtime is missing from $PREFIX (Kokkos was likely rebuilt without PETSc). Delete '$PREFIX/lib/petsc' and '$PREFIX/src-cache/petsc-3.25.1', then rerun to rebuild PETSc against the current Kokkos."

say "7 of 8: Configure and compile Frehg2"
export CMAKE_PREFIX_PATH="$PREFIX${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PREFIX/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# Start fresh: a CMake cache created using another compiler or MPI is unsafe.
rm -rf "$BUILD_ROOT"
TESTS="${FREHG_ENABLE_TESTS:-0}"
WERROR="${FREHG_WERROR:-ON}"
env CC="$REAL_CC" CXX="$REAL_CXX" cmake -S "$ROOT_DIR" -B "$BUILD_ROOT" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" \
  -DFREHG_ENABLE_TESTS="$TESTS" -DFREHG_WERROR="$WERROR"
cmake --build "$BUILD_ROOT" --parallel "$JOBS"
EXE="$BUILD_ROOT/src/frehg"
[[ -x "$EXE" ]] || die "CMake finished but expected executable $EXE was not produced."
echo "SUCCESS: executable created: $EXE"

say "8 of 8: Validate and perform a short scheduled verification"
export OMP_NUM_THREADS=1
"$EXE" --validate "$ROOT_DIR/benchmarks/b1-sw/b1-sw.yaml"
# srun uses this SLURM allocation, not an unmanaged login-node process.
pushd "$ROOT_DIR/benchmarks/b1-sw" >/dev/null
srun --ntasks="${SLURM_NTASKS:-1}" --cpus-per-task=1 "$EXE" b1-sw.yaml
popd >/dev/null

echo ""
echo "================================================================="
echo "FREHG2 BUILD AND VERIFICATION COMPLETED SUCCESSFULLY"
echo "Executable: $EXE"
echo "Dependency prefix: $PREFIX"
echo "For future run jobs, load exactly the module list shown below and"
echo "set CMAKE_PREFIX_PATH/PATH/LD_LIBRARY_PATH to the prefix if needed:"
module list 2>&1 || true
echo "================================================================="
