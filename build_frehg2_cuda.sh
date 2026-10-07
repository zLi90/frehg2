#!/usr/bin/env bash
#SBATCH --job-name=frehg2-cuda-build
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=frehg2-cuda-build-%j.out
#SBATCH --error=frehg2-cuda-build-%j.err
#
# Frehg2 automated GPU build script (Kokkos CUDA backend) -- interactive or SLURM
#
# >>> "--partition=gpu" above is a PLACEHOLDER. Replace it with your centre's
# >>> GPU partition (list them: sinfo -o '%P %G %l' | grep -i gpu) or override
# >>> it when submitting: sbatch --partition=<gpu-partition> build_frehg2_cuda.sh
# >>> (some sites spell the GPU request --gpus=1 instead of --gres=gpu:1).
# >>> FREHG_GPU_PARTITION=<name> only fills in the instructions this script
# >>> prints; #SBATCH lines cannot read shell variables.
#
# WHAT THIS SCRIPT DOES
#   1. Refuses to run unless an NVIDIA GPU is visible (nvidia-smi): it must
#      run ON A GPU NODE, not on the login node or a CPU (partition intel) node.
#   2. Loads CMake and gcc modules. Defaults are this site's CPU stack
#      (cmake/3.31.6 gcc/13.2.0, later openmpi/4.1.6 hdf5/1.14.6).
#   3. Picks a CUDA toolkit. The CUDA module name is unknown, so candidates
#      (cuda/12.x|13.x, CUDA/, cudatoolkit/, nvhpc last) are tried best-first:
#      each must have nvcc >= 12.2 (Kokkos 5.1.1's floor), compile C++20 with
#      the loaded gcc, run a small CUDA runtime test program on this GPU under
#      this driver, and target this GPU's arch; the first that passes wins.
#      FREHG_MODULE_CUDA forces one. Kokkos_ARCH_* comes from the GPU's
#      compute capability (7.0 VOLTA70, 8.0 AMPERE80, 8.6 AMPERE86, 8.9 ADA89,
#      9.0 HOPPER90, ...). FREHG_CUDA_ARCH overrides.
#   4. Loads MPI and reports whether it is CUDA-aware and what that means for
#      frehg2's runtime.gpu_aware_mpi and PETSc's -use_gpu_aware_mpi.
#   5. Builds into a SEPARATE prefix (default $HOME/frehg-deps-cuda): yaml-cpp
#      0.8.0, Kokkos 5.1.1 (SHARED; CUDA + Serial + OpenMP), parallel HDF5
#      1.14.6 if no module provides one, and PETSc 3.25.1 configured with CUDA,
#      --with-kokkos-dir=<that same Kokkos>, --download-kokkos-kernels, hypre
#      built for CUDA, --with-openmp and the --with-make-np OOM guard.
#   6. Builds frehg2 through Kokkos' nvcc_wrapper with -DFREHG_BACKEND=cuda in
#      build-cuda/ (binary build-cuda/src/frehg). The CPU build (build/ and
#      $HOME/frehg-deps) is never touched.
#   7. Validates benchmarks/b1-sw and runs it on 1 rank / 1 GPU (in a copy
#      under build-cuda/), checks the "Kokkos backend Cuda" banner and the run
#      record, runs a secondary b1-sw with hypre BoomerAMG on the device, and
#      prints the run-time environment.
#
# STATUS: WHAT IS NOT VERIFIED (read before relying on it)
#   - This script has NEVER been executed. It was written on a machine with no
#     GPU and no access to the cluster; only "bash -n" and shellcheck were run.
#   - Not known when written: the GPU partition name, the CUDA module name and
#     version, the GPU model, whether GPU nodes have Internet access, whether
#     openmpi/4.1.6 is CUDA-aware, whether hdf5/1.14.6 is parallel, and whether
#     Open MPI's hwloc works on the GPU nodes (it crashed on partition intel
#     on 2026-09-01; step 10 retries once with the workaround). Hint: frehg2's
#     first nvcc compile (2026-08-27) was on the owner's GPU cluster with CUDA
#     13.0 + gcc 13.2 and L40 GPUs (ADA89). If that is this cluster, expect a
#     CUDA 13.0 module and ADA89; the script still detects both.
#   - Never compiled together anywhere: Kokkos 5.1.1 (CUDA) + the Kokkos
#     Kernels 5.1.0 that PETSc 3.25.1 downloads. Kokkos Kernels 5.1.0's CMake
#     accepts Kokkos 5.0.2-5.1.99, and build_frehg2_slurm.sh uses the same
#     5.1.1/5.1.0 pair for the host backend, but the cuda-compile CI lane
#     (the only CUDA build that has compiled frehg2) used Kokkos 5.1.0 + Kokkos
#     Kernels 5.1.0. If Kokkos Kernels fails to compile, wipe the prefix (see
#     the message) and rerun with FREHG_KOKKOS_VERSION=5.1.0 (exact lockstep).
#   - Frehg2's GPU lane itself is EXPERIMENTAL (docs/agents/build.md section
#     7, docs/developer-guide/gpu-acceptance.md): compile/link-verified in CI,
#     device execution unverified until the owner's p6 acceptance bundle
#     (scripts/gpu_acceptance.sh) passes. A green run here shows the binary
#     runs on the GPU. It does NOT validate GPU physics.
#
# HOW TO USE IT (interactive GPU node; this is the normal route)
#   A. This file sits at the TOP LEVEL of the Frehg2 repository, beside
#      CMakeLists.txt and build_frehg2_slurm.sh.
#   B. Only if GPU nodes have no Internet: on the LOGIN node run once
#         bash build_frehg2_cuda.sh --fetch-only
#      (downloads every source archive into the prefix's src-cache).
#   C. Start tmux/screen on the login node (an ssh drop kills an srun shell),
#      then get a shell on a GPU node, for example:
#         srun --partition=<gpu-partition> --gres=gpu:1 --cpus-per-task=8 \
#              --mem=64G --time=06:00:00 --pty bash
#   D. On the GPU node, from the repository top level:
#         nvidia-smi                       # must list a GPU
#         bash build_frehg2_cuda.sh 2>&1 | tee frehg2-cuda-build.log
#   E. On success the executable is build-cuda/src/frehg and the script prints
#      the exact modules and environment that run jobs need.
#   Batch alternative: sbatch --partition=<gpu-partition> build_frehg2_cuda.sh
#   Re-running is safe: finished dependencies in the prefix are reused, and
#   build-cuda/ is rebuilt from scratch every time.
#   Expect 1.5-3 h on 8 cores the first time. Kokkos Kernels and hypre take
#   longest under nvcc.
#
# OPTIONAL OVERRIDES (environment variables)
#   FREHG_PREFIX=$HOME/frehg-deps-cuda   dependency prefix (NOT the CPU one).
#                                        One prefix = one MPI: to build against
#                                        another MPI module, give it a new
#                                        prefix (the script stops otherwise)
#   FREHG_BUILD_NAME=build-cuda          frehg2 build directory inside the repo
#                                        (build-<name>; rebuilt from scratch);
#                                        a second name keeps the first binary
#   FREHG_MODULE_CMAKE=cmake/3.31.6      exact module names; without them the
#   FREHG_MODULE_COMPILER=gcc/13.2.0     site defaults shown are used if
#   FREHG_MODULE_CUDA=cuda/12.6          available, otherwise discovered.
#   FREHG_MODULE_MPI=openmpi/4.1.6
#   FREHG_MODULE_HDF5=hdf5/1.14.6
#   FREHG_CUDA_ARCH=ADA89                Kokkos arch name, or 8.9 / 89
#   FREHG_GPU_PARTITION=<name>           used only in printed instructions
#   FREHG_JOBS=8                         compile parallelism
#   FREHG_PETSC_MAKE_NP=2                PETSc package-build parallelism (the
#                                        default caps it at ~8 GB RAM per job)
#   FREHG_MPI_GPU_AWARE=auto|0|1         override the CUDA-aware-MPI detection
#   FREHG_GPU_AWARE_MPI=OFF|ON           frehg2 compile-time default behind
#                                        runtime.gpu_aware_mpi: auto
#   FREHG_KOKKOS_VERSION=5.1.1           5.1.0 = PETSc's Kokkos Kernels pin
#   FREHG_PETSC_EXTRA="--download-umpire"  extra PETSc configure arguments
#   FREHG_ENABLE_TESTS=1                 also configure network-fetched tests
#   FREHG_WERROR=OFF                     only if a new nvcc/gcc gives warnings
#
# SOURCES BEHIND THE DEVICE SETTINGS (beyond build_frehg2_slurm.sh)
#   frehg2: CMakeLists.txt (FREHG_BACKEND check, PETSC_HAVE_KOKKOS requirement,
#     LTO scrub for nvcc), docs/agents/build.md s.7, docs/user-guide/
#     installation.md, docs/developer-guide/gpu-acceptance.md,
#     .github/workflows/cuda-compile.yml, scripts/ci_install_deps.sh
#     (FREHG_CUDA=1 recipe), src/core/LinearSystem.cpp (aijkokkos forced on
#     device), src/core/PetscSession.cpp (startup banner), src/io/RunRecord.cpp.
#   Kokkos 5.1.1 source (https://github.com/kokkos/kokkos/tree/5.1.1):
#     cmake/kokkos_compiler_id.cmake (KOKKOS_NVCC_MINIMUM 12.2.0),
#     cmake/kokkos_enable_options.cmake (CUDA_LAMBDA, CUDA_CONSTEXPR,
#     IMPL_CUDA_MALLOC_ASYNC, RDC needs static libs), cmake/kokkos_arch.cmake
#     (arch names), core/src/impl/Kokkos_Core.cpp (rank->GPU mapping).
#   PETSc 3.25.1 source (https://gitlab.com/petsc/petsc/-/tree/v3.25.1):
#     config/BuildSystem/config/packages/{kokkos,kokkos-kernels,hypre,CUDA,
#     MPI}.py, src/sys/objects/init.c and
#     src/sys/objects/device/impls/cupm/cupmdevice.cxx (-use_gpu_aware_mpi).
#   https://petsc.org/release/install/install/ (GPU section)
#   https://hypre.readthedocs.io/en/latest/ch-misc.html (GPU build options)
#   https://www.open-mpi.org/faq/?category=runcuda (ompi_info CUDA check)
#   https://docs.nvidia.com/cuda/archive/13.0.1/cuda-toolkit-release-notes/
#     (CUDA 13 dropped Maxwell/Pascal/Volta, i.e. compute capability < 7.5)
#   https://docs.nvidia.com/cuda/cuda-installation-guide-linux/ (host gcc table)

set -Eeuo pipefail
shopt -s nullglob
on_error() {
  local rc=$? line="$1"
  echo "ERROR: failed at line $line (exit $rc). Read the log above; no successful build was claimed." >&2
  exit "$rc"
}
trap 'on_error "$LINENO"' ERR

# Under sbatch, SLURM runs a spooled copy of this file, so its directory is not
# the repository. Fall back to the submission directory in that case.
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -f "$ROOT_DIR/CMakeLists.txt" && -n "${SLURM_SUBMIT_DIR:-}" && -f "$SLURM_SUBMIT_DIR/CMakeLists.txt" ]]; then
  ROOT_DIR="$SLURM_SUBMIT_DIR"
fi
PREFIX="${FREHG_PREFIX:-$HOME/frehg-deps-cuda}"
PREFIX="${PREFIX%/}"
SRC_CACHE="$PREFIX/src-cache"
BUILD_NAME="${FREHG_BUILD_NAME:-build-cuda}"
# It is removed and rebuilt below, so only a plain build-<name> directory in
# the repository is accepted (never the CPU build's build/).
[[ "$BUILD_NAME" =~ ^build-[A-Za-z0-9._-]+$ ]] \
  || { echo "ERROR: FREHG_BUILD_NAME='$BUILD_NAME' must look like build-<name> (letters, digits, . _ -), e.g. build-cuda-gpumpi." >&2; exit 1; }
BUILD_ROOT="$ROOT_DIR/$BUILD_NAME"
EXE="$BUILD_ROOT/src/frehg"
JOBS="${FREHG_JOBS:-${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-4}}}}"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || JOBS=4
GPU_PARTITION_HINT="${FREHG_GPU_PARTITION:-<gpu-partition>}"

# Version pins (identical to build_frehg2_slurm.sh unless overridden).
YAML_VERSION=0.8.0
KOKKOS_VERSION="${FREHG_KOKKOS_VERSION:-5.1.1}"
HDF5_VERSION=1.14.6
PETSC_VERSION=3.25.1
KOKKOS_NVCC_MIN=12.2.0   # Kokkos 5.1.1 cmake/kokkos_compiler_id.cmake
# What PETSc 3.25.1's own --download-* would fetch (its packages/*.py). The
# script fetches them itself so a GPU node without Internet can still build,
# after --fetch-only has run on the login node.
PETSC_KK_VERSION=5.1.0
PETSC_HYPRE_COMMIT=2395097204558c0f65e110302e055000c24238ff   # hypre v3.1.0 + PR 1487
PETSC_F2C=f2cblaslapack-3.8.0.q2

# This site's CPU-build modules (nueces/submitfrehg-*.sh). Explicit
# FREHG_MODULE_* values win, then these if available, then discovery.
SITE_MODULE_CMAKE=cmake/3.31.6
SITE_MODULE_COMPILER=gcc/13.2.0
SITE_MODULE_MPI=openmpi/4.1.6
SITE_MODULE_HDF5=hdf5/1.14.6
FREHG_MODULE_CMAKE="${FREHG_MODULE_CMAKE:-}"
FREHG_MODULE_COMPILER="${FREHG_MODULE_COMPILER:-}"
FREHG_MODULE_CUDA="${FREHG_MODULE_CUDA:-}"
FREHG_MODULE_MPI="${FREHG_MODULE_MPI:-}"
FREHG_MODULE_HDF5="${FREHG_MODULE_HDF5:-}"

say() { printf '\n========== %s ==========\n' "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }
need_command() { command -v "$1" >/dev/null 2>&1 || die "Required command '$1' is unavailable."; }
version_ge() { # true if $1 >= $2; works for numeric dotted versions
  [[ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" == "$2" ]]
}
ver_mm() { # major.minor of a dotted version
  local IFS=. a b rest
  read -r a b rest <<<"$1"
  printf '%s.%s\n' "$a" "${b:-0}"
}
first_existing() {
  local f
  for f in "$@"; do [[ -e "$f" ]] && { printf '%s\n' "$f"; return 0; }; done
  return 1
}
usage() {
  awk 'NR == 1 || /^#SBATCH/ { next } /^#/ { sub(/^# ?/, ""); print; next } { exit }' "${BASH_SOURCE[0]}"
}

FETCH_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --fetch-only) FETCH_ONLY=1 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument '$arg' (supported: --fetch-only, --help)." ;;
  esac
done

[[ -f "$ROOT_DIR/CMakeLists.txt" ]] || die "Run this script from the top-level Frehg2 repository directory. CMakeLists.txt was not found beside this script."
[[ "$PREFIX" != "$HOME/frehg-deps" ]] || die "FREHG_PREFIX=$PREFIX is the CPU build's prefix. The CUDA stack must live in its own prefix (default \$HOME/frehg-deps-cuda) so the CPU build is not clobbered."
mkdir -p "$PREFIX" "$SRC_CACHE"

# ---- download helpers (used by --fetch-only and by the build) --------------
fetch_file() { # fetch_file <url> <destination>; partial downloads never count
  local url="$1" dest="$2"
  [[ -s "$dest" ]] && return 0
  need_command curl
  echo "Downloading $url"
  if ! curl --fail --location --retry 3 --connect-timeout 30 --output "$dest.part" "$url"; then
    rm -f "$dest.part"
    die "Download failed: $url. If this node has no Internet access, run 'bash build_frehg2_cuda.sh --fetch-only' on a login node first (same FREHG_PREFIX), then rerun here."
  fi
  mv "$dest.part" "$dest"
}
fetch_tarball() { # fetch_tarball <url> <archive> <directory>
  local url="$1" archive="$2" directory="$3"
  [[ -d "$directory" ]] && return 0
  fetch_file "$url" "$archive"
  # Extract into a scratch dir and move the single top-level directory to the
  # expected name (GitLab tag archives name it <project>-<tag>-<sha>).
  local scratch="$SRC_CACHE/.extract-$$" entry entries=()
  rm -rf "$scratch"; mkdir -p "$scratch"
  tar -xf "$archive" -C "$scratch"
  while IFS= read -r entry; do entries+=("$entry"); done \
    < <(find "$scratch" -mindepth 1 -maxdepth 1 ! -name '.DS_Store')
  if [[ ${#entries[@]} -eq 1 && -d "${entries[0]}" ]]; then
    mv "${entries[0]}" "$directory"
    rm -rf "$scratch"
  else
    mv "$scratch" "$directory"
  fi
  [[ -d "$directory" ]] || die "Could not extract $archive into $directory."
}
cmake_build_install() {
  local source="$1" build="$2"; shift 2
  cmake -S "$source" -B "$build" -G "Unix Makefiles" "$@"
  cmake --build "$build" --parallel "$JOBS"
  cmake --install "$build"
}

URL_YAML="https://github.com/jbeder/yaml-cpp/archive/refs/tags/${YAML_VERSION}.tar.gz"
URL_KOKKOS="https://github.com/kokkos/kokkos/archive/refs/tags/${KOKKOS_VERSION}.tar.gz"
URL_HDF5="https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_${HDF5_VERSION}.tar.gz"
URL_PETSC="https://gitlab.com/petsc/petsc/-/archive/v${PETSC_VERSION}/petsc-v${PETSC_VERSION}.tar.gz"
URL_KK="https://github.com/kokkos/kokkos-kernels/archive/${PETSC_KK_VERSION}.tar.gz"
URL_HYPRE="https://github.com/hypre-space/hypre/archive/${PETSC_HYPRE_COMMIT}.tar.gz"
URL_F2C="https://web.cels.anl.gov/projects/petsc/download/externalpackages/${PETSC_F2C}.tar.gz"
ARC_YAML="$SRC_CACHE/yaml-cpp-${YAML_VERSION}.tar.gz"
ARC_KOKKOS="$SRC_CACHE/kokkos-${KOKKOS_VERSION}.tar.gz"
ARC_HDF5="$SRC_CACHE/hdf5-${HDF5_VERSION}.tar.gz"
ARC_PETSC="$SRC_CACHE/petsc-${PETSC_VERSION}.tar.gz"
PKG_KK="$SRC_CACHE/petsc-pkg-kokkos-kernels-${PETSC_KK_VERSION}.tar.gz"
PKG_HYPRE="$SRC_CACHE/petsc-pkg-hypre-${PETSC_HYPRE_COMMIT}.tar.gz"
PKG_F2C="$SRC_CACHE/petsc-pkg-${PETSC_F2C}.tar.gz"

if [[ "$FETCH_ONLY" -eq 1 ]]; then
  say "Fetch-only mode: downloading all sources into $SRC_CACHE"
  fetch_file "$URL_YAML" "$ARC_YAML"
  fetch_file "$URL_KOKKOS" "$ARC_KOKKOS"
  fetch_file "$URL_HDF5" "$ARC_HDF5"
  fetch_file "$URL_PETSC" "$ARC_PETSC"
  fetch_file "$URL_KK" "$PKG_KK"
  fetch_file "$URL_HYPRE" "$PKG_HYPRE"
  fetch_file "$URL_F2C" "$PKG_F2C"
  echo "All sources are cached. Now get a GPU node and run: bash build_frehg2_cuda.sh"
  exit 0
fi

say "1 of 10: Record environment and verify a GPU is visible"
printf 'Host: %s\nDate: %s\nRepository: %s\nPrefix: %s\nBuild dir: %s\nSLURM job: %s\nJobs: %s\n' \
  "$(hostname)" "$(date '+%Y-%m-%dT%H:%M:%S%z')" "$ROOT_DIR" "$PREFIX" "$BUILD_ROOT" \
  "${SLURM_JOB_ID:-not inside a SLURM allocation}" "$JOBS"
no_gpu_abort() {
  cat >&2 <<EOF
ERROR: no NVIDIA GPU is visible on $(hostname): $1
This looks like a login node or a CPU-only node, or a GPU job that did not
request a GPU. The CUDA build must run ON A GPU NODE: Kokkos' arch, the
driver/toolkit checks and the final GPU verification run all need the device.

Interactive (recommended; start tmux/screen on the login node first):
  srun --partition=${GPU_PARTITION_HINT} --gres=gpu:1 --cpus-per-task=8 --mem=64G --time=06:00:00 --pty bash
  cd "$ROOT_DIR" && bash build_frehg2_cuda.sh 2>&1 | tee frehg2-cuda-build.log
Batch:
  cd "$ROOT_DIR" && sbatch --partition=${GPU_PARTITION_HINT} build_frehg2_cuda.sh
List the GPU partitions:  sinfo -o '%P %G %l' | grep -i gpu
If GPU nodes have no Internet access, first run on this login node:
  bash build_frehg2_cuda.sh --fetch-only
EOF
  exit 1
}
SMI="$(command -v nvidia-smi || true)"
[[ -z "$SMI" && -x /usr/bin/nvidia-smi ]] && SMI=/usr/bin/nvidia-smi
[[ -n "$SMI" ]] || no_gpu_abort "nvidia-smi was not found"
GPU_LIST="$("$SMI" -L 2>&1 || true)"
[[ "$GPU_LIST" =~ GPU\ [0-9]+: ]] || no_gpu_abort "nvidia-smi lists no GPU: $GPU_LIST"
echo "$GPU_LIST"
"$SMI" --query-gpu=index,name,memory.total,driver_version --format=csv,noheader 2>/dev/null || true
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
SMI_OUT="$("$SMI" 2>/dev/null || true)"
DRIVER_CUDA=""
re_driver='CUDA Version: *([0-9]+\.[0-9]+)'
[[ "$SMI_OUT" =~ $re_driver ]] && DRIVER_CUDA="${BASH_REMATCH[1]}"
echo "Highest CUDA version the driver supports: ${DRIVER_CUDA:-unknown}"
CC_LIST="$("$SMI" --query-gpu=compute_cap --format=csv,noheader 2>/dev/null || true)"
GPU_CC="${CC_LIST%%$'\n'*}"; GPU_CC="${GPU_CC//[[:space:]]/}"
[[ "$GPU_CC" =~ ^[0-9]+\.[0-9]+$ ]] || GPU_CC=""
if [[ -n "$GPU_CC" ]]; then
  [[ "$(sort -u <<<"$CC_LIST" | sed '/^[[:space:]]*$/d' | wc -l)" -le 1 ]] \
    || warn "GPUs with different compute capabilities are visible ($(tr '\n' ' ' <<<"$CC_LIST")); building for the first one ($GPU_CC). Set FREHG_CUDA_ARCH to choose."
  echo "Compute capability (nvidia-smi): $GPU_CC"
else
  echo "nvidia-smi did not report compute_cap (old driver?); a CUDA probe will determine it in step 3."
fi
cc_num() { local IFS=. a b; read -r a b <<<"$1"; printf '%s\n' "$(( 10 * a + b ))"; }
GPU_CC_NUM=""; [[ -n "$GPU_CC" ]] && GPU_CC_NUM="$(cc_num "$GPU_CC")"

say "2 of 10: Load the CMake and compiler modules"
if ! type module >/dev/null 2>&1; then
  # A few clusters require sourcing this only in non-interactive jobs.
  # shellcheck disable=SC1091
  [[ -r /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  # shellcheck disable=SC1091
  [[ -r /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
fi
type module >/dev/null 2>&1 || die "The environment-modules command is unavailable. Ask your HPC support team which shell initialization file enables modules on GPU nodes."
module purge

# Return the first complete module name matching a regular expression. The
# module command commonly writes its list to stderr, hence 2>&1. Re-queried
# every call so Lmod hierarchies (openmpi under gcc, ...) are honoured.
find_module() {
  local regex="$1" x
  while IFS= read -r x; do
    x="${x%%(*}"; x="${x// /}"
    [[ "$x" =~ ^$regex$ ]] && { printf '%s\n' "$x"; return 0; }
  done < <(module -t avail 2>&1 | sed '/^$/d' | sort -Vr)
  return 1
}
module_listed() { # exact module name present in "module avail"?
  local want="$1" x
  while IFS= read -r x; do
    x="${x%%(*}"; x="${x// /}"
    [[ "$x" == "$want" ]] && return 0
  done < <(module -t avail "$want" 2>&1)
  return 1
}
LOADED_MODULES=()
load_module() { # load_module <VARIABLE> <requested> <site-default> <regex>
  local variable="$1" requested="$2" site_default="$3" regex="$4" selected=""
  if [[ -n "$requested" ]]; then
    selected="$requested"
  elif [[ -n "$site_default" ]] && module_listed "$site_default"; then
    selected="$site_default"
  else
    selected="$(find_module "$regex" || true)"
  fi
  [[ -n "$selected" ]] || die "Could not find a suitable $variable module automatically. Run: module spider cmake gcc cuda openmpi hdf5 ; then rerun with FREHG_MODULE_${variable}=<exact-module-name>."
  echo "Loading $variable module: $selected"
  module load "$selected" || die "Could not load $variable module '$selected'. Supply its exact full module name through FREHG_MODULE_${variable}."
  LOADED_MODULES+=("$selected")
}

load_module CMAKE "$FREHG_MODULE_CMAKE" "$SITE_MODULE_CMAKE" '(cmake|CMake)/[0-9].*'
load_module COMPILER "$FREHG_MODULE_COMPILER" "$SITE_MODULE_COMPILER" '(gcc|gnu|GCC)/[0-9].*'
need_command cmake
CMAKE_VERSION="$(cmake --version | awk 'NR==1{print $3}')"
version_ge "$CMAKE_VERSION" 3.23 || die "CMake $CMAKE_VERSION is too old; Frehg2 needs CMake >= 3.23. Load a newer CMake module."
# Real compiler paths, never mpicc/mpicxx, for Frehg2/Kokkos (docs/agents/
# build.md invariant 2). Under CUDA the real g++ is nvcc_wrapper's host
# compiler (NVCC_WRAPPER_DEFAULT_COMPILER); nvcc_wrapper is CMake's CXX.
if command -v g++ >/dev/null 2>&1 && [[ "$(g++ --version 2>&1)" =~ (GCC|g\+\+) ]]; then
  REAL_CXX="$(command -v g++)"; REAL_CC="$(command -v gcc)"
else
  die "No real gcc/g++ was available after loading the compiler module (nvcc on Linux needs a GNU host compiler here)."
fi
COMPILER_VERSION="$("$REAL_CXX" --version 2>&1)"; COMPILER_VERSION="${COMPILER_VERSION%%$'\n'*}"
echo "CMake: $CMAKE_VERSION"
echo "Real C compiler: $REAL_CC"
echo "Real C++ compiler (nvcc host compiler): $REAL_CXX ($COMPILER_VERSION)"

say "3 of 10: Select a working CUDA toolkit and the Kokkos GPU arch"
KNOWN_ARCHS="MAXWELL50 MAXWELL52 MAXWELL53 PASCAL60 PASCAL61 VOLTA70 VOLTA72 TURING75 AMPERE80 AMPERE86 AMPERE87 ADA89 HOPPER90 BLACKWELL100 BLACKWELL103 BLACKWELL120 BLACKWELL121"
cc_to_kokkos_arch() { # NVIDIA archs Kokkos 5.1.1 knows (cmake/kokkos_arch.cmake)
  case "$1" in
    5.0) echo MAXWELL50 ;; 5.2) echo MAXWELL52 ;; 5.3) echo MAXWELL53 ;;
    6.0) echo PASCAL60 ;; 6.1) echo PASCAL61 ;;
    7.0) echo VOLTA70 ;; 7.2) echo VOLTA72 ;; 7.5) echo TURING75 ;;
    8.0) echo AMPERE80 ;; 8.6) echo AMPERE86 ;; 8.7) echo AMPERE87 ;; 8.9) echo ADA89 ;;
    9.0) echo HOPPER90 ;;
    10.0) echo BLACKWELL100 ;; 10.3) echo BLACKWELL103 ;;
    12.0) echo BLACKWELL120 ;; 12.1) echo BLACKWELL121 ;;
    *) return 1 ;;
  esac
}
KOKKOS_ARCH=""; ARCH_NUM=""
if [[ -n "${FREHG_CUDA_ARCH:-}" ]]; then
  REQ="$(tr '[:lower:]' '[:upper:]' <<<"$FREHG_CUDA_ARCH")"
  if [[ "$REQ" =~ ^[0-9]+\.[0-9]+$ ]]; then
    KOKKOS_ARCH="$(cc_to_kokkos_arch "$REQ" || true)"
  elif [[ "$REQ" =~ ^[0-9]{2,3}$ ]]; then
    KOKKOS_ARCH="$(cc_to_kokkos_arch "${REQ:0:${#REQ}-1}.${REQ: -1}" || true)"
  else
    KOKKOS_ARCH="$REQ"
  fi
  [[ -n "$KOKKOS_ARCH" && " $KNOWN_ARCHS " == *" $KOKKOS_ARCH "* ]] \
    || die "FREHG_CUDA_ARCH='$FREHG_CUDA_ARCH' is not a Kokkos 5.1.1 NVIDIA arch. Use one of: $KNOWN_ARCHS (or a compute capability such as 8.9)."
  echo "Kokkos arch (FREHG_CUDA_ARCH): $KOKKOS_ARCH"
elif [[ -n "$GPU_CC" ]]; then
  KOKKOS_ARCH="$(cc_to_kokkos_arch "$GPU_CC" || true)"
  [[ -n "$KOKKOS_ARCH" ]] || die "Compute capability $GPU_CC has no Kokkos 5.1.1 arch name. Set FREHG_CUDA_ARCH to one of: $KNOWN_ARCHS."
  echo "Kokkos arch (from compute capability $GPU_CC): $KOKKOS_ARCH"
fi
[[ -n "$KOKKOS_ARCH" ]] && ARCH_NUM="${KOKKOS_ARCH//[!0-9]/}"   # ADA89 -> 89 (PETSc/hypre form)

# CUDA module candidates, best first: the newest the driver supports
# (major.minor <= driver CUDA), then newer same-major ones closest to the
# driver (CUDA minor-version compatibility), then newer majors (need NVIDIA's
# forward-compatibility package). Skipped: nvcc below Kokkos' 12.2 floor, and
# CUDA >= 13 for compute capability < 7.5 (CUDA 13 dropped Maxwell/Pascal/
# Volta). NVIDIA HPC SDK modules come last (the -nompi flavour first).
cuda_module_candidates() {
  local re='^(cuda|CUDA|cudatoolkit|cuda-toolkit|CUDAcore|nvidia/cuda)/([0-9]+(\.[0-9]+)*)$' x mm seen=" " i
  local fit=() minor=() major=()
  while IFS= read -r x; do
    x="${x%%(*}"; x="${x// /}"
    [[ "$x" =~ $re ]] || continue
    [[ "$seen" == *" $x "* ]] && continue
    seen+="$x "
    mm="$(ver_mm "${BASH_REMATCH[2]}")"
    version_ge "$mm" "$(ver_mm "$KOKKOS_NVCC_MIN")" || continue
    if [[ -n "$ARCH_NUM" ]] && (( ARCH_NUM < 75 && ${mm%%.*} >= 13 )); then continue; fi
    if [[ -z "$DRIVER_CUDA" ]] || version_ge "$DRIVER_CUDA" "$mm"; then
      fit+=("$x")
    elif [[ "${mm%%.*}" == "${DRIVER_CUDA%%.*}" ]]; then
      minor+=("$x")
    else
      major+=("$x")
    fi
  done < <(module -t avail 2>&1 | sed '/^$/d' | sort -Vr)
  for x in ${fit[@]+"${fit[@]}"}; do printf '%s\n' "$x"; done
  for (( i = ${#minor[@]} - 1; i >= 0; i-- )); do printf '%s\n' "${minor[i]}"; done
  for (( i = ${#major[@]} - 1; i >= 0; i-- )); do printf '%s\n' "${major[i]}"; done
  find_module 'nvhpc-nompi/[0-9].*' || true
  find_module '(nvhpc|NVHPC)/[0-9].*' || true
}

# Probe: compiles C++20 with the candidate nvcc + host g++ (catches an
# unsupported host gcc, see NVIDIA's host-compiler table), runs it on the GPU
# with that toolkit's runtime (catches toolkit/driver mismatch) and checks the
# target arch, all before hours of builds.
PROBE_DIR="$SRC_CACHE/.cuda-probe-$$"
rm -rf "$PROBE_DIR"; mkdir -p "$PROBE_DIR"
cat > "$PROBE_DIR/probe.cu" <<'EOF'
#include <cstdio>
#include <cuda_runtime.h>
int main() {
  int n = 0, drv = 0, rt = 0;
  const cudaError_t e = cudaGetDeviceCount(&n);
  cudaDriverGetVersion(&drv);
  cudaRuntimeGetVersion(&rt);
  if (e != cudaSuccess || n < 1) {
    std::printf("PROBE_ERROR %s (driver API %d, runtime %d)\n", cudaGetErrorString(e), drv, rt);
    return 1;
  }
  cudaDeviceProp p{};
  cudaGetDeviceProperties(&p, 0);
  std::printf("PROBE_OK cc=%d.%d devices=%d driver=%d runtime=%d name=%s\n",
              p.major, p.minor, n, drv, rt, p.name);
  return 0;
}
EOF
echo "int main() { return 0; }" > "$PROBE_DIR/arch.cu"
re_nvcc='release ([0-9]+\.[0-9]+)'
CUDA_FAIL=""
try_cuda() { # try_cuda <label>: is the nvcc now on PATH usable here?
  local label="$1" nvcc_path real dir d out
  CUDA_FAIL=""
  nvcc_path="$(command -v nvcc || true)"
  [[ -n "$nvcc_path" ]] || { CUDA_FAIL="$label: no nvcc on PATH"; return 1; }
  out="$("$nvcc_path" --version 2>&1 || true)"
  [[ "$out" =~ $re_nvcc ]] || { CUDA_FAIL="$label: cannot parse 'nvcc --version'"; return 1; }
  NVCC_VERSION="${BASH_REMATCH[1]}"; NVCC_MAJOR="${NVCC_VERSION%%.*}"
  # Kokkos 5.x needs C++20 device code; its configure hard-fails below nvcc
  # 12.2 (cmake/kokkos_compiler_id.cmake; the CI lane hit this with CUDA 12.0).
  version_ge "$NVCC_VERSION" "$(ver_mm "$KOKKOS_NVCC_MIN")" \
    || { CUDA_FAIL="$label: CUDA $NVCC_VERSION is too old (Kokkos $KOKKOS_VERSION requires nvcc >= $KOKKOS_NVCC_MIN for C++20)"; return 1; }
  if [[ -n "$ARCH_NUM" ]] && (( NVCC_MAJOR >= 13 && ARCH_NUM < 75 )); then
    CUDA_FAIL="$label: CUDA $NVCC_VERSION cannot target sm_$ARCH_NUM (CUDA 13 removed compute capability < 7.5)"; return 1
  fi
  # Toolkit root: the directory holding bin/nvcc and include/cuda_runtime.h.
  real="$(readlink -f "$nvcc_path")"; dir="$(dirname "$(dirname "$real")")"
  if [[ ! -f "$dir/include/cuda_runtime.h" ]]; then
    for d in "${CUDA_HOME:-}" "${CUDA_PATH:-}" "${CUDA_ROOT:-}" "${NVHPC_CUDA_HOME:-}" "$dir/../cuda"; do
      if [[ -n "$d" && -x "$d/bin/nvcc" && -f "$d/include/cuda_runtime.h" ]]; then dir="$(cd "$d" && pwd -P)"; break; fi
    done
  fi
  [[ -f "$dir/include/cuda_runtime.h" ]] || { CUDA_FAIL="$label: cannot locate the toolkit root of $real (no include/cuda_runtime.h)"; return 1; }
  rm -f "$PROBE_DIR/probe" "$PROBE_DIR/arch.o"
  # Commands whose failure is an expected outcome run directly in the "if"
  # with output to a log: inside $( ) the inherited ERR trap (set -E) would
  # fire and print a false "failed at line" message.
  if ! "$dir/bin/nvcc" -ccbin "$REAL_CXX" -std=c++20 -o "$PROBE_DIR/probe" "$PROBE_DIR/probe.cu" >"$PROBE_DIR/log" 2>&1; then
    CUDA_FAIL="$label: CUDA $NVCC_VERSION cannot compile C++20 with host $COMPILER_VERSION: $(head -n 3 "$PROBE_DIR/log" | tr '\n' ' ')"; return 1
  fi
  if ! "$PROBE_DIR/probe" >"$PROBE_DIR/log" 2>&1; then
    CUDA_FAIL="$label: a CUDA $NVCC_VERSION program cannot use this GPU (driver supports CUDA ${DRIVER_CUDA:-?}): $(head -n 3 "$PROBE_DIR/log" | tr '\n' ' ')"; return 1
  fi
  PROBE_OUT="$(head -n 1 "$PROBE_DIR/log")"
  if [[ -n "$ARCH_NUM" ]] && ! "$dir/bin/nvcc" -ccbin "$REAL_CXX" -arch="sm_$ARCH_NUM" -c -o "$PROBE_DIR/arch.o" "$PROBE_DIR/arch.cu" >"$PROBE_DIR/log" 2>&1; then
    CUDA_FAIL="$label: CUDA $NVCC_VERSION does not support sm_$ARCH_NUM ($KOKKOS_ARCH)"; return 1
  fi
  CUDA_DIR="$dir"
  return 0
}

CUDA_OK=0; CUDA_REJECTS=(); CUDA_CANDIDATES=()
if [[ -n "$FREHG_MODULE_CUDA" ]]; then
  CUDA_CANDIDATES=("$FREHG_MODULE_CUDA")
else
  while IFS= read -r x; do [[ -n "$x" ]] && CUDA_CANDIDATES+=("$x"); done < <(cuda_module_candidates)
  echo "CUDA module candidates (best first): ${CUDA_CANDIDATES[*]:-none}"
fi
for cand in ${CUDA_CANDIDATES[@]+"${CUDA_CANDIDATES[@]}"}; do
  echo "Trying CUDA module: $cand"
  if ! module load "$cand"; then
    CUDA_REJECTS+=("$cand: module load failed"); continue
  fi
  if try_cuda "$cand"; then
    LOADED_MODULES+=("$cand"); CUDA_OK=1
    [[ "$cand" == *nvhpc* ]] && warn "Using $cand (NVIDIA HPC SDK) for CUDA. If it bundles an MPI, make sure the 'MPI compiler wrapper' in step 4 is the site Open MPI."
    break
  fi
  echo "  rejected: $CUDA_FAIL"; CUDA_REJECTS+=("$CUDA_FAIL")
  module unload "$cand" || true
done
if [[ "$CUDA_OK" -eq 0 && -z "$FREHG_MODULE_CUDA" ]]; then
  if command -v nvcc >/dev/null 2>&1; then
    echo "Trying the nvcc already on PATH: $(command -v nvcc)"
    if try_cuda "nvcc on PATH"; then CUDA_OK=1; else CUDA_REJECTS+=("$CUDA_FAIL"); fi
  fi
  if [[ "$CUDA_OK" -eq 0 && -x /usr/local/cuda/bin/nvcc ]]; then
    echo "Trying /usr/local/cuda"
    export PATH="/usr/local/cuda/bin:$PATH"
    if try_cuda "/usr/local/cuda"; then CUDA_OK=1; else CUDA_REJECTS+=("$CUDA_FAIL"); fi
  fi
fi
if [[ "$CUDA_OK" -eq 0 ]]; then
  [[ ${#CUDA_REJECTS[@]} -gt 0 ]] || CUDA_REJECTS=("no CUDA module, no nvcc on PATH, no /usr/local/cuda")
  printf '  - %s\n' "${CUDA_REJECTS[@]}" >&2
  die "No usable CUDA toolkit (reasons above). Run: module spider cuda ; then rerun with FREHG_MODULE_CUDA=<exact-module-name>. It needs nvcc >= $KOKKOS_NVCC_MIN, must support host $COMPILER_VERSION (else change FREHG_MODULE_COMPILER) and must run under this driver (CUDA ${DRIVER_CUDA:-?}). Host-compiler table: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/"
fi
echo "nvcc: $CUDA_DIR/bin/nvcc (CUDA $NVCC_VERSION)"
echo "CUDA probe: $PROBE_OUT"
if [[ -n "$DRIVER_CUDA" ]] && ! version_ge "$DRIVER_CUDA" "$NVCC_VERSION"; then
  warn "CUDA toolkit $NVCC_VERSION is newer than the driver's CUDA $DRIVER_CUDA. It ran the probe (CUDA compatibility), but newer-than-driver features can still fail at run time; the step-10 GPU run is the real test."
fi
re_probe_cc='cc=([0-9]+\.[0-9]+)'
if [[ -z "$KOKKOS_ARCH" ]]; then
  [[ "$PROBE_OUT" =~ $re_probe_cc ]] || die "Could not determine the GPU's compute capability. Set FREHG_CUDA_ARCH (e.g. AMPERE80, ADA89, HOPPER90)."
  GPU_CC="${BASH_REMATCH[1]}"; GPU_CC_NUM="$(cc_num "$GPU_CC")"
  KOKKOS_ARCH="$(cc_to_kokkos_arch "$GPU_CC" || true)"
  [[ -n "$KOKKOS_ARCH" ]] || die "Compute capability $GPU_CC has no Kokkos 5.1.1 arch name. Set FREHG_CUDA_ARCH to one of: $KNOWN_ARCHS."
  ARCH_NUM="${KOKKOS_ARCH//[!0-9]/}"
  echo "Kokkos arch (from the probe's compute capability $GPU_CC): $KOKKOS_ARCH"
  if (( NVCC_MAJOR >= 13 && ARCH_NUM < 75 )); then
    die "CUDA $NVCC_VERSION cannot target sm_$ARCH_NUM: CUDA 13 removed Maxwell/Pascal/Volta (compute capability < 7.5). Load a CUDA 12.x module (FREHG_MODULE_CUDA=cuda/12.x)."
  fi
  "$CUDA_DIR/bin/nvcc" -ccbin "$REAL_CXX" -arch="sm_$ARCH_NUM" -c -o "$PROBE_DIR/arch.o" "$PROBE_DIR/arch.cu" \
    || die "nvcc $NVCC_VERSION does not support sm_$ARCH_NUM ($KOKKOS_ARCH). Load a newer CUDA module."
fi
if [[ -n "$GPU_CC_NUM" && "$ARCH_NUM" != "$GPU_CC_NUM" ]]; then
  warn "Building for sm_$ARCH_NUM but this GPU is sm_$GPU_CC_NUM; the verification run in step 10 may fail on this node."
fi
rm -rf "$PROBE_DIR"
# Pin every consumer (nvcc_wrapper reads CUDA_ROOT; PETSc gets --with-cuda-dir;
# CMake gets CUDAToolkit_ROOT) to this one toolkit.
export CUDA_ROOT="$CUDA_DIR" CUDA_HOME="$CUDA_DIR"
export PATH="$CUDA_DIR/bin:$PATH"
echo "CUDA toolkit root: $CUDA_DIR"
# PETSc 3.25.1's CUDA package hard-requires NVML (nvml.h, nvmlInit_v2).
NVML_H="$(find -L "$CUDA_DIR" -maxdepth 4 -name nvml.h -print -quit 2>/dev/null || true)"
[[ -n "$NVML_H" ]] || warn "nvml.h was not found under $CUDA_DIR. PETSc 3.25.1 requires the NVML headers (part of the full CUDA toolkit); its configure will fail without them. Ask support for a complete CUDA toolkit module."
export NVCC_WRAPPER_DEFAULT_COMPILER="$REAL_CXX"
echo "PETSc/hypre CUDA arch: $ARCH_NUM ; NVCC_WRAPPER_DEFAULT_COMPILER=$NVCC_WRAPPER_DEFAULT_COMPILER"

say "4 of 10: Load MPI and detect whether it is GPU-aware (CUDA-aware)"
# MPI after CUDA: on Lmod hierarchies a CUDA-aware MPI often lives under the
# CUDA module, and loading MPI last keeps its mpicc first on PATH (an nvhpc
# module can carry its own HPC-X MPI).
load_module MPI "$FREHG_MODULE_MPI" "$SITE_MODULE_MPI" '(openmpi|OpenMPI|mpi/openmpi|mpich|MPICH|intel-mpi|impi)/[0-9].*'
need_command mpicc
# A CUDA-aware MPI module often loads the CUDA toolkit it was built with as a
# prerequisite (openmpi/4.1.6-gpu loads cuda/12.4 on the A800 cluster). That
# would swap nvcc on PATH under the toolkit selected above and mix two CUDA
# versions in one build. Require the MPI module to leave the selection alone.
NVCC_AFTER_MPI="$(readlink -f "$(command -v nvcc 2>/dev/null)" 2>/dev/null || true)"
if [[ -n "$NVCC_AFTER_MPI" && "$NVCC_AFTER_MPI" != "$(readlink -f "$CUDA_DIR/bin/nvcc")" ]]; then
  die "Loading the MPI module put another CUDA toolkit first on PATH (nvcc is now $NVCC_AFTER_MPI; the build selected $CUDA_DIR). The MPI module most likely loads the CUDA it was built with ('module show <mpi-module>' lists it). Build with that same CUDA: rerun with FREHG_MODULE_CUDA=<its CUDA module>, e.g. FREHG_MODULE_CUDA=cuda/12.4 for openmpi/4.1.6-gpu."
fi
need_command mpicxx
echo "MPI compiler wrapper: $(command -v mpicxx)"
MPI_VERSION_OUT="$(mpiexec --version 2>&1 || true)"
printf '%s\n' "${MPI_VERSION_OUT%%$'\n'*}"
# Open MPI: ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
# (https://www.open-mpi.org/faq/?category=runcuda). Captured first: a grep -q
# that exits early would SIGPIPE ompi_info and, under pipefail, read as false.
MPI_GPU_AWARE=0
MPI_GPU_AWARE_HOW="not determined (not Open MPI or MPICH)"
OMPI_INFO="$(dirname "$(command -v mpicc)")/ompi_info"
[[ -x "$OMPI_INFO" ]] || OMPI_INFO="$(command -v ompi_info || true)"
case "${FREHG_MPI_GPU_AWARE:-auto}" in
  1|0) MPI_GPU_AWARE="$FREHG_MPI_GPU_AWARE"; MPI_GPU_AWARE_HOW="FREHG_MPI_GPU_AWARE=$FREHG_MPI_GPU_AWARE (override)" ;;
  auto)
    if [[ -n "$OMPI_INFO" ]]; then
      OMPI_ALL="$("$OMPI_INFO" --parsable --all 2>/dev/null || true)"
      CUDA_LINE="$(grep -m1 'mpi_built_with_cuda_support:value' <<<"$OMPI_ALL" || true)"
      echo "ompi_info: ${CUDA_LINE:-<no mpi_built_with_cuda_support line>}"
      OMPI_PLAIN="$("$OMPI_INFO" 2>/dev/null || true)"
      EXT_LINE="$(grep -m1 'MPI extensions' <<<"$OMPI_PLAIN" || true)"
      [[ -n "$EXT_LINE" ]] && echo "ompi_info:${EXT_LINE}   (PETSc's own detection reads this line)"
      if [[ "$CUDA_LINE" == *":value:true" ]]; then
        MPI_GPU_AWARE=1; MPI_GPU_AWARE_HOW="Open MPI built with CUDA support"
      else
        MPI_GPU_AWARE_HOW="Open MPI built WITHOUT CUDA support"
      fi
    elif command -v mpichversion >/dev/null 2>&1; then
      MPICH_OUT="$(mpichversion 2>&1 || true)"
      if grep -Eqi -- '--with-cuda|--enable-gpu|yaksa.*cuda' <<<"$MPICH_OUT"; then
        MPI_GPU_AWARE=1; MPI_GPU_AWARE_HOW="MPICH configured with CUDA"
      else
        MPI_GPU_AWARE_HOW="MPICH without visible CUDA configure options"
      fi
    fi ;;
  *) die "FREHG_MPI_GPU_AWARE must be auto, 0 or 1." ;;
esac
echo "GPU-aware MPI: $([[ $MPI_GPU_AWARE -eq 1 ]] && echo yes || echo no) -- $MPI_GPU_AWARE_HOW"
if [[ "$MPI_GPU_AWARE" -eq 1 ]]; then
  cat <<EOF
  -> frehg2: runtime.gpu_aware_mpi may be 'on' for multi-GPU runs (device
     buffers go straight to MPI; experimental, exercised only by p6 item 4).
     'auto' (the default) follows the compile-time FREHG_GPU_AWARE_MPI, which
     this build sets to ${FREHG_GPU_AWARE_MPI:-OFF}. Single-GPU runs are unaffected.
  -> PETSc keeps its default -use_gpu_aware_mpi 1 and checks it at run time.
     Some Open MPI builds still need OMPI_MCA_pml=ucx for CUDA buffers.
EOF
else
  CUDA_MPI_CANDIDATES="$(module -t avail 2>&1 | grep -Ei '((openmpi|ompi).*(cuda|gpu))|hpcx|hpc-x|mvapich2-gdr' | tr '\n' ' ' || true)"
  cat <<EOF
  -> Not fatal. Single-GPU (1-rank) runs are unaffected. For multi-GPU runs:
     frehg2: keep runtime.gpu_aware_mpi at 'auto' or 'off' (halos staged
       through host mirrors). Never 'on' with this MPI.
     PETSc: its default -use_gpu_aware_mpi 1 makes PETSc ABORT at GPU init
       ("your MPI is not GPU-aware"), even on 1 rank. Every run must set
         export PETSC_OPTIONS="-use_gpu_aware_mpi 0"
       or pass -use_gpu_aware_mpi 0 after the YAML on the frehg command line.
       solver.petsc_options_file is too late: PETSc reads this flag inside
       PetscInitialize, and frehg2 inserts that file afterwards.
     hypre: PETSc's CUDA hypre build defaults to --enable-gpu-aware-mpi.
       This script appends --disable-gpu-aware-mpi so hypre stages through
       the host.
     For a CUDA-aware MPI, ask support. Candidate modules here: ${CUDA_MPI_CANDIDATES:-none found}
     Switching MPI means rebuilding this whole prefix against it.
EOF
fi

say "5 of 10: Load and verify a parallel HDF5 module if one exists"
if [[ -n "$FREHG_MODULE_HDF5" ]]; then
  echo "Loading requested HDF5 module: $FREHG_MODULE_HDF5"
  module load "$FREHG_MODULE_HDF5" || die "Could not load HDF5 module '$FREHG_MODULE_HDF5'."
  LOADED_MODULES+=("$FREHG_MODULE_HDF5")
else
  HDF_CANDIDATE=""
  if module_listed "$SITE_MODULE_HDF5"; then
    HDF_CANDIDATE="$SITE_MODULE_HDF5"
  else
    HDF_CANDIDATE="$(find_module '(hdf5|HDF5|parallel-hdf5|hdf5-parallel)/[0-9].*' || true)"
  fi
  if [[ -n "$HDF_CANDIDATE" ]]; then
    echo "Trying HDF5 module: $HDF_CANDIDATE"
    if module load "$HDF_CANDIDATE"; then LOADED_MODULES+=("$HDF_CANDIDATE"); else HDF_CANDIDATE=""; fi
  else
    echo "No HDF5 module found; parallel HDF5 will be built below."
  fi
fi
USE_MODULE_HDF5=0
# A CMake-built HDF5 ships only h5cc (which wraps mpicc when parallel); the
# autotools build additionally ships h5pcc. Accept either wrapper.
HDF5_WRAPPER="$(command -v h5pcc || command -v h5cc || true)"
HDF5_SHOWCONFIG=""
[[ -n "$HDF5_WRAPPER" ]] && HDF5_SHOWCONFIG="$("$HDF5_WRAPPER" -showconfig 2>/dev/null || true)"
if grep -Eqi 'Parallel HDF5:[[:space:]]+yes' <<<"$HDF5_SHOWCONFIG"; then
  USE_MODULE_HDF5=1
  echo "Usable parallel HDF5 detected: $HDF5_WRAPPER"
else
  echo "A usable parallel HDF5 was not detected; building HDF5 $HDF5_VERSION into $PREFIX."
  # Drop an auto-tried serial HDF5 module again so its paths cannot shadow the
  # parallel build and the printed run-time module list stays accurate.
  if [[ -z "$FREHG_MODULE_HDF5" && -n "${HDF_CANDIDATE:-}" ]]; then
    echo "Unloading the non-parallel HDF5 module $HDF_CANDIDATE."
    module unload "$HDF_CANDIDATE" || true
    KEPT=()
    for m in "${LOADED_MODULES[@]}"; do [[ "$m" == "$HDF_CANDIDATE" ]] || KEPT+=("$m"); done
    LOADED_MODULES=("${KEPT[@]}")
  fi
fi

say "6 of 10: Build yaml-cpp and the CUDA Kokkos in $PREFIX"
# Command to reset this (CUDA-only) prefix but keep downloaded archives.
WIPE_HINT="rm -rf '$PREFIX'/{bin,include,lib,lib64,share} '$SRC_CACHE'/build-* '$SRC_CACHE/petsc-$PETSC_VERSION'"
# One prefix holds parallel HDF5, PETSc and hypre linked against ONE MPI. A
# rerun with another MPI module would otherwise reuse them and mix two MPI
# builds in one process. The first build records its mpicc; prefixes built
# before this check are identified by the mpicc PETSc recorded.
MPICC_NOW="$(readlink -f "$(command -v mpicc)" 2>/dev/null || command -v mpicc)"
MPI_STAMP="$PREFIX/.frehg-mpicc"
MPICC_THEN=""
if [[ -s "$MPI_STAMP" ]]; then
  MPICC_THEN="$(<"$MPI_STAMP")"
else
  PV_OLD="$(first_existing "$PREFIX/lib/petsc/conf/petscvariables" "$PREFIX/lib64/petsc/conf/petscvariables" || true)"
  [[ -n "$PV_OLD" ]] && MPICC_THEN="$(sed -n 's/^PCC[[:space:]]*=[[:space:]]*\([^[:space:]]*\).*/\1/p' "$PV_OLD" | head -n 1)"
  [[ -n "$MPICC_THEN" ]] && MPICC_THEN="$(readlink -f "$MPICC_THEN" 2>/dev/null || echo "$MPICC_THEN")"
fi
if [[ -n "$MPICC_THEN" && "$MPICC_THEN" != "$MPICC_NOW" ]]; then
  die "$PREFIX was built with the MPI at $MPICC_THEN, but this run loaded $MPICC_NOW. HDF5, PETSc and hypre in a prefix are linked against one MPI. Build the new MPI into its own prefix (this keeps the existing build), e.g.
    FREHG_PREFIX=\$HOME/frehg-deps-cuda-gpumpi FREHG_BUILD_NAME=build-cuda-gpumpi FREHG_MODULE_MPI=<module> bash build_frehg2_cuda.sh
  or reset this prefix to rebuild it against the new MPI: $WIPE_HINT"
fi
echo "$MPICC_NOW" > "$MPI_STAMP"
echo "MPI for this prefix: $MPICC_NOW"
YAML_SRC="$SRC_CACHE/yaml-cpp-${YAML_VERSION}"
fetch_tarball "$URL_YAML" "$ARC_YAML" "$YAML_SRC"
if ! first_existing "$PREFIX/lib/cmake/yaml-cpp/yaml-cpp-config.cmake" "$PREFIX/lib64/cmake/yaml-cpp/yaml-cpp-config.cmake" >/dev/null; then
  # yaml-cpp has no device code: build it with the plain host compiler.
  cmake_build_install "$YAML_SRC" "$SRC_CACHE/build-yaml-cpp-${YAML_VERSION}" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_C_COMPILER="$REAL_CC" -DCMAKE_CXX_COMPILER="$REAL_CXX" \
    -DYAML_BUILD_SHARED_LIBS=ON -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
else
  echo "yaml-cpp already installed in $PREFIX; reusing it."
fi

KOKKOS_SRC="$SRC_CACHE/kokkos-${KOKKOS_VERSION}"
fetch_tarball "$URL_KOKKOS" "$ARC_KOKKOS" "$KOKKOS_SRC"
KOKKOS_CONFIG="$(first_existing "$PREFIX/lib/cmake/Kokkos/KokkosConfig.cmake" "$PREFIX/lib64/cmake/Kokkos/KokkosConfig.cmake" || true)"
KOKKOS_CFG_H="$PREFIX/include/KokkosCore_config.h"
if [[ -z "$KOKKOS_CONFIG" ]]; then
  # BUILD_SHARED_LIBS=ON is load-bearing: frehg links Kokkos and so does
  # PETSc's downloaded kokkos-kernels. A static Kokkos core is absorbed
  # into BOTH the frehg executable and libkokkoskernels, giving two copies
  # of Kokkos's runtime singleton in one process -- it initializes twice
  # and the first VecKokkos access segfaults. One shared libkokkoscore keeps
  # a single runtime for frehg2 AND PETSc (docs/agents/build.md invariant 7).
  # Shared libraries also force relocatable device code OFF, which Kokkos
  # requires anyway ("Relocatable device code requires static libraries").
  # Device flags (Kokkos 5.1.1 cmake/kokkos_enable_options.cmake):
  #   CUDA_LAMBDA=ON  deprecated and ON by default with CUDA, but PETSc's
  #                   kokkos.py refuses a Kokkos without KOKKOS_ENABLE_CUDA_LAMBDA
  #                   in KokkosCore_config.h, so it is pinned explicitly.
  #   CUDA_CONSTEXPR  left at its default (OFF for nvcc): not needed by Kokkos,
  #                   Kokkos Kernels 5.1.0, PETSc or frehg2 (the cuda-compile CI
  #                   lane builds frehg2 + tests without it).
  #   IMPL_CUDA_MALLOC_ASYNC=OFF  what PETSc sets when it builds Kokkos itself:
  #                   cudaMallocAsync interferes with CUDA-aware MPI.
  # nvcc_wrapper from the Kokkos source is the CXX compiler (gpu-acceptance.md
  # "CMAKE_CXX_COMPILER=<kokkos>/bin/nvcc_wrapper"; ci_install_deps.sh).
  rm -rf "$SRC_CACHE/build-kokkos-${KOKKOS_VERSION}"
  cmake_build_install "$KOKKOS_SRC" "$SRC_CACHE/build-kokkos-${KOKKOS_VERSION}" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
    -DCMAKE_C_COMPILER="$REAL_CC" -DCMAKE_CXX_COMPILER="$KOKKOS_SRC/bin/nvcc_wrapper" \
    -DCUDAToolkit_ROOT="$CUDA_DIR" \
    -DBUILD_SHARED_LIBS=ON -DCMAKE_CXX_STANDARD=20 \
    -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON \
    "-DKokkos_ARCH_${KOKKOS_ARCH}=ON" \
    -DKokkos_ENABLE_CUDA_LAMBDA=ON \
    -DKokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF \
    -DKokkos_ENABLE_TESTS=OFF
  KOKKOS_CONFIG="$(first_existing "$PREFIX/lib/cmake/Kokkos/KokkosConfig.cmake" "$PREFIX/lib64/cmake/Kokkos/KokkosConfig.cmake" || true)"
  [[ -n "$KOKKOS_CONFIG" ]] || die "Kokkos installed but KokkosConfig.cmake is missing from $PREFIX."
else
  echo "Kokkos already installed in $PREFIX; checking it matches this GPU build."
fi
# Verify (fresh or reused) Kokkos: version, CUDA backend, arch, lambda, shared.
KOKKOS_FOUND_VERSION="$(sed -n 's/.*set(PACKAGE_VERSION "\([0-9.]*\)").*/\1/p' "$(dirname "$KOKKOS_CONFIG")/KokkosConfigVersion.cmake" 2>/dev/null || true)"
[[ -f "$KOKKOS_CFG_H" ]] || die "Kokkos header $KOKKOS_CFG_H is missing. Reset the CUDA prefix: $WIPE_HINT"
grep -Eq '^#define KOKKOS_ENABLE_CUDA([[:space:]]|$)' "$KOKKOS_CFG_H" \
  || die "The Kokkos in $PREFIX has no CUDA backend (a CPU Kokkos?). Refusing to reuse it. Use a fresh FREHG_PREFIX or reset it: $WIPE_HINT"
FOUND_ARCH="$(grep -Eo '^#define KOKKOS_ARCH_(MAXWELL|PASCAL|VOLTA|TURING|AMPERE|ADA|HOPPER|BLACKWELL)[0-9]+' "$KOKKOS_CFG_H" | sed 's/^#define KOKKOS_ARCH_//' | tr '\n' ' ' || true)"
[[ " $FOUND_ARCH " == *" $KOKKOS_ARCH "* ]] \
  || die "The Kokkos in $PREFIX was built for GPU arch '${FOUND_ARCH:-none}', not $KOKKOS_ARCH. For a second GPU type use FREHG_PREFIX=\$HOME/frehg-deps-cuda-$KOKKOS_ARCH; otherwise reset: $WIPE_HINT"
grep -Eq '^#define KOKKOS_ENABLE_CUDA_LAMBDA' "$KOKKOS_CFG_H" \
  || die "Kokkos lacks KOKKOS_ENABLE_CUDA_LAMBDA, which PETSc requires. Reset the prefix: $WIPE_HINT"
[[ -z "$KOKKOS_FOUND_VERSION" || "$KOKKOS_FOUND_VERSION" == "$KOKKOS_VERSION" ]] \
  || die "Kokkos $KOKKOS_FOUND_VERSION is installed in $PREFIX but $KOKKOS_VERSION was requested. Reset the prefix: $WIPE_HINT"
first_existing "$PREFIX"/lib/libkokkoscore.so "$PREFIX"/lib64/libkokkoscore.so >/dev/null \
  || die "No shared libkokkoscore.so in $PREFIX (static Kokkos => two runtimes => segfault). Reset the prefix: $WIPE_HINT"
[[ -x "$PREFIX/bin/nvcc_wrapper" ]] || die "Kokkos did not install $PREFIX/bin/nvcc_wrapper."
echo "Kokkos ${KOKKOS_FOUND_VERSION:-$KOKKOS_VERSION}: CUDA + $KOKKOS_ARCH + lambda, shared, nvcc_wrapper at $PREFIX/bin/nvcc_wrapper"

say "7 of 10: Obtain a matching parallel HDF5"
if [[ "$USE_MODULE_HDF5" -eq 0 ]]; then
  HDF5_SRC="$SRC_CACHE/hdf5-hdf5_${HDF5_VERSION}"
  fetch_tarball "$URL_HDF5" "$ARC_HDF5" "$HDF5_SRC"
  if ! grep -q '#define H5_HAVE_PARALLEL 1' "$PREFIX/include/H5pubconf.h" 2>/dev/null; then
    cmake_build_install "$HDF5_SRC" "$SRC_CACHE/build-hdf5-${HDF5_VERSION}" \
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" \
      -DCMAKE_C_COMPILER="$(command -v mpicc)" \
      -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_TOOLS=ON \
      -DHDF5_BUILD_FORTRAN=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TESTING=OFF \
      -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF -DHDF5_ENABLE_SZIP_SUPPORT=OFF
  else
    echo "User-built parallel HDF5 already installed in $PREFIX; reusing it."
  fi
  export PATH="$PREFIX/bin:$PATH"
  export HDF5_ROOT="$PREFIX"
  grep -q '#define H5_HAVE_PARALLEL 1' "$PREFIX/include/H5pubconf.h" 2>/dev/null \
    || die "HDF5 is not parallel. Frehg2 requires an MPI-enabled HDF5 built with the currently loaded MPI."
  echo "Parallel HDF5 confirmed (source build): $PREFIX"
else
  echo "Using module-provided parallel HDF5: $HDF5_WRAPPER"
fi

say "8 of 10: Build PETSc (CUDA + shared Kokkos + Kokkos Kernels + CUDA hypre)"
PETSC_SRC="$SRC_CACHE/petsc-${PETSC_VERSION}"
fetch_tarball "$URL_PETSC" "$ARC_PETSC" "$PETSC_SRC"
# Memory-aware cap for PETSc's own downloaded-package builds (Kokkos Kernels
# ETI units peak at several GB each, more under nvcc): ~8 GB per job, never
# above $JOBS (same rule as scripts/ci_install_deps.sh). SLURM's allocation
# is used when known, because /proc/meminfo shows the whole node.
MEM_GB=0
if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" =~ ^[0-9]+$ ]]; then
  MEM_GB=$(( SLURM_MEM_PER_NODE / 1024 ))
elif [[ -n "${SLURM_MEM_PER_CPU:-}" && "${SLURM_MEM_PER_CPU}" =~ ^[0-9]+$ ]]; then
  MEM_GB=$(( SLURM_MEM_PER_CPU * ${SLURM_CPUS_PER_TASK:-1} / 1024 ))
elif [[ -r /proc/meminfo ]]; then
  MEM_GB="$(awk '/MemAvailable/ {print int($2/1024/1024)}' /proc/meminfo)"
fi
if [[ "$MEM_GB" -gt 0 ]]; then
  MEM_NP=$(( MEM_GB / 8 )); (( MEM_NP < 1 )) && MEM_NP=1
  PETSC_MAKE_NP=$(( MEM_NP < JOBS ? MEM_NP : JOBS ))
else
  PETSC_MAKE_NP="$JOBS"
fi
PETSC_MAKE_NP="${FREHG_PETSC_MAKE_NP:-$PETSC_MAKE_NP}"
echo "PETSc --with-make-np=$PETSC_MAKE_NP (memory ${MEM_GB} GB, jobs $JOBS)"
PETSC_CONF_H="$PREFIX/include/petscconf.h"
if ! first_existing "$PREFIX/lib/petsc/conf/petscvariables" "$PREFIX/lib64/petsc/conf/petscvariables" >/dev/null; then
  # The three PETSc-pinned packages are given as local tarballs (identical to
  # what PETSc's --download-* would fetch) so no network is needed here once
  # --fetch-only has run.
  fetch_file "$URL_KK" "$PKG_KK"
  fetch_file "$URL_HYPRE" "$PKG_HYPRE"
  fetch_file "$URL_F2C" "$PKG_F2C"
  # PETSc correctly uses MPI wrappers for ITS build (not Frehg2's compiler).
  #  --with-cuda --with-cuda-dir --with-cuda-arch: the device lane (gpu-
  #     acceptance.md "--with-cuda=1"; PETSc's arch is the numeric sm, e.g. 89,
  #     and is handed on to hypre's --with-gpu-arch).
  #  --with-kokkos-dir=$PREFIX: PETSc reuses the ONE shared CUDA Kokkos built
  #     above (no --download-kokkos), so PETSc and frehg2 share one Kokkos
  #     runtime; --download-kokkos-kernels builds Kokkos Kernels against it.
  #     With CUDA found, PETSc's aijkokkos/VecKokkos run on the device -- the
  #     path frehg2 forces on device builds (solver.*.mat_type -> aijkokkos).
  #  --download-hypre: with CUDA found, PETSc's hypre.py configures hypre
  #     --with-cuda --with-gpu-arch=<arch> CUDA_HOME=... --enable-gpu-aware-mpi
  #     (hypre's own default is off: hypre.readthedocs.io ch-misc "GPU build
  #     options"). Without a CUDA-aware MPI, --disable-gpu-aware-mpi is
  #     appended via --download-hypre-configure-arguments (appended last, so it
  #     wins in hypre's autoconf). Umpire is not downloaded: hypre warns "not
  #     recommended" (performance only); FREHG_PETSC_EXTRA=--download-umpire
  #     adds it, but its build needs git submodules (git on the node).
  #  --with-openmp=1: Kokkos has an OpenMP host space; same as the CPU script
  #     and the cuda-compile CI lane (CUDA hypre ignores OpenMP).
  #  --with-make-np caps PETSc's downloaded-package builds. Without it PETSc
  #     uses every core and kokkos-kernels OOM-kills cc1plus ("make -j136").
  PETSC_ARGS=(
    --prefix="$PREFIX"
    --with-cc="$(command -v mpicc)" --with-cxx="$(command -v mpicxx)" --with-fc=0
    --with-debugging=0
    --download-f2cblaslapack="$PKG_F2C"
    --download-hypre="$PKG_HYPRE"
    --with-kokkos-dir="$PREFIX"
    --download-kokkos-kernels="$PKG_KK"
    --with-openmp=1
    --with-cuda=1 --with-cuda-dir="$CUDA_DIR" --with-cuda-arch="$ARCH_NUM"
    --with-make-np="$PETSC_MAKE_NP"
  )
  # hypre's extra configure arguments (appended last by PETSc, so they win):
  #  - --disable-gpu-aware-mpi without a CUDA-aware MPI (see above);
  #  - LIBS=-ldl: hypre's CUDA objects pull in NVIDIA's header-only NVTX v3
  #    through CUB/Thrust, whose loader calls dlopen/dlclose, and hypre links
  #    libHYPRE.so with -Wl,-z,defs (no undefined symbols). Before glibc 2.34
  #    (RHEL/CentOS 7-8) those live in libdl, which hypre's link line lacks:
  #    "undefined reference to `dlclose'" (first seen on the A800 cluster).
  #    Harmless on newer glibc, where libdl is an empty stub that still links.
  #    PETSc passes hypre no LIBS of its own outside cross builds.
  HYPRE_EXTRA=( LIBS=-ldl )
  [[ "$MPI_GPU_AWARE" -eq 1 ]] || HYPRE_EXTRA=( --disable-gpu-aware-mpi "${HYPRE_EXTRA[@]}" )
  PETSC_ARGS+=( --download-hypre-configure-arguments="${HYPRE_EXTRA[*]}" )
  read -r -a PETSC_USER_EXTRA <<<"${FREHG_PETSC_EXTRA:-}"
  PETSC_ARGS+=( ${PETSC_USER_EXTRA[@]+"${PETSC_USER_EXTRA[@]}"} )
  dump_petsc_log() { # surface the real error buried in PETSc's logs
    local log="$1" n
    [[ -f "$log" ]] || return 0
    # configure.log is full of EXPECTED failures (PETSc's probe programs:
    # MPI vendor macros, Intel/ARM flags, absent Kokkos backends), so a plain
    # grep for "error:" shows noise. A failed downloaded package is reported
    # once as "Error running <step> on <PACKAGE>: <its stdout+stderr>"; show
    # the matching lines of that block first.
    n="$(grep -nE 'Error running (configure|make|make; make install|cmake) on ' "$log" | head -n 1 | cut -d: -f1)"
    if [[ -n "$n" ]]; then
      echo "======== $log: the failing package step (from line $n) ========" >&2
      sed -n "${n}p" "$log" | cut -c1-300 >&2
      tail -n +"$n" "$log" \
        | grep -nE -B2 -A2 'error|fatal|cannot find|undefined reference|No such file|not found|Killed|out of memory|Error [0-9]' \
        | cut -c1-400 | head -n 120 >&2 || true
      echo "======== the same block, last 40 lines before PETSc's summary ========" >&2
      tail -n +"$n" "$log" | grep -n -m 1 'UNABLE to CONFIGURE' >/dev/null \
        && tail -n +"$n" "$log" | sed '/UNABLE to CONFIGURE/q' | tail -n 40 | cut -c1-400 >&2
    else
      echo "======== $log: compiler errors / OOM markers ========" >&2
      grep -nE 'error:|fatal error|Killed|cannot allocate|out of memory|virtual memory exhausted' "$log" | tail -n 60 >&2 || true
    fi
    echo "======== $log: last 80 lines ========" >&2
    tail -n 80 "$log" >&2 || true
  }
  pushd "$PETSC_SRC" >/dev/null
  unset PETSC_DIR PETSC_ARCH   # never let a site PETSc module steer this build
  export PETSC_DIR="$PETSC_SRC" PETSC_ARCH=arch-frehg2-cuda
  echo "PETSc configure: ./configure ${PETSC_ARGS[*]}"
  ./configure "${PETSC_ARGS[@]}" \
    || { dump_petsc_log "$PETSC_SRC/configure.log"; die "PETSc configure failed (log: $PETSC_SRC/configure.log). If Kokkos Kernels failed to compile against Kokkos $KOKKOS_VERSION, reset ($WIPE_HINT) and rerun with FREHG_KOKKOS_VERSION=$PETSC_KK_VERSION."; }
  make -j "$JOBS" all || { dump_petsc_log "$PETSC_SRC/$PETSC_ARCH/lib/petsc/conf/make.log"; die "PETSc make failed."; }
  make install
  popd >/dev/null
  unset PETSC_DIR PETSC_ARCH
else
  echo "PETSc already installed in $PREFIX; reusing it after the checks below."
fi
[[ -f "$PETSC_CONF_H" ]] || die "$PETSC_CONF_H is missing. Reset the prefix: $WIPE_HINT"
petsc_has() { grep -Eq "^#define PETSC_$1 1([[:space:]]|\$)" "$PETSC_CONF_H"; }
for feature in HAVE_CUDA HAVE_KOKKOS HAVE_KOKKOS_KERNELS HAVE_HYPRE HAVE_HYPRE_DEVICE; do
  petsc_has "$feature" || die "PETSc in $PREFIX lacks PETSC_$feature (needed: CUDA, Kokkos, Kokkos Kernels for aijkokkos, hypre with device support for solver amg). Reset the prefix and rerun: $WIPE_HINT"
  echo "petscconf.h: PETSC_$feature = yes"
done
PETSC_MIN_ARCH="$(sed -n 's/^#define PETSC_PKG_CUDA_MIN_ARCH[[:space:]]*"\{0,1\}\([0-9]*\).*/\1/p' "$PETSC_CONF_H")"
if [[ -n "$PETSC_MIN_ARCH" && "$PETSC_MIN_ARCH" != "$ARCH_NUM" ]]; then
  die "PETSc in $PREFIX was configured for sm_$PETSC_MIN_ARCH, not sm_$ARCH_NUM. Reset the prefix: $WIPE_HINT"
fi
echo "petscconf.h: PETSC_PKG_CUDA_MIN_ARCH = ${PETSC_MIN_ARCH:-<not recorded>}"
echo "petscconf.h: PETSC_HAVE_MPI_GPU_AWARE = $(petsc_has HAVE_MPI_GPU_AWARE && echo yes || echo no) (PETSc's own MPI detection)"
# A reused prefix must have been built with THIS toolkit: a PETSc/Kokkos built
# against another CUDA (say 12.x vs 13.x) would load two libcudart versions.
PETSC_VARS="$(first_existing "$PREFIX/lib/petsc/conf/petscvariables" "$PREFIX/lib64/petsc/conf/petscvariables" || true)"
PETSC_CUDAC="$(sed -n 's/^CUDAC[[:space:]]*=[[:space:]]*\([^[:space:]]*\).*/\1/p' "$PETSC_VARS" 2>/dev/null | head -n 1)"
if [[ "$PETSC_CUDAC" == /* ]]; then
  [[ "$(readlink -f "$PETSC_CUDAC" 2>/dev/null || echo "$PETSC_CUDAC")" == "$(readlink -f "$CUDA_DIR/bin/nvcc")" ]] \
    || die "PETSc in $PREFIX was built with $PETSC_CUDAC, but this run selected $CUDA_DIR/bin/nvcc. Load the same CUDA module again (FREHG_MODULE_CUDA=...), or reset the prefix to rebuild against this one: $WIPE_HINT"
  echo "petscvariables: CUDAC = $PETSC_CUDAC (matches the selected toolkit)"
fi
HYPRE_CFG="$PREFIX/include/HYPRE_config.h"
HYPRE_AWARE=0
if [[ -f "$HYPRE_CFG" ]]; then
  grep -Eq '^#define HYPRE_USING_CUDA 1' "$HYPRE_CFG" \
    || die "hypre in $PREFIX was not built with CUDA (HYPRE_USING_CUDA missing). Reset the prefix: $WIPE_HINT"
  grep -Eq '^#define HYPRE_USING_GPU_AWARE_MPI 1' "$HYPRE_CFG" && HYPRE_AWARE=1
  echo "HYPRE_config.h: HYPRE_USING_CUDA = yes, HYPRE_USING_GPU_AWARE_MPI = $([[ $HYPRE_AWARE -eq 1 ]] && echo yes || echo no)"
  [[ "$HYPRE_AWARE" -le "$MPI_GPU_AWARE" ]] \
    || warn "hypre passes device buffers to MPI but the MPI is not CUDA-aware: multi-rank 'amg' runs will crash. Reset the prefix and rerun to rebuild hypre with --disable-gpu-aware-mpi."
else
  warn "HYPRE_config.h not found in $PREFIX/include; relying on PETSC_HAVE_HYPRE_DEVICE above."
fi
# A reused PETSc must still have its Kokkos Kernels runtime: petscconf.h
# advertises KOKKOS_KERNELS even after the library is gone ("Kokkos rebuilt
# but PETSc reused" -> dangling libpetsc, cryptic loader abort at first run).
first_existing "$PREFIX"/lib/libkokkoskernels.so "$PREFIX"/lib64/libkokkoskernels.so >/dev/null \
  || die "PETSc's Kokkos Kernels runtime is missing from $PREFIX (Kokkos was likely rebuilt without PETSc). Reset the prefix: $WIPE_HINT"

say "9 of 10: Configure and compile Frehg2 in $BUILD_NAME/"
export CMAKE_PREFIX_PATH="$PREFIX${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PREFIX/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# Start fresh: a CMake cache from another compiler, MPI or backend is unsafe.
rm -rf "$BUILD_ROOT"
TESTS="${FREHG_ENABLE_TESTS:-0}"
WERROR="${FREHG_WERROR:-ON}"
GPU_AWARE_DEFAULT="${FREHG_GPU_AWARE_MPI:-OFF}"
# Frehg2's CMakeLists requires HDF5_IS_PARALLEL, which only CMake's
# module-mode FindHDF5 sets; force module mode and prefer the parallel build.
HDF5_CMAKE_ARGS=( -DHDF5_NO_FIND_PACKAGE_CONFIG_FILE=TRUE -DHDF5_PREFER_PARALLEL=TRUE )
[[ "$USE_MODULE_HDF5" -eq 0 ]] && HDF5_CMAKE_ARGS+=( -DHDF5_ROOT="$PREFIX" )
# Device lane exactly as docs/agents/build.md s.7 and cuda-compile.yml:
# CMAKE_CXX_COMPILER=<prefix>/bin/nvcc_wrapper (host g++ via
# NVCC_WRAPPER_DEFAULT_COMPILER, never mpicxx) and FREHG_BACKEND=cuda, which
# makes CMake verify the Kokkos has CUDA and the PETSc has PETSC_HAVE_KOKKOS.
env CC="$REAL_CC" CXX="$PREFIX/bin/nvcc_wrapper" cmake -S "$ROOT_DIR" -B "$BUILD_ROOT" \
  -DCMAKE_CXX_COMPILER="$PREFIX/bin/nvcc_wrapper" \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" \
  -DCUDAToolkit_ROOT="$CUDA_DIR" \
  "${HDF5_CMAKE_ARGS[@]}" \
  -DFREHG_BACKEND=cuda -DFREHG_GPU_AWARE_MPI="$GPU_AWARE_DEFAULT" \
  -DFREHG_ENABLE_TESTS="$TESTS" -DFREHG_WERROR="$WERROR"
cmake --build "$BUILD_ROOT" --parallel "$JOBS"
[[ -x "$EXE" ]] || die "CMake finished but expected executable $EXE was not produced."
# Link evidence: every library resolves and exactly ONE libkokkoscore (ours).
LDD_OUT="$(ldd "$EXE" 2>&1 || true)"
# ldd also lists optional ELF filter entries (DT_AUXILIARY / DT_FILTER) as
# "=> not found". A compiler flag like -fvisibility=hidden that reaches GNU ld
# directly (via nvcc) is read as "-f visibility=hidden": the auxiliary filter
# "visibility=hidden". The loader skips a missing auxiliary filter, so only
# DT_NEEDED entries are real missing libraries; readelf tells them apart.
MISSING="$(awk '/=> not found/ {print $1}' <<<"$LDD_OUT" | sort -u)"
if [[ -n "$MISSING" ]]; then
  ELF_OBJS=( "$EXE" )
  while IFS= read -r lib; do ELF_OBJS+=( "$lib" ); done \
    < <(awk '/=> \// {print $3}' <<<"$LDD_OUT" | sort -u)
  HAVE_READELF=0
  command -v readelf >/dev/null 2>&1 && HAVE_READELF=1
  REAL_MISSING=()
  while IFS= read -r name; do
    [[ -n "$name" ]] || continue
    kind=unknown; owner=""
    if [[ "$HAVE_READELF" -eq 1 ]]; then
      for obj in "${ELF_OBJS[@]}"; do
        dyn="$(readelf -d "$obj" 2>/dev/null || true)"
        if grep -qF "Shared library: [$name]" <<<"$dyn"; then kind=needed; owner="$obj"; break; fi
        if grep -qF "Auxiliary library: [$name]" <<<"$dyn" || grep -qF "Filter library: [$name]" <<<"$dyn"; then
          kind=auxiliary; owner="$obj"
        fi
      done
    fi
    if [[ "$kind" == auxiliary ]]; then
      warn "ldd lists '$name => not found', but it is an optional ELF filter entry (DT_AUXILIARY/DT_FILTER) in $owner, not a library dependency; the loader ignores it. Continuing."
    else
      REAL_MISSING+=( "$name${owner:+ (needed by $owner)}" )
    fi
  done <<<"$MISSING"
  if [[ ${#REAL_MISSING[@]} -gt 0 ]]; then
    printf '  missing: %s\n' "${REAL_MISSING[@]}" >&2
    [[ "$HAVE_READELF" -eq 1 ]] || warn "readelf is unavailable, so filter entries could not be told apart from missing libraries."
    die "$EXE has unresolved shared libraries (see above). Check LD_LIBRARY_PATH and the loaded modules."
  fi
fi
KOKKOS_CORES="$(awk '/libkokkoscore/ {print $3}' <<<"$LDD_OUT" | sort -u)"
[[ -n "$KOKKOS_CORES" && "$(wc -l <<<"$KOKKOS_CORES")" -eq 1 ]] \
  || die "Expected exactly one shared libkokkoscore in $EXE, found: ${KOKKOS_CORES:-none} (two Kokkos runtimes => segfault)."
[[ "$(readlink -f "$KOKKOS_CORES")" == "$(readlink -f "$PREFIX")"/* ]] \
  || die "$EXE loads libkokkoscore from $KOKKOS_CORES, not from $PREFIX (a different Kokkos on LD_LIBRARY_PATH?)."
echo "Linked libraries (selection):"
grep -E 'kokkos|cudart|cublas|cusparse|petsc|HYPRE|libmpi|hdf5|yaml' <<<"$LDD_OUT" || true
echo "SUCCESS: executable created: $EXE"

say "10 of 10: Validate and run b1-sw on the GPU (1 rank)"
export OMP_NUM_THREADS=1
# Always launch through mpirun, even for --validate and 1 rank: this site's
# Open MPI has no SLURM PMI, and a bare singleton inside an srun shell is
# taken for an srun direct launch ("OMPI was not built with SLURM's PMI
# support"). One local rank needs no srun step.
MPI_LAUNCHER="$(command -v mpirun || command -v mpiexec || true)"
[[ -n "$MPI_LAUNCHER" ]] || die "No mpirun/mpiexec found."
PETSC_FORCED_NO_GPU_AWARE=0
if [[ "$MPI_GPU_AWARE" -eq 0 && "${PETSC_OPTIONS:-}" != *use_gpu_aware_mpi* ]]; then
  export PETSC_OPTIONS="${PETSC_OPTIONS:+$PETSC_OPTIONS }-use_gpu_aware_mpi 0"
  PETSC_FORCED_NO_GPU_AWARE=1
fi
echo "PETSC_OPTIONS=${PETSC_OPTIONS:-<unset>}"
# Work in copies under build-cuda/ so benchmarks/b1-sw/out (the CPU script's
# verification output) is never overwritten.
VERIFY_DIR="$BUILD_ROOT/verify-b1-sw"
rm -rf "$VERIFY_DIR"; cp -R "$ROOT_DIR/benchmarks/b1-sw" "$VERIFY_DIR"; rm -rf "$VERIFY_DIR/out"
run_case() { # run_case <log> <frehg args...>; 1 rank, output to screen and log
  local log="$1"; shift
  "$MPI_LAUNCHER" -np 1 "$EXE" "$@" 2>&1 | tee "$log"
}
pushd "$VERIFY_DIR" >/dev/null
# The first launch doubles as an mpirun health check. On partition intel this
# Open MPI segfaulted inside hwloc (2026-09-01); if that happens here, retry
# once with the known workaround and keep it for the remaining runs.
HWLOC_WORKAROUND=0
if ! run_case validate.log --validate b1-sw.yaml; then
  if [[ -n "$OMPI_INFO" ]] && grep -Eqi 'hwloc|rtc_hwloc|opal_hwloc|Segmentation fault|signal 11' validate.log; then
    warn "mpirun crashed in/around hwloc; retrying with OMPI_MCA_rtc=^hwloc OMPI_MCA_hwloc_base_binding_policy=none."
    export OMPI_MCA_rtc='^hwloc' OMPI_MCA_hwloc_base_binding_policy=none
    HWLOC_WORKAROUND=1
    run_case validate.log --validate b1-sw.yaml || die "--validate failed even with the hwloc workaround; see $VERIFY_DIR/validate.log."
  else
    die "--validate failed; see $VERIFY_DIR/validate.log."
  fi
fi
grep -q '^VALID:' validate.log || die "--validate did not print VALID:; see $VERIFY_DIR/validate.log."
# -device_view (a PETSc option, passed after the YAML) prints the GPU PETSc
# uses; it forces eager device init, which also runs PETSc's GPU-aware-MPI check.
if ! run_case run-gpu.log b1-sw.yaml -device_view; then
  if grep -q 'MPI is not GPU-aware' run-gpu.log && [[ "${PETSC_OPTIONS:-}" != *use_gpu_aware_mpi* ]]; then
    warn "The MPI reported CUDA support, but PETSc's run-time check says it is NOT GPU-aware here. Retrying with -use_gpu_aware_mpi 0; treat this MPI as not GPU-aware (frehg2 runtime.gpu_aware_mpi: off)."
    export PETSC_OPTIONS="${PETSC_OPTIONS:+$PETSC_OPTIONS }-use_gpu_aware_mpi 0"
    PETSC_FORCED_NO_GPU_AWARE=1; MPI_GPU_AWARE=0
    rm -rf out
    run_case run-gpu.log b1-sw.yaml -device_view \
      || die "GPU verification run failed again; see $VERIFY_DIR/run-gpu.log."
  else
    die "GPU verification run failed; see $VERIFY_DIR/run-gpu.log."
  fi
fi
grep -q 'Kokkos backend Cuda' run-gpu.log \
  || die "The run did not report 'Kokkos backend Cuda' at startup: the binary is not on the CUDA backend."
grep -q 'run complete' run-gpu.log || die "The run did not print 'run complete'; see $VERIFY_DIR/run-gpu.log."
RECORD="$VERIFY_DIR/out/run-record.yaml"
[[ -f "$RECORD" ]] || die "No run record at $RECORD."
grep -Eq 'kokkos_backend: *Cuda' "$RECORD" || die "Run record does not report kokkos_backend: Cuda."
grep -Eq 'mat_type: *(seq|mpi)aijkokkos' "$RECORD" \
  || die "Run record does not show a resolved aijkokkos solve (expected seqaijkokkos on 1 rank)."
popd >/dev/null
echo ""
echo "Evidence the device backend is active (from $VERIFY_DIR):"
grep -E 'Kokkos backend|device build \(|\[0\] name:|Compute capability:' "$VERIFY_DIR/run-gpu.log" || true
grep -E 'kokkos_backend:|mat_type: *(seq|mpi)aijkokkos' "$RECORD" || true

# Secondary, NON-fatal check: the same case with solver.surface.preconditioner
# amg, i.e. PETSc's CUDA hypre BoomerAMG on device (frehg2 pins PMIS
# coarsening on device builds). A failure here is reported, not fatal.
AMG_STATUS="not run"
AMG_DIR="$BUILD_ROOT/verify-b1-sw-amg"
rm -rf "$AMG_DIR"; cp -R "$ROOT_DIR/benchmarks/b1-sw" "$AMG_DIR"; rm -rf "$AMG_DIR/out"
printf '\nsolver:\n  surface:\n    preconditioner: amg\n' >> "$AMG_DIR/b1-sw.yaml"
pushd "$AMG_DIR" >/dev/null
if "$MPI_LAUNCHER" -np 1 "$EXE" --validate b1-sw.yaml >validate.log 2>&1 \
   && run_case run-gpu-amg.log b1-sw.yaml >/dev/null \
   && grep -q 'run complete' run-gpu-amg.log \
   && grep -Eq 'amg_coarsen_type: *PMIS' out/run-record.yaml; then
  AMG_STATUS="PASS (hypre BoomerAMG ran on the device; PMIS coarsening recorded)"
else
  AMG_STATUS="FAILED -- see $AMG_DIR/run-gpu-amg.log (build is still usable with the default bjacobi-icc solver)"
fi
popd >/dev/null
echo "hypre-on-GPU check: $AMG_STATUS"

LOAD_LINE="module purge; module load ${LOADED_MODULES[*]}"
echo ""
echo "================================================================="
echo "FREHG2 CUDA BUILD AND 1-GPU VERIFICATION COMPLETED SUCCESSFULLY"
echo "Executable:        $EXE"
echo "Dependency prefix: $PREFIX"
echo "Backend:           Kokkos Cuda ($KOKKOS_ARCH, sm_$ARCH_NUM), CUDA $NVCC_VERSION, PETSc $PETSC_VERSION (CUDA+Kokkos+hypre-CUDA)"
echo "GPU-aware MPI:     $([[ $MPI_GPU_AWARE -eq 1 ]] && echo yes || echo no) ($MPI_GPU_AWARE_HOW)"
echo "hypre on GPU:      $AMG_STATUS"
if [[ "$HYPRE_AWARE" -eq 1 && "$MPI_GPU_AWARE" -eq 0 ]]; then
  echo "WARNING: hypre was built to pass device buffers to MPI, but this MPI is not"
  echo "  GPU-aware: multi-rank 'amg' runs will crash. Reset the prefix and rerun with"
  echo "  FREHG_MPI_GPU_AWARE=0 to rebuild hypre with --disable-gpu-aware-mpi:"
  echo "  $WIPE_HINT"
fi
echo "STATUS: frehg2's GPU lane is EXPERIMENTAL. This run shows the binary"
echo "executes on the GPU; it does not validate GPU physics. Compare against"
echo "CPU results (or run scripts/gpu_acceptance.sh) before trusting them."
echo "-----------------------------------------------------------------"
echo "RUN-TIME ENVIRONMENT (put this in every GPU job script):"
echo "  #SBATCH --partition=${GPU_PARTITION_HINT}"
echo "  #SBATCH --gres=gpu:N           # N GPUs"
echo "  #SBATCH --ntasks=N             # one MPI rank per GPU"
echo "  $LOAD_LINE"
echo "  export LD_LIBRARY_PATH=$PREFIX/lib:$PREFIX/lib64:$CUDA_DIR/lib64:\$LD_LIBRARY_PATH"
echo "  export OMP_NUM_THREADS=1        # host threads; device kernels unaffected"
if [[ "$PETSC_FORCED_NO_GPU_AWARE" -eq 1 || "$MPI_GPU_AWARE" -eq 0 ]]; then
  echo "  export PETSC_OPTIONS=\"-use_gpu_aware_mpi 0\"   # REQUIRED: this MPI is not GPU-aware"
fi
if [[ "$HWLOC_WORKAROUND" -eq 1 ]]; then
  echo "  export OMPI_MCA_rtc=^hwloc OMPI_MCA_hwloc_base_binding_policy=none   # REQUIRED: mpirun crashes in hwloc here"
fi
echo "  cd <case-dir>"
echo "  mpirun -np 1 $EXE <case>.yaml            # single GPU"
echo "Multi-GPU (experimental):"
echo "  - one rank per GPU: mpirun -np N $EXE <case>.yaml"
echo "    Kokkos binds a rank to visible GPU (node-local rank % #GPUs, read from"
echo "    OMPI_COMM_WORLD_LOCAL_RANK) and PETSc to (rank % #GPUs). These agree on one"
echo "    node. Across nodes, place exactly #GPUs-per-node ranks per node"
echo "    (mpirun --map-by ppr:<gpus-per-node>:node). To pin one GPU per rank"
echo "    explicitly when SLURM numbers the job's GPUs from 0:"
echo "      mpirun -np N bash -c 'export CUDA_VISIBLE_DEVICES=\$OMPI_COMM_WORLD_LOCAL_RANK; exec $EXE <case>.yaml'"
if [[ "$MPI_GPU_AWARE" -eq 1 ]]; then
  echo "  - halo exchange: in the case YAML set   runtime: {gpu_aware_mpi: on}"
  echo "    to pass device buffers to MPI ('off' stages through the host and is safe)."
else
  echo "  - halo exchange: keep runtime.gpu_aware_mpi at auto/off (this MPI is not GPU-aware)."
fi
echo "Solvers: device builds force solver.*.mat_type to aijkokkos; 'amg' uses CUDA"
echo "  hypre (PMIS coarsening on device). Small grids are launch-latency bound"
echo "  (the docs suggest ~1e6 cells before a GPU pays off)."
if [[ "$HWLOC_WORKAROUND" -eq 0 ]]; then
  echo "If mpirun ever segfaults in hwloc on these nodes (seen on partition intel):"
  echo "  export OMPI_MCA_rtc=^hwloc OMPI_MCA_hwloc_base_binding_policy=none"
fi
echo "Owner GPU acceptance bundle (needs the CPU build in build/ from the same commit):"
echo "  scripts/gpu_acceptance.sh --frehg-gpu $EXE --frehg-cpu $ROOT_DIR/build/src/frehg \\"
echo "    --mpiexec \$(command -v mpiexec) --gpus <N> --work gpu-acceptance"
echo "Modules loaded for this build:"
module list 2>&1 || true
echo "================================================================="
