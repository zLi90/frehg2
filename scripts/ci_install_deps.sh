#!/usr/bin/env bash
# ci_install_deps.sh — build the from-source dependencies (Kokkos, PETSc)
# into a prefix for CI. System packages (MPI, parallel HDF5, yaml-cpp,
# doxygen) come from apt in the workflow. The prefix is cached by the
# workflows keyed on the version pins below.

set -euo pipefail
PREFIX="${1:-$HOME/.frehg-deps}"
# Kokkos and Kokkos Kernels are released in lockstep and MUST be the identical
# version, compiled against each other. This version therefore has to equal the
# kokkos/kokkos-kernels gitcommit that the target PETSc release pins (PETSc
# 3.25.1 pins 5.1.0 in config/BuildSystem/config/packages/kokkos*.py). Bumping
# PETSc means re-checking that pin. See the Kokkos Kernels build below.
KOKKOS_VERSION="${KOKKOS_VERSION:-5.1.0}"
PETSC_VERSION="${PETSC_VERSION:-3.25.1}"
JOBS="$(getconf _NPROCESSORS_ONLN)"

# Memory-aware compile-parallelism cap for the explicit-template-instantiation
# (ETI) heavyweights: our own Kokkos Kernels build below, and hypre inside
# PETSc's configure (--with-make-np, otherwise PETSc auto-detects every core
# and runs e.g. "-j4"). Each Kokkos Kernels ETI translation unit peaks at
# several GB of compiler memory; on a memory-limited box (a 16 GB CI runner, a
# laptop) an unbounded -j OOM-kills cc1plus. Cap it ~8 GB per concurrent
# compile (a heavy ETI TU at -O2 can approach that), never above $JOBS, at
# least 1. Override with PETSC_MAKE_NP for a machine with plenty of RAM/core.
if [ -r /proc/meminfo ]; then
  MEM_GB="$(awk '/MemTotal/ {print int($2/1024/1024)}' /proc/meminfo)"
elif command -v sysctl >/dev/null 2>&1; then
  MEM_GB="$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1024 / 1024 / 1024 ))"
else
  MEM_GB=0
fi
MEM_NP="$(( MEM_GB / 8 ))"
[ "$MEM_NP" -lt 1 ] && MEM_NP=1
PETSC_MAKE_NP="${PETSC_MAKE_NP:-$(( MEM_NP < JOBS ? MEM_NP : JOBS ))}"
echo "using JOBS=$JOBS, PETSC_MAKE_NP=$PETSC_MAKE_NP (detected ${MEM_GB} GB RAM)"

if [ -f "$PREFIX/.complete" ]; then
  echo "dependencies already present at $PREFIX (cache hit)"
  exit 0
fi
mkdir -p "$PREFIX/src"

# ---- Kokkos (OpenMP + Serial host backends; +CUDA when FREHG_CUDA=1) ------
cd "$PREFIX/src"
curl -fsSL -o kokkos.tar.gz \
  "https://github.com/kokkos/kokkos/archive/refs/tags/${KOKKOS_VERSION}.tar.gz"
tar xf kokkos.tar.gz
# The CUDA compile-only lane (plan §10 P5) builds Kokkos with the CUDA
# backend through nvcc_wrapper; no GPU is needed to compile. The device
# architecture is a compile target only (FREHG_CUDA_ARCH, default AMPERE80).
KOKKOS_EXTRA=()
if [ "${FREHG_CUDA:-0}" = "1" ]; then
  KOKKOS_EXTRA+=(
    -DKokkos_ENABLE_CUDA=ON
    "-DKokkos_ARCH_${FREHG_CUDA_ARCH:-AMPERE80}=ON"
    "-DCMAKE_CXX_COMPILER=$PREFIX/src/kokkos-${KOKKOS_VERSION}/bin/nvcc_wrapper")
fi
# BUILD_SHARED_LIBS=ON is load-bearing (A8 invariant 7): frehg links Kokkos and
# so does PETSc (via the Kokkos Kernels built below, linked into libpetsc.so). A
# static Kokkos core is absorbed into BOTH, giving two copies of Kokkos's
# runtime singleton in one process -- it initializes twice and the first
# VecKokkos access segfaults. A single shared libkokkoscore keeps one runtime.
# It also produces -fPIC objects, without which the static Kokkos Kernels
# archive cannot be linked into the shared libpetsc.so (R_X86_64_TPOFF32 TLS
# relocation error on the ETI translation units).
cmake -S "kokkos-${KOKKOS_VERSION}" -B kokkos-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_CXX_STANDARD=20 \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  "${KOKKOS_EXTRA[@]}"
cmake --build kokkos-build -j "$JOBS"
cmake --install kokkos-build

# ---- Kokkos Kernels (SAME version as Kokkos, built against it) -------------
# PETSc's --download-kokkos-kernels is NOT usable here: it git-checks-out the
# version PETSc pins (5.1.0) and compiles it against whatever Kokkos we point
# --with-kokkos-dir at. Because A8 invariant 7 forces PETSc to reuse the single
# Kokkos we built above (a second Kokkos runtime would double-init the
# singleton and segfault), a Kokkos/Kokkos-Kernels version skew there is a
# lockstep violation and fails to compile. So we build Kokkos Kernels ourselves
# at ${KOKKOS_VERSION} against that same Kokkos and hand PETSc a pre-built one
# via --with-kokkos-kernels-dir below. Bonus: any compile error lands on CI
# stdout instead of being buried in PETSc's configure.log.
cd "$PREFIX/src"
curl -fsSL -o kokkos-kernels.tar.gz \
  "https://github.com/kokkos/kokkos-kernels/archive/refs/tags/${KOKKOS_VERSION}.tar.gz"
tar xf kokkos-kernels.tar.gz
KK_EXTRA=()
if [ "${FREHG_CUDA:-0}" = "1" ]; then
  # Match the Kokkos build's device compiler; the CUDA backend is inherited
  # from the installed Kokkos, so no Kokkos_ENABLE_CUDA flag is needed here.
  KK_EXTRA+=("-DCMAKE_CXX_COMPILER=$PREFIX/src/kokkos-${KOKKOS_VERSION}/bin/nvcc_wrapper")
fi
# Host BLAS/LAPACK TPLs off: PETSc supplies its own (f2cblaslapack); letting
# Kokkos Kernels link a different system BLAS invites ABI/link skew. ETI is left
# at the Kokkos Kernels defaults (double / LayoutLeft), which is what PETSc's
# aijkokkos backend uses. Build parallelism is the memory-aware cap (ETI TUs
# peak at several GB each at -O2 -- the exact OOM hazard flagged above).
cmake -S "kokkos-kernels-${KOKKOS_VERSION}" -B kokkos-kernels-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DBUILD_SHARED_LIBS=ON \
  -DCMAKE_CXX_STANDARD=20 \
  -DKokkos_ROOT="$PREFIX" \
  -DKokkosKernels_ENABLE_TPL_BLAS=OFF \
  -DKokkosKernels_ENABLE_TPL_LAPACK=OFF \
  "${KK_EXTRA[@]}"
cmake --build kokkos-kernels-build -j "$PETSC_MAKE_NP"
cmake --install kokkos-kernels-build

# ---- PETSc (MPI + BLAS/LAPACK + hypre BoomerAMG + Kokkos backend, no
# Fortran; v2 plan §2.2 and §2B.2 B3) -----------------------------------------
# PETSc reuses the single Kokkos AND the Kokkos Kernels installed above (OpenMP
# host backend on CPU lanes, CUDA on the FREHG_CUDA=1 lane), which is what makes
# -mat_type aijkokkos / solver.*.mat_type=aijkokkos real. The CUDA lane adds
# --with-cuda; compiling PETSc's CUDA support needs no GPU.
cd "$PREFIX/src"
curl -fsSL -o petsc.tar.gz \
  "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz"
tar xf petsc.tar.gz
cd "petsc-${PETSC_VERSION}"
PETSC_EXTRA=()
if [ "${FREHG_CUDA:-0}" = "1" ]; then
  PETSC_EXTRA+=(--with-cuda=1)
fi
# On any configure failure, surface the real error: PETSc redirects its
# downloaded-package build output (hypre, ...) into configure.log, so a failed
# package build otherwise shows only a generic "Error running make on <PKG>" on
# stdout. Dump the tail and the actual compiler diagnostics so CI logs
# distinguish an OOM ("Killed" / "cc1plus") from a genuine compile error
# without needing the runner. (Kokkos Kernels is built above, not here, so its
# errors already reach stdout directly.)
dump_petsc_log() {
  local log="$PWD/configure.log"
  [ -f "$log" ] || return 0
  echo "======== configure.log: compiler errors / OOM markers ========"
  grep -nE 'error:|fatal error|Killed|cannot allocate|out of memory|virtual memory exhausted' \
    "$log" | tail -n 60 || true
  echo "======== configure.log: last 100 lines ========"
  tail -n 100 "$log" || true
}
trap 'rc=$?; if [ "$rc" -ne 0 ]; then dump_petsc_log; fi' ERR
./configure --prefix="$PREFIX" \
  --with-fc=0 \
  --with-debugging=0 \
  --download-f2cblaslapack \
  --download-hypre \
  --with-kokkos-dir="$PREFIX" \
  --with-kokkos-kernels-dir="$PREFIX" \
  --with-openmp=1 \
  --with-make-np="$PETSC_MAKE_NP" \
  "${PETSC_EXTRA[@]}" \
  COPTFLAGS=-O2 CXXOPTFLAGS=-O2
make -j "$JOBS" all
make install
trap - ERR

touch "$PREFIX/.complete"
echo "dependencies installed at $PREFIX"
