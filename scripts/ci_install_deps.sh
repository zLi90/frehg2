#!/usr/bin/env bash
# ci_install_deps.sh — build the from-source dependencies (Kokkos, PETSc)
# into a prefix for CI. System packages (MPI, parallel HDF5, yaml-cpp,
# doxygen) come from apt in the workflow. The prefix is cached by the
# workflows keyed on the version pins below.

set -euo pipefail
PREFIX="${1:-$HOME/.frehg-deps}"
KOKKOS_VERSION="${KOKKOS_VERSION:-5.1.1}"
PETSC_VERSION="${PETSC_VERSION:-3.25.1}"
JOBS="$(getconf _NPROCESSORS_ONLN)"

# Parallelism for PETSc's OWN downloaded-package builds (kokkos-kernels,
# hypre) is separate from the "make -j $JOBS" that builds PETSc's sources:
# PETSc's configure otherwise auto-detects every core and runs e.g. "-j4" on
# kokkos-kernels, whose explicit-template-instantiation (ETI) translation
# units each peak at several GB of compiler memory. On a memory-limited box
# (a 16 GB CI runner, a laptop) that OOM-kills cc1plus -- the build dies with
# a bare "Error running make on KOKKOS-KERNELS". Cap it memory-aware:
# ~8 GB per concurrent compile (a heavy ETI TU at -O2 can approach that),
# never above $JOBS, at least 1. Override with PETSC_MAKE_NP for a machine
# with plenty of RAM per core.
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
cmake -S "kokkos-${KOKKOS_VERSION}" -B kokkos-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_CXX_STANDARD=20 \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  "${KOKKOS_EXTRA[@]}"
cmake --build kokkos-build -j "$JOBS"
cmake --install kokkos-build

# ---- PETSc (MPI + BLAS/LAPACK + hypre BoomerAMG + Kokkos backend, no
# Fortran; v2 plan §2.2 and §2B.2 B3) -----------------------------------------
# PETSc reuses the Kokkos installed above (OpenMP host backend on CPU lanes,
# CUDA on the FREHG_CUDA=1 lane) and builds Kokkos Kernels against it, which
# is what makes -mat_type aijkokkos / solver.*.mat_type=aijkokkos real. The
# CUDA lane adds --with-cuda; compiling PETSc's CUDA support needs no GPU.
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
# downloaded-package build output (kokkos-kernels, hypre, ...) into
# configure.log, so a failed package build otherwise shows only the generic
# "Error running make on KOKKOS-KERNELS" on stdout. Dump the tail and the
# actual compiler diagnostics so CI logs distinguish an OOM ("Killed" /
# "cc1plus") from a genuine compile error without needing the runner.
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
  --download-kokkos-kernels \
  --with-openmp=1 \
  --with-make-np="$PETSC_MAKE_NP" \
  "${PETSC_EXTRA[@]}" \
  COPTFLAGS=-O2 CXXOPTFLAGS=-O2
make -j "$JOBS" all
make install
trap - ERR

touch "$PREFIX/.complete"
echo "dependencies installed at $PREFIX"
