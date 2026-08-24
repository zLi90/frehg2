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

# ---- PETSc (minimal: MPI + BLAS/LAPACK, no Fortran) ------------------------
cd "$PREFIX/src"
curl -fsSL -o petsc.tar.gz \
  "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PETSC_VERSION}.tar.gz"
tar xf petsc.tar.gz
cd "petsc-${PETSC_VERSION}"
./configure --prefix="$PREFIX" \
  --with-fc=0 \
  --with-debugging=0 \
  --download-f2cblaslapack \
  COPTFLAGS=-O2 CXXOPTFLAGS=-O2
make -j "$JOBS" all
make install

touch "$PREFIX/.complete"
echo "dependencies installed at $PREFIX"
