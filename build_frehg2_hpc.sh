#!/usr/bin/env bash
# Frehg2 one-command private build using this HPC's existing HDF5 module.
#
# RUN THIS ONE COMMAND from the Frehg2 top-level directory:
#   bash build_frehg2_once_modules.sh
#
# What this script does automatically:
#   - loads the compatible HPC modules listed below;
#   - uses the HPC's HDF5/1.14.6 module (does NOT download/build HDF5);
#   - shallow-clones pinned yaml-cpp, Kokkos and PETSc source releases;
#   - builds them privately under ./frehg-deps;
#   - compiles build/src/frehg;
#   - validates the b1-sw input; and
#   - creates, submits, waits for, and checks a short SLURM MPI verification.
#
# Nothing is installed under /share/apps or another system/public package
# location. The private project-local prefix is ./frehg-deps.
#
# The four exact modules below came from this HPC's module output. HDF5/1.14.6
# requires OpenMPI/4.1.6, so all MPI-dependent software is deliberately built
# and run with OpenMPI/4.1.6. Do not change just one of these modules.

set -Eeuo pipefail
umask 077
ROOT_DIR="$PWD"
PREFIX="$ROOT_DIR/frehg-deps"
CACHE="$PREFIX/src-cache"
JOBS=2
STAMP="$(date +%Y%m%d-%H%M%S)"
LOG="$ROOT_DIR/frehg2-build-$STAMP.log"
trap 'status=$?; echo; echo "BUILD FAILED at line $LINENO (exit $status). Full log: $LOG"; exit "$status"' ERR
exec > >(tee -a "$LOG") 2>&1

say() { printf '\n==================== %s ====================\n' "$*"; }
die() { echo "ERROR: $*"; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "Required command '$1' is unavailable after loading modules."; }

[[ -f "$ROOT_DIR/CMakeLists.txt" ]] || die "Run this script from the Frehg2 top-level directory containing CMakeLists.txt."
[[ -d "$ROOT_DIR/benchmarks/b1-sw" ]] || die "benchmarks/b1-sw is missing; this is not a complete Frehg2 checkout."
mkdir -p "$PREFIX" "$CACHE"
chmod 700 "$PREFIX" "$CACHE" || true

say "Frehg2 private one-command build"
echo "Repository: $ROOT_DIR"
echo "Private dependency directory: $PREFIX"
echo "Build log: $LOG"
echo "Compilation jobs: $JOBS (safe login-node setting)"

say "1/7: load the matched compiler, MPI, and HDF5 modules"
if ! type module >/dev/null 2>&1; then
  [[ -r /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  [[ -r /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
fi
type module >/dev/null 2>&1 || die "The 'module' command is not available in this shell."
module purge
module load cmake/3.31.6
module load gcc/13.2.0
module load openmpi/4.1.6
module load hdf5/1.14.6
need cmake; need gcc; need g++; need mpicc; need mpicxx; need curl; need sbatch; need squeue
printf 'CMake: %s\nC++ compiler: %s\nMPI C wrapper: %s\n' "$(cmake --version | awk 'NR==1 {print $3}')" "$(g++ --version | head -1)" "$(command -v mpicc)"
module list 2>&1 || true

say "2/7: obtain source releases (git clone, then tarball fallback)"
# clone_source GIT_URL TAG DESTINATION TARBALL_URL TARBALL_FILE
clone_source() {
  local git_url="$1" tag="$2" dest="$3" tar_url="$4" tar_file="$5"
  if [[ -d "$dest" ]]; then echo "Reusing source tree: $dest"; return 0; fi
  if command -v git >/dev/null 2>&1; then
    echo "git clone --depth 1 --branch $tag $git_url"
    if git clone --depth 1 --branch "$tag" "$git_url" "$dest"; then return 0; fi
    echo "Git clone failed; using the identical pinned release archive."
    rm -rf "$dest"
  fi
  [[ -f "$tar_file" ]] || curl --fail --location --retry 3 --output "$tar_file" "$tar_url"
  local extracted="$(tar -tf "$tar_file" | head -1 | cut -d/ -f1)"
  tar -xf "$tar_file" -C "$CACHE"
  [[ -d "$CACHE/$extracted" ]] || die "Could not locate extracted sources from $tar_file."
  mv "$CACHE/$extracted" "$dest"
}
clone_source https://github.com/jbeder/yaml-cpp.git 0.8.0 "$CACHE/yaml-cpp-0.8.0" https://github.com/jbeder/yaml-cpp/archive/refs/tags/0.8.0.tar.gz "$CACHE/yaml-cpp-0.8.0.tar.gz"
clone_source https://github.com/kokkos/kokkos.git 5.1.1 "$CACHE/kokkos-5.1.1" https://github.com/kokkos/kokkos/archive/refs/tags/5.1.1.tar.gz "$CACHE/kokkos-5.1.1.tar.gz"
clone_source https://gitlab.com/petsc/petsc.git v3.25.1 "$CACHE/petsc-3.25.1" https://gitlab.com/petsc/petsc/-/archive/v3.25.1/petsc-v3.25.1.tar.gz "$CACHE/petsc-3.25.1.tar.gz"

cmake_install() {
  local source="$1" build="$2"; shift 2
  cmake -S "$source" -B "$build" -G "Unix Makefiles" "$@"
  cmake --build "$build" --parallel "$JOBS"
  cmake --install "$build"
}

say "3/7: build yaml-cpp and Kokkos under ./frehg-deps"
if [[ ! -f "$PREFIX/lib64/cmake/yaml-cpp/yaml-cpp-config.cmake" && ! -f "$PREFIX/lib/cmake/yaml-cpp/yaml-cpp-config.cmake" ]]; then
  rm -rf "$CACHE/build-yaml-cpp-0.8.0"
  cmake_install "$CACHE/yaml-cpp-0.8.0" "$CACHE/build-yaml-cpp-0.8.0" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_CXX_COMPILER="$(command -v g++)" -DYAML_BUILD_SHARED_LIBS=ON -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
else echo "yaml-cpp already installed; reusing it."; fi
if [[ ! -f "$PREFIX/lib64/cmake/Kokkos/KokkosConfig.cmake" && ! -f "$PREFIX/lib/cmake/Kokkos/KokkosConfig.cmake" ]]; then
  rm -rf "$CACHE/build-kokkos-5.1.1"
  cmake_install "$CACHE/kokkos-5.1.1" "$CACHE/build-kokkos-5.1.1" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_CXX_COMPILER="$(command -v g++)" -DBUILD_SHARED_LIBS=ON -DCMAKE_CXX_STANDARD=20 -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_TESTS=OFF
else echo "Kokkos already installed; reusing it."; fi

say "4/7: build PETSc under ./frehg-deps with OpenMPI 4.1.6"
if [[ ! -f "$PREFIX/lib/petsc/conf/petscvariables" && ! -f "$PREFIX/lib64/petsc/conf/petscvariables" ]]; then
  pushd "$CACHE/petsc-3.25.1" >/dev/null
  # Kokkos-enabled PETSc (v2 plan §2B.2 B3): enables mat_type aijkokkos.
  # On a GPU node, add --with-cuda=1 (or --with-hip=1) and build the step-3
  # Kokkos with the matching device backend first — everything else is
  # unchanged (v2 plan §2B.2; docs/developer-guide/gpu-acceptance.md).
  ./configure --prefix="$PREFIX" --with-cc="$(command -v mpicc)" --with-cxx="$(command -v mpicxx)" --with-fc=0 --with-debugging=0 --download-f2cblaslapack --download-hypre --with-kokkos-dir="$PREFIX" --download-kokkos-kernels --with-openmp=1
  make -j "$JOBS" all
  make install
  popd >/dev/null
else echo "PETSc already installed; reusing it."; fi

say "5/7: configure and compile Frehg2"
export CMAKE_PREFIX_PATH="$PREFIX${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
export PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PREFIX/lib64/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
rm -rf "$ROOT_DIR/build"
# Do NOT replace gcc/g++ below with mpicc/mpicxx. CMake finds MPI itself.
env CC="$(command -v gcc)" CXX="$(command -v g++)" cmake -S "$ROOT_DIR" -B "$ROOT_DIR/build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH" -DHDF5_PREFER_PARALLEL=TRUE -DFREHG_ENABLE_TESTS=OFF -DFREHG_WERROR=ON
cmake --build "$ROOT_DIR/build" --parallel "$JOBS"
EXE="$ROOT_DIR/build/src/frehg"
[[ -x "$EXE" ]] || die "The expected executable was not created: $EXE"
echo "SUCCESS: executable created: $EXE"

say "6/7: validate a Frehg2 benchmark input"
export OMP_NUM_THREADS=1
"$EXE" --validate "$ROOT_DIR/benchmarks/b1-sw/b1-sw.yaml"

say "7/7: submit and wait for a SLURM MPI verification"
VERIFY="$ROOT_DIR/run_frehg2_verification.slurm"
cat > "$VERIFY" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=frehg2-verify
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=$ROOT_DIR/frehg2-verify-%j.out
#SBATCH --error=$ROOT_DIR/frehg2-verify-%j.err
set -Eeuo pipefail
module purge
module load cmake/3.31.6
module load gcc/13.2.0
module load openmpi/4.1.6
module load hdf5/1.14.6
export PATH="$PREFIX/bin:\$PATH"
export LD_LIBRARY_PATH="$PREFIX/lib:$PREFIX/lib64\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1
cd "$ROOT_DIR/benchmarks/b1-sw"
srun --ntasks=2 --cpus-per-task=1 "$EXE" b1-sw.yaml
EOF
chmod 700 "$VERIFY"
JOB_ID="$(sbatch --parsable "$VERIFY" | cut -d';' -f1)"
echo "Submitted SLURM verification job: $JOB_ID. Waiting for it to finish..."
while squeue -h -j "$JOB_ID" 2>/dev/null | grep -q .; do sleep 10; done
sleep 3
OUT="$ROOT_DIR/frehg2-verify-$JOB_ID.out"
ERR="$ROOT_DIR/frehg2-verify-$JOB_ID.err"
[[ -f "$OUT" ]] && cat "$OUT"
[[ ! -s "$ERR" ]] || { echo "----- SLURM error file -----"; cat "$ERR"; die "SLURM verification wrote errors."; }
STATE="$(sacct -n -X -j "$JOB_ID" --format=State 2>/dev/null | awk 'NF {print $1; exit}' || true)"
[[ -z "$STATE" || "$STATE" == COMPLETED ]] || die "SLURM reports verification job state: $STATE"
chmod -R go-rwx "$PREFIX" || true

echo
echo "================================================================"
echo "ALL FREHG2 BUILD AND VERIFICATION TASKS COMPLETED SUCCESSFULLY"
echo "Executable: $EXE"
echo "Private project dependency directory: $PREFIX"
echo "Build log: $LOG"
echo "================================================================"
