# Build guide for AI assistants

**Purpose:** everything an LLM needs to write automated install / build /
compile scripts for Frehg2 on a fresh system — desktop or HPC, Linux or
macOS, CPU or (compile-only) GPU — without reading the source code.
Facts here were established during the gated P0–P5 development and are
enforced by the repository's own CI scripts; where a fact lives in an
executable file, that file is named so your script can defer to it.

## 1. What you are building

One executable, **`build/src/frehg`**, from five static C++20 libraries.
Build system: CMake ≥ 3.23 with plain `find_package` discovery — no
submodules, no vendored dependencies, no network access at configure
time (except the test suite, which fetches GoogleTest; disable tests to
build fully offline).

## 2. Dependency matrix

| Dependency | Minimum | Pinned known-good | Notes |
|---|---|---|---|
| CMake | 3.23 | — | |
| C++ compiler | gcc ≥ 12 or clang ≥ 15 | gcc 13–15, Apple clang 15 tested | must support `-std=c++20` |
| MPI | any MPI-3 | MPICH 4.3, Open MPI | C++ bindings required (`find_package(MPI COMPONENTS CXX)`) |
| Kokkos | 5.1 | **5.1.1** | Serial + OpenMP host backends; CUDA optional |
| PETSc | 3.20 | **3.25.1** | KSP only; `--with-fc=0` (no Fortran) is fine |
| HDF5 | 1.10, **parallel** | 1.14.x | must be built `--enable-parallel` against the *same* MPI |
| yaml-cpp | 0.7 | 0.7–0.8 | |

Python (for the test harness and post-processing only): `python3` with
`yaml`, `h5py`, `numpy`.

## 3. Invariants — get these right or the build fails confusingly

1. **Every dependency must be built against the same MPI.** A serial
   HDF5 will `find_package` successfully but Frehg2's configure aborts
   with an explicit `HDF5_IS_PARALLEL` error; mixed MPIs fail at output
   time. On Ubuntu, install the `-mpich-` *or* the `-openmpi-` flavor of
   HDF5 consistently with the MPI you install.
2. **Never set `CMAKE_CXX_COMPILER=mpicxx`** (or any MPI wrapper). The
   wrapper injects plain `-I` include paths that defeat CMake's
   SYSTEM-header treatment, and the strict `-Wconversion -Werror` build
   then fails *inside Kokkos headers*. Set `CXX`/`CC` to the real
   compiler; CMake adds MPI flags itself via `find_package(MPI)`.
3. **Kokkos and PETSc must be consumed by the same compiler family that
   built them.** If the user has gcc-built deps, build Frehg2 with gcc.
4. **`CMAKE_PREFIX_PATH` points at the dependency install prefix** (the
   directory containing `lib/cmake/Kokkos/KokkosConfig.cmake`), never at
   Frehg2's `build/` directory.
5. The default build is **zero-warning with `-Werror`**
   (`-Wall -Wextra -Wpedantic -Wshadow -Wconversion`). On a compiler
   version outside the tested set, a spurious warning can stop the
   build: the escape hatch is `-DFREHG_WERROR=OFF` (report the warning
   upstream rather than patching source).
6. macOS: gcc on Apple Silicon ships **no sanitizer runtimes** — a
   `FREHG_SANITIZE=ON` build there must use Apple clang (and clang-built
   Kokkos/yaml-cpp). Plain builds work with either compiler.

## 4. Install routes (pick by system)

### Route A — package manager + pinned source builds (recommended)

System packages provide compiler, CMake, MPI, parallel HDF5, yaml-cpp;
the repository's canonical script builds pinned Kokkos + PETSc into a
user-chosen prefix (idempotent; drops a `.complete` marker):

```bash
# Debian/Ubuntu
sudo apt install cmake g++ mpich libmpich-dev libhdf5-mpich-dev libyaml-cpp-dev
# macOS (Homebrew) — hdf5-mpi is the parallel HDF5
brew install cmake open-mpi hdf5-mpi yaml-cpp

scripts/ci_install_deps.sh $HOME/frehg-deps      # builds Kokkos 5.1.1 + PETSc 3.25.1
```

Script interface: prefix as `$1` (default `~/.frehg-deps`); env
overrides `KOKKOS_VERSION`, `PETSC_VERSION`; `FREHG_CUDA=1` adds the
CUDA backend (§7); needs `curl`, internet, and ~10–30 min. PETSc builds
its own BLAS/LAPACK (`--download-f2cblaslapack`), so none is needed from
the system.

### Route B — everything from system packages / Spack / modules

If the distribution or an HPC module system already provides Kokkos
≥ 5.1 and PETSc ≥ 3.20 (recent Spack and Fedora do), skip the script.
Load the modules so the packages land on `CMAKE_PREFIX_PATH`, and honor
invariant 3 (module compiler = build compiler). On Cray systems prefer
`CC`-the-real-compiler plus MPICH from the environment over the `CC`
wrapper (invariant 2 applies to Cray wrappers too).

### Route C — fully manual

Follow Route A's script as a recipe (it is short and readable): Kokkos
with `-DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_SERIAL=ON
-DCMAKE_CXX_STANDARD=20`; PETSc with `--with-fc=0 --with-debugging=0`.

## 5. Configure, build, verify

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/frehg-deps
cmake --build build -j
```

Where the dependency compiler is not the system default (e.g. macOS with
gcc-built deps at a custom prefix):

```bash
env CXX=g++-15 CC=gcc-15 CMAKE_PREFIX_PATH=/path/to/deps cmake -B build -S .
```

CMake options:

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `RelWithDebInfo` | use `Release` for production |
| `FREHG_WERROR` | `ON` | warnings are errors (escape hatch: `OFF`) |
| `FREHG_ENABLE_TESTS` | `ON` | test suite (fetches GoogleTest at configure; `OFF` for offline/production-only builds) |
| `FREHG_SANITIZE` | `OFF` | ASan+UBSan instrumented build |
| `FREHG_GPU_AWARE_MPI` | `OFF` | compile-time default for GPU-aware MPI staging (runtime-overridable) |
| `FREHG_LEGACY_BENCHMARKS` | `../legacy/benchmarks` | path to regression goldens — **not distributed**; see §6 |

**Verification (no external data needed):**

```bash
ctest --test-dir build -L 'unit|mpi' --output-on-failure   # must be 100 % pass
OMP_NUM_THREADS=1 build/src/frehg --validate benchmarks/b1-sw/b1-sw.yaml
cd benchmarks/b1-sw && OMP_NUM_THREADS=1 ../../build/src/frehg b1-sw.yaml  # ~0.3 s, ends "run complete"
```

For a deeper acceptance test, run cases from `validation/` — each README
states the expected physical result and wall time (seconds to minutes
for most).

## 6. What an outside user CANNOT run

The `regression`/`regression_nightly` ctest labels compare against the
legacy model's golden outputs, which are deliberately **not
distributed**. Every entry point (`scripts/ci_build_and_test.sh`, the CI
workflows, the sanitizer lane) detects their absence and falls back
loudly to the self-contained subset — the restart-determinism and
rank-invariance gates run fine without goldens. Do not treat the skip
message as a broken install, and do not try to fabricate goldens.

## 7. GPU (CUDA) builds — compile-only support

Frehg2 is **GPU-ready, GPU-unvalidated**: the CUDA backend compiles
warning-free (CI-gated) but no benchmark has been validated on a device.
Scripts may offer a CUDA build as *experimental*, with that caveat
stated. Recipe:

```bash
export CC=gcc-12 CXX=g++-12          # host gcc must satisfy the local nvcc's support matrix
FREHG_CUDA=1 FREHG_CUDA_ARCH=AMPERE80 scripts/ci_install_deps.sh $HOME/frehg-deps-cuda
cmake -B build-cuda -DCMAKE_CXX_COMPILER=$HOME/frehg-deps-cuda/bin/nvcc_wrapper \
      -DCMAKE_PREFIX_PATH=$HOME/frehg-deps-cuda
cmake --build build-cuda -j
```

Facts: Kokkos CUDA builds go through its `nvcc_wrapper` (installed into
the prefix); `FREHG_CUDA_ARCH` is a Kokkos arch name (`AMPERE80`,
`HOPPER90`, `VOLTA70`, …) — match the user's GPU; nvcc constrains the
host gcc version (e.g. CUDA 12.0 needs gcc ≤ 12), so pick the host
compiler from `nvcc --version`; no GPU is needed to *compile*. Set
`NVCC_WRAPPER_DEFAULT_COMPILER` to the chosen host g++. Anyone
validating GPU *execution* should run the benchmark gates and compare
against the documented CPU records before trusting results.

## 8. Runtime environment for the scripts you generate

- `OMP_NUM_THREADS=1` for small/benchmark-scale grids (launch-latency
  bound); threads only pay off on large per-rank subdomains.
- MPI: `mpirun -np N build/src/frehg case.yaml`; decomposition is
  automatic (or pinned via `domain.decomposition`).
- `FI_PROVIDER=tcp` if MPICH+libfabric hangs in `MPI_Finalize`
  (post-run, results unaffected — seen on macOS).
- Extra arguments after the YAML path go to PETSc (prefixes `fs_` for
  the free-surface solver, `gw_` for the subsurface).

## 9. Troubleshooting map

| Symptom | Cause → fix |
|---|---|
| configure: "requires a parallel (MPI) HDF5" | serial HDF5 found → install the MPI flavor matching your MPI (invariant 1) |
| CMake can't find Kokkos/PETSc/HDF5 | wrong/empty prefix → check `lib/cmake/Kokkos/KokkosConfig.cmake` exists under `CMAKE_PREFIX_PATH` |
| `-Werror` failures inside Kokkos headers | CXX is an MPI wrapper (invariant 2) |
| `-Werror` failure in Frehg2 source on an untested compiler | `-DFREHG_WERROR=OFF`, report the warning |
| link errors mixing libstdc++/libc++ | deps and Frehg2 built by different compiler families (invariant 3) |
| MPI/HDF5 errors at output time | stack built against mixed MPIs → rebuild deps against one MPI |
| tiny run is very slow | default thread count → `OMP_NUM_THREADS=1` |
| hang at exit, 100 % CPU | `FI_PROVIDER=tcp` |
| regression tests fail/missing | goldens not distributed (§6) — expected |

## 10. Reference implementations

The CI workflows are complete, working install-and-build scripts to
crib from: `.github/workflows/build.yml` (Ubuntu gcc + clang),
`openmp.yml` (threaded run lane), `sanitize.yml` (ASan/UBSan),
`cuda-compile.yml` (CUDA compile-only). Human-facing detail:
`docs/user-guide/installation.md` and `docs/user-guide/troubleshooting.md`.

## 11. Remote HPC workflows (ssh/scp)

Facts that change the scripts you generate on a cluster:

- **No sudo.** Skip every `apt`/`brew` line. Compiler, CMake, MPI, and
  (usually) parallel HDF5 come from the module system — inspect
  `module avail cmake gcc openmpi mpich hdf5 cuda` and load an HDF5
  variant built against the loaded MPI (often named `hdf5/<ver>-mpi` or
  `hdf5-parallel`). **yaml-cpp is rarely a module**: build it from
  source into the same user prefix (plain CMake,
  `-DYAML_BUILD_SHARED_LIBS=ON`, ~1 minute, no root).
- `scripts/ci_install_deps.sh $HOME/frehg-deps` needs **no root** but
  does need **internet** — run it on the **login node** (compute nodes
  often have none). Building Frehg2 itself is also normally done on the
  login node; configure with `-DFREHG_ENABLE_TESTS=OFF` if the login
  node blocks the GoogleTest download, or pre-fetch on a connected
  machine.
- **Runs go through the scheduler** (SLURM/PBS/LSF). Generate a job
  script alongside the build script, and make it `module load` the
  *exact same modules* the build used (invariants 1–3 apply per-session:
  a job that loads a different MPI than the build will fail
  confusingly). A minimal SLURM template:

  ```bash
  #!/bin/bash
  #SBATCH -N 1 -n 8 -t 02:00:00
  module load gcc/13 openmpi/4.1 hdf5/1.14-mpi cmake     # match the build!
  export OMP_NUM_THREADS=1
  cd $SLURM_SUBMIT_DIR/mycase
  srun ../build/src/frehg mycase.yaml
  ```

- Verification on the cluster: `ctest -L unit` runs on the login node;
  the `mpi`-label tests and real cases belong in a (short) job. GPU
  nodes: the CUDA build compiles on the login node without a GPU (§7),
  but remember the GPU-unvalidated caveat before running physics there.

Three ways to drive all this with an LLM:

1. **Coding agent on your workstation, driving the cluster over ssh**
   (works with any cluster, no software installed remotely): clone the
   repository locally so the agent can read the guides, set up
   **passwordless ssh** (keys + `ControlMaster` — the agent cannot
   answer password prompts), then ask e.g.:

   > Read AGENTS.md and docs/agents/build.md. The cluster is reachable
   > as `ssh hpc1`. Inspect it (`module avail`, compilers, scheduler),
   > then write, upload (scp), and execute scripts that install the
   > dependencies into ~/frehg-deps, build frehg2, run the unit tests,
   > and submit a one-node verification job. Show me the job output.

   The agent gathers facts with `ssh hpc1 '...'`, iterates the scripts
   against real error output, and `scp`s case files up / results back.
2. **Coding agent on the login node** — if the center permits it and
   login nodes have outbound HTTPS: clone the repository there, start
   the agent in the repository root, and work as if local. Respect
   login-node etiquette (build with a modest `-j`; never run
   simulations outside the scheduler).
3. **Chat-only LLM + scp loop** — no agent anywhere: on the cluster run
   an information-gathering one-liner and paste its output into the
   chat together with this file:

   ```bash
   (hostname; cat /etc/os-release | head -2; module avail 2>&1 | head -50; \
    which sbatch qsub; gcc --version | head -1; nvcc --version 2>/dev/null | tail -1) | tee sysinfo.txt
   ```

   The LLM returns an install/build script and a job script; `scp` them
   up, run, and paste any failure output back for the next iteration.
   The verification recipe (§5) is the loop's exit condition.
