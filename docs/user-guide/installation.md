# Installation

This page is the **manual** route: you run the commands yourself.
Alternatively, an AI assistant can write and run tailored
install/build scripts for your system from the condensed
[build guide for AI assistants](../agents/build.md) — the README's
*"Getting started — Path 2"* section shows exactly what to hand the LLM
and example prompts. Both routes end at the same verified build.

## Dependencies

| Dependency | Version | Notes |
|---|---|---|
| CMake | ≥ 3.23 | build system |
| C++20 compiler | gcc ≥ 12 or clang ≥ 15 | `-std=c++20` |
| MPI | any MPI-3 (MPICH, Open MPI) | must provide the C++ bindings |
| Kokkos | ≥ 5.1 | Serial + OpenMP host backends (or a GPU backend) |
| PETSc | ≥ 3.20 | KSP linear solver; built `--with-fc=0` is fine |
| HDF5 | ≥ 1.10, **parallel** | must be built `--enable-parallel` against your MPI |
| yaml-cpp | ≥ 0.7 | configuration parsing |

All dependencies must be built against the **same MPI**. A parallel HDF5 and
a serial HDF5 will both `find_package`, but only the parallel one links; if
you get MPI errors at output time, that is the usual cause.

Frehg2 does **not** bundle or auto-download its dependencies. You install them
once into a **dependency prefix** — a directory you create solely to hold the
built libraries — and then point the build at that prefix. The prefix is a
separate directory from the source tree and from the `build/` directory;
pointing `CMAKE_PREFIX_PATH` at your `build/` folder will *not* work (that is
where Frehg2 compiles *to*, not where the libraries live). Throughout this
guide the prefix is written as `$HOME/frehg-deps` — substitute any path you
like, but use the *same* one for install and for the build.

## Option A — install the heavy dependencies from source (recommended)

Kokkos and PETSc are the two that are simplest to pin by building them into a
single prefix. The repository ships a script that does exactly this. **Run it
first, before configuring the build**, and give it your dependency prefix as
the argument:

```bash
scripts/ci_install_deps.sh $HOME/frehg-deps
```

It downloads and installs Kokkos (OpenMP + Serial backends) and PETSc
(MPI + bundled BLAS/LAPACK, no Fortran) into `$HOME/frehg-deps`. It expects
MPI, parallel HDF5, and yaml-cpp to come from your system package manager
(install those first):

```bash
# Debian/Ubuntu
sudo apt install cmake g++ mpich libhdf5-mpich-dev libyaml-cpp-dev doxygen

# macOS (Homebrew) — note Homebrew's hdf5-mpi provides parallel HDF5
brew install cmake open-mpi hdf5-mpi yaml-cpp doxygen
```

The script pins tested versions (Kokkos 5.1.1, PETSc 3.25.1). Override them
with environment variables if you need to:

```bash
KOKKOS_VERSION=5.1.1 PETSC_VERSION=3.25.1 scripts/ci_install_deps.sh $HOME/frehg-deps
```

When it finishes it prints `dependencies installed at $HOME/frehg-deps` and
drops a `.complete` marker so re-running it is a fast no-op.

## Option B — use system packages for everything

If your distribution ships Kokkos ≥ 5.1 and PETSc ≥ 3.20 (recent Spack,
Fedora, or an HPC module system usually do), install them and skip the
script. Load the modules / activate the Spack environment so the packages
are on `CMAKE_PREFIX_PATH`, then in the build step below you can drop the
`-DCMAKE_PREFIX_PATH` argument entirely.

## A note on compiler wrappers

Do **not** set `CMAKE_CXX_COMPILER=mpicxx`. Kokkos and PETSc that were built
with a plain compiler (e.g. `g++-15`) must be consumed with that same
compiler; point CMake at MPI through `find_package(MPI)` (automatic) rather
than the wrapper. Set `CXX`/`CC` to the real compiler and let CMake add the
MPI flags. For example, with gcc-built dependencies installed under
`$HOME/frehg-deps`:

```bash
env CXX=g++-15 CC=gcc-15 \
    PATH=$HOME/frehg-deps/bin:$PATH \
    CMAKE_PREFIX_PATH=$HOME/frehg-deps \
    cmake -B build -S .
```

## Building Frehg2

With the dependencies installed, configure and build, pointing
`CMAKE_PREFIX_PATH` at the **same dependency prefix** you installed into —
not the `build/` directory:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=$HOME/frehg-deps
cmake --build build -j
```

(If you used Option B and your packages are already on the default search
path, omit `-DCMAKE_PREFIX_PATH`.) The executable is written to
**`build/src/frehg`**.

If `find_package` still reports it cannot find Kokkos, PETSc, or HDF5, the
prefix is wrong or empty — confirm the install step ran and that
`$HOME/frehg-deps/lib/cmake/Kokkos/KokkosConfig.cmake` exists.

Useful CMake options (all default sensibly):

| Option | Default | Effect |
|---|---|---|
| `CMAKE_BUILD_TYPE` | `RelWithDebInfo` | use `Release` for production runs |
| `FREHG_WERROR` | `ON` | warnings are errors (`-Wall -Wextra -Wpedantic -Wshadow -Wconversion`) |
| `FREHG_ENABLE_TESTS` | `ON` | build the GoogleTest / MPI / regression suites |
| `FREHG_SANITIZE` | `OFF` | Address + UndefinedBehavior sanitizers (debugging) |
| `FREHG_GPU_AWARE_MPI` | `OFF` | compile-time default for GPU-aware MPI staging (runtime-overridable) |
| `FREHG_BACKEND` | `auto` | requested execution backend (`serial`/`openmp`/`cuda`/`hip`), verified against the found Kokkos at configure time (v2 plan §2B.2) |
| `CMAKE_PREFIX_PATH` | — | where to find Kokkos/PETSc/HDF5/yaml-cpp |

For a fast production binary:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/frehg-deps
cmake --build build -j
```

## Choosing a parallel lane (v2 plan §2B.2)

One source tree covers every deployment; the execution backend is a *build*
choice (which Kokkos you link) and the linear-algebra backend follows from it:

| Lane | Kokkos backend | `solver.*.mat_type` | Launch | Notes |
|---|---|---|---|---|
| laptop serial | Serial | `aij` | 1 rank | the debug/validation lane; the default |
| laptop OpenMP | OpenMP | `aijkokkos` | 1 rank × N threads | needs ≥ ~1e5 cells to pay off |
| laptop/HPC MPI | Serial or OpenMP | `aij` | N ranks | best CPU lane at benchmark scale (measured) |
| HPC hybrid | OpenMP | `aijkokkos` | ranks × threads, ~1 rank/NUMA domain | the p3 gate records the winning ratio |
| single GPU | CUDA/HIP | forced `aijkokkos` | 1 rank | ≥ ~1e6 cells to pay off |
| multi-GPU | CUDA/HIP | forced `aijkokkos` | 1 rank per GPU + `runtime.gpu_aware_mpi` | see `developer-guide/gpu-acceptance.md` |

`aijkokkos` requires a Kokkos-enabled PETSc (`build_frehg2_local.sh` and
`scripts/ci_install_deps.sh` build one: `--with-kokkos-dir
--download-kokkos-kernels --with-openmp`). On device builds the Kokkos types
are forced and a configure-time check rejects a PETSc without Kokkos support.

## GPU builds (experimental — compile-verified, not device-executed by CI)

All physics kernels are Kokkos and contain no backend `#ifdef`s; building
against a CUDA-enabled Kokkos and a Kokkos+CUDA PETSc
(`scripts/ci_install_deps.sh` with `FREHG_CUDA=1`, compiled through
`nvcc_wrapper`, `-DFREHG_BACKEND=cuda`) compiles **and links, tests
included**, warning-free — enforced by the `cuda-compile` CI lane (gate p4).
The memory-space discipline the CPU lanes cannot see is covered statically
(gate p5: forbidden-pattern scan, compile-time backend invariant with a
negative test, debug-build memtype assertions). GPU lanes carry the
**experimental** status until a passing p6 acceptance bundle
(`scripts/gpu_acceptance.sh`, run on a GPU machine; criteria in
`developer-guide/gpu-acceptance.md`) is returned — v2 plan §2B.4.

## Verifying the build

Run the unit and MPI test labels — they should all pass:

```bash
ctest --test-dir build -L 'unit|mpi' --output-on-failure
```

The single-command per-PR gate (strict build, unit + MPI + fast
regressions, forbidden-pattern scan, parameter-docs lockstep, Doxygen)
is:

```bash
scripts/ci_build_and_test.sh
```

On macOS with gcc-built dependencies, prefix it with the compiler environment from the
note above. The fast regression gates (b1–b4, the coupled b5
restart-determinism and rank-invariance checks) additionally need the
legacy goldens and digitized reference curves; they are pointed at via the
CMake cache variable `FREHG_LEGACY_BENCHMARKS` and are never committed.
The long b5/b6 gate runs carry the separate `regression_nightly` ctest
label and are excluded from the per-PR set — run them deliberately, on
capable hardware, with `ctest --test-dir build -L regression_nightly`.

Smoke-test the executable itself against a shipped case:

```bash
OMP_NUM_THREADS=1 build/src/frehg --validate benchmarks/b1-sw/b1-sw.yaml
```

`--validate` runs the full schema and cross-field check and prints the
effective configuration without running the solver. (Why
`OMP_NUM_THREADS=1`: see
[Running a simulation](running.md#threads-and-performance).)
