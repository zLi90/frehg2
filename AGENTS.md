# AGENTS.md — Frehg2 guide for AI assistants

This file is the entry point for an AI assistant (chat LLM or coding
agent) helping a user install, run, or analyze **Frehg2** — a coupled
semi-implicit 2D shallow-water / 3D variably-saturated groundwater /
solute-transport model (C++20, Kokkos, MPI, PETSc, parallel HDF5,
strict-schema YAML configuration). You do **not** need to read the
source code for the tasks below; each task has a condensed guide that
carries the full contract. (Human readers: the README's *"Getting
started — Path 2"* section shows how to hand these files to your LLM,
with example prompts.)

## Task routing

| The user wants to… | Read (paste into the LLM) |
|---|---|
| install dependencies, build, and compile Frehg2 on their system | `docs/agents/build.md` |
| create a new simulation case (YAML config + input files) | `docs/agents/case-setup.md` **plus** `docs/user-guide/parameters.md` (the complete, CI-checked key reference) |
| write Python scripts to read results and make plots | `docs/agents/postprocessing.md` |

For humans, the full manual is the mkdocs site under `docs/` (user
guide, theory reference with legacy provenance, developer guide) and the
top-level `README.md`.

## Universal facts (apply to every task)

- **Units are SI everywhere** (meters, seconds, m/s); elevations are
  absolute (one datum for bed, water surface, and heads).
- **The executable is `build/src/frehg`**; a run is
  `frehg <case>.yaml`, with output paths resolved against the current
  working directory (run from the case directory).
- **Validate before running.** `frehg --validate <case>.yaml` checks the
  full schema (unknown keys rejected *with nearest-key suggestions*),
  types, ranges, cross-field rules, and input-file existence — iterate
  until it prints `VALID:`. This loop is the backbone of reliable
  AI-generated configs.
- **Small grids run fastest single-threaded**: prefix runs with
  `OMP_NUM_THREADS=1` for anything under ~1M cells (kernel-launch
  latency dominates otherwise). On large per-rank subdomains threads pay
  off, and with a Kokkos-aware PETSc the linear *solve* threads too
  (v2 Q3; `solver.<sys>.mat_type: aijkokkos`). Parallel runs:
  `mpirun -np N ...`; results are rank-invariant to round-off.
- If a run **hangs at exit** with ~100 % CPU under MPICH/libfabric, set
  `FI_PROVIDER=tcp` (known platform issue; results are unaffected).
- **GPU status**: all kernels are Kokkos with no backend `#ifdef`s and
  the CUDA lane compiles+links in CI, but the GPU path ships
  **experimental** — CPU-physics-verified, device execution unverified
  pending the owner GPU-acceptance bundle. Never promise validated GPU
  results.
- Ready-made worked examples live in `benchmarks/` (the six validation
  gates) and `validation/` (19 extended cases); every case directory is
  runnable and its README states the expected result. Start new work by
  copying the nearest example, not from a blank file.

## Repository shape (what matters for these tasks)

```
CMakeLists.txt          CMake ≥ 3.23, C++20; options documented in docs/agents/build.md
src/                    the model (you do not need to read it)
scripts/                canonical install/CI scripts (ci_install_deps.sh builds Kokkos+PETSc)
benchmarks/, validation/  runnable example cases with READMEs and provenance
docs/user-guide/        installation, configuration, parameters.md, output reference
docs/agents/            the condensed task guides this file routes to
tests/                  ctest suite; `ctest -L 'unit|mpi'` needs no external data
```

## What to collect from the user's system before writing scripts

Operating system and version; package manager (apt/dnf/brew/spack) or
HPC module system; existing compiler and version (`g++ --version`,
`clang++ --version`); existing MPI (`mpiexec --version`); whether a GPU
and CUDA toolkit are present (`nvcc --version`); core count and RAM;
whether the user has sudo. On an HPC additionally: the scheduler
(SLURM/PBS/LSF), whether login nodes have internet access, and the
relevant `module avail` output — `docs/agents/build.md` §11 covers the
three remote (ssh/scp) workflows. That guide turns these answers into a
concrete install path.
