# Frehg2

Frehg2 is a production-grade rewrite of Frehg 1.0 — a coupled semi-implicit
2D shallow-water / 3D mixed-form Richards / solute-transport model — in
C++20 with Kokkos (CPU/GPU-portable), MPI, PETSc, yaml-cpp configuration,
and single-file parallel HDF5 output. The core physics is preserved exactly
from the validated legacy code; each of the six benchmarks is a blocking,
quantitative validation gate. The binding specification is
[`docs/developer-guide/FREHG2_UPGRADE_PLAN.md`](docs/developer-guide/FREHG2_UPGRADE_PLAN.md).

**What Frehg2 simulates**

- **surface water** — the θ-scheme 2D shallow-water solver: rainfall and
  evaporation, wind stress, Manning or Chézy friction, wetting/drying, and
  an implicit free surface (PETSc CG);
- **groundwater** — the mass-conservative PCA mixed-form Richards solver
  (Li et al. 2021): van Genuchten soils, adaptive subsurface stepping,
  head/flux/free-drainage/hydrostatic boundary conditions,
  terrain-following or partial-cell meshes;
- **coupled surface–subsurface flow** — both modules exchanging
  infiltration and seepage each step (wet cells drive the subsurface as a
  ponded-head boundary limited to the water actually available; seepage
  returns to the surface), in lockstep (`sync`) or subcycled mode;
- **scalar transport** — one scalar (e.g. salinity) advected with upwind or
  TVD-superbee fluxes, diffused on the surface, dispersed anisotropically
  in the subsurface, exchanged through the seepage, and optionally fed back
  into the subsurface density/viscosity (baroclinic saltwater intrusion).

Every mode has HDF5 field output, mass-audit tables (including the scalar
budget), point monitors, and checkpoint/restart (restart reproduces the
uninterrupted run bitwise). Frehg2 is **GPU-ready but GPU-unvalidated**:
all kernels are Kokkos with no backend `#ifdef`s and the CUDA backend is
compile-gated in CI, but the released version has been validated on CPU
backends only.

## Getting started — two paths

Frehg2 is designed to be set up either **by hand** (Path 1) or **through
an AI assistant** (Path 2): the repository ships condensed context files
that let an LLM install, build, configure, and post-process Frehg2 for
you, on your specific system, without reading the source code. Both
paths end at the same place — a built `build/src/frehg` verified by the
test suite — so pick whichever suits you and switch freely.

### Path 1 — build it yourself

Install the dependencies (MPI, parallel HDF5, yaml-cpp from your package
manager; Kokkos + PETSc via the pinned script), then build:

```bash
scripts/ci_install_deps.sh $HOME/frehg-deps
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$HOME/frehg-deps
cmake --build build -j
```

Validate and run a shipped case:

```bash
build/src/frehg --validate benchmarks/b1-sw/b1-sw.yaml
cd benchmarks/b1-sw && OMP_NUM_THREADS=1 ../../build/src/frehg b1-sw.yaml
```

Run the test suite:

```bash
ctest --test-dir build -L 'unit|mpi' --output-on-failure
```

Details, options, per-platform notes, and troubleshooting: the
[installation guide](docs/user-guide/installation.md).

### Path 2 — let an AI assistant do it

Three task-scoped guide files carry everything an LLM needs — the
dependency matrix, the build invariants and their failure modes, the
configuration conventions and pitfalls, and the output-file contract —
so you never pay for the model to read the code base:

| Task | Give the LLM | It produces |
|---|---|---|
| install + build + verify | [AGENTS.md](AGENTS.md) + [docs/agents/build.md](docs/agents/build.md) | install/build scripts tailored to your OS, package manager, compiler, and (optionally) GPU |
| set up a simulation case | [docs/agents/case-setup.md](docs/agents/case-setup.md) **+** [docs/user-guide/parameters.md](docs/user-guide/parameters.md) | a validated YAML configuration and its input data files |
| analyze results / make plots | [docs/agents/postprocessing.md](docs/agents/postprocessing.md) | Python scripts reading the HDF5 output (maps, hydrographs, budget checks) |

**With a coding agent** (Claude Code, Cursor, Copilot Workspace, …):
clone the repository, start the agent in the repository root, and ask in
plain language — agents that follow the `AGENTS.md` convention find the
guides themselves; otherwise start your request with *"read AGENTS.md
first."* For example:

> Read AGENTS.md, inspect this machine, and write and run scripts that
> install the dependencies, build frehg2, and verify the build.

**With a chat-only LLM** (no file access): paste the guide file(s) from
the table above into the conversation together with your request and
your system details. Example prompts:

> Here is frehg2's build guide [paste `docs/agents/build.md`]. I am on
> Ubuntu 24.04 with sudo, gcc 13, no GPU. Write me a single script that
> installs everything, builds frehg2, and runs the verification steps.

> Here are frehg2's case-setup guide and parameter reference [paste
> both]. Create a case: rain of 20 mm/h for 1 h on a 100 m × 50 m plane
> with 2 % slope and a 1 m loam soil column underneath, draining east.

> Here is frehg2's post-processing guide [paste
> `docs/agents/postprocessing.md`]. Write a script that plots the outlet
> hydrograph and checks the mass balance of `out/output.h5`.

Whatever the LLM produces, the model's own gates are your safety net:
configurations must pass `frehg --validate` (strict schema, actionable
errors), builds must pass `ctest -L 'unit|mpi'`, and every shipped
example case's README states its expected result. Insist the assistant
runs those checks — the guides instruct it to.

**On a remote HPC (ssh/scp)?** The same files work three ways — a coding
agent on your workstation driving the cluster over ssh, an agent on the
login node, or a chat-only LLM iterating through an scp loop — including
module-system installs (no sudo) and scheduler job scripts. See
[docs/agents/build.md](docs/agents/build.md) §11.

## Documentation

The full manual lives under `docs/` and builds into a site with
`python3 -m mkdocs build` (mkdocs-material). The API reference is Doxygen
(`doxygen docs/Doxyfile` → `docs/api/html`).

| I want to… | Read |
|---|---|
| install and build | [user-guide/installation.md](docs/user-guide/installation.md) |
| set up a run | [user-guide/configuration.md](docs/user-guide/configuration.md), [user-guide/parameters.md](docs/user-guide/parameters.md) (every key, CI-checked) |
| prepare input files | [user-guide/input-data.md](docs/user-guide/input-data.md) |
| run and analyze | [user-guide/running.md](docs/user-guide/running.md), [user-guide/output.md](docs/user-guide/output.md) |
| checkpoint and restart | [user-guide/restart.md](docs/user-guide/restart.md) |
| reproduce the validation | [user-guide/benchmarks.md](docs/user-guide/benchmarks.md) |
| build a case from scratch | [user-guide/worked-example.md](docs/user-guide/worked-example.md) |
| port a legacy Frehg case | [user-guide/migration.md](docs/user-guide/migration.md) |
| fix a problem | [user-guide/troubleshooting.md](docs/user-guide/troubleshooting.md) |
| **use an AI assistant** to write install scripts, set up cases, or make plots | *Getting started — Path 2* above; the guide files live at [AGENTS.md](AGENTS.md) and [docs/agents/](docs/agents/) |
| understand the numerics | [docs/theory/](docs/theory/) (equations with legacy `file:line` provenance) |
| modify the code | [docs/developer-guide/](docs/developer-guide/) (architecture, grid/indexing, halo contract, tutorials, testing/CI, performance) |

## Benchmarks

Runnable example cases, one directory each under `benchmarks/`; each README
documents data provenance and configuration decisions, and
[user-guide/benchmarks.md](docs/user-guide/benchmarks.md) walks through
running and gating each one.

| Case | Mode | What it is |
|---|---|---|
| `b1-sw` | surface | tilted-plane rainfall runoff (Maxwell et al. 2014; legacy golden gate) |
| `b4-govindaraju` | surface | sloping-plane runoff hydrograph vs a kinematic-wave reference (Chézy friction) |
| `b2-gw` | groundwater | 1D infiltration column vs the Warrick (1971) analytical solution |
| `b3-kirkland` | groundwater | 2D layered-soil variably saturated flow (Kirkland et al. 1992) |
| `b5-vcatchment` | coupled | tilted-V catchment vs the Kollet et al. (2017) intercomparison envelope, rain and no-rain scenarios |
| `b6-kuan` | coupled + transport | Kuan tidal saltwater-intrusion sandbox, steady (`-ss`) and tidal (`-td`) variants vs the legacy goldens and the published experiment |

## Repository map

| Path | Contents |
|---|---|
| `src/core` | configuration, grid/decomposition, halo exchange, PETSc linear systems, time series, timers, logging |
| `src/swe` | the surface-water module (momentum, free surface, wet/dry, sources) |
| `src/gw` | the groundwater module (predictor/corrector, reallocation, terrain mesh, baroclinic ratios) |
| `src/transport` | the scalar-transport module (advection limiters, dispersion, scalar exchange) |
| `src/coupling` | the surface–subsurface coupler and exchange |
| `src/driver` | the time loops, output/monitor/checkpoint mediation, restart |
| `src/bc` | polygon regions and rasterized boundary conditions |
| `src/io` | parallel HDF5 output, monitors, checkpoints, gridded-input readers |
| `benchmarks` | the six benchmark configurations and their input data |
| `validation` | extended validation suites: 13 SERGHEI-benchmark ports and 6 classic transport benchmarks ([validation/README.md](validation/README.md)) |
| `tests` | GoogleTest unit suites, MPI drivers (1/2/4 ranks), regression harness |
| `scripts` | CI gates: build+test, forbidden scan, sanitizer/OpenMP lanes, scaling study, docs lockstep |
| `docs` | the manual (mkdocs), theory reference, developer guide, Doxyfile |

## Status

**v1.0.0** — all six benchmark gates green; phases P0 (foundation),
P1 (surface water), P2 (groundwater), P3 (coupling), P4 (transport +
density coupling), and P5 (hardening, performance, documentation, release)
complete. The per-phase Definitions of Done and handoff reports are under
[`docs/developer-guide/`](docs/developer-guide/).

## License and citation

BSD 3-Clause ([`LICENSE`](LICENSE)). If you use Frehg2 in research, cite
the software and the model papers — see [`CITATION.cff`](CITATION.cff):

- Li, Z., Hodges, B.R. (2019). Model instability and channel connectivity
  for 2D coastal marsh simulations. *Environ Fluid Mech* 19, 1309–1338.
  doi:10.1007/s10652-018-9623-7
- Li, Z., Hodges, B.R. (2019). Modeling subgrid-scale topographic effects
  on shallow marsh hydrodynamics and salinity transport. *Adv. Water
  Resour.* 129, 1–15. doi:10.1016/j.advwatres.2019.05.004
- Li, Z., Özgen-Xian, I., Maina, F.Z. (2021). A mass-conservative
  predictor-corrector solution to the 1D Richards equation with adaptive
  time control. *J. Hydrol.* 592, 125809. doi:10.1016/j.jhydrol.2020.125809
- Li, Z., Hodges, B.R., Shen, X. (2023). Modeling hypersalinity caused by
  evaporation and surface-subsurface exchange in a coastal marsh — the
  Frehg 1.0 model description and application. *J. Hydrol.* 618, 129268.
  doi:10.1016/j.jhydrol.2023.129268
