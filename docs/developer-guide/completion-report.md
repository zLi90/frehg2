# Frehg2 upgrade — completion report (v1.0.0)

**Project:** rewrite Frehg 1.0 — a validated but unreleasable research
code (~8,300 lines of serial C, LASPack, hand-rolled MPI, ASCII I/O, no
tests) — as Frehg2, a production-grade, publicly releasable C++20 /
Kokkos / MPI / CMake / PETSc / YAML / HDF5 model, preserving the core
physics exactly and gating every step on quantitative benchmark
validation. The binding specification is
[FREHG2_UPGRADE_PLAN.md](FREHG2_UPGRADE_PLAN.md); this report is the
project-level summary of what was delivered.

**Outcome:** released as **v1.0.0** with all six benchmark gates green in
a single verification pipeline, zero-warning strict builds on gcc and
clang, a sanitizer-clean test suite, bitwise-deterministic restart in
every mode, rank-invariant parallel results, complete user / developer /
theory documentation, and two extended validation suites (19 additional
published-reference cases) beyond the gates.

---

## 1. Delivery timeline

| Phase | Scope | Completed | Gates closed |
|---|---|---|---|
| P0 | Foundation: build system, config schema + validation, grid/decomposition, halo exchange, PETSc COO linear systems, polygon BCs, parallel HDF5 + checkpoint, CI scaffolding | 2026-07-18 | — (infrastructure criteria) |
| P1 | Surface-water module (θ-scheme SWE), driver, monitors | 2026-07-25 | b1, b4 |
| P2 | Groundwater module (PCA Richards, van Genuchten, adaptive stepping, reallocation) | 2026-07-27 | b2, b3 |
| P3 | Surface–subsurface coupling (sync + subcycled), coupled audits, coupled restart | 2026-08-11 | b5 |
| P4 | Scalar transport (upwind/superbee, dispersion tensor), baroclinic density coupling | 2026-08-21 | b6 (ss + td) |
| P5 | Hardening (sanitizer matrix, OpenMP/CUDA lanes), performance, documentation, release | 2026-08-22 | release pipeline |

Each phase closed with a verified Definition of Done (`dod-P*.md`) and a
handoff report (`report-P*.md`); all previous phases' criteria were
re-verified at every close. Deviations from the plan were recorded as 23
written amendments (A1–A23) in the plan's amendment log, each with the
measurement that forced it.

## 2. The validated model

Frehg2 simulates 2D semi-implicit shallow-water flow (Manning or Chézy
friction, wetting/drying, rain/evaporation, wind), 3D variably saturated
groundwater flow (the mass-conservative PCA predictor–corrector of Li et
al. 2020 with van Genuchten retention), their two-way coupling
(infiltration/seepage exchange, lockstep or subcycled), and scalar
transport with optional baroclinic density feedback — one executable,
five statically layered libraries, ~14.5k lines of C++20 with all
physics in Kokkos kernels (no backend `#ifdef`s), MPI-decomposed over a
2D rank grid, PETSc CG linear solves, strict-schema YAML configuration,
and single-file parallel HDF5 output with bitwise checkpoint/restart.

### Benchmark gate record (plan §9; achieved vs allowed)

| Gate | Reference | Achieved (allowed) |
|---|---|---|
| **b1** tilted-plane rainfall runoff | legacy goldens, 11 output times | worst field ratio 0.231 of tolerance; mass error 2.7e-11 of rainfall (1e-3) |
| **b2** infiltration column | legacy goldens + Warrick (1971) | worst field ratio 0.40; wetting front within 2.0 % of Warrick (5 %) |
| **b3** layered-soil infiltration | Kirkland et al. (1992) contours | contour max 0.160/0.170 m (0.2), RMS 0.073/0.088 m (0.1); mass 0.3 % (0.5 %) |
| **b4** overland hydrograph (Chézy) | kinematic-wave reference | rel-L2 0.082 (0.15); peak +0.3 % (±10 %); volume −0.10 % of rain (5 %) |
| **b5** tilted-V coupled catchment | Kollet et al. (2017) 4-model envelope | rain/sync: 0/48 discharge and 0/54 ponding points outside the envelope; peaks +2.6 %/−1.5 % (±15 %) — A15/A17 shortened-horizon record |
| **b6** Kuan tidal saltwater intrusion | laboratory experiment + legacy goldens | interface MAE 0.0334 m (ss) / 0.0389 m (td) vs caps 0.107 m — the td model beats the legacy golden's own 0.0430 m; salt mass +5.5/+8.8 % (10 %); salinity bounded [0, 35] emergent |

The final release pipeline (`scripts/ci_release_gate.sh`) re-ran all six
green on the release commit.

### Extended validation (beyond the gates)

Two independently authored suites under `validation/` (see
`validation/README.md` in the repository):

- **the `swe-*`/`re-*`/`swere-*` cases** — 13 ports of the published SERGHEI
  verification suite (SWE, Richards, and coupled SWE-RE *GMD* papers),
  with a feasibility analysis of the full suite, the SERGHEI→Frehg2
  boundary-condition translation table, and measured results/costs per
  case.
- **the `transport-*` cases** — 6 classic solute-transport benchmarks
  (surface tracer advection, Ogata–Banks, Henry, Goswami–Clement, Elder,
  Kuan) forming a complexity ladder to fully coupled variable-density
  transport, with the transport module's three load-bearing conventions
  documented where they shape case design.

## 3. Fidelity: what was preserved, fixed, and dropped

- **Preserved exactly** (plan §3.1, Appendix B provenance; per-module
  equation tables in `docs/theory/`): the θ-scheme with its CFL-damped
  advection and thin-layer drag switch, the PCA predictor–corrector with
  the `aev` cutoff and adaptive stepping, the flux-based coupling, the
  transport schemes with the legacy ledger structure, and every
  scheme-defining numerical constant.
- **Legacy defects resolved, not carried** — each measured, derived, and
  amendment-logged: the seepage unit ambiguity (§5.7; conservation
  proof + per-step audit), the infiltration-limit porosity factor (A9),
  the subcycled exchange accumulator (A9), rain silently skipped on dry
  cells in coupled runs (A11), the uninitialized surface
  exchange-diffusion coefficient (`smap->dz`, A19), rank-dependent
  superbee guards and unexchanged diagonal ghosts (A19), silent
  continuation on solver divergence (now fatal), and the legacy
  below-bed clamp measured into the mass audit instead of hidden (A16).
- **Deliberately dropped** (owner-approved, plan §3.2): the diffusive
  wave, subgrid topography, Newton and modified-Picard groundwater
  schemes, modified van Genuchten, the aerodynamic evaporation model,
  and every dead code path — all recorded with rationale in
  `docs/theory/removed-features.md`.

## 4. Quality infrastructure

- **38 ctest entries**: 146 GoogleTest unit tests, MPI drivers at 1/2/4
  ranks (bitwise halo exchange), 14 per-PR regression gates (b1–b4,
  restart determinism for every mode, rank-invariance lanes), and the
  nightly-class b5/b6 gates; ~5.8k lines of test code.
- **Restart is bitwise** in all four regimes (surface, groundwater,
  coupled, coupled+transport), enforced per-PR.
- **Rank invariance** proven at machine precision in strict solver mode
  and gated at data-derived bounds in production mode (A5/A8/A14).
- **Sanitizers**: ASan+UBSan clean over every per-PR label in full plus
  shortened runs of every nightly-class configuration (A21).
- **OpenMP lane**: the threaded-backend suite including a bitwise
  coupled+transport restart at 4 threads.
- **Static gates**: zero-warning builds (gcc + clang, `-Werror` with
  `-Wconversion`), the forbidden-pattern scan (no TODOs, no warning
  suppression, no raw allocation, no UVM, dropped-feature tripwire), the
  schema↔docs lockstep check, Doxygen with undocumented-public-API as
  error, and `mkdocs build --strict`.
- **CI workflows** for all of the above plus a CUDA compile-only lane —
  authored and pinned, but never executed (the development repository
  had no remote; see §6).

## 5. Performance

Full report: [performance.md](performance.md). Headlines (Apple M3,
4 P + 4 E cores): 79.0 % strong-scaling efficiency at 4 ranks on the
16×-b5 gate case (recorded min-over-repeats artifact; best observed
84.8 %; amendment A23 documents the fanless machine's measurement
variance), near-ideal 2-rank efficiency, and a per-module breakdown
showing the coupled b5 regime is subsurface-bound (~76 % PCA kernels,
~24 % CG solve). Benchmark-scale grids run fastest at
`OMP_NUM_THREADS=1` (kernel-launch bound).

## 6. Known limitations and deferred items

- **GPU-ready, GPU-unvalidated** (A22, owner dispensation): all kernels
  are Kokkos and the CUDA compile-only CI lane is authored, but no GPU
  toolchain existed locally — the lane has never compiled and no gate
  has executed on a device. First order of business on GPU hardware.
- **CI has never run on real runners** (no remote existed); the
  golden-dependent regression steps also need a runner-side legacy
  checkout that no workflow provisions. Watch the first push.
- The legacy benchmark **goldens are not committed** (plan §8.3); the
  regression gates require a local legacy checkout via
  `FREHG_LEGACY_BENCHMARKS`.
- Model-scope limits recorded where they surfaced: subsurface scalar
  injection is y+-only, the density coefficients are compile-time
  constants, dispersion acts on volumetric face fluxes
  (`validation/README.md`); no surface infiltration source
  term, Darcy–Weisbach friction, or periodic BCs
  (`validation/README.md`); masked-column coupled paths and
  transport-under-subcycled-coupling are implemented and unit-covered
  but not benchmark-gated (report-P4/P5).

## 7. Reproducing the record

```bash
scripts/ci_build_and_test.sh          # per-PR gate: build, unit+mpi, b1-b4, docs
scripts/ci_release_gate.sh            # all six benchmark gates in one pipeline
scripts/run_sanitizers.sh --full      # the A21 sanitizer matrix
scripts/run_openmp_lane.sh build 4    # the threaded lane
scripts/run_scaling.py --case synthetic --ranks 1 4 --repeats 3 ...  # A23 scaling gate
```

Per-phase verification commands live in `dod-P0.md` … `dod-P5.md`; the
decisions and their reasons in `report-P0.md` … `report-P5.md` and the
plan's amendment log.
