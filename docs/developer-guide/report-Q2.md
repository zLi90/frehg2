# Q2 phase report — Run provenance and timing log

**Scope:** v2 development plan §2A (Q2). **Closed:** 2026-09-11.
**Companion:** [dod-Q2.md](dod-Q2.md) (verification commands).

## 1. What was delivered

Every simulation now leaves a complete, machine-readable account of itself:
`run-record.yaml` beside the HDF5 output, holding seven sections —
provenance, the resolved configuration, modules, boundary conditions,
the timer tree, solver telemetry, and mass-audit closure (the §2A.2 design,
implemented as specified; format documented in the user guide's output
reference). Key properties:

- **Truthful under failure.** The record is rewritten at every output
  flush (write-then-rename, so a crash cannot truncate the previous copy)
  and finalized after `run()` returns; a killed run leaves
  `finished: false` with the last completed state.
- **The resolved config is the record, not the input.** A new serializer
  (`resolvedConfigYaml`) is the exact inverse of extraction; `frehg
  --resolve <input>` prints it (sentinel-fenced — the logger and PETSc
  share stdout). Gate r1 verifies both directions on every gate run: the
  embedded config revalidates against the schema, and re-resolving the
  original input reproduces it exactly.
- **One source of truth per number.** The timer tree comes from a new
  structured `Timer::merged()` that the stdout table is now rendered
  from; the solver section is built by the same lambda that prints the
  `solver summary` lines; the closure section reads the same rank-0
  accumulators the mass-audit monitors record.
- **Provenance hashes.** The input YAML and every BC series file carry a
  SHA-256 (self-contained FIPS 180-4 implementation, ~90 lines — no new
  dependency); restarted runs record their parent checkpoint and resume
  time (verified in the b1-restart lane's restarted record).

**Timing coverage completed.** New sections: `init`, `io/output`,
`io/checkpoint`, `io/run_record`, `monitors`, `swe/begin_step`, and `halo`
— message passing, timed once at the single funnel (`exchangeSelected`)
every public exchange goes through. With the pre-existing swe/gw/transport/
solve sections, the tree now separates every category the owner asked for:
shallow-water solver, groundwater solver, transport, linear-solve
setup/solve, message passing, initialization, I/O, and monitors.

## 2. Gate record

| Gate | Result |
|---|---|
| r1 (record on every gate run: presence, schema, launch, round-trip) | **PASS across the full per-PR sweep, 40/40** (2026-09-11) — all four regimes, restart (with parentage), 1/2/4 ranks, amg/gamg lanes |
| r1 negative battery (`unit.run_record_checker`) | **PASS** — valid record accepted; all 15 corruptions rejected |
| r2 (`regression.r2_timer_coverage`) | **PASS** — coverage **100 %** at 1 rank (≥ 90 %); halo 0.139 s / io nonzero at 4 ranks; record write **0.01–0.03 %** of the loop (≤ 1 %) |
| Physics untouched | bitwise restart + both rank-invariance modes green in the same sweep |
| Sanitizer lane (clang ASan+UBSan, unit + mpi labels) | clean (`run_sanitizers.sh: clean`) |

## 3. Decisions worth recording

- **`--resolve` output is sentinel-fenced** rather than clean stdout: the
  logger banner and PETSc's unused-options warnings share the stream, and
  suppressing them globally would change unrelated behavior. The checker
  extracts between the fences.
- **The r2 coverage rule is section-aware, not slash-counting:** frame
  names may embed slashes ("swe/free_surface" is one frame), so top-level
  cover is computed against the actual section set (a section counts if no
  other section is its ancestor). Measured cover is 100 % because the
  three loop bodies are now fully tiled (begin_step + free_surface +
  velocity + transport + monitors for the surface loop; groundwater +
  monitors for the gw loop; coupled_step + transport + monitors coupled).
- **`init` covers the pre-loop work inside `run()`** (restart restore or
  first outputs). Module/grid construction happens before `run()` and is
  visible as `provenance.wall_seconds` minus the `simulation` timer — kept
  out of the loop timer deliberately so `simulation` still means "the time
  loop" everywhere (A23 comparability).
- **`omp_threads` reports `OMP_NUM_THREADS` when set, else Kokkos
  host concurrency** — the launch-geometry number the harness can assert
  against; Kokkos's own thread pool is what it reflects when unpinned.
- The forbidden-pattern scan rejected `snprintf` in the SHA-256 hex
  rendering; rewritten with `<iomanip>` (the scan is the arbiter, plan
  §11.3 — no suppression).
- **Sanitizer-lane hang fixed in passing:** the first Q2 sanitizer run
  wedged for 8+ hours inside `MPI_Finalize` (libfabric's sockets provider
  spinning in `sock_cq_sreadfrom` — the known A19 failure). The regression
  ctest entries pin `FI_PROVIDER=tcp`, but `run_sanitizers.sh`'s own
  unit/mpi labels did not; the script now exports the pin for every label.
  The Q1 lane had dodged this by scheduling luck ("intermittently wedges",
  per the A19 record), not by immunity.

## 4. Notes for Q3+ (evaporation onward)

- New physics phases get their timing and provenance **for free**: any new
  `Timer::Scoped` section appears in the record automatically, and r1
  already validates every new gate case's record the moment the case is
  wired into the harness.
- The resolved-config serializer must be extended **in the same PR** as
  any schema change (new keys ⇒ new emit lines), or r1's round-trip fails
  immediately — this is deliberate: the record cannot silently fall behind
  the schema. The failure message names the differing keys (unified diff).
- `run_scaling.py` reads records now; if a future mode needs a metric that
  is not in the record, add it to the record rather than parsing stdout.
- The v1 OpenMP-lane scripts and CI workflows need no change: the record
  is written unconditionally and adds ≤ 0.03 % to the loop.
