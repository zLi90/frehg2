# Report — Stage-3 phase T3 (periodic lateral topology)

Date: 2026-10-28. Companion: `dod-T3.md`. Plan §7/§13.B; owner decision O3 (axis-generic, both axes gated). Amendment S3-A6 records the first-green adjudications and the coupled-lane redesign.

## 1. What landed

| Piece | Where |
|---|---|
| `domain.periodic: {x, y}` schema + cross-field rule (`nx >= 3`), resolved-config echo, run-record `topology:` line | `Config.*`, `ConfigSchema.cpp`, `RunRecord.cpp` (W3.1) |
| Wrapped neighbour ranks (`rankWest_ = pi_ > 0 ? rank_ - 1 : (periodicX_ ? rank_ + (px_ - 1) : MPI_PROC_NULL)` and the three mirrors) — the halo, gid-exchange and matrix legs follow with **no change** to `HaloExchanger` (paired tags 200 + dir / 200 + (dir ^ 1) disambiguate self- and 2-rank wraps) | `core/Grid.cpp` (W3.2) |
| `Grid::EdgeTopology` (KOKKOS_INLINE predicates `west/east/south/north/outsideX/outsideY`, all false on a periodic axis) replacing every raw global-index edge test | `core/Grid.hpp` (W3.2) |
| The 58-site conversion: Predictor (18 predicate sites), Corrector (8), FreeSurface (8, incl. the closed-edge leg fold and the boundary audit), Baroclinic (4), RichardsSolver (4), WetDry (4, stage reconstruction behind topo), SurfaceTransport (4), SubsurfaceTransport (4), TerrainMetric (2, the W3.6 wrapped seam-face slope), SurfaceSolver (1) | the §7 site table (W3.3) |
| Periodic-axis BC rule: side faces never attach on a periodic axis; all-skipped ⇒ fatal with the §13.B message; per-side + skipped counts in the run-record BC summary | `bc/BoundarySet.cpp` (W3.4) |
| Seam diagnostic (bed/soil/Manning/IC rasters): seam jump vs largest interior jump, warning + run-record `topology.warnings` | `Simulation.cpp` (W3.5) |
| Dihedral battery flag swap (axis-swapping transforms map periodic-x onto periodic-y) + constant wind vectors rotated with the domain (the coord_op's linear part) | `run_regression.py` (W3.7) |
| Forbidden rule "raw global-edge test (use EdgeTopology)" — zero hits on the converted tree | `scripts/check_forbidden.sh` (W3.8) |
| Gates pd1–pd8 + pd-restart, `gate_periodic.py` (discrete wind–friction balance, roll/rank/ledger/linear-head/zero-flux/restart/exchange checks), `stage3-pd.yaml`, `unit.pd_gates.negative`; `mpi.halo.n{1,2,4}` gain the periodic far-edge halo case | `benchmarks/pd-periodic/`, `tests/regression/`, `scripts/test_pd_gates.py` |

## 2. Conversion neutrality (the default-bitwise spine)

The conversion is pure refactoring for non-periodic domains: `Ext`-free
predicates fold to the v2.0 comparisons. Proof by execution: 10/10 rl/ob/
bitwise gates green across the T2→T3 binary boundary inside the T2 close
pipeline, plus `default_bitwise` PASS on the converted tree (T0 hashes,
§2.6 declared exceptions only).

## 3. The S3-A6 adjudication (what is a seam defect vs the scheme's own noise)

The branch-free lanes hold strict bounds: pd2-gw roll 1e-13, pd3-gw rank
5e-13, pd-restart flow fields **1.9e-15 absolute**, pd5 qx at the denormal
class, a terrain-following pd2-gw control at 6e-15. The surface lanes
deviate at their measured branch-amplification classes (superbee/wet-dry/
wind flips seeded by 1-ulp solver summation-order changes; exponential
growth uu 3.2e-10 @ 60 s → 7.3e-8 @ 600 s, spread over all columns, NOT
seam-localized), and the **walled control** — the same case with periodic
removed and x walls added — measures strictly LARGER deviations on every
field (restart flow 1.3e-13 vs 1.9e-15; rank uu 6.6e-7 vs 7.0e-8). The
periodic machinery introduces nothing; the bounds now carry the measured
field classes with ≥ 14× headroom while defect-class seam errors fail by
2–8 orders (battery-proven).

## 4. The coupled lane (two failed designs, one diagnostic tool)

The **1e-12-twin control** (perturb one bed cell by 1e-12, rerun, compare)
separates trajectory sensitivity from seam defects. It rejected two
pd2-coupled designs: a marginally-drying row over a vg_n = 1.23 soil
(chaotic, Lyapunov ~40 s — the twin reproduced the roll plateau exactly:
eta 3e-4, uu 1e-2, head 2.0 m) and a vadose gap under a wet surface (a
singular wetting-front event retains O(1e-2) subsurface-tracer memory of
any seed). The adopted saturated tidal basin (Ss = 1e-3 loading, external
closure, consistent salinity both ways, terrain-following) measures roll
deltas uu 1.5e-7 / head 2.1e-11 / c_gw 2.0e-10 with a real exchange signal
(max |Δc_gw| = 0.284 psu, gated at the 0.05 anti-vacuity floor), and its
rank lanes (1/2/4 and 2 × 2) hold the same classes. The twin-control
method is recorded here for reuse: **author identity gates only on
dynamically stable trajectories, and prove stability with the twin before
adjudicating any bound.**

## 5. pd8 (axis generality, O3)

Two batteries through the 8-transform machinery with the W3.7 flag swap:
the rl5 vadose gw base (1e-14 class vs the 1e-12 strict bounds) and the
new surface base (periodic torus bump, oblique wind rotated with the
domain, tracer slab across the seam). The surface battery compares the
scalar fields (eta/depth/tracer — the swe-orient convention; uu/vv are
staggered vector components the scalar push cannot express): transpose at
the pd2 field classes (measured eta 6.5e-11, tracer 3.6e-7), the six
edge-swapping transforms at the x-orient swe class (measured ≤ 1.2e-4
vs 1e-2).

## 6. Tier P / Tier N (S3-A1 scope)

Tier P: regression at `-j8` minus `perf_baseline` + `perf_baseline` solo +
unit/mpi/stage3 — see `records/T3/timings.md`. Tier N (S3-A1 — T3 touched
surface, groundwater and transport): b5 24 h record (the record running
inside the T2 close pipeline executes on the T3-converted binary, so one
record serves both phases — noted in both timings files), b6 ss + td,
g5 × 2, g8. Details and durations: `stage3/records/T3/timings.md`.

## 7. Deviations

S3-A6 (per-field-class bounds for pd2/pd3/pd-restart, the coupled-lane
redesign, the pd8 surface scalar convention — each entry carries its
measured datum). The harness `read_raster2d_file` now accepts
whitespace-separated multi-value lines (the documented input-data.md
format). No v2.0 tolerance changed; no zero_gradient or legacy-gate
behaviour touched.
