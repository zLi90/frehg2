# b1-sw — tilted-plane rainfall runoff (Maxwell et al. 2014, MAXP1)

Rainfall runoff on a tilted plane: a 1 × 10 strip of 80 m cells
(`follow_terrain: true`, bed from `input/bathymetry.dat`), Manning
n = 0.019, a rain-pulse series, dt = 5 s to t = 18000 s (3600 steps).
Surface-water module only, with no configured boundary conditions: the
domain drains through the *preserved legacy edge behavior* (boundary
faces see the pre-source ghost stage, so rain induces a small outward
boundary-face flow each step — report-P1 §"fidelity findings" item 3),
which the mass audit counts as boundary outflow. The case is the MAXP1
tilted-plane problem of the Maxwell et al. (2014) integrated hydrologic
model intercomparison, and it is the P1 gate's primary legacy-fidelity
benchmark: the same algorithm as legacy Frehg, so the comparison is
element-wise and tight.

Two configuration details preserve legacy behavior the goldens embed:

1. **`rainfall.exclude` covers the outlet row j = 9** (plan amendment
   A3): legacy `evaprain` hard-coded "no rain on the last global row"
   (`shallowwater.c:596`), a b1-specific rule the golden outputs contain.
   Cases without the quirk simply omit the key.
2. The monitor at (i = 0, j = 4) reproduces the legacy monitor location,
   giving the per-step hydrograph the intercomparison plots use.

`input.v1.yaml` is the archival v1 draft; the committed `b1-sw.yaml` is
its `tools/migrate_yaml_v1_to_v2.py` output augmented from the legacy
input file (provenance in the YAML header), and the migration test
asserts that embedding on every CI run.

Gate (plan §9 b1, run by `tests/regression/run_regression.py b1`):
element-wise against the legacy ASCII goldens at 11 output times —
`eta`/`depth` within max(5e-4 m, 1 %), `uu`/`vv` within max(1e-3 m/s,
2 %) — plus cumulative mass-balance error ≤ 0.1 % of total rainfall.
The goldens live in the development-side `legacy/benchmarks/` archive
(not distributed with this repository), which also carries the digitized
ParFlow hydrograph for visual comparison.

Gate results at P1 (allowed in parentheses): worst achieved/allowed
field ratios — depth 0.231, vv 0.131, eta 0.063 (1.0); cumulative
mass-balance error 2.7e-11 of rainfall volume (1e-3). Runtime: well
under a second serial (`OMP_NUM_THREADS=1`).
