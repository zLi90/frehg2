# Benchmark walkthroughs

The six benchmarks are the model's validation gates: each is a published or
analytically known problem, each has a quantitative pass criterion from the
[upgrade plan §9](../developer-guide/FREHG2_UPGRADE_PLAN.md#9-benchmark-validation-gates-and-tolerances),
and a release is only made with all six green. Every case ships as a
runnable configuration under `benchmarks/` with its own README documenting
data provenance and configuration decisions.

Beyond the gates, the repository also ships two **extended validation
suites** under `validation/` — 13 ports of the published SERGHEI benchmark
suite and 6 classic solute-transport benchmarks (surface tracer advection, Ogata–Banks, Henry,
Goswami–Clement, Elder, Kuan) — each case self-contained with its measured
result and run command (`validation/README.md`). They are evidence and
worked examples, not CI gates.

Reference data (legacy ASCII goldens, digitized published curves) is
**never committed**; the regression harness reads it from a sibling legacy
checkout pointed at by the CMake cache variable `FREHG_LEGACY_BENCHMARKS`.

## Running a gate

Each gate is a ctest entry that stages the case into a scratch directory,
runs it, and applies the plan §9 criteria, printing achieved-vs-allowed for
every metric:

```bash
ctest --test-dir build -R 'regression.b1$' --output-on-failure   # any of b1..b4
ctest --test-dir build -L regression_nightly                     # the long b5/b6 runs
```

You can also run any benchmark directly and analyze the HDF5 output
yourself (see the [output reference](output.md)):

```bash
cd benchmarks/b1-sw && OMP_NUM_THREADS=1 /path/to/build/src/frehg b1-sw.yaml
```

## b1 — tilted-plane rainfall runoff (surface water)

The Maxwell et al. (2014) intercomparison tilted plane: a 1×10 column of
80 m cells with terrain-following bed, Manning n = 0.019, and a rain pulse.
Gate: element-wise agreement with the legacy Frehg outputs at 11 output
times (`eta`, `depth` within max(5e-4 m, 1 %); velocities within
max(1e-3 m/s, 2 %)) plus cumulative mass balance ≤ 0.1 % of rainfall.

*Record:* worst achieved/allowed ratio 0.231 (depth); mass-balance error
2.7e-11 of rainfall. Plot `depth` at the outlet monitor to see the
rising/falling hydrograph limb against the digitized ParFlow curve in the
legacy reference directory.

## b2 — vertical infiltration column vs Warrick (groundwater)

A 1×1×100 column (dz = 1 cm) infiltrating toward a Warrick (1971)
analytical solution; exercises the PCA predictor–corrector, the adaptive
subsurface step, and the `aev` saturation cutoff. Gate: element-wise
agreement with the legacy outputs (head within max(1e-3 m, 1 %), θ within
0.005) and the wetting-front depth within 5 % of Warrick at each output
time (and ≤ 1.2× the legacy code's own error).

*Record:* worst field ratio 0.40 (head); wetting front within 2.0 % at
11700 s. Plot the θ(z) profiles at the output times against
`warrick_water_content_profile.csv`.

## b3 — layered-soil infiltration (Kirkland)

The Kirkland et al. (1992) 2D infiltration problem: 50×1×30 cells, two soil
types, a fixed-flux strip on top. Gate: every digitized h = 0 and h = −400
pressure contour point within 0.2 m (2 cells) of the simulated contour,
RMS ≤ 0.1 m, internal mass balance ≤ 0.5 %.

*Record:* h = 0 contour max 0.160 m / RMS 0.073 m; h = −400 max 0.170 m /
RMS 0.088 m; mass balance 0.3 %. Contour `hydraulic_head` at the final time
over the digitized points (the legacy `makeplot.py` reads Frehg2's HDF5
directly).

## b4 — overland-flow hydrograph (Govindaraju / Chézy)

A 200×10 sloping plane against a kinematic-wave reference, run with the
**Chézy** friction law (C = 1.767) added by plan §5.8. Gate: outlet
hydrograph rel-L2 ≤ 0.15, peak within ±10 %, time-to-peak within ±10 %,
outflow volume within 5 % of rainfall and 3 % of the reference.

*Record:* rel-L2 0.082, peak +0.3 %, time-to-99 %-peak −7.9 %, volume
−0.10 % of rainfall / −0.22 % of reference. Plot outlet discharge from the
mass-audit table (the conserved source of boundary outflow) against
`ReferenceData/outflow.txt`.

## b5 — tilted-V catchment (coupled)

The Kollet et al. (2017) tilted-V intercomparison: 101×55×25 cells, rain
and no-rain scenarios, gated against the envelope of the published
ParFlow/CATHY/HGS/Cast3M curves (discharge and ponding storage in the
envelope band, peaks within ±15 % of the envelope mean, recession volume
within 10 %). Both coupling modes are defined in the suite; the four
full-horizon runs (120 sim-hours; wall-hours each) carry the
`regression_nightly` label.

*Record (amendments A15/A17):* the gated record is rain/sync on the
shortened 24 h horizon — discharge 0/48 reference times outside the
envelope, peak +2.6 %, integral +2.4 %; ponding 0/54 outside, peak −1.5 %,
integral −2.5 %; budgets closed with the legacy below-bed clamp measured at
1.6 % of rain. The rain/subcycled lane at the committed `time.dt = 5 s` is
a measured, diagnosed FAIL (fixed-step CFL excursions; the config README
records the ~2 s recommendation). The case README documents every
configuration decision, including the free-outfall outlet strip and the
envelope-reconstructed roughness.

## b6 — Kuan tidal saltwater intrusion (coupled + transport)

The Kuan et al. laboratory sandbox: a 1×68(×18) tank with tidal stage,
salt transport, and baroclinic density feedback, in steady (`ss`) and tidal
(`td`) variants. Primary gate (experiment): 50 %-isohaline interface
position MAE ≤ 15 % of tank height and ≤ 1.5× the legacy golden's own
error; tidally averaged salt mass within 10 % of the golden; salinity
bounded in [0, 35] everywhere. Secondary gate: loose element-wise
agreement with the legacy goldens (which used the dropped Newton scheme).

*Record:* interface MAE 0.0334 m (ss) and 0.0389 m (td) — the td model
agrees with the experiment *better than the legacy golden itself*
(0.0430 m); salt mass +5.5 %/+8.8 %; salinity exactly [0, 35] emergent.
The td configuration prescribes the golden's sea-surface salinity as an
explicit `scalar_value` condition (amendment A20; the case README has the
provenance). Plot the 50 %-isohaline over the experimental points from
`/transport/concentration` at the tidally averaged state.
