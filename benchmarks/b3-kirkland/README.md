# b3-kirkland — 2D layered-soil infiltration (Kirkland et al. 1992)

Groundwater-only gate case for phase P2 (plan §9): a 5 m x 3 m vertical
slice with a Berino loamy fine sand trench (x in [1, 4] m, top 1 m and below
2 m depth) embedded in Glendale clay loam, infiltrated through a fixed-flux
strip (5.787e-6 m/s over x in [0.99, 4.01]) for 24 h. Authored from the
SERGHEI-style reference inputs (`legacy/benchmarks/b3-kirkland/input/`);
there is no legacy-Frehd golden for this case.

## Gate

`tests/regression/run_regression.py b3` (ctest `regression.b3`): every
digitized h = 0 and h = -400 contour point within 0.2 m (2 cell widths) of
the simulated contour at t = 86400 s, RMS over all points within 0.1 m, and
internal mass balance within 0.5 % of the injected volume (plan §9). The
digitized points are parsed from the legacy `makeplot.py` (goldens are never
committed); `makeplot.py` also plots the comparison.

## Adjudications recorded at the P2 gate

- **Contour units.** The P0 report left open whether the digitized h = 0 /
  h = -400 contours are meters. The gate passes with the levels read in the
  model's own head units (h = -400 on the simulated field) — the digitized
  points sit on the simulated contours to well under a cell — so the
  SERGHEI reference and Frehg2 share the same convention and no rescaling
  applies.
- **`specific_storage: 0`.** The SERGHEI-style inputs carry no storage term
  (Kirkland is an unsaturated problem). A nonzero Ss feeds the PCA
  compressibility/consistency-restore feedback on the very dry,
  flat-retention Glendale flank and runs it dry (P2 finding;
  `docs/theory/groundwater.md`, "preserved quirks").
- **`reallocation_surplus: redistribute`** (amendment A7). The saturation
  bulb perched on the clay layer is starved if the post-allocation surplus
  of saturation-adjacent cells is discarded (the legacy default embedded in
  the b2 golden): about a quarter of the injected volume vanished and the
  h = 0 contour never formed. With the redistribute mode the budget closes
  to 0.3 % and both contours pass with margin.

Measured at the gate (2026-07-26): h = 0 max distance 0.160 m, RMS 0.073 m;
h = -400 max 0.170 m, RMS 0.088 m; mass balance 0.3 %.
