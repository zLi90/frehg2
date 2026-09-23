# g5 reference digitization record (v2 plan §6.4)

Source: Geng, X. & Boufadel, M. C. (2015), *Numerical modeling of water flow
and salt transport in bare saline soil subjected to evaporation*, Journal of
Hydrology 524:427-438. Code-to-code reference (MARUN); no dataset was
published — the curves exist only as figures (v2 plan §3.4).

Produced by `tools/digitize_geng2015.py` on 2026-09-23 against the paper PDF
(not committed; SHA-256 recorded below). Re-running the script against the
same PDF reproduces these CSVs byte-identically.

## Method: exact vector extraction, not pixel picking

The §6.4 protocol budgeted for raster digitization (~20-40 points per curve).
Inspection of the PDF's content streams showed the gate figures are **vector
line art**, so the extraction is exact up to axis calibration:

| Figure | Vector form | Extraction |
|---|---|---|
| Fig. 3 ER curve | stroked polyline, 18,033 vertices in 3 chained subpaths, x strictly monotone | vertices taken verbatim; resampled to a 0.1 h grid |
| Fig. 4 profiles | 26 cubic Béziers per curve, each uniquely identified by its (dash pattern, gray, stroke width) triple matching the legend 1:1 | Béziers sampled at 25 points per segment; clipped to the axes box; interpolated to a 5 mm elevation grid |
| Fig. 9b profiles | same construction; base case K₀ = 9e-4 m/s = the two solid curves (black 20 h, gray 50 h) | as Fig. 4 |
| Fig. 7 | **raster** (embedded image) | not digitized — g5(iv) is an inequality between two model runs; the figure fixes only its sign and scale (V2-A14) |

Axis calibration: a linear map least-squares fitted to the printed tick-label
centers of each axis. Fit residuals (max over labels): Fig. 3 x 0.057 pt,
y 0.004 pt; Fig. 4 x 0.008 pt, y 0.005 pt; Fig. 9b x 0.006 pt, y 0.005 pt.
(1 pt ≈ 0.35 mm on the page; the curves are stroked 0.7-2.3 pt wide, so the
calibration error is far inside the stroke width.)

## Self-checks (plan §6.3 reference sanity — all enforced by the script)

1. **Closed-form anchor:** digitized E(0.1 h) = 1.4474e-7 m/s vs the Table-1
   bulk-aerodynamic closed form 1.470e-7 m/s (V2-A13 derivation) — ratio
   0.985, required within 5 %.
2. **Known-value curve:** the Fig. 4 "0 h" profile is the exactly saturated
   initial condition; the digitized curve reads mean moisture ratio 0.9997
   (required 1.000 ± 0.01).
3. **Deep tails:** every Fig. 4 profile holds moisture ratio 1.000 below
   z = 1.70 m; both Fig. 9b profiles hold the initial 25.0 g/L there
   (required ± 2 g/L).
4. **g5(ii) peak precondition:** the digitized 50 h salinity peak is
   106.8 g/L > 60 g/L, so the near-surface-peak criterion is readable off
   the reference.

Key digitized scalars (the V2-A13 g5(i) reference values):
E(10 h) = 3.317e-8 m/s; E(50 h)/E(0) = 4.24e-2; Fig. 4 surface moisture
0.207 / 0.157 / 0.134 / 0.114 at 10 / 20 / 30 / 50 h.

## Known limitation

The Fig. 4 legend lists six curves; only five exist as vector paths in the
PDF (0/10/20/30/50 h). The **40 h** curve has no vector content (checked
exhaustively over all 26-segment strokes on the page) and is left empty in
`fig4_moisture_ratio.csv`. g5 gates 20 h and 50 h, so nothing is lost.

## Files

- `fig3_evaporation_rate.csv` — time_h (0.1-50, step 0.1), evaporation_m_per_s
- `fig4_moisture_ratio.csv` — elevation_m (1.60-2.00, step 0.005), moisture ratio at 0/10/20/30/50 h
- `fig9b_salinity.csv` — elevation_m (same grid), salinity [g/L] at 20/50 h (base case)
- `overlay_fig3.png`, `overlay_fig4.png`, `overlay_fig9b.png` — digitized
  points re-plotted over the original figure rendering: **the owner check.**
  Look for the open circles riding their curves; sign-off is recorded in the
  Q4 DoD.

PDF provenance: `references/Geng_and_Boufadel_JH2015.pdf` (local, not
committed), SHA-256
`d3af5ef38a3b132ce29d229a777ee1d128e6bb0a9117373b900d7267bd003717`.
