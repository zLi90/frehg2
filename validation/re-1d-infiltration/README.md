# re-1d-infiltration — Warrick et al. (1985) 1D infiltration

Groundwater-only (Richards) validation case converted from the SERGHEI-RE
reference deck. A homogeneous, initially dry 1 m soil column is infiltrated
from the top by a thin ponding layer, and the moisture/head profiles are
compared against the Warrick, Lomen & Islas (1985) analytical (quasi-linear)
solution for 1D infiltration into a dry column.

## Geometry and soil

- 1 x 1 x 100 vertical column; `dx = dy = 1.0 m`, `dz = 0.01 m`, 1.0 m tall
  (`subsurface.input`: height 1.0, dz_base 0.01, ndepth 100; `dem.input`
  cellsize 1.0). Flat box, `bottom_elevation = 0`.
- Homogeneous van Genuchten soil (`vg.input`): Ks = 2.89e-6 m/s
  (isotropic), theta_s = 0.33, theta_r = 0.0, alpha = 1.43 /m, n = 1.56,
  air-entry value 0.
- `specific_storage = 0` (the deck carries no storage term; unsaturated
  problem — same choice as b3-kirkland).
- `reallocation_surplus = redistribute` (mass-conserving surplus handling,
  so the infiltrated volume is preserved for the analytical comparison).

## Initial and boundary conditions

- Initial condition (`initialMode 3`, `theta.input` is uniform 0.03):
  `moisture: {constant: 0.03}` — a dry start.
- Top (`gwbc` dir 6, bctype 2, bcvals 1e-5): Dirichlet head of 1e-5 m
  (a thin ponding layer ~ 0) on `groundwater_top`.
- Bottom: no-flux (default; the deck sets no bottom BC).

## Reference solution

Warrick, Lomen & Islas (1985), *Water Resources Research* — analytical
quasi-linear solution for constant-head infiltration into a semi-infinite
dry column. Compare the simulated water-content / pressure-head profiles at
the output times against the analytical front position and shape.

## Scheme note

The source deck uses `gw_scheme 1` (PCA), which maps to frehg2's sole
groundwater scheme, `pca`.

## Source

`serghei-re-tests/case1-1d-infiltration/input/` (SERGHEI supplement; not shipped in this repository — see `validation/README.md`)
(parameters, subsurface, vg, gwbc, soilID, theta, polygontop, dem `.input`).
This is a feasible, clean conversion (no geometry approximation).

## How to run

```
cd validation/re-1d-infiltration
OMP_NUM_THREADS=1 ../../build/src/frehg re-1d-infiltration.yaml
```

Validate the config without running:

```
OMP_NUM_THREADS=1 ../../build/src/frehg --validate re-1d-infiltration.yaml
```

## Compute cost

Groundwater costs ~2.0e-6 s per cell-step (`OMP_NUM_THREADS=1`).
Cell-steps ~ nx*ny*nz * (t_end / dt_max) = 100 * (46800 / 0.4) = 1.17e7,
so ~23 s of wall time. dt_max is conservative (adaptive steps often larger),
so treat this as an upper bound. **Recommendation: LOCAL.**

A coarser-controller variant with `dt_max = 10 s` (mentioned in the deck
context) drops this to ~1 s; the 0.4 s ceiling is kept here to resolve the
sharp wetting front for the analytical comparison.
