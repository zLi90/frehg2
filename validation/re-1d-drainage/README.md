# re-1d-drainage — Abeele (1984) crushed-tuff gravity drainage

Groundwater-only (Richards) validation case converted from the SERGHEI-RE
reference deck. An initially fully saturated 6 m column of crushed Bandelier
tuff drains under gravity through an open bottom. The transient
moisture/head profiles are compared against the Abeele (1984) crushed-tuff
drainage experiment (and companion Richards-equation solutions).

## Geometry and soil

- 1 x 1 x 100 vertical column; `dx = dy = 1.0 m`, `dz = 0.06 m`, 6.0 m tall
  (`subsurface.input`: height 6.0, dz_base 0.06, ndepth 100). Flat box,
  `bottom_elevation = 0`.
- Homogeneous van Genuchten soil (`vg.input`, same VG parameters as case1):
  Ks = 2.89e-6 m/s, theta_s = 0.33, theta_r = 0.0, alpha = 1.43 /m,
  n = 1.56, air-entry value 0.
- `specific_storage = 0`; `reallocation_surplus = redistribute`
  (mass-conserving surplus handling as the column desaturates).

## Initial and boundary conditions

- Initial condition (`initialMode 1`, fully saturated):
  `moisture: {constant: 0.33}` (= theta_s), which seeds a hydrostatic head
  from the bed.
- Bottom (`gwbc` dir 5, bctype 2, bcvals -1e-5): Dirichlet head of -1e-5 m
  (~ 0, an open free-outflow datum) on `groundwater_bottom`.
- Top: no-flux (default; the deck sets no top BC).

## Reference solution

Abeele, W.V. (1984), *Hydraulic testing of crushed Bandelier Tuff*,
LA-10037-MS — column drainage of initially saturated crushed tuff. Compare
the simulated drainage-front descent and moisture-profile evolution.

## Duration note (deck vs spec)

The shipped deck runs `simLength = 8640000 s = 100 days` (the paper's
drainage horizon), **not** the 1-day (86400 s) value assumed in the
conversion spec. The deck is ground truth, so `t_end = 8640000 s` is used
here; `dt_max = 1200 s` (from `subsurface.input`).

## Scheme note

The source deck uses `gw_scheme 2` (Picard), which maps to frehg2's sole
groundwater scheme, `pca`.

## Source

`serghei-re-tests/case2-1d-drainage/input/` (SERGHEI supplement; not shipped in this repository — see `validation/README.md`)
(parameters, subsurface, vg, gwbc, soilID, polygonbot, dem `.input`).
Feasible, clean conversion (no geometry approximation).

## How to run

```
cd validation/re-1d-drainage
OMP_NUM_THREADS=1 ../../build/src/frehg re-1d-drainage.yaml
```

Validate only: append `--validate` before the yaml path.

## Compute cost

Groundwater costs ~2.0e-6 s per cell-step (`OMP_NUM_THREADS=1`).
Cell-steps ~ nx*ny*nz * (t_end / dt_max) = 100 * (8640000 / 1200) = 7.2e5,
so ~1.5 s of wall time. dt_max is conservative (adaptive steps often
larger), so treat this as an upper bound. **Recommendation: LOCAL.**
