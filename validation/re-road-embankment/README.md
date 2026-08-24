# re-road-embankment — Chavez-Negrete et al. (2018) road-embankment infiltration

Groundwater-only (Richards) validation case converted from the SERGHEI-RE
reference deck. Rainfall-driven infiltration through the crest of a road
embankment advances a wetting front into an initially dry fill sitting above
a regional water table. The simulated pressure-head / wetting-front field is
compared against the GFDM (generalized finite-difference method) reference of
Chavez-Negrete, Domínguez-Mota & Santana-Quinteros (2018).

## Geometry and soil

- 90 x 1 x 50 vertical x-z slice; `dx = dy = 0.1 m` (`dem.input` cellsize),
  `dz = 0.1 m`, 5 m tall (`subsurface.input`: height 5, dz_base 0.1,
  ndepth 50). Domain is 9.0 m wide.
- Homogeneous van Genuchten soil (`vg.input`): Ks = 2.89e-6 m/s,
  theta_s = 0.37, theta_r = 0.04, alpha = 0.873 /m, n = 1.57, air-entry 0.
- `specific_storage = 0`; `reallocation_surplus = redistribute`.

## Geometry caveat

The SERGHEI deck uses a **sloped / terrain-following embankment mesh** — the
surface in `dem.input` is flat at 0.0 m for x < 3 m then drops to -3.0 m.
frehg2's rectilinear grid approximates it as a **flat box**
(`follow_terrain: false`, `bottom_elevation: -5.0`, z in [-5, 0] m with the
flat-approach surface at z = 0).

## Initial and boundary conditions

- Initial condition (`initialMode 4`, water-table file): `wt.input` is a
  uniform -5.0 m water-table elevation -> `water_table: {constant: -5.0}`,
  a flat table at the box bottom with the column starting unsaturated above.
- Top (`gwbc` dir 6, bctype 2, bcvals 1e-5): Dirichlet head of 1e-5 m
  (thin water layer) over the infiltrating strip `polygontop`
  (x in [1.99, 9.01] m) on `groundwater_top`.
- Bottom (`gwbc` dir 5, bctype 2): Dirichlet head on `groundwater_bottom`.
  **Judgment call:** `headbot.input` is hydrostatic from the original water
  table — over the sloped bottom it reads 0 -> 3 m, and
  `headbot(i) = -5.0 - z_bottom(i)` recovers one water-table elevation,
  eta = -5.0 m. In the flat box the bottom face sits at absolute z = -5.0,
  where that table gives a uniform pressure head of 0, so the BC is a
  `constant: 0.0` head. frehg2 forbids a hydrostatic value on
  `groundwater_bottom` (side-only) and offers no per-column value form, so
  the x-variation of `headbot` (a pure terrain effect) is approximated away
  here — consistent with the flat-box geometry caveat above.

## Reference solution

Chavez-Negrete, C., Dominguez-Mota, F.J. & Santana-Quinteros, D. (2018),
*Numerical solution of Richards' equation of water flow by generalized finite
differences*, Computers and Geotechnics 101 — road-embankment infiltration.
Compare the simulated pressure-head field / wetting-front position.

## Scheme note

The source deck uses `gw_scheme 2` (Picard), which maps to frehg2's sole
groundwater scheme, `pca`.

## Source

`serghei-re-tests/case5-2d-infiltration-into-a-road-embankment/input/` (SERGHEI supplement; not shipped in this repository — see `validation/README.md`)
(parameters, subsurface, vg, gwbc, soilID, wt, headbot, polygon* `.input`;
dem `.input`). `t_end = 129600 s = 1.5 d` (`parameters.input`),
`dt_max = 10 s`.

## How to run

```
cd validation/re-road-embankment
OMP_NUM_THREADS=1 ../../build/src/frehg re-road-embankment.yaml
```

Validate only: append `--validate` before the yaml path.

## Compute cost

Groundwater costs ~2.0e-6 s per cell-step (`OMP_NUM_THREADS=1`).
Cell-steps ~ nx*ny*nz * (t_end / dt_max) = 4500 * (129600 / 10) = 5.83e7,
so ~2 min of wall time. dt_max is conservative (adaptive steps often
larger), so treat this as an upper bound. **Recommendation: LOCAL.**
