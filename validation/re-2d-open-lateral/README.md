# re-2d-open-lateral — 2D infiltration with open lateral boundaries

Groundwater-only (Richards) validation case converted from the SERGHEI-RE
reference deck. A 4 km-wide, 15 m-deep aquifer slice receives a multi-year
atmospheric infiltration/evaporation signal on top while its two lateral
faces are held at fixed (Dirichlet) heads, so the water table fluctuates and
drains laterally. The simulated water-table position is compared against the
Hydrus-2D reference solution shipped with the deck (`HYDRUS2D.csv`,
`hydrus-time.csv`).

## Geometry and soil

- 40 x 1 x 60 vertical x-z slice; `dx = dy = 100 m` (`dem.input` cellsize),
  `dz = 0.25 m`, 15 m tall (`subsurface.input`: height 15, dz_base 0.25,
  ndepth 60). Domain is 4000 m wide.
- Homogeneous van Genuchten soil (`vg.input`): Ks = 5.78e-4 m/s,
  theta_s = 0.45, theta_r = 0.1, alpha = 1.65 /m, n = 2.0, air-entry 0.
- `specific_storage = 0`; `reallocation_surplus = redistribute`.

## Geometry caveat

The SERGHEI deck uses an **inclined / terrain-following mesh** — the surface
in `dem.input` drops from 3.9 m to 0.0 m across the slice. frehg2's
rectilinear grid approximates this as a **flat box** (`follow_terrain:
false`, `bottom_elevation: 0`, z in [0, 15] m). Head fields and lateral
water-table elevations are expressed in this flat-box z frame.

## Initial and boundary conditions

- Initial condition (`initialMode 2`, head file): the exact 3D pressure-head
  field from `head.input`, reproduced as `input/head.dat`
  (`(j*nx+i)*nz+k`, k=0=bottom; source z-rows reversed). It is hydrostatic
  below a water table sloping from eta ~ 7.0 m (west) to ~ 0.9 m (east),
  capped at -1.25 m suction above the table.
- Top (`gwbc` dir 6, bctype 6): atmospheric time-varying flux on
  `groundwater_top`, from `atmosphere.input` -> `input/atmosphere.dat`
  (`t[s]  flux[m/s]`; already SI, negative = infiltration). Units verified:
  source values are m/s (e.g. -4.7e-9 to -1.97e-8), used unchanged.
- West face (`gwbc` dir 2) and east face (`gwbc` dir 1): bctype 2 Dirichlet
  head on `groundwater_side`. `headleft.input` / `headright.input` are fixed
  vertical profiles that are hydrostatic below the water table (constant
  -1.25 m suction above). **Judgment call:** frehg2 head BCs offer no
  fixed-vertical-profile value form, so each is represented as a hydrostatic
  head from its water-table elevation (`eta = 7.005 m` west, `0.905 m`
  east). This reproduces the saturated-zone head exactly and only
  approximates the -1.25 m suction cap in the (near-inert) unsaturated zone.

## Reference solution

Hydrus-2D water-table time series (`HYDRUS2D.csv` / `hydrus-time.csv` in the
source case). Compare the frehg2 water-table elevation vs time against these.

## Scheme note

`gw_scheme 1` (PCA) -> `pca`.

## Source

`serghei-re-tests/case4-2d-infiltration-with-open-lateral-boundaries/input-uniform-dz/` (SERGHEI supplement; not shipped in this repository — see `validation/README.md`)
(parameters, subsurface, vg, gwbc, soilID, head, headleft, headright,
atmosphere, polygon* `.input`; dem `.input`).

## Duration note (deck vs spec)

`t_end = 157680000 s` (5 yr, `parameters.input`). `dt_max = 300 s` from
`subsurface.input` — **not** the 3600 s assumed in the conversion spec.
The deck is ground truth, so 300 s is used.

## How to run

```
cd validation/re-2d-open-lateral
OMP_NUM_THREADS=1 ../../build/src/frehg re-2d-open-lateral.yaml
```

Validate only: append `--validate` before the yaml path.

## Compute cost

Groundwater costs ~2.0e-6 s per cell-step (`OMP_NUM_THREADS=1`).
Cell-steps ~ nx*ny*nz * (t_end / dt_max) = 2400 * (157680000 / 300) =
1.26e9, so ~42 min of wall time. dt_max is conservative (adaptive steps
often larger), so treat this as an upper bound. **Recommendation: HPC**
(or a long local/overnight run) since the estimate exceeds the 10-minute
local threshold.

Note: relaxing the controller to `dt_max = 3600 s` (the value the spec
assumed) drops the estimate to ~3.5 min (LOCAL) at the cost of coarser
time resolution; the deck's 300 s ceiling is retained here for fidelity.
