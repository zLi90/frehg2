# ogata-banks-column — 1-D advection–dispersion breakthrough (saturated column)

A **subsurface** solute-transport validation case (groundwater + transport, no
surface flow, no density). It reproduces the canonical **Ogata & Banks (1961)**
problem: continuous injection of a conservative tracer into a semi-infinite
saturated column carrying steady uniform seepage, which has the closed-form
complementary-error-function (erfc) solution.

## The problem

A 1.0 m column resolved along `y` (`nx=1, ny=100, nz=1`, `dy = 0.01 m`), fully
saturated homogeneous sand (`Ks = 1e-4 m/s`, `θs = 0.35`). A pressure-head drop
of 1.0 m across the column drives a steady Darcy flux `q = Ks·dH/L = 1e-4 m/s`;
with `θs = 0.35` the pore velocity is `v = q/θs = 2.86e-4 m/s`, so the mean
tracer arrival time at the far end is `t50 = L/v ≈ 3500 s`. At `t = 0` the
inflow-side concentration steps to `c0 = 1`; the far-end breakthrough curve and
the interior concentration profile are compared against the Ogata–Banks erfc
solution.

## Reference solution

Ogata, A. & Banks, R. B. (1961), *A solution of the differential equation of
longitudinal dispersion in porous media*, USGS Professional Paper **411-A**
([pubs.usgs.gov/pp/0411a/report.pdf](https://pubs.usgs.gov/pp/0411a/report.pdf)):

```
c(x,t)/c0 = ½ erfc[(x − v t)/(2√(D t))] + ½ exp(v x/D) erfc[(x + v t)/(2√(D t))]
```

(the second term is negligible here, giving the simple `½ erfc` front).

## Measured result (OMP_NUM_THREADS=1)

| quantity | model | reference | error |
|---|---|---|---|
| breakthrough `t50` at the outlet | 3431 s | `L/v` = 3500 s | −2.0 % |
| interior profile vs erfc (`t = 1750 s`) | — | `½ erfc`, `D_eff = 2.5e-6 m²/s` | RMS 0.004 |

The advective arrival time matches `L/v` to 2 %, and the front shape overlays the
analytical erfc almost exactly (e.g. at the front `xi ≈ 0.485 m`, model 0.573 vs
erfc 0.564). Wall time is **< 1 s**.

## Important frehg2 setup constraints (learned building this case)

- **The inflow must be on the y+ (north) side.** The subsurface transport module
  wires a prescribed side-ghost concentration only to the y+ boundary
  (`src/transport/SubsurfaceTransport.cpp` `sideCodeYp`, legacy `s_yp`); there is
  no injection path on y-, x-, or x+. An inflow placed on any other side carries
  **zero** tracer (the water enters but the scalar does not). This mirrors
  `b6-kuan`, whose salt enters from its seaward y+ boundary. Flow here therefore
  runs north → south, and breakthrough is read at the south (`j = 0`) end.
- **Dispersion quirk.** frehg2 applies the dispersion term to *volumetric* face
  fluxes, not Darcy velocities (the face-area factor is not divided out). So the
  configured `dispersion.longitudinal = 0.1` is **not** the physical dispersivity
  — it is calibrated to the face-flux scale to produce a realistic erfc spread.
  The quantitatively exact, quirk-independent check is the advective `t50`; the
  effective `D` is recovered by fitting the front (`D_eff ≈ 2.5e-6 m²/s`).
- The column is fully saturated by a positive initial pressure head
  (`head: {constant: 0.3}`, `θ = θs = 0.35` everywhere), so the transport runs on
  a clean steady flow field.

## Cost & recommendation

100 cells, ~700 groundwater substeps → ~10⁵ cell-steps. **Run locally** (< 1 s).

## Run

```bash
cd validation/transport-ogata-banks-column
OMP_NUM_THREADS=1 ../../build/src/frehg ogata-banks-column.yaml
```

Output in `out/`. Compare `transport/concentration` at the south end against
`t50 = L/v`, and the interior profile against the `½ erfc` solution.
