# surface-tracer-advection — passive tracer in uniform open-channel flow

A **surface** solute-transport validation case (SWE + transport, no groundwater,
no density). It is the surface analogue of the classic 1-D advection–dispersion
verification: a passive tracer carried by a steady, uniform open-channel flow.
Its job is to validate frehg2's **surface advection schemes** — first-order
`upwind` vs. second-order TVD `superbee` — against the exact answer.

Two YAMLs share this directory and `input/`, differing only in the advection
scheme:

| file | scheme | what it shows |
|---|---|---|
| `advect-upwind.yaml`   | `upwind`   | strongly smeared front (numerical diffusion) |
| `advect-superbee.yaml` | `superbee` | sharp front, minimal smearing |

## The problem

A 1000 m × 1 m channel (`nx = 200`, `dx = 5 m`) on a mild 0.1 % slope. With
Manning `n = 0.03` the normal depth is `h = 0.5 m` and the uniform velocity is
`u = 0.664 m/s` (unit discharge `q = 0.332 m²/s`). The bed drops linearly from
0.9975 m at the west inlet to 0.0025 m at the east outlet
(`input/dem_channel.dat`); the initial stage is `bed + 0.5 m`
(`input/eta_channel.dat`), matched to the inflow discharge so the flow is
uniform from `t = 0`.

At `t = 0` the inflow concentration steps from 0 → 1. With no physical diffusion
(`surface_diffusivity ≈ 1e-10`), the exact solution is a **step front at
`x = u·t`** (pure advection). Any spreading of the front is therefore purely
*numerical* diffusion introduced by the scheme.

## Reference solution

- Pure advection: a sharp front translating at `u`; front position `x_f = u·t`.
- With dispersion `D`, the front is the Ogata & Banks (1961) complementary-error
  -function profile `c(x,t) = ½ erfc[(x − u·t)/(2√(Dt))]`.
  Ogata, A. & Banks, R. B. (1961), *A solution of the differential equation of
  longitudinal dispersion in porous media*, USGS Professional Paper **411-A**.
  (The same 1-D ADE governs a uniform surface flow; see also Fischer et al.,
  *Mixing in Inland and Coastal Waters*, Academic Press, 1979, ch. 5.)

## Measured result (OMP_NUM_THREADS=1)

At `t = 1000 s` the exact front is at `u·t = 664 m`. Both schemes conserve mass
exactly (`c_max = 1.000`) and put the `c = 0.5` front at **657.5 m** — within one
cell of exact. They differ only in front sharpness:

| scheme | front (c=0.5) | 10–90 % smear width | wall time |
|---|---|---|---|
| upwind   | 657.5 m | 140 m (28 cells) | ~1 s |
| superbee | 657.5 m |  15 m ( 3 cells) | <1 s |

Superbee keeps the front ~9× sharper — the expected TVD-vs-donor-cell contrast.

## Notes on the frehg2 setup

- The tracer is injected with a `scalar_value` condition placed on the **same
  polygon as the `discharge` inflow**: per the P4 semantics
  (`docs/user-guide/parameters.md`) a `scalar_value` member that overlaps
  a `discharge` condition injects that inflow's concentration (legacy
  `s_inflow`), rather than acting as a wet-cell Dirichlet.
- Depth 0.5 m sits well above the `thin_layer_depth` friction-regularization band
  (0.1 m), so the normal-depth flow is clean (no cm-scale offset like the shallow
  SERGHEI steady-plane cases).
- `dt = 1 s` gives an advective Courant number `u·dt/dx ≈ 0.13`; the scalar
  advection is explicit, the free surface semi-implicit.
- `bounds.max: 1.0` clamps the scalar to its physical ceiling (the step value).

## Cost & recommendation

200 cells × 1000 steps ≈ 2×10⁵ cell-steps. **Run locally** — each variant is
~1 s. (Rule of thumb for this collection: > 10 min ⇒ prefer HPC.)

## Run

```bash
cd validation/transport-surface-tracer-advection
OMP_NUM_THREADS=1 ../../build/src/frehg advect-upwind.yaml
OMP_NUM_THREADS=1 ../../build/src/frehg advect-superbee.yaml
```

Outputs land in `out-upwind/` and `out-superbee/`. Compare
`transport/concentration_surface` at the final step against `x_f = u·t = 664 m`.
