# elder-free-convection — density-driven free convection (Elder problem)

A **variable-density subsurface** validation case (groundwater + transport +
density coupling, no surface module, **no imposed through-flow**). The Elder
problem is the hardest of the classic density benchmarks: dense fluid overlying
lighter fluid drives buoyant convection *with no applied flow at all* — the motion
is created entirely by the density field feeding back into the Darcy solve. It is
the definitive test of a code's baroclinic coupling, and it is famously
mesh-/scheme-sensitive (the plume count and whether the centre rises or sinks
*bifurcates* between codes), so it validates the coupling **mechanism**, not a
single reference number.

## ⚠️ Structural limitation — read before using

The canonical Elder problem imposes a **sustained fixed-concentration (Dirichlet)
on the top boundary** over the central half of the box. frehg2's subsurface
transport wires a prescribed scalar ghost **only on the y+ lateral side**
(`SubsurfaceTransport.cpp` `sideCodeYp`) — there is **no top-face (z)
scalar-injection path**. The sustained top source therefore **cannot** be
reproduced. This case instead seeds the salt as an **initial condition** in the
top-centre row (`input/salt_ic.dat`) and lets it convect: a **decaying-source**
free-convection demonstration, not the steady sustained-source Elder. Read the
result as a qualitative check that frehg2 produces gravitationally-driven
fingering with the right character, not a quantitative isochlor match.

## The problem

Full Elder box **600 m (y) × 150 m (z)**, `nx=1, ny=60 (dy=10 m), nz=30
(dz=5 m)`, closed (impermeable sides/base, no-flow lid). Voss & Souza (1987)
parameters:

| Elder parameter | value | frehg2 mapping |
|---|---|---|
| permeability / `Ks` | k=4.845e-13 m² → `Ks≈4.75e-6 m/s` | `soil.ksx/ksy/ksz` |
| porosity | 0.1 | `soil.theta_s` |
| density ratio `Δρ/ρ` | 0.2 (Ra ≈ 400) | salt surrogate `s_max = 269` |
| diffusion `D` | 3.565e-6 m²/s | `transport.dispersion.molecular` |

`s_max = 0.2 / 7.44e-4 = 269` under frehg2's hardwired `β_ρ` (`r_ρ = 1 + β_ρ·s`).
Salt is seeded at `s_max` in the top-centre row over the central 300 m
(`input/salt_ic.dat`, flat `(j·nx+i)·nz+k` order, k=0 = top).

## Reference solution

- Elder, J. W. (1967), *Transient convection in a porous medium*, **J. Fluid
  Mech. 27(3)**, 609–623.
- Voss, C. I. & Souza, W. R. (1987), *Variable density flow and solute transport
  simulation…*, **Water Resour. Res. 23(10)**, 1851–1866,
  [doi:10.1029/WR023i010p01851](https://doi.org/10.1029/WR023i010p01851).

There is no unique reference field (the problem bifurcates); codes are compared on
the *pattern* — symmetric descending plumes on the flanks with central
up-/down-welling — and the salt penetration depth versus time.

## Measured result (OMP_NUM_THREADS=1)

| time | salinity range | max penetration | pattern |
|---|---|---|---|
| 0 yr | [0, 269] | 2 m (seed) | top-centre slug |
| 5 yr | [0, 33] | 52 m | descending |
| 10 yr | [0, 23] | 58 m of 150 | **symmetric two-lobe cell** |

At 10 years the salt has convected to ~58 m depth and the lateral profile at depth
shows a **symmetric two-lobe structure** (central minimum at the box centre,
flanking maxima) — the hallmark Elder two-plume convection cell with central
upwelling. This confirms frehg2's buoyancy coupling drives the correct
free-convection instability. The peak salinity decays (269 → 23) because the
source is not sustained (see the limitation above), so the plumes weaken over time
rather than reaching a steady Elder pattern. Wall time **~110 s** for 10
simulated years.

## Caveats

- **No sustained top source** (structural, above) → decaying-source variant only.
- **β_ρ and the viscosity ratio are compile-time constants**, so the 0.2 contrast
  is set purely through `s_max`, and the constant-viscosity assumption Elder makes
  may not be matched by frehg2's fixed viscosity ratio.
- **Dispersion quirk** (term on volumetric face fluxes) applies.
- **Bifurcating / mesh-sensitive**: plume count and symmetry depend on mesh and
  scheme; do not expect reproducibility across resolutions.

## Cost & recommendation

1800 cells, density coupled; a 10-year decaying-source run is **~110 s — runnable
locally** (this variant). The **canonical sustained-source Elder** — a finer mesh
(≥ 44×44), a proper top-face Dirichlet, and a 20-year horizon to resolve the
steady plume field — is a **research / HPC** endeavour and, as noted, is not
faithfully expressible in the current frehg2 transport BC set. **Recommend HPC**
for any quantitative Elder study; use this local case only as a coupling
smoke-test.

## Run

```bash
cd validation/transport-elder-free-convection
OMP_NUM_THREADS=1 FI_PROVIDER=tcp ../../build/src/frehg elder-free-convection.yaml
```

Output in `out/`. Reshape `transport/concentration` to `[ny=60, nz=30]` (k=0 top);
track the penetration depth and the lateral two-lobe structure over time.
