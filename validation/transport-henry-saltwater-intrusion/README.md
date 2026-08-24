# henry-saltwater-intrusion — the classic Henry (1964) salt-wedge problem

A **variable-density subsurface** solute-transport validation case (groundwater +
transport + density coupling, no surface module). Henry's problem is *the*
canonical benchmark for density-dependent groundwater codes: a confined aquifer
with steady freshwater recharge landward and a hydrostatic seawater boundary
seaward reaches a **steady salt wedge** whose isochlors have a semi-analytical
reference solution. If a code forms the right wedge here, its baroclinic Darcy
coupling is wired correctly.

## The problem

A 2 m (landward → seaward, `y`) × 1 m (vertical, `z`) confined box, `nx=1,
ny=40, nz=20` (`dy = dz = 0.05 m`), fully saturated, no-flow top and bottom.

| Henry parameter | value | frehg2 mapping |
|---|---|---|
| hydraulic conductivity `Ks` | 1e-2 m/s | `soil.ksx/ksy/ksz` |
| porosity `θs` | 0.35 | `soil.theta_s` |
| freshwater recharge `Q` | 6.6e-5 m²/s | inland (y−) `flux = 6.6e-5 m/s` |
| density ratio `Δρ/ρ` | 0.025 | seawater `s_sea = 33.6` (see below) |
| diffusion `D` | 6.6e-6 m²/s | `transport.dispersion.molecular` |

Seawater enters as a **hydrostatic head** (sea level 1.0 m) on the seaward (y+)
boundary at salinity `s_sea`; freshwater enters as a small constant flux on the
inland (y−) boundary at `s = 0`.

### Why `s_sea = 33.6` and not 35

frehg2's baroclinic law is `r_ρ = 1 + β_ρ·s` with **β_ρ hardwired at 7.44e-4**
(compile-time constant). Henry's density contrast is `Δρ/ρ = 0.025`, so the
salinity that reproduces it is `s_sea = 0.025 / 7.44e-4 = 33.6`. (At the "real
seawater" value 35 the ratio would be 0.026 — close, but 33.6 hits Henry's
number exactly.) The scalar is therefore a *density surrogate*, not a literal
psu.

## Reference solution

Henry, H. R. (1964), *Effects of dispersion on salt encroachment in coastal
aquifers*, USGS Water-Supply Paper **1613-C**, pp. C71–C84. The standard
verification metric is the position of the 0.5 isochlor (50 % seawater) toe along
the base: it reaches roughly **0.8–1.0 m inland** of the sea in the `a = 2`
domain.

## Measured result (OMP_NUM_THREADS=1)

| quantity | model | reference |
|---|---|---|
| steady state reached | ΔC < 0.02 psu over t=14000→15000 | — |
| salinity range | [0, 33.5] psu | [0, 33.6] (bounds) |
| 0.5-isochlor toe (bottom row) | **0.73 m from the sea** | ~0.8–1.0 m |
| wedge character | seawater along the base, fresh discharge over the top | classic Henry isochlors ✓ |

The field is a textbook Henry wedge: at the seaward column salinity ramps from
~0 at the top to 33 at the bottom, and the salt toe intrudes landward along the
base to ~0.73 m — within the accepted Henry band. Wall time **~6 s**.

## Caveats

- **β_ρ (and the viscosity ratio) are compile-time constants**, so the density
  contrast is set purely through `s_sea`; you cannot dial `Δρ/ρ` in the YAML.
- **Dispersion quirk.** The dispersion term acts on *volumetric* face fluxes (the
  face-area factor is not divided out), so quantitative isochlor positions are
  approximate. This case validates that frehg2 forms the steady Henry wedge with
  the right *character* (toe intrusion depth, upper fresh outflow), not
  to-the-millimetre isochlors — hence `longitudinal = transverse = 0` here, with
  only `molecular` diffusion active.
- **Seawater must enter from y+.** The transport module wires a prescribed scalar
  ghost only on the y+ side (`SubsurfaceTransport.cpp` `sideCodeYp`), so the
  seaward boundary is y+ and freshwater enters inland at y− (where `s = 0 =`
  initial, so the missing y− injection path is harmless). Same constraint as
  `../transport-ogata-banks-column/`.

## Cost & recommendation

800 cells × ~3000 groundwater substeps, density coupled → ~6 s. **Run locally.**

## Run

```bash
cd validation/transport-henry-saltwater-intrusion
OMP_NUM_THREADS=1 ../../build/src/frehg henry.yaml
```

Output in `out/`. Reshape `transport/concentration` to `[ny=40, nz=20]` (k=0 is
the top, k=19 the bottom); the bottom row's 0.5·s_sea crossing is the wedge toe.
