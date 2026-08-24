# goswami-clement-swi — laboratory salt-wedge intrusion (Goswami & Clement 2007)

A **variable-density subsurface** solute-transport validation case (groundwater +
transport + density coupling, no surface module). Where `../transport-henry-saltwater-intrusion/`
is the semi-analytical benchmark, Goswami & Clement is the standard **laboratory**
benchmark for density-dependent groundwater codes: a homogeneous glass-bead flow
tank with a freshwater reservoir landward and a saltwater reservoir seaward,
reaching a **steady salt wedge** whose toe was measured by image analysis. It is
the case people cite when validating SEAWAT / SUTRA / FEFLOW against a *real*
experiment rather than a formula.

## The problem

The 0.53 m (landward → seaward, `y`) × 0.26 m (vertical, `z`) flow chamber,
`nx=1, ny=53, nz=26` (`dy = dz = 0.01 m`), homogeneous, fully saturated
(confined — both reservoir heads sit above the tank top).

| Goswami–Clement parameter | value | frehg2 mapping |
|---|---|---|
| hydraulic conductivity `Ks` | 1.09e-2 m/s | `soil.ksx/ksy/ksz` |
| porosity `θs` | 0.385 | `soil.theta_s` |
| density ratio `Δρ/ρ` | ~0.026 (NaCl, ρₛ≈1026) | seawater `s_sea = 35` |
| boundaries | constant-head reservoirs | seaward (y+) & inland (y−) `head` |

The saltwater reservoir (seaward, y+) and freshwater reservoir (inland, y−) are
set as constant hydrostatic heads; the small inland-over-seaward head excess
(here **13.5 mm**) drives net fresh outflow while density drives the intruding
wedge. `s_sea = 35` is the density surrogate for `Δρ/ρ ≈ 0.026` under frehg2's
hardwired `β_ρ = 7.44e-4` (`r_ρ = 1 + β_ρ·s`).

## Reference solution

Goswami, R. R. & Clement, T. P. (2007), *Laboratory-scale investigation of
saltwater intrusion dynamics*, **Water Resour. Res. 43(4)**, W04418,
[doi:10.1029/2006WR005151](https://doi.org/10.1029/2006WR005151). The measured
metric is the **steady 0.5-isochlor toe length** from the saltwater boundary; the
steady runs (SS1–SS3) sit in the **~0.24–0.26 m** range for their intermediate
head settings.

## Measured result (OMP_NUM_THREADS=1)

| quantity | model | reference |
|---|---|---|
| steady state | ΔC = 0 over the last 1000 s | — |
| salinity range | [0, 35.0] psu (fully fresh landward) | [0, 35] |
| 0.5-isochlor toe (base) | **0.285 m from the sea** | ~0.24–0.26 m |
| wedge character | seawater along the base, fresh discharge over the top | classic SWI wedge ✓ |

The steady wedge is well-formed and its toe lands ~1 cell beyond the experimental
band — a good match given the caveats below.

## Caveats

- **Extreme toe sensitivity to boundary heads.** This is a documented feature of
  the experiment itself: in this model a **sub-millimetre** change in the inland
  head swings the toe by ~5 cm (e.g. Δh = 15 mm → toe 0.205 m; Δh = 13.5 mm →
  0.285 m). Do not read the toe to better than ~±2 cells. The head excess here
  was tuned to bracket the experimental range, not fitted to a single cell.
- **β_ρ and the viscosity ratio are compile-time constants**, so the density
  contrast is set purely through `s_sea` (same as Henry).
- **Dispersion quirk.** The dispersion term acts on volumetric face fluxes, so
  quantitative isochlor spread is approximate; `longitudinal = transverse = 0`
  here with only `molecular = 1e-6` active. This case validates the steady wedge
  *character* and toe *order*, not a to-the-mm isochlor.
- **Seawater enters from y+** (the only side the transport module injects a scalar
  ghost on, `SubsurfaceTransport.cpp` `sideCodeYp`); freshwater enters inland at
  y− where `s = 0`. Same constraint as `../transport-ogata-banks-column/`.

## Cost & recommendation

1378 cells, density coupled, ~8000 groundwater substeps → **~27 s**. **Run
locally.**

## Run

```bash
cd validation/transport-goswami-clement-swi
OMP_NUM_THREADS=1 FI_PROVIDER=tcp ../../build/src/frehg goswami-clement-swi.yaml
```

Output in `out/`. Reshape `transport/concentration` to `[ny=53, nz=26]` (k=0 top,
k=25 bottom); the bottom row's 0.5·s_sea crossing is the wedge toe, measured from
the seaward (y+) end.
