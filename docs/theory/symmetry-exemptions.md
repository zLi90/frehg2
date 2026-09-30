# Symmetry exemptions (v2 plan §8.1)

The dihedral orientation batteries (`regression.swe_orient`,
`regression.gw_orient`, `regression.transport_orient`,
`regression.wind_orient`, `regression.heat_orient`) run each base case
through the 8 dihedral transforms of the square — identity, `transpose`,
and the six *mixed* transforms (`flipx`, `flipy`, `rot180`, `rot90`,
`rot270`, `antitranspose`) — and require the inverse-transformed output to
match the identity run. This page is the plan §8.1 **documented-asymmetry
exemption table**: the only mechanism by which a symmetry deviation may be
excused. Every entry is a legacy-faithful behavior, waived *per term, with
provenance*; an over-tolerance deviation, or a deviation without a row
here, is a bug. New v2 physics (evaporation, heat exchange, wind stress,
Cd laws) has **no exemptions** — it is exactly symmetric by construction,
and its batteries gate at the strict bound.

## Why `transpose` is strict and the mixed class is not

`transpose` maps x− ↔ y− and x+ ↔ y+, preserving each domain edge's
minus/plus class, so it exercises scheme/indexing orientation symmetry
without touching any of the preserved legacy minus/plus differences below.
It therefore gates at the strict rank-invariance floor (1e-12, strict
solver mode). The six mixed transforms map at least one minus edge onto a
plus edge; for the surface module the preserved legacy edge arithmetic
*differs by design* between the classes, so mixed-class agreement has a
measured floor above round-off. The mixed bound
(`tests/regression/tolerances/x-orient.yaml`) is calibrated above that
floor and far below any V2-A11-class side defect (~0.1 m).

## Exemption table

| # | Waived asymmetry | Provenance (legacy) | Frehg2 site | Batteries it floors |
|---|---|---|---|---|
| 1 | **Closed-edge explicit flux**: the RHS divergence reads the ghost-side E at west/south faces (zero) but the owned-cell E at east/north faces; the implicit legs are folded on all four edges. | `shallowwater_rhs` (legacy `shallowwater.c:367-393` fold; the E-read asymmetry is the preserved v1 scheme) | `FreeSurface.cpp` RHS assembly | swe, transport, heat (mixed class) |
| 2 | **Pre-source ghost stage at boundary faces**: physical-edge ghost η refreshes only inside the free-surface phase, so rain/evaporation induce a small outward boundary-face velocity at live east/north default faces each step (a flat "closed" basin leaks slowly there, not at west/south). | legacy call order `solve.c:51-98` | `SurfaceSolver` phase order | swe, transport (mixed class, forced cases) |
| 3 | **Outflow coefficient, plus-reuse vs minus-recompute**: an east/north transmissive outlet reuses the assembled interior face coefficient (`Sxp`/`Syp`); a west/south outlet recomputes `g·dt²·A²·D/V` from the halo-column face geometry (post-V2-A11 the true boundary face). The two evaluations agree only to the geometry/rounding floor. | legacy `shallowwater.c` outflow row; V2-A11 rewrote the minus form | `FreeSurface.cpp::applyOutflowCorrections` | swe (mixed class; `regression.outflow_staircase` gates both forms against the same absolute analytic bands) |
| 4 | **Velocity-condition edge slots**: the wall/velocity correction restores west/south edge slots through a reconstruct-and-restore path whose rounding differs from the plus-side direct write (bitwise-restart-preserving, V1 P3). | `enforceVeloBc` restore order (P3) | `WetDry.cpp::enforceVeloBc` | swe (mixed class, walled cases) |
| 5 | **Eddy-viscosity face areas**: difX uses `Asx(i)` on both x-faces and difY uses `Asy(j)` on both y-faces (one-sided face-area reuse), so mirrored shear layers see the area of the opposite face. | legacy `shallowwater.c:170-190` | `Momentum.cpp` (`momentumSource`) | swe (mixed class, sheared cases) |
| 6 | **Subsurface hydrostatic side-ghost head**: coupled live ghost head is `(bed − z_c)·r_face + depth·r_face` at plus sides but `(bed − z_c)·r_face + depth` at minus sides (the density factor multiplies the depth term on one class only). | legacy `enforce_head_bc:769-788` (`:777` vs `:787`), preserved by A18 | `gw/Predictor.cpp::ghostHead` | gw/transport coupled saline cases (dormant in the gw battery base: uncoupled, uniform density) |
| 7 | **y+-only scalar side-ghost admission (salinity)**: the legacy limiter admits prescribed side ghosts into the wet-stencil min/max on y+ ONLY (the b6 sea rule); salinity Dirichlets on x−/x+/y− are limiter-clipped. Temperature admits all four sides (V2-A18). | legacy `scalar.c:397-400` | `SubsurfaceTransport.cpp` limiter ghost admission (salinity spec) | transport (mixed class, side-value cases; dormant in the battery base — no side scalar BC) |
| 8 | **BoundarySet corner-face spill**: a side polygon spanning a domain corner hands its condition to *every* domain-edge face of the corner cell. Load-bearing for b6's 1-wide tank (V2-A18); makes a corner-adjacent BC's footprint direction-dependent under rotation. | v1 BoundarySet rasterization, pinned by b6 | `bc/BoundarySet.cpp` | all batteries (bases avoid corner-spanning polygons by design) |
| 9 | **Uncoupled top-face advection**: superbee runs pass the donor value through head-condition tops while upwind runs zero the inflow side — a z-axis (not dihedral) asymmetry between schemes, listed for completeness. | legacy `scalar.c:653-697`, derived branch outcomes | `SubsurfaceTransport.cpp` top-face branch | none (z is never transformed) |

Rows 6–8 are *dormant by design* in the battery bases (uncoupled,
uniform-density, no side scalar values, no corner-spanning polygons), so
they do not contribute to the measured floors; they are listed because the
batteries must not be "fixed" into exercising them without adding the
matching exemption bound first. Rows 1–5 set the surface mixed-class
floor.

## Measured floors (release record, v2.0.0, 2026-09-26, macOS arm64)

Measured by the batteries on the release tree (strict solver mode; each
inverse-transformed variant against the identity run; max over the listed
fields):

| Battery | transpose class (measured / bound) | mixed class (measured / bound) | Mixed-class floor set by |
|---|---|---|---|
| swe (`eta`, `depth`) | 4.7e-15 / 1e-12 | 3.478e-3 / 1e-2 | rows 1–5; per-axis and additive: x-edge swaps (flipx, rot270) 2.227e-3, y-edge swaps (flipy, rot90) 1.336e-3, both (rot180, antitranspose) 3.478e-3 |
| gw (`hydraulic_head`, `water_content`) | 1.3e-15 / 1e-12 | 3.75e-14 / 1e-12 | none — the groundwater side arithmetic is exactly symmetric (both classes gate strict) |
| transport (`eta`, `concentration_surface`) | 2.1e-15 / 1e-6 | eta 3.09e-10, concentration 2.17e-9 / 1e-6 | x-edge swaps only (y-swaps are at machine precision): the conveyor's live-edge terms (row 1) shift eta at 3e-10 and the limiter/slope-select branches amplify that into the scalar at ~7× |
| wind (`eta`) | ≤ 1.1e-15 m / 1e-7 | ≤ 1.1e-15 m / 1e-7 | none — Q6 physics, no exemptions |
| heat (`t_subsurface` ≤ 1.5e-13, `head` ≤ 1.8e-15) | / 1e-12 | / 1e-12 | none — Q5 physics, no exemptions (V2-A18 landed all-side admission for temperature) |

The transport bound is 1e-6 in both classes although the transpose class
measures machine precision: the class distinction protects the *edge*
arithmetic, while the limiter's minmod/slope-select branches can flip on
round-off-level differences in either class and amplify them to
O(local flux·dt); the shared bound absorbs that mechanism wherever it
lands, three orders above the measured mixed floor and six below the
initial-condition contrast.

The g9 wind battery's exact symmetry (≤ 1.1e-15 m against a 1e-7 bound)
is the constructive proof that the v2-physics "no exemptions" rule is
attainable: forcing terms written once against face geometry, without
per-side branches, are symmetric to round-off.
