# Surface-water module (frehg::swe)

The Frehg2 SWE module preserves legacy Frehg's semi-implicit θ-scheme
exactly (upgrade plan §3.1). This page is the P1 provenance table required
by plan §11.1 rule 7 and §12: every preserved equation names the legacy
`file:line` it reproduces and the Frehg2 function that implements it, then
the deliberate deviations and correctness fixes are listed with rationale.

## Scheme summary

Per step (legacy `solve.c:51` and `:97`; `driver::Simulation::run` →
`SurfaceSolver::beginStep/solveFreeSurface/updateVelocity`):

1. save η<sup>n</sup>, evaluate forcings at the end-of-step time;
2. explicit momentum predictor E = (u + Δt(diffusion − damped advection)
   [+ wind]) · D with the point-implicit drag factor
   D = 1/(1 + ½ Δt C_D |u| A_zface/V_face);
3. implicit 5-point free-surface system with `coef = g Δt²`
   (SPD; PETSc CG, options prefix `fs_`);
4. wetting limiter, rainfall/evaporation on η;
5. depth and face-geometry update, drag update, velocity update with
   wet/dry limiters, uy/vx interpolation.

## Equation provenance

| Piece | Legacy provenance | Frehg2 implementation |
|---|---|---|
| Upwind advection with the 0.5–0.7 CFL damping ramp | `shallowwater.c:151-169` | `Momentum.cpp` (`momentumSource`), `SweFormulas.hpp` (`cflDampedAdvection`) |
| Central eddy viscosity (incl. the legacy face-area asymmetry: difX uses Asx(i) on both x-faces, difY uses Asy(j) on both y-faces) | `shallowwater.c:170-190` | `Momentum.cpp` (`momentumSource`) |
| Point-implicit drag factor | `shallowwater.c:192-202` | `SweFormulas.hpp` (`pointImplicitFactor`) |
| Manning drag with the thin-layer exponent switch (2/3 below hD, 1/3 above) | `update_drag_coef`, `shallowwater.c:977-994` | `SweFormulas.hpp` (`manningDrag`), `Momentum.cpp` (`updateDragCoef`) |
| Chezy drag C_D = g/C² (plan §5.8 additive extension; no exponent switch — that is Manning-specific) | — (addition) | `SweFormulas.hpp` (`chezyDrag`) |
| Quadratic wind stress with thin-layer attenuation | `wind_source`, `shallowwater.c:258-292` | `SweFormulas.hpp` (`windStress`, `windStressAttenuated`), `Momentum.cpp` |
| Free-surface right-hand side | `shallowwater_rhs`, `shallowwater.c:296-317` | `FreeSurface.cpp` (`assembleRhs`) |
| Matrix coefficients, dry-cell regularization, closed-edge folds | `shallowwater_mat_coeff`, `shallowwater.c:320-393` | `FreeSurface.cpp` (`assembleCoefficients`) |
| Prescribed-stage (tide) rows | `shallowwater.c:396-411` | `FreeSurface.cpp` (`fillAndSolve`) |
| η clamp, stage enforcement, edge ghosts | `enforce_surf_bc`, `shallowwater.c:502-540` | `FreeSurface.cpp` (`enforceSurfBc`) |
| One-cell wetting limiter + sub-threshold drying | `cfl_limiter`, `shallowwater.c:543-574` | `WetDry.cpp` (`cflLimiter`) |
| Rain/evaporation on η with the final clamp | `evaprain`, `shallowwater.c:577-636` | `SurfaceSources.cpp` (`evapRain`) |
| Velocity update (incl. the double drag factor, see below) | `update_velocity`, `shallowwater.c:719-794` | `Momentum.cpp` (`updateVelocityField`) |
| Velocity limiters (vanished faces, dry-cell outflow, wetting flags) | `shallowwater.c:796-833` | `WetDry.cpp` (`applyVelocityLimiters`) |
| Velocity edge ghosts + stage-boundary velocity correction | `enforce_velo_bc`, `shallowwater.c:891-958` | `WetDry.cpp` (`enforceVeloBc`) |
| uy/vx four-point interpolation | `interp_velocity`, `shallowwater.c:962-974` | `Momentum.cpp` (`interpolateVelocity`) |
| Center/face depths (higher of two surfaces over higher of two bottoms) | `update_depth`, `initialize.c:918-995` | `WetDry.cpp` (`updateDepth`) |
| Cell volumes and face areas/volumes (non-subgrid branch) | `shallowwater.c:1046-1099`, `initialize.c:584-611` | `WetDry.cpp` (`updateGeometry`) |
| Bathymetry offset (see below) | `read_bathymetry`, `initialize.c:142-147` | `SurfaceSolver.cpp` (`readBathymetry`) |
| Elevation/velocity initial state | `ic_surface`, `initialize.c:457-655` | `SurfaceSolver.cpp` (`applyInitialConditions`) |

Named constants preserved with their legacy values (plan §3.1 item 5):
the 0.5–0.7 damping band, air density 1.225 kg/m³ and water density
998 kg/m³ in the wind stress (`SweFormulas.hpp`), linear-solve rtol 1e-8
(`LinearSystem`).

## Behaviors that look like bugs but are preserved scheme

These are load-bearing quirks: the b1 goldens embed them, so they are part
of the algorithm, not defects.

- **The elevation offset is preserved.** Legacy lifts all elevations by
  −min(bottom) so the bed is non-negative (`initialize.c:142-147`) and
  subtracts the offset on output. Under the legacy dry-cell row
  (`Sct = dx dy`, RHS `η·dx dy`, `shallowwater.c:349-361`) the shift was
  algorithmically significant — a dry cell's stage was pulled toward its
  neighbors' *absolute* η, so wetting depended on the frame, and b1's
  outlet row never wets without the lift. P3 replaced that row with the
  translation-invariant zero-depth continuity closure (amendment A13; the
  correctness-fix table below) — b1 re-gates green with unchanged
  tolerances either way — and the offset remains as the shared elevation
  frame convention of both modules (`SurfaceSolver::elevationOffset`,
  `TerrainMetric::elevationOffset`).
- **The stored velocity applies the drag factor twice.** `update_velocity`
  multiplies the full expression by D although E already carries one factor
  (`shallowwater.c:745-746`, the "ignore drag inversion" comment). The
  velocities that transport mass in the continuity equation are the
  matrix-consistent ones; the *stored* uu/vv systematically understate the
  conserved face fluxes. Consequence: discharge measured as `uu·depth`
  underestimates the mass-consistent outflow (b4 measures the hydrograph
  from the volume budget instead; `tests/regression/tolerances/`).
- **Boundary faces see the pre-source ghost stage.** Physical-edge ghost η
  is refreshed only inside the free-surface phase; `update_velocity` runs
  after rain/evaporation without a refresh (legacy call order,
  `solve.c:51-98`), so rain induces a small outward boundary-face velocity
  each step. A flat "closed" basin therefore leaks slowly at its east/north
  edges. Rank-invariant: interface halos are exchanged fresh, matching what
  a single rank computes at interior faces.
- **Rain skips a configured exclusion region.** Legacy hardcoded "no rain
  on the last global row" (`shallowwater.c:596`) — a b1-specific outlet-row
  rule. Frehg2 expresses it as `surface_water.rainfall.exclude` (b1's
  configuration sets it; uniform-rain cases omit it).
- **West/south closed edges carry no explicit flux, east/north edges do.**
  The RHS divergence reads the ghost-side E at west/south faces (zero) but
  the owned-cell E at east/north faces (legacy asymmetry in
  `shallowwater_rhs`); the implicit legs are folded on all four edges.

- **The point-implicit drag factor surge/stalls at drag numbers above ~1.**
  D = 1/(1 + ½ Δt C_D |u^n|/h) linearizes the implicit drag at the lagged
  speed; when the equilibrium drag number ½ Δt C_D |u_eq|/h exceeds ~1 the
  fixed-point iteration is a 2-cycle (alternating surge and stall), and the
  pair-averaged conveyance collapses far below the Manning balance
  (measured at the b5 outlet, n = 2 against a bed-level reservoir: ~60×
  under-conveyance, an artificial 0.92 m impoundment, catastrophic release
  when the forcing stops). This is the preserved legacy scheme — every
  gated regime is either low-drag-number (b1–b4) or configured out of the
  regime (b5's free-outfall strip, amendment A16). A self-consistent drag
  solve would remove the limitation; it is deferred because it would
  perturb the b1/b4 goldens for no gated benefit.

- **The below-bed clamp creates volume, now measured.** `enforceSurfBc`
  lifts η onto the bed with no compensating flux (`shallowwater.c:507-511`)
  — preserved, since removing it changes wetting. In b1–b4 regimes the
  created volume is negligible; b5's hard-drawn films make it ~1.6 % of
  the rain volume, so the step audit measures it (`SurfaceStepAudit::
  clampVolume`, the mass_audit `clamped` column, amendment A16) and budget
  identities include it: closure is exact bookkeeping with the legacy
  defect reported as data.

## Deviations from legacy (deliberate, with rationale)

- **One global free-surface solve.** Legacy solved each rank's block
  separately with the neighbors' previous-step η as Dirichlet data
  (`shallowwater.c:371-392`), which is decomposition-dependent by
  construction. Frehg2 assembles one distributed SPD system (plan §5.1);
  single-rank results are unchanged, and rank invariance (plan §8.2) becomes
  achievable.
- **Prescribed-stage columns are eliminated.** Legacy kept the neighbor→
  stage-cell coupling, making the matrix unsymmetric under CG. Frehg2 folds
  the known stage into the right-hand side; the linear solution is
  identical and the matrix stays SPD.
- **kind `outflow` (plan §5.6 amendment).** SERGHEI-style free outflow for
  b4: the ghost stage at member domain-edge faces is extrapolated down the
  continued bed slope, so the boundary face keeps the interior momentum
  balance. Decided by the P1 b4 gate exactly as the plan anticipated: the
  provisional stage-sink mapping keeps the outlet column dry, and a plain
  open boundary retains ~half the rain.
- **Wind direction interpolates linearly.** Legacy held the direction
  piecewise-constant at the previous sample while interpolating the speed
  (`solve.c:177-188`). No benchmark exercises wind; the difference is
  documented rather than reproduced.
- **Monitors and outputs**: HDF5 per plan §7 instead of ASCII; the mass
  audit table `/monitor/mass_audit` records the domain volume and
  cumulative rain/evaporation/outflow/boundary-inflow volumes every step.

## Correctness fixes (legacy defects of the plan §2.1 hazard class)

| Defect | Legacy location | Fix |
|---|---|---|
| uy/vx interpolation index arithmetic wraps to the far end of the row at i = 0 and reads out of bounds at j = 0 | `interp_velocity`, `shallowwater.c:969-972` (`icjP[ii]-1`, `iPjc[ii]-nx`) | intended four-point stencil with ghost/corner values (`Momentum.cpp`); identical on b1 (nx = 1 reduces both to the same average); b4's j = 0 row previously read heap garbage |
| discharge (inflow) divides by the per-rank member count — a region spanning R ranks injects R× the flow | `shallowwater.c:313`, `get_BC_location` | divide by the region's global member count (`FreeSurface.cpp::assembleRhs`) |
| north-edge stage-velocity correction rank test uses `>` for `>=`, silently skipping the first top-row rank (including single-rank runs) | `enforce_velo_bc`, `shallowwater.c:952` | corrected global-edge test (`WetDry.cpp::enforceVeloBc`) |
| ghost flow rates Fu/Fv read uninitialized memory (Fu is never initialized; `initialize.c:649` zeroes Fv twice) | `init_Data`/`ic_surface` | all fields zero-initialized; physical-edge ghost Fu/Fv are defined as 0 |
| dry-cell rows hold η with live off-diagonal legs whenever a face velocity is nonzero and the RHS overwritten to η·dx·dy: the held row reads η = η_old + (S/A)·η_neighbor — a frame-dependent pull on the neighbors' *absolute* stage that mints volume out of the datum (legacy's accidental wetting mechanism) and, in the velocity-zero case, leaves the wet neighbor's leg unmatched (asymmetric matrix, CG stalls). ~1e-4 of a cell area on b1/b4; fatal at b5's g·dt² scale (found at P3, amendment A13) | `shallowwater_mat_coeff`, `shallowwater.c:349-361` | dry cells close as zero-depth continuity rows: free-surface area = cell area (∂V/∂η at the bed), legs live, standard RHS. The higher-of-two-bottoms face depths vanish whenever the neighbors' stage sits below the dry cell's bed, so live legs only carry flux *into* the dry cell — implicit, conservative, symmetric, translation-invariant wetting (`FreeSurface.cpp::assembleCoefficients`). b1/b4 re-gated green with unchanged tolerances (b1 worst achieved/allowed 0.231, the P1 value) |
| drag coefficient written only where Vs > 0: persistently dry cells keep the initial zero forever, and a face between a dry cell and a wet neighbor carries water with *no* drag — the point-implicit update accumulates gravity unopposed (b5's dry outlet-adjacent faces reached km/s; found at P3, amendment A13) | `update_drag_coef`, `shallowwater.c:984` | dry cells evaluate the same drag law at their deepest adjacent face depth (the water actually flowing across them: an outfall reservoir sees the upstream depth, a hairline film face sees the film), floored at min_depth; wet-cell values are bit-identical (`Momentum.cpp::updateDragCoef`) |
| the volume audit counts only the east/north domain edges (the west/south ghost-side E is identically zero under the preserved legacy asymmetry, so nothing is missed today — but any later boundary kind that fills the ghost momentum would leak unaudited) | `Frehg2 P1` (legacy had no audit) | all four domain edges counted (`FreeSurface.cpp::accumulateBoundaryFluxes`) |

Dead legacy code found during the port is recorded in
[removed-features.md](removed-features.md).
