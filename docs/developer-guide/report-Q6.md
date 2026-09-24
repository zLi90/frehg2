# Q6 Report — Wind stress

**Phase:** v2 Q6 (v2 plan §5)
**Completed:** 2026-09-24
**Machine:** macOS arm64 (M3), gcc-16.2, MPICH 4.3, PETSc 3.25.1,
Kokkos 5.1.1 OpenMP host backend.

## What the phase delivered

The v1 wind path — implemented in P1 and never exercised by any gate
(the completion report's own hard lesson, and the reason §8.3 exists) —
is now validated and generalized. Selectable Cd(U₁₀) laws (Garratt,
Smith & Banke, Wu, Large & Pond, each capped; the uncapped legacy
constant Cw stays the default), grid-frame (u10, v10) component forcing
alongside the legacy compass pair, and a documented, unit-tested
shortest-arc convention for direction-series interpolation across the
360°→0° wrap (the legacy linear-in-degrees interpolation took the long
way round; the path had never run). All per-step evaluation lives in
`swe/WindForcing.hpp`, unit-testable standalone; the momentum kernel
consumes three scalars and is otherwise untouched.

## Gate results

| gate | measured | bound |
|---|---|---|
| g9(a) nonlinear setup vs closed form (TELEMAC case-10 replica) | 0.22 % | 0.5 % |
| g9(b) linear setup, Garratt law via components, Cd recomputed offline | 0.17 % | 1 % |
| g9(c) sloping-bottom profile RMS vs quadrature (h_min ≈ 0.096 m) | 0.40 % of range | 1 % |
| g10 seiche period vs Merian (Cr = 0.50) | 0.06 % | 2 % |
| 8-orientation symmetry battery | ≤ 1.1e-15 m | 1e-7 m |
| Cd law tables / wrap / component equivalence | exact to 1e-12 | unit battery |

## Findings

1. **The v1 wind path is correct.** Its first exercise ever (g9(a))
   converges to the analytic setup at first order in dx
   (1.70/0.87/0.44/0.22 % at dx = 20/10/5/2.5 m) — the bias is the
   face-depth evaluation in the stress term, a discretization property,
   vanishing under refinement. Five phases after it was written, the
   authored-unexercised debt is paid.
2. **Steady wind setup is the ideal spin-up-friendly gate**: the steady
   state has u = 0 identically, so mild Manning friction damps the
   transient without touching what is gated, and a quasi-static wind
   ramp avoids exciting the seiche at all (V2-A16 makes the plan's
   "mild relaxation" concrete). Steadiness is asserted in-gate.
3. **g10 had to move basins (V2-A16):** releasing the g9(a) state gives
   a 40 %-of-depth nonlinear seiche where the linear Merian period does
   not apply at 2 %; the g9(b) release (0.25 % of depth) measures
   1010.3 s vs 1009.7 — 0.06 %.
4. The g9(c) shallow tip exposes two known scheme properties, both
   recorded rather than gated: the first-order shallow-cell deviation
   (max-norm 2.4 % at the gated dx = 5 while the RMS is 0.40 %) and the
   P3 thin-layer drag limit cycle (~5e-4 m span at h ≈ 0.1 m).

## Cost

g9 ≈ 20 s, g10 ≈ 4 s, orientation battery ≈ 5 s — all per-PR. The gates
add three fast tests to the standing per-PR tier.

## Carried forward

Per-cell wind fields (future extension); the rejected-for-gating
hindcast references stay rejected (§5.3); Kraus & Militello (1999)
remains the optional time-varying-wind stretch case.
