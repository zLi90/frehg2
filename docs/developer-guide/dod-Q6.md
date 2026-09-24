# Definition of Done — Q6 (Wind stress)

Per v2 development plan §5 and §9, under the v1 plan's §11.4 DoD
conventions. Every item names its verification; commands were run green
on 2026-09-24 on macOS arm64 (gcc-16.2, MPICH 4.3, PETSc 3.25.1, Kokkos
5.1.1 OpenMP host backend; `FI_PROVIDER=tcp`).

## Gates (blocking, §9 Q6 row)

- [x] **g9 — steady wind setup** (per-PR, `regression.g9`): PASS, all
      three variants. (a) the TELEMAC case-10 replica (nonlinear,
      0.815 m setup in 2 m of water) against the closed form
      h² = C + 2τx/(ρg) with volume-conserving C: 0.22 % (bound 0.5 %)
      at the §6.3 convergence-study resolution dx = 2.5 m — **the first
      exercise ever of the v1 wind path (authored P1, ungated since):
      it is correct**, with first-order discretization convergence
      (1.70/0.87/0.44/0.22 % at dx = 20/10/5/2.5). (b) linear regime,
      Garratt law at U₁₀ = 15 m/s through the (u10, v10) component form,
      Cd recomputed independently in python: 0.17 % (bound 1 %). (c)
      sloping bottom (−0.12 → −2 m; h_min ≈ 0.096 m grazing the wet/dry
      machinery) vs the RK4+bisection quadrature: profile RMS 0.40 % of
      the setup range (bound 1 %; max deviation recorded — it
      concentrates in the shallow-tip cell and converges first order),
      shallow-tip thin-layer drag limit cycle recorded.
- [x] **g10 — seiche relaxation** (per-PR, `regression.g10`): PASS.
      Fundamental period from zero crossings 1010.3 s vs Merian
      1009.7 s (0.06 %, bound 2 %) at Cr = 0.50, released by a
      forcing-series step to zero (also gating the series edge
      handling); frictionless amplitude decay per period recorded
      (θ-scheme dissipation). Basis: the g9(b) basin per V2-A16 (the
      g9(a) state is a 40 %-of-depth nonlinear seiche where Merian does
      not apply at 2 %).
- [x] **Cd unit battery** (per-PR, in `unit.all`): hand tables for
      every law at U₁₀ ∈ {0, 5, 10, 20, 30, 40, 60} m/s to 1e-12
      relative, including caps and the Large & Pond 11 m/s breakpoint
      (the published 5e-6 jump kept faithful), the constant law's
      no-cap legacy behavior, the direction-wrap interpolation
      (350°→10° through 0°, never 180°), and component-vs-compass
      equivalence.
- [x] **x-battery, 8-orientation setup case**
      (`regression.wind_orient`, per-PR): the same square-basin setup
      under 8 compass directions maps onto itself under the matching
      rotations/reflections to ≤ 1.1e-15 m (tolerance 1e-7) — the wind
      path is exactly symmetric.
- [x] **Gate-first (§6.1):** authored red in commit `45e4255` ((b)/(c),
      g10, orientation "capability absent"; (a) green — it gates
      pre-existing behavior); §6.3 battery `unit.g910_gates.negative`
      (26 expectations) green throughout, each closed form checked two
      independent ways.
- [x] **b1–b6 unchanged (wind off):** full per-PR tier green
      post-capability. The momentum refactor moved the per-step wind
      evaluation into `WindForcing` (beginStep) with identical
      arithmetic on the legacy path; wind stays off in every b-case.
- [x] **s-matrix:** the wind sample is three per-step host scalars into
      the existing momentum kernel — no new parallel structure; the
      standing s-gates ride the nightly runner unchanged.

## Deliverables (§5.2) — state

- [x] Selectable Cd(U₁₀): `constant` (uncapped legacy Cw, default),
      `garratt`, `smith-banke`, `wu`, `large-pond`, capped at
      `wind.cap` (default 3.5e-3, the ADCIRC practice). Schema
      cross-rules: `cd` belongs to the constant law, `cap` to the U₁₀
      laws (and the resolve emitter mirrors them — the r1 round-trip).
- [x] Component wind forcing: `wind.u10`/`wind.v10` (grid frame,
      wrap-free, `north_angle` not applied) exclusive with the legacy
      compass `speed`/`direction` pair.
- [x] Direction-series wrap convention: sampled directions become unit
      vectors at load, interpolate linearly, angle recovered by atan2 —
      shortest-arc (chord) interpolation, documented (deviation from
      constant rate O(Δ²/8); a 180° step between samples is degenerate)
      and unit-tested across the 360°→0° wrap. The legacy
      constant-direction arithmetic (truncated π literal) preserved
      bitwise.
- [x] Thin-layer attenuation unchanged (legacy-faithful; both gate
      basins deep relative to hD except the g9(c) tip, which stays
      above its configured hD).
- [x] `swe/WindForcing.hpp`: host-side per-step evaluation, unit-tested
      standalone; the momentum kernel consumes three scalars.

## Amendments this phase

V2-A16 (g10 basis moved to g9(b); ramp+friction spin-up made concrete;
the g9(a) and g9(c) convergence studies archived).

## Known limitations carried forward

1. Wind remains spatially uniform per step (per-cell wind fields are a
   future extension; the plan's Q6 scope is per-component *series*).
2. The g9(c) shallow-tip max deviation (first-order, recorded) and the
   thin-layer drag limit cycle there (P3 finding, recorded) are scheme
   properties, not wind defects.
3. Rejected-for-gating references (Okeechobee/Erie hindcasts, Csanady
   gyres) stay rejected per §5.3; Kraus & Militello remains the
   optional time-varying-wind stretch case.
