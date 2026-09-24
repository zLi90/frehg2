# g9/g10 — wind setup and seiche gates (v2 Q6, plan §5.3)

Per-PR gates, closed forms computed in-script (`tests/regression/
gate_wind.py`; harness subcommands `g9`, `g10`, `wind-orient`;
tolerances in `tolerances/g9-wind.yaml`; §6.3 negative battery in
`scripts/test_g910_gates.py`).

| Case | Reference | Notes |
|---|---|---|
| `g9a-setup` | h(x)² = C + 2(τ/ρw)x/g, volume-conserving C (TELEMAC-2D validation case 10 replica: 500 × 100 m, 2 m, τ/ρw = 1.2615e-3·U², U = 5 m/s) | constant-Cd law — **the v1 wind path's first gate ever**; end-to-end setup within 0.5 % |
| `g9b-linear` | same closed form in the linear regime (0.25 % of depth) | Garratt Cd(U₁₀) at 15 m/s through the (u10, v10) component form; Cd recomputed independently in python |
| `g9c-slope` | RK4 + bisection quadrature of dη/dx = (τ/ρw)/(g(η−b(x))) | Garratt at 12 m/s; bed −0.12 → −2.0 m; h_min ≈ 0.096 m grazes the wet/dry machinery without crossing; RMS-gated (see below) |
| `g10-seiche` | Merian T = 2L/√(gH) = 1009.7 s | wind series ramps, holds, steps to zero; frictionless; period from zero crossings within 2 % at Cr = 0.50; decay/period recorded |
| `orient-base` | the run set itself | the §8.1/§9 x-battery row: 8 compass orientations of the same square-basin setup must map onto each other under the matching rotation/reflection |

**Spin-up ("mild relaxation", the plan's own words):** the wind ramps
quasi-statically (over ≥ 20 seiche periods) and the setup cases carry a
mild Manning n — legitimate because the gated steady state has u = 0, so
friction cannot alter it; it only damps the transient. Steadiness is
asserted in-gate (trailing-window span), so a still-sloshing basin
cannot gate. g10 is frictionless (its decay is θ-scheme dissipation,
recorded not gated) and rides the g9(b) basin, not g9(a) — V2-A16: the
(a) state is a 0.8 m setup in 2 m of water, a 40 %-of-depth nonlinear
seiche where the linear Merian period does not apply at 2 %.

**§6.3 convergence study (g9c, archived here):** the deviation from the
quadrature concentrates in the shallow-tip cell and converges first
order — max deviation 5.5 / 2.4 / 1.1 % of the setup range at
dx = 10 / 5 / 2.5 m, profile RMS 0.40 / 0.15 % at dx = 5 / 2.5. The gate
therefore runs dx = 5 with the **RMS** as the gated 1 % norm (2.5×
headroom) and the max deviation recorded. Steadiness gates at the deep
end: the shallow tip carries the P3-documented thin-layer drag limit
cycle (~5e-4 m span), legacy scheme behavior recorded, not gated.

**§6.3 convergence study (g9a, archived here):** the discrete steady
setup converges to the closed form at first order in dx — rel. error
1.70 % / 0.87 % / 0.44 % / 0.22 % at dx = 20 / 10 / 5 / 2.5 m (the bias
is the face-depth evaluation in the momentum stress term, not a wind
defect: it shrinks uniformly under refinement). The gate runs dx = 2.5,
where the 0.5 % bound holds with 2.3× headroom.

Closed basins: every case closes its edges with the explicit
`kind: velocity` value-0 condition (the g4 bring-up finding — default
edges are legacy-transmissive).
