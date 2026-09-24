# Definition of Done — Q5 (Temperature transport, surface + subsurface)

Per v2 development plan §4 and §9, under the v1 plan's §11.4 DoD
conventions. Every item names its verification; commands were run green
on 2026-09-24 on macOS arm64 (gcc-16.2, MPICH 4.3, PETSc 3.25.1, Kokkos
5.1.1 OpenMP host backend; `FI_PROVIDER=tcp`, `OMP_NUM_THREADS=1` for
serial gate runs). Realization decisions: plan amendments V2-A17
(design) and V2-A18 (the orientation-battery catch).

## Gates (blocking, §9 Q5 row)

- [x] **g6 — subsurface heat, analytic** (per-PR, `regression.g6`):
      PASS, all four parts, each at the MEASURED Darcy flux (V2-A17;
      window-asserted, landed 1e-17–1e-19 from target).
      (a) B&P steady Péclet sweep at measured Pe ∈ {−5, −1, ~0, +1, +5}:
      worst max-norm 0.087 K of the 0.1 K bound (|Pe| = 5 — the exact
      discrete-operator residual ln(1+Pe_c)/Pe_c; Pe = 0 measures
      3.4e-10 K). t_end re-derived at bring-up: Pe ≠ 0 excites the n = 1
      advective-diffusive mode (τ₁ = L²/(α(π² + Pe²/4)) = 1.36e7 s at
      |Pe| = 1) that the Pe = 0 odd symmetry annihilates — 8e7 s ≥
      5.9 τ₁ everywhere (the 2e7 s first cut left a 0.83 K transient).
      (b) Ogata–Banks heat step, OGS parameter set, vertical graded
      column (dz₀ = 0.125, stretch 1.02, nz = 145): superbee
      0.04–0.38 % of the step at ~1/5/10/20/50 m (bound 1 %) for
      t ≥ 2e6 s; upwind recorded 0.34–4.9 %, not gated. The gating
      window + mesh bound the pinned-inlet-cell half-thickness
      ambiguity (errors halved with dz₀ — representation, not scheme).
      (c) Stallman diel damping, q_z ∈ {0, ±2e-6 m/s}: amplitude ratios
      within 0.04–3.4 % (bound 5 %), phase lags within 11–339 s (bound
      600 s) at z = 0.08/0.18/0.28 m relative to the pinned top cell's
      own fit; the ±q amplitude separation (~30 % at the deepest gated
      depth) carries the direction discrimination.
      (c, coupled) the composed variant (plan §4.2): eta-pinned pond,
      diel wet-cell temperature reset, heat entering through the
      coupled seepage + conduction exchange — 0.32–0.89 % amplitude,
      85–331 s phase, indistinguishable from the prescribed-BC physics.
- [x] **g7 — surface heat exchange, analytic** (per-PR,
      `regression.g7`): PASS. (a) Edinger relaxation: max
      |T − exact| = 0.0079 K of the 0.05 K bound — the explicit-Euler
      bias exactly; per-step heat-ledger closure 3.8e-16 (bound 1e-8)
      on `/monitor/temperature_audit`. (a2) full bulk mode settles at
      the offline root T_e = 20.550 °C to < 1e-3 °C (bound 0.05;
      trailing span 1.1e-6 °C) — the python chain mirrors
      `atm::netHeatFlux` term for term, with hand-computed spot values
      pinned in the §6.3 battery. (b) channel plume: measured
      u = 0.5000 (bound ±1 %), steady L2 = 0.01 % (bound 1 %),
      mid-channel transient 1.75 % of the step (bound 2 %) against the
      closed-form step response convolved with the measured inlet
      history. §6.3 rule-3 study in the case comments: dx = 500
      (front < 2 cells) measured 6.2 %, dx = 100 (~9 cells) gates.
- [x] **g8 — thermal convection onset** (nightly, recommended;
      `regression.g8`): PASS. Ra_eff = K β_T ΔT H_eff/α_e with
      H_eff = H − dz (pinned-cell centres, V2-A17): 28.5 decays the
      seeded critical mode to 2.0e-6 of seed (bound 0.1) with trailing
      Nu = 1.0000; 52.25 convects at trailing Nu_min = 1.53 (bound
      1.05), mode saturated at 2.8 K. Seed = the admissible m = 2
      cosine mode (the drafted sine seed was orthogonal to it — §6.3
      catch, orthogonality now a battery expectation; generator
      committed).
- [x] **x-battery for heat** (per-PR): `regression.heat_orient` — the
      8-dihedral battery on the 2D thermal case (side Dirichlet
      visiting all four sides, thermal density coupling active,
      strict-mode solver): temperature ≤ 1.5e-13, head ≤ 1.8e-15
      against the §8.1 1e-12 tolerance. **First x-battery execution to
      catch a defect** (V2-A18): the legacy y+-only limiter side-ghost
      admission clipped side thermal Dirichlets on the other three
      sides (compounded by the documented corner-face spill) —
      temperature now admits on all four sides; salinity keeps the
      b6-pinned y+ rule byte-for-byte.
- [x] **s-matrix for the new kernels**:
      `regression.heat_rank_invariance.{strict,default}` at 1/2/4
      ranks: 4.3e-15 (bound 1e-12) / 4.5e-8 (bound 1e-4);
      `regression.heat_restart` two-scalar restart determinism ≤ 1e-12
      (exact-class). Risk-register item "restart lanes extended to the
      two-scalar checkpoint layout" — the `t_*` fields mirror `s_*`, a
      salinity-only checkpoint stays byte-identical to v1.
- [x] **Gate-first (§6.1):** authored red in `b67f9aa` — all six
      regression entries failing "capability absent" against the
      pre-Q5 schema; `unit.g678_gates.negative` (49 expectations, every
      closed form two independent ways, every criterion demonstrably
      able to fail, the dihedral array/coordinate ops proven
      consistent) green throughout.

## Also required (§9 Q5 row)

- [x] **b6 unchanged with temperature off** — the ScalarSpec refactor
      is spec-keyed with exact-neutral arithmetic on the salinity path:
      per-PR tier (b1–b4 exact-bound gates, b1/b2/b6 restart
      determinism at 1e-12, b6-gw smoke, rank-invariance lanes) green
      post-capability; full-horizon nightly b6 ss + td PASS (this
      machine, close day).
- [x] **g5 cross-check (§4.3)** — `scripts/cross_check_g5_temperature.py`:
      g5 rerun with the full temperature module active (transport of
      the uniform 20 °C field, thermal conduction/retardation, β_T =
      2e-4 armed with (T − T₀) ≡ 0) against the Q4-mode baseline:
      **bitwise identical** on water_content, hydraulic_head, and
      concentration at every output (max |Δ| = 0.0), temperature holds
      20 °C to rounding — stronger than the plan's "statistically
      identical".
- [x] **Density-coefficient config migration is golden-neutral** —
      `beta_saline`/`beta_saline_viscosity` defaults equal the legacy
      compile-time constants (same literals → same doubles); with no
      temperature attached the r_ρ arithmetic is unchanged
      (`1 + β_s·s − β_T·0` never evaluated; guarded reads). Verified by
      the b-tier above; b6 runs the coupling live.
- [x] **p1/p5 for the new kernels (V2-A6)** — `check_forbidden.sh`
      clean (no backend branches, MemSpace-view idiom; one
      nvcc-hazardous lambda-capturing-lambda flattened before the CUDA
      lane); the kernels add no PETSc interaction, and the aij solver
      lanes of the per-PR tier are green.

## Sweeps at close (2026-09-24)

- [x] Full per-PR tier: **68/69**, the one failure being
      `regression.perf_baseline` under the concurrent `-j4` tier — the
      documented A23 M3 timing lottery (it failed identically in the
      PRE-Q5 baseline tier and passes standalone before and after:
      118–124 s). Every gate, restart, rank-invariance, validate, unit,
      and mpi entry green.
- [x] Sanitizers (ASan+UBSan build): g7a (equilibrium exchange), g7a2
      (bulk chain), g6c short (pins + series), g6c-coupled short
      (seepage heat), orient-base (side Dirichlet + thermal density),
      and the 179-test unit suite — zero reports.
- [x] Docs: `docs/theory/temperature.md` (nav'd), output.md
      (temperature group + `temperature_audit` columns), parameters.md
      (16 new keys; `unit.parameter_docs` lockstep green), case
      READMEs with the §6.3 studies; `mkdocs build --strict` clean.
- [x] r1 lockstep: every new schema key extends `resolvedConfigYaml`
      mirroring the module cross-rules; every gate run passed the
      harness's `check_run_record` round-trip hook.

## Open items leaving Q5

- §8.1/§8.2 x-gate backfill (pre-existing; V2-A18 seeds the exemption
  table: salinity y+ admission, corner-face spill), V2-A11 west/south
  `outflow` first x-gate deliverable, owner p6 GPU bundle, owner §6.4
  digitization sign-off. Remaining phase: **Q7 (release)**.
