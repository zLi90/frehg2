# pd — periodic-topology gates (Stage-3 plan T3)

Per-PR gates for `domain.periodic: {x, y}` (Stage-3 plan §7/§13.B; owner
decision O3 — axis-generic, both axes gated). Harness subcommands `pd1`,
`pd2`, `pd3`, `pd4`, `pd5`, `pd6`, `pd8` in
`tests/regression/run_regression.py`; pure references in
`tests/regression/gate_periodic.py`; tolerances in
`tolerances/stage3-pd.yaml`; §6.3 battery in `scripts/test_pd_gates.py`.

| Gate | Case | Reference | Asserts |
|---|---|---|---|
| pd1 | `pd1-alongshore.yaml` (+ the harness-generated transposed lane: periodic y, v10) | the **discrete wind–friction balance** recomputed in-script (Garratt Cd independent, relative-wind stress, the point-implicit double drag factor), bisection-solved; authoring self-check: the fixed point closes to 1e-12 and the dt→0 limit recovers the independent continuum closed form (0.07446 vs 0.07455 m/s at dt = 2 s) | strict lane: η spread, u spread/value vs the balance, trailing-window steadiness; default-solver lane: the bjacobi-ICC translation asymmetry recorded (the T10c 1e-6 class), not gated |
| pd2 | `pd2-surface.yaml` (periodic bed bump + angled wind + tracer slab ON the seam), `pd2-gw.yaml` (exactly-periodic quantized K, side water tables + recharge), `pd2-coupled.yaml` (saturated tidal basin, external closure, consistent salinity both ways, Ss-driven seepage exchange — see the S3-A6 design notes in the file header) | the exact translation identity: run(roll(inputs, k)) == roll(run(inputs), k) for k ∈ {3, 6}, strict solver | per-field-class bounds (S3-A6): gw lane `roll_abs` 1e-11; surface/coupled flow `roll_surface_abs`, tracers `roll_tracer_abs` (the measured branch-amplification classes); anti-vacuity: the tracer crosses the seam (pd2-surface) and the tidal exchange moves salt into the soil (`coupled_exchange_floor`, pd2-coupled) |
| pd3 | pd2 bases at `mpi_nx` = 1 (self-wrap), 2 (same rank both sides), 4; the coupled lane adds 2 × 2 | the n = 1 run | per-field-class rank invariance (S3-A6): gw 1e-12 local (V2-A9 1e-9 on CI), surface/coupled at the pd2 field classes; `FREHG_EXPECT_RANKS` guards the world size |
| pd4 | `pd4-tidal-periodic.yaml` (ob3 made periodic-x, closure external) + unit `SweModule.PeriodicSeamCarriesNoBoundaryFlux` | ledger identity | the surface ledger closes with ONLY the y⁺ stage exchange booked; the unit test asserts `boundaryOutflow == 0.0` exactly on a periodic channel with uniform flow |
| pd5 | `pd5-confined.yaml`, lanes uniform-K and periodic column-varying K | confined Darcy: head linear in y, uniform in x/z; qx ≡ 0 | `check_linear_head` (x/z spread, vanishing second differences, drop direction), `check_zero_x_flux` vs the qy scale; the K-varying lane must match the uniform lane (no x gradient exists) |
| pd6 | schema/setup negatives | — | nx < 3 rejected; an x-edge-only BC under periodic x fatal at setup with the §13.B message; the pd2 seam warning fires on a deliberately non-periodic bump lane; y-periodic + y-edge eta rejected |
| pd7 | `regression.s3.default_bitwise` | T0 hashes | `periodic` absent ⇒ bitwise |
| pd8 | two batteries with the W3.7 flag swap (rot90 maps periodic-x onto periodic-y): the `rl-realloc` vadose square (gw lane, x-orient gw classes) and `pd8-orient-base.yaml` (periodic torus bump, oblique wind rotated with the domain, tracer slab, y walls; scalar fields — eta/depth/tracer — per the swe-orient convention) | run identity | gw: x-orient strict classes; surface: pd2 field classes on `transpose`, the x-orient swe mixed class on the six edge-swapping transforms (`stage3-pd.yaml` pd8) |
| restart | `s3-restart` lane pd2-surface (checkpoint mid-run) | the project restart standard (S3-A6) | flow fields ≤ 1e-12 ABSOLUTE (measured ulp class 1.9e-15); tracer at its measured branch class (`tracer_rel` 1e-4, measured 4.0e-6) — the walled control of the same case measures LARGER deviations on every field |

Steadiness is asserted in-gate (pd1); every check is provably able to fail
(`unit.pd_gates.negative`), including a 1e-9 seam glitch, a seam booked as
boundary flux, and a 1e-12 m/s wrapped-leg x-flux leak.

**Red evidence (gate-first, plan §2.3):** `domain: unknown key 'periodic'`
on the T2 binary for all five cases — `stage3/records/T3/red/`. The
`pd2-coupled.yaml` and `pd8-orient-base.yaml` lanes were authored after the
capability landed (S3-A6); their red evidence is the same schema rejection
(the `domain.periodic` key is what every pd case fails on pre-W3.1).
