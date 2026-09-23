# g5 — Geng & Boufadel (2015) bare-soil salinization (v2 Q4, plan §3.3/§3.4)

Nightly-class code-to-code gate against MARUN (Geng & Boufadel, *J.
Hydrology* 524:427-438, 2015): bulk-aerodynamic evaporation with humidity
feedback, α₁ moisture limiting, evaporative salt concentration under the
zero-total-flux (Cauchy) top condition, and density feedback. 50 × 2 m
x-z slice, Table-1 parameters, 50 h.

- `g5-geng2015.yaml` — the gated run (density coupling on).
- `g5-geng2015-nodensity.yaml` — the β = 0 control for g5(iv) (V2-A14:
  at 50 h the β run's 30 g/L plume edge must reach ≥ 0.05 m deeper).
- `reference/` — digitized MARUN curves + `DIGITIZATION.md` (the §6.4
  protocol record) + the overlay plots for the owner's visual check.

Harness: `tests/regression/run_regression.py g5` (criteria per V2-A13,
V2-A14 and V2-A15; tolerances in `tolerances/g5-geng2015.yaml`; §6.3
negative battery in `scripts/test_g45_gates.py`, which also asserts the
digitized reference passes its own criteria).

**V2-A15 (read it before touching the criteria):** the paper's figures
are mutually inconsistent — Fig. 4/9's profiles hold 5-8x more moisture
deficit than Fig. 3's evaporation supplies, with interval increments in
the constant ratio ~4.28 = R_air(1)/R_air(5): the profile figures
evidently come from a ~4.3x harder-forced run than Table 1 states. The
gate therefore anchors on the self-consistent subset (the Table-1 closed
form, Fig. 3's t->0 value, monotone decay, the salinization peak, the
salt-mass bound, and the V2-A14 density inequality) and computes the
Fig. 4/9b profile RMS as a *recorded* comparison, not a criterion. The
E(t) decay-rate mismatch (~3x at 10 h) is the internodal-conductivity
difference between the legacy-pinned upstream-K faces and MARUN's
starving surface node — insensitive to mesh (2 mm top cells move it 4 %)
and out of scope to "fix" (b2/b3 goldens pin the scheme).

Observable contract the capability must satisfy (checked by the gate):
`/monitor/gw_mass_audit` gains a cumulative `evaporation` column
(appended, index 7), from which the gate reads the actual evaporation-rate
history over the x ∈ [1, 49] m zone.
