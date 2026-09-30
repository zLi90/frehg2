# Atmosphere module (frehg::atm)

The `atmosphere` configuration block and `src/atm/` supply one
bulk-aerodynamic (mass-transfer) chain with several consumers: open-water
evaporation and soil evaporation (v2 Q4, plan §3.2), and the surface
latent/sensible heat terms (v2 Q5, plan §4.1). This chapter consolidates
the shared physics; the consumer-side behavior lives with its module —
[surface water § evaporation modes](surface-water.md#evaporation-modes-v2-q4-plan-32)
and [temperature § surface exchange](temperature.md).

**Provenance.** Geng & Boufadel (2015) Eqs. (1)–(6) — the formulation
chain whose Table-1 values the V2-A13 derivation reproduces exactly
(q_a = 2.896e-3 vs the stated 2.9e-3; E(0) = 1.470e-7 m/s vs the digitized
1.442e-7). The aerodynamic resistance is the Liu et al. (2006) fit named
by the plan. The legacy aerodynamic model hardwired these constants (v1
removed it for that, not for its physics); in v2 every input is
configuration. All formulas are pure host/device functions
(`BulkAerodynamic.hpp`); the per-step evaluation of the configured series
is host-side (`MetForcing.hpp`), so kernels consume sampled scalars.

## The chain

Saturated vapor pressure (Tetens, [kPa], T in °C):

$$e_{sat}(T) = 0.6108\, \exp\!\left(\frac{17.27\,T}{T + 237.3}\right)$$

Saturated specific humidity at pressure \(P\) [kPa]:

$$q_{sat}(T, P) = \frac{0.622\, e_{sat}}{P - 0.376\, e_{sat}}$$

Aerodynamic resistance (Liu et al. 2006 fit, U [m/s]):

$$R_{air}(U) = 94.909\, U^{-0.9036}$$

Dry-air density from the ideal gas law:
\(\rho_a(T, P) = 1000\,P / (R_d\,(T + 273.15))\), \(R_d = 287.05\).

Evaporative water flux (Mahfouf & Noilhan form; [m/s] of liquid water,
positive = evaporation):

$$E = \frac{\rho_a}{R_{air}\,\rho_w}\,\big(q_g - q_a\big),
  \qquad q_g = \alpha_1\, q_{sat}(T_s, P)$$

- **Open water**: \(\alpha_1 = 1\) (the potential rate).
- **Soil** (`groundwater.evaporation`): \(\alpha_1(w_g) = \min(1,\;
  1.8\,w_g/(w_g + 0.30))\) (Barton / Lee & Pielke, paper Eq. 6) at the
  surface-cell volumetric water content.
- **Negative results (condensation) are returned as-is** — clamping to
  zero would hide a humidity-gradient sign error from the mass audits
  (V2-A13 decision). The surface consumer deposits the condensate.

## Q5 surface net heat flux

Positive into the water (evaluated per wet surface cell):

$$Q_{net} = Q_{sw} + \varepsilon\,\big(LW_{in} - \sigma T_K^4\big)
  - Q_{lat} - Q_{sens}$$

with \(\varepsilon = 0.97\) (Kirchhoff both ways: the same emissivity
absorbs incoming longwave and emits), \(Q_{lat} = \rho_w\, L_v(T)\, E\)
on the Q4 `evaporationRate` chain **verbatim** — including its air density
evaluated at the *water* temperature (the Q4 convention; do not "fix" it,
the g7(a2) offline root-finder mirrors it term for term) —
\(L_v(T) = 2.501\times 10^6 - 2370\,T\) [J/kg], and
\(Q_{sens} = \rho_a(T_{air})\, c_{pa}\, (T - T_{air}) / R_{air}\),
\(c_{pa} = 1005\).

Both bulk terms evaluate at the **local water temperature**, never at
`atmosphere.surface_temperature` (V2-A17) — that key exists only for the
Q4 evaporation consumers (it is required exactly when one is configured).

## Forcing block conventions

- Every input (`air_temperature`, `pressure`, `wind_speed`,
  `shortwave`/`longwave_in` for heat, and exactly one humidity form) is a
  constant or a time series, linearly interpolated at the step start.
- `relative_humidity` is folded at sample time:
  \(q_a = RH \cdot q_{sat}(T_{air}, P)\) — consumers only ever see a
  specific humidity.
- `atmosphere.wind_speed` is floored at `atmosphere.wind_speed_floor`
  (default 0.5 m/s, the GLM still-air convention) for **all** bulk
  consumers: \(R_{air}(U \to 0) \to \infty\) is a removable singularity
  the floor regularizes.
- `atmosphere.wind_speed` (scalar, bulk transfer) is deliberately
  separate from `surface_water.wind` (vector, momentum): the Q5 review
  decided against unifying a scalar and a vector schema surface. Point
  both at the same series when physical consistency matters.

## Validation

The chain is gated by g4 (analytic drawdown + evaporative concentration,
`regression.g4`), g5 (Geng & Boufadel 2015 laboratory salinization,
`regression.g5.geng2015` and its transposed slice), and g7
(equilibrium relaxation + the measured-inlet-convolved channel plume,
`regression.g7`); unit hand-value tables cover every formula above
(`unit.all`, suites `BulkAerodynamic` and `MetForcing`).
