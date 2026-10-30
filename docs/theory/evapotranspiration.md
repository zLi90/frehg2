# Evapotranspiration (Stage-3 T5)

The `evapotranspiration` block replaces both legacy evaporation keys
(`surface_water.evaporation`, `groundwater.evaporation`; the schema
rejects mixing the forms — the legacy paths themselves are untouched and
stay golden-pinned). One potential mode feeds two channels:

- the **open-water channel** removes `min(E_w·dt, depth)` from every wet
  surface cell per step (a negative potential condenses onto water);
- the **soil channel** removes the capped soil flux `qVap` through each
  exposed column's top face per groundwater substep.

## Potential modes and aerodynamic resistance

`potential: bulk` is the Geng & Boufadel (2015) bulk-aerodynamic chain —
at `r_a = liu2006`, `r_s = 0` it is the v2.0 `atm::evaporationRate`
expression tree term for term (`src/atm/EtFormulas.hpp`; the et0 unit
gate asserts it bitwise). `potential: prescribed` reads a series/constant
rate. `potential: penman_monteith` is the T7 energy-budget form (next
section).

`atmosphere.aerodynamic_resistance.law` selects `r_a` for the ET
consumers only (the Q4/Q5 bulk consumers stay pinned to `liu2006`):
`liu2006` = 94.909 U^-0.9036; `fao56` = 208/u2 with
u2 = U·4.87/ln(67.8 z_ref − 5.42) (FAO-56 eq. 47); `log_profile` =
ln((z−d)/z0m)·ln((z−d)/z0h)/(κ² U).

## The Penman–Monteith potential (Stage-3 T7)

`potential: penman_monteith` evaluates the combination equation in
**specific-humidity form**, sharing the bulk chain's functions
(`src/atm/PenmanMonteith.hpp`) so the two modes are commensurable:

$$
\lambda \rho_w E_x(r_s) \;=\;
\frac{s\,A \;+\; \lambda \rho_a\,\bigl(q_{sat}(T_{air}) - q_a\bigr)/r_a}
     {s + 1 + r_s/r_a},
\qquad s = \frac{\lambda}{c_{pa}}\,\frac{dq_{sat}}{dT}\Big|_{T_{air}}
$$

with λ = `latentHeatOfVaporization(T_air)`, ρ_a = `airDensity(T_air, P)`,
c_pa = 1005 J/kg/K, dq_sat/dT the **analytic** Tetens derivative
(`dqsatdT`; unit-pinned against a central difference), and the available
energy

$$
A = R_n - G, \quad
R_n = (1-\alpha)\,R_s + \varepsilon\,(LW_{in} - \sigma T_K^4), \quad
G = f_G\,R_n
$$

at ε = 0.97 (`kWaterEmissivity`) and T_K = T_air + 273.15. R_s is the
**incident** shortwave `atmosphere.shortwave_incoming` (distinct from
the Q5 `shortwave`, which stays the absorbed SW of the heat exchange);
`longwave_in` is shared. The albedo α and ground-heat fraction f_G are
**per channel**: open water takes `energy.albedo_water` with the `water`
fraction (r_s = 0); the soil channel takes the top soil type's `albedo`
with the `soil` fraction (r_s = 0, then the `beta_lee_pielke`/`none`
limiter — `alpha_barton` acts inside the bulk humidity difference and is
schema-rejected with PM); the canopy takes the vegetation type's
`albedo` with the `canopy` fraction and r_s = r_c. Negative rates
(condensation-dominated) are returned as-is, matching the bulk chain's
V2-A13 decision (the canopy clamp T_pot ≥ 0 still applies).

**Derivation** (the standard elimination, in q): the surface energy
balance A = λρ_a(q_sat(T_s) − q_a)/r_a' + ρ_a c_pa(T_s − T_air)/r_a with
r_a' = r_a + r_s eliminates the unknown T_s by linearizing
q_sat(T_s) ≈ q_sat(T_air) + (dq_sat/dT)(T_s − T_air). **Commensurability**:
at the surface temperature T_s* that closes the balance, the PM rate and
the bulk rate E(T_s*) agree with an absolute gap that is second order in
(T_s* − T_air) — both rates vanish linearly, the secant-vs-tangent
linearization error enters at first order relative (measured slopes 2.03
and 0.98; `AtmEt.PmCommensurateWithBulkAtEnergyBalanceTemperature`,
S3-A10). The shipped bulk evaluates ρ_a at T_s, a recorded (ungated)
convention difference.

**FAO-56 convention difference**: FAO's reference chain (Chapter 4) is
the vapor-pressure form with the psychrometric constant
γ = c_pa P/(0.622 λ) and tabulated radiation conventions. The q-form
with the model's constants reproduces FAO's hourly worked example
(Example 19, N'Diaye 14.00–15.00 h, ET₀ = 0.63 mm/h) to a measured
convention gap δ_conv = 0.016 mm/h (2.6 %): the et0 part-4 battery pins
every intermediate of the transcribed examples
(`benchmarks/et-evapotranspiration/reference/fao56-examples.yaml`) to
the published precision, and gate et10 holds the model inside
δ_conv × 1.2 + print rounding of FAO's value on a grass column
configured to FAO's reference surface (LAI 2.88, r_c = 69.4 s/m,
albedo 0.23, G = 0.1 R_n).

PM requires the atmosphere block with `shortwave_incoming` and
`longwave_in`; it is rejected with `modules.temperature` (the heat
module's latent flux uses the bulk rate at the simulated temperature —
the combination is ungated).

## Soil-channel limiters, the supply cap, θ_fc

`alpha_barton` (bulk only) puts α₁(w_g) inside the humidity difference
with w_g the v2.0 surface extrapolation of the top two θ values (the
lines are copied verbatim from the Q4 path; et1 holds the chains to the
fma-contraction class). `beta_lee_pielke` multiplies the potential by
β = 1 for θ ≥ θ_fc, else 0.25(1 − cos(πθ/θ_fc))² (Lee & Pielke 1992);
θ_fc is `soil.types[].theta_fc`, defaulting to θ(h = −3.3 m) from the
van Genuchten curve. `none` applies the raw potential.

The **supply cap** bounds the applied flux by the top cell's physical
store per substep: evaporation takes at most `(θ − θ_r)·dz/dtg`;
condensation fills at most the pore room below the legacy 0.9999 θ_s
mark. `unrealised = demand − applied` is booked, never silently lost
(the et7 gate).

## Time-basis complementarity (no double counting)

Per cell, the surface channel covers the wet fraction of each step
(`wetFrac = min(1, depth/(E_w·dt))` when evaporating; its integral is the
`wet_time` clock) and the soil channel covers the complement: each
groundwater substep [a, b] reconstructs the wet clock at b by linear
growth, `W_b = wetClockStepStart + wetFrac·(b − stepStart)`, takes
`f_wet = clamp((W_b − wetClockLast)/dtg, 0, 1)`, marks the clock consumed,
and applies `f_dry = 1 − f_wet` to the soil demand. By construction
`wet_time + soil_evaporation_time ≤ t` per cell always (the et3/et9
gates); the v2.0 combination (surface bulk + `groundwater_top` flux)
double counts across the subcycled window lag — et3's real negative
measures the violation.

A coupled wet column (top-exchange mode ≠ 0) leaves evaporation entirely
to the surface channel: `qVap = 0` there.

## The D1/D2 fixes on the ET path

The applied flux `qVap` is read **guard-free** on both sides of the
substep — the predictor's RHS and the corrector's flux both take it as
is (the cap already bounded it physically). The v2.0 path guarded the
predictor with `θⁿ < θ_s` (skipping the sink on a saturated column) but
the corrector only with `θ > θ_r` — the D1 mismatch the et1b lane
measures. The coupled bookkeeping strips exactly `qVap·dtg` from the
exchange accumulator (any sign, no depth test) and books exactly that
volume as the channel's evaporation — not the whole top-face flux (D2).

## Salt: evaporation removes water only

- **Surface (O4, `transport.surface_dryout: store`)**: a drying cell's
  remaining scalar mass moves into the per-cell **salt store**
  [psu·m³], released into the water at the first wet step (after the
  limiter/bounds/dilution so the released mass is never clipped). The
  v2.0 `delete` books it into `surf_adjust` and loses it (the et8c
  negative). The transport audit's `Δ(surf_mass + surf_store)` closes on
  the unchanged RHS; the instantaneous `surf_store` column and the
  `salt_store` output expose the store.
- **Subsurface**: the vapour passes the top face with no scalar (the
  scalar-cauchy relation applied to ET columns), and the top-cell
  limiter admits the exact concentration factor from the window vapour
  volume, `f = (V_gflux + vap)/V_gflux` (condensation mirrors on the
  lower bound). This replaces the post-step-depth `coupledEvapTop` flag
  on ET columns (D5); non-ET columns keep the legacy branch.
- **Held seepage (D9)**: exfiltrated water's salt leaves the subsurface
  as the water does — each step's NEW exfiltration debits `exfNew·A·s_top`
  whether the exchange applied or is held below `min_depth`; a held
  cell's share lands in the salt store (the cell is dry), and the later
  application carries no additional exchange salt (its salt is already in
  the store, which the rewet step releases).
- Rain/evaporation dilution uses the actual per-cell volumes
  `(rainAdded − evapRemoved)·A` — rain and evaporation together, masks
  honoured (the D8 fix; legacy chose one or the other per cell).

## Vegetation and transpiration (Stage-3 T6)

Static per-cell vegetation (`vegetation.types` + an id map, 0 = bare;
optional per-cell LAI override) partitions the potential through
Beer–Lambert transmission τ = exp(−k·LAI):

- **ground rates** (open water under canopy, exposed soil) are scaled by
  τ — exactly 1.0 on bare cells, so vegetation-absent runs are bitwise;
- **potential transpiration** T_pot = (1 − τ)·E_x(r_s = r_c) with the
  big-leaf canopy resistance r_c = r_l/(0.5·LAI); bulk evaluates the
  canopy at T_s with r_c inside the chain, prescribed uses E_pot and
  ignores r_s; a condensing canopy does not transpire (T_pot ≥ 0).

Transpiration is a **non-compensated volumetric sink** over the root
zone: per cell S_k = min(α_F(hⁿ_k)·b_k·T_pot/Δz_k, (θⁿ_k − θ_r)/dtg)
[1/s], where b_k is the closed-form root-density integral over the
cell's depth interval (uniform 1/L or linear 2(L−d)/L², renormalised on
columns shallower than L with a warning) measured from the terrain-
following/masked column's top face, and α_F is the Feddes (1978)
piecewise-linear stress (h1 ≥ h2 ≥ h3 > h4; **h1/h2 omitted for wetland
species** — α = 1 at any h ≥ h3, so transpiration continues while
inundated, owner decision O7). The sink enters the predictor RHS and the
corrector numerator with one shared sign (a saturated confined cell
answers dh = −S·dtg/Ss with θ pinned at θ_s — the W6.4 unit gate), and
runs regardless of surface water: no depth test anywhere in the chain.

**Salt**: the uptake removes water only. The sink's window volume joins
the transport flux-volume subtraction (`vgflux −= sinkWin`) and the
limiter admits the exact per-cell factor f = (V + sinkWin)/V (combined
with the vapour factor on ET top cells), so c·θ·V is conserved per cell
and root-zone salinity rises in proportion to the uptake.

## Audits and outputs

`mass_audit.evaporation` keeps the open-water removal. `gw_mass_audit`
appends `soil_evaporation` and `soil_evaporation_unrealised` [m³]
(`soil_evaporation` is a sub-account of `boundary_in`) and, with
vegetation, `transpiration` and `transpiration_unrealised` — the sink
leaves through no boundary face, so the closure identity becomes
ΔV = boundary_in − ss_storage + realloc − vloss − transpiration (the
et4 gates assert both forms). The `/atmosphere` output group carries the interval-mean
rates (`evaporation_rate`, `soil_evaporation_rate`,
`soil_evaporation_unrealised_rate`; the T4 IntervalRates mechanism) and
the cumulative clocks (`wet_time`, `soil_evaporation_time`); T6 adds
`transpiration_rate`, `transpiration_unrealised_rate` and the 3-D
`groundwater: root_uptake` (interval-mean sink, 1/s — the registry's
Field3 leg).
Checkpoint: `et_wet_clock`, `et_wet_clock_last`, `et_soil_time`,
`s_store` (all prognostic; the et restart lanes hold them bitwise).
Run-record closure: `open_water_evaporation_m3`, `soil_evaporation_m3`,
`soil_evaporation_unrealised_m3`, `salt_store_final`.

## Deviations from v2.0 removed on the ET path

D1 predictor/corrector guard mismatch; D2 audit summing the whole top
flux; D4 no supply cap; D5 the post-step-depth transport flag; D6
dry-out salt deletion; D8 the rain/evap either/or dilution; D9
held-seepage salt; the subcycled window-lag double count. D3
(condensation/strip threshold) is subsumed by the cap + the guard-free
strip; D7 (solubility) stays out of scope. The legacy keys reproduce
v2.0 byte for byte (`default_bitwise`).
