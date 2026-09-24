# Temperature (v2 Q5)

Temperature is the model's second registered scalar (v2 plan §4.1,
realization V2-A17): it rides the same finite-volume transport
machinery as salinity — the advection schemes, the local min/max
limiter, the decomposition-invariant far-neighbor staging, the ledger
architecture — instantiated a second time with the thermal physics
switched in. Everything below is the delta against
[Transport](transport.md); anything not mentioned behaves identically.

## Governing forms

Subsurface (per column cell, water flux **q** from the completed
Richards substep):

$$\frac{\partial}{\partial t}\big[(\theta + \kappa)\,T\big]
  + \nabla\cdot(\mathbf{q}\,T)
  = \nabla\cdot\big(\mathbf{D}_T \nabla T\big),
  \qquad \kappa = (1-\theta_s)\,\frac{(\rho c)_s}{(\rho c)_w}$$

dividing the volumetric heat equation by \((\rho c)_w\): the **thermal
retardation** \(\kappa\) (the SEAWAT temperature-as-species form) joins
the moisture basis in the scalar mass, the flux volume, and the
owned-mass reduction, so the thermal front travels at
\(v_T = q\,(\rho c)_w/(\rho c)_b\) with
\((\rho c)_b = \theta(\rho c)_w + (1-\theta_s)(\rho c)_s\).

The dispersion tensor's molecular slot carries the **effective thermal
diffusivity** \(\alpha_e = \lambda_{\mathrm{eff}}/(\rho c)_w\) — a bulk
property used as-is, *not* multiplied by \(\theta_s\) the way solute
molecular diffusion is (conduction passes through grains and water
alike; `temperature.thermal_conductivity` is the saturated-bulk
\(\lambda_{\mathrm{eff}}\), held constant in v2.0 even in unsaturated
cells). The mechanical terms remain available as thermal dispersivity
(`temperature.dispersivity`, default 0).

Surface (depth-integrated, per wet cell):

$$\frac{\partial (h\,T)}{\partial t} + \nabla\cdot(\mathbf{u}\,h\,T)
  = \nabla\cdot(K\,h\,\nabla T) + \frac{Q_{\mathrm{net}}}{(\rho c)_w}$$

with the same superbee/upwind advection on the pre-correction flow
rates and `temperature.surface_diffusivity` in the diffusion slot.

## Where temperature deviates from the salinity code path

All of these are keyed on the scalar spec (V2-A17); the salinity
instance is byte-identical to v1.

1. **Uncoupled top faces advect the donor value in both directions.**
   Heat travels with the water: an evaporative top flux removes heat at
   the cell's own temperature (T is invariant under pure mass removal —
   exact when the transport step equals the subsurface substep), and
   recharge under a flux/head top enters at the zero-gradient ghost
   value (the cell itself). The salinity-specific zero-advection and
   evaporative-concentration machinery (`scalar_cauchy`, the legacy
   allowance, the exact in-step factor) never applies. To prescribe an
   *inflow* temperature at a top or bottom face, pin the boundary cell
   (below).
2. **Rain/evaporation dilution is skipped.** Rain enters and
   evaporation leaves at the cell's own temperature — the volume change
   carries no temperature change. The rain-heat approximation this
   makes is *measured* into the surface anchor ledger column, not
   hidden. (Latent-heat cooling is a flux term, not a dilution term.)
3. **Cell-pinning Dirichlet rows.** `scalar_value` with
   `scalar: temperature` on `groundwater_top`/`groundwater_bottom`
   re-imposes the member columns' top/bottom cells after every update —
   the b6 tide rule applied vertically (temperature-only in v2.0; the
   pin delta is booked into `subs_boundary`). The effective isothermal
   plane sits at the pinned-cell *centre*: analytic gates measure
   domain lengths between pinned centres (g6a's \(L-\mathrm{d}z\),
   g8's \(H_{\mathrm{eff}} = H-\mathrm{d}z\)), and the pinned cell's
   half-thickness is a boundary-position ambiguity that decays with
   the front width (the g6b gating window).
4. **Side-ghost limiter admission on all four sides.** Legacy admits a
   prescribed side-ghost value into the limiter window on y+ only (the
   b6 sea-side rule, golden-pinned for salinity). For temperature the
   admission applies on every side with a non-no-flux code — the
   y+-only rule is orientation-asymmetric, and the heat 8-orientation
   battery caught exactly that (a side thermal Dirichlet was
   limiter-admitted on y+ but clipped on the other three sides).
   Related sharp edge, unchanged in v2.0: a side polygon that spans a
   domain corner hands the condition *every* domain-edge face of the
   corner cell (`BoundarySet` semantics, which b6's one-cell-wide tank
   relies on) — with the all-sides admission this is
   orientation-equivariant for temperature.
5. **Open default bounds.** `temperature.bounds` min/max both default
   open (salinity's lower bound defaults to 0); dry surface cells hold
   temperature 0 by the same convention that dry cells hold no scalar.

## Surface heat exchange

`temperature.surface_exchange` adds a per-wet-cell source inside the
update pass, *after* the limiter (an atmospheric flux must be able to
cross the local advective extrema — relaxation toward an equilibrium
below the wet-stencil minimum would otherwise stall), audited in the
`surf_atmos` ledger column:

- `mode: equilibrium` — the Edinger (1968) linearization
  \(Q = -K_e\,(T - T_e)\) (gate g7a).
- `mode: bulk` — \(Q_{\mathrm{net}} = Q_{sw} + \varepsilon\,(LW_{in} -
  \sigma T_K^4) - Q_{lat} - Q_{sens}\) with \(\varepsilon = 0.97\),
  prescribed absorbed shortwave and incident longwave
  (`atmosphere.shortwave`/`longwave_in`), and the latent/sensible terms
  on the shared Q4 bulk-aerodynamic chain (Tetens \(q_{sat}\), Liu
  \(R_{air}\) with the `wind_speed_floor` still-air bound,
  \(L_v(T) = 2.501{\times}10^6 - 2370\,T\)) evaluated at the **local
  water temperature** — never at `atmosphere.surface_temperature`,
  which only the Q4 evaporation-volume consumers read. The g7a2 gate
  root-finds the zero of this exact chain offline and requires the
  basin to settle there within 0.05 °C.

The coupled surface–subsurface heat exchange rides the seepage
(advective at the upwind temperature) plus the two-point conductive
term with the thermal tensor's top coefficient — the composed coupled
Stallman gate (g6c-coupled) verifies that path against the same
dispersion-relation roots as the prescribed-BC form.

## Density coupling

With `groundwater.density_coupling` enabled the cell ratios become

$$r_\rho = 1 + \beta_s\,s - \beta_T\,(T - T_0),
  \qquad r_\mu = \frac{1}{1 + \beta_{s\mu}\,s}$$

with all coefficients configuration (defaults = the legacy compile-time
constants; `thermal_expansion` \(\beta_T\) defaults 0). \(\mu(T)\) is
out of scope in v2.0 — \(r_\mu\) stays salinity-only, recorded. The
Horton–Rogers–Lapwood gate (g8) brackets the analytic convection
threshold \(Ra_c = 4\pi^2\) of exactly this coupling in the Darcy
momentum: \(Ra_{\mathrm{eff}} = K \beta_T \Delta T H_{\mathrm{eff}} /
\alpha_e = 28.5\) must stay conductive and \(52.25\) must convect.

## Ledger

An active temperature module writes `/monitor/temperature_audit`: the
transport-audit identity on the temperature scalar with the extra
`surf_atmos` column, in K·m³ (surface) and K·m³ on the
\(\theta+\kappa\) basis (subsurface):

    Δ(surf_heat) = exchange + surf_source + surf_boundary
                 + surf_adjust + surf_anchor + surf_atmos

closing to rounding per step — the g7a energy-ledger criterion gates it
at 10⁻⁸ of the basin heat content.
