#!/usr/bin/env python3
"""gate_heat.py — pure closed forms and check functions for the v2 Q5
temperature gates g6 (subsurface heat: Bredehoeft & Papadopulos steady
Peclet sweep / Ogata-Banks heat breakthrough / Stallman diel damping),
g7 (surface heat exchange: Edinger equilibrium relaxation / channel
thermal plume), and g8 (Horton-Rogers-Lapwood convection onset), plan
§4.2.

Every reference is computed in-script (plan §6.2: fully self-contained).
The §6.3 rule-1 cross-checks live in scripts/test_g678_gates.py — each
closed form is verified two independent ways (hand values, limiting
cases, residual substitution into the governing ODE/PDE). The HDF5
plumbing lives in run_regression.py.

Conventions: temperatures in degC (differences in K), thermal properties
volumetric ((rho c) in J/m^3/K), alpha denotes a thermal diffusivity
[m^2/s]; v_t denotes the THERMAL FRONT velocity (the water flux scaled
by (rho c)_w / (rho c)_bulk).
"""

from __future__ import annotations

import math

import numpy as np

WATER_HEAT_CAPACITY = 4.184e6  # (rho c)_w [J/m^3/K] (the g6(c) value)


# ---------------------------------------------------------------------------
# Measured-flux windows (V2-A17): the g6 gates drive the column through
# head values but never trust that arithmetic — they read the achieved
# Darcy flux from the output, assert it landed within a window of the
# target (so the sweep really spans the intended range), and evaluate
# the closed forms at the MEASURED value.
# ---------------------------------------------------------------------------

def check_flux_window(measured: float, target: float, window_rel: float,
                      floor_abs: float, label: str):
    """|measured - target| <= max(window_rel |target|, floor_abs). The
    absolute floor keeps the zero-flux sub-cases assertable."""
    limit = max(window_rel * abs(target), floor_abs)
    err = abs(measured - target)
    ok = err <= limit
    return ok, [f"{label}: measured flux {measured:.4e} vs target {target:.4e} "
                f"m/s (|d| = {err:.2e}, allowed {limit:.2e}) "
                f"{'ok' if ok else 'FAIL'}"]


# ---------------------------------------------------------------------------
# g6(a) — Bredehoeft & Papadopulos (1965) steady profile.
# ---------------------------------------------------------------------------

def bp_profile(xi, peclet: float):
    """Steady temperature between fixed ends, as the fraction
    f = (T - T(0)) / (T(L) - T(0)) at relative position xi = x/L:
    f = (exp(Pe xi) - 1) / (exp(Pe) - 1), linear in the Pe -> 0 limit.
    Pe = v_t L / alpha, signed with the flow direction (+ from the T(0)
    end toward the T(L) end)."""
    xi = np.asarray(xi, dtype=float)
    if abs(peclet) < 1.0e-12:
        return xi.copy()
    # expm1 keeps small-Pe precision.
    return np.expm1(peclet * xi) / math.expm1(peclet)


# ---------------------------------------------------------------------------
# g6(b) — Ogata & Banks (1961), heat form (the OGS benchmark).
# ---------------------------------------------------------------------------

def _exp_erfc(a: float, b: float) -> float:
    """exp(a) * erfc(b), overflow-safe (the second Ogata-Banks term has
    a ~ v x / alpha up to ~70 while erfc(b) underflows): for large b use
    the asymptotic erfc(b) = exp(-b^2)/(b sqrt(pi)) (1 - 1/(2 b^2) +
    3/(4 b^4))."""
    if b < 25.0:
        return math.exp(a) * math.erfc(b)
    if b > 0.0:
        series = 1.0 - 0.5 / (b * b) + 0.75 / (b ** 4)
        return math.exp(a - b * b) * series / (b * math.sqrt(math.pi))
    return math.exp(a) * math.erfc(b)


def ogata_banks(x: float, t: float, front_velocity: float, alpha: float,
                t_initial: float, t_step: float) -> float:
    """Temperature at distance x, time t for a step change t_step at the
    inlet of a semi-infinite column initially at t_initial:
    T = T0 + dT/2 [erfc((x - v t)/(2 sqrt(a t)))
                   + exp(v x / a) erfc((x + v t)/(2 sqrt(a t)))],
    v = thermal front velocity, a = thermal diffusivity."""
    if t <= 0.0:
        return t_initial
    denom = 2.0 * math.sqrt(alpha * t)
    term1 = math.erfc((x - front_velocity * t) / denom)
    term2 = _exp_erfc(front_velocity * x / alpha,
                      (x + front_velocity * t) / denom)
    return t_initial + 0.5 * (t_step - t_initial) * (term1 + term2)


# ---------------------------------------------------------------------------
# g6(c) — Stallman (1965) damped diel sinusoid under vertical flow.
# ---------------------------------------------------------------------------

def stallman_coefficients(front_velocity: float, alpha: float,
                          omega: float) -> tuple[float, float]:
    """(a, b) of T = A0 exp(-a z) sin(omega t - b z) satisfying
    dT/dt + v_t dT/dz = alpha d2T/dz2 (z positive downward, v_t positive
    downward). Substituting exp(-a z + i(omega t - b z)) gives
        alpha (a^2 - b^2) + v_t a = 0
        2 alpha a b + v_t b - omega = 0
    solved by damped Newton from the no-flow seed a = b =
    sqrt(omega / 2 alpha). The §6.3 checks: the v_t = 0 limit recovers
    the classical damping depth, and the residuals vanish at the root."""
    a = b = math.sqrt(omega / (2.0 * alpha))
    for _ in range(200):
        f1 = alpha * (a * a - b * b) + front_velocity * a
        f2 = 2.0 * alpha * a * b + front_velocity * b - omega
        j11 = 2.0 * alpha * a + front_velocity
        j12 = -2.0 * alpha * b
        j21 = 2.0 * alpha * b
        j22 = 2.0 * alpha * a + front_velocity
        det = j11 * j22 - j12 * j21
        da = (f1 * j22 - f2 * j12) / det
        db = (j11 * f2 - j21 * f1) / det
        a -= da
        b -= db
        if abs(da) + abs(db) < 1.0e-14 * (abs(a) + abs(b)):
            break
    return a, b


def stallman_amplitude_phase(z: float, front_velocity: float, alpha: float,
                             period: float) -> tuple[float, float]:
    """(amplitude ratio A(z)/A0, phase lag [s]) at depth z."""
    omega = 2.0 * math.pi / period
    a, b = stallman_coefficients(front_velocity, alpha, omega)
    return math.exp(-a * z), b * z / omega


def fit_sinusoid(times, values, period: float) -> tuple[float, float, float]:
    """(amplitude, phase [s], mean) least-squares fit of
    m + A sin(omega (t - lag)) to a sampled signal — the analysis the
    Stallman gate applies identically to the model and to itself."""
    t = np.asarray(times, dtype=float)
    v = np.asarray(values, dtype=float)
    omega = 2.0 * math.pi / period
    design = np.column_stack([np.sin(omega * t), np.cos(omega * t),
                              np.ones_like(t)])
    coef, *_ = np.linalg.lstsq(design, v, rcond=None)
    s, c, mean = coef
    amplitude = float(math.hypot(s, c))
    # m + A sin(w t + phi), phi = atan2(c, s); lag = -phi/omega.
    lag = float(-math.atan2(c, s) / omega)
    return amplitude, lag, float(mean)


# ---------------------------------------------------------------------------
# g7(a) — Edinger et al. (1968) equilibrium-temperature relaxation.
# ---------------------------------------------------------------------------

def edinger_relaxation(t, t0: float, t_equilibrium: float, k_e: float,
                       depth: float,
                       rho_c: float = WATER_HEAT_CAPACITY):
    """T(t) = T_e + (T0 - T_e) exp(-K_e t / (rho c_p h))."""
    tau = rho_c * depth / k_e
    return t_equilibrium + (t0 - t_equilibrium) * np.exp(-np.asarray(t, dtype=float) / tau)


# ---------------------------------------------------------------------------
# g7(a) stage 2 — the offline bulk-formula equilibrium (V2-A17). This is
# the SAME formula chain as src/atm/BulkAerodynamic.hpp plus the Q5 heat
# terms, computed independently in python. §6.3 rule-1 spot values in
# scripts/test_g678_gates.py: Q_net(20 C) and Q_net(22 C) for the g7a2
# forcing checked against hand-evaluated terms; the root bracketed by
# their sign change.
# ---------------------------------------------------------------------------

STEFAN_BOLTZMANN = 5.670374419e-8   # [W/m^2/K^4]
EMISSIVITY = 0.97                   # water surface (Q_lw = eps (LW_in - sigma T^4))
AIR_HEAT_CAPACITY = 1005.0          # c_pa [J/kg/K]
GAS_CONSTANT_DRY_AIR = 287.05       # [J/kg/K]
WATER_DENSITY = 1000.0              # [kg/m^3]


def tetens(temp_c: float) -> float:
    """Saturated vapor pressure [kPa] (Tetens, the atm-module form)."""
    return 0.6108 * math.exp(17.27 * temp_c / (temp_c + 237.3))


def q_sat(temp_c: float, pressure_kpa: float) -> float:
    esat = tetens(temp_c)
    return 0.622 * esat / (pressure_kpa - 0.376 * esat)


def r_air(wind_ms: float) -> float:
    """Liu et al. (2006) aerodynamic resistance [s/m]."""
    return 94.909 * wind_ms ** -0.9036


def air_density(temp_c: float, pressure_kpa: float) -> float:
    return pressure_kpa * 1000.0 / (GAS_CONSTANT_DRY_AIR * (temp_c + 273.15))


def latent_heat(temp_c: float) -> float:
    """L_v(T) [J/kg]."""
    return 2.501e6 - 2370.0 * temp_c


def net_heat_flux(t_water: float, t_air: float, pressure_kpa: float,
                  wind_ms: float, q_air: float, shortwave: float,
                  longwave_in: float) -> float:
    """Q_net [W/m^2, positive into the water] =
    Q_sw + eps (LW_in - sigma T_K^4) - Q_lat - Q_sens. The latent term
    reuses the Q4 atm-module chain verbatim — including its air density
    at the WATER temperature (atm::evaporationRate evaluates rho_a at
    its surface-temperature argument); the sensible term uses rho_a at
    the air temperature (its own, air-side quantity). Both evaluate at
    the LOCAL water temperature (V2-A17), never at
    atmosphere.surface_temperature."""
    t_k = t_water + 273.15
    lw = EMISSIVITY * (longwave_in - STEFAN_BOLTZMANN * t_k ** 4)
    evap = (air_density(t_water, pressure_kpa) / r_air(wind_ms) *
            (q_sat(t_water, pressure_kpa) - q_air) / WATER_DENSITY)
    q_lat = WATER_DENSITY * latent_heat(t_water) * evap
    q_sens = (air_density(t_air, pressure_kpa) * AIR_HEAT_CAPACITY *
              (t_water - t_air) / r_air(wind_ms))
    return shortwave + lw - q_lat - q_sens


def equilibrium_root(t_air: float, pressure_kpa: float, wind_ms: float,
                     q_air: float, shortwave: float, longwave_in: float,
                     lo: float = -5.0, hi: float = 60.0) -> float:
    """Bisection root of Q_net(T) = 0 (Q_net is strictly decreasing in T)."""
    flo = net_heat_flux(lo, t_air, pressure_kpa, wind_ms, q_air, shortwave, longwave_in)
    fhi = net_heat_flux(hi, t_air, pressure_kpa, wind_ms, q_air, shortwave, longwave_in)
    if not (flo > 0.0 > fhi):
        raise ValueError(f"equilibrium root not bracketed: Q({lo})={flo}, Q({hi})={fhi}")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = net_heat_flux(mid, t_air, pressure_kpa, wind_ms, q_air, shortwave, longwave_in)
        if fm > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-10:
            break
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# g7(b) — channel thermal plume: steady and transient ADE + first-order
# relaxation toward ambient T_e (van Genuchten's A1 solution family).
# ---------------------------------------------------------------------------

def channel_steady(x, velocity: float, dispersion: float, decay: float):
    """Steady excess-temperature fraction under u dT/dx = D d2T/dx2 - k T:
    T/T_in = exp(x (u - sqrt(u^2 + 4 k D)) / (2 D))."""
    w = math.sqrt(velocity * velocity + 4.0 * decay * dispersion)
    return np.exp(np.asarray(x, dtype=float) * (velocity - w) / (2.0 * dispersion))


def channel_transient(x: float, t: float, velocity: float, dispersion: float,
                      decay: float) -> float:
    """Excess-temperature fraction for a step release at x = 0, t = 0
    (ADE + decay, van Genuchten A1):
    T/T_in = 1/2 [exp(x(u-w)/2D) erfc((x - w t)/(2 sqrt(D t)))
                + exp(x(u+w)/2D) erfc((x + w t)/(2 sqrt(D t)))],
    w = sqrt(u^2 + 4 k D)."""
    if t <= 0.0:
        return 0.0
    w = math.sqrt(velocity * velocity + 4.0 * decay * dispersion)
    denom = 2.0 * math.sqrt(dispersion * t)
    term1 = _exp_erfc(x * (velocity - w) / (2.0 * dispersion),
                      (x - w * t) / denom)
    term2 = _exp_erfc(x * (velocity + w) / (2.0 * dispersion),
                      (x + w * t) / denom)
    return 0.5 * (term1 + term2)


# ---------------------------------------------------------------------------
# g8 — Horton-Rogers-Lapwood onset.
# ---------------------------------------------------------------------------

RAYLEIGH_CRITICAL = 4.0 * math.pi * math.pi  # 39.478


def rayleigh_number(conductivity: float, beta_t: float, delta_t: float,
                    height: float, alpha_effective: float) -> float:
    """Darcy-Rayleigh number in hydraulic-conductivity form:
    Ra = K beta_T dT H / alpha_e, alpha_e = lambda_bulk / (rho c)_w —
    the buoyant Darcy flux scale K beta dT against the conductive scale
    alpha_e / H. The gate evaluates it at H_eff between the pinned-cell
    centres (V2-A17)."""
    return conductivity * beta_t * delta_t * height / alpha_effective


def hrl_mode_amplitude(temp_xz: np.ndarray, x: np.ndarray, depth: np.ndarray,
                       z_top_pin: float, h_eff: float, width: float) -> float:
    """Amplitude of the seeded critical mode cos(pi x / H) sin(pi zeta)
    in a temperature slice (nx, nz): project T - <T>_x (the horizontal
    mean removes the conduction profile without prescribing it) onto the
    normalized mode. zeta = (depth - z_top_pin)/H_eff."""
    zeta = np.clip((np.asarray(depth, dtype=float) - z_top_pin) / h_eff, 0.0, 1.0)
    mode = (np.cos(math.pi * np.asarray(x, dtype=float) / (width / 2.0))[:, None] *
            np.sin(math.pi * zeta)[None, :])
    anomaly = temp_xz - np.mean(temp_xz, axis=0, keepdims=True)
    return float(np.sum(anomaly * mode) / np.sum(mode * mode))


def hrl_nusselt(temp_xz: np.ndarray, qz_xz: np.ndarray, k_face: int,
                dz: float, conductivity_w_per_mk: float, rho_c_w: float,
                delta_t: float, h_eff: float) -> float:
    """Nusselt number at the horizontal face below layer k_face: the
    x-averaged advective + conductive vertical heat flux over the
    conductive reference lambda dT / H_eff. qz_xz is the per-area Darcy
    flux at that face (positive up, the output 'qz' at plane k_face+1
    stored per cell as the lower-face flux); the face temperature is the
    two-cell mean (an analysis choice — the 1.05 threshold carries the
    margin)."""
    t_face = 0.5 * (temp_xz[:, k_face] + temp_xz[:, k_face + 1])
    # Positive-up flux of a bottom-heated slab transports heat up where
    # q and T anomalies correlate; the conductive part uses the discrete
    # gradient across the face.
    advective = rho_c_w * qz_xz * t_face
    conductive = conductivity_w_per_mk * (temp_xz[:, k_face + 1] - temp_xz[:, k_face]) / dz
    reference = conductivity_w_per_mk * delta_t / h_eff
    return float(np.mean(advective + conductive) / reference)


# ---------------------------------------------------------------------------
# Heat x-battery (plan §8.1 + V2-A17): the 8 dihedral variants of the
# square base case. Each entry pairs the ARRAY op (mapping a base (ny,
# nx) layer onto the variant's layer) with the COORDINATE op (mapping a
# base (x, y) point onto the variant's frame) — the harness builds the
# variant IC rasters with the array op and the variant BC polygons with
# the coordinate op, and the gate compares array_op(base_field) against
# the variant's field, so no inverse is ever formed and the two ops'
# consistency is provable on a synthetic field
# (scripts/test_g678_gates.py).
# ---------------------------------------------------------------------------

def dihedral_table(length: float):
    """name -> (array_op, coord_op) for a square domain of side length."""
    ll = length
    return {
        "id": (lambda f: f, lambda x, y: (x, y)),
        "flipx": (lambda f: f[:, ::-1], lambda x, y: (ll - x, y)),
        "flipy": (lambda f: f[::-1, :], lambda x, y: (x, ll - y)),
        "rot180": (lambda f: f[::-1, ::-1], lambda x, y: (ll - x, ll - y)),
        "transpose": (lambda f: f.T, lambda x, y: (y, x)),
        "antitranspose": (lambda f: f.T[::-1, ::-1], lambda x, y: (ll - y, ll - x)),
        "rot90": (lambda f: np.rot90(f, k=-1), lambda x, y: (ll - y, x)),
        "rot270": (lambda f: np.rot90(f, k=1), lambda x, y: (y, ll - x)),
    }


def check_heat_orientations(fields: dict, length: float, tol_abs: float,
                            label: str):
    """fields: name -> (ny, nx, nz) array per dihedral variant. Each
    variant must equal the base field pushed forward per layer by its
    array op, within tol_abs."""
    msgs = []
    ok = True
    base = fields["id"]
    for name, (array_op, _) in dihedral_table(length).items():
        if name == "id":
            continue
        pushed = np.stack([array_op(base[:, :, k])
                           for k in range(base.shape[2])], axis=2)
        dev = float(np.max(np.abs(pushed - fields[name])))
        line = (f"{label} orientation {name}: max |delta| = {dev:.3e} "
                f"(allowed {tol_abs:.1e})")
        if dev > tol_abs:
            ok = False
            line += " FAIL"
        msgs.append(line)
    return ok, msgs


# ---------------------------------------------------------------------------
# Checks.
# ---------------------------------------------------------------------------

def check_bp(peclets, max_errors, delta_t: float, tol_fraction: float):
    """g6(a): per-Peclet max-norm profile error <= tol_fraction of dT."""
    msgs = []
    ok = True
    limit = tol_fraction * abs(delta_t)
    for pe, err in zip(peclets, max_errors):
        line = f"B&P Pe={pe:+g}: max |T - closed form| = {err:.3e} K (allowed {limit:.2e})"
        if err > limit:
            ok = False
            line += " FAIL"
        msgs.append(line)
    return ok, msgs


def check_breakthrough(points, rel_errors, tol_rel: float):
    """g6(b): per-observation-point max relative error (of the step
    magnitude) over the breakthrough history."""
    msgs = []
    ok = True
    for x, err in zip(points, rel_errors):
        line = f"Ogata-Banks x={x:g} m: max rel error {err:.2%} (allowed {tol_rel:.0%})"
        if err > tol_rel:
            ok = False
            line += " FAIL"
        msgs.append(line)
    return ok, msgs


def check_stallman(depths, amp_model, amp_ref, lag_model, lag_ref,
                   amp_tol_rel: float, lag_tol_s: float, label: str):
    """g6(c): amplitude ratios within amp_tol_rel (relative to the
    reference ratio) and phase lags within lag_tol_s, gated separately
    (amplitude catches dissipation, phase catches dispersion error)."""
    msgs = []
    ok = True
    for z, am, ar, lm, lr in zip(depths, amp_model, amp_ref, lag_model, lag_ref):
        rel = abs(am - ar) / ar
        line = (f"{label} z={z:g} m: amplitude ratio {am:.4f} vs {ar:.4f} "
                f"(rel {rel:.2%}, allowed {amp_tol_rel:.0%})")
        if rel > amp_tol_rel:
            ok = False
            line += " FAIL"
        msgs.append(line)
        dlag = abs(lm - lr)
        line = (f"{label} z={z:g} m: phase lag {lm:.0f} s vs {lr:.0f} s "
                f"(|d| = {dlag:.0f} s, allowed {lag_tol_s:.0f})")
        if dlag > lag_tol_s:
            ok = False
            line += " FAIL"
        msgs.append(line)
    return ok, msgs


def check_relaxation(times, t_model, t0: float, t_e: float, k_e: float,
                     depth: float, tol_rel: float):
    """g7(a): T(t) within tol_rel of the exact exponential, relative to
    the initial excess |T0 - Te|."""
    ref = edinger_relaxation(times, t0, t_e, k_e, depth)
    err = float(np.max(np.abs(np.asarray(t_model, dtype=float) - ref)))
    limit = tol_rel * abs(t0 - t_e)
    ok = err <= limit
    return ok, [f"Edinger relaxation: max |T - exact| = {err:.3e} K "
                f"(allowed {limit:.2e}) {'ok' if ok else 'FAIL'}"]


def check_energy_ledger(residuals, tol_per_step: float):
    """g7(a): the per-step heat-ledger closure residual (the temperature
    analogue of the scalar-mass audit identity), relative units."""
    worst = float(np.max(np.abs(np.asarray(residuals, dtype=float))))
    ok = worst <= tol_per_step
    return ok, [f"energy ledger: worst per-step residual {worst:.2e} "
                f"(allowed {tol_per_step:.1e}) {'ok' if ok else 'FAIL'}"]


def check_equilibrium_consistency(t_model_final: float, t_e_offline: float,
                                  tol_abs: float):
    """g7(a) stage 2: the full bulk-formula run must settle at the T_e
    found offline by root-finding the net-flux zero."""
    err = abs(t_model_final - t_e_offline)
    ok = err <= tol_abs
    return ok, [f"bulk-mode equilibrium: settled {t_model_final:.3f} C vs "
                f"offline root {t_e_offline:.3f} C (|d| = {err:.3f} K, "
                f"allowed {tol_abs:g}) {'ok' if ok else 'FAIL'}"]


def check_channel(x, steady_model, steady_ref, l2_tol: float,
                  transient_model, transient_ref, transient_tol: float,
                  step_k: float):
    """g7(b): steady profile L2 (relative to the inlet excess) and the
    mid-channel transient breakthrough max error (of the step)."""
    msgs = []
    ok = True
    sm = np.asarray(steady_model, dtype=float)
    sr = np.asarray(steady_ref, dtype=float)
    l2 = float(np.sqrt(np.mean((sm - sr) ** 2))) / abs(step_k)
    line = f"channel steady L2: {l2:.2%} of the inlet excess (allowed {l2_tol:.0%})"
    if l2 > l2_tol:
        ok = False
        line += " FAIL"
    msgs.append(line)
    tm = np.asarray(transient_model, dtype=float)
    tr = np.asarray(transient_ref, dtype=float)
    err = float(np.max(np.abs(tm - tr))) / abs(step_k)
    line = (f"channel transient (mid-channel): max {err:.2%} of the step "
            f"(allowed {transient_tol:.0%})")
    if err > transient_tol:
        ok = False
        line += " FAIL"
    msgs.append(line)
    return ok, msgs


def check_hrl(perturbation_history_sub, nusselt_history_super,
              decay_factor: float, nusselt_min: float):
    """g8: at Ra = 28.5 the seeded-mode amplitude (signed projection)
    must decay by at least decay_factor from its seeded value; at
    Ra = 52.25 the Nusselt number must exceed nusselt_min persistently
    (over the trailing half of the history)."""
    msgs = []
    ok = True
    sub = np.abs(np.asarray(perturbation_history_sub, dtype=float))
    ratio = float(sub[-1] / sub[0]) if sub[0] > 0 else float("nan")
    line = (f"HRL Ra=28.5: |mode amplitude| ratio end/seed = {ratio:.3e} "
            f"(must decay below {decay_factor:g})")
    if not ratio < decay_factor:
        ok = False
        line += " FAIL"
    msgs.append(line)
    nu = np.asarray(nusselt_history_super, dtype=float)
    tail = nu[nu.size // 2:]
    line = (f"HRL Ra=52.25: trailing Nusselt min {float(np.min(tail)):.3f} "
            f"(must exceed {nusselt_min:g} persistently)")
    if not np.all(tail > nusselt_min):
        ok = False
        line += " FAIL"
    msgs.append(line)
    return ok, msgs
