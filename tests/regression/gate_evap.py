#!/usr/bin/env python3
"""gate_evap.py — pure check functions for the v2 Q4 gates g4 (analytic
drawdown / evaporative concentration) and g5 (Geng & Boufadel 2015
bare-soil salinization), plan §3.3 as amended by V2-A13 and V2-A14.

Every function takes plain arrays and returns (ok, messages); the HDF5
plumbing lives in run_regression.py. Keeping the criteria pure lets the
§6.3 negative battery (test_g45_gates.py) prove each check can fail on
synthetic data without running the model.
"""

from __future__ import annotations

import math

import numpy as np

# ---------------------------------------------------------------------------
# The bulk-aerodynamic closed form (plan §3.2; the g4(c) offline reference
# and the g5(i) E(0) anchor of V2-A13).
# ---------------------------------------------------------------------------

R_DRY_AIR = 287.05     # [J/(kg K)]
RHO_WATER = 1000.0     # [kg/m^3]


def tetens_esat_kpa(temp_c: float) -> float:
    """Saturated vapor pressure [kPa] (Tetens; paper Eq. 5)."""
    return 0.6108 * math.exp(17.27 * temp_c / (temp_c + 237.3))


def saturated_specific_humidity(temp_c: float, pressure_kpa: float) -> float:
    """q_sat [-] from e_sat (paper Eq. 4)."""
    esat = tetens_esat_kpa(temp_c)
    return 0.622 * esat / (pressure_kpa - 0.376 * esat)


def aerodynamic_resistance(wind_m_s: float) -> float:
    """R_air [s/m] (Liu et al. 2006 fit; plan §3.2)."""
    return 94.909 * wind_m_s ** (-0.9036)


def bulk_evaporation_rate(surface_temp_c: float, pressure_kpa: float,
                          wind_m_s: float, q_air: float,
                          alpha1: float = 1.0) -> float:
    """Potential (alpha1 = 1) or humidity-limited evaporation rate [m/s]
    (Mahfouf & Noilhan form, plan §3.2). Positive = evaporation."""
    q_ground = alpha1 * saturated_specific_humidity(surface_temp_c, pressure_kpa)
    rho_air = pressure_kpa * 1000.0 / (R_DRY_AIR * (surface_temp_c + 273.15))
    return (rho_air / aerodynamic_resistance(wind_m_s)) * (q_ground - q_air) / RHO_WATER


def soil_relative_humidity(water_content: float) -> float:
    """alpha_1(w_g) (paper Eq. 6, Barton / Lee & Pielke)."""
    return min(1.0, 1.8 * water_content / (water_content + 0.30))


# ---------------------------------------------------------------------------
# g4 checks.
# ---------------------------------------------------------------------------

def check_drawdown(times, eta_max_abs_err, eta0, rate, total_drawdown,
                   tol_fraction):
    """g4(a)/(c) drawdown: at each output time, the max-over-cells absolute
    deviation of eta from eta0 - rate*t must stay within tol_fraction of
    the total drawdown. `eta_max_abs_err` is that per-time deviation."""
    msgs = []
    ok = True
    limit = tol_fraction * total_drawdown
    for t, err in zip(times, eta_max_abs_err):
        line = f"drawdown t={t:g}: |eta - analytic| = {err:.3e} (allowed {limit:.3e})"
        if err > limit:
            ok = False
            line += " FAIL"
        msgs.append(line)
    return ok, msgs


def check_concentration(times, s_values, depths, s0, depth0, tol_rel,
                        mass_values, tol_mass_rel):
    """g4(b)/(c) evaporative concentration: s(t) = s0*V0/V(t) to tol_rel,
    and total salt mass constant to tol_mass_rel of its initial value."""
    msgs = []
    ok = True
    for t, s, d in zip(times, s_values, depths):
        if d <= 0.0:
            continue
        s_ref = s0 * depth0 / d
        rel = abs(s - s_ref) / s_ref
        line = f"concentration t={t:g}: s = {s:.6f}, analytic {s_ref:.6f} (rel {rel:.2e})"
        if rel > tol_rel:
            ok = False
            line += " FAIL"
        msgs.append(line)
    mass0 = mass_values[0]
    drift = float(np.max(np.abs(np.asarray(mass_values) - mass0))) / abs(mass0)
    line = f"salt mass: max drift {drift:.2e} of initial (allowed {tol_mass_rel:.1e})"
    if drift > tol_mass_rel:
        ok = False
        line += " FAIL"
    msgs.append(line)
    return ok, msgs


def check_closure(delta_volume, rain, evap, outflow, bc_inflow, reference,
                  tol_fraction):
    """Mass-audit closure |dV + outflow - rain + evap - bc| <= tol*reference.
    g4(d)'s teeth: pre-Q4, the evaporation ledger keeps accumulating the
    potential rate after the basin dries, so the identity breaks by the
    unaudited shortfall."""
    error = delta_volume + outflow - rain + evap - bc_inflow
    limit = tol_fraction * reference
    ok = abs(error) <= limit
    msg = (f"closure: |{error:.4e}| m^3 (allowed {limit:.3e}) "
           f"{'ok' if ok else 'FAIL'}")
    return ok, [msg]


def check_positivity_and_dry(depth_min_per_time, depth_final_max, min_depth):
    """g4(d): depth never negative at any output; fully dry at the end."""
    msgs = []
    ok = True
    worst = float(np.min(depth_min_per_time))
    line = f"positivity: min depth over outputs {worst:.3e} m"
    if worst < -1.0e-12:
        ok = False
        line += " FAIL"
    msgs.append(line)
    line = f"dry-out: final max depth {depth_final_max:.3e} m (must be <= {min_depth:g})"
    if depth_final_max > min_depth:
        ok = False
        line += " FAIL"
    msgs.append(line)
    return ok, msgs


# ---------------------------------------------------------------------------
# g5 checks (V2-A13 criteria for (i); plan §3.3 for (ii)/(iii); V2-A14
# for (iv)).
# ---------------------------------------------------------------------------

def check_rate_series(t_hours, e_m_s, e0_closed_form, e10_reference,
                      surface_saturation_50h, decay_bound=0.5,
                      saturation_band=(0.0915, 0.5)):
    """g5(i), the V2-A15 form. Gated: E(0) within 5 % of the Table-1 closed
    form; the rate monotone decreasing after the stage-1 plateau (t > 3 h);
    E(50h)/E(0) <= decay_bound (a missing alpha_1 limiter holds the ratio
    near 1); the extrapolated surface saturation at 50 h inside
    saturation_band (floor = the closed-form alpha_1 equilibrium — crossing
    it implies a humidity-gradient sign error). Recorded, not gated: the
    Fig.-3 comparison at 10 h (V2-A15: the E(t) decay rate is an
    internodal-conductivity artifact in both codes)."""
    t = np.asarray(t_hours, dtype=float)
    e = np.asarray(e_m_s, dtype=float)
    msgs = []
    ok = True

    e0 = float(e[0])
    ratio0 = e0 / e0_closed_form
    line = f"E(0) = {e0:.4e} m/s vs closed form {e0_closed_form:.3e} (ratio {ratio0:.3f})"
    if not 0.95 <= ratio0 <= 1.05:
        ok = False
        line += " FAIL"
    msgs.append(line)

    tail = e[t > 3.0]
    drops = np.diff(tail) <= 1.0e-12
    frac = float(np.mean(drops)) if drops.size else 0.0
    line = f"monotone decay after 3 h: {frac:.1%} of intervals decreasing"
    if frac < 0.99:
        ok = False
        line += " FAIL"
    msgs.append(line)

    e50 = float(np.interp(50.0, t, e))
    decay = e50 / e0
    line = f"E(50h)/E(0) = {decay:.4f} (allowed <= {decay_bound:g}; measured headroom V2-A15)"
    if decay > decay_bound:
        ok = False
        line += " FAIL"
    msgs.append(line)

    lo, hi = saturation_band
    line = (f"extrapolated surface saturation at 50h = {surface_saturation_50h:.4f} "
            f"(band [{lo:g}, {hi:g}]; floor = closed-form alpha1 equilibrium)")
    if not lo <= surface_saturation_50h <= hi:
        ok = False
        line += " FAIL"
    msgs.append(line)

    e10 = float(np.interp(10.0, t, e))
    msgs.append(f"recorded (not gated): E(10h) = {e10:.4e} m/s vs digitized Fig. 3 "
                f"{e10_reference:.3e} (factor {e10 / e10_reference:.2f}; V2-A15)")
    return ok, msgs


def _rms_on_overlap(z_grid, model, reference):
    """RMS of model - reference where both are defined; count returned."""
    m = np.asarray(model, dtype=float)
    r = np.asarray(reference, dtype=float)
    mask = ~(np.isnan(m) | np.isnan(r))
    n = int(np.sum(mask))
    if n == 0:
        return float("nan"), 0
    return float(np.sqrt(np.mean((m[mask] - r[mask]) ** 2))), n


def record_moisture_profiles(z_grid, model20, model50, ref20, ref50):
    """g5(ii) moisture RMS vs digitized Fig. 4 — RECORDED, not gated
    (V2-A15: the profile figures are mass-inconsistent with the Table-1
    forcing by the identified ~4.3x factor)."""
    msgs = []
    for label, model, ref in (("20h", model20, ref20), ("50h", model50, ref50)):
        rms, n = _rms_on_overlap(z_grid, model, ref)
        msgs.append(f"recorded (not gated): moisture {label} RMS {rms:.4f} over "
                    f"{n} levels vs digitized Fig. 4 (V2-A15)")
    return True, msgs


def check_salinity_profiles(z_grid, model20, model50, ref20, ref50,
                            peak_min_g_l, peak_above_z):
    """g5(ii), the V2-A15 form. Gated: the near-surface salinization
    signature — 50 h peak > peak_min_g_l above z = peak_above_z m (a
    missing scalar_cauchy condition caps it near the initial 25 g/L).
    Recorded, not gated: the profile RMS vs digitized Fig. 9b."""
    msgs = []
    ok = True
    for label, model, ref in (("20h", model20, ref20), ("50h", model50, ref50)):
        rms, n = _rms_on_overlap(z_grid, model, ref)
        msgs.append(f"recorded (not gated): salinity {label} RMS {rms:.2f} g/L over "
                    f"{n} levels vs digitized Fig. 9b (V2-A15)")
    z = np.asarray(z_grid, dtype=float)
    m50 = np.asarray(model50, dtype=float)
    near = ~np.isnan(m50) & (z > peak_above_z)
    peak = float(np.max(m50[near])) if np.any(near) else float("nan")
    line = (f"50h near-surface peak: {peak:.1f} g/L above z = {peak_above_z:g} m "
            f"(required > {peak_min_g_l:g})")
    if not peak > peak_min_g_l:
        ok = False
        line += " FAIL"
    msgs.append(line)
    return ok, msgs


def check_salt_mass(mass_series, tol_fraction):
    """g5(iii): total salt mass conserved to tol_fraction (the paper's own
    MARUN budget bound, 4 %; do not gate tighter than the reference)."""
    mass = np.asarray(mass_series, dtype=float)
    drift = float(np.max(np.abs(mass - mass[0]))) / abs(float(mass[0]))
    ok = drift <= tol_fraction
    msg = (f"salt mass: max drift {drift:.4f} of initial "
           f"(allowed {tol_fraction:g}) {'ok' if ok else 'FAIL'}")
    return ok, [msg]


def plume_edge_min_elevation(z_centers, mean_salinity_xz, edge_g_l):
    """Deepest cell-centre elevation with salinity >= edge_g_l anywhere in
    the slice (the Fig. 7 '30 g/L plume edge', fingers included).
    `mean_salinity_xz` is (nx, nz) — the y-averaged x-z slice."""
    s = np.asarray(mean_salinity_xz, dtype=float)
    z = np.asarray(z_centers, dtype=float)
    hit = np.any(s >= edge_g_l, axis=0)
    if not np.any(hit):
        return float("nan")
    return float(np.min(z[hit]))


def check_density_contrast(zmin_beta, zmin_control, min_margin):
    """g5(iv), the V2-A14 direction: at 50 h the beta = 7.44e-4 run's
    30 g/L plume edge must reach at least min_margin deeper than the
    beta = 0 control's, and both runs must have a plume at all."""
    msgs = []
    ok = True
    if math.isnan(zmin_beta) or math.isnan(zmin_control):
        return False, ["plume edge: missing in one of the runs FAIL"]
    margin = zmin_control - zmin_beta
    line = (f"density contrast: edge at z = {zmin_beta:.3f} m (beta run) vs "
            f"{zmin_control:.3f} m (beta = 0); margin {margin:.3f} m "
            f"(required >= {min_margin:g})")
    if margin < min_margin:
        ok = False
        line += " FAIL"
    msgs.append(line)
    return ok, msgs
