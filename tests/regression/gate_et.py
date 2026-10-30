#!/usr/bin/env python3
"""gate_et.py — pure references and checks for the Stage-3 T5
evapotranspiration gates et1–et9 (Stage-3 plan §9).

The bulk-potential reference mirrors src/atm/EtFormulas.hpp term by term
(Tetens, the Geng/Boufadel specific-humidity and air-density chain, the
Liu 2006 / FAO-56 / log-profile aerodynamic resistances); the §6.3
battery (scripts/test_et_gates.py) proves it against an independent
second evaluation and every checker's ability to fail. The HDF5
plumbing lives in run_regression.py.
"""

from __future__ import annotations

import math

import numpy as np

GAS_CONSTANT_DRY_AIR = 287.05  # J/(kg K)   (EtFormulas kGasConstantDryAir)
WATER_DENSITY = 1000.0         # kg/m^3     (EtFormulas kWaterDensity)
VON_KARMAN = 0.41


def tetens(temp_c: float) -> float:
    return 0.6108 * math.exp(17.27 * temp_c / (temp_c + 237.3))


def saturated_specific_humidity(temp_c: float, pressure_kpa: float) -> float:
    esat = tetens(temp_c)
    return 0.622 * esat / (pressure_kpa - 0.376 * esat)


def air_density(temp_c: float, pressure_kpa: float) -> float:
    return pressure_kpa * 1000.0 / (GAS_CONSTANT_DRY_AIR * (temp_c + 273.15))


def aero_resistance(law: str, wind: float, ref_height: float = 2.0,
                    z0m: float = 0.01, z0h: float = 0.001,
                    displacement: float = 0.0) -> float:
    """r_a [s/m]: liu2006 (the v2.0 fit), fao56 (eq. 47, 208/u2 with the
    measurement-height wind conversion), log_profile."""
    if law == "liu2006":
        return 94.909 * wind ** -0.9036
    if law == "fao56":
        u2 = wind * 4.87 / math.log(67.8 * ref_height - 5.42)
        return 208.0 / u2
    if law == "log_profile":
        z = ref_height - displacement
        return (math.log(z / z0m) * math.log(z / z0h)
                / (VON_KARMAN * VON_KARMAN * wind))
    raise ValueError(law)


def soil_relative_humidity(water_content: float) -> float:
    """alpha_1(w_g) (Barton / Lee & Pielke; BulkAerodynamic.hpp)."""
    return min(1.0, 1.8 * water_content / (water_content + 0.30))


def potential_bulk(surface_temp_c: float, pressure_kpa: float, r_a: float,
                   r_s: float, q_air: float, alpha1: float = 1.0) -> float:
    """E [m/s] = rho_a (alpha1 q_sat(T_s) - q_a) / (rho_w (r_a + r_s)),
    ordered exactly as EtFormulas::potentialBulk (at r_s = 0 it is the
    v2.0 atm::evaporationRate expression tree)."""
    q_ground = alpha1 * saturated_specific_humidity(surface_temp_c, pressure_kpa)
    return (air_density(surface_temp_c, pressure_kpa) / (r_a + r_s)
            * (q_ground - q_air) / WATER_DENSITY)


def beta_lee_pielke(theta: float, theta_fc: float) -> float:
    if theta >= theta_fc:
        return 1.0
    return 0.25 * (1.0 - math.cos(math.pi * theta / theta_fc)) ** 2


def supply_cap(e: float, theta: float, theta_r: float, theta_s: float,
               dz: float, dtg: float) -> float:
    """Evaporation capped by the drainable store above theta_r;
    condensation by the pore room below the 0.9999 theta_s mark."""
    if e > 0.0:
        return min(e, max(0.0, theta - theta_r) * dz / dtg)
    return max(e, -max(0.0, 0.9999 * theta_s - theta) * dz / dtg)


def theta_fc_default(theta_r: float, theta_s: float, vg_alpha: float,
                     vg_n: float) -> float:
    """theta(h = -3.3 m) from van Genuchten (the SoilType default)."""
    m = 1.0 - 1.0 / vg_n
    se = (1.0 + (vg_alpha * 3.3) ** vg_n) ** (-m)
    return theta_r + se * (theta_s - theta_r)


def met_specific_humidity(air_temp_c: float, pressure_kpa: float,
                          relative_humidity: float) -> float:
    """q_a from RH (the MetForcing fold)."""
    return relative_humidity * saturated_specific_humidity(air_temp_c,
                                                           pressure_kpa)


# ---------------------------------------------------------------------------
# Checks.
# ---------------------------------------------------------------------------

def check_pond_phase(times, depths, cum_open, h0: float, e_w: float,
                     cum_exchange, area: float, tol_rel: float):
    """et3 (1): during the pond phase (t < t_d), cumulative open-water
    evaporation = E_w t and depth + cum_exchange/A = h0 - E_w t."""
    t_d = h0 / e_w
    ok = True
    parts = []
    worst_e = worst_d = 0.0
    n = 0
    for t, d, ce, cx in zip(times, depths, cum_open, cum_exchange):
        if t <= 0.0 or t >= t_d:
            continue
        n += 1
        worst_e = max(worst_e, abs(ce / area - e_w * t) / (e_w * t))
        # Surface ledger during the pond phase: depth = h0 - E_w t +
        # exchange/A (exchange = the audited cumulative seepage, positive
        # INTO the surface; the clamp is 0 until dry-out).
        worst_d = max(worst_d, abs(d - (h0 - e_w * t) - cx / area) / h0)
    if n == 0 or worst_e > tol_rel or worst_d > tol_rel:
        ok = False
    parts.append(f"pond phase ({n} outputs < t_d = {t_d:.1f} s): "
                 f"cum-open err {worst_e:.3e}, depth err {worst_d:.3e} "
                 f"(<= {tol_rel:.0e})")
    return ok, "; ".join(parts)


def check_no_double_count(wet_time: np.ndarray, soil_time: np.ndarray,
                          t: float):
    """et3 (2) / et9: per cell, wet_time + soil_evaporation_time never
    exceeds elapsed time."""
    total = np.nan_to_num(wet_time) + np.nan_to_num(soil_time)
    worst = float(np.max(total)) - t
    ok = worst <= 1.0e-9
    return ok, (f"time-basis complementarity: max(wet+soil) - t = "
                f"{worst:.3e} s (<= 1e-9)")


def check_combined_rate(cum_open: float, cum_soil: float, e_w: float,
                        t: float, area: float, tol_rel: float):
    """et3 (3): open + soil = E_w t while alpha_1 = 1."""
    err = abs((cum_open + cum_soil) / area - e_w * t) / (e_w * t)
    ok = err <= tol_rel
    return ok, (f"combined-rate identity at t = {t:g}: err {err:.3e} "
                f"(<= {tol_rel:.0e})")


def check_cell_salt(c: np.ndarray, theta: np.ndarray, c0: np.ndarray,
                    theta0: np.ndarray, tol_rel: float):
    """et8a: per cell, c theta = c0 theta0 (uniform dz cancels)."""
    m = np.isfinite(c) & np.isfinite(c0)
    now = c[m] * theta[m]
    ref = c0[m] * theta0[m]
    worst = float(np.max(np.abs(now - ref) / np.maximum(ref, 1e-300)))
    ok = worst <= tol_rel
    return ok, f"per-cell c*theta identity: worst rel err {worst:.3e} (<= {tol_rel:.0e})"


def check_pond_salt(s: float, h: float, s0h0: float, tol_rel: float,
                    label: str):
    """et8b: s h = s0 h0 whenever wet (evaporation removes water only,
    rain adds salt-free volume, the store carries the dry spell)."""
    err = abs(s * h - s0h0) / s0h0
    ok = err <= tol_rel
    return ok, f"{label}: s*h = {s * h:.6e} vs {s0h0:.6e} (rel err {err:.3e})"


def check_conserved(total: float, ref: float, tol_rel: float, label: str):
    err = abs(total - ref) / max(abs(ref), 1e-300)
    ok = err <= tol_rel
    return ok, f"{label}: {total:.10e} vs {ref:.10e} (rel err {err:.3e})"


def check_cap_trace(theta_min: float, theta_r: float, unrealised: float,
                    demand: float, applied: float, vloss: float,
                    tol_abs: float):
    """et7: the floor holds, the books balance, nothing is clamped."""
    ok = True
    parts = []
    if theta_min < theta_r - 1.0e-14:
        ok = False
    parts.append(f"min theta - theta_r = {theta_min - theta_r:.3e} (>= -1e-14)")
    err = abs(unrealised - (demand - applied))
    if err > tol_abs * max(demand, 1.0):
        ok = False
    parts.append(f"unrealised - (demand-applied) = {err:.3e}")
    if abs(vloss) > 1.0e-14:
        ok = False
    parts.append(f"vloss = {vloss:.3e} (== 0)")
    return ok, "; ".join(parts)


# ---------------------------------------------------------------------------
# Stage-3 T7 (plan §11): the FAO-56 reference chain (the e-form, FAO's own
# constants — et0 part 4 checks it against the transcribed Chapter 4
# examples) and the model's q-form PM (et10's 1e-10 reference).
# ---------------------------------------------------------------------------

STEFAN_DAY = 4.903e-9      # [MJ m-2 K-4 day-1] (FAO)
STEFAN_HOUR = 4.903e-9 / 24.0  # = 2.043e-10 [MJ m-2 K-4 h-1] (p. 74)


def fao_pressure(z: float) -> float:
    """FAO Eq. 7 [kPa]."""
    return 101.3 * ((293.0 - 0.0065 * z) / 293.0) ** 5.26


def fao_gamma(p_kpa: float) -> float:
    """FAO Eq. 8 [kPa/C]."""
    return 0.665e-3 * p_kpa


def fao_e0(t_c: float) -> float:
    """FAO Eq. 11 [kPa] (the same Tetens as the model's chain)."""
    return 0.6108 * math.exp(17.27 * t_c / (t_c + 237.3))


def fao_delta(t_c: float) -> float:
    """FAO Eq. 13 [kPa/C]."""
    return 4098.0 * fao_e0(t_c) / (t_c + 237.3) ** 2


def fao56_daily(t_max: float, t_min: float, rh_max: float, rh_min: float,
                u2: float, z: float, ra_mj: float, big_n: float,
                n_sun: float) -> dict:
    """The Example-18 chain (daily, G = 0, albedo 0.23, Angstrom
    0.25/0.50). Returns every intermediate the transcription prints."""
    out = {}
    out["pressure_kpa"] = fao_pressure(z)
    t_mean = 0.5 * (t_max + t_min)
    out["t_mean_c"] = t_mean
    out["delta_kpa_per_c"] = fao_delta(t_mean)
    out["gamma_kpa_per_c"] = fao_gamma(out["pressure_kpa"])
    out["one_plus_034u2"] = 1.0 + 0.34 * u2
    denom = out["delta_kpa_per_c"] + out["gamma_kpa_per_c"] * out["one_plus_034u2"]
    out["delta_term"] = out["delta_kpa_per_c"] / denom
    out["gamma_term"] = out["gamma_kpa_per_c"] / denom
    out["wind_factor"] = 900.0 / (t_mean + 273.0) * u2
    out["e0_tmax_kpa"] = fao_e0(t_max)
    out["e0_tmin_kpa"] = fao_e0(t_min)
    out["es_kpa"] = 0.5 * (out["e0_tmax_kpa"] + out["e0_tmin_kpa"])
    out["ea_kpa"] = 0.5 * (out["e0_tmin_kpa"] * rh_max / 100.0 +
                           out["e0_tmax_kpa"] * rh_min / 100.0)   # Eq. 17
    out["deficit_kpa"] = out["es_kpa"] - out["ea_kpa"]
    out["n_over_bign"] = n_sun / big_n
    out["rs_mj"] = (0.25 + 0.50 * out["n_over_bign"]) * ra_mj     # Eq. 35
    out["rso_mj"] = (0.75 + 2.0e-5 * z) * ra_mj                   # Eq. 37
    out["rs_over_rso"] = out["rs_mj"] / out["rso_mj"]
    out["rns_mj"] = 0.77 * out["rs_mj"]                           # Eq. 38
    sig_max = STEFAN_DAY * (t_max + 273.16) ** 4
    sig_min = STEFAN_DAY * (t_min + 273.16) ** 4
    out["sigma_tmax4_mj"] = sig_max
    out["sigma_tmin4_mj"] = sig_min
    out["sigma_t4_mean_mj"] = 0.5 * (sig_max + sig_min)
    out["emissivity_term"] = 0.34 - 0.14 * math.sqrt(out["ea_kpa"])
    out["cloud_term"] = 1.35 * out["rs_over_rso"] - 0.35
    out["rnl_mj"] = (out["sigma_t4_mean_mj"] * out["emissivity_term"] *
                     out["cloud_term"])                           # Eq. 39
    out["rn_mj"] = out["rns_mj"] - out["rnl_mj"]                  # Eq. 40
    out["g_mj"] = 0.0                                             # Eq. 42
    out["radiation_mm"] = 0.408 * (out["rn_mj"] - out["g_mj"])
    out["radiation_term_mm"] = out["radiation_mm"] * out["delta_term"]
    out["wind_term_mm"] = (out["gamma_term"] * out["wind_factor"] *
                           out["deficit_kpa"])
    out["et0_mm_per_day"] = out["radiation_term_mm"] + out["wind_term_mm"]
    return out


def fao56_hourly(t_hr: float, rh: float, u2: float, rs_mj: float, z: float,
                 ra_mj: float, rs_over_rso_assumed=None, daytime=True) -> dict:
    """The Example-19 chain (hourly, Eq. 53; G = 0.1 Rn day / 0.5 Rn
    night; sigma hourly)."""
    out = {}
    out["pressure_kpa"] = fao_pressure(z)
    out["delta_kpa_per_c"] = fao_delta(t_hr)
    out["gamma_kpa_per_c"] = fao_gamma(out["pressure_kpa"])
    out["e0_kpa"] = fao_e0(t_hr)
    out["ea_kpa"] = out["e0_kpa"] * rh / 100.0                    # Eq. 54
    out["deficit_kpa"] = out["e0_kpa"] - out["ea_kpa"]
    out["rso_mj_per_h"] = (0.75 + 2.0e-5 * z) * ra_mj if ra_mj > 0 else 0.0
    out["rns_mj_per_h"] = 0.77 * rs_mj
    out["sigma_tk4_mj_per_h"] = STEFAN_HOUR * (t_hr + 273.16) ** 4
    out["emissivity_term"] = 0.34 - 0.14 * math.sqrt(out["ea_kpa"])
    rr = (rs_over_rso_assumed if rs_over_rso_assumed is not None
          else rs_mj / out["rso_mj_per_h"])
    out["rs_over_rso"] = rr
    out["cloud_term"] = 1.35 * rr - 0.35
    out["rnl_mj_per_h"] = (out["sigma_tk4_mj_per_h"] * out["emissivity_term"] *
                           out["cloud_term"])
    out["rn_mj_per_h"] = out["rns_mj_per_h"] - out["rnl_mj_per_h"]
    out["g_mj_per_h"] = (0.1 if daytime else 0.5) * out["rn_mj_per_h"]
    out["rn_minus_g_mj"] = out["rn_mj_per_h"] - out["g_mj_per_h"]
    out["radiation_mm"] = 0.408 * out["rn_minus_g_mj"]
    denom = (out["delta_kpa_per_c"] +
             out["gamma_kpa_per_c"] * (1.0 + 0.34 * u2))
    out["radiation_term_mm"] = (out["radiation_mm"] *
                                out["delta_kpa_per_c"] / denom)
    out["wind_term_mm"] = (out["gamma_kpa_per_c"] / denom * 37.0 /
                           (t_hr + 273.0) * u2 * out["deficit_kpa"])
    out["et0_mm_per_h"] = out["radiation_term_mm"] + out["wind_term_mm"]
    return out


# --- The model's q-form (W7.2, mirrored term for term) ---------------------

AIR_HEAT_CAPACITY = 1005.0       # c_pa [J/kg/K] (BulkAerodynamic.hpp)
WATER_EMISSIVITY = 0.97          # kWaterEmissivity
STEFAN_SI = 5.670374419e-8       # [W m-2 K-4]


def latent_heat(t_c: float) -> float:
    return 2.501e6 - 2370.0 * t_c


def dqsat_dt(t_c: float, p_kpa: float) -> float:
    """Analytic d(q_sat)/dT of the model's Tetens chain [1/K]."""
    e0 = tetens(t_c)
    de = e0 * 17.27 * 237.3 / (t_c + 237.3) ** 2
    denom = p_kpa - 0.376 * e0
    return 0.622 * de * (p_kpa / (denom * denom))


def pm_q_form(t_air_c: float, p_kpa: float, q_air: float, rs_wm2: float,
              lw_in_wm2: float, albedo: float, g_frac: float, r_a: float,
              r_s: float) -> float:
    """The W7.2 formula [m/s of liquid water]:
      lambda rho_w E = [s A + lambda rho_a (q_sat - q_a)/r_a]
                       / [s + 1 + r_s/r_a]
    with s = lambda dq_sat/dT / c_pa at T_air and A = (1-albedo) R_s +
    eps (LW_in - sigma T_K^4) - G, G = g_frac * A_net-radiation."""
    lam = latent_heat(t_air_c)
    s = lam * dqsat_dt(t_air_c, p_kpa) / AIR_HEAT_CAPACITY
    rho_a = air_density(t_air_c, p_kpa)
    t_k = t_air_c + 273.15
    rn = ((1.0 - albedo) * rs_wm2 +
          WATER_EMISSIVITY * (lw_in_wm2 - STEFAN_SI * t_k ** 4))
    a = rn - g_frac * rn
    q_def = saturated_specific_humidity(t_air_c, p_kpa) - q_air
    num = s * a + lam * rho_a * q_def / r_a
    den = s + 1.0 + r_s / r_a
    return num / den / (lam * WATER_DENSITY)

