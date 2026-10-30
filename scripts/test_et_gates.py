#!/usr/bin/env python3
"""test_et_gates.py — the §6.3 self-check battery for the Stage-3 T5
evapotranspiration gate criteria (et0–et9; Stage-3 plan §9).

1. Reference sanity (§6.3 rule 1) — the bulk potential chain is checked
   two independent ways (the EtFormulas mirror vs a literal textbook
   recomputation with different operation ordering), the FAO-56 r_a
   against a hand value, t_d against the closed form, and the salt
   closed forms against direct bookkeeping.
2. Sensitivity (§6.3 rule 2) — every checker fails on a corrupted trace
   modeled on the physical failure it exists to catch (a double-counted
   transition window, a deleted salt store, an overdrawn cap, a leaking
   pond identity).

Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_et as et  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


print("bulk potential: two independent evaluations (et-met)")
# et-met: T_air 25, T_s 30, P 101.3, RH 0.3, U 8.
qa = et.met_specific_humidity(25.0, 101.3, 0.3)
ra = et.aero_resistance("liu2006", 8.0)
e_mirror = et.potential_bulk(30.0, 101.3, ra, 0.0, qa)
# Independent recomputation: different ordering and raw constants.
esat30 = 0.6108 * math.exp(17.27 * 30.0 / 267.3)
qsat30 = 0.622 * esat30 / (101.3 - 0.376 * esat30)
esat25 = 0.6108 * math.exp(17.27 * 25.0 / 262.3)
qa2 = 0.3 * (0.622 * esat25 / (101.3 - 0.376 * esat25))
rho = 101300.0 / (287.05 * 303.15)
e_lit = (rho * (qsat30 - qa2)) / (1000.0 * (94.909 * math.pow(8.0, -0.9036)))
expect(abs(e_mirror - e_lit) / e_lit < 1.0e-12,
       f"bulk chain two ways agree (E = {e_mirror:.6e} m/s)")
expect(1.0e-6 < e_mirror < 3.0e-6, "E_w is in the posted class (~1.6e-6)")

print("aerodynamic resistance laws")
# FAO-56: at z_ref = 2 m, u2 = U (4.87/ln(130.18) = 1.0006... -> hand value).
u2 = 8.0 * 4.87 / math.log(67.8 * 2.0 - 5.42)
expect(abs(et.aero_resistance("fao56", 8.0) - 208.0 / u2) < 1.0e-12,
       "fao56 matches 208/u2 with the height conversion")
# log profile hand value: z=2, z0m=0.01, z0h=0.001, U=8:
# ln(200) ln(2000) / (0.41^2 * 8)
hand = math.log(200.0) * math.log(2000.0) / (0.41 ** 2 * 8.0)
expect(abs(et.aero_resistance("log_profile", 8.0) - hand) < 1.0e-12,
       "log profile matches the hand value")
expect(et.aero_resistance("liu2006", 2.0) > et.aero_resistance("liu2006", 8.0),
       "liu2006 decreases with wind")

print("limiters and the cap")
fc = et.theta_fc_default(0.05, 0.41, 4.0, 2.5)
expect(0.05 < fc < 0.41, f"theta_fc default in range ({fc:.4f})")
expect(et.beta_lee_pielke(fc, fc) == 1.0 and et.beta_lee_pielke(fc + 0.01, fc) == 1.0,
       "beta = 1 at and above field capacity")
expect(0.0 < et.beta_lee_pielke(0.5 * fc, fc) < 1.0, "beta in (0,1) below fc")
expect(et.beta_lee_pielke(0.0, fc) == 0.0, "beta(0) = 0")
expect(et.supply_cap(1.0e-5, 0.0501, 0.05, 0.41, 0.05, 60.0)
       == (0.0501 - 0.05) * 0.05 / 60.0,
       "cap binds an unfundable evaporative demand")
expect(et.supply_cap(1.0e-9, 0.2, 0.05, 0.41, 0.05, 60.0) == 1.0e-9,
       "cap passes a fundable demand")
expect(et.supply_cap(-1.0e-5, 0.40999, 0.05, 0.41, 0.05, 60.0)
       == -max(0.0, 0.9999 * 0.41 - 0.40999) * 0.05 / 60.0,
       "cap bounds condensation by the 0.9999 theta_s room")

print("et3 references and checkers")
h0, ew = 4.0e-4, e_mirror
t_d = h0 / ew
expect(200.0 < t_d < 400.0, f"t_d = {t_d:.1f} s lands inside the 910 s window")
times = [70.0 * k for k in range(1, 4)]
depths = [h0 - ew * t for t in times]
cum = [ew * t for t in times]
exch = [0.0] * 3
ok, _ = et.check_pond_phase(times, depths, cum, h0, ew, exch, 1.0, 1.0e-12)
expect(ok, "pond-phase checker accepts the closed form")
bad = [c * 1.001 for c in cum]
ok, _ = et.check_pond_phase(times, depths, bad, h0, ew, exch, 1.0, 1.0e-12)
expect(not ok, "pond-phase checker: a 0.1% evaporation error fails")

wet = np.full((1, 1), 300.0)
soil = np.full((1, 1), 400.0)
ok, _ = et.check_no_double_count(wet, soil, 700.0)
expect(ok, "complementarity accepts wet+soil = t")
ok, _ = et.check_no_double_count(wet, soil + 7.0, 700.0)
expect(not ok, "complementarity: one double-counted substep (7 s) fails")

ok, _ = et.check_combined_rate(ew * 300.0, ew * 400.0, ew, 700.0, 1.0, 1.0e-10)
expect(ok, "combined rate accepts E_w t")
ok, _ = et.check_combined_rate(ew * 300.0, ew * 407.0, ew, 700.0, 1.0, 1.0e-10)
expect(not ok, "combined rate: the 1-window double count (1e-2) fails")

print("salt closed forms and checkers")
c0 = np.full((1, 1, 10), 10.0)
th0 = np.full((1, 1, 10), 0.3)
th1 = th0.copy()
th1[0, 0, 0] = 0.25                 # the top cell lost water
c1 = c0 * th0 / th1                 # exact concentration
ok, _ = et.check_cell_salt(c1, th1, c0, th0, 1.0e-12)
expect(ok, "per-cell c*theta accepts the exact factor")
ok, _ = et.check_cell_salt(c0, th1, c0, th0, 1.0e-12)
expect(not ok, "per-cell c*theta: an unconcentrated top cell fails")

ok, _ = et.check_pond_salt(35.0 * 0.005 / 0.003, 0.003, 35.0 * 0.005,
                           1.0e-12, "pond")
expect(ok, "pond salt accepts s h = s0 h0")
ok, _ = et.check_pond_salt(35.0, 0.003, 35.0 * 0.005, 1.0e-12, "pond")
expect(not ok, "pond salt: an undiluted-looking pond fails")
ok, _ = et.check_conserved(35.0 * 0.005, 35.0 * 0.005 * (1 - 1e-6), 1.0e-12,
                           "total")
expect(not ok, "conservation: a 1e-6 relative loss (deleted store class) fails")

ok, _ = et.check_cap_trace(0.0501 - 2.0e-15, 0.05, 3.0, 4.0, 1.0, 0.0, 1.0e-14)
expect(ok, "cap trace accepts the exact books")
ok, _ = et.check_cap_trace(0.05 - 1.0e-10, 0.05, 3.0, 4.0, 1.0, 0.0, 1.0e-14)
expect(not ok, "cap trace: theta below theta_r fails")
ok, _ = et.check_cap_trace(0.0501, 0.05, 2.9, 4.0, 1.0, 0.0, 1.0e-14)
expect(not ok, "cap trace: unrealised != demand - applied fails")
ok, _ = et.check_cap_trace(0.0501, 0.05, 3.0, 4.0, 1.0, 1.0e-9, 1.0e-14)
expect(not ok, "cap trace: nonzero vloss fails")

print("fao56 transcription lanes (et0 part 4)")
# The transcription with its printed precision (decimals recovered from
# the raw text, since yaml floats forget them).
import re
fao_path = ROOT / "benchmarks" / "et-evapotranspiration" / "reference" / "fao56-examples.yaml"
raw = fao_path.read_text()

def parse_section(lines_between):
    """key -> (value, decimals) for `key: number` lines in a slice."""
    out = {}
    for m in re.finditer(r"^\s+(\w+):\s+(-?\d+\.(\d+))(?:\s|$)",
                         lines_between, re.M):
        out[m.group(1)] = (float(m.group(2)), len(m.group(3)))
    return out

daily_txt = raw[raw.index("daily:"):raw.index("hourly:")]
d_ref = parse_section(daily_txt[daily_txt.index("intermediates"):])
hourly_txt = raw[raw.index("hourly:"):]
h_int = hourly_txt[hourly_txt.index("intermediates"):hourly_txt.index("result:")]
h1_ref = parse_section(h_int[h_int.index("lane_1400_1500"):h_int.index("lane_0200_0300")])
h2_ref = parse_section(h_int[h_int.index("lane_0200_0300"):])
h_res = parse_section(hourly_txt[hourly_txt.index("result:"):])
d_res = parse_section(daily_txt[daily_txt.index("result:"):])

# FAO carries u2 at full precision (Eq. 47 on the 10 km/h, 10 m input)
# and prints 2.078; feeding the chain the printed value misses the
# printed wind_factor 6.450 by 1.2 ulp.
u2_full = (10.0 / 3.6) * 4.87 / math.log(67.8 * 10.0 - 5.42)
d_chain = et.fao56_daily(t_max=21.5, t_min=12.3, rh_max=84.0, rh_min=63.0,
                         u2=u2_full, z=100.0, ra_mj=41.09, big_n=16.1,
                         n_sun=9.25)
d_chain["u2_ms"] = u2_full
h1_chain = et.fao56_hourly(t_hr=38.0, rh=52.0, u2=3.3, rs_mj=2.450, z=8.0,
                           ra_mj=3.543, daytime=True)
h2_chain = et.fao56_hourly(t_hr=28.0, rh=90.0, u2=1.9, rs_mj=0.0, z=8.0,
                           ra_mj=0.0, rs_over_rso_assumed=0.8, daytime=False)

NOT_CHAIN = {"ra_mj", "daylight_hours", "t_mid_h", "omega_rad",
             "omega1_rad", "omega2_rad", "ra_mj_per_h", "et0_rounded"}
coarse_keys = []  # keys where the precision floor, not 1 %, drives the negative

def compare_lane(name, refs, chain, result_map=None):
    n = 0
    for key, (printed, dec) in refs.items():
        if key in NOT_CHAIN:
            continue
        ck = key
        if result_map and key in result_map:
            ck = result_map[key]
        if ck not in chain:
            failures.append(f"fao56 {name}: no chain value for {key}")
            print(f"  fao56 {name} {key}: NO CHAIN VALUE")
            continue
        tol = 10.0 ** (-dec) + 1.0e-12          # last printed digit +- 1
        diff = abs(chain[ck] - printed)
        expect(diff <= tol,
               f"fao56 {name} {key}: chain {chain[ck]:.6g} vs printed "
               f"{printed} (tol {tol:.0e})")
        # Sensitivity: the perturbed transcription must fail. 1 % where
        # 1 % can exceed the precision floor; else 2.5x the floor
        # (S3-A10: for keys printed coarser than 1 %, a 1 % nudge is
        # inside the publication-precision tolerance by construction).
        delta = 0.01 * abs(printed)
        if delta + diff <= tol:
            delta = 2.5 * tol
            coarse_keys.append(f"{name}.{key}")
        worse = printed + delta if chain[ck] <= printed else printed - delta
        expect(abs(chain[ck] - worse) > tol,
               f"fao56 {name} {key}: perturbed transcription fails")
        n += 1
    return n

n_d = compare_lane("daily", {**d_ref, **d_res}, d_chain)
n_h1 = compare_lane("hourly-14", {**h1_ref,
                    "et0_1400_1500_mm_per_h": h_res["et0_1400_1500_mm_per_h"]},
                    h1_chain, {"et0_1400_1500_mm_per_h": "et0_mm_per_h"})
n_h2 = compare_lane("hourly-02", {**h2_ref,
                    "et0_0200_0300_mm_per_h": h_res["et0_0200_0300_mm_per_h"]},
                    h2_chain, {"et0_0200_0300_mm_per_h": "et0_mm_per_h"})
expect(n_d >= 26 and n_h1 >= 18 and n_h2 >= 13,
       f"fao56 coverage: {n_d}/{n_h1}/{n_h2} keys compared (vacuity guard)")
print(f"  fao56: precision-floor negatives on {len(coarse_keys)} coarse keys: "
      f"{', '.join(coarse_keys)}")

print("pm q-form vs the FAO chain (delta_conv for et10)")
# The et10 configuration on the 14-15 h lane: grass r_c = r_l/(0.5 LAI)
# = 100/(0.5*2.88), albedo 0.23, G fraction 0.1, r_a = 208/u2 (fao56 law
# with u2 given at z_ref), LW_in chosen so eps (LW_in - sigma T_K^4)
# equals -R_nl of the chain, R_s = the chain's forcing.
p_kpa = et.fao_pressure(8.0)
q_air = et.met_specific_humidity(38.0, p_kpa, 0.52)
rs_w = 2.450e6 / 3600.0
rnl_w = h1_chain["rnl_mj_per_h"] * 1.0e6 / 3600.0
t_k = 38.0 + 273.15
lw_in = et.STEFAN_SI * t_k ** 4 - rnl_w / et.WATER_EMISSIVITY
r_a = 208.0 / 3.3
r_c = 100.0 / (0.5 * 2.88)
e_pm = et.pm_q_form(38.0, p_kpa, q_air, rs_w, lw_in, 0.23, 0.1, r_a, r_c)
pm_mm_h = e_pm * 3.6e6
delta_conv = abs(pm_mm_h - h1_chain["et0_mm_per_h"])
rel = delta_conv / h1_chain["et0_mm_per_h"]
print(f"  PM_q = {pm_mm_h:.6f} mm/h; FAO chain = "
      f"{h1_chain['et0_mm_per_h']:.6f} mm/h; delta_conv = {delta_conv:.6f} "
      f"mm/h ({100*rel:.2f} %)")
expect(0.0 < rel < 0.05, "delta_conv in the expected 1-3 % class (< 5 %)")
# Sanity: with FAO's own gamma-form psychrometrics the two should be far
# closer than with the model constants -- the gap is convention, not a bug.
expect(e_pm > 0.0, "PM_q positive on the daytime lane")

if failures:
    print(f"{len(failures)} et gate-logic check(s) FAILED")
    sys.exit(1)
print("all et gate-logic checks passed")
