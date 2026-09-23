#!/usr/bin/env python3
"""test_g45_gates.py — the §6.3 self-check battery for the v2 Q4 gate
criteria (g4/g5, plan §3.3 as amended by V2-A13/V2-A14).

Two obligations, both discharged on synthetic data (the capability the
gates gate does not exist yet — this battery is green while the gates
themselves are red, the same split as scripts/test_p_gates.py):

1. Reference sanity (§6.3 rule 1): the closed forms are checked against
   independently computed spot values, and the digitized reference curves
   themselves must PASS the criteria derived from them (a reference that
   fails its own gate is the V2-A13 defect).
2. Sensitivity (§6.3 rule 2): every criterion demonstrably FAILS when fed
   a corrupted input modeled on the physical failure it exists to catch
   (limiter missing -> rate does not decay; wrong alpha_1 -> wrong E;
   profile offset; salt leak; no density-driven deepening).

Exits nonzero on the first violated expectation.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_evap  # noqa: E402

REF = ROOT / "benchmarks" / "g5-geng2015" / "reference"

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    status = "ok" if condition else "FAIL"
    print(f"  {label}: {status}")
    if not condition:
        failures.append(label)


def read_csv(path: pathlib.Path) -> dict[str, np.ndarray]:
    rows = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")]
    header = rows[0].split(",")
    data = np.array([[float(v) if v else np.nan for v in line.split(",")]
                     for line in rows[1:]])
    return {name: data[:, k] for k, name in enumerate(header)}


# ---------------------------------------------------------------------------
print("closed forms (two independent spot values each, §6.3 rule 1):")
# Tetens at 20 C: 0.6108*exp(17.27*20/257.3) = 0.6108*exp(1.342383...)
#   hand value 2.3383 kPa; and at 25 C: 0.6108*exp(431.75/262.3) = 3.1678 kPa.
expect(abs(gate_evap.tetens_esat_kpa(20.0) - 2.3383) < 5.0e-4, "e_sat(20C) = 2.3383 kPa")
expect(abs(gate_evap.tetens_esat_kpa(25.0) - 3.1678) < 5.0e-4, "e_sat(25C) = 3.1678 kPa")
# q_sat(20C, 101.325 kPa) = 0.622*2.3383/(101.325 - 0.376*2.3383) = 0.014480;
# Table 1's own cross-check: q_a = 20% of q_sat = 2.896e-3 ("2.9e-3").
qsat = gate_evap.saturated_specific_humidity(20.0, 101.325)
expect(abs(qsat - 0.014480) < 2.0e-6, "q_sat(20C) = 0.014480")
expect(abs(0.2 * qsat - 2.9e-3) < 1.0e-5, "q_a = 0.2 q_sat matches Table 1's 2.9e-3")
# R_air(1 m/s) = 94.909 s/m (the fit's own coefficient); R_air(2) = 94.909/2^0.9036.
expect(abs(gate_evap.aerodynamic_resistance(1.0) - 94.909) < 1.0e-9, "R_air(1) = 94.909")
expect(abs(gate_evap.aerodynamic_resistance(2.0) - 50.734) < 5.0e-3, "R_air(2) = 50.734")
# E(0): rho_a = 101325/(287.05*293.15) = 1.2041 kg/m^3;
#   (1.2041/94.909)*(0.014480 - 0.0028959)/1000 = 1.4696e-7 m/s.
e0 = gate_evap.bulk_evaporation_rate(20.0, 101.325, 1.0, 0.2 * qsat)
expect(abs(e0 - 1.4696e-7) < 2.0e-10, "E(0) = 1.4696e-7 m/s (Table-1 chain)")
# alpha_1: equilibrium inversion w_g = 0.0375 -> 0.2; saturated 0.41 -> 1.
expect(abs(gate_evap.soil_relative_humidity(0.0375) - 0.2) < 1.0e-3,
       "alpha_1(0.0375) = 0.2 (the V2-A13 equilibrium anchor)")
expect(gate_evap.soil_relative_humidity(0.41) == 1.0, "alpha_1(0.41) capped at 1")

# ---------------------------------------------------------------------------
print("g4 drawdown:")
t = np.arange(0.0, 50001.0, 5000.0)
rate = 1.0e-6
exact_err = np.zeros_like(t)
ok, _ = gate_evap.check_drawdown(t, exact_err, 0.0, rate, 0.05, 1.0e-6)
expect(ok, "exact drawdown passes")
# A 0.1 % rate error accumulates to 5e-5 m at the end: >> 5e-8 allowance.
bad_err = np.abs((-1.001e-6 * t) - (-rate * t))
ok, _ = gate_evap.check_drawdown(t, bad_err, 0.0, rate, 0.05, 1.0e-6)
expect(not ok, "0.1% rate error fails")

print("g4 concentration:")
depth = 0.1 - rate * t
s_exact = 10.0 * 0.1 / depth
mass = np.full(t.size, 12.0)
ok, _ = gate_evap.check_concentration(t, s_exact, depth, 10.0, 0.1, 1.0e-4,
                                      mass, 1.0e-10)
expect(ok, "exact concentration passes")
# The pre-Q4 throttle analogue: concentration pinned at s0.
ok, _ = gate_evap.check_concentration(t, np.full(t.size, 10.0), depth, 10.0,
                                      0.1, 1.0e-4, mass, 1.0e-10)
expect(not ok, "throttled (constant) concentration fails")
mass_leak = mass * (1.0 - 1.0e-6 * np.arange(t.size))
ok, _ = gate_evap.check_concentration(t, s_exact, depth, 10.0, 0.1, 1.0e-4,
                                      mass_leak, 1.0e-10)
expect(not ok, "salt-mass leak fails")

print("g4 closure and dry-out:")
ok, _ = gate_evap.check_closure(-0.24, 0.0, 0.24, 0.0, 0.0, 0.24, 1.0e-8)
expect(ok, "closed ledger passes")
# The pre-Q4 g4(d) state: potential evap booked for the dry half of the run.
ok, _ = gate_evap.check_closure(-0.24, 0.0, 0.48, 0.0, 0.0, 0.24, 1.0e-8)
expect(not ok, "unaudited dry-out shortfall fails")
ok, _ = gate_evap.check_positivity_and_dry([0.02, 0.001, 0.0], 0.0, 1.0e-8)
expect(ok, "nonnegative + dry passes")
ok, _ = gate_evap.check_positivity_and_dry([0.02, -1.0e-6, 0.0], 0.0, 1.0e-8)
expect(not ok, "negative depth fails")
ok, _ = gate_evap.check_positivity_and_dry([0.02, 0.01, 0.005], 0.005, 1.0e-8)
expect(not ok, "not drying out fails")

# ---------------------------------------------------------------------------
print("g5(i) rate criteria (V2-A15) — the digitized reference passes its own gate:")
fig3 = read_csv(REF / "fig3_evaporation_rate.csv")
ok, msgs = gate_evap.check_rate_series(fig3["time_h"], fig3["evaporation_m_per_s"],
                                       1.470e-7, 3.317e-8, 0.109)
for msg in msgs:
    print(f"    {msg}")
expect(ok, "digitized Fig. 3 + equilibrium-band surface value pass")
# Limiter-missing analogue: E stays at its potential value -> the decay
# bound (0.5) and monotonicity both fail.
flat = np.full(fig3["time_h"].size, 1.47e-7)
ok, _ = gate_evap.check_rate_series(fig3["time_h"], flat, 1.470e-7, 3.317e-8, 0.109)
expect(not ok, "undecaying (no alpha_1 limiter) rate fails")
# Wrong-formula analogue: E(0) off by 10 %.
ok, _ = gate_evap.check_rate_series(fig3["time_h"],
                                    0.9 * fig3["evaporation_m_per_s"],
                                    1.470e-7, 3.317e-8, 0.109)
expect(not ok, "10% E(0) formula error fails")
# Sign-error analogue: the surface saturation crosses the closed-form
# alpha_1 equilibrium floor (condensation would have to be feeding it).
ok, _ = gate_evap.check_rate_series(fig3["time_h"], fig3["evaporation_m_per_s"],
                                    1.470e-7, 3.317e-8, 0.05)
expect(not ok, "surface saturation below the equilibrium floor fails")
# Non-monotone analogue: an oscillating tail fails.
wob = fig3["evaporation_m_per_s"] * (1.0 + 0.05 * np.sign(np.sin(fig3["time_h"] * 2.0)))
ok, _ = gate_evap.check_rate_series(fig3["time_h"], wob, 1.470e-7, 3.317e-8, 0.109)
expect(not ok, "oscillating rate history fails monotonicity")

print("g5(ii) — the salinization peak gates, the profile RMS records (V2-A15):")
fig4 = read_csv(REF / "fig4_moisture_ratio.csv")
fig9 = read_csv(REF / "fig9b_salinity.csv")
z = fig4["elevation_m"]
ok, msgs = gate_evap.record_moisture_profiles(z, fig4["moisture_20h"],
                                              fig4["moisture_50h"],
                                              fig4["moisture_20h"],
                                              fig4["moisture_50h"])
expect(ok and all("not gated" in m for m in msgs),
       "moisture RMS is recorded, never gated")
ok, _ = gate_evap.check_salinity_profiles(z, fig9["salinity_20h_gL"],
                                          fig9["salinity_50h_gL"],
                                          fig9["salinity_20h_gL"],
                                          fig9["salinity_50h_gL"],
                                          60.0, 1.9)
expect(ok, "digitized salinity passes (the >60 g/L peak present)")
ok, _ = gate_evap.check_salinity_profiles(z, fig9["salinity_20h_gL"],
                                          np.minimum(fig9["salinity_50h_gL"], 55.0),
                                          fig9["salinity_20h_gL"],
                                          fig9["salinity_50h_gL"],
                                          60.0, 1.9)
expect(not ok, "missing near-surface peak (no scalar_cauchy analogue) fails")

print("g5(iii) salt mass:")
ok, _ = gate_evap.check_salt_mass([100.0, 99.0, 98.5, 99.5], 0.04)
expect(ok, "1.5% drift passes (bound 4%)")
ok, _ = gate_evap.check_salt_mass([100.0, 97.0, 94.0], 0.04)
expect(not ok, "6% drift fails")

print("g5(iv) density contrast (V2-A14 direction):")
zc = np.linspace(1.99, 1.61, 20)
plume = np.full((30, 20), 25.0)
plume[:, zc >= 1.86] = 40.0        # control: edge at ~1.86 m
plume_deep = plume.copy()
plume_deep[3, zc >= 1.70] = 40.0   # beta run: one finger to ~1.70 m
zb = gate_evap.plume_edge_min_elevation(zc, plume_deep, 30.0)
zctl = gate_evap.plume_edge_min_elevation(zc, plume, 30.0)
ok, _ = gate_evap.check_density_contrast(zb, zctl, 0.05)
expect(ok, "finger-deepened beta run vs flat control passes")
ok, _ = gate_evap.check_density_contrast(zctl, zctl, 0.05)
expect(not ok, "equal plume depths fail (the density effect absent)")
ok, _ = gate_evap.check_density_contrast(float("nan"), zctl, 0.05)
expect(not ok, "missing plume fails")

# ---------------------------------------------------------------------------
if failures:
    print(f"test_g45_gates: {len(failures)} FAILED expectation(s)")
    for name in failures:
        print(f"  - {name}")
    sys.exit(1)
print("test_g45_gates: OK")
