#!/usr/bin/env python3
"""test_g910_gates.py — the §6.3 self-check battery for the v2 Q6 wind
gate criteria (g9/g10 + the 8-orientation battery, plan §5.3).

1. Reference sanity (§6.3 rule 1) — each closed form checked two
   independent ways: the g9(a) analytic-integral volume constraint vs a
   numeric trapezoid quadrature; the g9(c) ODE quadrature vs the flat-
   bottom closed form; Merian and every Cd law vs hand values.
2. Sensitivity (§6.3 rule 2) — every check fails on a corrupted input
   modeled on the physical failure it exists to catch.

Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_wind  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


# ---------------------------------------------------------------------------
print("Cd(U10) laws vs hand values (published formulas):")
expect(abs(gate_wind.cd_garratt(15.0) - 1.755e-3) < 1e-15, "garratt(15) = 1.755e-3")
expect(abs(gate_wind.cd_garratt(45.0) - 3.5e-3) < 1e-15, "garratt(45) caps at 3.5e-3")
expect(abs(gate_wind.cd_smith_banke(10.0) - 1.29e-3) < 1e-15, "smith-banke(10) = 1.29e-3")
expect(abs(gate_wind.cd_wu(0.0) - 0.8e-3) < 1e-15, "wu(0) = 0.8e-3")
expect(abs(gate_wind.cd_large_pond(10.0) - 1.2e-3) < 1e-15, "large-pond(10) = 1.2e-3")
expect(abs(gate_wind.cd_large_pond(20.0) - 1.79e-3) < 1e-15, "large-pond(20) = 1.79e-3")
expect(abs(gate_wind.cd_large_pond(60.0) - 3.5e-3) < 1e-15, "large-pond(60) caps")
# kinematic stress: 1.225 * 1.755e-3 * 225 / 998 = 4.8469e-4 (hand).
expect(abs(gate_wind.kinematic_stress(1.755e-3, 15.0) - 4.8469e-4) < 1e-8,
       "tau'/rho_w hand value at Garratt(15)")

print("g9(a) closed form (two independent volume evaluations):")
taup = gate_wind.kinematic_stress(1.0277551, 5.0)
x = np.linspace(0.0, 500.0, 20001)
h = gate_wind.setup_flat(taup, 500.0, 2.0, x)
expect(abs(np.trapz(h, x) - 2.0 * 500.0) < 1e-3,
       "volume conserved (numeric quadrature agrees with the analytic integral)")
expect(abs(h[0] - 1.5643) < 2e-4, "h(0) = 1.5643 m (hand value)")
expect(abs(h[-1] - 2.3795) < 2e-4, "h(500) = 2.3795 m (hand value)")

print("g9(c) quadrature cross-check (flat bed reproduces the closed form):")
xq, eta_q = gate_wind.setup_slope(taup, 500.0, lambda xx: np.full_like(
    np.asarray(xx, dtype=float), -2.0), n=4000)
h_flat = gate_wind.setup_flat(taup, 500.0, 2.0, xq)
expect(float(np.max(np.abs((eta_q + 2.0) - h_flat))) < 1e-4,
       "quadrature vs closed form max dev < 1e-4 m")

print("Merian: T = 2*5000/sqrt(9.81*10) = 1009.66 s (hand):")
expect(abs(gate_wind.merian_period(5000.0, 10.0) - 1009.66) < 0.05, "Merian hand value")

# ---------------------------------------------------------------------------
print("sensitivity (each check fails on its physical failure):")
ok, _ = gate_wind.check_setup(0.8152, 0.8152, 0.005, "syn")
expect(ok, "exact setup passes")
ok, _ = gate_wind.check_setup(0.8152 * 1.02, 0.8152, 0.005, "syn")
expect(not ok, "2% setup error fails (0.5% bound)")

t = np.arange(0.0, 20000.0, 1.0)
ok, _ = gate_wind.check_steady(t, np.full(t.size, 0.4), 0.1, 2e-3, "syn")
expect(ok, "flat tail passes steadiness")
ok, _ = gate_wind.check_steady(t, 0.4 + 0.01 * np.sin(t / 200.0), 0.1, 2e-3, "syn")
expect(not ok, "sloshing tail fails steadiness")

ref = np.linspace(-0.4, 0.4, 50)
ok, _ = gate_wind.check_profile(ref, ref, 0.8, 0.01, "syn")
expect(ok, "exact profile passes")
ok, _ = gate_wind.check_profile(ref + 0.02, ref, 0.8, 0.01, "syn")
expect(not ok, "2.5%-of-range profile offset fails")
ok, _ = gate_wind.check_profile_rms(ref, ref, 0.8, 0.01, "syn")
expect(ok, "exact profile passes the RMS norm")
ok, _ = gate_wind.check_profile_rms(ref + 0.02, ref, 0.8, 0.01, "syn")
expect(not ok, "2.5%-of-range offset fails the RMS norm")

merian = gate_wind.merian_period(5000.0, 10.0)
tt = np.arange(30000.0, 38200.0, 5.0)
sig = 0.012 * np.cos(2 * np.pi * (tt - 30000.0) / merian) * np.exp(-(tt - 30000.0) / 40000.0)
period, amps = gate_wind.seiche_metrics(tt, sig, 30500.0)
ok, msgs = gate_wind.check_seiche(period, merian, 0.02, amps, 0.5)
for m in msgs:
    print(f"    {m}")
expect(ok, "synthetic Merian oscillation passes")
period5, _ = gate_wind.seiche_metrics(
    tt, 0.012 * np.cos(2 * np.pi * (tt - 30000.0) / (merian * 1.05)), 30500.0)
ok, _ = gate_wind.check_seiche(period5, merian, 0.02, None, 0.5)
expect(not ok, "5% period error fails (2% bound)")
ok, _ = gate_wind.check_seiche(float("nan"), merian, 0.02, None, 0.5)
expect(not ok, "no oscillation (dead forcing-series handling) fails")

print("orientation battery:")
rng = np.random.default_rng(7)
base = np.add.outer(np.zeros(20), np.linspace(0.0, 3.6e-3, 20))  # E-wind setup
diag = np.add.outer(np.linspace(0.0, 2.5e-3, 20), np.linspace(0.0, 2.5e-3, 20))
fields = {
    "E": base, "N": base.T.copy(), "W": base[:, ::-1].copy(),
    "S": base[:, ::-1].T.copy(),
    "NE": diag, "NW": diag[:, ::-1].copy(), "SW": diag[::-1, ::-1].copy(),
    "SE": diag[::-1, :].copy(),
}
ok, _ = gate_wind.check_orientations(fields, 1e-12)
expect(ok, "consistent orientation set passes")
bad = dict(fields)
bad["W"] = fields["W"] + 1e-6 * rng.standard_normal((20, 20))
ok, _ = gate_wind.check_orientations(bad, 1e-7)
expect(not ok, "a perturbed orientation fails")

# ---------------------------------------------------------------------------
if failures:
    print(f"test_g910_gates: {len(failures)} FAILED expectation(s)")
    for name in failures:
        print(f"  - {name}")
    sys.exit(1)
print("test_g910_gates: OK")
