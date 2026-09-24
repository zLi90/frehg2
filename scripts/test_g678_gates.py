#!/usr/bin/env python3
"""test_g678_gates.py — the §6.3 self-check battery for the v2 Q5
temperature gate criteria (g6/g7/g8, plan §4.2).

1. Reference sanity (§6.3 rule 1) — every closed form checked two
   independent ways: residual substitution into the governing equation,
   limiting-case recovery (Pe -> 0 linear; v_t = 0 classical damping
   depth; k -> 0 plain Ogata-Banks), and hand values.
2. Sensitivity (§6.3 rule 2) — every criterion fails on a corrupted
   input modeled on the physical failure it exists to catch (numerical
   dissipation -> amplitude loss; dispersion error -> phase error;
   broken exchange -> wrong T_e; dead density coupling -> no
   convection).

Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_heat  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


# ---------------------------------------------------------------------------
print("g6(a) Bredehoeft & Papadopulos:")
xi = np.linspace(0.0, 1.0, 101)
# Limit check: Pe -> 0 is linear.
expect(float(np.max(np.abs(gate_heat.bp_profile(xi, 1.0e-13) - xi))) < 1e-9,
       "Pe -> 0 recovers the linear profile")
# Residual check: alpha f'' - v f' = 0 with Pe = vL/alpha -> f'' = Pe f'
# (finite differences on the closed form).
for pe in (-5.0, -1.0, 1.0, 5.0):
    f = gate_heat.bp_profile(xi, pe)
    d1 = np.gradient(f, xi)
    d2 = np.gradient(d1, xi)
    resid = np.max(np.abs(d2[5:-5] - pe * d1[5:-5])) / np.max(np.abs(d2[5:-5]))
    expect(resid < 5e-3, f"Pe={pe:+g} closed form satisfies f'' = Pe f'")
# Hand value: Pe=5 at xi=0.5: (e^2.5-1)/(e^5-1) = 11.1825/147.413 = 0.07585.
expect(abs(float(gate_heat.bp_profile(np.array([0.5]), 5.0)[0]) - 0.075853) < 1e-5,
       "Pe=5 midpoint hand value 0.07585")

print("g6(b) Ogata-Banks heat form (OGS parameter set):")
alpha, vt = 1.1e-6, 1.5e-6
# Limit: x = 0 gives the full step immediately; t -> inf gives the step.
expect(abs(gate_heat.ogata_banks(0.0, 1.0, vt, alpha, 300.0, 330.0) - 330.0) < 1e-9,
       "inlet holds the step")
expect(abs(gate_heat.ogata_banks(1.0, 5e9, vt, alpha, 300.0, 330.0) - 330.0) < 1e-6,
       "long-time limit reaches the step")
# Hand value at the front: x = v t exactly -> erfc(0)/2 + small term.
t_at = 20.0 / vt  # front at x = 20 m
val = gate_heat.ogata_banks(20.0, t_at, vt, alpha, 300.0, 330.0)
# erfc(0)/2 of the step plus the second-term tail: hand value 316.59.
expect(abs(val - 316.592) < 0.05, "front midpoint hand value 316.59")
# Overflow guard: x = 50 m early time must not overflow and stay ~T0.
early = gate_heat.ogata_banks(50.0, 1e5, vt, alpha, 300.0, 330.0)
expect(300.0 <= early < 300.0001, "far field early time stays at T0 (guarded)")
# PDE residual on a grid (dT/dt + v dT/dx = a d2T/dx2).
xs = np.linspace(5.0, 30.0, 60)
tq = 1.5e7
h = 1.0e-3
T = np.array([gate_heat.ogata_banks(x, tq, vt, alpha, 0.0, 1.0) for x in xs])
Tt = np.array([(gate_heat.ogata_banks(x, tq * (1 + h), vt, alpha, 0.0, 1.0) -
                gate_heat.ogata_banks(x, tq * (1 - h), vt, alpha, 0.0, 1.0)) /
               (2 * tq * h) for x in xs])
Tx = np.gradient(T, xs)
Txx = np.gradient(Tx, xs)
resid = np.max(np.abs(Tt[5:-5] + vt * Tx[5:-5] - alpha * Txx[5:-5]))
scale = np.max(np.abs(vt * Tx))
expect(resid < 5e-3 * scale, "Ogata-Banks satisfies the ADE (residual check)")

print("g6(c) Stallman:")
omega = 2.0 * math.pi / 86400.0
alpha_c = 2.0 / 2.96e6  # lambda / (rho c)_bulk
a0, b0 = gate_heat.stallman_coefficients(0.0, alpha_c, omega)
classic = math.sqrt(omega / (2.0 * alpha_c))
expect(abs(a0 - classic) < 1e-12 and abs(b0 - classic) < 1e-12,
       "v_t = 0 recovers the classical damping depth")
for vt_c in (5e-6 * 4.184e6 / 2.96e6, -5e-6 * 4.184e6 / 2.96e6):
    a, b = gate_heat.stallman_coefficients(vt_c, alpha_c, omega)
    r1 = alpha_c * (a * a - b * b) + vt_c * a
    r2 = 2.0 * alpha_c * a * b + vt_c * b - omega
    expect(abs(r1) < 1e-12 * omega and abs(r2) < 1e-12 * omega,
           f"dispersion-relation residuals vanish (v_t = {vt_c:+.2e})")
# Downward flow carries heat down: smaller damping (larger amplitude);
# the PHASE is even in v (a(-v) = a(+v) + v/alpha, b identical) — the
# published reason the Hatch (2006) phase method cannot resolve flow
# DIRECTION, so the g6(c) amplitude criterion carries the direction
# discrimination and the phase criterion the speed.
amp_down, lag_down = gate_heat.stallman_amplitude_phase(0.2, +7.07e-6, alpha_c, 86400.0)
amp_up, lag_up = gate_heat.stallman_amplitude_phase(0.2, -7.07e-6, alpha_c, 86400.0)
_, lag_still = gate_heat.stallman_amplitude_phase(0.2, 0.0, alpha_c, 86400.0)
expect(amp_down > amp_up, "downward flow deepens the diel signal")
expect(abs(lag_down - lag_up) < 1e-9, "phase lag is even in v (direction-blind)")
expect(lag_down < lag_still, "flow of either sign shortens the lag")
a_up, _ = gate_heat.stallman_coefficients(-7.07e-6, alpha_c, omega)
a_dn, _ = gate_heat.stallman_coefficients(+7.07e-6, alpha_c, omega)
expect(abs((a_up - a_dn) - 7.07e-6 / alpha_c) < 1e-9,
       "a(-v) = a(+v) + v/alpha (the exact antisymmetry)")
# The fit machinery recovers a known sinusoid exactly.
t = np.arange(0.0, 3 * 86400.0, 600.0)
sig = 12.0 + 3.2 * np.sin(omega * (t - 4000.0))
amp, lag, mean = gate_heat.fit_sinusoid(t, sig, 86400.0)
expect(abs(amp - 3.2) < 1e-9 and abs(lag - 4000.0) < 1e-6 and abs(mean - 12.0) < 1e-9,
       "sinusoid fit recovers (A, lag, mean) exactly")

print("g7(a) Edinger:")
# Hand value: tau = 4.184e6 * 1 / 30 = 1.3947e5 s; at t = tau, excess/e.
tau = 4.184e6 / 30.0
val = float(gate_heat.edinger_relaxation(np.array([tau]), 30.0, 20.0, 30.0, 1.0)[0])
expect(abs(val - (20.0 + 10.0 / math.e)) < 1e-9, "T(tau) = Te + (T0-Te)/e")
# ODE residual: dT/dt = -K/(rho c h) (T - Te).
ts = np.linspace(0.0, 5 * tau, 4000)
T = gate_heat.edinger_relaxation(ts, 30.0, 20.0, 30.0, 1.0)
dT = np.gradient(T, ts)
resid = np.max(np.abs(dT[2:-2] + (30.0 / 4.184e6) * (T[2:-2] - 20.0)))
expect(resid < 1e-4 * float(np.max(np.abs(dT))),
       "Edinger closed form satisfies its ODE (interior FD residual)")

print("g7(b) channel plume:")
u, D, k = 0.5, 5.0, 25.0 / 4.184e6
# k -> 0 steady limit is flat.
expect(float(np.max(np.abs(gate_heat.channel_steady(np.linspace(0, 1e4, 50), u, D, 0.0) - 1.0))) < 1e-12,
       "k = 0 steady limit is uniform")
# ODE residual: u f' = D f'' - k f on the closed form.
xs = np.linspace(0.0, 3.0e5, 4000)
f = gate_heat.channel_steady(xs, u, D, k)
f1 = np.gradient(f, xs)
f2 = np.gradient(f1, xs)
resid = np.max(np.abs(u * f1[5:-5] - D * f2[5:-5] + k * f[5:-5])) / k
expect(resid < 5e-3, "channel steady form satisfies u f' = D f'' - k f")
# Transient long-time limit reaches the steady profile.
xm = 4.2e4
expect(abs(gate_heat.channel_transient(xm, 5e6, u, D, k) -
           float(gate_heat.channel_steady(np.array([xm]), u, D, k)[0])) < 1e-9,
       "transient long-time limit = steady profile")
# k = 0 transient reduces to plain Ogata-Banks (one harness, two configs).
ob = gate_heat.ogata_banks(xm, 9e4, u, D, 0.0, 1.0)
expect(abs(gate_heat.channel_transient(xm, 9e4, u, D, 0.0) - ob) < 1e-12,
       "k = 0 transient reduces to Ogata-Banks")

print("g8 Rayleigh:")
expect(abs(gate_heat.RAYLEIGH_CRITICAL - 39.478) < 1e-3, "Ra_c = 4 pi^2")
# Dimensional check: K [m/s] * beta dT [-] * H [m] / alpha [m2/s] = [-].
ra = gate_heat.rayleigh_number(1.0e-4, 2.0e-4, 10.0, 1.0, 5.0e-7)
expect(abs(ra - 0.4) < 1e-12, "Ra assembly hand value")

# ---------------------------------------------------------------------------
print("sensitivity (each criterion fails on its physical failure):")
ok, _ = gate_heat.check_bp([5.0], [0.05], 10.0, 0.01)
expect(ok, "B&P within bound passes")
ok, _ = gate_heat.check_bp([5.0], [0.2], 10.0, 0.01)
expect(not ok, "2% B&P profile error fails (1% of dT bound)")

ok, _ = gate_heat.check_breakthrough([20.0], [0.005], 0.01)
expect(ok, "breakthrough within bound passes")
ok, _ = gate_heat.check_breakthrough([20.0], [0.03], 0.01)
expect(not ok, "3% breakthrough error fails")

depths = [0.1, 0.2]
ar = [gate_heat.stallman_amplitude_phase(z, 0.0, alpha_c, 86400.0)[0] for z in depths]
lr = [gate_heat.stallman_amplitude_phase(z, 0.0, alpha_c, 86400.0)[1] for z in depths]
ok, _ = gate_heat.check_stallman(depths, ar, ar, lr, lr, 0.05, 600.0, "syn")
expect(ok, "exact Stallman passes")
ok, _ = gate_heat.check_stallman(depths, [a * 0.9 for a in ar], ar, lr, lr,
                                 0.05, 600.0, "syn")
expect(not ok, "10% amplitude loss (numerical dissipation) fails")
ok, _ = gate_heat.check_stallman(depths, ar, ar, [l + 900 for l in lr], lr,
                                 0.05, 600.0, "syn")
expect(not ok, "15-minute phase error (dispersion) fails")

ts = np.linspace(0.0, 10 * 86400.0, 200)
exact = gate_heat.edinger_relaxation(ts, 30.0, 20.0, 30.0, 1.0)
ok, _ = gate_heat.check_relaxation(ts, exact, 30.0, 20.0, 30.0, 1.0, 0.005)
expect(ok, "exact relaxation passes")
ok, _ = gate_heat.check_relaxation(ts, gate_heat.edinger_relaxation(ts, 30.0, 20.0, 27.0, 1.0),
                                   30.0, 20.0, 30.0, 1.0, 0.005)
expect(not ok, "10% exchange-coefficient error fails")
ok, _ = gate_heat.check_energy_ledger([1e-9, -5e-9], 1e-8)
expect(ok, "closed energy ledger passes")
ok, _ = gate_heat.check_energy_ledger([1e-9, 5e-7], 1e-8)
expect(not ok, "leaking energy ledger fails")
ok, _ = gate_heat.check_equilibrium_consistency(21.32, 21.30, 0.05)
expect(ok, "bulk equilibrium within 0.05 K passes")
ok, _ = gate_heat.check_equilibrium_consistency(21.45, 21.30, 0.05)
expect(not ok, "0.15 K equilibrium inconsistency fails")

xs = np.linspace(0.0, 8.4e4, 160)
sr = gate_heat.channel_steady(xs, u, D, k)
tt = np.linspace(1e3, 2.5e5, 60)
tr = np.array([gate_heat.channel_transient(4.2e4, t_, u, D, k) for t_ in tt])
ok, _ = gate_heat.check_channel(xs, sr, sr, 0.01, tr, tr, 0.02, 1.0)
expect(ok, "exact channel passes")
ok, _ = gate_heat.check_channel(xs, sr * 0.97, sr, 0.01, tr, tr, 0.02, 1.0)
expect(not ok, "3% steady offset fails the 1% L2")
ok, _ = gate_heat.check_channel(xs, sr, sr, 0.01, np.clip(tr - 0.05, 0, 1), tr, 0.02, 1.0)
expect(not ok, "5%-of-step transient error fails")

ok, _ = gate_heat.check_hrl([1.0, 0.5, 0.05], [1.2, 1.3, 1.25, 1.3], 0.2, 1.05)
expect(ok, "decaying subcritical + convecting supercritical passes")
ok, _ = gate_heat.check_hrl([1.0, 1.1, 1.2], [1.2, 1.3, 1.25, 1.3], 0.2, 1.05)
expect(not ok, "growing subcritical perturbation fails")
ok, _ = gate_heat.check_hrl([1.0, 0.5, 0.05], [1.2, 1.3, 1.01, 1.02], 0.2, 1.05)
expect(not ok, "non-convecting supercritical (dead density coupling) fails")

# ---------------------------------------------------------------------------
print("measured-flux window (V2-A17):")
ok, _ = gate_heat.check_flux_window(1.01e-6, 1.0e-6, 0.02, 0.0, "syn")
expect(ok, "1% flux deviation passes the 2% window")
ok, _ = gate_heat.check_flux_window(1.05e-6, 1.0e-6, 0.02, 0.0, "syn")
expect(not ok, "5% flux deviation fails the 2% window")
ok, _ = gate_heat.check_flux_window(-1.0e-6, 1.0e-6, 0.02, 0.0, "syn")
expect(not ok, "sign-flipped flux (direction error) fails")
ok, _ = gate_heat.check_flux_window(1.0e-9, 0.0, 0.02, 2.0e-9, "syn")
expect(ok, "zero-flux sub-case passes inside the absolute floor")
ok, _ = gate_heat.check_flux_window(1.0e-8, 0.0, 0.02, 2.0e-9, "syn")
expect(not ok, "residual flow in the zero-flux sub-case fails")

print("g7(a2) bulk-formula chain (rule-1 spot values, hand-computed):")
# q_sat(25 C, 101.325 kPa): e_sat = 0.6108 exp(17.27*25/262.3) = 3.1677 kPa
# -> q_sat = 0.622*3.1677/(101.325 - 0.376*3.1677) = 0.0196769.
expect(abs(gate_heat.q_sat(25.0, 101.325) - 0.0196769) < 2e-6,
       "q_sat(25 C) hand value 0.0196769")
# R_air(2 m/s) = 94.909 * 2^-0.9036 = 50.73 s/m.
expect(abs(gate_heat.r_air(2.0) - 50.73) < 0.01, "R_air(2 m/s) hand value 50.73")
# rho_a(25 C) = 101325/(287.05*298.15) = 1.1839 kg/m3.
expect(abs(gate_heat.air_density(25.0, 101.325) - 1.1839) < 1e-3,
       "air density hand value 1.1839")
# The g7a2 forcing: Q_net(20) = +44.9 and Q_net(22) = -122.2 W/m2 by
# hand (LW: eps(350 - sigma*293.15^4) = -66.8; latent at rho_a(T_w),
# q_sat(20) = 0.014480: 155.9; sensible -117.3; SW 150).
q_air = 0.6 * gate_heat.q_sat(25.0, 101.325)
expect(abs(gate_heat.net_heat_flux(20.0, 25.0, 101.325, 2.0, q_air, 150.0, 350.0)
           - 44.9) < 1.0, "Q_net(20 C) hand value +44.9 W/m2")
expect(abs(gate_heat.net_heat_flux(22.0, 25.0, 101.325, 2.0, q_air, 150.0, 350.0)
           + 122.2) < 1.0, "Q_net(22 C) hand value -122.2 W/m2")
root = gate_heat.equilibrium_root(25.0, 101.325, 2.0, q_air, 150.0, 350.0)
expect(20.0 < root < 22.0, "equilibrium root lands in the hand-checked bracket")
expect(abs(gate_heat.net_heat_flux(root, 25.0, 101.325, 2.0, q_air, 150.0, 350.0))
       < 1e-6, "Q_net vanishes at the root (residual substitution)")
# Sensitivity: a 10% exchange error moves the settled temperature well
# past the 0.05 C consistency bound.
ok, _ = gate_heat.check_equilibrium_consistency(root + 0.2, root, 0.05)
expect(not ok, "0.2 C settled-temperature offset fails the consistency check")

print("g8 analysis machinery:")
nx, nz = 40, 20
dx = dz = 0.05
x_c = (np.arange(nx) + 0.5) * dx
depths = (np.arange(nz) + 0.5) * dz
z_pin, h_eff, width = 0.025, 0.95, 2.0
zeta = np.clip((depths - z_pin) / h_eff, 0.0, 1.0)
mode = np.cos(math.pi * x_c / (width / 2.0))[:, None] * np.sin(math.pi * zeta)[None, :]
conduction = 20.0 + 10.0 * zeta
synth = conduction[None, :] + 0.37 * mode
amp = gate_heat.hrl_mode_amplitude(synth, x_c, depths, z_pin, h_eff, width)
expect(abs(amp - 0.37) < 1e-9, "mode projector recovers a seeded 0.37 K amplitude")
# The ORIGINAL sine seed sin(pi x / H) is orthogonal to the admissible
# cosine mode — the projection is ~0 (the V2-A17 seed correction).
sine = conduction[None, :] + 0.5 * (np.sin(math.pi * x_c / (width / 2.0))[:, None] *
                                    np.sin(math.pi * zeta)[None, :])
amp = gate_heat.hrl_mode_amplitude(sine, x_c, depths, z_pin, h_eff, width)
expect(abs(amp) < 1e-9, "the rejected sine seed projects to zero (orthogonality)")
# Nusselt: pure conduction gives exactly 1 (q = 0, linear profile).
lam, rc_w = 2.092, 4.184e6
nu = gate_heat.hrl_nusselt(np.tile(conduction, (nx, 1)), np.zeros(nx),
                           nz // 2 - 1, dz, lam, rc_w, 10.0, h_eff)
expect(abs(nu - 1.0) < 1e-9, "conduction-only Nusselt = 1 exactly")
# A synthetic roll (q up where warm) must raise Nu above 1.
q_roll = 1.0e-6 * np.cos(math.pi * x_c / (width / 2.0))
warm = np.tile(conduction, (nx, 1)) + 0.5 * mode
nu = gate_heat.hrl_nusselt(warm, q_roll, nz // 2 - 1, dz, lam, rc_w, 10.0, h_eff)
expect(nu > 1.0, "correlated upwelling raises the Nusselt number")

print("heat orientation battery (array/coordinate consistency):")
length = 4.0
n = 8
d = length / n
centers = (np.arange(n) + 0.5) * d
base = (1.7 * centers[None, :] + 0.3 * centers[:, None] ** 2 +
        np.outer(np.sin(centers), np.cos(centers)))  # injective test field
for name, (array_op, coord_op) in gate_heat.dihedral_table(length).items():
    variant = array_op(base)
    worst = 0.0
    for j in range(n):
        for i in range(n):
            xv, yv = coord_op(centers[i], centers[j])
            iv = int(round(xv / d - 0.5))
            jv = int(round(yv / d - 0.5))
            worst = max(worst, abs(variant[jv, iv] - base[j, i]))
    expect(worst < 1e-12,
           f"{name}: array op and coordinate op agree on every cell")
fields = {name: np.stack([op(base)] * 3, axis=2)
          for name, (op, _) in gate_heat.dihedral_table(length).items()}
ok, _ = gate_heat.check_heat_orientations(fields, length, 1e-12, "syn")
expect(ok, "consistent dihedral fields pass")
fields["rot90"] = np.stack([gate_heat.dihedral_table(length)["rot270"][0](base)] * 3,
                           axis=2)
ok, _ = gate_heat.check_heat_orientations(fields, length, 1e-12, "syn")
expect(not ok, "a swapped rotation direction fails")

# ---------------------------------------------------------------------------
if failures:
    print(f"test_g678_gates: {len(failures)} FAILED expectation(s)")
    for name in failures:
        print(f"  - {name}")
    sys.exit(1)
print("test_g678_gates: OK")
