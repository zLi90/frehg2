#!/usr/bin/env python3
"""test_pd_gates.py — the §6.3 self-check battery for the Stage-3 T3
periodic-topology gate criteria (pd1–pd8; Stage-3 plan §7).

1. Reference sanity (§6.3 rule 1) — the pd1 discrete balance is checked
   two independent ways: the bisection solution satisfies the fixed-point
   relation to 1e-12, and it approaches the independent continuum closed
   form as dt -> 0 (the double drag factor's O(dt) bias vanishes).
2. Sensitivity (§6.3 rule 2) — every check fails on a corrupted input
   modeled on the physical failure it exists to catch (a seam glitch, a
   booked seam flux, an x-flux leak through the wrapped matrix legs).

Also runs the shared Stage-3 harness battery (plan §3.1).
Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))
sys.path.insert(0, str(ROOT / "scripts"))

import gate_periodic as pd  # noqa: E402
import test_s3_shared  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


print("pd1 discrete balance: two independent evaluations")
u10, h, n, dt = 10.0, 0.5, 0.05, 2.0
u_disc = pd.discrete_balance(u10, h, n, dt)
# Evaluation 1: the fixed-point relation holds at the bisection solution.
cd_w = pd.cd_garratt(u10)
cd_b = pd.manning_cd(n, h)
tau = pd.AIR_DENSITY * cd_w * (u10 - u_disc) ** 2
d = 1.0 / (1.0 + 0.5 * dt * cd_b * u_disc / h)
resid = (u_disc + dt * tau / (pd.WATER_DENSITY * h)) * d * d - u_disc
expect(abs(resid) < 1.0e-12, f"fixed-point residual vanishes ({resid:.2e})")
# Evaluation 2: dt -> 0 recovers the continuum closed form.
u_cont = pd.continuum_balance(u10, h, n)
u_small = pd.discrete_balance(u10, h, n, 1.0e-4)
print(f"    u_disc(dt=2) = {u_disc:.5f}, u_disc(dt->0) = {u_small:.5f}, "
      f"continuum = {u_cont:.5f} m/s")
expect(abs(u_small - u_cont) / u_cont < 1.0e-4,
       "dt -> 0 limit recovers the continuum closed form")
expect(abs(pd.cd_garratt(10.0) - 1.42e-3) < 1.0e-15, "garratt(10) = 1.42e-3 (hand)")
# The plan's §7 hand prediction at these parameters (assessment P8 class).
expect(0.06 < u_disc < 0.09, f"balance is the 0.07 m/s class ({u_disc:.4f})")

print("sensitivity: each check fails on its target corruption")
flat = np.zeros((4, 8))
expect(pd.check_flat(flat, 1.0e-12)[0], "flatness accepts a flat field")
tilt = flat.copy()
tilt[:, 4] += 1.0e-9
expect(not pd.check_flat(tilt, 1.0e-12)[0], "flatness: a 1e-9 m seam step fails")

u_field = np.full((4, 8), u_disc)
expect(pd.check_balance(u_field, u_disc, 1.0e-10)[0], "balance accepts itself")
expect(not pd.check_balance(u_field * 1.8, u_disc, 1.0e-10)[0],
       "balance: the walled-channel speedup (1.8x) fails")

expect(pd.check_steady(np.full(5, u_disc), 1.0e-13)[0], "steadiness accepts frozen")
expect(not pd.check_steady(np.array([u_disc, u_disc * (1 + 1e-9)]), 1.0e-13)[0],
       "steadiness: a 1e-9 relative drift fails")

rng = np.random.default_rng(3)
base = rng.normal(size=(6, 12))
rolled_output = np.roll(base, 3, axis=1)       # a correct periodic run
unrolled = np.roll(rolled_output, -3, axis=1)  # the gate's comparison
expect(pd.check_roll(base, unrolled, 3, 1.0e-11, "eta")[0],
       "roll identity accepts a translation-invariant pair")
glitch = unrolled.copy()
glitch[2, 0] += 1.0e-9                          # a seam-column defect
expect(not pd.check_roll(base, glitch, 3, 1.0e-11, "eta")[0],
       "roll identity: a 1e-9 seam glitch fails the gw bound")
# S3-A6 surface-lane bounds: the branch-amplified rounding class passes,
# a defect-class seam error still fails by orders of magnitude.
amplified = unrolled + rng.normal(scale=2.0e-8, size=unrolled.shape)
expect(pd.check_roll(base, amplified, 3, 1.0e-6, "uu")[0],
       "surface roll bound accepts the measured amplification class (7e-8)")
seam_defect = unrolled.copy()
seam_defect[:, 0] += 1.0e-4                      # defect-class seam column
expect(not pd.check_roll(base, seam_defect, 3, 1.0e-6, "uu")[0],
       "surface roll bound: a 1e-4 seam-column defect fails")
tracer_defect = unrolled.copy()
tracer_defect[:, 0] += 1.0e-2                    # defect-class tracer error
expect(not pd.check_roll(base, tracer_defect, 3, 1.0e-4, "tracer")[0],
       "tracer roll bound: a 1e-2 seam-column defect fails")

# S3-A6 restart checks: flow fields absolute, tracer relative.
full_f = rng.normal(scale=0.05, size=(6, 12))
expect(pd.check_restart_field(full_f, full_f + 1.9e-15, "uu",
                              tol_abs=1.0e-12)[0],
       "restart flow accepts the measured ulp class (1.9e-15)")
expect(not pd.check_restart_field(full_f, full_f + 1.0e-9, "uu",
                                  tol_abs=1.0e-12)[0],
       "restart flow: a 1e-9 m/s restart-order defect fails")
full_c = np.abs(rng.normal(scale=3.0, size=(6, 12))) + 1.0
expect(pd.check_restart_field(full_c, full_c * (1.0 + 4.0e-6), "tracer",
                              tol_rel=1.0e-4)[0],
       "restart tracer accepts the measured branch class (4e-6 rel)")
expect(not pd.check_restart_field(full_c, full_c * (1.0 + 1.0e-2), "tracer",
                                  tol_rel=1.0e-4)[0],
       "restart tracer: a 1e-2 rel defect fails")

frames = [np.zeros((6, 12)), np.zeros((6, 12))]
frames[1][:, 0] = 2.0                           # tracer arrived past the seam
expect(pd.check_seam_crossing(frames, 0, 1.0)[0],
       "seam crossing accepts an advected slab")
expect(not pd.check_seam_crossing([frames[0], frames[0]], 0, 1.0)[0],
       "seam crossing: a run whose tracer never crosses fails (anti-vacuity)")

c0 = np.full((8, 12, 6), 5.0)
c1 = c0.copy()
c1[7, :, 5] = 5.3                                # seepage brought salt in
expect(pd.check_coupled_exchange(c0, c1, 0.05)[0],
       "coupled exchange accepts the measured tidal signal (0.3 psu)")
expect(not pd.check_coupled_exchange(c0, c0, 0.05)[0],
       "coupled exchange: a hollow lane (no subsurface change) fails")

audit = np.zeros((25, 7))
audit[:, 0] = np.arange(25) * 300.0
audit[:, 4] = 10.0 * np.sin(2.0 * np.pi * audit[:, 0] / 3600.0)
audit[:, 1] = 4000.0 - audit[:, 4]
expect(pd.check_surface_ledger(audit, 1.0e-11)[0], "ledger accepts consistency")
seam_booked = audit.copy()
seam_booked[:, 4] += np.linspace(0.0, 1.0e-3, 25)  # a seam booked as boundary
expect(not pd.check_surface_ledger(seam_booked, 1.0e-11)[0],
       "ledger: a seam booked as boundary flux fails")

head = np.zeros((10, 6, 5))
for j in range(10):
    head[j] = 3.0 - 0.02 * j
expect(pd.check_linear_head(head, np.arange(10) + 0.5, 3.0, 2.8, 1.0e-10)[0],
       "linear head accepts the Darcy profile")
kink = head.copy()
kink[5] += 1.0e-8
expect(not pd.check_linear_head(kink, np.arange(10) + 0.5, 3.0, 2.8, 1.0e-10)[0],
       "linear head: a 1e-8 m kink fails")
xdep = head.copy()
xdep[:, 0, :] += 1.0e-9
expect(not pd.check_linear_head(xdep, np.arange(10) + 0.5, 3.0, 2.8, 1.0e-10)[0],
       "linear head: a 1e-9 m x dependence fails")

qx = np.zeros((10, 6, 5))
expect(pd.check_zero_x_flux(qx, 1.0e-7, 1.0e-14)[0], "zero x flux accepts zeros")
leak = qx.copy()
leak[3, 0, 2] = 1.0e-12                          # a wrapped-leg leak
expect(not pd.check_zero_x_flux(leak, 1.0e-7, 1.0e-14)[0],
       "zero x flux: a 1e-12 m/s wrapped-leg leak fails")

print("shared Stage-3 harness battery (plan §3.1):")
failures.extend(test_s3_shared.run_shared_checks())

if failures:
    print(f"{len(failures)} pd gate-logic check(s) FAILED")
    sys.exit(1)
print("all pd gate-logic checks passed")
