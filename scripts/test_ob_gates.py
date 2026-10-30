#!/usr/bin/env python3
"""test_ob_gates.py — the §6.3 self-check battery for the Stage-3 T2
open-boundary gate criteria (ob1/ob3/ob5; Stage-3 plan §6).

1. Reference sanity (§6.3 rule 1) — the ob1 discrete chain is checked two
   independent ways: its per-face steady residuals vanish at its own
   solution, and its external-lane velocity agrees with the continuum
   Bernoulli + Manning-friction balance (an independent closed form) to
   the discretisation class.
2. Sensitivity (§6.3 rule 2) — every check fails on a corrupted input
   modeled on the physical failure it exists to catch, including the
   factor-two weight-1 entrance (plan §3.4 of the assessment).

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

import gate_openbc as ob  # noqa: E402
import test_s3_shared  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


# The ob1 configuration (benchmarks/ob-openbc/ob1-channel.yaml).
OB1 = dict(nx=24, dx=1.0, dy=1.0, dt=0.25, bed=0.0, n_manning=0.01,
           eta_west=1.005)
ETA_EAST = 1.0

print("ob1 chain, lane A (external, still reservoir):")
chain_a = ob.Chain(closure="external", u_ext=0.0, **OB1)
q_a, eta_a, u_a = ob.chain_solve(chain_a, ETA_EAST)
res = np.max(np.abs(chain_a.residuals(q_a)[1:]))
expect(res < 1.0e-12, f"per-face steady residuals vanish (max {res:.2e})")
expect(np.all(np.diff(eta_a[1:]) <= 1.0e-15),
       "eta decreases monotonically down the channel")
# Independent evaluation 2: continuum Bernoulli + Manning friction.
# u^2/2g (entrance) + sum Sf dx (friction) = total drop; solve by fixed point.
u_bern = 0.3
for _ in range(100):
    h = 1.0
    sf = OB1["n_manning"]**2 * u_bern**2 / h**(4.0 / 3.0)
    u_bern = np.sqrt(max(2.0 * ob.GRAVITY * (0.005 - sf * (OB1["nx"] - 2)), 1e-12))
u_in_a = u_a[1]
dev = abs(u_in_a - u_bern) / u_bern
print(f"    chain u_in = {u_in_a:.4f} m/s; continuum Bernoulli+Manning = "
      f"{u_bern:.4f} m/s; rel dev {dev:.3f}")
expect(dev < 0.15,
       "chain agrees with the continuum closed form to the discretisation class")
head_a = u_in_a**2 / (2.0 * ob.GRAVITY)
ratio_a = (1.005 - eta_a[2]) / head_a
print(f"    chain entrance ratio R_A = {ratio_a:.4f} (the seed measured 1.01)")
expect(0.80 < ratio_a < 1.20, "entrance ratio is the Bernoulli head class")

print("ob1 chain, lane B (external, u_ext = 0.1):")
chain_b = ob.Chain(closure="external", u_ext=0.1, **OB1)
q_b, eta_b, u_b = ob.chain_solve(chain_b, ETA_EAST)
expect(q_b > q_a, "incoming momentum buys more discharge than the still reservoir")
expect(np.max(np.abs(chain_b.residuals(q_b)[1:])) < 1.0e-12,
       "lane B residuals vanish")

print("ob1 chain, lane C (zero_gradient):")
chain_c = ob.Chain(closure="zero_gradient", **OB1)
q_c, eta_c, u_c = ob.chain_solve(chain_c, ETA_EAST)
u_in_c = u_c[1]
ratio_c = (1.005 - eta_c[2]) / (u_in_c**2 / (2.0 * ob.GRAVITY))
print(f"    chain u_in: A {u_in_a:.4f}, B {u_b[1]:.4f}, C {u_in_c:.4f} m/s; "
      f"R_C = {ratio_c:.4f}")
expect(u_in_c >= 2.0 * u_in_a, "zero_gradient runs >= 2x the still-reservoir lane")
expect(ratio_c <= 0.2, "zero_gradient pays almost no entrance head")
cfl = u_in_c * OB1["dt"] / OB1["dx"]
expect(cfl < 0.5, f"lane C stays under the advection-damping ramp (CFL {cfl:.3f})")

print("sensitivity: each check fails on its target corruption")
ok, _ = ob.check_velocity_profile(u_a[1:], u_a[1:], 1.0e-6)
expect(ok, "velocity profile accepts the chain itself")
w1 = np.sqrt(ob.GRAVITY * 0.005)  # the weight-1 entrance: u = sqrt(g dh)
bad = u_a[1:].copy()
bad[:] = w1
expect(not ob.check_velocity_profile(bad, u_a[1:], 1.0e-6)[0],
       "the factor-two weight-1 entrance fails the velocity check")
expect(ob.check_steadiness(np.full(5, 0.31), 1.0e-9)[0],
       "steadiness accepts a frozen tail")
expect(not ob.check_steadiness(np.array([0.31, 0.31, 0.311]), 1.0e-9)[0],
       "steadiness: a 1e-3 m/s trailing drift fails")
expect(ob.check_entrance_ratio(1.005, float(eta_a[2]), u_in_a, ratio_a, 0.05)[0],
       "entrance ratio accepts the chain state")
expect(not ob.check_entrance_ratio(1.005, float(eta_c[2]), u_in_c, ratio_a, 0.05)[0],
       "entrance ratio: the zero_gradient state fails lane A's check")
expect(ob.check_lane_contrast(u_in_a, u_in_c, 2.0, ratio_a, ratio_c, 0.2)[0],
       "lane contrast accepts the chain pair")
expect(not ob.check_lane_contrast(u_in_a, 1.5 * u_in_a, 2.0, ratio_a, ratio_c, 0.2)[0],
       "lane contrast: a 1.5x-only speedup fails")

print("ob3 ledger and prism checks:")
audit = np.zeros((25, 7))
audit[:, 0] = np.arange(25) * 300.0
audit[:, 4] = 10.0 * np.sin(2.0 * np.pi * audit[:, 0] / 3600.0)  # outflow
audit[:, 1] = 4000.0 - audit[:, 4]  # volume consistent with exchange
expect(ob.check_surface_ledger(audit, 1.0e-12)[0],
       "ledger accepts a consistent table")
hole = audit.copy()
hole[-1, 1] += 1.0e-6
expect(not ob.check_surface_ledger(hole, 1.0e-12)[0],
       "ledger: a 1e-6 m^3 hole fails")
prisms = ob.prism_per_cycle(audit, 3600.0)
expect(len(prisms) == 2 and abs(prisms[0] - prisms[1]) < 1.0e-9,
       "prism extraction: two equal cycles from the sinusoid")
expect(ob.check_prism_agreement([10.0], [10.001], 1.0e-3, 2.0, 1.0e-4)[0],
       "prism agreement accepts an Fr^2-class difference")
expect(not ob.check_prism_agreement([10.0], [10.5], 1.0e-3, 2.0, 1.0e-4)[0],
       "prism agreement: a 5% exchange shift fails")

print("ob5 KE checks:")
ke = np.array([100.0, 80.0, 60.0, 50.0, 45.0, 44.0])
expect(ob.check_ke_decay(ke, 1, 1.0e-6)[0], "KE decay accepts a decaying series")
grow = ke.copy()
grow[4] = 70.0
expect(not ob.check_ke_decay(grow, 1, 1.0e-6)[0], "KE decay: a mid-run rise fails")
expect(ob.check_ke_influx(1.0e-4, 100.0, 1.0e-3)[0], "influx bound accepts quiet faces")
expect(not ob.check_ke_influx(1.0, 100.0, 1.0e-3)[0], "influx bound: 1% of KE0 fails")
expect(ob.check_influx_contrast(5.0, 0.1, 10.0, 1.0)[0],
       "influx contrast accepts a 50x import")
expect(not ob.check_influx_contrast(0.5, 0.1, 10.0, 1.0)[0],
       "influx contrast: a 5x-only import fails (anti-vacuity)")
vv_row = np.array([-0.2, -0.1, 0.1, 0.0])
influx = ob.stage_ke_influx(vv_row, np.full(4, 2.0), 1.0, 30.0)
hand = 0.5 * 998.0 * (0.2**3 + 0.1**3) * 2.0 * 1.0 * 30.0
expect(abs(influx - hand) < 1.0e-9, "stage KE influx matches the hand value")

print("shared Stage-3 harness battery (plan §3.1):")
failures.extend(test_s3_shared.run_shared_checks())

if failures:
    print(f"{len(failures)} ob gate-logic check(s) FAILED")
    sys.exit(1)
print("all ob gate-logic checks passed")
