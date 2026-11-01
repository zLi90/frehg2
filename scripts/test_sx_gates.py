#!/usr/bin/env python3
"""test_sx_gates.py — the §6.3 self-check battery for the Stage-3 T10
solve-consistent transport gates (sx0–sx5; plan §12B). Each checker
accepts exact synthetic input and fails on a corruption modelled on the
defect it exists to catch: a transport advecting with post-limiter flow
rates instead of the solve's face volumes (the anchor and the drying-face
salt loss), a uniform field concentrated above its sources, an extremum
minted at a front."""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_sx as sx  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


print("telescoping (sx0 / C1)")
rng = np.random.default_rng(10)
ny, nx, area = 4, 6, 2.0
qx = rng.normal(0.0, 1.0e-3, (ny, nx + 1))
qy = rng.normal(0.0, 1.0e-3, (ny + 1, nx))
d_eta = -((qx[:, 1:] - qx[:, :-1]) + (qy[1:, :] - qy[:-1, :])) / area
res = sx.telescope_residuals(d_eta, area, qx, qy)
expect(float(np.max(np.abs(res))) < 1.0e-18, "exact solve volumes telescope")
bad = qx.copy()
bad[2, 3] = 0.0  # a face the limiter zeroed AFTER the water moved
res = sx.telescope_residuals(d_eta, area, bad, qy)
expect(float(np.max(np.abs(res))) > 1.0e-5,
       "negative: a post-hoc zeroed face breaks the balance")

print("uniform (sx1 / J1)")
frames = [np.full(10, 35.0), np.full(10, 35.0)]
wet = [np.ones(10, bool), np.array([True] * 6 + [False] * 4)]
ok, _ = sx.check_uniform(frames, wet, 35.0, 1.0e-9, "u")
expect(ok, "a uniform field passes (dry cells ignored)")
bad = [frames[0], frames[1].copy()]
bad[1][2] = 35.035  # the consistent update's film concentration
ok, _ = sx.check_uniform(bad, wet, 35.0, 1.0e-9, "u")
expect(not ok, "negative: a 35.035 concentrated film fails")

print("anchor (sx1/sx3 / J2)")
ok, _ = sx.check_anchor(1.0e-9, 1.0e4, 1.0e-10, "a")
expect(ok, "a solve-residual-class anchor passes")
ok, _ = sx.check_anchor(28.4, 1.19e3, 1.0e-10, "a")
expect(not ok, "negative: the ws1n +28.4 psu m^3 anchor (2.4 % of M) fails")

print("extrema (sx4)")
frames = [np.array([35.0, 30.0, 5.0, 5.0]), np.array([34.2, 20.0, 8.0, 5.0])]
wet = [np.ones(4, bool)] * 2
ok, _ = sx.check_extrema(frames, wet, 5.0, 35.0, 1.0e-9)
expect(ok, "a front inside its range passes")
bad = [frames[0], np.array([35.4, 20.0, 8.0, 5.0])]
ok, _ = sx.check_extrema(bad, wet, 5.0, 35.0, 1.0e-9)
expect(not ok, "negative: a minted maximum fails")
bad = [frames[0], np.array([34.0, 20.0, 8.0, 4.3])]
ok, _ = sx.check_extrema(bad, wet, 5.0, 35.0, 1.0e-9)
expect(not ok, "negative: a minted minimum fails")

print("total (sx2/sx4)")
ok, _ = sx.check_total(np.array([0.525, 0.525 * (1 + 1e-14), 0.525]), 1.0e-12, "t")
expect(ok, "an exact total passes")
ok, _ = sx.check_total(np.array([0.525, 0.525 - 0.109, 0.416]), 1.0e-12, "t")
expect(not ok, "negative: the et8d -0.109 psu m^3 stranding fails")

if failures:
    print(f"{len(failures)} sx gate-logic check(s) FAILED")
    sys.exit(1)
print("all sx gate-logic checks passed")
