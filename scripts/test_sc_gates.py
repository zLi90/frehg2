#!/usr/bin/env python3
"""test_sc_gates.py — the §6.3 self-check battery for the Stage-3 T12
saturated-contact classification gates (sc1–sc2; plan §12D). Each checker
accepts exact synthetic input and fails on a corruption modelled on the
defect it exists to catch: a seam-column restore chain (the T11 legacy
sweep's 2.4e-5 at n = 5), a lane that never reaches the seam, a step
without restores, and a one-ulp thread-order difference."""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_sc as sc  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


rng = np.random.default_rng(12)
base = {
    "groundwater/water_content": 0.08 + 0.32 * rng.random(2 * 55 * 101 * 25),
    "groundwater/hydraulic_head": -2.0 + 5.0 * rng.random(2 * 55 * 101 * 25),
    "surface/depth": 1.0e-3 * rng.random(2 * 55 * 101),
}
same = {key: value.copy() for key, value in base.items()}

print("sc1 rank invariance")
expect(sc.check_rank_invariance(base, same, 1.0e-12, 2)[0], "identical fields pass")
rounding = {key: value * (1.0 + 4.0e-15) for key, value in base.items()}
expect(sc.check_rank_invariance(base, rounding, 1.0e-12, 3)[0], "rounding-level drift passes")
seam = {key: value.copy() for key, value in base.items()}
# One seam-column fringe cell left at theta_s - 1.1e-6 instead of theta_s
# (the T11 legacy n = 3 divergence at (j=0, i=34, k=10)).
seam["groundwater/water_content"][(0 * 101 + 34) * 25 + 10] -= 1.1e-6
ok, msg = sc.check_rank_invariance(base, seam, 1.0e-12, 3)
expect(not ok and "water_content" in msg, "a seam-cell restore miss fails and is named")
nan_base = {key: value.copy() for key, value in base.items()}
nan_base["surface/depth"][7] = np.nan
expect(sc.check_rank_invariance(nan_base, nan_base, 1.0e-12, 2)[0],
       "masked (NaN) entries are skipped")

print("sc1 anti-vacuity")
expect(sc.check_legacy_exercises_seam({2: 1.0e-14, 3: 2.76e-6, 5: 2.37e-5}, 1.0e-8)[0],
       "legacy lane with the seam defect passes the floor")
expect(not sc.check_legacy_exercises_seam({2: 1.0e-14, 3: 6.0e-15, 4: 5.7e-14}, 1.0e-8)[0],
       "a legacy lane that never reaches the seam fails")
expect(not sc.check_legacy_exercises_seam({}, 1.0e-8)[0], "an empty legacy lane fails")
expect(sc.check_restore_active(1.91e-2, 1.0e-3)[0], "measured restore volume passes")
expect(not sc.check_restore_active(0.0, 1.0e-3)[0], "a step without restores fails")
expect(not sc.check_restore_active(float("nan"), 1.0e-3)[0], "a missing closure value fails")

print("sc2 thread bitwise")
expect(sc.check_bitwise(base, same, "threads")[0], "identical fields pass")
ulp = {key: value.copy() for key, value in base.items()}
ulp["groundwater/hydraulic_head"][123] = np.nextafter(ulp["groundwater/hydraulic_head"][123],
                                                      np.inf)
ok, msg = sc.check_bitwise(base, ulp, "threads")
expect(not ok and "hydraulic_head" in msg, "a one-ulp difference fails and is named")
short = {key: value[:-1] for key, value in base.items()}
expect(not sc.check_bitwise(base, short, "threads")[0], "a shape mismatch fails")

if failures:
    print(f"sc gate battery: FAIL ({len(failures)}: {', '.join(failures)})")
    sys.exit(1)
print("sc gate battery: PASS")
