#!/usr/bin/env python3
"""test_dg_gates.py — the §6.3 self-check battery for the Stage-3 T4
diagnostics gate criteria (dg1/dg2; Stage-3 plan §8).

1. Reference sanity (§6.3 rule 1) — the dg1 rate-ledger identity is
   verified on a synthetic bookkeeping pair constructed from the SAME
   per-step gain numbers both the rate frames and the audit column would
   sum, so the identity holds exactly by construction.
2. Sensitivity (§6.3 rule 2) — every check fails on a corrupted input
   modeled on the physical failure it exists to catch: a rate frame that
   dropped one substep's exchange (a window-lag bug), a snapshot that
   actually carries the exchange (a vacuous V3 lane), a decimated row with
   one perturbed column (an accumulate/write split that changed the sums),
   a missing t_end row, and a schedule that lost the checkpoint-forced row.

Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_diagnostics as dg  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


print("dg1 rate-ledger identity: synthetic bookkeeping pair")
rng = np.random.default_rng(7)
area = 4.0
ny, nx = 4, 4
times = [0.0, 300.0, 600.0]
# Per-step gains over two 300 s intervals (300 steps each), both signs.
gains = [rng.normal(scale=1.0e-6, size=(300, ny, nx)),
         rng.normal(scale=1.0e-6, size=(300, ny, nx))]
rates = [np.zeros((ny, nx))] + [g.sum(axis=0) / 300.0 for g in gains]
audit_times = np.arange(0.0, 601.0)
seep = np.zeros(601)
for k, g in enumerate(gains):
    per_step = g.sum(axis=(1, 2)) * area
    seep[1 + 300 * k:301 + 300 * k] = seep[300 * k] + np.cumsum(per_step)
ok, _ = dg.check_rate_ledger(times, rates, audit_times, seep, area, 1.0e-12)
expect(ok, "identity accepts the exact bookkeeping pair")
dropped = [rates[0], rates[1].copy(), rates[2]]
dropped[1] -= gains[0][17] / 300.0          # one substep's exchange lost
ok, _ = dg.check_rate_ledger(times, dropped, audit_times, seep, area, 1.0e-12)
expect(not ok, "identity: a dropped-substep rate frame fails")

print("dg1 anti-vacuity (V3 reproduction)")
depth = [np.full((ny, nx), 0.02)] * 3
snap_zero = [np.zeros((ny, nx))] * 3
rate_live = [np.full((ny, nx), 1.0e-7)] * 3
ok, _ = dg.check_rate_antivacuity(snap_zero, rate_live, depth, 1.0e-5)
expect(ok, "anti-vacuity accepts zero snapshots with live rates")
snap_live = [np.full((ny, nx), 1.0e-7)] * 3
ok, _ = dg.check_rate_antivacuity(snap_live, rate_live, depth, 1.0e-5)
expect(not ok, "anti-vacuity: a snapshot that carries the exchange fails"
               " (vacuous lane)")

print("dg2 decimation: bitwise row subset + forced rows")
full = np.column_stack([np.arange(0.0, 11.0),
                        rng.normal(size=11), np.cumsum(rng.normal(size=11))])
dec = full[::5].copy()                       # rows at 0, 5, 10
ok, _ = dg.check_decimated_subset(full, dec, 10.0, "mass_audit")
expect(ok, "subset accepts bitwise decimated rows ending at t_end")
shifted = dec.copy()
shifted[1, 2] += 1.0e-15                     # one column off by 1 ulp class
ok, _ = dg.check_decimated_subset(full, shifted, 10.0, "mass_audit")
expect(not ok, "subset: a 1e-15 column perturbation fails (bitwise)")
ok, _ = dg.check_decimated_subset(full, dec[:-1], 10.0, "mass_audit")
expect(not ok, "subset: a missing t_end row fails")

expected = [0.0, 77.0, 154.0, 180.0, 231.0]
ok, _ = dg.check_row_times(np.array(expected), expected, "sched")
expect(ok, "schedule accepts the exact forced-row set")
ok, _ = dg.check_row_times(np.array([0.0, 77.0, 154.0, 231.0]), expected,
                           "sched")
expect(not ok, "schedule: a lost checkpoint-forced row fails")
ok, _ = dg.check_row_times(np.array(expected + [233.0]), expected, "sched")
expect(not ok, "schedule: an extra unscheduled row fails")

if failures:
    print(f"{len(failures)} dg gate-logic check(s) FAILED")
    sys.exit(1)
print("all dg gate-logic checks passed")
