#!/usr/bin/env python3
"""test_ws_gates.py — the §6.3 self-check battery for the Stage-3 T9 paired
water–salt store gates (ws0–ws5; plan §12A).

1. Reference sanity: the closed forms (the R5 release mix, the S3-A12
   et8b rewet closed form) against direct bookkeeping.
2. Sensitivity: every checker accepts exact synthetic books and fails on a
   corruption modelled on the defect it exists to catch — a release that
   divides by vBasis alone (T8's), an identity missing film_removed (the
   S3-A12 1.2–3.2 % ledger gap), a store holding salt without its water
   (the S3-A12 mechanism), a store debited by evaporation (R2).

Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_ws as ws  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


print("closed forms")
# R5: 2 mm film at 20 psu meets a store of 0.5 mm carrying 35-psu water.
area = 0.5
v, m = 2.0e-3 * area, 20.0 * 2.0e-3 * area
vs, ms = 0.5e-3 * area, 35.0 * 0.5e-3 * area
c = ws.mix_release(m, v, ms, vs)
expect(abs(c - (20.0 * 2.0 + 35.0 * 0.5) / 2.5) < 1.0e-12,
       f"R5 mix = depth-weighted mean ({c:.6f} psu)")
expect(min(20.0, 35.0) <= c <= max(20.0, 35.0),
       "R5 mix is a convex combination of the sources")
t8 = (m + ms) / v  # T8: divides by the flux basis alone
expect(abs(t8 - c) / c > 1.0e-3, "negative: a vBasis-only release differs")
# et8b rewet: 5 mm at 35 psu dried with a parked 4e-6 m remainder, 2 mm rain.
s = ws.pond_rewet_concentration(35.0, 5.0e-3, 1.0, 4.0e-6, 2.0e-3)
expect(abs(s - 35.0 * 5.0e-3 / (4.0e-6 + 2.0e-3)) < 1.0e-9,
       f"et8b rewet closed form ({s:.4f} psu)")
expect(abs(ws.pond_rewet_concentration(35.0, 5.0e-3, 1.0, 0.0, 2.0e-3)
           - 87.5) < 1.0e-9, "V_r = 0 recovers the T5 closed form (87.5 psu)")

print("water identity (ws2)")
rng = np.random.default_rng(7)
n = 50


def synthetic_books(coupled: bool, store: bool):
    cols = ws.audit_columns(coupled, store)
    rows = np.zeros((n, len(cols)))
    rows[:, cols["time"]] = np.arange(n) * 60.0
    inc = lambda: np.concatenate([[0.0], np.cumsum(rng.uniform(0, 1e-2, n - 1))])
    for name in ("rain", "evaporation", "boundary_outflow", "bc_inflow",
                 "clamped"):
        rows[:, cols[name]] = inc()
    if coupled:
        rows[:, cols["seepage"]] = inc() - inc()
        rows[:, cols["stage_deposit"]] = inc()
    rows[:, cols["film_removed"]] = inc()
    pending = rng.uniform(0, 1e-4, n) * (0.0 if store else 1.0)
    rows[:, cols["film_pending"]] = pending
    flow = (rows[:, cols["rain"]] - rows[:, cols["evaporation"]]
            - rows[:, cols["boundary_outflow"]] + rows[:, cols["bc_inflow"]]
            + rows[:, cols["clamped"]])
    if coupled:
        flow += rows[:, cols["seepage"]] - rows[:, cols["stage_deposit"]]
    v0 = 10.0
    if store:
        release = rows[:, cols["film_removed"]] * rng.uniform(0.2, 0.8)
        release = np.maximum.accumulate(release)
        store_v = rows[:, cols["film_removed"]] - release
        rows[:, cols["store_release"]] = release
        rows[:, cols["store_volume"]] = store_v
        rows[:, cols["volume"]] = v0 + flow - store_v
    else:
        rows[:, cols["volume"]] = (v0 + flow - rows[:, cols["film_removed"]]
                                   - pending + pending[0])
    return rows, cols


for coupled in (False, True):
    for store in (False, True):
        rows, cols = synthetic_books(coupled, store)
        tag = f"{'coupled' if coupled else 'surface'}/{'store' if store else 'discard'}"
        ok, msg = ws.check_water_identity(rows, cols, store, 1.0e-10)
        expect(ok, f"{tag}: exact books close ({msg})")
        if not store:
            bad = rows.copy()
            bad[:, cols["film_removed"]] = 0.0  # the S3-A12 unbooked removal
            ok, _ = ws.check_water_identity(bad, cols, store, 1.0e-10)
            expect(not ok, f"{tag}: negative — dropping film_removed fails")
            bad = rows.copy()
            bad[:, cols["film_pending"]] = 0.0  # the R8 pending films
            ok, _ = ws.check_water_identity(bad, cols, store, 1.0e-10)
            expect(not ok, f"{tag}: negative — dropping film_pending fails")
        else:
            bad = rows.copy()
            bad[:, cols["store_volume"]] = 0.0  # parked water deleted
            ok, _ = ws.check_water_identity(bad, cols, store, 1.0e-10)
            expect(not ok, f"{tag}: negative — a deleted store fails")
            bad = rows.copy()
            bad[n // 2:, cols["volume"]] += 1.0e-7  # one-off created film
            ok, _ = ws.check_water_identity(bad, cols, store, 1.0e-10)
            expect(not ok, f"{tag}: negative — a 1e-7 m^3 creation fails")

print("salt identity (ws2)")
audit = np.zeros((n, 12))
audit[:, 0] = np.arange(n) * 60.0
terms = np.cumsum(rng.uniform(-1e-3, 1e-3, (n, 7)), axis=0)
terms[0] = 0.0
audit[:, 4:11] = terms
audit[:, 11] = rng.uniform(0, 0.1, n)
audit[:, 2] = 100.0 + rng.uniform(-1, 1, n)
audit[:, 1] = 50.0 + terms.sum(axis=1) - (audit[:, 2] - audit[0, 2]) \
    - (audit[:, 11] - audit[0, 11])
ok, msg = ws.check_salt_identity(audit, True, 1.0e-10)
expect(ok, f"exact salt books close ({msg})")
bad = audit.copy()
bad[:, 11] *= 1.0e3  # the S3-A12 store growth with no source column
ok, _ = ws.check_salt_identity(bad, True, 1.0e-10)
expect(not ok, "negative: a store growing without a booked source fails")

print("bounds (ws1 / I2)")
ok, _ = ws.check_source_bound([34.9, 35.0, 35.0 * (1 + 5e-13)], 35.0, 1.0e-12,
                              "surface")
expect(ok, "source bound accepts c <= c_max (1 + 1e-12)")
ok, _ = ws.check_source_bound([35.0, 72850.0], 35.0, 1.0e-12, "surface")
expect(not ok, "negative: the S3-A12 72 850 psu film fails")
ok, _ = ws.check_source_bound([35.0 * (1 + 1e-10)], 35.0, 1.0e-12, "surface")
expect(not ok, "negative: a 1e-10 overshoot fails")
vs = np.array([1e-5, 2e-5, 0.0, 5e-6])
ms = np.array([3.5e-4, 7.0e-4, 0.0, 1.0e-4])
ok, _ = ws.check_store_bound(ms, vs, 35.0, 1.0e-12)
expect(ok, "store bound accepts M_s <= 35 V_s")
bad = ms.copy()
bad[2] = 1.76e-2  # salt without water (the zoom-A cell)
ok, _ = ws.check_store_bound(bad, vs, 35.0, 1.0e-12)
expect(not ok, "negative: salt on a cell with V_s = 0 fails")
bad = ms.copy()
bad[0] *= 57.0  # the 57x store growth
ok, _ = ws.check_store_bound(bad, vs, 35.0, 1.0e-12)
expect(not ok, "negative: a 57x concentrated store fails")

print("pairing (ws4)")
vs = np.array([1e-5, 0.0, 3e-6, 2.5e-5])
ok, _ = ws.check_paired(35.0 * vs, vs, 35.0, 1.0e-12)
expect(ok, "pairing accepts M_s = 35 V_s exactly")
bad = 35.0 * vs * (1 + 1e-9)
ok, _ = ws.check_paired(bad, vs, 35.0, 1.0e-12)
expect(not ok, "negative: a 1e-9 pairing error fails")
ok, msg = ws.check_paired(np.zeros(3), np.zeros(3), 35.0, 1.0e-12)
expect(not ok, "anti-vacuity: no parked cell fails")

print("store monotone while dry (ws3, R2/R7)")
series = np.array([[0.0], [4e-6], [4e-6], [4e-6]])
dry = np.array([[False], [True], [True], [True]])
ok, _ = ws.check_monotone_while_dry(series, dry)
expect(ok, "a store that only parks passes")
bad = np.array([[0.0], [4e-6], [3e-6], [2e-6]])  # evaporation debits it
ok, _ = ws.check_monotone_while_dry(bad, dry)
expect(not ok, "negative: a store debited by evaporation fails")

if failures:
    print(f"{len(failures)} ws gate-logic check(s) FAILED")
    sys.exit(1)
print("all ws gate-logic checks passed")
