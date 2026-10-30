#!/usr/bin/env python3
"""gate_diagnostics.py — pure checks for the Stage-3 T4 diagnostics gates
dg1/dg2/dg-restart (Stage-3 plan §8).

dg1 gates the DEFINING identity of an interval-mean rate output: the rate
frame written at output time t_k integrates back to the audit's cumulative
column over the same interval, sum(rate * A * (t_k - t_{k-1})) ==
delta(mass_audit.seepage), because both sides are sums of the very same
per-step gain numbers. dg2 gates monitor decimation as a pure row-subset
property (a decimated table's rows are bitwise rows of the every-step
table). The HDF5 plumbing lives in run_regression.py; the §6.3 battery
(scripts/test_dg_gates.py) proves every check's ability to fail.
"""

from __future__ import annotations

import numpy as np


def check_rate_ledger(times: list[float], rate_frames: list[np.ndarray],
                      audit_times: np.ndarray, audit_seepage: np.ndarray,
                      cell_area: float, tol_rel: float):
    """dg1: per output interval (t_{k-1}, t_k], the rate frame at t_k
    satisfies |sum(rate * A * dt) - delta(seepage)| <= tol * max(|delta|,
    1e-12 m^3). audit_seepage is the CUMULATIVE audit column sampled at
    audit_times (every-step rows; the output times are a subset)."""
    ok = True
    msgs = []
    worst = 0.0
    for k in range(1, len(times)):
        dt_out = times[k] - times[k - 1]
        integ = float(np.nansum(rate_frames[k])) * cell_area * dt_out
        rows_prev = np.nonzero(np.isclose(audit_times, times[k - 1],
                                          rtol=0.0, atol=1.0e-6))[0]
        rows_this = np.nonzero(np.isclose(audit_times, times[k],
                                          rtol=0.0, atol=1.0e-6))[0]
        if rows_prev.size == 0 or rows_this.size == 0:
            ok = False
            msgs.append(f"interval {k}: no audit row at the output time")
            continue
        delta = float(audit_seepage[rows_this[0]] - audit_seepage[rows_prev[0]])
        err = abs(integ - delta)
        bound = tol_rel * max(abs(delta), 1.0e-12)
        worst = max(worst, err / max(bound, 1.0e-300) * tol_rel)
        if err > bound:
            ok = False
            msgs.append(f"interval {k} ({times[k-1]:g}-{times[k]:g} s): "
                        f"integ {integ:.6e} vs delta {delta:.6e} "
                        f"(err {err:.3e} > {bound:.3e})")
    msgs.insert(0, f"rate-ledger identity over {len(times) - 1} intervals: "
                   f"worst err = {worst:.3e} (tol {tol_rel:.0e} relative)")
    return ok, msgs


def check_rate_antivacuity(seepage_frames: list[np.ndarray],
                           rate_frames: list[np.ndarray],
                           depth_frames: list[np.ndarray], min_depth: float):
    """dg1 anti-vacuity (the V3 reproduction): in >= 1 output frame the
    per-step snapshot `seepage` is identically 0 on every wet cell while
    the interval-mean `seepage_rate` is nonzero somewhere — the snapshot
    is a per-step sample and misses the subcycled exchange."""
    hits = 0
    for snap, rate, depth in zip(seepage_frames, rate_frames, depth_frames):
        wet = np.nan_to_num(depth) > min_depth
        if not wet.any():
            continue
        snap_zero_on_wet = not np.nan_to_num(snap)[wet].any()
        rate_active = bool(np.nan_to_num(rate).any())
        if snap_zero_on_wet and rate_active:
            hits += 1
    ok = hits >= 1
    return ok, (f"V3 reproduction: {hits} frame(s) with snapshot seepage == 0 "
                f"on every wet cell while seepage_rate != 0 (need >= 1)")


def check_decimated_subset(full_rows: np.ndarray, dec_rows: np.ndarray,
                           t_end: float, label: str):
    """dg2: every decimated row equals the every-step row at the same time
    BITWISE (all columns), and the decimated table's last row is at t_end.
    Column 0 is time."""
    by_time = {row[0]: row for row in full_rows}
    ok = True
    parts = []
    bad = 0
    for row in dec_rows:
        full = by_time.get(row[0])
        if full is None or not np.array_equal(row, full):
            bad += 1
            ok = False
    parts.append(f"{label}: {len(dec_rows)} decimated rows, {bad} not bitwise "
                 f"rows of the every-step table")
    if abs(float(dec_rows[-1][0]) - t_end) > 1.0e-9:
        ok = False
        parts.append(f"{label}: last row at {dec_rows[-1][0]:g} s, not t_end "
                     f"{t_end:g} s")
    return ok, "; ".join(parts)


def check_row_times(dec_times: np.ndarray, expected: list[float], label: str):
    """dg2 forced-row lane: the decimated row-time set is exactly the
    schedule (t0, the interval multiples, the checkpoint time, t_end)."""
    got = sorted(float(t) for t in dec_times)
    want = sorted(expected)
    ok = (len(got) == len(want)
          and all(abs(a - b) <= 1.0e-9 for a, b in zip(got, want)))
    extra = [t for t in got if not any(abs(t - w) <= 1.0e-9 for w in want)]
    missing = [w for w in want if not any(abs(w - t) <= 1.0e-9 for t in got)]
    return ok, (f"{label}: {len(got)} rows vs {len(want)} expected"
                + (f"; extra {extra[:4]}" if extra else "")
                + (f"; missing {missing[:4]}" if missing else ""))
