#!/usr/bin/env python3
"""gate_realloc.py — pure check functions for the Stage-3 T1 reallocation
gates rl2/rl3/rl4 (Stage-3 plan §5).

Every reference is an exact identity recomputed in-script from the audit
tables (plan §6.2 style: self-contained, no stored goldens). The HDF5
plumbing lives in run_regression.py; keeping the criteria pure lets the
negative battery (scripts/test_rl_gates.py) prove every check can fail on
corrupted synthetic input.

/monitor/gw_mass_audit columns (docs/user-guide/output.md): time, volume,
boundary_in, ss_storage, realloc, realloc_dropped, vloss[, ...]. All
cumulative except volume; deltas are taken against the first row.
"""

from __future__ import annotations

import numpy as np

# gw_mass_audit column indices.
TIME, VOLUME, BOUNDARY_IN, SS_STORAGE, REALLOC, REALLOC_DROPPED, VLOSS = range(7)
# transport_audit column indices (output.md).
TR_SUBS_MASS, TR_SUBS_ANCHOR = 2, 10


def deltas(audit: np.ndarray) -> np.ndarray:
    """Column deltas last-minus-first (the cumulative columns' run totals)."""
    return audit[-1] - audit[0]


def throughput(audit: np.ndarray) -> float:
    """T = sum |boundary_in increments| — the ledger scale of the rl bounds."""
    return float(np.sum(np.abs(np.diff(audit[:, BOUNDARY_IN]))))


def check_ledger(audit: np.ndarray, tol_rel: float, tol_abs: float = 1.0e-12):
    """Cumulative identity dV = d(boundary_in) - d(ss_storage) + d(realloc)
    - d(vloss). Bound max(tol_rel * T, tol_abs): the absolute floor covers
    low-throughput runs whose residual is rounding accumulated over the
    audit rows, not a leak (rl3's first execution measured 4e-13 m^3
    against T = 0.19 m^3 over 864 rows)."""
    d = deltas(audit)
    residual = d[VOLUME] - (d[BOUNDARY_IN] - d[SS_STORAGE] + d[REALLOC] - d[VLOSS])
    scale = max(throughput(audit), 1.0e-12)
    bound = max(tol_rel * scale, tol_abs)
    ok = abs(residual) <= bound
    return ok, (f"ledger residual {residual:.3e} m^3 vs max({tol_rel:.0e} * T, "
                f"{tol_abs:.0e}) (T = {scale:.3e} m^3)")


def check_realloc_exactly_zero(audit: np.ndarray):
    """off mode: nothing is transferred, created, or dropped — exact zeros,
    not small numbers (the plan's W1.2 contract)."""
    re = float(audit[-1, REALLOC] - audit[0, REALLOC])
    dr = float(audit[-1, REALLOC_DROPPED] - audit[0, REALLOC_DROPPED])
    ok = (re == 0.0) and (dr == 0.0)
    return ok, f"realloc = {re!r}, realloc_dropped = {dr!r} (must be exactly 0.0)"


def check_vloss(audit: np.ndarray, tol_rel: float):
    """off mode: the finalize clamp's gross volume stays a rounding-class
    fraction of the throughput."""
    d = deltas(audit)
    scale = max(throughput(audit), 1.0e-12)
    ok = abs(d[VLOSS]) <= tol_rel * scale
    return ok, f"|vloss| {abs(d[VLOSS]):.3e} m^3 vs {tol_rel:.0e} * T ({scale:.3e})"


def check_theta_bounds(wc: np.ndarray, theta_r_min: float, theta_s_max: float,
                       slack: float):
    """theta within the soil table's global [min theta_r, max theta_s] on
    active (non-NaN) cells."""
    active = wc[np.isfinite(wc)]
    lo = float(active.min())
    hi = float(active.max())
    ok = (lo >= theta_r_min - slack) and (hi <= theta_s_max + slack)
    return ok, (f"theta in [{lo:.6f}, {hi:.6f}] vs table "
                f"[{theta_r_min}, {theta_s_max}] (slack {slack:.0e})")


def check_percell_sums(created_sum: float, dropped_sum: float, vloss_sum: float,
                       audit: np.ndarray, tol_rel: float):
    """The per-cell cumulative diagnostics integrate to the audit columns:
    sum(created - dropped) = d(realloc), sum(dropped) = d(realloc_dropped),
    sum(vloss) = d(vloss) (plan W1.3; exact by construction)."""
    d = deltas(audit)
    gross = abs(d[REALLOC]) + 2.0 * abs(d[REALLOC_DROPPED]) + abs(d[VLOSS])
    scale = max(1.0, gross)
    checks = [
        ("created - dropped vs realloc", created_sum - dropped_sum, d[REALLOC]),
        ("dropped vs realloc_dropped", dropped_sum, d[REALLOC_DROPPED]),
        ("vloss cells vs vloss column", vloss_sum, d[VLOSS]),
    ]
    ok = True
    parts = []
    for name, a, b in checks:
        good = abs(a - b) <= tol_rel * scale
        ok = ok and good
        parts.append(f"{name}: {a:.6e} vs {b:.6e}{'' if good else ' MISMATCH'}")
    return ok, "; ".join(parts)


def check_substep_ratio(count_off: int, count_legacy: int, bound: float):
    """off must not slow the adaptive controller beyond the stated factor."""
    ratio = count_off / max(count_legacy, 1)
    ok = ratio <= bound
    return ok, f"substeps off/legacy = {count_off}/{count_legacy} = {ratio:.3f} (<= {bound})"


def check_legacy_churn(audit_legacy: np.ndarray, floor_rel: float):
    """Anti-vacuity (plan §2.3): the legacy lane must actually exercise
    reallocation. Gross churn = created + dropped = d(realloc) +
    2 d(realloc_dropped), floored at floor_rel of the lane's own
    throughput (rl2 authoring measurement: churn/T = 2.9; floor 0.1 gives
    29x headroom while still failing any case too quiet to reallocate)."""
    d = deltas(audit_legacy)
    churn = d[REALLOC] + 2.0 * d[REALLOC_DROPPED]
    floor = floor_rel * throughput(audit_legacy)
    ok = churn >= floor
    return ok, f"legacy gross churn {churn:.3e} m^3 vs floor {floor:.3e} m^3"


def check_steady(first: np.ndarray, last: np.ndarray, tol_abs: float):
    """Trailing-window steadiness: max |theta change| over the window."""
    diff = np.abs(last - first)
    span = float(np.nanmax(diff))
    ok = span <= tol_abs
    return ok, f"trailing-window max |d theta| = {span:.3e} (<= {tol_abs:.0e})"


def check_uniform_flux(qz: np.ndarray, q_ref: float, tol_rel: float):
    """Steady column: every active face carries the top flux."""
    active = qz[np.isfinite(qz)]
    worst = float(np.max(np.abs(active - q_ref))) / abs(q_ref)
    ok = worst <= tol_rel
    return ok, f"max |qz - q_top| / |q_top| = {worst:.3e} (<= {tol_rel:.0e})"


def check_anchor(tr_disabled: np.ndarray, tr_legacy: np.ndarray, s_max: float,
                 gw_disabled: np.ndarray, tol_rel: float, improvement: float,
                 legacy_floor_rel: float):
    """rl4: with nothing created or dropped the reallocation contribution to
    the subsurface re-anchor vanishes; what remains is the coupled-exchange
    basis plus the scalar on the finalize clamp (first-execution measurement:
    legacy -38.7, disabled +1.32 over M = 5.8e4 — a 29x reduction). Asserts
    (a) |d anchor_disabled| <= tol_rel * M_subs, (b) |d anchor_disabled| <=
    |d anchor_legacy| / improvement, and (c) anti-vacuity: |d anchor_legacy|
    >= legacy_floor_rel * M_subs — a case whose legacy anchor is already
    negligible proves nothing about the reallocation contribution."""
    d_dis = float(tr_disabled[-1, TR_SUBS_ANCHOR] - tr_disabled[0, TR_SUBS_ANCHOR])
    d_leg = float(tr_legacy[-1, TR_SUBS_ANCHOR] - tr_legacy[0, TR_SUBS_ANCHOR])
    m_subs = float(tr_disabled[-1, TR_SUBS_MASS])
    d_vloss = float(gw_disabled[-1, VLOSS] - gw_disabled[0, VLOSS])
    ok = (abs(d_dis) <= tol_rel * abs(m_subs) and
          abs(d_dis) <= abs(d_leg) / improvement and
          abs(d_leg) >= legacy_floor_rel * abs(m_subs))
    return ok, (f"|d anchor| disabled {abs(d_dis):.3e} vs legacy {abs(d_leg):.3e} "
                f"(M_subs {m_subs:.3e}, s_max {s_max:.3g}, bound "
                f"{tol_rel:.0e}*M and legacy/{improvement:g}, legacy floor "
                f"{legacy_floor_rel:.0e}*M; |d vloss| {abs(d_vloss):.3e})")
