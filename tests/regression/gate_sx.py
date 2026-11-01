"""gate_sx.py — pure reference and check functions for the Stage-3 T10
solve-consistent surface transport gates (plan §12B, sx0–sx5). No HDF5
I/O; scripts/test_sx_gates.py feeds corrupted synthetic inputs to every
checker (the §6.3 negative battery)."""

from __future__ import annotations

import numpy as np


def telescope_residuals(d_eta: np.ndarray, cell_area: float,
                        qx: np.ndarray, qy: np.ndarray) -> np.ndarray:
    """sx0/C1: per-cell cellArea * d_eta + div(Q) for interior cells.
    d_eta is (ny, nx); qx[(j, i)] is the eastward face volume between
    cells i and i+1 with a leading west slot (ny, nx+1); qy likewise
    (ny+1, nx). Zero to the solve residual on source-free cells."""
    div = (qx[:, 1:] - qx[:, :-1]) + (qy[1:, :] - qy[:-1, :])
    return cell_area * d_eta + div


def check_uniform(frames: list[np.ndarray], wet: list[np.ndarray],
                  c0: float, tol_abs: float, label: str) -> tuple[bool, str]:
    """J1: every wet-cell concentration equals c0 to tol_abs."""
    worst = 0.0
    n = 0
    for c, w in zip(frames, wet):
        if np.any(w):
            worst = max(worst, float(np.max(np.abs(c[w] - c0))))
            n += int(np.sum(w))
    ok = worst <= tol_abs
    return ok, (f"{label}: max |c - {c0}| over {n} wet samples = {worst:.3e} "
                f"(<= {tol_abs:.0e})")


def check_anchor(anchor_final: float, mass_scale: float, tol_rel: float,
                 label: str) -> tuple[bool, str]:
    """J2: the cumulative surf_anchor stays at the solve-residual class."""
    rel = abs(anchor_final) / mass_scale if mass_scale > 0.0 else abs(anchor_final)
    ok = rel <= tol_rel
    return ok, (f"{label}: |surf_anchor| = {abs(anchor_final):.3e} = {rel:.3e} "
                f"of M {mass_scale:.3e} (<= {tol_rel:.0e})")


def check_extrema(frames: list[np.ndarray], wet: list[np.ndarray],
                  lo: float, hi: float, tol: float) -> tuple[bool, str]:
    """sx4: wet-cell values stay inside the initial range (no new extrema)."""
    worst_hi = -np.inf
    worst_lo = np.inf
    for c, w in zip(frames, wet):
        if np.any(w):
            worst_hi = max(worst_hi, float(np.max(c[w])))
            worst_lo = min(worst_lo, float(np.min(c[w])))
    ok = worst_hi <= hi + tol and worst_lo >= lo - tol
    return ok, (f"wet extrema [{worst_lo:.6f}, {worst_hi:.6f}] inside "
                f"[{lo}, {hi}] (+/- {tol:.0e})")


def check_total(series: np.ndarray, tol_rel: float, label: str) -> tuple[bool, str]:
    """Mass conservation: the total (mass + stores) drift per row."""
    scale = max(abs(float(series[0])), 1.0e-30)
    worst = float(np.max(np.abs(series - series[0])))
    ok = worst <= tol_rel * scale
    return ok, (f"{label}: worst drift {worst / scale:.3e} of {scale:.4e} "
                f"(<= {tol_rel:.0e})")
