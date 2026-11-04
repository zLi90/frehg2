"""gate_pc.py — pure reference and check functions for the Stage-3 T11
predictor–corrector consistency gates (plan §12C, pc0–pc8). No HDF5 I/O;
scripts/test_pc_gates.py feeds corrupted synthetic inputs to every checker
(the §6.3 negative battery).

gw_mass_audit columns as in gate_realloc: time, volume, boundary_in,
ss_storage, realloc, realloc_dropped, vloss (cumulative except volume).
"""

from __future__ import annotations

import numpy as np

TIME, VOLUME, BOUNDARY_IN, SS_STORAGE, REALLOC, REALLOC_DROPPED, VLOSS = range(7)


def deltas(audit: np.ndarray) -> np.ndarray:
    return audit[-1] - audit[0]


def throughput(audit: np.ndarray) -> float:
    """T = sum |boundary_in increments| (the rl ledger scale)."""
    return float(np.sum(np.abs(np.diff(audit[:, BOUNDARY_IN]))))


def churn(audit: np.ndarray) -> float:
    """Gross restore churn = created + dropped = d(realloc) +
    2 d(realloc_dropped) (the rl2/S3-A4 definition)."""
    d = deltas(audit)
    return float(d[REALLOC] + 2.0 * d[REALLOC_DROPPED])


def check_churn_cap(audit: np.ndarray, cap_rel: float,
                    label: str = "churn") -> tuple[bool, str]:
    """pc1/pc2 (K1 at the lane scale): the saturated-contact restore books
    at most cap_rel of the lane's own throughput — the INVERSE of the old
    rl2 anti-vacuity floor, which encoded M1 (S3-A16). Red on T10: rl2
    churn/T = 2.91; the report's patched binary measured 0.006."""
    c = churn(audit)
    scale = max(throughput(audit), 1.0e-12)
    ratio = c / scale
    ok = ratio <= cap_rel
    return ok, (f"{label}: gross churn {c:.3e} m^3 = {ratio:.3e} of T = "
                f"{scale:.3e} (<= {cap_rel:.0e})")


def check_hetero_activity(n_sat_contrast_faces: int, floor_faces: int,
                          churn_gross: float) -> tuple[bool, str]:
    """pc1's replacement anti-vacuity: the lane holds >= floor saturated
    faces with a Ks ratio > 10 (the regime the old floor proxied), and the
    restore actually fires (churn strictly > 0 — it is the scheme's
    consistency mechanism, small but never exactly zero on a live lane)."""
    ok = (n_sat_contrast_faces >= floor_faces) and (churn_gross > 0.0)
    return ok, (f"anti-vacuity: {n_sat_contrast_faces} saturated Ks-contrast "
                f"faces (>= {floor_faces}), gross churn {churn_gross:.3e} (> 0)")


def count_sat_contrast_faces(ks: np.ndarray, wc: np.ndarray, wcs: np.ndarray,
                             ratio_floor: float = 10.0,
                             sat_frac: float = 0.9999) -> int:
    """Count interior faces (all three axes) whose flanks are both active,
    both saturated (theta >= sat_frac * theta_s), with a saturated-Ks ratio
    above ratio_floor. Arrays are (ny, nx, nz) with NaN on inactive cells."""
    n = 0
    sat = np.isfinite(wc) & (wc >= sat_frac * wcs)
    for axis in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[axis] = slice(None, -1)
        b[axis] = slice(1, None)
        ka, kb = ks[tuple(a)], ks[tuple(b)]
        both = sat[tuple(a)] & sat[tuple(b)] & np.isfinite(ka) & np.isfinite(kb)
        lo = np.minimum(ka, kb)
        hi = np.maximum(ka, kb)
        contrast = both & (lo > 0.0) & (hi / np.where(lo > 0.0, lo, 1.0) > ratio_floor)
        n += int(np.sum(contrast))
    return n


def series_resistance(dz: np.ndarray, ks: np.ndarray) -> float:
    """pc3: the exact one-dimensional series resistance Sigma(L/K) of a
    saturated stack whose boundary faces sit at half spacing — i.e. each
    cell contributes its own full thickness over its own Ks (the halves
    pair across interfaces and with the half-spacing edge faces)."""
    return float(np.sum(np.asarray(dz) / np.asarray(ks)))


def check_series_flux(q_measured: float, dh_total: float, dz: np.ndarray,
                      ks: np.ndarray, area: float, tol_rel: float,
                      label: str) -> tuple[bool, str]:
    """pc3 (harmonic lane): the steady face flux equals the exact series
    value dH_total / Sigma(L/K) times the face area."""
    q_exact = dh_total / series_resistance(dz, ks) * area
    rel = abs(q_measured - q_exact) / max(abs(q_exact), 1.0e-300)
    ok = rel <= tol_rel
    return ok, (f"{label}: flux {q_measured:.6e} vs exact {q_exact:.6e} m^3/s "
                f"(rel {rel:.3e} <= {tol_rel:.0e})")


def check_face_resistance(dh_face: float, gravity_term: float, q: float,
                          delta: float, k_exact: float, tol_rel: float,
                          label: str) -> tuple[bool, str]:
    """pc3 (per interior face): the measured face resistance (driving head
    over flux per unit area) equals delta / K_exact. The driving head is
    the psi drop plus the face's gravity term (delta for a fresh z face,
    0 for a flat x face) — self-contained in the interior, independent of
    the edge-leg conventions (S3-A17: the side head-edge legs have their
    own predictor/corrector split, measured 3x, OUT of pc3's scope)."""
    r_meas = abs(dh_face + gravity_term) / max(abs(q), 1.0e-300)
    r_exact = delta / k_exact
    rel = abs(r_meas - r_exact) / r_exact
    ok = rel <= tol_rel
    return ok, (f"{label}: face R {r_meas:.6e} vs exact {r_exact:.6e} s "
                f"(rel {rel:.3e} <= {tol_rel:.0e})")


def check_steady(q_frames: list[float], tol_rel: float) -> tuple[bool, str]:
    """pc3/pc6b: the last two sampled fluxes agree (the run reached steady
    state before the assert)."""
    a, b = q_frames[-2], q_frames[-1]
    rel = abs(a - b) / max(abs(b), 1.0e-300)
    ok = rel <= tol_rel
    return ok, f"steady: last two fluxes {a:.6e} / {b:.6e} (rel {rel:.3e} <= {tol_rel:.0e})"


def check_wall_flux(q_wall_max: float, leak_arith: float,
                    label: str = "wall") -> tuple[bool, str]:
    """pc4 (harmonic lane): no water crosses a Ks = 0 wall — the wall-face
    flux is exactly zero (the harmonic face conductance of a zero flank is
    zero by construction, not small). The arithmetic lane's leak is passed
    in for the record line only."""
    ok = q_wall_max == 0.0
    return ok, (f"{label}: max |q| through the Ks = 0 wall = {q_wall_max!r} "
                f"(must be exactly 0.0; arithmetic lane leaks {leak_arith:.3e} m^3/s)")


def check_front_advance(wc_profile: np.ndarray, k_interface: int,
                        theta_floor: float) -> tuple[bool, str]:
    """pc5 (anti-stall): the wetting front crossed into the low-K layer —
    the first cell below the interface wetted above the floor."""
    val = float(wc_profile[k_interface])
    ok = val > theta_floor
    return ok, (f"anti-stall: theta just below the K interface = {val:.4f} "
                f"(> {theta_floor})")


def check_column_restore(created: float, dropped: float, moved: float,
                         tol_rel: float, label: str) -> tuple[bool, str]:
    """pc6 (K1/K2 at the column scale): the restore books at most tol_rel
    of the volume the column actually exchanged (infiltrated or exfiltrated)."""
    gross = abs(created) + abs(dropped)
    scale = max(abs(moved), 1.0e-12)
    ratio = gross / scale
    ok = ratio <= tol_rel
    return ok, (f"{label}: restore churn {gross:.3e} m^3 = {ratio:.3e} of the "
                f"exchanged {scale:.3e} m^3 (<= {tol_rel:.0e})")


def check_ledger_constant(series: np.ndarray, tol_rel: float,
                          label: str) -> tuple[bool, str]:
    """pc7: a conserved ledger (Sigma s theta V, Sigma (1 + beta s) theta V)
    holds its first value to tol_rel over every frame."""
    s0 = float(series[0])
    worst = float(np.max(np.abs(series - s0)))
    rel = worst / max(abs(s0), 1.0e-300)
    ok = rel <= tol_rel
    return ok, (f"{label}: worst drift {worst:.3e} of {s0:.6e} "
                f"(rel {rel:.3e} <= {tol_rel:.0e})")


def ledgers(theta: np.ndarray, conc: np.ndarray, cell_vol: np.ndarray,
            beta: float) -> tuple[float, float, float]:
    """pc7 frame ledgers from the fields: (Sigma theta V, Sigma s theta V,
    Sigma (1 + beta s) theta V) over active (finite) cells."""
    m = np.isfinite(theta)
    tv = theta[m] * cell_vol[m]
    s = np.where(np.isfinite(conc[m]), conc[m], 0.0)
    return (float(np.sum(tv)), float(np.sum(s * tv)),
            float(np.sum((1.0 + beta * s) * tv)))


def check_rest(q_max: float, ks: float, face_area: float, tol_rel: float,
               label: str) -> tuple[bool, str]:
    """pc8 (S1): the slab stays at rest — the largest face flux after
    settling is tol_rel of the conductivity scale Ks * A_face. Red on T10
    (salty lane): the S1 force K (rho/rho0 - 1) sin(a) drives a persistent
    flux of ~2.6e-3 * Ks * A at 35 psu on the 10 % slope."""
    scale = ks * face_area
    rel = q_max / scale
    ok = rel <= tol_rel
    return ok, (f"{label}: max |q| = {q_max:.3e} m^3/s = {rel:.3e} of "
                f"Ks*A = {scale:.3e} (<= {tol_rel:.0e})")
