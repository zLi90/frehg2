#!/usr/bin/env python3
"""gate_periodic.py — pure references and checks for the Stage-3 T3
periodic-topology gates pd1–pd8 (Stage-3 plan §7).

The pd1 reference is the DISCRETE wind–friction balance of the code's own
steady relation (Momentum.cpp updateVelocityField with the legacy double
drag factor, SweFormulas.hpp relative-wind Garratt stress), solved by
bisection with every coefficient recomputed here (the g9b practice). The
pd2 reference is the exact translation identity roll(run(roll(inputs, k)),
-k) == run(inputs). The HDF5 plumbing lives in run_regression.py; the
§6.3 battery (scripts/test_pd_gates.py) proves the balance against an
independent fixed-point evaluation and every check's ability to fail.
"""

from __future__ import annotations

import numpy as np

GRAVITY = 9.81
AIR_DENSITY = 1.225     # SweFormulas.hpp kAirDensity
WATER_DENSITY = 998.0   # SweFormulas.hpp kWaterDensity (wind + momentum rho)


def cd_garratt(u10: float, cap: float = 3.5e-3) -> float:
    """Garratt (1977), independent of gate_wind (recomputed per plan §6.2)."""
    return min((0.75 + 0.067 * u10) * 1.0e-3, cap)


def manning_cd(n: float, depth: float, thin_layer: float = 0.1) -> float:
    expo = 2.0 / 3.0 if depth < thin_layer else 1.0 / 3.0
    return GRAVITY * n * n / depth**expo


def discrete_balance(u10: float, h: float, n: float, dt: float,
                     cap: float = 3.5e-3, thin_layer: float = 0.1) -> float:
    """The steady along-shore velocity of the DISCRETE scheme on a flat
    periodic channel: no pressure gradient, no advection (u uniform), so
    updateVelocityField's fixed point is
        u = (u + dt tau(u) / (rho_w h)) D(u)^2
    with tau(u) = rho_a Cd (W - u)^2 (the relative wind) and the
    point-implicit factor D = 1 / (1 + dt/2 CD |u| / h) applied twice (the
    preserved legacy double drag, Momentum.cpp). Solved by bisection on
    g(u) = rhs(u) - u, which is positive at 0 and negative at W."""
    cd_w = cd_garratt(u10, cap)
    cd_b = manning_cd(n, h, thin_layer)

    def rhs(u: float) -> float:
        tau = AIR_DENSITY * cd_w * (u10 - u) ** 2
        d = 1.0 / (1.0 + 0.5 * dt * cd_b * abs(u) / h)
        return (u + dt * tau / (WATER_DENSITY * h)) * d * d

    lo, hi = 0.0, u10
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rhs(mid) > mid:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def continuum_balance(u10: float, h: float, n: float,
                      cap: float = 3.5e-3) -> float:
    """The continuum closed form (recorded, not gated): g n^2 u^2 / h^(4/3)
    = rho_a Cd (W - u)^2 / (rho_w h), solved by bisection."""
    cd_w = cd_garratt(u10, cap)

    def f(u: float) -> float:
        return (GRAVITY * n * n * u * u / h**(4.0 / 3.0)
                - AIR_DENSITY * cd_w * (u10 - u) ** 2 / (WATER_DENSITY * h))

    lo, hi = 0.0, u10
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Checks.
# ---------------------------------------------------------------------------

def check_flat(eta: np.ndarray, tol_abs: float, label: str = "eta"):
    """Spread of a field that must be uniform (NaN-masked)."""
    active = eta[np.isfinite(eta)]
    span = float(active.max() - active.min())
    ok = span <= tol_abs
    return ok, f"{label} spread = {span:.3e} (<= {tol_abs:.0e})"


def check_balance(u: np.ndarray, u_ref: float, tol_rel: float):
    """Every along-shore face velocity equals the discrete balance."""
    active = u[np.isfinite(u)]
    worst = float(np.max(np.abs(active - u_ref))) / abs(u_ref)
    ok = worst <= tol_rel
    return ok, f"max |u - u_disc| / u_disc = {worst:.3e} (<= {tol_rel:.0e})"


def check_steady(series: np.ndarray, tol_abs: float):
    span = float(np.max(series) - np.min(series))
    ok = span <= tol_abs
    return ok, f"trailing-window span = {span:.3e} (<= {tol_abs:.0e})"


def check_roll(base: np.ndarray, rolled: np.ndarray, k: int, tol_abs: float,
               label: str):
    """roll(rolled_output, -k) == base_output along the periodic axis
    (axis -1 for (ny, nx) and (ny, nx, nz)-ordered arrays: callers pass
    arrays whose LAST-but-one axis is x for 3-D, so they pre-roll; here we
    just compare against the pre-unrolled array)."""
    a = np.nan_to_num(base, nan=-9.0e99)
    b = np.nan_to_num(rolled, nan=-9.0e99)
    worst = float(np.max(np.abs(a - b)))
    ok = worst <= tol_abs
    return ok, f"{label} roll-{k} max |delta| = {worst:.3e} (<= {tol_abs:.0e})"


def check_seam_crossing(tracer_frames: list[np.ndarray], seam_col: int,
                        floor: float):
    """Anti-vacuity: the tracer actually crosses the seam — the column just
    PAST the seam (wrap target) must receive at least `floor` somewhere
    during the run while having started clean."""
    start = float(np.nanmax(tracer_frames[0][:, seam_col]))
    peak = max(float(np.nanmax(f[:, seam_col])) for f in tracer_frames[1:])
    ok = (start < 0.5 * floor) and (peak >= floor)
    return ok, (f"seam-crossing tracer at wrap column {seam_col}: start "
                f"{start:.3e}, peak {peak:.3e} (floor {floor:.1e})")


def check_surface_ledger(audit: np.ndarray, tol_rel: float):
    """Uncoupled surface ledger (mass_audit: time, volume, rain, evap,
    boundary_outflow, bc_inflow, clamped[, stage_deposit])."""
    d = audit[-1] - audit[0]
    residual = d[1] + d[4] - d[2] + d[3] - d[5] - d[6]
    scale = max(float(np.sum(np.abs(np.diff(audit[:, 4])))), 1.0e-12)
    ok = abs(residual) <= max(tol_rel * scale, 1.0e-12)
    return ok, (f"surface ledger residual {residual:.3e} m^3 vs "
                f"max({tol_rel:.0e} * X, 1e-12) (X = {scale:.3e} m^3)")


def check_linear_head(head: np.ndarray, y: np.ndarray, h_south: float,
                      h_north: float, tol_abs: float):
    """pd5: the confined steady PRESSURE head is uniform in x at every
    (j, k) (the periodic axis carries no gradient), linear in y per layer
    (second differences vanish), hydrostatic in z (the per-layer profiles
    are parallel — the saturated column adds h = H - z_c per layer), and
    drops in the imposed direction. head is (ny, nx, nz), NaN outside."""
    ny, _, nz = head.shape
    ok = True
    parts = []
    worst_x = 0.0
    for j in range(ny):
        for k in range(nz):
            layer = head[j, :, k][np.isfinite(head[j, :, k])]
            if layer.size:
                worst_x = max(worst_x, float(layer.max() - layer.min()))
    if worst_x > tol_abs:
        ok = False
    parts.append(f"max x spread = {worst_x:.3e}")
    profs = np.array([[float(np.nanmean(head[j, :, k])) for j in range(ny)]
                      for k in range(nz)])
    curv = float(np.max(np.abs(np.diff(profs, 2, axis=1))))
    if curv > tol_abs:
        ok = False
    parts.append(f"max |d2h/dy2| = {curv:.3e}")
    # Hydrostatic parallelism: each layer's profile is the one above it
    # shifted by the constant dz (h = H - z_c), so the MEAN-CENTERED
    # profiles must coincide across layers.
    centered = profs - profs.mean(axis=1, keepdims=True)
    parallel = float(np.max(np.abs(centered - centered[0])))
    if parallel > tol_abs:
        ok = False
    parts.append(f"layer-shape parallelism = {parallel:.3e}")
    prof = profs.mean(axis=0)
    direction_ok = (prof[0] - prof[-1]) * (h_south - h_north) > 0.0
    ok = ok and direction_ok
    parts.append(f"profile {prof[0]:.4f} -> {prof[-1]:.4f} "
                 f"({'with' if direction_ok else 'AGAINST'} the head drop)")
    return ok, "; ".join(parts)


def check_zero_x_flux(qx: np.ndarray, qy_scale: float, tol_rel: float):
    """pd5: every x face carries nothing against the y-flux scale."""
    worst = float(np.nanmax(np.abs(qx)))
    ok = worst <= tol_rel * qy_scale
    return ok, (f"max |qx| = {worst:.3e} vs {tol_rel:.0e} * qy scale "
                f"({qy_scale:.3e})")


def check_coupled_exchange(c0: np.ndarray, c1: np.ndarray, floor: float):
    """pd2-coupled anti-vacuity: the tidal seepage exchange actually moved
    surface salt into the soil — the subsurface concentration must change
    by at least `floor` somewhere over the run."""
    mask = np.isfinite(c0) & np.isfinite(c1)
    moved = float(np.max(np.abs(c1[mask] - c0[mask]))) if mask.any() else 0.0
    ok = moved >= floor
    return ok, (f"coupled exchange: max |dc_gw| over the run = {moved:.3e} "
                f"(floor {floor:.1e})")


def check_restart_field(full: np.ndarray, part: np.ndarray, label: str,
                        tol_abs: float | None = None,
                        tol_rel: float | None = None):
    """pd-restart (S3-A6): the restarted field agrees with the continuous
    run to the project's restart class — absolute for the SI-scaled flow
    fields, relative to the field's own magnitude for the tracer. Exactly
    one of tol_abs / tol_rel is given."""
    mask = np.isfinite(full) & np.isfinite(part)
    worst = float(np.max(np.abs(full[mask] - part[mask]))) if mask.any() else 0.0
    if tol_abs is not None:
        ok = worst <= tol_abs
        return ok, f"{label}: max |delta| = {worst:.3e} (<= {tol_abs:.0e} abs)"
    scale = max(float(np.max(np.abs(full[mask]))) if mask.any() else 0.0,
                1.0e-300)
    rel = worst / scale
    ok = rel <= tol_rel
    return ok, (f"{label}: max |delta| = {worst:.3e}, rel {rel:.3e} "
                f"(<= {tol_rel:.0e} rel, scale {scale:.3e})")
