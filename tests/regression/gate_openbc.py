#!/usr/bin/env python3
"""gate_openbc.py — pure references and checks for the Stage-3 T2
open-boundary-closure gates ob1/ob3/ob5 (Stage-3 plan §6).

The ob1 reference is the DISCRETE steady face chain: the code's own update
formulas (Momentum.cpp advection + point-implicit drag with the legacy
double drag factor, WetDry.cpp max-of-two face geometry, the W2.3
conservative entrance across a reservoir face), marched face by face at a
trial discharge Q and shot on Q to hit the downstream stage. Every
coefficient is recomputed here from the configuration numbers — no model
output enters the reference (plan §6.2 style). The HDF5 plumbing lives in
run_regression.py; the §6.3 battery (scripts/test_ob_gates.py) proves the
chain's internal consistency, its agreement with the continuum Bernoulli
limit, and that every check fails on corrupted input.

Face/cell convention (1-based like the kernels): cell 1 is the west stage
cell, cells 2..nx-1 interior, cell nx the east stage cell; face f is the
east face of cell f (u > 0 eastward), so face 1 is the stage cell's east
face — the inflow face whose upstream stencil value the closure selects.
"""

from __future__ import annotations

import numpy as np

GRAVITY = 9.81
WATER_DENSITY = 998.0  # SweFormulas.hpp kWaterDensity (the momentum rho)


# ---------------------------------------------------------------------------
# The ob1 discrete steady chain.
# ---------------------------------------------------------------------------

def manning_cd(n: float, depth: float, thin_layer: float = 0.1) -> float:
    """SweFormulas.hpp manningDrag: g n^2 / h^(1/3) above the thin-layer
    depth, g n^2 / h^(2/3) below (ob1 stays above)."""
    expo = 2.0 / 3.0 if depth < thin_layer else 1.0 / 3.0
    return GRAVITY * n * n / depth**expo


class Chain:
    """One steady march at discharge Q; see chain_solve for the shot."""

    def __init__(self, nx: int, dx: float, dy: float, dt: float, bed: float,
                 n_manning: float, eta_west: float, closure: str,
                 u_ext: float = 0.0):
        if closure not in ("external", "zero_gradient"):
            raise ValueError(closure)
        self.nx, self.dx, self.dy, self.dt = nx, dx, dy, dt
        self.bed, self.n = bed, n_manning
        self.eta_west, self.closure, self.u_ext = eta_west, closure, u_ext

    def march(self, q: float):
        """March faces 1..nx-1; returns (eta_cells, u_faces) with
        eta_cells[1] = the west stage and eta_cells[nx] = the implied
        downstream stage (compared against the prescribed one by the shot).
        Face geometry is the code's: deptx = max(eta_l, eta_r) - bed,
        Vsx = (Vs_l + Vs_r)/2, Aszx = dx dy, CD from the WEST cell's depth
        (Momentum.cpp uses CDx(j, i) — cell i's coefficient at its east
        face). The per-face relation is the steady point of
        updateVelocityField: u = (u - dt advX) D^2 - g dt effh (eta_r -
        eta_l) D, solved for eta_r with a short fixed point on the
        eta_r-dependent geometry."""
        area = self.dx * self.dy
        eta = np.zeros(self.nx + 1)
        u = np.zeros(self.nx)  # u[f], f = 1..nx-1
        eta[1] = self.eta_west
        u_up = None  # the advection upstream value of the current face
        for f in range(1, self.nx):
            eta_l = eta[f]
            eta_r = eta_l  # fixed-point start
            for _ in range(8):
                deptx = max(eta_l, eta_r) - self.bed
                asx = deptx * self.dy
                h_l = eta_l - self.bed
                vs_l = h_l * area
                vs_r = (eta_r - self.bed) * area
                vsx = 0.5 * (vs_l + vs_r)
                effh = asx / vsx
                uf = q / asx
                cd = manning_cd(self.n, h_l)
                d = 1.0 / (1.0 + 0.5 * self.dt * cd * abs(uf) * (area / vsx))
                if f == 1:
                    if self.closure == "external":
                        # W2.3 conservative entrance across the reservoir
                        # face: d(u^2/2)/dx with the reservoir state.
                        adv = (uf * uf - self.u_ext * self.u_ext) / (2.0 * self.dx)
                    else:
                        # zero_gradient: the boundary face carries the
                        # mass-consistent reconstruction u0 = Q / Asx(j,0),
                        # deptx(j,0) = eta_1 - bed (zero-gradient ghost).
                        u0 = q / (h_l * self.dy)
                        adv = (uf / self.dx) * (uf - u0)
                else:
                    adv = (uf / self.dx) * (uf - u_up)
                # Steady updateVelocityField (legacy double drag factor):
                # u = (u - dt adv) D^2 - g dt effh (eta_r - eta_l) D.
                d_eta = ((uf - self.dt * adv) * d * d - uf) / (GRAVITY * self.dt * effh * d)
                eta_r = eta_l + d_eta
            eta[f + 1] = eta_r
            u[f] = q / asx
            u_up = u[f]
        return eta, u

    def residuals(self, q: float) -> np.ndarray:
        """Per-face steady-relation residuals at the march's own solution —
        the battery's internal-consistency check (must be ~0 by
        construction of the fixed point)."""
        eta, u = self.march(q)
        area = self.dx * self.dy
        res = np.zeros(self.nx)
        u_up = None
        for f in range(1, self.nx):
            eta_l, eta_r = eta[f], eta[f + 1]
            deptx = max(eta_l, eta_r) - self.bed
            asx = deptx * self.dy
            h_l = eta_l - self.bed
            vsx = 0.5 * ((eta_l - self.bed) + (eta_r - self.bed)) * area
            effh = asx / vsx
            uf = u[f]
            cd = manning_cd(self.n, h_l)
            d = 1.0 / (1.0 + 0.5 * self.dt * cd * abs(uf) * (area / vsx))
            if f == 1:
                if self.closure == "external":
                    adv = (uf * uf - self.u_ext * self.u_ext) / (2.0 * self.dx)
                else:
                    u0 = q / (h_l * self.dy)
                    adv = (uf / self.dx) * (uf - u0)
            else:
                adv = (uf / self.dx) * (uf - u_up)
            res[f] = uf - ((uf - self.dt * adv) * d * d
                           - GRAVITY * self.dt * effh * (eta_r - eta_l) * d)
            u_up = uf
        return res


def chain_solve(chain: Chain, eta_east: float, q_lo: float = 1.0e-9,
                tol: float = 1.0e-14, iters: int = 200):
    """Shoot on Q so the implied downstream stage equals eta_east.
    Returns (q, eta_cells, u_faces)."""
    def implied(q: float) -> float:
        eta, _ = chain.march(q)
        return eta[chain.nx] - eta_east

    q_hi = q_lo * 2.0
    while implied(q_hi) > 0.0:
        q_hi *= 2.0
        if q_hi > 1.0e6:
            raise RuntimeError("chain_solve: no bracketing discharge")
    lo, hi = q_lo, q_hi
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if implied(mid) > 0.0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= tol * max(hi, 1.0):
            break
    q = 0.5 * (lo + hi)
    eta, u = chain.march(q)
    return q, eta, u


# ---------------------------------------------------------------------------
# Checks (pure; fed by run_regression.py from the HDF5 outputs).
# ---------------------------------------------------------------------------

def check_velocity_profile(u_model: np.ndarray, u_chain: np.ndarray,
                           tol_rel: float):
    """Model face velocities (faces 1..nx-1) against the chain's."""
    scale = float(np.max(np.abs(u_chain)))
    worst = float(np.max(np.abs(u_model - u_chain))) / scale
    ok = worst <= tol_rel
    return ok, f"max |u - u_chain| / max|u_chain| = {worst:.3e} (<= {tol_rel:.0e})"


def check_steadiness(u_frames: np.ndarray, tol_abs: float):
    """Trailing-window span of the inflow-face velocity over the last
    frames (steadiness asserted in-gate, plan §6)."""
    span = float(np.max(u_frames) - np.min(u_frames))
    ok = span <= tol_abs
    return ok, f"trailing-window u span = {span:.3e} m/s (<= {tol_abs:.0e})"


def check_entrance_ratio(eta_stage: float, eta_first: float, u_in: float,
                         ratio_chain: float, tol_abs: float):
    """R = (eta_stage - eta_first_interior) / (u_in^2 / 2g) against the
    chain's own ratio (the seed test's 1.01 observable, now chain-exact)."""
    head = u_in * u_in / (2.0 * GRAVITY)
    ratio = (eta_stage - eta_first) / head
    ok = abs(ratio - ratio_chain) <= tol_abs
    return ok, (f"entrance ratio {ratio:.4f} vs chain {ratio_chain:.4f} "
                f"(|diff| <= {tol_abs})")


def check_lane_contrast(u_a: float, u_c: float, factor: float, r_a: float,
                        r_c: float, r_cap: float):
    """The real negative lane: zero_gradient must run >= factor x faster
    and pay <= r_cap of the entrance head (plan ob1 lane C)."""
    ok = (u_c >= factor * u_a) and (r_c <= r_cap)
    return ok, (f"u_C {u_c:.3f} vs {factor} x u_A {factor * u_a:.3f}; "
                f"R_C {r_c:.3f} (<= {r_cap}) [R_A {r_a:.3f}]")


# ob3 — uncoupled surface ledger (mass_audit columns: time, volume, rain,
# evaporation, boundary_outflow, bc_inflow, clamped).

def check_surface_ledger(audit: np.ndarray, tol_rel: float):
    d = audit[-1] - audit[0]
    residual = d[1] + d[4] - d[2] + d[3] - d[5] - d[6]
    scale = max(float(np.sum(np.abs(np.diff(audit[:, 4])))), 1.0e-12)
    ok = abs(residual) <= max(tol_rel * scale, 1.0e-12)
    return ok, (f"surface ledger residual {residual:.3e} m^3 vs "
                f"max({tol_rel:.0e} * X, 1e-12) (X = {scale:.3e} m^3)")


def prism_per_cycle(audit: np.ndarray, period: float) -> list[float]:
    """Gross stage exchange per tidal cycle = sum |d(boundary_outflow)|
    over each cycle window."""
    t = audit[:, 0]
    out = audit[:, 4]
    prisms = []
    t0 = t[0]
    while t0 + period <= t[-1] + 1.0e-6:
        mask = (t >= t0 - 1.0e-6) & (t <= t0 + period + 1.0e-6)
        prisms.append(float(np.sum(np.abs(np.diff(out[mask])))))
        t0 += period
    return prisms


def check_prism_agreement(prisms_a: list[float], prisms_b: list[float],
                          fr2_max: float, coef: float, floor_rel: float):
    """Between-lane tidal-prism agreement: |dP|/P <= max(coef * Fr^2,
    floor_rel) per cycle — the closure moves the exchange only at the
    momentum-flux order (plan P14)."""
    ok = True
    parts = []
    bound = max(coef * fr2_max, floor_rel)
    for n, (a, b) in enumerate(zip(prisms_a, prisms_b)):
        rel = abs(a - b) / max(a, 1.0e-12)
        good = rel <= bound
        ok = ok and good
        parts.append(f"cycle {n}: {rel:.3e}{'' if good else ' FAIL'}")
    return ok, (f"prism diff vs bound {bound:.3e} (Fr^2 {fr2_max:.3e}): "
                + "; ".join(parts))


# ob5 — kinetic energy and stage-face KE influx.

def domain_ke(uu: np.ndarray, vv: np.ndarray, depth: np.ndarray,
              dx: float, dy: float) -> float:
    """0.5 rho h |u|^2 summed over cells, with cell speeds from face
    averages (adequate for a monotonicity check)."""
    uc = 0.5 * (uu + np.roll(uu, 1, axis=1))
    vc = 0.5 * (vv + np.roll(vv, 1, axis=0))
    ke = 0.5 * WATER_DENSITY * depth * (uc * uc + vc * vc) * dx * dy
    return float(np.nansum(ke))


def stage_ke_influx(vv_row: np.ndarray, depth_row: np.ndarray, dx: float,
                    dt_frame: float) -> float:
    """Incoming KE through the north stage faces over one frame interval:
    0.5 rho v^2 |v| deptx dx for inflow (v < 0, southward), summed [J]."""
    v = vv_row
    inflow = np.where(v < 0.0, -v, 0.0)
    return float(np.nansum(0.5 * WATER_DENSITY * v * v * inflow
                           * depth_row * dx) * dt_frame)


def check_ke_decay(ke_series: np.ndarray, t_settle_frames: int, tol_rise: float):
    """KE monotone decreasing after the seed's initial adjustment, within a
    per-frame rise tolerance (rounding/exchange ripple)."""
    tail = ke_series[t_settle_frames:]
    rises = np.diff(tail)
    worst = float(np.max(rises)) if rises.size else 0.0
    ok = worst <= tol_rise * ke_series[0]
    return ok, (f"worst KE rise {worst:.3e} J per frame vs "
                f"{tol_rise:.0e} * KE0 ({ke_series[0]:.3e} J)")


def check_ke_influx(influx_ext: float, ke0: float, tol_rel: float):
    ok = influx_ext <= tol_rel * ke0
    return ok, (f"external stage KE influx {influx_ext:.3e} J vs "
                f"{tol_rel:.0e} * KE0 ({ke0:.3e} J)")


def check_influx_contrast(influx_zg: float, influx_ext: float, factor: float,
                          floor_j: float):
    """Anti-vacuity + the real negative: zero_gradient must import >=
    factor x the external lane's KE and >= an absolute floor (a case whose
    zero_gradient lane imports nothing proves nothing)."""
    ok = (influx_zg >= factor * max(influx_ext, 1.0e-30)) and (influx_zg >= floor_j)
    return ok, (f"zero_gradient influx {influx_zg:.3e} J vs {factor} x external "
                f"{influx_ext:.3e} J and floor {floor_j:.1e} J")
