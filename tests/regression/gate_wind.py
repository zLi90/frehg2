#!/usr/bin/env python3
"""gate_wind.py — pure closed forms and check functions for the v2 Q6
wind gates g9 (steady setup: nonlinear / linear / sloping-bottom) and
g10 (Merian seiche relaxation), plan §5.3.

Every reference is computed in-script (plan §6.2: the Q6 gates are fully
self-contained). The §6.3 rule-1 cross-checks live in
scripts/test_g910_gates.py: the analytic-integral and numeric-quadrature
volume constraints agree to 1e-9, and the sloping-bottom quadrature
reproduces the flat-bottom closed form. The HDF5 plumbing lives in
run_regression.py; keeping the criteria pure lets the negative battery
prove every check can fail on synthetic data.

Constants mirror src/swe/SweFormulas.hpp (kAirDensity 1.225,
kWaterDensity 998) — the gate recomputes tau/rho_w from the configured
Cd with the code's own densities so no case constant is duplicated.
"""

from __future__ import annotations

import numpy as np

AIR_DENSITY = 1.225     # SweFormulas.hpp kAirDensity
WATER_DENSITY = 998.0   # SweFormulas.hpp kWaterDensity
GRAVITY = 9.81


# ---------------------------------------------------------------------------
# Cd(U10) laws — independent python references for the unit battery and the
# g9(b) offline Cd (plan §5.2; published formulas).
# ---------------------------------------------------------------------------

def cd_garratt(u10: float, cap: float = 3.5e-3) -> float:
    """Garratt (1977): Cd = (0.75 + 0.067 U10) 1e-3, capped."""
    return min((0.75 + 0.067 * u10) * 1.0e-3, cap)


def cd_smith_banke(u10: float, cap: float = 3.5e-3) -> float:
    """Smith & Banke (1975): Cd = (0.63 + 0.066 U10) 1e-3, capped."""
    return min((0.63 + 0.066 * u10) * 1.0e-3, cap)


def cd_wu(u10: float, cap: float = 3.5e-3) -> float:
    """Wu (1982): Cd = (0.8 + 0.065 U10) 1e-3, capped."""
    return min((0.8 + 0.065 * u10) * 1.0e-3, cap)


def cd_large_pond(u10: float, cap: float = 3.5e-3) -> float:
    """Large & Pond (1981): 1.2e-3 below 11 m/s, (0.49 + 0.065 U10) 1e-3
    above, capped. The published law is piecewise and carries a 5e-6
    jump at the 11 m/s breakpoint; kept faithful."""
    cd = 1.2e-3 if u10 < 11.0 else (0.49 + 0.065 * u10) * 1.0e-3
    return min(cd, cap)


def kinematic_stress(cd: float, u10: float) -> float:
    """tau / rho_w [m^2/s^2] for wind speed u10 against still water."""
    return AIR_DENSITY * cd * u10 * u10 / WATER_DENSITY


# ---------------------------------------------------------------------------
# Closed forms.
# ---------------------------------------------------------------------------

def setup_flat(taup: float, length: float, rest_depth: float, x):
    """Steady setup over a flat bottom (Dean & Dalrymple Eq. 5.96 with
    volume conservation): h(x)^2 = C + 2 taup x / g, C from
    integral(h) = rest_depth * length (analytic integral + bisection).
    Returns h at the requested x locations."""
    a = 2.0 * taup / GRAVITY

    def volume(c):
        return (2.0 / (3.0 * a)) * ((c + a * length) ** 1.5 - c ** 1.5)

    lo, hi = 1.0e-12, (rest_depth + np.sqrt(a * length)) ** 2 + 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if volume(mid) > rest_depth * length:
            hi = mid
        else:
            lo = mid
    c = 0.5 * (lo + hi)
    return np.sqrt(c + a * np.asarray(x, dtype=float))


def setup_slope(taup: float, length: float, bed_at, n: int = 20000):
    """Steady setup over an arbitrary bed b(x) (the g9(c) quadrature):
    d eta/dx = taup / (g (eta - b(x))) integrated by RK4, with eta(0)
    bisected so the volume matches the rest volume. Returns (x, eta)."""
    x = np.linspace(0.0, length, n + 1)
    bed = bed_at(x)

    def march(eta0):
        eta = np.empty_like(x)
        eta[0] = eta0
        dxs = x[1] - x[0]

        def f(xx, ee):
            return taup / (GRAVITY * (ee - bed_at(xx)))

        for i in range(n):
            k1 = f(x[i], eta[i])
            k2 = f(x[i] + 0.5 * dxs, eta[i] + 0.5 * dxs * k1)
            k3 = f(x[i] + 0.5 * dxs, eta[i] + 0.5 * dxs * k2)
            k4 = f(x[i] + dxs, eta[i] + dxs * k3)
            eta[i + 1] = eta[i] + dxs * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        return eta

    rest_volume = np.trapz(-bed, x)
    lo = float(bed[0]) + 1.0e-4  # h(0) just above dry
    hi = float(-bed.min())       # far above any plausible setup
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        eta = march(mid)
        if np.trapz(eta - bed, x) > rest_volume:
            hi = mid
        else:
            lo = mid
    return x, march(0.5 * (lo + hi))


def merian_period(length: float, depth: float) -> float:
    """Fundamental closed-basin seiche period T = 2L/sqrt(gH)."""
    return 2.0 * length / np.sqrt(GRAVITY * depth)


# ---------------------------------------------------------------------------
# Checks.
# ---------------------------------------------------------------------------

def check_steady(times, values, window_fraction: float, tol_abs: float,
                 label: str):
    """The monitored signal must be steady: its span over the trailing
    window_fraction of the run stays within tol_abs."""
    t = np.asarray(times, dtype=float)
    v = np.asarray(values, dtype=float)
    tail = v[t >= t[-1] * (1.0 - window_fraction)]
    span = float(np.max(tail) - np.min(tail))
    ok = span <= tol_abs
    msg = (f"{label} steadiness: trailing span {span:.3e} m "
           f"(allowed {tol_abs:.1e}) {'ok' if ok else 'FAIL'}")
    return ok, [msg]


def check_setup(delta_model: float, delta_analytic: float, tol_rel: float,
                label: str):
    """End-to-end setup vs the closed form, relative to the analytic."""
    rel = abs(delta_model - delta_analytic) / abs(delta_analytic)
    ok = rel <= tol_rel
    msg = (f"{label}: end-to-end setup {delta_model:.5f} m vs analytic "
           f"{delta_analytic:.5f} m (rel {rel:.2%}, allowed {tol_rel:.1%}) "
           f"{'ok' if ok else 'FAIL'}")
    return ok, [msg]


def check_profile_rms(eta_model, eta_reference, signal_range: float,
                      tol_fraction: float, label: str):
    """Profile RMS vs the reference as a fraction of the end-to-end
    signal range (the g9(c) gated norm); the max deviation is reported
    alongside as a record — it concentrates in the shallow-tip cell and
    converges first order (see the case README)."""
    dev = np.asarray(eta_model, dtype=float) - np.asarray(eta_reference, dtype=float)
    rms = float(np.sqrt(np.mean(dev * dev)))
    limit = tol_fraction * signal_range
    ok = rms <= limit
    msgs = [(f"{label} profile RMS: {rms:.3e} m (allowed {limit:.3e} = "
             f"{tol_fraction:.0%} of the {signal_range:.4f} m range) "
             f"{'ok' if ok else 'FAIL'}"),
            (f"recorded (not gated): {label} max |eta - reference| = "
             f"{float(np.max(np.abs(dev))):.3e} m (the shallow-tip cell)")]
    return ok, msgs


def check_profile(eta_model, eta_reference, signal_range: float,
                  tol_fraction: float, label: str):
    """Max profile deviation vs the reference, as a fraction of the
    end-to-end signal range."""
    dev = float(np.max(np.abs(np.asarray(eta_model) - np.asarray(eta_reference))))
    limit = tol_fraction * signal_range
    ok = dev <= limit
    msg = (f"{label} profile: max |eta - reference| = {dev:.3e} m "
           f"(allowed {limit:.3e} = {tol_fraction:.0%} of the "
           f"{signal_range:.4f} m range) {'ok' if ok else 'FAIL'}")
    return ok, [msg]


def seiche_metrics(times, eta, t_release: float):
    """(period, amplitudes) from the free oscillation after t_release:
    the fundamental period is twice the mean interval between zero
    crossings of eta - mean(eta); amplitudes are successive |extrema|."""
    t = np.asarray(times, dtype=float)
    v = np.asarray(eta, dtype=float)
    sel = t > t_release
    t, v = t[sel], v[sel]
    v = v - np.mean(v)
    sign = np.sign(v)
    flips = np.where(np.diff(sign) != 0)[0]
    crossings = []
    for i in flips:
        # linear interpolation of the crossing time
        frac = v[i] / (v[i] - v[i + 1])
        crossings.append(t[i] + frac * (t[i + 1] - t[i]))
    crossings = np.asarray(crossings)
    if crossings.size < 4:
        return float("nan"), np.array([])
    period = 2.0 * float(np.mean(np.diff(crossings)))
    # peak amplitude between successive crossings
    amps = []
    for a, b in zip(crossings[:-1], crossings[1:]):
        seg = np.abs(v[(t > a) & (t < b)])
        if seg.size:
            amps.append(float(np.max(seg)))
    return period, np.asarray(amps)


def check_seiche(period_model: float, period_merian: float, tol_rel: float,
                 amplitudes, courant: float):
    """g10: fundamental period within tol_rel of Merian; the per-period
    amplitude decay is recorded (scheme dissipation), not gated."""
    msgs = []
    ok = True
    if not np.isfinite(period_model):
        return False, ["seiche: too few zero crossings to estimate a period FAIL"]
    rel = abs(period_model - period_merian) / period_merian
    line = (f"seiche period {period_model:.1f} s vs Merian {period_merian:.1f} s "
            f"(rel {rel:.2%}, allowed {tol_rel:.1%}; Cr = {courant:.2f})")
    if rel > tol_rel:
        ok = False
        line += " FAIL"
    msgs.append(line)
    if amplitudes is not None and len(amplitudes) >= 3:
        # decay per full period = ratio across two half-period amplitudes
        ratios = amplitudes[2:] / amplitudes[:-2]
        msgs.append(f"recorded (not gated): amplitude decay per period "
                    f"{float(1.0 - np.mean(ratios)):.2%} (theta-scheme dissipation)")
    return ok, msgs


# ---------------------------------------------------------------------------
# The 8-orientation battery (v2 §8.1/§9): each steady eta field must map
# onto the reference orientation's field under the corresponding
# rotation/reflection of the square basin.
# ---------------------------------------------------------------------------

ORIENTATIONS = {
    #   name: (u10, v10, transform mapping this run's field back onto E's)
    "E": (1.0, 0.0, lambda f: f),
    "N": (0.0, 1.0, lambda f: f.T),
    "W": (-1.0, 0.0, lambda f: f[:, ::-1]),
    "S": (0.0, -1.0, lambda f: f.T[:, ::-1]),
    "NE": (1.0, 1.0, lambda f: f),
    "NW": (-1.0, 1.0, lambda f: f[:, ::-1]),
    "SW": (-1.0, -1.0, lambda f: f[::-1, ::-1]),
    "SE": (1.0, -1.0, lambda f: f[::-1, :]),
}


def check_orientations(fields: dict, tol_abs: float):
    """fields: name -> steady eta array (ny, nx) on the square basin.
    Cardinal runs must map onto E's field, diagonal runs onto NE's, each
    under the ORIENTATIONS transform, within tol_abs. NE itself must be
    symmetric under transposition (its own mirror axis)."""
    msgs = []
    ok = True
    for group, ref in (("E", ("N", "W", "S")), ("NE", ("NW", "SW", "SE"))):
        base = fields[group]
        for name in ref:
            mapped = ORIENTATIONS[name][2](fields[name])
            dev = float(np.max(np.abs(mapped - base)))
            line = (f"orientation {name} vs {group}: max |delta eta| = "
                    f"{dev:.3e} m (allowed {tol_abs:.1e})")
            if dev > tol_abs:
                ok = False
                line += " FAIL"
            msgs.append(line)
    dev = float(np.max(np.abs(fields["NE"] - fields["NE"].T)))
    line = (f"orientation NE self-transpose symmetry: max |delta eta| = "
            f"{dev:.3e} m (allowed {tol_abs:.1e})")
    if dev > tol_abs:
        ok = False
        line += " FAIL"
    msgs.append(line)
    return ok, msgs
