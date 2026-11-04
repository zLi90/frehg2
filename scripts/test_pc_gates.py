#!/usr/bin/env python3
"""test_pc_gates.py — the §6.3 self-check battery for the Stage-3 T11
predictor–corrector consistency gates (pc1–pc8; plan §12C). Each checker
accepts exact synthetic input and fails on a corruption modelled on the
defect it exists to catch: the M1 one-sided-Ks churn (rl2's measured
2.91·T), an arithmetic-mean series flux, a leaking Ks = 0 wall, a stalled
wetting front, a drifting salt ledger, and the S1 along-slope creep."""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_pc as pc  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


def audit_rows(boundary_total: float, realloc: float, dropped: float,
               rows: int = 9) -> np.ndarray:
    """A synthetic gw_mass_audit with linear cumulative columns."""
    t = np.linspace(0.0, 1.0, rows)
    a = np.zeros((rows, 7))
    a[:, pc.TIME] = t * 3600.0
    a[:, pc.BOUNDARY_IN] = boundary_total * t
    a[:, pc.REALLOC] = realloc * t
    a[:, pc.REALLOC_DROPPED] = dropped * t
    a[:, pc.VOLUME] = 1.0 + a[:, pc.BOUNDARY_IN] + a[:, pc.REALLOC]
    return a


# --- pc1/pc2: the churn cap -------------------------------------------------
clean = audit_rows(0.5, realloc=1.0e-4, dropped=1.0e-4)   # churn/T = 6e-4
ok, msg = pc.check_churn_cap(clean, 0.01)
expect(ok, "churn cap accepts the consistent-class churn")
# The measured T10 red: rl2 throughput 0.497, churn 2.91 * T.
red = audit_rows(0.497, realloc=0.040, dropped=0.704)     # churn = 1.448 = 2.91 T
ok, msg = pc.check_churn_cap(red, 0.01)
expect(not ok, "churn cap REJECTS the rl2 one-sided-Ks churn (2.91 T)")

# --- pc1 anti-vacuity -------------------------------------------------------
ks = np.full((4, 4, 6), 1.0e-5)
ks[:2] = 1.0e-3                       # a 100x contrast plane
wc = np.full((4, 4, 6), 0.41)
wcs = np.full((4, 4, 6), 0.41)
n = pc.count_sat_contrast_faces(ks, wc, wcs)
expect(n > 0, "contrast counter sees the saturated Ks jump")
ok, _ = pc.check_hetero_activity(n, 4, churn_gross=1.0e-6)
expect(ok, "anti-vacuity accepts a live heterogeneous lane")
ok, _ = pc.check_hetero_activity(n, 4, churn_gross=0.0)
expect(not ok, "anti-vacuity REJECTS a lane whose restore never fires")
dry = np.full((4, 4, 6), 0.10)
ok, _ = pc.check_hetero_activity(pc.count_sat_contrast_faces(ks, dry, wcs), 4, 1.0e-6)
expect(not ok, "anti-vacuity REJECTS an unsaturated (vacuous) lane")

# --- pc3: the series flux ---------------------------------------------------
dz = np.full(8, 0.25)
layer_ks = np.where(np.arange(8) < 4, 1.0e-3, 1.0e-5)
exact = 3.0 / pc.series_resistance(dz, layer_ks)          # per unit area
ok, msg = pc.check_series_flux(exact, 3.0, dz, layer_ks, 1.0, 1.0e-10, "z")
expect(ok, "series check accepts the exact harmonic flux")
# The arithmetic-mean bias: the interface face carries (Ka+Kb)/2 over the
# intercell span instead of the two half-cells in series.
r_arith = pc.series_resistance(dz, layer_ks) \
    - 0.125 / 1.0e-3 - 0.125 / 1.0e-5 + 0.25 / (0.5 * (1.0e-3 + 1.0e-5))
q_arith = 3.0 / r_arith
ok, msg = pc.check_series_flux(q_arith, 3.0, dz, layer_ks, 1.0, 1.0e-10, "z")
expect(not ok, f"series check REJECTS the arithmetic flux (bias {q_arith/exact:.3f}x)")
ok, _ = pc.check_steady([exact, exact * (1 + 1e-12)], 1.0e-8)
expect(ok, "steady check accepts a settled run")
ok, _ = pc.check_steady([exact, exact * 1.01], 1.0e-8)
expect(not ok, "steady check REJECTS a still-moving run")

# --- pc3 per-face form --------------------------------------------------------
ok, _ = pc.check_face_resistance(0.25 / 1.980198e-5 * 1.0e-6, 0.0, 1.0e-6,
                                 0.25, 1.980198e-5, 1.0e-8, "iface")
expect(ok, "face-resistance check accepts the exact harmonic face")
ok, _ = pc.check_face_resistance(0.25 / 5.05e-4 * 1.0e-6, 0.0, 1.0e-6,
                                 0.25, 1.980198e-5, 1.0e-8, "iface")
expect(not ok, "face-resistance check REJECTS the arithmetic face (25x)")
# the z gravity term: a fresh face at rest (q -> 0 is guarded by |q| floor;
# a flowing face must include +delta)
ok, _ = pc.check_face_resistance(0.25 * (1.0 / 1.0e-3 * 5.0e-6 - 1.0), 0.25,
                                 5.0e-6 * 0.0 + 5.0e-6 * 0.0 + 5.0e-6, 0.25,
                                 1.0e-3, 1.0e-6, "zface")
expect(ok, "face-resistance check carries the z gravity term")

# --- pc4: the wall ----------------------------------------------------------
ok, _ = pc.check_wall_flux(0.0, leak_arith=3.1e-5)
expect(ok, "wall check accepts the sealed harmonic wall")
ok, _ = pc.check_wall_flux(3.1e-5, leak_arith=3.1e-5)
expect(not ok, "wall check REJECTS the arithmetic half-K leak")

# --- pc5: the front ---------------------------------------------------------
profile = np.array([0.41, 0.41, 0.41, 0.41, 0.31, 0.12, 0.08, 0.08])
ok, _ = pc.check_front_advance(profile, 4, 0.20)
expect(ok, "anti-stall accepts a front inside the low-K layer")
stalled = profile.copy()
stalled[4:] = 0.08
ok, _ = pc.check_front_advance(stalled, 4, 0.20)
expect(not ok, "anti-stall REJECTS a stalled front")

# --- pc6: the column restore ------------------------------------------------
ok, _ = pc.check_column_restore(1.0e-14, 0.0, 1.0e-2, 1.0e-10, "pc6a")
expect(ok, "column restore accepts the consistent class")
ok, _ = pc.check_column_restore(2.0e-4, 1.0e-4, 1.0e-2, 1.0e-10, "pc6a")
expect(not ok, "column restore REJECTS the M2/M3 mismatch class")

# --- pc7: the ledgers -------------------------------------------------------
theta = np.full((8, 1, 8), 0.41)
vol = np.full((8, 1, 8), 0.25 * 1.0 * 0.25)
conc = np.zeros((8, 1, 8))
conc[:4] = 35.0
tv0, sv0, mv0 = pc.ledgers(theta, conc, vol, 7.44e-4)
series = np.array([sv0, sv0 * (1 + 2e-13), sv0 * (1 - 1e-13)])
ok, _ = pc.check_ledger_constant(series, 1.0e-12, "salt")
expect(ok, "ledger check accepts a conserved salt total")
ok, _ = pc.check_ledger_constant(np.array([sv0, sv0 * 1.001]), 1.0e-12, "salt")
expect(not ok, "ledger check REJECTS a 1e-3 salt drift")
expect(abs(mv0 - (tv0 + 7.44e-4 * sv0)) < 1.0e-12 * mv0,
       "fluid-mass ledger composes from theta-V and s-theta-V")

# --- pc8: rest --------------------------------------------------------------
ok, _ = pc.check_rest(1.0e-18, 1.0e-4, 0.25, 1.0e-9, "salty slab")
expect(ok, "rest check accepts a resting slab")
# The derived S1 creep: K (rho/rho0 - 1) sin(a) * A at 35 psu, 10 % slope.
creep = 1.0e-4 * (7.44e-4 * 35.0) * 0.0995 * 0.25
ok, msg = pc.check_rest(creep, 1.0e-4, 0.25, 1.0e-9, "salty slab")
expect(not ok, f"rest check REJECTS the S1 creep ({msg.split('= ')[2].split(' ')[0]} of Ks*A)")

print()
if failures:
    print(f"test_pc_gates: {len(failures)} FAILURES")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)
print("test_pc_gates: all checks behave (accept exact, reject corrupted)")
