"""gate_sc.py — pure check functions for the Stage-3 T12 saturated-contact
classification gates (plan §12D, sc1–sc2). No HDF5 I/O;
scripts/test_sc_gates.py feeds corrupted synthetic inputs to every checker
(the §6.3 negative battery).

A field set is a dict mapping "group/variable" to the concatenation of that
field's snapshots over the run's output times (the b5 rank lane's metric:
max |a - b| over the field's own max |a|).
"""

from __future__ import annotations

import numpy as np


def field_rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    """max |a - b| / max |a| over the finite entries of a (the b5
    rank-invariance metric, run_regression.gate_rank_invariance_b5)."""
    finite = np.isfinite(a)
    if not finite.any():
        return 0.0
    scale = max(float(np.abs(a[finite]).max()), 1.0e-12)
    return float(np.abs(a[finite] - b[finite]).max()) / scale


def worst_field(base: dict[str, np.ndarray],
                other: dict[str, np.ndarray]) -> tuple[float, str]:
    worst, name = 0.0, ""
    for key, a in base.items():
        d = field_rel_diff(a, other[key])
        if d > worst or not name:
            worst, name = d, key
    return worst, name


def check_rank_invariance(base: dict[str, np.ndarray], other: dict[str, np.ndarray],
                          bound: float, ranks: int) -> tuple[bool, str]:
    """sc1: under the head rule the one-step fields at n ranks match n = 1
    to rounding (the T11 legacy sweep: 3.7e-7 at n = 2 on x86, 2.4e-5 at
    n = 5 on arm64 — a seam-column restore chain)."""
    worst, name = worst_field(base, other)
    return worst <= bound, (f"n=1 vs n={ranks}: max rel diff {worst:.3e} ({name}; "
                            f"<= {bound:.0e})")


def check_legacy_exercises_seam(legacy_worst: dict[int, float],
                                floor: float) -> tuple[bool, str]:
    """sc1 anti-vacuity: the same scenario under the legacy in-place
    classification departs from n = 1 by more than `floor` at some rank
    count — the gate's scenario reaches the defect it guards."""
    if not legacy_worst:
        return False, "legacy lane: no rank counts measured"
    n, worst = max(legacy_worst.items(), key=lambda item: item[1])
    return worst > floor, (f"legacy lane: worst rel diff {worst:.3e} at n={n} "
                           f"(> {floor:.0e} required)")


def check_restore_active(realloc_m3: float, floor: float) -> tuple[bool, str]:
    """sc1 anti-vacuity: the saturated-contact restore fires in the step
    (a scenario without fringe cells passes vacuously)."""
    return (np.isfinite(realloc_m3) and abs(realloc_m3) > floor,
            f"restore volume {realloc_m3:.3e} m^3 (|.| > {floor:.0e} required)")


def check_bitwise(base: dict[str, np.ndarray], other: dict[str, np.ndarray],
                  label: str) -> tuple[bool, str]:
    """sc2: the threaded twin reproduces every field bitwise (the head
    classification reads only the pre-sweep snapshot, so no thread can
    observe another's writes)."""
    differing = [key for key, a in base.items()
                 if a.shape != other[key].shape
                 or not np.array_equal(a, other[key], equal_nan=True)]
    if differing:
        return False, f"{label}: fields differ: {', '.join(differing)}"
    return True, f"{label}: {len(base)} fields bitwise"
