"""gate_rs.py — pure check functions for the Stage-3 T13 side-face
right-hand-side scatter gates (plan §12E, rs1). No HDF5 I/O;
scripts/test_rs_gates.py feeds corrupted synthetic inputs to every checker
(the §6.3 negative battery). rs0 is compiled (tests/unit/test_swe.cpp,
SweRs.*).

A side-face list is the sequence of (j, i, face) entries BoundarySet's
appendSideFaces builds on one rank: owned cells in (j, i) order, inside the
polygon by their centre (edges and vertices count as inside), one entry per
non-periodic domain-edge face in the order x-, x+, y-, y+. A corner cell
holds two consecutive entries, a cell on an nx = 1 or ny = 1 axis two or
three.
"""

from __future__ import annotations

import numpy as np


def _on_segment(px: float, py: float, a: tuple, b: tuple) -> bool:
    (ax, ay), (bx, by) = a, b
    cross = (bx - ax) * (py - ay) - (by - ay) * (px - ax)
    if cross != 0.0:
        return False
    return min(ax, bx) <= px <= max(ax, bx) and min(ay, by) <= py <= max(ay, by)


def polygon_contains(polygon: list, x: float, y: float) -> bool:
    """Ray-cast containment with boundary points inside (src/bc/Polygon.hpp)."""
    n = len(polygon)
    for k in range(n):
        if _on_segment(x, y, polygon[k], polygon[(k + 1) % n]):
            return True
    inside = False
    for k in range(n):
        (xi, yi), (xj, yj) = polygon[k], polygon[(k + 1) % n]
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi) + xi:
            inside = not inside
    return inside


def side_face_list(nx: int, ny: int, dx: float, dy: float, polygon: list,
                   periodic_x: bool = False, periodic_y: bool = False) -> list:
    """Serial replica of appendSideFaces over the whole grid (one rank)."""
    entries = []
    for j in range(ny):
        for i in range(nx):
            if not polygon_contains(polygon, (i + 0.5) * dx, (j + 0.5) * dy):
                continue
            if i == 0 and not periodic_x:
                entries.append((j, i, "x-"))
            if i == nx - 1 and not periodic_x:
                entries.append((j, i, "x+"))
            if j == 0 and not periodic_y:
                entries.append((j, i, "y-"))
            if j == ny - 1 and not periodic_y:
                entries.append((j, i, "y+"))
    return entries


def static_chunk_starts(entries: int, threads: int) -> list[int]:
    """Entry indices at which a Kokkos OpenMP RangePolicy hands the launch to
    the next thread. For entries < 40 * threads the auto chunk size is 1
    (Kokkos_ExecPolicy.hpp set_auto_chunk_size) and the static partition
    gives every thread ceil(entries / threads) consecutive entries
    (HostThreadTeamData::set_work_partition)."""
    if entries <= 0 or threads <= 1:
        return []
    chunk = 1
    while chunk * 100 * threads < entries:
        chunk *= 2
    if chunk < 128:
        chunk = 1
        while chunk * 40 * threads < entries and chunk < 128:
            chunk *= 2
    chunks = -(-entries // chunk)
    part = -(-chunks // threads)
    return [k * part * chunk for k in range(1, threads) if k * part * chunk < entries]


def straddling_pairs(entries: list, threads: int) -> list[tuple[int, int]]:
    """(m - 1, m) for every pair of entries on one cell that the static
    partition hands to two different threads (the race's precondition)."""
    starts = set(static_chunk_starts(len(entries), threads))
    return [(m - 1, m) for m in range(1, len(entries))
            if m in starts and entries[m][:2] == entries[m - 1][:2]]


def check_straddle(entries: list, threads: int) -> tuple[bool, str]:
    """rs1 anti-vacuity: at the lane's thread count a duplicated cell's
    entries straddle a chunk boundary, so two threads can update one Srhs
    slot at once. (At 4 and 8 threads every shipped corner pair falls inside
    one chunk; that lane cannot lose an update.)"""
    pairs = straddling_pairs(entries, threads)
    cells = sorted({entries[a][:2] for a, _ in pairs})
    return bool(pairs), (f"{threads} threads, {len(entries)} entries: "
                         f"{len(pairs)} duplicated-cell pair(s) straddle a chunk "
                         f"boundary {cells}")


def check_fields_bitwise(serial: dict[str, np.ndarray], threaded: dict[str, np.ndarray],
                         label: str) -> tuple[bool, str]:
    """rs1: every surface field and point monitor of the threaded run equals
    the 1-thread run bitwise (the corner wall entries' second contribution is
    exactly 0, so the atomic order cannot round differently)."""
    missing = [key for key in serial if key not in threaded]
    differing = [key for key, a in serial.items() if key in threaded
                 and (a.shape != threaded[key].shape
                      or not np.array_equal(a, threaded[key], equal_nan=True))]
    if missing or differing:
        return False, (f"{label}: {len(differing)} dataset(s) differ "
                       f"{differing[:4]}, {len(missing)} missing {missing[:4]}")
    return True, f"{label}: {len(serial)} datasets bitwise"


def check_audit_close(serial: np.ndarray, threaded: np.ndarray, bound: float,
                      label: str) -> tuple[bool, str]:
    """rs1: the mass-audit table of the threaded run matches the 1-thread
    run to reduction order (its sums run in partition order). A lost
    Srhs update books the face volume but never moves it, so the audit's
    volume column departs from the booked terms."""
    if serial.shape != threaded.shape:
        return False, f"{label}: audit shape {threaded.shape} != {serial.shape}"
    worst = float(np.max(np.abs(serial - threaded))) if serial.size else 0.0
    return bool(np.isfinite(worst) and worst <= bound), (
        f"{label}: mass_audit max |d| {worst:.3e} (<= {bound:.0e})")
