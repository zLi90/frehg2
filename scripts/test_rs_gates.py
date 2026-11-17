#!/usr/bin/env python3
"""test_rs_gates.py — the §6.3 self-check battery for the Stage-3 T13
side-face scatter gate (rs1; plan §12E). Each checker accepts exact
synthetic input and fails on a corruption modelled on the defect it exists
to catch: a lost corner update (the T12 binary's g4a at 6 threads: depth
3.8e-9 m, audit 4.5e-8 m^3), a lane whose partition keeps every corner pair
on one thread (4 and 8 threads on g4a), and a one-ulp thread-order
difference. The list replica and the partition model are checked against
the counts the C++ builder produces (SweRs.SideFaceListsRepeatCornerAndAxisCells)
and against hand-derived Kokkos partitions."""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import gate_rs as rs  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


def cell_counts(entries: list) -> list[tuple[tuple[int, int], int]]:
    counts: list[tuple[tuple[int, int], int]] = []
    for j, i, _ in entries:
        if counts and counts[-1][0] == (j, i):
            counts[-1] = ((j, i), counts[-1][1] + 1)
        else:
            counts.append(((j, i), 1))
    return counts


whole = [(-1.0, -1.0), (99.0, -1.0), (99.0, 99.0), (-1.0, 99.0)]

print("side-face list replica (same shapes as SweRs.SideFaceListsRepeatCornerAndAxisCells)")
box = rs.side_face_list(3, 2, 1.0, 1.0, whole)
expect(cell_counts(box) == [((0, 0), 2), ((0, 1), 1), ((0, 2), 2),
                            ((1, 0), 2), ((1, 1), 1), ((1, 2), 2)],
       "3 x 2: corners twice, edge cells once")
column = rs.side_face_list(1, 3, 1.0, 1.0, whole)
expect(cell_counts(column) == [((0, 0), 3), ((1, 0), 2), ((2, 0), 3)],
       "1 x 3: both x faces everywhere, a y face at each end")
expect([f for _, _, f in box[:2]] == ["x-", "y-"], "face order x-, x+, y-, y+")
expect(rs.side_face_list(3, 2, 1.0, 1.0, whole, periodic_x=True)
       == [(0, 0, "y-"), (0, 1, "y-"), (0, 2, "y-"), (1, 0, "y+"), (1, 1, "y+"), (1, 2, "y+")],
       "a periodic x axis attaches no x faces")
strip = [(4.0, -0.1), (5.1, -0.1), (5.1, 9.0), (4.0, 9.0)]
expect(cell_counts(rs.side_face_list(5, 3, 1.0, 1.0, strip))
       == [((0, 4), 2), ((1, 4), 1), ((2, 4), 2)], "an east strip: its two corners twice")
expect(rs.polygon_contains([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)], 1.0, 0.5),
       "a centre on the polygon edge counts as inside")

print("Kokkos static partition model")
expect(rs.static_chunk_starts(14, 4) == [4, 8, 12], "14 entries / 4 threads: 4 each")
expect(rs.static_chunk_starts(14, 6) == [3, 6, 9, 12], "14 entries / 6 threads: 3 each, 5 used")
expect(rs.static_chunk_starts(14, 1) == [], "one thread: no boundary")
expect(rs.static_chunk_starts(1000, 4) == [256, 512, 768],
       "1000 entries / 4 threads: chunk 8 (auto), 32 chunks per thread")

print("rs1 anti-vacuity (straddle)")
g4a = rs.side_face_list(4, 3, 1.0, 1.0, whole)
expect(len(g4a) == 14, "g4a: 14 entries")
for threads, want in ((3, True), (4, False), (5, True), (6, True), (7, False), (8, False)):
    expect(rs.check_straddle(g4a, threads)[0] is want,
           f"g4a at {threads} threads {'straddles' if want else 'keeps every pair on one thread'}")
expect(rs.straddling_pairs(g4a, 6) == [(8, 9)], "g4a at 6 threads: corner (2, 0), entries 8/9")

print("rs1 fields bitwise")
rng = np.random.default_rng(13)
serial = {f"surface/depth/{t}": 1.0e-2 * rng.random(12) for t in range(0, 50001, 5000)}
serial["monitor/point_a"] = rng.random((11, 3))
same = {key: value.copy() for key, value in serial.items()}
expect(rs.check_fields_bitwise(serial, same, "threads")[0], "identical datasets pass")
lost = {key: value.copy() for key, value in serial.items()}
lost["surface/depth/40000"][8] += 3.8e-9  # the T12 binary's 6-thread corner loss
ok, msg = rs.check_fields_bitwise(serial, lost, "threads")
expect(not ok and "depth/40000" in msg, "a lost corner update fails and is named")
ulp = {key: value.copy() for key, value in serial.items()}
ulp["monitor/point_a"][4, 1] = np.nextafter(ulp["monitor/point_a"][4, 1], np.inf)
expect(not rs.check_fields_bitwise(serial, ulp, "threads")[0], "a one-ulp difference fails")
gone = {key: value for key, value in serial.items() if key != "monitor/point_a"}
expect(not rs.check_fields_bitwise(serial, gone, "threads")[0], "a missing dataset fails")

print("rs1 audit")
audit = np.cumsum(rng.random((5001, 10)), axis=0)
expect(rs.check_audit_close(audit, audit + 4.4e-16, 1.0e-12, "threads")[0],
       "reduction-order drift passes")
lossy = audit.copy()
lossy[4000:, 1] -= 4.5e-8
expect(not rs.check_audit_close(audit, lossy, 1.0e-12, "threads")[0],
       "the T12 binary's audit departure fails")
expect(not rs.check_audit_close(audit, audit[:-1], 1.0e-12, "threads")[0],
       "a short audit table fails")
nan = audit.copy()
nan[3, 3] = np.nan
expect(not rs.check_audit_close(audit, nan, 1.0e-12, "threads")[0], "a NaN audit fails")

if failures:
    print(f"rs gate battery: FAIL ({len(failures)}: {', '.join(failures)})")
    sys.exit(1)
print("rs gate battery: PASS")
