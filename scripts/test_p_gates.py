#!/usr/bin/env python3
"""Negative tests for the p2/p3 gate logic in run_scaling.py (v2 plan §6.3).

Each p-gate assertion must provably be able to fail: these tests feed
gate_p2_solver_threads / gate_p3_hybrid synthetic result sets containing
exactly one defect each and assert the gate returns nonzero, plus one clean
set per gate asserting it passes. final_volume is monkeypatched so no HDF5
files are needed.

Run: python3 scripts/test_p_gates.py   (also a ctest unit entry)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_scaling  # noqa: E402

FAILURES = []


def check(name: str, condition: bool) -> None:
    print(f"  {'ok  ' if condition else 'FAIL'} {name}")
    if not condition:
        FAILURES.append(name)


def result(ranks: int, threads: int, iters: float, solve_s: float,
           sim_s: float, volume: float) -> dict:
    return {
        "ranks": ranks, "threads": threads, "wall_s": sim_s,
        "timers": {"simulation": {"min": sim_s, "mean": sim_s,
                                  "max": sim_s, "count": 1}},
        "solver": {"gw": {"solves": 30, "iters_mean": iters, "iters_max": iters,
                          "rebuilds": 2, "retries": 0, "setup_s": 0.1,
                          "solve_s": solve_s}},
        "output": f"vol:{volume!r}",  # consumed by the patched final_volume
        "record": "unused",
    }


def patched_final_volume(output: Path) -> float:
    return float(str(output).split("vol:", 1)[1])


run_scaling.final_volume = patched_final_volume


def p2(results):
    return run_scaling.gate_p2_solver_threads(results)


def p3(results):
    return run_scaling.gate_p3_hybrid(results, 4)


def main() -> int:
    print("p2 negative battery:")
    clean = [result(1, t, 15.0, 10.0 / t, 100.0 / t, 5.0) for t in (1, 2, 4)]
    check("clean sweep passes", p2(clean) == 0)

    slow = [result(1, 1, 15.0, 10.0, 100.0, 5.0),
            result(1, 2, 15.0, 9.0, 60.0, 5.0),
            result(1, 4, 15.0, 11.0, 40.0, 5.0)]  # 4-thread solve SLOWER
    check("non-threading solve fails", p2(slow) != 0)

    drift = [result(1, 1, 15.0, 10.0, 100.0, 5.0),
             result(1, 2, 15.0, 6.0, 60.0, 5.0),
             result(1, 4, 16.0, 4.0, 40.0, 5.0)]  # iters 15 -> 16 = 6.7 %
    check("iteration drift > 2 % fails", p2(drift) != 0)

    vol = [result(1, 1, 15.0, 10.0, 100.0, 5.0),
           result(1, 4, 15.0, 4.0, 40.0, 5.0 * (1 + 1e-6))]  # volume off
    check("volume mismatch fails", p2(vol) != 0)

    nobase = [result(1, 4, 15.0, 4.0, 40.0, 5.0)]
    check("missing 1-thread baseline is an error", p2(nobase) == 2)

    print("p3 negative battery:")
    clean = [result(4, 1, 15.0, 3.0, 30.0, 5.0),
             result(2, 2, 15.5, 3.5, 35.0, 5.0),
             result(1, 4, 14.8, 4.0, 42.0, 5.0)]
    check("clean placements pass", p3(clean) == 0)

    vol = [result(4, 1, 15.0, 3.0, 30.0, 5.0),
           result(2, 2, 15.0, 3.5, 35.0, 5.0 * (1 + 1e-6))]
    check("placement volume mismatch fails", p3(vol) != 0)

    iters = [result(4, 1, 15.0, 3.0, 30.0, 5.0),
             result(1, 4, 17.0, 4.0, 42.0, 5.0)]  # 13.3 % > 10 %
    check("placement iteration drift > 10 % fails", p3(iters) != 0)

    if FAILURES:
        print(f"\ntest_p_gates: {len(FAILURES)} FAILURE(S): {FAILURES}")
        return 1
    print("\ntest_p_gates: all negative tests behave")
    return 0


if __name__ == "__main__":
    sys.exit(main())
