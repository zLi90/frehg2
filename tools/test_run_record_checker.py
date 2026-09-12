#!/usr/bin/env python3
"""Negative tests for tools/check_run_record.py (v2 plan §6.3: a gate that
cannot fail cannot pass). Builds synthetic records and asserts the checker
accepts the valid one and rejects each corruption.

Usage: test_run_record_checker.py --checker tools/check_run_record.py
"""

from __future__ import annotations

import argparse
import copy
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

VALID = {
    "provenance": {
        "frehg_version": "1.0.0", "git_sha": "deadbeef", "build_type": "Release",
        "hostname": "test", "mpi_ranks": 4, "decomposition": [2, 2],
        "omp_threads": 1, "kokkos_backend": "OpenMP",
        "start_time": "2026-09-10T00:00:00Z", "end_time": "2026-09-10T00:01:00Z",
        "wall_seconds": 60.0, "input_file": "case.yaml",
        "input_sha256": "0" * 64, "finished": True,
    },
    "configuration": {"simulation": {"id": "t"}},
    "modules": {"surface_water": True, "groundwater": False, "transport": False,
                "coupling_mode": "none"},
    "boundary_conditions": [
        {"name": "tide", "target": "surface", "kind": "eta",
         "value": {"form": "constant", "constant": 1.0}, "global_cells": 42},
    ],
    "timers": {"simulation": {"count": 1, "min_s": 1.0, "mean_s": 1.0, "max_s": 1.0}},
    "solver": {},
    "closure": {},
}


def run_checker(checker: Path, record: dict, extra: list[str]) -> int:
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "run-record.yaml"
        path.write_text(yaml.safe_dump(record))
        r = subprocess.run([sys.executable, str(checker), str(path)] + extra,
                           capture_output=True, text=True)
        return r.returncode


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checker", type=Path, required=True)
    args = parser.parse_args()

    failures = []

    def expect(name: str, record: dict, extra: list[str], want_ok: bool) -> None:
        code = run_checker(args.checker, record, extra)
        ok = (code == 0) == want_ok
        print(f"  {name}: {'ok' if ok else 'UNEXPECTED (exit ' + str(code) + ')'}")
        if not ok:
            failures.append(name)

    expect("valid record accepted", VALID, ["--ranks", "4", "--threads", "1"], True)

    # One corruption per structural rule.
    for key in ("provenance", "configuration", "modules", "boundary_conditions",
                "timers", "solver", "closure"):
        broken = copy.deepcopy(VALID)
        del broken[key]
        expect(f"missing '{key}' rejected", broken, [], False)

    broken = copy.deepcopy(VALID)
    del broken["provenance"]["input_sha256"]
    expect("missing provenance.input_sha256 rejected", broken, [], False)

    broken = copy.deepcopy(VALID)
    del broken["timers"]["simulation"]
    broken["timers"]["other"] = {"count": 1, "min_s": 0.0, "mean_s": 0.0, "max_s": 0.0}
    expect("missing 'simulation' timer rejected", broken, [], False)

    broken = copy.deepcopy(VALID)
    del broken["timers"]["simulation"]["mean_s"]
    expect("timer missing a field rejected", broken, [], False)

    broken = copy.deepcopy(VALID)
    del broken["boundary_conditions"][0]["global_cells"]
    expect("BC entry missing a field rejected", broken, [], False)

    expect("rank mismatch rejected", VALID, ["--ranks", "2"], False)
    expect("thread mismatch rejected", VALID, ["--ranks", "4", "--threads", "8"], False)

    broken = copy.deepcopy(VALID)
    broken["provenance"]["decomposition"] = [3, 2]
    expect("decomposition/ranks mismatch rejected", broken, ["--ranks", "4"], False)

    broken = copy.deepcopy(VALID)
    broken["provenance"]["finished"] = False
    expect("unfinished record rejected (default)", broken, [], False)
    expect("unfinished record accepted with --allow-partial", broken,
           ["--allow-partial"], True)

    if failures:
        print("run-record checker negative tests: FAIL:", ", ".join(failures))
        return 1
    print("run-record checker negative tests: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
