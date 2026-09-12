#!/usr/bin/env python3
"""check_run_record.py — validate a Frehg2 run record (v2 plan §2A, gate r1).

Checks, in order:
  1. presence: the record file exists and parses as YAML;
  2. schema: every required top-level key and required subkey is present
     with the right type;
  3. launch match: ranks/threads/decomposition equal the values the caller
     observed at launch (--ranks/--threads when given);
  4. round-trip: the embedded resolved configuration, re-serialized to YAML,
     revalidates against the frehg schema (`frehg --validate`), and
     re-resolving the ORIGINAL input yields a resolved configuration equal
     to the embedded one (checked via `frehg --resolve`, which prints the
     resolved-config YAML for any valid input);
  5. finalized: for a completed run (the default), `provenance.finished`
     is true and `provenance.end_time` is set; --allow-partial accepts a
     truthful partial record instead.

Exit 0 when everything passes; nonzero with a message per failure.

Usage:
  check_run_record.py RECORD.yaml [--frehg BIN] [--input CONFIG.yaml]
      [--ranks N] [--threads N] [--allow-partial]

Round-trip checks (4) need both --frehg and --input; without them the
structural checks (1-3, 5) still run — the regression harness always
passes both.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

# Required top-level keys and, per key, required subkeys (v2 plan §2A.2).
REQUIRED = {
    "provenance": [
        "frehg_version", "git_sha", "build_type", "hostname",
        "mpi_ranks", "decomposition", "omp_threads", "kokkos_backend",
        "start_time", "input_file", "input_sha256", "finished",
    ],
    "configuration": [],       # the resolved config tree (validated by round-trip)
    "modules": ["surface_water", "groundwater", "transport", "coupling_mode"],
    "boundary_conditions": [],  # list (possibly empty)
    "timers": [],              # map path -> {count, min_s, mean_s, max_s}
    "solver": [],              # per-system telemetry (may be empty pre-solve)
    "closure": [],             # final audit rows (empty for t_end == t_start)
}


def fail(msg: str) -> None:
    print(f"check_run_record: FAIL: {msg}")
    sys.exit(1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    parser.add_argument("--frehg", type=Path, default=None)
    parser.add_argument("--input", type=Path, default=None,
                        help="the original input YAML the run was launched with")
    parser.add_argument("--ranks", type=int, default=None)
    parser.add_argument("--threads", type=int, default=None)
    parser.add_argument("--allow-partial", action="store_true")
    args = parser.parse_args()

    # 1. Presence.
    if not args.record.exists():
        fail(f"record file missing: {args.record}")
    try:
        record = yaml.safe_load(args.record.read_text())
    except yaml.YAMLError as e:
        fail(f"record does not parse as YAML: {e}")
    if not isinstance(record, dict):
        fail("record is not a YAML mapping")

    # 2. Schema.
    for key, subkeys in REQUIRED.items():
        if key not in record:
            fail(f"missing top-level key '{key}'")
        for sub in subkeys:
            if not isinstance(record[key], dict) or sub not in record[key]:
                fail(f"missing key '{key}.{sub}'")
    if not isinstance(record["boundary_conditions"], list):
        fail("'boundary_conditions' is not a list")
    if not isinstance(record["timers"], dict) or not record["timers"]:
        fail("'timers' is empty — no section was recorded")
    for path, entry in record["timers"].items():
        for field in ("count", "min_s", "mean_s", "max_s"):
            if not isinstance(entry, dict) or field not in entry:
                fail(f"timer '{path}' missing field '{field}'")
    if not args.allow_partial and "simulation" not in record["timers"]:
        fail("'timers' has no 'simulation' section (the time loop)")
    for bc in record["boundary_conditions"]:
        for field in ("name", "target", "kind", "value", "global_cells"):
            if field not in bc:
                fail(f"boundary condition entry missing field '{field}'")

    prov = record["provenance"]

    # 3. Launch match.
    if args.ranks is not None and int(prov["mpi_ranks"]) != args.ranks:
        fail(f"provenance.mpi_ranks = {prov['mpi_ranks']}, launched with {args.ranks}")
    if args.threads is not None and int(prov["omp_threads"]) != args.threads:
        fail(f"provenance.omp_threads = {prov['omp_threads']}, launched with {args.threads}")
    decomp = prov["decomposition"]
    if not (isinstance(decomp, list) and len(decomp) == 2):
        fail(f"provenance.decomposition is not [px, py]: {decomp}")
    if args.ranks is not None and int(decomp[0]) * int(decomp[1]) != args.ranks:
        fail(f"decomposition {decomp} does not multiply to {args.ranks} ranks")

    # 5 (before 4: cheap). Finalization.
    if not args.allow_partial:
        if not prov.get("finished", False):
            fail("provenance.finished is not true for a completed run")
        if not prov.get("end_time"):
            fail("provenance.end_time missing for a completed run")
        if "wall_seconds" not in prov:
            fail("provenance.wall_seconds missing for a completed run")

    # 4. Round-trip.
    if args.frehg is not None and args.input is not None:
        with tempfile.TemporaryDirectory() as tmp:
            embedded = Path(tmp) / "embedded.yaml"
            embedded.write_text(yaml.safe_dump(record["configuration"],
                                               sort_keys=False))
            # 4a. The embedded resolved configuration re-validates. Resolved
            # configs reference input files by path; those paths were
            # relative to the original config's directory, so validate from
            # there via a sibling temp file.
            side = args.input.parent / ".rr-roundtrip-embedded.yaml"
            side.write_text(embedded.read_text())
            try:
                r = subprocess.run([str(args.frehg), "--validate", str(side)],
                                   capture_output=True, text=True)
                if r.returncode != 0:
                    fail("embedded configuration does not revalidate:\n" +
                         r.stdout[-2000:])
                # 4b. Re-resolving the original input reproduces the embedded
                # resolved configuration exactly.
                r = subprocess.run([str(args.frehg), "--resolve", str(args.input)],
                                   capture_output=True, text=True)
                if r.returncode != 0:
                    fail("frehg --resolve failed on the original input:\n" +
                         r.stdout[-2000:])
                # The logger and PETSc share stdout; the YAML is fenced by
                # sentinels.
                begin = "--- FREHG RESOLVED CONFIG BEGIN ---"
                end = "--- FREHG RESOLVED CONFIG END ---"
                if begin not in r.stdout or end not in r.stdout:
                    fail("frehg --resolve output lacks the sentinel fence")
                yaml_text = r.stdout.split(begin, 1)[1].split(end, 1)[0]
                resolved = yaml.safe_load(yaml_text)
                if resolved != record["configuration"]:
                    import difflib
                    a = yaml.safe_dump(record["configuration"], sort_keys=True)
                    b = yaml.safe_dump(resolved, sort_keys=True)
                    diff = "\n".join(difflib.unified_diff(
                        a.splitlines(), b.splitlines(),
                        "record.configuration", "frehg --resolve", lineterm=""))
                    fail("embedded configuration != re-resolved input:\n" +
                         diff[:3000])
            finally:
                side.unlink(missing_ok=True)

    print(f"check_run_record: OK ({args.record})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
