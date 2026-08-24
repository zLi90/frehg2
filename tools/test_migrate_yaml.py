#!/usr/bin/env python3
"""Test for migrate_yaml_v1_to_v2.py (run by CTest as unit.migrate_yaml).

Checks, for each archival v1 draft (b1-sw, b2-gw):
  1. migration succeeds and applies the header rename table (spot checks);
  2. the migrated output passes `frehg --validate` when anchored next to the
     benchmark's input/ directory;
  3. every leaf key produced by migration is present with an equal value in
     the committed v2 benchmark configuration (the committed configs are the
     migrated files plus additions sourced from the legacy `input` files,
     documented in each benchmark README).
"""

import argparse
import pathlib
import shutil
import subprocess
import sys
import tempfile

import yaml

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import migrate_yaml_v1_to_v2 as migrate_tool  # noqa: E402


FAILURES = []


def check(ok, what):
    if not ok:
        FAILURES.append(what)
        print(f"FAIL: {what}")


def leaf_items(node, prefix=""):
    """Yield (path, value) for every scalar leaf of a nested dict/list."""
    if isinstance(node, dict):
        for key, value in node.items():
            yield from leaf_items(value, f"{prefix}.{key}" if prefix else key)
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from leaf_items(value, f"{prefix}[{index}]")
    else:
        yield prefix, node


def lookup(node, path):
    """Resolve a leaf path produced by leaf_items; raises KeyError if absent."""
    current = node
    for part in path.replace("]", "").replace("[", ".").split("."):
        if isinstance(current, list):
            current = current[int(part)]
        else:
            current = current[part]
    return current


def run_case(repo, frehg, case_dir, v1_name, committed_name):
    case = repo / "benchmarks" / case_dir
    v1_path = case / v1_name
    with open(v1_path, "r", encoding="utf-8") as handle:
        v1 = yaml.safe_load(handle)

    v2 = migrate_tool.migrate(v1)
    notes = v2.get("_migration_notes", [])
    text = migrate_tool.render(dict(v2))
    v2.pop("_migration_notes", None)

    # 1. Rename-table spot checks.
    check("domain" in v2 and "grid" not in v2, f"{case_dir}: grid -> domain")
    check(v2["time"].get("t_end") == v1["time"]["Tend"], f"{case_dir}: Tend -> t_end")
    check(v2["time"].get("output_interval") == v1["time"]["dt_out"],
          f"{case_dir}: dt_out -> output_interval")
    check("max_step" not in v2.get("time", {}), f"{case_dir}: max_step dropped")
    check(any("max_step" in n for n in notes), f"{case_dir}: max_step drop is noted")
    check(v2["output"]["filename"].endswith("/output.h5"),
          f"{case_dir}: io.dir -> output.filename")
    if v1["modules"].get("solute") is not None:
        check(v2["modules"]["transport"] == bool(v1["modules"]["solute"]),
              f"{case_dir}: solute -> transport")

    # 2. Migrated output validates against the frozen schema.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = pathlib.Path(tmp)
        migrated = tmp_path / "migrated.yaml"
        migrated.write_text(text, encoding="utf-8")
        input_dir = case / "input"
        if input_dir.is_dir():
            shutil.copytree(input_dir, tmp_path / "input")
        result = subprocess.run([str(frehg), "--validate", str(migrated)],
                                capture_output=True, text=True)
        check(result.returncode == 0,
              f"{case_dir}: migrated output passes --validate\n{result.stdout}")

    # 3. The committed v2 config embeds the migration (values may be
    # augmented from the legacy `input` file where it is more specific than
    # the draft, so compare only keys whose provenance is the v1 draft; the
    # per-case exceptions are listed here with their legacy source).
    exceptions = {
        "b1-sw": {"simulation.title"},
        "b2-gw": {
            "simulation.title",
            # The v1 draft disagreed with the committed legacy input
            # (b2-gw/input); fidelity follows the legacy input:
            "groundwater.timestep.dt_init",   # legacy dt_min = 0.0001
            "groundwater.timestep.dt_min",    # legacy dt_min = 0.0001
            "groundwater.timestep.dt_max",    # legacy dt_max = 2.0
            "initial_conditions.groundwater.moisture.constant",  # legacy 0.033
            "domain.bottom_elevation.constant",  # land surface datum at 0
        },
    }
    with open(case / committed_name, "r", encoding="utf-8") as handle:
        committed = yaml.safe_load(handle)
    for path, value in leaf_items(v2):
        if path in exceptions.get(case_dir, set()):
            continue
        try:
            committed_value = lookup(committed, path)
        except (KeyError, IndexError, TypeError):
            check(False, f"{case_dir}: committed config is missing migrated key {path}")
            continue
        check(committed_value == value,
              f"{case_dir}: committed {path} = {committed_value!r} != migrated {value!r}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", required=True, type=pathlib.Path)
    parser.add_argument("--frehg", required=True, type=pathlib.Path)
    args = parser.parse_args()

    run_case(args.repo, args.frehg, "b1-sw", "input.v1.yaml", "b1-sw.yaml")
    run_case(args.repo, args.frehg, "b2-gw", "input.v1.yaml", "b2-gw.yaml")

    # Unknown v1 keys must be a hard error, never silently dropped.
    try:
        migrate_tool.migrate({"grid": {"nx": 1, "ny": 1, "nz": 1, "dx": 1, "dy": 1,
                                       "dz": 1, "mystery_key": 7}})
        check(False, "unknown v1 key must raise MigrationError")
    except migrate_tool.MigrationError:
        pass

    if FAILURES:
        print(f"test_migrate_yaml: FAILED ({len(FAILURES)})")
        return 1
    print("test_migrate_yaml: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
