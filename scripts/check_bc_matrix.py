#!/usr/bin/env python3
"""check_bc_matrix.py — the plan §8.2 BC-matrix lockstep check (Q7).

Validates tests/coverage/bc_matrix.csv against three contracts:

1. COMPLETENESS: every cell of the canonical BC kind x target x coupling
   enumeration (mirroring the ConfigSchema acceptance matrix) has at least
   one row, and every row names a canonical cell (typo guard).
2. CROSS-REFERENCE: every test id a `tested` or `schema-rejected` row cites
   exists — as a gtest in tests/unit/*.cpp, as a ctest name registered in
   tests/CMakeLists.txt (simple nested-foreach generators are expanded), or
   as an in-tree validation-case record (validation/<case>).
3. RELEASE RULE: `limitation` rows fail (plan §8.2: the category must be
   empty at the Q7 release); pass --allow-limitations only mid-phase.

Schema drift alarm: the `kind` and `target` vocabularies are re-read from
src/core/ConfigSchema.cpp; a kind or target the canonical table does not
know fails the check, so a new BC kind cannot land without extending the
matrix (the same lockstep contract as check_parameter_docs.py).

Composition rows: a cell may cite several ids joined by ';' when its
coverage is composed (e.g. a side-sweep gate for the side axis plus a
coupled unit test for the exchange term); the CSV's header comments
document each composition. The checker validates every cited id.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

SIDES = ("west", "east", "south", "north")
COUPLINGS = ("uncoupled", "coupled")

# --------------------------------------------------------------------------
# Canonical enumeration. Kind families split where the code path splits
# (BcValueConfig::Form -> GwBcCode): head-constant vs head-hydrostatic vs
# flux vs flux-gravity; scalar_value per scalar spec.
# --------------------------------------------------------------------------


def canonical_cells() -> set[tuple[str, str, str]]:
    cells: set[tuple[str, str, str]] = set()
    surface_kinds = ("eta", "discharge", "velocity", "outflow",
                     "scalar_value-salinity", "scalar_value-temperature")
    for kind in surface_kinds:
        for side in SIDES:
            for coupled in COUPLINGS:
                cells.add((kind, f"surface-{side}-edge", coupled))
    gw_side_kinds = ("head-constant", "head-hydrostatic", "flux",
                     "scalar_value-salinity", "scalar_value-temperature")
    for kind in gw_side_kinds:
        for side in SIDES:
            for coupled in COUPLINGS:
                cells.add((kind, f"groundwater_side-{side}", coupled))
    for coupled in COUPLINGS:
        cells.add(("head-constant", "groundwater_top", coupled))
        cells.add(("flux", "groundwater_top", coupled))
        cells.add(("scalar_value-temperature", "groundwater_top", coupled))
        cells.add(("head-constant", "groundwater_bottom", coupled))
        cells.add(("flux", "groundwater_bottom", coupled))
        cells.add(("flux-gravity", "groundwater_bottom", coupled))
        cells.add(("scalar_value-temperature", "groundwater_bottom", coupled))
        cells.add(("scalar_cauchy", "groundwater_top", coupled))
    # Whole-family schema rejections ride an `any` coupling axis.
    for kind, target in (
            ("scalar_value-salinity", "groundwater_top"),
            ("scalar_value-salinity", "groundwater_bottom"),
            ("scalar_cauchy", "groundwater_side"),
            ("scalar_cauchy", "groundwater_bottom"),
            ("scalar_cauchy", "surface"),
            ("eta", "groundwater"), ("discharge", "groundwater"),
            ("velocity", "groundwater"), ("outflow", "groundwater"),
            ("head-constant", "surface"), ("flux", "surface"),
            ("flux-gravity", "off-bottom"), ("head-hydrostatic", "off-side")):
        cells.add((kind, target, "any"))
    return cells


# --------------------------------------------------------------------------
# The test-name universe.
# --------------------------------------------------------------------------


def gtest_names() -> set[str]:
    names = set()
    pattern = re.compile(r"TEST(?:_F)?\(\s*(\w+)\s*,\s*(\w+)\s*\)")
    for path in sorted((ROOT / "tests" / "unit").glob("*.cpp")):
        for suite, test in pattern.findall(path.read_text(encoding="utf-8")):
            names.add(f"{suite}.{test}")
    for path in sorted((ROOT / "tests" / "mpi").glob("*.cpp")):
        for suite, test in pattern.findall(path.read_text(encoding="utf-8")):
            names.add(f"{suite}.{test}")
    return names


def ctest_names() -> set[str]:
    """add_test(NAME ...) names from tests/CMakeLists.txt, expanding the
    file's own simple literal-list foreach loops (nested included)."""
    text = (ROOT / "tests" / "CMakeLists.txt").read_text(encoding="utf-8")
    loops: list[tuple[str, list[str]]] = []
    for var, items in re.findall(r"foreach\((\w+)((?:\s+[^)\s]+)+)\)", text):
        values = items.split()
        loops.append((var, values))
    # Compound list entries like "rank-invariance;b1" bind via list(GET);
    # expand them into their parts for substitution purposes.
    names = set()
    for match in re.findall(r"add_test\(NAME\s+([^\s)]+)", text):
        pending = {match}
        for var, values in loops:
            token = "${" + var + "}"
            expanded = set()
            for name in pending:
                if token in name:
                    for value in values:
                        for part in value.split(";"):
                            expanded.add(name.replace(token, part))
                else:
                    expanded.add(name)
            pending = expanded
        unresolved = {n for n in pending if "${" in n}
        names |= pending - unresolved
        for name in unresolved:
            # Variables bound by list(GET ...) (the rank-invariance pairs):
            # accept the prefix up to the first unresolved variable.
            names.add(name.split("${", 1)[0])
    return names


def validation_records() -> set[str]:
    return {f"validation/{p.name}" for p in (ROOT / "validation").iterdir()
            if p.is_dir()}


def id_exists(test_id: str, universe: dict[str, set[str]]) -> bool:
    if test_id in universe["gtest"] or test_id in universe["validation"]:
        return True
    if test_id in universe["ctest"]:
        return True
    # Generated ctest families registered via a prefix (see ctest_names).
    return any(test_id.startswith(prefix) and prefix
               for prefix in universe["ctest_prefixes"])


# --------------------------------------------------------------------------
# Schema drift alarm.
# --------------------------------------------------------------------------


def schema_vocabulary() -> tuple[set[str], set[str]]:
    text = (ROOT / "src" / "core" / "ConfigSchema.cpp").read_text(encoding="utf-8")
    kind_match = re.search(
        r'\{"kind",\s*enumeration\(true,\s*\{([^}]*)\}', text, re.S)
    target_match = re.search(
        r'\{"target",\s*enumeration\(true,\s*\{([^}]*)\}', text, re.S)
    if not kind_match or not target_match:
        raise SystemExit("check_bc_matrix: cannot locate the kind/target "
                         "vocabularies in ConfigSchema.cpp — update the "
                         "extraction regexes alongside the schema")
    kinds = set(re.findall(r'"(\w+)"', kind_match.group(1)))
    targets = set(re.findall(r'"(\w+)"', target_match.group(1)))
    return kinds, targets


CANONICAL_KIND_FAMILIES = {
    "eta": "eta", "discharge": "discharge", "velocity": "velocity",
    "outflow": "outflow", "head": "head-constant|head-hydrostatic",
    "flux": "flux|flux-gravity",
    "scalar_value": "scalar_value-salinity|scalar_value-temperature",
    "scalar_cauchy": "scalar_cauchy",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-limitations", action="store_true",
                        help="mid-phase mode: limitation rows warn instead "
                             "of failing (forbidden at the Q7 release)")
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []

    kinds, targets = schema_vocabulary()
    for kind in kinds:
        if kind not in CANONICAL_KIND_FAMILIES:
            errors.append(f"schema kind '{kind}' is unknown to the canonical "
                          "matrix — extend canonical_cells() and the CSV")
    expected_targets = {"surface", "groundwater_top", "groundwater_bottom",
                        "groundwater_side"}
    for target in targets:
        if target not in expected_targets:
            errors.append(f"schema target '{target}' is unknown to the "
                          "canonical matrix — extend canonical_cells()")

    universe = {
        "gtest": gtest_names(),
        "ctest": ctest_names(),
        "validation": validation_records(),
    }
    universe["ctest_prefixes"] = {n for n in universe["ctest"]
                                  if n.endswith(".")}

    csv_path = ROOT / "tests" / "coverage" / "bc_matrix.csv"
    rows = []
    with open(csv_path, encoding="utf-8") as handle:
        reader = csv.reader(line for line in handle
                            if line.strip() and not line.startswith("#"))
        header = next(reader)
        if header != ["kind", "target", "coupled", "disposition",
                      "test_or_record"]:
            errors.append(f"unexpected CSV header: {header}")
        rows = [row for row in reader if row]

    cells = canonical_cells()
    seen: set[tuple[str, str, str]] = set()
    for lineno, row in enumerate(rows, start=2):
        if len(row) != 5:
            errors.append(f"row {lineno}: expected 5 columns, got {len(row)}")
            continue
        kind, target, coupled, disposition, tests = (col.strip() for col in row)
        cell = (kind, target, coupled)
        if cell not in cells:
            errors.append(f"row {lineno}: ({kind}, {target}, {coupled}) is "
                          "not a canonical cell (typo, or extend "
                          "canonical_cells())")
            continue
        seen.add(cell)
        if disposition == "limitation":
            message = (f"row {lineno}: LIMITATION cell ({kind}, {target}, "
                       f"{coupled}) — {tests}")
            (warnings if args.allow_limitations else errors).append(message)
            continue
        if disposition not in ("tested", "schema-rejected"):
            errors.append(f"row {lineno}: unknown disposition '{disposition}'")
            continue
        for test_id in (t.strip() for t in tests.split(";")):
            if not test_id:
                errors.append(f"row {lineno}: empty test id")
            elif not id_exists(test_id, universe):
                errors.append(f"row {lineno}: test id '{test_id}' not found "
                              "in tests/unit, tests/mpi, tests/CMakeLists.txt "
                              "or validation/")

    for cell in sorted(cells - seen):
        errors.append(f"canonical cell has no row: {cell}")

    for message in warnings:
        print(f"WARNING: {message}")
    if errors:
        for message in errors:
            print(f"ERROR: {message}")
        print(f"check_bc_matrix: FAIL ({len(errors)} error(s))")
        return 1
    print(f"check_bc_matrix: OK ({len(rows)} rows, {len(cells)} canonical "
          f"cells, {len(universe['gtest'])} gtests, "
          f"{len(universe['ctest'])} ctest names)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
