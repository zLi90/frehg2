#!/usr/bin/env python3
"""check_feature_coverage.py — the plan §8.3 feature-interaction lockstep
check (Q7). Validates docs/developer-guide/feature-coverage.md:

1. Every table row's coverage cell is non-empty and free of gap markers
   (GAP / TODO / unexercised / authored-unexercised).
2. Every test-like token a cell cites exists — same universe as
   check_bc_matrix.py (gtests, expanded ctest names, validation cases).
3. Every row of the GPU section carries an explicit `experimental` or
   `supported` status (the §8.3 release rule for device rows).

A feature or interaction removed from a release must be removed from the
table too — an empty or missing cell fails the build, so coverage decay is
loud, not silent (the v1 wind lesson).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from check_bc_matrix import (ROOT, ctest_names, gtest_names,  # noqa: E402
                             id_exists, validation_records)

DOC = ROOT / "docs" / "developer-guide" / "feature-coverage.md"
GAP_MARKERS = ("GAP", "TODO", "unexercised", "n/a")
TOKEN = re.compile(
    r"((?:regression|scaling|unit|validate|mpi)\.[\w.\-]+"
    r"|\b[A-Z]\w+\.[A-Z]\w+"
    r"|validation/[\w\-]+"
    r"|\.github/workflows/[\w\-]+\.yml"
    r"|scripts/[\w\-]+\.(?:sh|py))")


def main() -> int:
    errors: list[str] = []
    universe = {
        "gtest": gtest_names(),
        "ctest": ctest_names(),
        "validation": validation_records(),
    }
    universe["ctest_prefixes"] = {n for n in universe["ctest"]
                                  if n.endswith(".")}

    text = DOC.read_text(encoding="utf-8")
    section = ""
    rows = 0
    for lineno, line in enumerate(text.splitlines(), start=1):
        if line.startswith("## "):
            section = line[3:].strip()
            continue
        if not line.startswith("|") or set(line.replace("|", "").strip()) <= {"-"}:
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 2 or cells[0] in ("Feature (config surface)",
                                          "Interaction", "Row"):
            continue
        rows += 1
        gpu_row = section.startswith("3.") and len(cells) >= 3
        label = cells[0]
        coverage = cells[2] if gpu_row else cells[1]
        if gpu_row:
            status = cells[1].lower()
            if status not in ("experimental", "supported"):
                errors.append(f"line {lineno}: GPU row '{label}' status "
                              f"'{cells[1]}' must be experimental or supported")
        if not coverage:
            errors.append(f"line {lineno}: '{label}' has an empty coverage cell")
            continue
        for marker in GAP_MARKERS:
            if re.search(rf"\b{re.escape(marker)}\b", coverage):
                errors.append(f"line {lineno}: '{label}' carries the gap "
                              f"marker '{marker}' — forbidden in a release")
        tokens = TOKEN.findall(coverage)
        file_like = [t for t in tokens if "/" in t and
                     not t.startswith("validation/")]
        for token in file_like:
            if not (ROOT / token).exists():
                errors.append(f"line {lineno}: '{label}' cites missing file "
                              f"'{token}'")
        test_tokens = [t for t in tokens if t not in file_like]
        if not test_tokens and not file_like and not gpu_row:
            # GPU rows are governed by their explicit status (the §8.3
            # device rule); evidence tokens are validated when present.
            errors.append(f"line {lineno}: '{label}' cites no recognizable "
                          "test, gate, workflow, or case token")
        for token in test_tokens:
            if not id_exists(token, universe):
                errors.append(f"line {lineno}: '{label}' cites unknown test "
                              f"'{token}'")

    if rows < 20:
        errors.append(f"only {rows} table rows parsed — the table structure "
                      "changed; update this checker with it")
    if errors:
        for message in errors:
            print(f"ERROR: {message}")
        print(f"check_feature_coverage: FAIL ({len(errors)} error(s))")
        return 1
    print(f"check_feature_coverage: OK ({rows} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
