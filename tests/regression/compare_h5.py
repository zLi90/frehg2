#!/usr/bin/env python3
"""Element-wise comparison of Frehg2 HDF5 outputs against a golden file.

Implements the plan §8.3 rule: a cell passes when
``|model - ref| <= max(abs_floor, rel * |ref|)`` for its variable, applied
per cell, per variable, per output time. A per-case tolerance file names the
variables, their floors, and the allowed exceedance fraction (0 for legacy
goldens unless stated). Achieved-vs-allowed errors are always printed so
tolerances can later be tightened from data.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml


def compare(model_file: Path, golden_file: Path, tolerances: dict) -> bool:
    """Compare /<group>/<var>/<t> datasets; returns True when all pass.

    The dataset group defaults to "surface"; groundwater cases set
    ``group: groundwater`` in their tolerance file.
    """
    ok = True
    group = str(tolerances.get("group", "surface"))
    allowed_fraction = float(tolerances.get("allowed_exceedance_fraction", 0.0))
    with h5py.File(model_file, "r") as model, h5py.File(golden_file, "r") as golden:
        for var, spec in tolerances["variables"].items():
            abs_floor = float(spec["abs_floor"])
            rel = float(spec["rel"])
            gvar = golden.get(f"/{group}/{var}")
            mvar = model.get(f"/{group}/{var}")
            if gvar is None:
                print(f"  {var}: no golden datasets, skipped")
                continue
            if mvar is None:
                print(f"  {var}: MISSING from model output")
                ok = False
                continue
            worst_ratio = 0.0
            worst_time = None
            for t in sorted(gvar.keys(), key=int):
                if t not in mvar:
                    print(f"  {var}@{t}: MISSING from model output")
                    ok = False
                    continue
                ref = gvar[t][:]
                mod = mvar[t][:]
                if ref.shape != mod.shape:
                    print(f"  {var}@{t}: shape {mod.shape} != golden {ref.shape}")
                    ok = False
                    continue
                allowed = np.maximum(abs_floor, rel * np.abs(ref))
                ratio = np.abs(mod - ref) / allowed
                exceed = float(np.mean(ratio > 1.0))
                if ratio.max() > worst_ratio:
                    worst_ratio = float(ratio.max())
                    worst_time = t
                if exceed > allowed_fraction:
                    print(
                        f"  {var}@{t}: FAIL max|d|={np.abs(mod - ref).max():.3e} "
                        f"(allowed max({abs_floor:.1e}, {rel:.0%}|ref|)), "
                        f"{exceed:.2%} of cells exceed"
                    )
                    ok = False
            print(
                f"  {var}: worst achieved/allowed = {worst_ratio:.3f}"
                + (f" at t={worst_time}" if worst_time else "")
            )
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("golden", type=Path)
    parser.add_argument("tolerances", type=Path)
    args = parser.parse_args()
    with open(args.tolerances, encoding="utf-8") as handle:
        tolerances = yaml.safe_load(handle)
    print(f"comparing {args.model} vs {args.golden}")
    return 0 if compare(args.model, args.golden, tolerances) else 1


if __name__ == "__main__":
    sys.exit(main())
