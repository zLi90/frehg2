#!/usr/bin/env python3
"""Convert legacy Frehg ASCII golden outputs into the Frehg2 HDF5 layout.

The legacy model wrote one ASCII file per variable per output time
(``<var>_<t_seconds>``, one value per line, ``j*NX + i`` order for surface
fields). This tool packs a directory of those files into a single HDF5 file
whose datasets mirror the plan §7 layout (``/surface/<var>/<t>``), so the
regression harness compares goldens and model output with the same reader.

Variable name mapping (legacy file stem -> Frehg2 dataset):
  surf -> eta, depth -> depth, uu -> uu, vv -> vv, seepage -> seepage

Goldens are fetched from ``legacy/benchmarks/`` and are never committed
(plan §8.3).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import h5py
import numpy as np

LEGACY_TO_FREHG2 = {
    "surf": "eta",
    "depth": "depth",
    "uu": "uu",
    "vv": "vv",
    "seepage": "seepage",
}

LEGACY_TO_FREHG2_3D = {
    "head": "hydraulic_head",
    "moisture": "water_content",
    "qx": "qx",
    "qy": "qy",
    "qz": "qz",
}


def convert(golden_dir: Path, output_file: Path, variables: list[str]) -> int:
    """Pack ``<var>_<t>`` files into ``/surface/<var>/<t>`` datasets.

    Returns the number of datasets written.
    """
    stems = {v: k for k, v in LEGACY_TO_FREHG2.items()}
    pattern = re.compile(r"^([a-z]+)_(\d+)$")
    written = 0
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_file, "w") as out:
        for path in sorted(golden_dir.iterdir()):
            match = pattern.match(path.name)
            if not match:
                continue
            stem, t = match.group(1), int(match.group(2))
            var = LEGACY_TO_FREHG2.get(stem)
            if var is None or var not in variables:
                continue
            values = np.loadtxt(path, dtype=np.float64, ndmin=1)
            out.create_dataset(f"/surface/{var}/{t}", data=values)
            written += 1
        out.attrs["source"] = str(golden_dir)
    if written == 0:
        print(f"error: no golden files matched in {golden_dir}", file=sys.stderr)
        return 0
    _ = stems  # documented mapping; the inverse is not needed programmatically
    return written


def convert3d(golden_dir: Path, output_file: Path, variables: list[str]) -> int:
    """Pack 3D ``<var>_<t>`` files into ``/groundwater/<var>/<t>`` datasets.

    Legacy 3D ASCII files are already flat in the plan §7 order
    ``(j*NX + i)*NZ + k`` (reorder_subsurf); values transfer verbatim.
    Returns the number of datasets written.
    """
    pattern = re.compile(r"^([a-z]+)_(\d+)$")
    written = 0
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(output_file, "w") as out:
        for path in sorted(golden_dir.iterdir()):
            match = pattern.match(path.name)
            if not match:
                continue
            stem, t = match.group(1), int(match.group(2))
            var = LEGACY_TO_FREHG2_3D.get(stem)
            if var is None or var not in variables:
                continue
            values = np.loadtxt(path, dtype=np.float64, ndmin=1)
            out.create_dataset(f"/groundwater/{var}/{t}", data=values)
            written += 1
        out.attrs["source"] = str(golden_dir)
    if written == 0:
        print(f"error: no golden files matched in {golden_dir}", file=sys.stderr)
    return written


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("golden_dir", type=Path, help="directory of <var>_<t> ASCII files")
    parser.add_argument("output_file", type=Path, help="HDF5 file to create")
    parser.add_argument(
        "--variables",
        default="eta,depth,uu,vv",
        help="comma-separated Frehg2 variable names to convert",
    )
    args = parser.parse_args()
    written = convert(args.golden_dir, args.output_file, args.variables.split(","))
    if written == 0:
        return 1
    print(f"wrote {written} datasets to {args.output_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
