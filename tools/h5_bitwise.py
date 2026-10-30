#!/usr/bin/env python3
"""h5_bitwise.py — per-dataset SHA-256 hashing and bitwise comparison of
Frehg2 HDF5 outputs (Stage-3 plan §3.3, the default-bitwise gate).

  h5_bitwise.py hash OUTPUT.h5 [--out HASHES.json] [--case NAME]
  h5_bitwise.py compare OUTPUT.h5 HASHES.json --case NAME
      [--allow GLOB]... [--require-toolchain]

`hash` walks /surface /groundwater /transport /temperature /monitor and
records the SHA-256 of each dataset's raw bytes under the case name,
plus a toolchain fingerprint. With --out naming an existing file, the
case is merged into it. `compare` re-hashes the model file and requires
every recorded dataset to match and no recorded dataset to be missing;
a dataset matching an --allow glob (fnmatch on the dataset path) is
exempt from comparison in every direction (new, changed, or missing) —
the declared-exception channel of plan §2.6, for audit tables that gain
appended columns. A toolchain mismatch exits 77 (ctest SKIP_RETURN_CODE) because
bitwise reproduction is only meaningful on the recording toolchain.
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

GROUPS = ("surface", "groundwater", "transport", "temperature", "monitor")


def toolchain_fingerprint() -> dict:
    compiler = "/opt/homebrew/bin/g++-16"
    try:
        version = subprocess.run([compiler, "--version"], capture_output=True,
                                 text=True).stdout.splitlines()[0]
    except OSError:
        version = "unavailable"
    return {"compiler": compiler, "compiler_version": version,
            "machine": platform.machine(), "system": platform.system()}


def dataset_hashes(path: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    with h5py.File(path, "r") as handle:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset) and name.split("/")[0] in GROUPS:
                data = obj[()]
                hashes[name] = hashlib.sha256(
                    np.ascontiguousarray(data).tobytes()).hexdigest()
        handle.visititems(visit)
    if not hashes:
        raise SystemExit(f"h5_bitwise: no datasets under {GROUPS} in {path}")
    return hashes


def cmd_hash(args: argparse.Namespace) -> int:
    doc = {"toolchain": toolchain_fingerprint(), "cases": {}}
    if args.out and args.out.exists():
        doc = json.loads(args.out.read_text())
    doc["cases"][args.case] = dataset_hashes(args.model)
    text = json.dumps(doc, indent=1, sort_keys=True)
    if args.out:
        args.out.write_text(text)
        print(f"h5_bitwise: recorded {len(doc['cases'][args.case])} datasets "
              f"for '{args.case}' -> {args.out}")
    else:
        print(text)
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    doc = json.loads(args.reference.read_text())
    if args.require_toolchain and doc["toolchain"] != toolchain_fingerprint():
        print("h5_bitwise: toolchain differs from the recording toolchain; "
              "bitwise comparison skipped (exit 77)")
        print(f"  recorded: {doc['toolchain']}")
        print(f"  current:  {toolchain_fingerprint()}")
        return 77
    if args.case not in doc["cases"]:
        print(f"h5_bitwise: FAIL: case '{args.case}' not in {args.reference}")
        return 1
    recorded = doc["cases"][args.case]
    model = dataset_hashes(args.model)
    ok = True

    def allowed(name: str) -> bool:
        return any(fnmatch.fnmatch(name, glob) for glob in args.allow)

    for name, digest in sorted(recorded.items()):
        if name not in model:
            print(f"  {'missing (allowed)' if allowed(name) else 'MISSING'} {name}")
            ok = ok and allowed(name)
        elif model[name] != digest:
            print(f"  {'differs (allowed)' if allowed(name) else 'DIFFERS'} {name}")
            ok = ok and allowed(name)
    for name in sorted(set(model) - set(recorded)):
        print(f"  {'new (allowed)' if allowed(name) else 'NEW (not allowed)'} {name}")
        ok = ok and allowed(name)
    print(f"h5_bitwise [{args.case}]: {len(recorded)} recorded datasets:",
          "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main_from(argv: list[str] | None = None) -> int:
    """Entry point callable with an explicit argv (the regression harness
    imports this module instead of shelling out)."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p_hash = sub.add_parser("hash")
    p_hash.add_argument("model", type=Path)
    p_hash.add_argument("--out", type=Path, default=None)
    p_hash.add_argument("--case", required=True)
    p_cmp = sub.add_parser("compare")
    p_cmp.add_argument("model", type=Path)
    p_cmp.add_argument("reference", type=Path)
    p_cmp.add_argument("--case", required=True)
    p_cmp.add_argument("--allow", action="append", default=[])
    p_cmp.add_argument("--require-toolchain", action="store_true")
    args = parser.parse_args(argv)
    return cmd_hash(args) if args.cmd == "hash" else cmd_compare(args)


if __name__ == "__main__":
    sys.exit(main_from())
