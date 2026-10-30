#!/usr/bin/env python3
"""test_s3_shared.py — self-checks of the Stage-3 shared harness helpers
(plan §3.1/§3.2/§3.3): h5_bitwise hash/compare round-trip and sensitivity,
and the FREHG_STAGE3_OVERRIDES application and module gating. Imported by
every scripts/test_<x>_gates.py battery until that family's own checks
land (gate-first: the shared helpers must be provably able to fail before
any gate relies on them).

Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import json
import os
import pathlib
import sys
import tempfile

import h5py
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tests" / "regression"))

import h5_bitwise  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


def run_shared_checks() -> list[str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = pathlib.Path(tmpdir)
        model = tmp / "out.h5"
        with h5py.File(model, "w") as f:
            f.create_dataset("surface/eta/100", data=np.arange(12.0))
            f.create_dataset("monitor/mass_audit", data=np.ones((3, 5)))
        ref = tmp / "ref.json"

        print("h5_bitwise round trip and sensitivity:")
        expect(h5_bitwise.main_from(
            ["hash", str(model), "--case", "t", "--out", str(ref)]) == 0,
            "hash records")
        expect(h5_bitwise.main_from(
            ["compare", str(model), str(ref), "--case", "t"]) == 0,
            "identical file passes")
        with h5py.File(model, "a") as f:
            f["surface/eta/100"][3] += 1.0e-15
        expect(h5_bitwise.main_from(
            ["compare", str(model), str(ref), "--case", "t"]) == 1,
            "a 1e-15 perturbation fails (bitwise, not approximate)")
        expect(h5_bitwise.main_from(
            ["compare", str(model), str(ref), "--case", "t",
             "--allow", "surface/eta/*"]) == 0,
            "the allow glob exempts the perturbed dataset")
        with h5py.File(model, "a") as f:
            f.create_dataset("monitor/new_audit", data=np.zeros(3))
        expect(h5_bitwise.main_from(
            ["compare", str(model), str(ref), "--case", "t",
             "--allow", "surface/eta/*"]) == 1,
            "an unallowed new dataset fails")
        doc = json.loads(ref.read_text())
        doc["toolchain"]["compiler_version"] = "other"
        ref.write_text(json.dumps(doc))
        expect(h5_bitwise.main_from(
            ["compare", str(model), str(ref), "--case", "t",
             "--require-toolchain"]) == 77,
            "a toolchain mismatch exits 77 (skip)")

    print("FREHG_STAGE3_OVERRIDES application and module gating:")
    os.environ["FREHG_STAGE3_OVERRIDES"] = \
        "groundwater.reallocation=disabled;transport.surface_dryout=store"
    # re-import with the env set (run_regression reads it at import)
    for mod in [m for m in list(sys.modules) if m == "run_regression"]:
        del sys.modules[mod]
    import run_regression  # noqa: E402
    doc_gw = {"modules": {"groundwater": True, "transport": False}}
    applied = run_regression.apply_stage3_overrides(doc_gw)
    expect(applied == ["groundwater.reallocation=disabled"],
           "groundwater override applies, transport override gated off")
    expect(doc_gw["groundwater"]["reallocation"] == "disabled",
           "the dotted path lands in the document")
    doc_sw = {"modules": {"surface_water": True}}
    expect(run_regression.apply_stage3_overrides(doc_sw) == [],
           "no module, no override")
    del os.environ["FREHG_STAGE3_OVERRIDES"]
    return failures


if __name__ == "__main__":
    run_shared_checks()
    if failures:
        print(f"{len(failures)} shared-harness check(s) FAILED")
        sys.exit(1)
    print("all shared-harness checks passed")
