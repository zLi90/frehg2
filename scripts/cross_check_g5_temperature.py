#!/usr/bin/env python3
"""cross_check_g5_temperature.py — the v2 plan §4.3 follow-on cross-check
(NOT a gate): rerun g5 with the full temperature module active at the
paper's constant 20 C and require the salinity/flow results to be
statistically identical to the Q4-mode run.

With T uniform at the reference temperature the thermal density term is
exactly -beta_T * 0.0, temperature transport moves a uniform field, and
the temperature instance touches none of the flow/salinity state — so
the comparison is made BITWISE (max |delta| = 0 expected on every
salinity/flow dataset), which is stronger than "statistically
identical". Any nonzero delta means the temperature module leaks into a
run that should not feel it (the accidental-activation regression this
check exists to catch, in the spirit of "b1 stays golden with modules
off").

Usage: cross_check_g5_temperature.py --frehg BIN --repo ROOT --work DIR
"""

from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys

import h5py
import numpy as np
import yaml


def run(frehg: pathlib.Path, config: pathlib.Path, tag: str) -> pathlib.Path:
    workdir = config.parent
    log = workdir / f"run-{tag}.log"
    with open(log, "w", encoding="utf-8") as handle:
        result = subprocess.run([str(frehg), str(config)], cwd=workdir,
                                stdout=handle, stderr=subprocess.STDOUT, check=False)
    if result.returncode != 0:
        print(f"error: {tag} run exited {result.returncode}; see {log}")
        sys.exit(1)
    out = workdir / "out" / "output.h5"
    dst = workdir / "out" / f"{tag}.h5"
    shutil.move(out, dst)
    return dst


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frehg", type=pathlib.Path, required=True)
    parser.add_argument("--repo", type=pathlib.Path, required=True)
    parser.add_argument("--work", type=pathlib.Path, required=True)
    args = parser.parse_args()

    case = args.work / "g5-geng2015"
    if case.exists():
        shutil.rmtree(case)
    shutil.copytree(args.repo / "benchmarks" / "g5-geng2015", case)

    baseline = run(args.frehg, case / "g5-geng2015.yaml", "q4-baseline")

    with open(case / "g5-geng2015.yaml", encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    doc["simulation"]["id"] = "g5-geng2015-temperature"
    doc["modules"]["temperature"] = True
    # The full T-dependent module: transport of the (uniform 20 C) field,
    # thermal conduction/retardation, and the thermal density term armed
    # with a nonzero beta_T whose (T - T0) factor is exactly zero.
    doc["temperature"] = {
        "scheme": {"advection": "upwind"},
        "thermal_conductivity": 2.0,
        "heat_capacity_water": 4.184e6,
        "heat_capacity_solid": 2.0548e6,
        "dispersivity": {"longitudinal": 0.0, "transverse": 0.0},
    }
    doc["groundwater"]["density_coupling"]["thermal_expansion"] = 2.0e-4
    doc["groundwater"]["density_coupling"]["reference_temperature"] = 20.0
    doc["initial_conditions"]["temperature"] = {"groundwater": {"constant": 20.0}}
    doc["output"]["variables"]["temperature"] = ["temperature"]
    variant = case / "g5-temperature.yaml"
    with open(variant, "w", encoding="utf-8") as handle:
        yaml.safe_dump(doc, handle, sort_keys=False)

    with_temp = run(args.frehg, variant, "q5-temperature")

    ok = True
    with h5py.File(baseline, "r") as a, h5py.File(with_temp, "r") as b:
        for group in ("/groundwater/water_content", "/groundwater/hydraulic_head",
                      "/transport/concentration"):
            worst = 0.0
            for t in a[group]:
                va = a[f"{group}/{t}"][:]
                vb = b[f"{group}/{t}"][:]
                worst = max(worst, float(np.nanmax(np.abs(va - vb))))
            line = f"  {group}: max |delta| over all outputs = {worst:.3e}"
            if worst != 0.0:
                ok = False
                line += " NONZERO"
            print(line)
        tmax = sorted(int(k) for k in b["/temperature/temperature"])[-1]
        temp = b[f"/temperature/temperature/{tmax}"][:]
        drift = float(np.nanmax(np.abs(temp[~np.isnan(temp)] - 20.0)))
        print(f"  temperature stays 20 C: max |T - 20| at t_end = {drift:.3e} K")
        if drift > 1.0e-9:
            ok = False
    print("g5 temperature cross-check (plan §4.3):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
