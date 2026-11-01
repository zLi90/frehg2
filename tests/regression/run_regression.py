#!/usr/bin/env python3
"""Frehg2 regression driver (plan §8.3): runs a benchmark case in a scratch
directory and applies its plan §9 gate.

Subcommands:
  b1                element-wise fields vs the legacy ASCII goldens
                    (converted on the fly; goldens are never committed) plus
                    the cumulative mass-balance gate
  b2                element-wise subsurface fields vs the legacy ASCII
                    goldens plus the Warrick wetting-front gate (plan §9 b2)
  b3                digitized Kirkland contour gate + internal mass balance
                    (plan §9 b3; reference points parsed from the legacy
                    makeplot.py)
  b4                outlet-hydrograph metric gate vs the kinematic-wave
                    reference
  b5                coupled tilted-V envelope gate (plan §9 b5): outlet
                    discharge and ponding storage against the
                    PF/CATHY/HGS/Cast3M curves, selected by --scenario
                    (rain | norain) and --coupling (sync | subcycled)
  b1-restart        run -> checkpoint at half time -> restart; the restarted
                    run must match the uninterrupted one to 1e-12 relative
  b2-restart        the groundwater restart determinism gate (checkpoint
                    keyed to the crossed whole-second boundary, exact state
                    time in the header)
  b5-restart        coupled restart determinism (SWE + GW + coupler state:
                    the seepage accumulator and the adaptive clock)
  rank-invariance   fields at 1/2/4 ranks agree to 1e-12 relative in strict
                    mode (rank-invariant jacobi preconditioning, tightened
                    tolerances) and 1e-4 in the default bjacobi/icc mode
                    (plan §8.2, amendment A5). The b2 variant runs the
                    benchmark column replicated to 4x4 so 2- and 4-rank
                    decompositions exist (amendment A8). The b5 variant runs
                    the coupled case through t = 600 s (plan §8.2); its
                    strict one-step bound is 1e-9 (V2-A9: the coupled
                    threshold density admits one platform-rounding branch
                    flip inside even a single step).
  smoke-b5          sanitizer path-coverage run of a b5 scenario x coupling
                    combination on a shortened horizon (plan §10 P5,
                    amendment A21); clean instrumented exit is the check,
                    the plan §9 metrics are not applied
  smoke-b6          sanitizer path-coverage run of a b6 variant on a
                    shortened horizon (same purpose as smoke-b5)
  g4                v2 Q4 analytic drawdown + evaporative concentration
                    (plan §3.3, four sub-criteria (a)-(d); per-PR;
                    authored failing per §6.1 — (c)/(d) turn green with
                    the Q4 capability)
  g5                v2 Q4 Geng & Boufadel (2015) bare-soil salinization,
                    code-to-code vs the digitized MARUN figures
                    (nightly-class; V2-A13/V2-A14 criteria; authored
                    failing per §6.1)
  g9                v2 Q6 steady wind setup (a nonlinear TELEMAC replica /
                    b linear Garratt / c sloping-bottom quadrature;
                    per-PR, closed forms in-script)
  g10               v2 Q6 seiche relaxation vs the Merian period (per-PR)
  wind-orient       v2 Q6 8-orientation wind-setup symmetry battery
                    (the §8.1/§9 x-battery row; per-PR)

Every run happens in a copy of the benchmark case directory so relative
input paths keep working and the repository stays clean.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import math
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml

HERE = Path(__file__).resolve().parent

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "tools"))

import ascii_golden_to_h5  # noqa: E402
import compare_h5  # noqa: E402

# Strict rank-invariance mode: a rank-invariant preconditioner and
# tolerances tight enough that every solve lands at machine precision, so
# rank-count differences stay at rounding level (plan §8.2).
STRICT_PETSC_OPTIONS = ["-fs_pc_type", "jacobi", "-fs_ksp_rtol", "1e-13", "-fs_ksp_atol", "1e-16"]
STRICT_PETSC_OPTIONS_GW = ["-gw_pc_type", "jacobi", "-gw_ksp_rtol", "1e-13",
                           "-gw_ksp_atol", "1e-16"]

# g1 solver-invariance gate (v2 plan §2.3): when --solver is given, every
# staged configuration gets a solver block selecting that preconditioner for
# both systems, and the gate criteria are unchanged — the solver choice must
# not change the physics. Set from args in main(); staging applies it.
SOLVER_OVERRIDE: str | None = None

# p1 backend-invariance gate (v2 plan §2B.3): when --mat-type is given, the
# staged solver block additionally selects the PETSc matrix/vector backend
# ("aij" host, "aijkokkos" Kokkos Kernels). Same contract as SOLVER_OVERRIDE:
# gate criteria unchanged — the linear-algebra backend must not change the
# physics. This is the CPU rehearsal of the exact code path a GPU build takes.
MAT_TYPE_OVERRIDE: str | None = None

# Stage-3 override channel (Stage-3 plan §3.2): FREHG_STAGE3_OVERRIDES is a
# ';'-separated list of dotted.key=value pairs applied to every staged YAML
# whose modules make the key meaningful. check_run_record() then passes each
# as --expect, so an override that never reached a run fails that run's r1
# hook (anti-vacuity: the V2-A9 lesson).
STAGE3_OVERRIDES: list[tuple[str, str]] = [
    (item.partition("=")[0], item.partition("=")[2])
    for item in os.environ.get("FREHG_STAGE3_OVERRIDES", "").split(";") if item
]


def override_applies(doc: dict, dotted: str) -> bool:
    """A groundwater.* override needs modules.groundwater; transport.* needs
    modules.transport; anything else applies unconditionally."""
    root = dotted.split(".")[0]
    if root in ("groundwater", "transport"):
        return bool(doc.get("modules", {}).get(root, False))
    return True


def apply_stage3_overrides(doc: dict) -> list[str]:
    """Apply the active overrides to one config document; returns the dotted
    keys actually applied (for the --expect hook)."""
    applied = []
    for dotted, value in STAGE3_OVERRIDES:
        if not override_applies(doc, dotted):
            continue
        node = doc
        parts = dotted.split(".")
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        parsed: object = value
        if value in ("true", "false"):
            parsed = (value == "true")
        else:
            try:
                parsed = int(value)
            except ValueError:
                try:
                    parsed = float(value)
                except ValueError:
                    pass
        node[parts[-1]] = parsed
        applied.append(f"{dotted}={value}")
    return applied


def apply_solver_override(case_dir: Path) -> None:
    """Write the SOLVER_OVERRIDE / MAT_TYPE_OVERRIDE and any Stage-3
    overrides (FREHG_STAGE3_OVERRIDES) into every staged YAML."""
    if SOLVER_OVERRIDE is None and MAT_TYPE_OVERRIDE is None and not STAGE3_OVERRIDES:
        return
    for config in sorted(case_dir.glob("*.yaml")):
        with open(config, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        if SOLVER_OVERRIDE is not None or MAT_TYPE_OVERRIDE is not None:
            solver = doc.setdefault("solver", {})
            for system in ("surface", "groundwater"):
                entry = solver.setdefault(system, {})
                if SOLVER_OVERRIDE is not None:
                    entry["preconditioner"] = SOLVER_OVERRIDE
                if MAT_TYPE_OVERRIDE is not None:
                    entry["mat_type"] = MAT_TYPE_OVERRIDE
            print(f"  solver override: {config.name} -> "
                  f"preconditioner {SOLVER_OVERRIDE or '(default)'}"
                  f" mat_type {MAT_TYPE_OVERRIDE or '(default)'}")
        applied = apply_stage3_overrides(doc)
        for item in applied:
            print(f"  stage3 override: {config.name} -> {item}")
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)


def run_case(frehg: Path, config: Path, workdir: Path, mpiexec: Path, ranks: int,
             extra_args: list[str] | None = None) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(mpiexec), "-n", str(ranks), str(frehg), str(config)] + (extra_args or [])
    log = workdir / f"run_n{ranks}.log"
    with open(log, "w", encoding="utf-8") as handle:
        result = subprocess.run(cmd, cwd=workdir, stdout=handle, stderr=subprocess.STDOUT,
                                check=False)
    if result.returncode != 0:
        print(f"error: frehg exited {result.returncode}; log tail ({log}):")
        print("\n".join(log.read_text().splitlines()[-25:]))
        raise SystemExit(1)
    check_run_record(frehg, config, ranks)


def check_run_record(frehg: Path, config: Path, ranks: int) -> None:
    """r1 gate hook (v2 plan §2A.3): every successful gate run must leave a
    valid run record beside its HDF5 output — present, schema-valid, launch-
    matched, and configuration-round-tripping."""
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    output = Path(doc["output"]["filename"])
    record = (config.parent / output).parent / "run-record.yaml"
    checker = HERE.parent.parent / "tools" / "check_run_record.py"
    cmd = [sys.executable, str(checker), str(record),
           "--frehg", str(frehg), "--input", str(config), "--ranks", str(ranks)]
    threads = os.environ.get("OMP_NUM_THREADS")
    if threads:
        cmd += ["--threads", threads]
    for dotted, value in STAGE3_OVERRIDES:
        if override_applies(doc, dotted):
            cmd += ["--expect", f"{dotted}={value}"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout, end="")
    if result.returncode != 0:
        print(f"r1 run-record gate: FAIL for {config.name}")
        raise SystemExit(1)


def stage_case(repo: Path, case: str, workdir: Path) -> Path:
    """Copy a benchmark (or validation) case directory into the scratch area."""
    source = repo / "benchmarks" / case
    if not source.exists():
        source = repo / "validation" / case
    target = workdir / case
    if target.exists():
        shutil.rmtree(target)
    # ignore=: a case run in place (validation/ cases are runnable from the
    # source tree) leaves an out/ behind; staging must not carry stale results
    # into the gate.
    shutil.copytree(source, target, ignore=shutil.ignore_patterns("out"))
    apply_solver_override(target)
    return target


def mass_balance_error(output: Path) -> tuple[float, float]:
    """(closure error, total rain) in m^3 from the mass-audit table."""
    with h5py.File(output, "r") as handle:
        table = handle["/monitor/mass_audit"][:]
    volume, rain, evap, outflow, bc = (table[:, k] for k in range(1, 6))
    error = (volume[-1] - volume[0]) + outflow[-1] - rain[-1] + evap[-1] - bc[-1]
    return float(error), float(rain[-1])


def gate_b1(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "b1-sw", args.work)
    run_case(args.frehg, case_dir / "b1-sw.yaml", case_dir, args.mpiexec, args.ranks)

    golden_src = args.legacy / "b1-sw" / "out"
    golden_h5 = args.work / "b1-golden.h5"
    written = ascii_golden_to_h5.convert(golden_src, golden_h5, ["eta", "depth", "uu", "vv"])
    if written == 0:
        return 1
    with open(HERE / "tolerances" / "b1-sw.yaml", encoding="utf-8") as handle:
        tolerances = yaml.safe_load(handle)
    print(f"b1 gate: comparing against {written} golden datasets")
    ok = compare_h5.compare(case_dir / "out" / "output.h5", golden_h5, tolerances)

    error, rain = mass_balance_error(case_dir / "out" / "output.h5")
    limit = tolerances["mass_balance"]["max_error_fraction_of_rain"]
    print(f"  mass balance: |error| = {abs(error):.4e} m^3 of {rain:.1f} m^3 rain "
          f"({abs(error) / rain:.2e}, allowed {limit:.1e})")
    if abs(error) > limit * rain:
        print("  mass balance: FAIL")
        ok = False
    print("b1 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def wetting_front_depth(profile: np.ndarray, dz: float, threshold: float) -> float:
    """Depth of the first top-down theta crossing of ``threshold`` [m < 0]."""
    z = -(np.arange(profile.size) + 0.5) * dz
    below = np.where(profile < threshold)[0]
    if below.size == 0 or below[0] == 0:
        return float("nan")
    k = int(below[0])
    frac = (threshold - profile[k - 1]) / (profile[k] - profile[k - 1])
    return float(z[k - 1] + frac * (z[k] - z[k - 1]))


def gate_b2(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "b2-gw", args.work)
    run_case(args.frehg, case_dir / "b2-gw.yaml", case_dir, args.mpiexec, args.ranks)

    golden_src = args.legacy / "b2-gw" / "out"
    golden_h5 = args.work / "b2-golden.h5"
    written = ascii_golden_to_h5.convert3d(golden_src, golden_h5,
                                           ["hydraulic_head", "water_content"])
    if written == 0:
        return 1
    with open(HERE / "tolerances" / "b2-gw.yaml", encoding="utf-8") as handle:
        tolerances = yaml.safe_load(handle)
    print(f"b2 gate: comparing against {written} golden datasets")
    ok = compare_h5.compare(case_dir / "out" / "output.h5", golden_h5, tolerances)

    # Wetting-front gate vs the Warrick analytical reference (plan §9 b2).
    spec = tolerances["wetting_front"]
    threshold = float(spec["threshold"])
    reference = {}
    csv_path = args.legacy / "b2-gw" / "reference" / "warrick_water_content_profile.csv"
    with open(csv_path, encoding="utf-8") as handle:
        header = handle.readline().strip().split(",")
        for line in handle:
            row = dict(zip(header, line.strip().split(",")))
            if abs(float(row["water_content"]) - threshold) < 1.0e-9:
                reference[int(row["time_s"])] = float(row["z_m"])
    with h5py.File(case_dir / "out" / "output.h5", "r") as model, \
         h5py.File(golden_h5, "r") as golden:
        for t, z_ref in sorted(reference.items()):
            front_m = wetting_front_depth(model[f"/groundwater/water_content/{t}"][:],
                                          0.01, threshold)
            front_g = wetting_front_depth(golden[f"/groundwater/water_content/{t}"][:],
                                          0.01, threshold)
            err_m = abs(front_m - z_ref)
            err_g = abs(front_g - z_ref)
            rel = err_m / abs(z_ref)
            ratio = err_m / max(err_g, 1.0e-6)
            passed = rel <= spec["max_rel_error_vs_warrick"] and \
                ratio <= spec["max_ratio_vs_legacy_error"]
            print(f"  front@{t}: model {front_m:.4f} m, Warrick {z_ref:.4f} m "
                  f"(rel {rel:.2%} of {spec['max_rel_error_vs_warrick']:.0%}; "
                  f"{ratio:.2f}x legacy error, allowed "
                  f"{spec['max_ratio_vs_legacy_error']:.1f}x) {'ok' if passed else 'FAIL'}")
            ok = ok and passed
    print("b2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def parse_kirkland_reference(makeplot: Path) -> dict:
    """Extract the digitized ref0/ref400 point arrays from the legacy
    b3 makeplot.py (the golden data lives there; goldens are never
    committed)."""
    import ast

    tree = ast.parse(makeplot.read_text())
    points = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and \
                isinstance(node.targets[0], ast.Name) and \
                node.targets[0].id in ("ref0", "ref400"):
            values = ast.literal_eval(node.value)
            points[node.targets[0].id] = np.asarray(values, dtype=float).reshape(-1, 2)
    if set(points) != {"ref0", "ref400"}:
        raise SystemExit(f"error: ref0/ref400 not found in {makeplot}")
    return points


def contour_point_cloud(field: np.ndarray, x: np.ndarray, z: np.ndarray,
                        level: float) -> np.ndarray:
    """Level-set crossings of field[i, k] along grid edges."""
    pts = []
    nxc, nzc = field.shape
    for i in range(nxc):
        for k in range(nzc - 1):
            a, b = field[i, k] - level, field[i, k + 1] - level
            if a * b < 0:
                frac = a / (a - b)
                pts.append((x[i], z[k] + frac * (z[k + 1] - z[k])))
    for k in range(nzc):
        for i in range(nxc - 1):
            a, b = field[i, k] - level, field[i + 1, k] - level
            if a * b < 0:
                frac = a / (a - b)
                pts.append((x[i] + frac * (x[i + 1] - x[i]), z[k]))
    return np.asarray(pts)


def gate_b3(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "b3-kirkland", args.work)
    run_case(args.frehg, case_dir / "b3-kirkland.yaml", case_dir, args.mpiexec, args.ranks)

    with open(HERE / "tolerances" / "b3-kirkland.yaml", encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)
    reference = parse_kirkland_reference(args.legacy / "b3-kirkland" / "makeplot.py")

    nx, nz, dcell = 50, 30, 0.1
    x = (np.arange(nx) + 0.5) * dcell
    z = nz * dcell - (np.arange(nz) + 0.5) * dcell  # height above the box bottom
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        head = handle["/groundwater/hydraulic_head/86400"][:].reshape(nx, nz)
        audit = handle["/monitor/gw_mass_audit"][:]

    for level, ref in ((0.0, reference["ref0"]), (-400.0, reference["ref400"])):
        cloud = contour_point_cloud(head, x, z, level)
        if cloud.size == 0:
            print(f"  h={level:g}: FAIL — no contour in the simulated field "
                  f"(head range [{head.min():.3g}, {head.max():.3g}])")
            ok = False
            continue
        dist = np.sqrt(((ref[:, None, :] - cloud[None, :, :]) ** 2).sum(-1)).min(1)
        max_d, rms = float(dist.max()), float(np.sqrt((dist ** 2).mean()))
        passed = max_d <= spec["contours"]["max_distance_m"] and \
            rms <= spec["contours"]["max_rms_m"]
        print(f"  h={level:g}: max dist {max_d:.3f} m "
              f"(allowed {spec['contours']['max_distance_m']}), rms {rms:.3f} m "
              f"(allowed {spec['contours']['max_rms_m']}) {'ok' if passed else 'FAIL'}")
        ok = ok and passed

    # Internal mass balance: total moisture change vs the strip inflow
    # recorded by the audit (columns: t, volume, boundary_in, ss_storage,
    # realloc, realloc_dropped, vloss).
    d_volume = float(audit[-1, 1] - audit[0, 1])
    inflow = float(audit[-1, 2])
    error = abs(d_volume - inflow) / inflow
    limit = spec["mass_balance"]["max_error_fraction_of_inflow"]
    print(f"  mass balance: |dV - inflow| = {error:.2%} of inflow (allowed {limit:.1%})")
    ok = ok and error <= limit
    print("b3 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_b2_restart(args: argparse.Namespace) -> int:
    # Uninterrupted run with a checkpoint keyed at half time (the adaptive
    # dtg run crosses the boundary mid-step; the checkpoint stores the
    # exact state time and dtg in its header).
    case_dir = stage_case(args.repo, "b2-gw", args.work)
    config = case_dir / "b2-gw.yaml"

    def add_checkpoint(doc):
        doc["output"]["checkpoint"] = {"interval": 23400}

    rewrite_config(config, add_checkpoint)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"

    restart_dir = stage_case(args.repo, "b2-gw", args.work / "restart")
    restart_config = restart_dir / "b2-gw.yaml"

    def add_restart(doc):
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()), "time": 23400}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)

    ok = True
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for var in ["hydraulic_head", "water_content"]:
            for t in ["35100", "46800"]:
                a = full[f"/groundwater/{var}/{t}"][:]
                b = part[f"/groundwater/{var}/{t}"][:]
                scale = max(float(np.abs(a).max()), 1.0e-12)
                diff = float(np.abs(a - b).max()) / scale
                print(f"  restart {var}@{t}: max rel diff = {diff:.3e}")
                if diff > 1.0e-12:
                    ok = False
    print("b2 restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rank_invariance_b2(args: argparse.Namespace) -> int:
    strict = args.mode == "strict"
    extra = STRICT_PETSC_OPTIONS_GW if strict else []
    # Strict mode proves assembly/physics rank invariance at 1e-12 on the
    # first half of the run (through t = 23400 s, the ponded front crossing
    # ~40 cells with the adaptive dtg at its ceiling). Past that, discrete
    # saturation-front cell crossings echo the rank-layout rounding of the
    # PETSc reductions — the threshold-crossing amplification amendment A5
    # measured for b1's wet/dry fronts — and the full horizon measures
    # front chaos, not assembly correctness (measured: <= 3.1e-14 at
    # 23400 s, 9.4e-4 at 46800 s even at machine-precision solves). The
    # default bjacobi/icc lane runs the full horizon against the
    # data-derived amendment-A8 bound.
    limit = 1.0e-12 if strict else 5.0e-3
    outputs = {}
    for ranks in (1, 2, 4):
        case_dir = stage_case(args.repo, "b2-gw", args.work / f"n{ranks}")
        config = case_dir / "b2-gw.yaml"

        # b2's 1x1 column cannot be decomposed; replicate it to 4x4
        # identical columns (lateral conductivity is zero under
        # use_full3d = false, so the physics per column is unchanged) —
        # amendment A8. The ponded-top region widens with the domain so
        # every column keeps the b2 boundary condition.
        def widen(doc):
            doc["domain"]["nx"] = 4
            doc["domain"]["ny"] = 4
            for bc in doc["boundary_conditions"]:
                bc["region"]["polygon"] = [[-0.1, -0.1], [4.1, -0.1],
                                           [4.1, 4.1], [-0.1, 4.1]]
            if strict:
                doc["time"]["t_end"] = 23400

        rewrite_config(config, widen)
        run_case(args.frehg, config, case_dir, args.mpiexec, ranks, extra)
        outputs[ranks] = case_dir / "out" / "output.h5"

    ok = True
    with h5py.File(outputs[1], "r") as base:
        for ranks in (2, 4):
            with h5py.File(outputs[ranks], "r") as other:
                worst = 0.0
                for var in ["hydraulic_head", "water_content"]:
                    for t in base[f"/groundwater/{var}"]:
                        a = base[f"/groundwater/{var}/{t}"][:]
                        b = other[f"/groundwater/{var}/{t}"][:]
                        scale = max(float(np.abs(a).max()), 1.0e-12)
                        worst = max(worst, float(np.abs(a - b).max()) / scale)
                print(f"  n=1 vs n={ranks} [{args.mode}]: max rel diff = {worst:.3e} "
                      f"(allowed {limit:.0e})")
                if worst > limit:
                    ok = False
    print(f"b2 rank invariance ({args.mode}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_b6_gw_smoke(args: argparse.Namespace) -> int:
    """The plan §10 P2 exit item: an early groundwater-only b6 run compared
    against the syncV4 steady-state golden columns. This is a *recorded*
    comparison, not a pass/fail gate — it sizes the P4 tolerance risk
    (PCA vs the golden's Newton scheme, plus the missing surface coupling
    and salinity in a groundwater-only run) before those tolerances freeze.
    The measured deviations are archived in dod-P2.md / report-P2.md."""
    case_dir = stage_case(args.repo, "b6-kuan", args.work)
    config = case_dir / "b6-kuan-ss.yaml"

    def gw_only(doc):
        doc["modules"] = {"surface_water": False, "groundwater": True, "transport": False}
        doc["groundwater"]["density_coupling"] = {"enabled": False}
        for section in ("surface_water", "transport"):
            doc.pop(section, None)
        doc["initial_conditions"] = {
            "groundwater": doc["initial_conditions"]["groundwater"]}
        doc["boundary_conditions"] = [
            bc for bc in doc["boundary_conditions"]
            if bc["target"].startswith("groundwater") and bc["kind"] != "scalar_value"]
        doc["output"]["variables"] = {
            "groundwater": ["hydraulic_head", "water_content"]}
        doc["output"].pop("monitors", None)
        doc["output"]["filename"] = "out/gw-smoke.h5"

    rewrite_config(config, gw_only)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)

    golden_dir = args.legacy / "b6-kuan" / "out-ss-syncV4"
    t_final = 36000
    head_g = np.loadtxt(golden_dir / f"head_{t_final}")
    moisture_g = np.loadtxt(golden_dir / f"moisture_{t_final}")
    with h5py.File(case_dir / "out" / "gw-smoke.h5", "r") as handle:
        head_m = handle[f"/groundwater/hydraulic_head/{t_final}"][:]
        moisture_m = handle[f"/groundwater/water_content/{t_final}"][:]

    finite = np.isfinite(head_m)
    dh = np.abs(head_m - head_g)[finite]
    dth = np.abs(moisture_m - moisture_g)[finite]
    saturated = (moisture_g >= 0.4599)[finite.nonzero()[0]] if finite.all() else \
        (moisture_g[finite] >= 0.4599)
    print(f"b6 gw-only smoke vs syncV4 golden at t = {t_final} s "
          f"({int(finite.sum())} active cells):")
    print(f"  |dh|: max {dh.max():.4f} m, mean {dh.mean():.4f} m; "
          f"saturated zone max {dh[saturated].max():.4f} m")
    print(f"  |dtheta|: max {dth.max():.4f}, mean {dth.mean():.4f}; "
          f"saturated zone max {dth[saturated].max():.4f}")
    print("  (context: the golden is a coupled Newton run with salinity; the")
    print("   P4 b6 tolerances are head <= max(0.01 m, 5 %), salinity gates")
    print("   separate — plan §9. Recorded, not gated, at P2.)")
    return 0


# ---------------------------------------------------------------------------
# b6-kuan (plan §9, P4 blocking gate; both variants)
# ---------------------------------------------------------------------------

B6_DIM = (68, 18)          # (ny, nz)
B6_BOTZ = -0.4             # legacy botZ [m]
B6_DY = 0.05               # [m]
B6_AVG_WINDOW = (30000, 36000, 600)  # tidal-average window (script_plotall.py)


def b6_geometry(legacy_dir: Path):
    """Fine-grid geometry of the legacy plotting script (script_plotall.py).

    Returns (bath, dzf, nzf, actv): the per-column bed elevations, the fine
    layer thickness, the fine layer count, and per-column active fine-layer
    counts."""
    bath = np.genfromtxt(legacy_dir / "b6-kuan" / "ex3_kuan_input" / "bath")
    dzf = (bath.min() - B6_BOTZ) / B6_DIM[1]
    nzf = int(np.ceil((bath.max() - B6_BOTZ) / dzf))
    actv = np.ceil((bath - B6_BOTZ) / dzf).astype(int)
    return bath, dzf, nzf, actv


def b6_column_to_fine(field: np.ndarray, actv: np.ndarray, nzf: int) -> np.ndarray:
    """The legacy c2f mapping: interpolate each column's terrain-following
    cells onto the uniform fine grid, flipped so index 0 is the fine-grid
    top (script_plotall.py c2f + flipud)."""
    nz, ny = field.shape
    out = np.nan * np.ones((nzf, ny))
    xc = np.arange(0, nz, 1)
    for j in range(ny):
        xf = np.linspace(0, nz, actv[j])
        colf = np.interp(xf, xc, field[:, j])
        out[:actv[j], j] = np.flip(colf)
    return np.flipud(out)


def b6_average_salt_model(output: Path) -> np.ndarray:
    """Tidally averaged subsurface salinity (nz, ny) from the model output."""
    t0, t1, dt = B6_AVG_WINDOW
    acc = np.zeros(B6_DIM[::-1])
    count = 0
    with h5py.File(output, "r") as handle:
        for t in range(t0, t1 + 1, dt):
            flat = handle[f"/transport/concentration/{t}"][:]
            acc += np.reshape(flat, (B6_DIM[1], B6_DIM[0]), order="F")
            count += 1
    return acc / count


def b6_average_salt_golden(golden_dir: Path) -> np.ndarray:
    t0, t1, dt = B6_AVG_WINDOW
    acc = np.zeros(B6_DIM[::-1])
    count = 0
    for t in range(t0, t1 + 1, dt):
        flat = np.genfromtxt(golden_dir / f"scalar_subs1__{t}")
        acc += np.reshape(flat, (B6_DIM[1], B6_DIM[0]), order="F")
        count += 1
    return acc / count


def b6_interface_mae(saltFine: np.ndarray, expPoints: np.ndarray, dzf: float,
                     actv: np.ndarray) -> float:
    """MAE [m] of the 50 %-isohaline position against the experimental points.

    For each experimental point (y [m], z [m] with the 0.73 m tank datum of
    script_plotall.py): collect every fine-grid depth where the column's
    normalized salinity crosses 0.5 (the tidal variant's columns are
    non-monotone — the upper saline plume and the deep wedge give two
    crossings, and the experimental points trace both) and take the crossing
    nearest the point; a column without any crossing clamps to its bottom.
    The identical extraction runs on the model and the golden, so the
    1.5x-golden criterion compares like with like."""
    errors = []
    for y, z in expPoints:
        j = int(round(y / B6_DY))
        j = min(max(j, 0), saltFine.shape[1] - 1)
        col = saltFine[:, j] / 35.0
        # Fine index 0 is the fine-grid TOP after the flipud in
        # b6_column_to_fine; valid entries span the column's active layers
        # at the bottom of the fine grid.
        idx = np.where(np.isfinite(col))[0]
        crossings = []
        for a, b in zip(idx[:-1], idx[1:]):
            if b != a + 1:
                continue
            lo, hi = col[a], col[b]
            if (lo - 0.5) * (hi - 0.5) <= 0.0 and lo != hi:
                crossings.append(a + (0.5 - lo) / (hi - lo))
        if not crossings:
            crossings = [float(idx[-1])]  # interface absent: clamp to the bottom
        zExpIndex = (0.73 - z) / dzf  # the plotting script's datum mapping
        errors.append(min(abs(c - zExpIndex) for c in crossings) * dzf)
    del actv  # geometry enters through the NaN mask of saltFine
    return float(np.mean(errors))


def b6_salt_mass(salt_series, moisture_series, bath: np.ndarray) -> float:
    """Tidally averaged subsurface salt mass sum(s theta V) [psu m^3]; the
    scaled terrain mesh gives dz3d = (bath - botZ)/nz per column."""
    dz = (bath - B6_BOTZ) / B6_DIM[1]  # per column
    az = 1.0 * B6_DY
    total = 0.0
    for salt, theta in zip(salt_series, moisture_series):
        # salt/theta are (nz, ny)
        total += float(np.nansum(salt * theta * (az * dz)[None, :]))
    return total / len(salt_series)


def apply_ob8_override(case_dir: Path) -> None:
    """Stage-3 T2 ob8 (nightly record): FREHG_STAGE3_OB8=external sets
    inflow: {closure: external} on every eta condition of the staged b6
    configurations. Unlike the generic §3.2 channel this edits a BC list
    entry, so it has its own hook; the b6 criteria are unchanged (the
    experiment bounds gate, the legacy-golden residuals are recorded)."""
    if os.environ.get("FREHG_STAGE3_OB8", "") != "external":
        return
    for config in sorted(case_dir.glob("*.yaml")):
        with open(config, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        touched = False
        for bc in doc.get("boundary_conditions", []):
            if bc.get("kind") == "eta":
                bc["inflow"] = {"closure": "external"}
                touched = True
        if touched:
            with open(config, "w", encoding="utf-8") as handle:
                yaml.safe_dump(doc, handle, sort_keys=False)
            print(f"  ob8 override: {config.name} -> inflow closure external")


def gate_b6(args: argparse.Namespace) -> int:
    variant = args.variant
    case_dir = args.work / "b6-kuan"
    output = case_dir / f"out-{variant}" / "output.h5"
    if args.reuse_output and output.exists():
        print(f"b6 gate [{variant}]: reusing {output}")
    else:
        case_dir = stage_case(args.repo, "b6-kuan", args.work)
        apply_ob8_override(case_dir)
        config = case_dir / f"b6-kuan-{variant}.yaml"
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)

    with open(HERE / "tolerances" / "b6-kuan.yaml", encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)
    golden_dir = args.legacy / "b6-kuan" / f"out-{variant}-syncV4"
    exp_file = args.legacy / "b6-kuan" / "ex3_kuan_input" / (
        "Kuan_Exp_NoTide" if variant == "ss" else "Kuan_Exp_Tide")
    bath, dzf, nzf, actv = b6_geometry(args.legacy)
    ok = True

    # Stage-3 T2 ob8 (plan §6): under FREHG_STAGE3_OB8=external the run
    # deliberately departs from the LEGACY GOLDEN, which embeds the
    # zero_gradient stage-deposit kick (0.90 m^3 of the 2.31 m^3 b6-td
    # seepage lands on the stage row — dod-Q7's unclosed term). The ob8
    # record therefore gates the EXPERIMENT criteria (interface MAE caps,
    # concentration bounds) and RECORDS every golden-referenced residual.
    # First execution measured: td golden head 0.08 % outside (one beach
    # cell, 0.0196 m vs the 0.015 m floor), salt mass 9.88 % (vs 8.79 %
    # zero_gradient), interface MAE 0.0381 m vs the golden's own 0.0430 m
    # — agreement with the experiment IMPROVES under external.
    ob8 = os.environ.get("FREHG_STAGE3_OB8", "") == "external"

    def golden_gate(bad: bool) -> bool:
        if ob8 and bad:
            print(f"b6-{variant} ob8: golden residual above the zero_gradient "
                  "bound — recorded, not gated (the golden embeds the "
                  "stage-deposit kick)")
            return False
        return bad

    # --- Golden field comparison at quasi-steady state (t = 36000 s).
    t_final = 36000
    head_g = np.loadtxt(golden_dir / f"head_{t_final}")
    salt_g = np.loadtxt(golden_dir / f"scalar_subs1__{t_final}")
    with h5py.File(output, "r") as handle:
        head_m = handle[f"/groundwater/hydraulic_head/{t_final}"][:]
        salt_m = handle[f"/transport/concentration/{t_final}"][:]
    finite = np.isfinite(head_m)
    g = spec["golden"]
    dh = np.abs(head_m - head_g)[finite]
    h_allow = np.maximum(g["head_abs_floor"], g["head_rel"] * np.abs(head_g[finite]))
    h_bad = float(np.mean(dh > h_allow))
    print(f"b6-{variant} golden head: |dh| max {dh.max():.4f} m, mean {dh.mean():.4f} m, "
          f"{100.0 * h_bad:.2f} % of cells outside (allowed 0 %)")
    if golden_gate(h_bad > 0.0):
        ok = False
    ds = np.abs(salt_m - salt_g)[finite]
    s_allow = np.maximum(g["salt_abs_floor"], g["salt_rel"] * np.abs(salt_g[finite]))
    s_bad = float(np.mean(ds > s_allow))
    print(f"b6-{variant} golden salinity: |ds| max {ds.max():.3f} psu, mean {ds.mean():.3f}, "
          f"{100.0 * s_bad:.2f} % of cells outside "
          f"(allowed {100.0 * g['salt_exceed_fraction']:.0f} %)")
    if golden_gate(s_bad > g["salt_exceed_fraction"]):
        ok = False

    # --- Primary gate: the 50 %-isohaline vs the Kuan experiment.
    exp_points = np.genfromtxt(exp_file, delimiter=",")
    salt_avg_m = b6_average_salt_model(output)
    salt_avg_g = b6_average_salt_golden(golden_dir)
    fine_m = b6_column_to_fine(salt_avg_m, actv, nzf)
    fine_g = b6_column_to_fine(salt_avg_g, actv, nzf)
    mae_m = b6_interface_mae(fine_m, exp_points, dzf, actv)
    mae_g = b6_interface_mae(fine_g, exp_points, dzf, actv)
    height = nzf * dzf
    e = spec["experiment"]
    mae_cap = e["interface_mae_of_height"] * height
    ratio_cap = e["interface_mae_golden_ratio"] * mae_g
    print(f"b6-{variant} interface MAE vs experiment: model {mae_m:.4f} m, "
          f"golden {mae_g:.4f} m (caps: {mae_cap:.4f} m = "
          f"{100.0 * e['interface_mae_of_height']:.0f} % of {height:.3f} m tank, "
          f"{ratio_cap:.4f} m = {e['interface_mae_golden_ratio']:.1f}x golden)")
    if mae_m > mae_cap or mae_m > ratio_cap:
        ok = False

    # --- Tidally averaged salt mass vs the golden.
    t0, t1, dt = B6_AVG_WINDOW
    salts_m, thetas_m, salts_g, thetas_g = [], [], [], []
    with h5py.File(output, "r") as handle:
        for t in range(t0, t1 + 1, dt):
            salts_m.append(np.reshape(handle[f"/transport/concentration/{t}"][:],
                                      (B6_DIM[1], B6_DIM[0]), order="F"))
            thetas_m.append(np.reshape(handle[f"/groundwater/water_content/{t}"][:],
                                       (B6_DIM[1], B6_DIM[0]), order="F"))
            salts_g.append(np.reshape(np.genfromtxt(golden_dir / f"scalar_subs1__{t}"),
                                      (B6_DIM[1], B6_DIM[0]), order="F"))
            thetas_g.append(np.reshape(np.genfromtxt(golden_dir / f"moisture_{t}"),
                                       (B6_DIM[1], B6_DIM[0]), order="F"))
    mass_m = b6_salt_mass(salts_m, thetas_m, bath)
    mass_g = b6_salt_mass(salts_g, thetas_g, bath)
    rel = abs(mass_m - mass_g) / abs(mass_g)
    print(f"b6-{variant} tidally averaged salt mass: model {mass_m:.3f}, golden {mass_g:.3f} "
          f"psu m^3 ({100.0 * rel:.2f} % off; allowed "
          f"{100.0 * spec['salt_mass_rel']:.0f} %)")
    if golden_gate(rel > spec["salt_mass_rel"]):
        ok = False

    # --- Salinity bounded in [0, s_boundary] at every output (plan §10 P4).
    b = spec["bounds"]
    worst_lo, worst_hi = 0.0, 0.0
    with h5py.File(output, "r") as handle:
        for t in handle["/transport/concentration"]:
            s3 = handle[f"/transport/concentration/{t}"][:]
            s2 = handle[f"/transport/concentration_surface/{t}"][:]
            for s in (s3[np.isfinite(s3)], s2[np.isfinite(s2)]):
                worst_lo = min(worst_lo, float(s.min()))
                worst_hi = max(worst_hi, float(s.max()))
    print(f"b6-{variant} salinity range over all outputs: [{worst_lo:.3e}, {worst_hi:.6f}] "
          f"(bounds [{b['min']}, {b['max']}] + {b['slack']:.0e})")
    if worst_lo < b["min"] - b["slack"] or worst_hi > b["max"] + b["slack"]:
        ok = False

    print(f"b6-{variant} gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_b6_restart(args: argparse.Namespace) -> int:
    # A shortened coupled+transport run with a checkpoint at half time; the
    # restart must match the uninterrupted run bitwise on every output —
    # the P3 determinism criterion extended over the P4 transport state
    # (scalars, the flow-rate snapshots, the carried top-cell dispersion).
    case_dir = stage_case(args.repo, "b6-kuan", args.work)
    config = case_dir / "b6-kuan-ss.yaml"

    def shorten(doc):
        doc["time"]["t_end"] = 240
        doc["time"]["output_interval"] = 60
        doc["output"]["checkpoint"] = {"interval": 120}
        doc["output"]["filename"] = "out/full.h5"

    rewrite_config(config, shorten)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "full.h5"

    restart_dir = stage_case(args.repo, "b6-kuan", args.work / "restart")
    restart_config = restart_dir / "b6-kuan-ss.yaml"

    def add_restart(doc):
        shorten(doc)
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()), "time": 120}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)

    ok = True
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        checks = [("surface", v) for v in ("eta", "depth", "uu", "vv", "seepage")] + \
                 [("groundwater", v) for v in ("hydraulic_head", "water_content")] + \
                 [("transport", v) for v in ("concentration", "concentration_surface")]
        for group, var in checks:
            for t in ["180", "240"]:
                a = full[f"/{group}/{var}/{t}"][:]
                b = part[f"/{group}/{var}/{t}"][:]
                finite = np.isfinite(a)
                scale = max(float(np.abs(a[finite]).max()), 1.0e-12)
                diff = float(np.abs(a[finite] - b[finite]).max()) / scale
                print(f"  restart {var}@{t}: max rel diff = {diff:.3e}")
                if diff > 1.0e-12:
                    ok = False
    print("b6 restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_smoke_b5(args: argparse.Namespace) -> int:
    # Sanitizer path-coverage run (plan §10 P5, amendment A21): the four
    # nightly-class b5 envelope runs are hours long even uninstrumented, so
    # the sanitizer matrix executes every scenario x coupling combination on
    # a shortened horizon instead. The plan §9 gate metrics are NOT applied
    # here (the full-horizon runs own them); the sanitizers themselves are
    # the check — any finding aborts the run (halt_on_error) and fails this
    # subcommand through the nonzero exit.
    case_dir = stage_case(args.repo, "b5-vcatchment", args.work)
    config_name = "b5-vcatchment.yaml" if args.scenario == "rain" else "b5-vcatchment-norain.yaml"
    config = case_dir / config_name
    horizon = float(args.t_end) if args.t_end is not None else 1800.0

    def shorten(doc):
        if args.coupling == "subcycled":
            doc["coupling"]["mode"] = "subcycled"
        doc["time"]["t_end"] = horizon
        doc["time"]["output_interval"] = min(float(doc["time"]["output_interval"]), horizon)
        doc["output"]["checkpoint"] = {"interval": horizon / 2}

    rewrite_config(config, shorten)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    print(f"b5 smoke [{args.scenario}, {args.coupling}, t_end={horizon:g} s]: "
          "PASS (clean instrumented run)")
    return 0


def gate_smoke_b6(args: argparse.Namespace) -> int:
    # Sanitizer path-coverage run for the two b6 variants (plan §10 P5,
    # amendment A21): full transport + baroclinic + (td) tidal-boundary and
    # sea-surface-salinity paths on a shortened horizon. As with smoke-b5,
    # the sanitizers are the check, not the plan §9 metrics.
    case_dir = stage_case(args.repo, "b6-kuan", args.work)
    config = case_dir / f"b6-kuan-{args.variant}.yaml"
    horizon = float(args.t_end) if args.t_end is not None else 1800.0

    def shorten(doc):
        doc["time"]["t_end"] = horizon
        doc["time"]["output_interval"] = min(float(doc["time"]["output_interval"]), horizon)
        doc["output"]["checkpoint"] = {"interval": horizon / 2}

    rewrite_config(config, shorten)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    print(f"b6 smoke [{args.variant}, t_end={horizon:g} s]: PASS (clean instrumented run)")
    return 0


def gate_smoke_g5(args: argparse.Namespace) -> int:
    # Sanitizer path-coverage run for the nightly-class g5 case (plan §9 Q7
    # row: the sanitizer matrix covers the new labels' code paths): bulk
    # evaporation + exact evaporative concentration + baroclinic density on
    # the Geng slice, shortened. Same contract as smoke-b5/b6 — the
    # sanitizers are the check, not the g5 metrics.
    case_dir = stage_case(args.repo, "g5-geng2015", args.work)
    config = case_dir / "g5-geng2015.yaml"
    horizon = float(args.t_end) if args.t_end is not None else 1800.0

    def shorten(doc):
        doc["time"]["t_end"] = horizon
        doc["time"]["output_interval"] = min(float(doc["time"]["output_interval"]), horizon)
        doc["output"]["checkpoint"] = {"interval": horizon / 2}

    rewrite_config(config, shorten)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    print(f"g5 smoke [t_end={horizon:g} s]: PASS (clean instrumented run)")
    return 0


def gate_smoke_g8(args: argparse.Namespace) -> int:
    # Sanitizer path-coverage run for the nightly-class g8 case: subsurface
    # temperature transport + thermal density feedback at the supercritical
    # conductivity (the convecting path), shortened.
    case_dir = stage_case(args.repo, "g8-hrl", args.work)
    config = case_dir / "g8-hrl.yaml"
    horizon = float(args.t_end) if args.t_end is not None else 1800.0

    def shorten(doc):
        for soil in doc["soil"]["types"]:
            soil["ksx"] = soil["ksy"] = soil["ksz"] = 1.375e-2  # supercritical
        doc["time"]["t_end"] = horizon
        doc["time"]["output_interval"] = min(float(doc["time"]["output_interval"]), horizon)
        doc["output"]["checkpoint"] = {"interval": horizon / 2}

    rewrite_config(config, shorten)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    print(f"g8 smoke [supercritical, t_end={horizon:g} s]: PASS (clean instrumented run)")
    return 0


def hydrograph_from_audit(output: Path, halfwidth: float):
    """Times (2 s cadence) and outflow rate from the mass-audit table."""
    with h5py.File(output, "r") as handle:
        table = handle["/monitor/mass_audit"][:]
    t, cum = table[:, 0], table[:, 4]

    def rate(x: np.ndarray) -> np.ndarray:
        hi = np.interp(x + halfwidth, t, cum)
        lo = np.interp(x - halfwidth, t, cum)
        return (hi - lo) / (2.0 * halfwidth)

    sample_t = np.arange(2.0, float(t[-1]) - 1.0, 2.0)
    return sample_t, rate(sample_t), rate, float(cum[-1])


def gate_b4(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "b4-govindaraju", args.work)
    run_case(args.frehg, case_dir / "b4-govindaraju.yaml", case_dir, args.mpiexec, args.ranks)

    with open(HERE / "tolerances" / "b4-govindaraju.yaml", encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)["hydrograph"]
    nx, ny, dx = 200, 10, 0.109725
    area = nx * ny * dx * dx

    reference = np.loadtxt(args.legacy / "b4-govindaraju" / "ReferenceData" / "outflow.txt",
                           skiprows=1)
    ref_t, ref_q = reference[:, 0], reference[:, 1]

    sample_t, sample_q, rate, total_out = hydrograph_from_audit(
        case_dir / "out" / "output.h5", spec["smoothing_halfwidth_seconds"])
    sample_q = sample_q / area
    model_at_ref = rate(ref_t) / area

    rel_l2 = float(np.sqrt(((model_at_ref - ref_q) ** 2).sum()) / np.sqrt((ref_q ** 2).sum()))
    peak_m, peak_r = float(sample_q.max()), float(ref_q.max())
    frac = spec["time_to_peak_fraction"]
    ttp_m = float(sample_t[np.argmax(sample_q >= frac * peak_m)])
    ttp_r = float(ref_t[np.argmax(ref_q >= frac * peak_r)])

    rain = np.loadtxt(args.repo / "benchmarks" / "b4-govindaraju" / "input" / "rain.dat")
    rain_vol = 0.0
    t_end = float(sample_t[-1])
    for k in range(len(rain) - 1):
        left, right = float(rain[k, 0]), min(float(rain[k + 1, 0]), t_end)
        if left >= t_end:
            break
        rain_vol += 0.5 * (rain[k, 1] + rain[k + 1, 1]) * (right - left)
    rain_vol *= area
    ref_vol = float(np.trapz(ref_q * area, ref_t))

    checks = [
        ("rel_L2", rel_l2, spec["rel_l2_max"], rel_l2 <= spec["rel_l2_max"]),
        ("peak", peak_m / peak_r - 1.0, spec["peak_rel_tol"],
         abs(peak_m / peak_r - 1.0) <= spec["peak_rel_tol"]),
        (f"time-to-{frac:.0%}-peak", ttp_m / ttp_r - 1.0, spec["time_to_peak_rel_tol"],
         abs(ttp_m / ttp_r - 1.0) <= spec["time_to_peak_rel_tol"]),
        ("volume-vs-rain", total_out / rain_vol - 1.0, spec["volume_vs_rain_tol"],
         abs(total_out / rain_vol - 1.0) <= spec["volume_vs_rain_tol"]),
        ("volume-vs-reference", total_out / ref_vol - 1.0, spec["volume_vs_reference_tol"],
         abs(total_out / ref_vol - 1.0) <= spec["volume_vs_reference_tol"]),
    ]
    ok = True
    for name, achieved, allowed, passed in checks:
        print(f"  {name}: achieved {achieved:+.4f} (allowed {allowed:.2f}) "
              f"{'ok' if passed else 'FAIL'}")
        ok = ok and passed
    print("b4 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def rewrite_config(path: Path, transform) -> None:
    with open(path, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    transform(doc)
    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(doc, handle, sort_keys=False)


# ---------------------------------------------------------------------------
# b5 (plan §9): coupled tilted-V envelope gate.
# ---------------------------------------------------------------------------

# The reference CSVs digitize the published intercomparison curves for the
# FULL tilted-V; the model runs the symmetric half (one hillslope plus half
# the channel), exactly as the SERGHEI reference run did, so model volumes
# and discharges double before comparison (plot_discharge.py: "2.0 * data").
B5_SYMMETRY_FACTOR = 2.0
B5_RAIN_END_SECONDS = 72000.0


def b5_reference_curves(legacy: Path, kind: str, scenario: str) -> dict:
    """Digitized reference curves {model: (t_hours, value)} for 'discharge'
    [m^3/h] or 'ponding' [m^3] and the given scenario."""
    names = {
        "discharge": [f"PF-{scenario}.csv", f"CATHY-{scenario}.csv", f"HGS-{scenario}.csv",
                      f"cast3m-{scenario}.csv"],
        "ponding": [f"PF-ponding-{scenario}.csv", f"CATHY-ponding-{scenario}.csv",
                    f"HGS-ponding-{scenario}.csv", f"cast3m-ponding-{scenario}.csv"],
    }[kind]
    curves = {}
    for name in names:
        path = legacy / "b5-vcatchment" / kind / name
        if not path.exists():
            continue  # cast3m has no rain-scenario discharge curve
        data = np.loadtxt(path, delimiter=",")
        order = np.argsort(data[:, 0])
        curves[name.split("-")[0]] = (data[order, 0], np.clip(data[order, 1], 0.0, None))
    return curves


def b5_model_series(output: Path) -> tuple:
    """(t_hours, discharge m^3/h, ponding m^3) of the full V from the
    mass-audit table (the P1 finding: the audit is the mass-consistent
    discharge source; stored velocities under-read)."""
    with h5py.File(output, "r") as handle:
        table = handle["/monitor/mass_audit"][:]
    t = table[:, 0]
    cum_out = table[:, 4]
    volume = table[:, 1]

    # Smooth the per-step outflow into an hourly-scale rate (the reference
    # curves resolve hours; the audit resolves single adaptive steps).
    half = 900.0  # seconds
    sample_t = np.arange(half, float(t[-1]) - half, 600.0)
    hi = np.interp(sample_t + half, t, cum_out)
    lo = np.interp(sample_t - half, t, cum_out)
    q = (hi - lo) / (2.0 * half) * 3600.0 * B5_SYMMETRY_FACTOR  # m^3/h, full V
    pond = np.interp(sample_t, t, volume) * B5_SYMMETRY_FACTOR
    return sample_t / 3600.0, q, pond


def b5_envelope_check(name: str, t_model, v_model, curves: dict, spec: dict) -> bool:
    """Plan §9 b5: at every reference time, the model inside
    [env_min - margin, env_max + margin] with margin = 0.1 Qpeak floored at
    the per-kind minimum (the digitized no-rain ponding values sit at the
    digitization noise floor — a fraction-of-peak band below physical
    meaning would gate digitization noise); peak within 15 % of the
    envelope-mean peak (skipped below the floor for the same reason);
    integrated volume within 10 %."""
    tmax_model = float(t_model[-1])
    spans = [(c[0][0], c[0][-1]) for c in curves.values()]
    t_lo = max(s[0] for s in spans)
    t_hi = min(min(s[1] for s in spans), tmax_model)

    # The envelope at each model's own digitized times (restricted to the
    # window every model covers).
    ref_times = np.unique(np.concatenate([c[0] for c in curves.values()]))
    ref_times = ref_times[(ref_times >= t_lo) & (ref_times <= t_hi)]
    env = np.array([[np.interp(tt, c[0], c[1]) for c in curves.values()] for tt in ref_times])
    env_min = env.min(axis=1)
    env_max = env.max(axis=1)
    env_mean = env.mean(axis=1)
    peak_mean = float(env_mean.max())
    floor = float(spec.get("min_margin", {}).get(name, 0.0))
    margin = max(spec["envelope_margin_of_peak"] * peak_mean, floor)

    model_at_ref = np.interp(ref_times, t_model, v_model)
    low = env_min - margin
    high = env_max + margin
    bad = (model_at_ref < low) | (model_at_ref > high)
    worst = float(np.max(np.maximum(low - model_at_ref, model_at_ref - high) / max(peak_mean, 1e-12)))
    print(f"  {name}: {int(bad.sum())}/{len(ref_times)} reference times outside the envelope "
          f"(margin {spec['envelope_margin_of_peak']:.0%} of peak {peak_mean:.3g}; "
          f"worst overshoot {worst:+.2%} of peak)")
    ok = not bad.any()

    peak_model = float(np.max(np.interp(np.linspace(t_lo, t_hi, 2000), t_model, v_model)))
    if peak_mean <= floor:
        print(f"  {name} peak: envelope mean {peak_mean:.3g} below the {floor:.3g} floor — "
              f"model {peak_model:.3g} gated by the envelope band only")
    else:
        peak_err = peak_model / peak_mean - 1.0
        peak_ok = abs(peak_err) <= spec["peak_rel_tol"]
        print(f"  {name} peak: model {peak_model:.3g} vs envelope mean {peak_mean:.3g} "
              f"({peak_err:+.1%}, allowed {spec['peak_rel_tol']:.0%}) {'ok' if peak_ok else 'FAIL'}")
        ok = ok and peak_ok

    grid = np.linspace(t_lo, t_hi, 2000)
    vol_model = float(np.trapz(np.interp(grid, t_model, v_model), grid))
    vol_mean = float(np.trapz(np.interp(grid, ref_times, env_mean), grid))
    if abs(vol_mean) <= floor * (t_hi - t_lo):
        print(f"  {name} integral over [{t_lo:.1f}, {t_hi:.1f}] h below the floor — "
              f"model {vol_model:.4g} gated by the envelope band only")
        return ok
    vol_err = vol_model / max(vol_mean, 1e-12) - 1.0
    vol_ok = abs(vol_err) <= spec["volume_rel_tol"]
    print(f"  {name} integral over [{t_lo:.1f}, {t_hi:.1f}] h: model {vol_model:.4g} vs "
          f"envelope mean {vol_mean:.4g} ({vol_err:+.1%}, allowed {spec['volume_rel_tol']:.0%}) "
          f"{'ok' if vol_ok else 'FAIL'}")
    return ok and vol_ok


def gate_b5(args: argparse.Namespace) -> int:
    scenario = args.scenario
    coupling = args.coupling
    output = args.work / "b5-vcatchment" / "out" / "output.h5"
    if args.reuse_output and output.exists():
        # Re-apply the envelope checks to an existing run (the simulations
        # are hours long; the checks are seconds).
        print(f"b5 gate [{scenario}, {coupling}]: reusing {output}")
        case_dir = args.work / "b5-vcatchment"
    else:
        case_dir = stage_case(args.repo, "b5-vcatchment", args.work)
        config_name = "b5-vcatchment.yaml" if scenario == "rain" else "b5-vcatchment-norain.yaml"
        config = case_dir / config_name

        def adjust(doc):
            if coupling == "subcycled":
                doc["coupling"]["mode"] = "subcycled"
            if args.t_end is not None:
                # Short-horizon variant: the envelope checks clip to the
                # covered window automatically. The full-horizon runs remain
                # the plan §9 gate (peak/recession need the whole curve);
                # a 6 h rain window covers onset and ramp in ~1/12 the wall
                # time, a >= 24 h norain window the seepage onset.
                doc["time"]["t_end"] = float(args.t_end)
                doc["time"]["output_interval"] = min(
                    doc["time"]["output_interval"], 3600)
                doc["output"]["checkpoint"] = {"interval": 0}

        rewrite_config(config, adjust)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)

    with open(HERE / "tolerances" / "b5-vcatchment.yaml", encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)["envelope"]
    t_model, q_model, pond_model = b5_model_series(case_dir / "out" / "output.h5")

    ok = True
    for kind, series in (("discharge", q_model), ("ponding", pond_model)):
        curves = b5_reference_curves(args.legacy, kind, scenario)
        if not curves:
            print(f"error: no reference curves for {kind}/{scenario}")
            return 1
        print(f"b5 [{scenario}, {coupling}] {kind} vs {sorted(curves)}:")
        ok = b5_envelope_check(kind, t_model, series, curves, spec) and ok

    # The coupled budgets stay closed over the full horizon (plan §10 P3):
    # surface closure to the b1-class bound, subsurface to rounding.
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        ma = handle["/monitor/mass_audit"][:]
        ga = handle["/monitor/gw_mass_audit"][:]
    vol, rain, evap, out, bcin, seep, clamped = ma[-1, 1:8]
    # The below-bed clamp is measured into the audit (column 'clamped'), so
    # the closure identity holds with the legacy defect visible: it is
    # reported here as a defect size, not hidden in the residual.
    surf_resid = vol - (rain - evap - out + bcin + seep + clamped)
    gw_resid = (ga[-1, 1] - ga[0, 1]) - (ga[-1, 2] - ga[-1, 3] + ga[-1, 4] - ga[-1, 6])
    scale = max(rain, ga[0, 1])
    print(f"  budgets: surface residual {surf_resid:.3e} m^3, subsurface residual "
          f"{gw_resid:.3e} m^3 (allowed {spec['budget_rel_tol']:.1e} of {scale:.3g} m^3); "
          f"below-bed clamp creation {clamped:.3g} m^3 ({clamped / max(rain, 1e-12):.2%} of rain)")
    budget_ok = abs(surf_resid) <= spec["budget_rel_tol"] * scale and \
        abs(gw_resid) <= spec["budget_rel_tol"] * scale
    ok = ok and budget_ok

    print(f"b5 gate [{scenario}, {coupling}]:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_b5_restart(args: argparse.Namespace) -> int:
    # A shortened coupled run (rain scenario, committed sync mode) with a
    # checkpoint at half time; the restart must match the uninterrupted run
    # to 1e-12 on every prognostic field (plan §10 P3). The coupled restart
    # state adds the seepage accumulator, qss, and the adaptive-clock
    # scalars to the P1/P2 sets.
    case_dir = stage_case(args.repo, "b5-vcatchment", args.work)
    config = case_dir / "b5-vcatchment.yaml"

    def shorten(doc):
        doc["time"]["t_end"] = 7200
        doc["time"]["output_interval"] = 1800
        doc["output"]["checkpoint"] = {"interval": 3600}

    rewrite_config(config, shorten)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"

    restart_dir = stage_case(args.repo, "b5-vcatchment", args.work / "restart")
    restart_config = restart_dir / "b5-vcatchment.yaml"

    def add_restart(doc):
        shorten(doc)
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()), "time": 3600}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)

    ok = True
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        checks = [("surface", v) for v in ("eta", "depth", "uu", "vv", "seepage")] + \
                 [("groundwater", v) for v in ("hydraulic_head", "water_content")]
        for group, var in checks:
            for t in ["5400", "7200"]:
                a = full[f"/{group}/{var}/{t}"][:]
                b = part[f"/{group}/{var}/{t}"][:]
                finite = np.isfinite(a)
                scale = max(float(np.abs(a[finite]).max()), 1.0e-12)
                diff = float(np.abs(a[finite] - b[finite]).max()) / scale
                print(f"  restart {var}@{t}: max rel diff = {diff:.3e}")
                if diff > 1.0e-12:
                    ok = False
    print("b5 restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rank_invariance_b5(args: argparse.Namespace) -> int:
    strict = args.mode == "strict"
    # Trajectory comparison is unavailable for the coupled case (amendment
    # A14): the adaptive controller forks on a single theta_s-boundary
    # cell's ULP (measured: the dt sequences separate at step 3, 7.8 s vs
    # 1.84 s, from a 2e-6 water-content difference at one water-table
    # cell), and even fixed-dtg runs drift to 2e-1 relative within 600 s
    # through the exchange's binary thresholds — the A5/A8
    # threshold-amplification class at coupled density. The strict lane
    # therefore proves the coupled assembly and exchange operators are
    # rank-invariant over one step (bound 1e-9 per V2-A9 — rounding-level
    # on the dev machine, one platform-rounding threshold flip of footprint
    # 2.2e-10 on the x86-64 CI runner at n=4); the
    # default lane gates the rank-robust bulk observables over the plan
    # §8.2 600 s window and prints the per-field number for the record.
    extra = (STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW) if strict else []
    outputs = {}
    for ranks in (1, 2, 4):
        case_dir = stage_case(args.repo, "b5-vcatchment", args.work / f"n{ranks}")
        config = case_dir / "b5-vcatchment.yaml"

        def shorten(doc):
            doc["time"]["t_end"] = 5 if strict else 600
            doc["time"]["output_interval"] = 5 if strict else 300
            doc["output"]["checkpoint"] = {"interval": 0}

        rewrite_config(config, shorten)
        run_case(args.frehg, config, case_dir, args.mpiexec, ranks, extra)
        outputs[ranks] = case_dir / "out" / "output.h5"

    ok = True
    with h5py.File(outputs[1], "r") as base:
        base_ma = base["/monitor/mass_audit"][:]
        base_ga = base["/monitor/gw_mass_audit"][:]
        for ranks in (2, 4):
            with h5py.File(outputs[ranks], "r") as other:
                worst = 0.0
                per_field: dict[str, float] = {}
                for group, var in [("surface", "eta"), ("surface", "depth"),
                                   ("surface", "seepage"), ("groundwater", "hydraulic_head"),
                                   ("groundwater", "water_content")]:
                    field_worst = 0.0
                    for t in base[f"/{group}/{var}"]:
                        a = base[f"/{group}/{var}/{t}"][:]
                        b = other[f"/{group}/{var}/{t}"][:]
                        finite = np.isfinite(a)
                        scale = max(float(np.abs(a[finite]).max()), 1.0e-12)
                        field_worst = max(
                            field_worst,
                            float(np.abs(a[finite] - b[finite]).max()) / scale)
                    per_field[f"{group}/{var}"] = field_worst
                    worst = max(worst, field_worst)
                if strict:
                    # Bound 1e-9, re-derived 2026-09-20 (V2-A9) from the first
                    # true multi-rank CI execution (the runner's MPICH launched
                    # singletons until the noble PMI fix, so every earlier CI
                    # "pass" compared serial runs). The dev-machine record is
                    # rounding-level (9e-15 A14 / 5.7e-14 dod-P3, arm64 + Apple
                    # libm); the x86-64 glibc/f2cblaslapack runner measures a
                    # stable 2.189e-10 at n=4 (n=2: 8.5e-15) — one discrete
                    # branch of the A14 threshold class flips inside the single
                    # step under that platform's rounding, and A14 already
                    # establishes no solver tolerance removes the class. 1e-9
                    # is ~4.5x above the measured flip footprint and still
                    # orders below real assembly/exchange rank-dependence,
                    # which shows across whole halo lines (>= ~1e-6), not as
                    # an isolated-cell footprint. The per-field breakdown is
                    # printed so any future breach names the flipped field.
                    limit = 1.0e-9
                    print(f"  n=1 vs n={ranks} [strict, one step]: max rel field diff = "
                          f"{worst:.3e} (allowed {limit:.0e}, V2-A9)")
                    for name, val in sorted(per_field.items()):
                        print(f"    {name}: {val:.3e}")
                    ok = ok and worst <= limit
                    continue
                ma = other["/monitor/mass_audit"][:]
                ga = other["/monitor/gw_mass_audit"][:]
                # The surface holds only thin films over this window
                # (~0.7 m^3), so its difference is scaled by the exchanged
                # volume (~93 m^3), not by itself: front-cell wet/dry
                # states flip under any rounding change (measured 9.7e-2
                # self-relative under the sanitizer binary, 7e-4 against
                # the exchange).
                exchanged = max(abs(base_ma[-1, 6]), 1.0e-12)
                # Default-mode bounds are data-derived per solver (A5/A8/A14
                # discipline). The bjacobi-icc numbers are the A14 record.
                # AMG-class preconditioners are not rank-invariant by
                # construction (BoomerAMG/GAMG coarsening depends on the
                # decomposition), and the solver *choice alone* at a fixed
                # decomposition already moves these observables past the
                # bjacobi bounds (measured at Q1: subsurface 5.9e-5,
                # exchanged 4.4e-3 between bjacobi and amg at n=1), so the
                # amg/gamg bounds are derived from the Q1 measurements
                # (subsurface 1.8e-4, exchanged 8.1e-3, surface 1.4e-3 at
                # 600 s) with ~2x headroom — amendment V2-A4.
                if SOLVER_OVERRIDE in ("amg", "gamg"):
                    bounds = {"subsurface": 4.0e-4, "exchanged": 1.6e-2,
                              "surface": 3.0e-3}
                else:
                    bounds = {"subsurface": 3.0e-5, "exchanged": 2.0e-3,
                              "surface": 2.0e-3}
                checks = [
                    ("subsurface volume", base_ga[-1, 1], ga[-1, 1],
                     abs(base_ga[-1, 1]), bounds["subsurface"]),
                    ("exchanged volume", base_ma[-1, 6], ma[-1, 6], exchanged,
                     bounds["exchanged"]),
                    ("surface volume", base_ma[-1, 1], ma[-1, 1], exchanged,
                     bounds["surface"]),
                ]
                print(f"  n=1 vs n={ranks} [default, 600 s]: max rel field diff = {worst:.3e} "
                      "(recorded, not gated — trajectory chaos, amendment A14)")
                for name, va, vb, scale, bound in checks:
                    rel = abs(va - vb) / max(scale, 1.0e-12)
                    passed = rel <= bound
                    print(f"    {name}: {rel:.3e} (allowed {bound:.0e}) "
                          f"{'ok' if passed else 'FAIL'}")
                    ok = ok and passed
    print(f"b5 rank invariance ({args.mode}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# v2 Q4 gates (plan §3.3, V2-A13/V2-A14): g4 analytic drawdown + evaporative
# concentration (per-PR), g5 Geng & Boufadel 2015 bare-soil salinization
# (nightly-class). Authored gate-first per §6.1: g4(c)/(d) and all of g5
# FAIL against pre-Q4 code — (c)/g5 because the schema rejects the Q4
# target keys ("capability absent"), (d) because the pre-Q4 evaporation
# ledger books the potential rate after dry-out. The pure check functions
# live in gate_evap.py so the §6.3 negative battery
# (scripts/test_g45_gates.py) can prove each criterion able to fail.
# ---------------------------------------------------------------------------
import gate_evap  # noqa: E402

# ---------------------------------------------------------------------------
# v2 Q6 gates (plan §5.3): g9 steady wind setup (a nonlinear TELEMAC
# replica / b linear Garratt / c sloping-bottom quadrature), g10 Merian
# seiche relaxation, and the §8.1/§9 8-orientation symmetry battery.
# Authored gate-first per §6.1: (b)/(c)/orient use the Q6 target schema
# (law/cap/u10/v10) and fail "capability absent" until it lands; (a) and
# g10's basis run on configuration the schema already accepts — their
# first execution is also the first time the v1 wind path is exercised
# at all (completion report §6: implemented, ungated).
# ---------------------------------------------------------------------------
import gate_wind  # noqa: E402

# ---------------------------------------------------------------------------
# v2 Q5 gates (plan §4.2, realization V2-A17): g6 subsurface heat
# (B&P Peclet sweep / Ogata-Banks breakthrough / Stallman diel damping +
# the composed coupled variant), g7 surface heat exchange (Edinger
# relaxation + ledger / bulk equilibrium consistency / channel plume),
# g8 Horton-Rogers-Lapwood convection onset, the heat 8-orientation
# battery, the two-scalar restart-determinism lane, and the heat
# rank-invariance lane. Authored gate-first per §6.1: every case uses
# the Q5 target schema (modules.temperature) and fails "capability
# absent" until the capability lands. The closed forms and check
# functions live in gate_heat.py; the §6.3 negative battery is
# scripts/test_g678_gates.py.
# ---------------------------------------------------------------------------
import gate_heat  # noqa: E402

CAPABILITY_ABSENT_Q5 = ("capability absent: the schema rejects the Q5 "
                        "temperature keys (expected until the Q5 capability "
                        "lands — plan §6.1)")


def monitor_series(handle, name: str):
    table = handle[f"/monitor/{name}"][:]
    return table[:, 0], table[:, 1]


def surface_eta(handle, t: int, ny: int, nx: int) -> np.ndarray:
    return handle[f"/surface/eta/{t}"][:].reshape(ny, nx)


def gate_g9(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g9-wind", args.work)
    with open(HERE / "tolerances" / "g9-wind.yaml", encoding="utf-8") as handle:
        tol = yaml.safe_load(handle)
    ok = True

    def report(sub: str, sub_ok: bool, msgs: list[str]) -> bool:
        for msg in msgs:
            print(f"  {msg}")
        print(f"g9({sub}):", "ok" if sub_ok else "FAIL")
        return sub_ok

    # (a) nonlinear setup, constant-Cd law (the v1 wind path's first gate).
    config_a = case_dir / "g9a-setup.yaml"
    with open(config_a, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    run_case(args.frehg, config_a, case_dir, args.mpiexec, args.ranks)
    out_a = case_dir / "out" / "g9a.h5"
    shutil.move(case_dir / "out" / "output.h5", out_a)
    cd = float(doc["surface_water"]["wind"]["cd"])
    u10 = 5.0  # the ramp's held value
    taup = gate_wind.kinematic_stress(cd, u10)
    nx, ny, dx = 200, 5, 2.5  # the §6.3 convergence-study resolution
    x_c = (np.arange(nx) + 0.5) * dx
    h_ref = gate_wind.setup_flat(taup, nx * dx, 2.0, x_c)
    with h5py.File(out_a, "r") as handle:
        t_end = int(doc["time"]["t_end"])
        eta = surface_eta(handle, t_end, ny, nx)[ny // 2]
        tm, vm = monitor_series(handle, "eta_east")
    delta_ref = float(h_ref[-1] - h_ref[0])
    s_ok, s_msgs = gate_wind.check_steady(
        tm, vm, 0.1, float(tol["g9a"]["steady_tol_m"]), "g9a")
    d_ok, d_msgs = gate_wind.check_setup(
        float(eta[-1] - eta[0]), delta_ref, float(tol["g9a"]["setup_tol_rel"]), "g9a")
    p_ok, p_msgs = gate_wind.check_profile(
        eta, h_ref - 2.0, delta_ref, float(tol["g9a"]["profile_tol_fraction"]), "g9a")
    ok &= report("a", s_ok and d_ok and p_ok, s_msgs + d_msgs + p_msgs)

    # (b) linear regime, Garratt law through the component form.
    config_b = case_dir / "g9b-linear.yaml"
    if not validate_config(args.frehg, config_b):
        ok &= report("b", False,
                     ["capability absent: the schema rejects the Q6 wind "
                      "law/cap/u10/v10 keys (expected until the Q6 "
                      "capability lands — plan §6.1)"])
    else:
        run_case(args.frehg, config_b, case_dir, args.mpiexec, args.ranks)
        out_b = case_dir / "out" / "g9b.h5"
        shutil.move(case_dir / "out" / "output.h5", out_b)
        u10 = 15.0
        cd = gate_wind.cd_garratt(u10)  # independent python reference
        taup = gate_wind.kinematic_stress(cd, u10)
        nx, dx = 50, 100.0
        x_c = (np.arange(nx) + 0.5) * dx
        h_ref = gate_wind.setup_flat(taup, nx * dx, 10.0, x_c)
        with h5py.File(out_b, "r") as handle:
            eta = surface_eta(handle, 40000, 1, nx)[0]
            tm, vm = monitor_series(handle, "eta_east")
        delta_ref = float(h_ref[-1] - h_ref[0])
        s_ok, s_msgs = gate_wind.check_steady(
            tm, vm, 0.1, float(tol["g9b"]["steady_tol_m"]), "g9b")
        d_ok, d_msgs = gate_wind.check_setup(
            float(eta[-1] - eta[0]), delta_ref, float(tol["g9b"]["setup_tol_rel"]),
            "g9b (Garratt Cd computed offline)")
        ok &= report("b", s_ok and d_ok, s_msgs + d_msgs)

    # (c) sloping bottom vs the quadrature reference.
    config_c = case_dir / "g9c-slope.yaml"
    if not validate_config(args.frehg, config_c):
        ok &= report("c", False,
                     ["capability absent: the schema rejects the Q6 wind "
                      "law/cap/u10/v10 keys (expected until the Q6 "
                      "capability lands — plan §6.1)"])
    else:
        run_case(args.frehg, config_c, case_dir, args.mpiexec, args.ranks)
        out_c = case_dir / "out" / "g9c.h5"
        shutil.move(case_dir / "out" / "output.h5", out_c)
        u10 = 12.0
        taup = gate_wind.kinematic_stress(gate_wind.cd_garratt(u10), u10)
        nx, dx = 200, 5.0  # the convergence-study resolution (case README)

        x_c = (np.arange(nx) + 0.5) * dx

        def bed_at(xx):
            return -(0.12 + 1.88 * np.asarray(xx, dtype=float) / 1000.0)

        xq, eta_q = gate_wind.setup_slope(taup, nx * dx, bed_at)
        eta_ref = np.interp(x_c, xq, eta_q)
        with h5py.File(out_c, "r") as handle:
            eta = surface_eta(handle, 40000, 1, nx)[0]
            # Steadiness gates at the deep end; the shallow tip carries the
            # P3-documented thin-layer drag limit cycle, recorded below.
            tm, vm = monitor_series(handle, "eta_east")
            tw, vw = monitor_series(handle, "eta_west")
        rng = float(eta_ref[-1] - eta_ref[0])
        h_min = float(np.min(eta - bed_at(x_c)))
        s_ok, s_msgs = gate_wind.check_steady(
            tm, vm, 0.1, float(tol["g9c"]["steady_tol_m"]), "g9c (deep end)")
        tip = vw[tw >= tw[-1] * 0.9]
        s_msgs.append(f"recorded (not gated): shallow-tip trailing span "
                      f"{float(np.max(tip) - np.min(tip)):.2e} m (the thin-layer "
                      f"drag limit cycle)")
        p_ok, p_msgs = gate_wind.check_profile_rms(
            eta, eta_ref, rng, float(tol["g9c"]["profile_tol_fraction"]), "g9c")
        graze = h_min < float(tol["g9c"]["graze_max_depth_m"])
        g_msgs = [f"g9c grazing: h_min = {h_min:.4f} m "
                  f"(must stay wet and < {tol['g9c']['graze_max_depth_m']} m) "
                  f"{'ok' if graze and h_min > 0 else 'FAIL'}"]
        ok &= report("c", s_ok and p_ok and graze and h_min > 0,
                     s_msgs + p_msgs + g_msgs)

    print("g9 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_g10(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g9-wind", args.work)
    with open(HERE / "tolerances" / "g9-wind.yaml", encoding="utf-8") as handle:
        tol = yaml.safe_load(handle)
    config = case_dir / "g10-seiche.yaml"
    if not validate_config(args.frehg, config):
        print("g10 gate: capability absent — the schema rejects the Q6 wind "
              "law/cap/u10/v10 keys (expected until the Q6 capability lands)")
        print("g10 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        tm, vm = monitor_series(handle, "eta_east")
    period, amps = gate_wind.seiche_metrics(tm, vm, float(tol["g10"]["release_time_s"]))
    merian = gate_wind.merian_period(5000.0, 10.0)
    courant = float(np.sqrt(gate_wind.GRAVITY * 10.0) * 5.0 / 100.0)
    ok, msgs = gate_wind.check_seiche(period, merian,
                                      float(tol["g10"]["period_tol_rel"]),
                                      amps, courant)
    for msg in msgs:
        print(f"  {msg}")
    print("g10 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_wind_orient(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g9-wind", args.work)
    with open(HERE / "tolerances" / "g9-wind.yaml", encoding="utf-8") as handle:
        tol = yaml.safe_load(handle)
    base = case_dir / "orient-base.yaml"
    if not validate_config(args.frehg, base):
        print("wind-orient battery: capability absent — the schema rejects "
              "the Q6 wind law/cap/u10/v10 keys (expected until the Q6 "
              "capability lands)")
        print("wind-orient battery: FAIL")
        return 1
    u_mag = 10.0
    fields = {}
    for name, (ux, uy, _) in gate_wind.ORIENTATIONS.items():
        config = case_dir / f"orient-{name}.yaml"
        with open(base, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        doc["simulation"]["id"] = f"g9-orient-{name}"
        doc["surface_water"]["wind"]["u10"] = {"constant": ux * u_mag}
        doc["surface_water"]["wind"]["v10"] = {"constant": uy * u_mag}
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            fields[name] = surface_eta(handle, 6000, 20, 20)
        (case_dir / "out" / "output.h5").unlink()
    ok, msgs = gate_wind.check_orientations(fields, float(tol["orient"]["tol_abs_m"]))
    for msg in msgs:
        print(f"  {msg}")
    print("wind-orient battery:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# outflow-staircase — the plan §8.2 BC-matrix cells `Outflow x {W,S,E,N}`
# (amendment V2-A11): a transmissive outlet on any edge must drain to
# Manning normal depth and conserve volume, in both feeding directions.
# The west pair is validation/swe-outflow-staircase as committed; the
# south, east and north variants are generated here by transposing and/or
# mirroring the strip (the defect had one line per MINUS direction, and
# the plus directions anchor the correct behaviour). The criteria are the
# case's own check.py — absolute analytic bands (Manning normal depth,
# volume budget, clamp, first-flow latency), deliberately robust to the
# preserved legacy plus/minus edge-arithmetic differences that the §8.1
# batteries measure (docs/theory/symmetry-exemptions.md).
# ---------------------------------------------------------------------------

# name -> (transposed?, reversed bed?, coord transform on (x, y), outlet
# ravel index, neighbour ravel index); L = the strip length in metres.
STAIRCASE_VARIANTS = {
    "south": (True, False, lambda x, y, L: (y, x), 0, 1),
    "east": (False, True, lambda x, y, L: (L - x, y), 9, 8),
    "north": (True, True, lambda x, y, L: (y, L - x), 9, 8),
}


def staircase_variants(case_dir: Path) -> list[tuple[Path, int, int]]:
    """Generate the south/east/north twins of each west-edge staircase
    configuration. Returns (config, outlet, neighbour) tuples."""
    lines = (case_dir / "input" / "dem.dat").read_text().splitlines()
    values = lines[6].split()  # one data row after the 6-line header
    length = float(len(values))
    configs = []
    for variant, (transposed, reverse, op, outlet, neighbor) in \
            STAIRCASE_VARIANTS.items():
        ordered = list(reversed(values)) if reverse else values
        dem = case_dir / "input" / f"dem-{variant}.dat"
        with open(dem, "w", encoding="utf-8") as handle:
            if transposed:
                # first data row is j = 0 (input-data.md; not ESRI north-up)
                handle.write(f"ncols 1\nnrows {len(values)}\nxllcorner 0.0\n"
                             "yllcorner 0.0\ncellsize 1.0\nnodata_value -9999\n")
                handle.write("\n".join(ordered) + "\n")
            else:
                handle.write(f"ncols {len(values)}\nnrows 1\nxllcorner 0.0\n"
                             "yllcorner 0.0\ncellsize 1.0\nnodata_value -9999\n")
                handle.write(" ".join(ordered) + "\n")
        for name in ("upslope-fed", "outlet-fed"):
            with open(case_dir / f"{name}.yaml", encoding="utf-8") as handle:
                doc = yaml.safe_load(handle)
            doc["simulation"]["id"] += f"-{variant}"
            if transposed:
                doc["domain"]["nx"], doc["domain"]["ny"] = (doc["domain"]["ny"],
                                                            doc["domain"]["nx"])
            doc["domain"]["bottom_elevation"] = {"file": f"input/dem-{variant}.dat"}
            for bc in doc["boundary_conditions"]:
                bc["region"]["polygon"] = [
                    list(op(float(p[0]), float(p[1]), length))
                    for p in bc["region"]["polygon"]]
            monitor = {0: (0, 0), 9: ((9, 0) if transposed else (0, 9))}[outlet]
            doc["output"]["monitors"][0]["j"], doc["output"]["monitors"][0]["i"] = monitor
            doc["output"]["filename"] = f"out/{name}-{variant}.h5"
            config = case_dir / f"{name}-{variant}.yaml"
            with open(config, "w", encoding="utf-8") as handle:
                yaml.safe_dump(doc, handle, sort_keys=False)
            configs.append((config, outlet, neighbor))
    return configs


def gate_outflow_staircase(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "swe-outflow-staircase", args.work)
    runs = [(case_dir / "upslope-fed.yaml", 0, 1),
            (case_dir / "outlet-fed.yaml", 0, 1)]
    runs += staircase_variants(case_dir)
    for config, _, _ in runs:
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    spec = importlib.util.spec_from_file_location("staircase_check",
                                                  case_dir / "check.py")
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    fails = []
    for config, outlet, neighbor in runs:
        stem = config.stem
        edge = stem.rsplit("-", 1)[-1]
        edge = edge if edge in STAIRCASE_VARIANTS else "west"
        fails += checker.check(f"{stem} (outflow x {edge})",
                               f"out/{stem}.h5", outlet, neighbor)
    for msg in fails:
        print(f"  FAIL: {msg}")
    print("outflow-staircase x-gate (Outflow x {W,S,E,N}):",
          "PASS" if not fails else f"FAIL ({len(fails)})")
    return 0 if not fails else 1


def output_times(handle, group: str) -> list[int]:
    """Sorted integer-second output times of an HDF5 field group."""
    return sorted(int(k) for k in handle[group].keys())


def validate_config(frehg: Path, config: Path) -> bool:
    """True when `frehg --validate` accepts the configuration."""
    result = subprocess.run([str(frehg), "--validate", str(config)],
                            cwd=config.parent, capture_output=True, text=True,
                            check=False)
    return result.returncode == 0


def g4_closure_error(output: Path) -> float:
    """The full surface closure identity including the audited clamp column
    (v2 Q4: the evaprain clamp's signed volume is measured into 'clamped',
    so dV + outflow - rain + evap - bc - clamped closes to rounding)."""
    with h5py.File(output, "r") as handle:
        table = handle["/monitor/mass_audit"][:]
    volume, rain, evap, outflow, bc = (table[:, k] for k in range(1, 6))
    clamped = table[:, 6] if table.shape[1] > 6 else np.zeros_like(rain)
    return float((volume[-1] - volume[0]) + outflow[-1] - rain[-1] + evap[-1] -
                 bc[-1] - clamped[-1])


def g4_drawdown_errors(output: Path, eta0: float, rate: float):
    """(times, max-over-cells |eta - (eta0 - rate t)| per time)."""
    times, errs = [], []
    with h5py.File(output, "r") as handle:
        for t in output_times(handle, "/surface/eta"):
            field = handle[f"/surface/eta/{t}"][:]
            errs.append(float(np.max(np.abs(field - (eta0 - rate * t)))))
            times.append(t)
    return times, errs


def g4_concentration_series(output: Path):
    """(times, max-cell concentration, max-cell depth) per output, plus the
    per-step surf_mass column of the transport audit."""
    times, s_vals, depths = [], [], []
    with h5py.File(output, "r") as handle:
        for t in output_times(handle, "/transport/concentration_surface"):
            s = handle[f"/transport/concentration_surface/{t}"][:]
            d = handle[f"/surface/depth/{t}"][:]
            times.append(t)
            s_vals.append(float(np.max(s)))
            depths.append(float(np.max(d)))
        mass = handle["/monitor/transport_audit"][:, 1]
    return times, s_vals, depths, mass


def run_and_rename(args: argparse.Namespace, case_dir: Path, config: str) -> Path:
    """Run one g4 sub-case and move its output aside (the sub-cases share
    out/output.h5)."""
    run_case(args.frehg, case_dir / config, case_dir, args.mpiexec, args.ranks)
    src = case_dir / "out" / "output.h5"
    dst = case_dir / "out" / (Path(config).stem + ".h5")
    shutil.move(src, dst)
    return dst


def gate_g4(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g4-evap", args.work)
    with open(HERE / "tolerances" / "g4-evap.yaml", encoding="utf-8") as handle:
        tol = yaml.safe_load(handle)
    ok = True

    def report(sub: str, sub_ok: bool, msgs: list[str]) -> bool:
        for msg in msgs:
            print(f"  {msg}")
        print(f"g4({sub}):", "ok" if sub_ok else "FAIL")
        return sub_ok

    # (a) prescribed-rate drawdown against the closed form.
    out_a = run_and_rename(args, case_dir, "g4a-drawdown.yaml")
    rate = 1.0e-6
    total = rate * 50000.0
    times, errs = g4_drawdown_errors(out_a, 0.0, rate)
    sub_ok, msgs = gate_evap.check_drawdown(
        times, errs, 0.0, rate, total, tol["drawdown"]["tol_fraction_of_total"])
    # Closure via the audit identity (reference: total evaporated volume).
    error = g4_closure_error(out_a)
    c_ok = abs(error) <= tol["closure"]["tol_fraction"] * (total * 12.0)
    c_msgs = [f"closure: |{error:.4e}| m^3 of {total * 12.0:.3f} m^3 evaporated "
              f"{'ok' if c_ok else 'FAIL'}"]
    ok &= report("a", sub_ok and c_ok, msgs + c_msgs)

    # (b) evaporative concentration against s0 V0 / V(t).
    out_b = run_and_rename(args, case_dir, "g4b-concentration.yaml")
    times, s_vals, depths, mass = g4_concentration_series(out_b)
    sub_ok, msgs = gate_evap.check_concentration(
        times, s_vals, depths, 10.0, 0.1, tol["concentration"]["tol_rel"],
        mass, tol["concentration"]["tol_mass_rel"])
    ok &= report("b", sub_ok, msgs)

    # (c) bulk mode under constant met forcing: E computable offline.
    config_c = case_dir / "g4c-bulk.yaml"
    if not validate_config(args.frehg, config_c):
        ok &= report("c", False,
                     ["capability absent: the schema rejects the Q4 "
                      "atmosphere/bulk keys (expected until the Q4 "
                      "capability PRs land — plan §6.1)"])
    else:
        with open(config_c, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        atm = doc["atmosphere"]
        e_bulk = gate_evap.bulk_evaporation_rate(
            atm["surface_temperature"]["constant"], atm["pressure"]["constant"],
            atm["wind_speed"]["constant"], atm["specific_humidity"]["constant"])
        print(f"  offline bulk rate E = {e_bulk:.4e} m/s")
        out_c = run_and_rename(args, case_dir, "g4c-bulk.yaml")
        t_end = doc["time"]["t_end"]
        times, errs = g4_drawdown_errors(out_c, 0.0, e_bulk)
        d_ok, d_msgs = gate_evap.check_drawdown(
            times, errs, 0.0, e_bulk, e_bulk * t_end,
            tol["drawdown"]["tol_fraction_of_total"])
        times, s_vals, depths, mass = g4_concentration_series(out_c)
        s_ok, s_msgs = gate_evap.check_concentration(
            times, s_vals, depths, 10.0, 0.1, tol["concentration"]["tol_rel"],
            mass, tol["concentration"]["tol_mass_rel"])
        ok &= report("c", d_ok and s_ok, d_msgs + s_msgs)

    # (d) drain-to-dry: positivity + the audited shortfall.
    out_d = run_and_rename(args, case_dir, "g4d-drain.yaml")
    with h5py.File(out_d, "r") as handle:
        d_times = output_times(handle, "/surface/depth")
        depth_min = [float(np.min(handle[f"/surface/depth/{t}"][:])) for t in d_times]
        depth_final = float(np.max(handle[f"/surface/depth/{d_times[-1]}"][:]))
    p_ok, p_msgs = gate_evap.check_positivity_and_dry(depth_min, depth_final, 1.0e-8)
    error = g4_closure_error(out_d)
    v0 = 0.02 * 12.0  # initial volume [m^3]
    c_ok = abs(error) <= tol["closure"]["tol_fraction"] * v0
    if c_ok:
        c_msgs = [f"closure: |{error:.4e}| m^3 of the {v0:.2f} m^3 initial volume ok"]
    else:
        c_msgs = [f"closure: |{error:.4e}| m^3 of the {v0:.2f} m^3 initial volume "
                  "FAIL (pre-Q4: the evaporation ledger books the potential "
                  "rate after dry-out)"]
    ok &= report("d", p_ok and c_ok, p_msgs + c_msgs)

    print("g4 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def read_reference_csv(path: Path) -> dict[str, np.ndarray]:
    """Column-name -> array from a digitized-reference CSV (comments '#';
    empty cells -> NaN)."""
    with open(path, encoding="utf-8") as handle:
        rows = [line.strip() for line in handle
                if line.strip() and not line.startswith("#")]
    header = rows[0].split(",")
    data = np.array([[float(v) if v else np.nan for v in line.split(",")]
                     for line in rows[1:]])
    return {name: data[:, k] for k, name in enumerate(header)}


def g5_mean_profile(handle, group: str, t: int, nz: int) -> np.ndarray:
    """Horizontal (whole-domain) average of a 3D field at time t -> (nz,)."""
    field = handle[f"{group}/{t}"][:].reshape(-1, nz)
    return field.mean(axis=0)


def transpose_g5_case(case_dir: Path) -> None:
    """Rewrite the staged g5 YAMLs as the y-z slice (the §8.1/§9 'transposed
    g5 slice' x-battery row): swap nx/ny (square cells) and every polygon's
    (x, y). The gate metrics are orientation-agnostic; finger positions are
    instability-set and deliberately not compared (plan §3.4)."""
    for config in sorted(case_dir.glob("*.yaml")):
        with open(config, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        dom = doc["domain"]
        dom["nx"], dom["ny"] = dom["ny"], dom["nx"]
        dom["dx"], dom["dy"] = dom["dy"], dom["dx"]

        def swap(polygon):
            return [[pt[1], pt[0]] for pt in polygon]

        if "evaporation" in doc.get("groundwater", {}):
            region = doc["groundwater"]["evaporation"]["region"]
            region["polygon"] = swap(region["polygon"])
        for bc in doc.get("boundary_conditions", []):
            bc["region"]["polygon"] = swap(bc["region"]["polygon"])
        doc["simulation"]["id"] += "-transposed"
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        print(f"  transposed: {config.name}")


def gate_g5(args: argparse.Namespace) -> int:
    reuse = args.reuse_output and (args.work / "g5-geng2015" / "out").exists()
    if reuse:
        case_dir = args.work / "g5-geng2015"
        print("g5 gate: --reuse-output — re-applying the metrics to the "
              "existing runs")
    else:
        case_dir = stage_case(args.repo, "g5-geng2015", args.work)
        if args.transposed:
            transpose_g5_case(case_dir)
    with open(HERE / "tolerances" / "g5-geng2015.yaml", encoding="utf-8") as handle:
        tol = yaml.safe_load(handle)
    config = case_dir / "g5-geng2015.yaml"
    if not validate_config(args.frehg, config):
        print("g5 gate: capability absent — the schema rejects the Q4 "
              "atmosphere/evaporation keys (expected until the Q4 capability "
              "PRs land — plan §6.1)")
        print("g5 gate: FAIL")
        return 1

    ok = True
    outputs = {}
    for name in ("g5-geng2015", "g5-geng2015-nodensity"):
        outputs[name] = case_dir / "out" / f"{name}.h5"
        if reuse and outputs[name].exists():
            continue
        run_case(args.frehg, case_dir / f"{name}.yaml", case_dir, args.mpiexec,
                 args.ranks)
        shutil.move(case_dir / "out" / "output.h5", outputs[name])

    ref_dir = args.repo / "benchmarks" / "g5-geng2015" / "reference"
    fig3 = read_reference_csv(ref_dir / "fig3_evaporation_rate.csv")
    fig4 = read_reference_csv(ref_dir / "fig4_moisture_ratio.csv")
    fig9 = read_reference_csv(ref_dir / "fig9b_salinity.csv")
    porosity = float(tol["porosity"])
    zone_area = float(tol["evaporation_zone_area_m2"])
    z_ref = fig4["elevation_m"]

    with open(case_dir / "g5-geng2015.yaml", encoding="utf-8") as handle:
        nz = int(yaml.safe_load(handle)["domain"]["nz"])

    with h5py.File(outputs["g5-geng2015"], "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
        if audit.shape[1] < 8:
            print("g5 gate: capability absent — gw_mass_audit has no "
                  "evaporation column (Q4 audit extension pending)")
            print("g5 gate: FAIL")
            return 1
        t_s = audit[:, 0]
        cumulative = audit[:, 7]  # cumulative evaporated volume [m^3]
        rate = np.diff(cumulative) / np.diff(t_s) / zone_area
        t_mid_h = 0.5 * (t_s[1:] + t_s[:-1]) / 3600.0

        # frehg's subsurface datum is the land surface at 0 (cells at
        # negative z, the b2 convention); the paper's elevations put the
        # surface at 2.0 m. 3D datasets are stored flat in the §7
        # (j nx + i) nz + k order, so reshape(-1, nz) yields one row per
        # column.
        z_model = 2.0 + handle["/groundwater/zcell/0"][:].reshape(-1, nz)[0]
        theta20 = g5_mean_profile(handle, "/groundwater/water_content", 72000, nz)
        theta50 = g5_mean_profile(handle, "/groundwater/water_content", 180000, nz)
        salt20 = g5_mean_profile(handle, "/transport/concentration", 72000, nz)
        salt50 = g5_mean_profile(handle, "/transport/concentration", 180000, nz)
        # Extrapolated surface saturation over the evaporation zone (the
        # V2-A15 anchor): 1.5 theta0 - 0.5 theta1 at the surface face,
        # columns x in [1, 49] m -> global i in [10, 490).
        wc50 = handle["/groundwater/water_content/180000"][:].reshape(-1, nz)
        w_face = 1.5 * wc50[10:490, 0] - 0.5 * wc50[10:490, 1]
        s_surf_50h = float(np.mean(np.clip(w_face, 0.0, porosity))) / porosity
        salt_mass = handle["/monitor/transport_audit"][:, 2]
        s50_full = handle["/transport/concentration/180000"][:].reshape(-1, nz)

    # (i) rate criteria (V2-A13).
    sub_ok, msgs = gate_evap.check_rate_series(
        t_mid_h, rate, float(tol["rate"]["e0_closed_form"]),
        float(tol["rate"]["e10_reference"]), s_surf_50h,
        decay_bound=float(tol["rate"]["decay_bound"]),
        saturation_band=(float(tol["rate"]["saturation_min"]),
                         float(tol["rate"]["saturation_max"])))
    for msg in msgs:
        print(f"  {msg}")
    ok &= sub_ok

    # (ii) profiles vs the digitized references (model interpolated onto
    # the reference elevation grid; zcell descends with k, so re-sort).
    order = np.argsort(z_model)

    def on_ref(profile: np.ndarray) -> np.ndarray:
        return np.interp(z_ref, z_model[order], profile[order])

    sub_ok, msgs = gate_evap.record_moisture_profiles(
        z_ref, on_ref(theta20) / porosity, on_ref(theta50) / porosity,
        fig4["moisture_20h"], fig4["moisture_50h"])
    for msg in msgs:
        print(f"  {msg}")
    ok &= sub_ok
    sub_ok, msgs = gate_evap.check_salinity_profiles(
        z_ref, on_ref(salt20), on_ref(salt50),
        fig9["salinity_20h_gL"], fig9["salinity_50h_gL"],
        float(tol["salinity"]["peak_min_g_l"]),
        float(tol["salinity"]["peak_above_z"]))
    for msg in msgs:
        print(f"  {msg}")
    ok &= sub_ok

    # (iii) salt-mass conservation (the paper's own 4 % bound).
    sub_ok, msgs = gate_evap.check_salt_mass(salt_mass,
                                             float(tol["salt_mass"]["tol_fraction"]))
    for msg in msgs:
        print(f"  {msg}")
    ok &= sub_ok

    # (iv) density contrast at 50 h (V2-A14 direction).
    with h5py.File(outputs["g5-geng2015-nodensity"], "r") as handle:
        s50_control = handle["/transport/concentration/180000"][:].reshape(-1, nz)
    edge = float(tol["density"]["edge_g_l"])
    zmin_beta = gate_evap.plume_edge_min_elevation(z_model, s50_full, edge)
    zmin_ctrl = gate_evap.plume_edge_min_elevation(z_model, s50_control, edge)
    sub_ok, msgs = gate_evap.check_density_contrast(
        zmin_beta, zmin_ctrl, float(tol["density"]["min_margin_m"]))
    for msg in msgs:
        print(f"  {msg}")
    ok &= sub_ok

    print("g5 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_b1_restart(args: argparse.Namespace) -> int:
    # Uninterrupted run with a checkpoint at half time.
    case_dir = stage_case(args.repo, "b1-sw", args.work)
    config = case_dir / "b1-sw.yaml"

    def add_checkpoint(doc):
        doc["output"]["checkpoint"] = {"interval": 9000}

    rewrite_config(config, add_checkpoint)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"

    # Restarted run from t = 9000 into a separate output file.
    restart_dir = stage_case(args.repo, "b1-sw", args.work / "restart")
    restart_config = restart_dir / "b1-sw.yaml"

    def add_restart(doc):
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()), "time": 9000}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)

    ok = True
    # Bitwise-class determinism (1e-12) on the default solver and on
    # aijkokkos (measured exactly 0.0 at one thread). Under an amg/gamg
    # override the straight and restarted runs rebuild their reused
    # hierarchies on different solve counters, so the fields carry a
    # deterministic offset at the 1e-8 solve-tolerance scale (measured
    # 6.9e-8 max on b1, Q7): the g1/AMG+restart lane bound is 1e-6 —
    # 14x the measured offset and far below any restart-state corruption.
    limit = 1.0e-6 if SOLVER_OVERRIDE in ("amg", "gamg") else 1.0e-12
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for var in ["eta", "depth", "uu", "vv"]:
            for t in ["10800", "14400", "18000"]:
                a = full[f"/surface/{var}/{t}"][:]
                b = part[f"/surface/{var}/{t}"][:]
                scale = max(float(np.abs(a).max()), 1.0e-12)
                diff = float(np.abs(a - b).max()) / scale
                print(f"  restart {var}@{t}: max rel diff = {diff:.3e} "
                      f"(allowed {limit:.0e})")
                if diff > limit:
                    ok = False
    print("b1 restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rank_invariance(args: argparse.Namespace) -> int:
    strict = args.mode == "strict"
    extra = STRICT_PETSC_OPTIONS if strict else []
    # Strict mode proves rank invariance of assembly and physics (plan
    # §8.2): rank-invariant preconditioning at machine-precision tolerances
    # leaves only rounding, 1e-12. The default bjacobi/icc lane bounds the
    # preconditioner-induced drift at production tolerances: measured
    # 5.4e-5 over b1's 3600 steps (per-solve rtol-level differences echoed
    # through wet/dry threshold crossings) — plan amendment A5 sets the
    # bound to 1e-4 from that data.
    limit = 1.0e-12 if strict else 1.0e-4
    outputs = {}
    for ranks in (1, 2, 4):
        case_dir = stage_case(args.repo, "b1-sw", args.work / f"n{ranks}")
        run_case(args.frehg, case_dir / "b1-sw.yaml", case_dir, args.mpiexec, ranks, extra)
        outputs[ranks] = case_dir / "out" / "output.h5"

    ok = True
    with h5py.File(outputs[1], "r") as base:
        for ranks in (2, 4):
            with h5py.File(outputs[ranks], "r") as other:
                worst = 0.0
                for var in ["eta", "depth", "uu", "vv"]:
                    for t in base[f"/surface/{var}"]:
                        a = base[f"/surface/{var}/{t}"][:]
                        b = other[f"/surface/{var}/{t}"][:]
                        scale = max(float(np.abs(a).max()), 1.0e-12)
                        worst = max(worst, float(np.abs(a - b).max()) / scale)
                print(f"  n=1 vs n={ranks} [{args.mode}]: max rel diff = {worst:.3e} "
                      f"(allowed {limit:.0e})")
                if worst > limit:
                    ok = False
    print(f"b1 rank invariance ({args.mode}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# v2 Q5 gate implementations (plan §4.2 + V2-A17).
# ---------------------------------------------------------------------------

def heat_tolerances() -> dict:
    with open(HERE / "tolerances" / "g6-g8-heat.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def subsurface_column(handle, group: str, t: int, nz: int) -> np.ndarray:
    """A (nz,) column from a flat nx = ny = 1 3D dataset at time t."""
    return handle[f"{group}/{t}"][:].reshape(nz)


def subsurface_slice(handle, group: str, t: int, nx: int, nz: int) -> np.ndarray:
    """A (nx, nz) x-z slice from a flat ny = 1 3D dataset at time t."""
    return handle[f"{group}/{t}"][:].reshape(nx, nz)


def measured_qz_down(handle, t: int, nz: int) -> float:
    """Downward-positive Darcy flux [m/s], mid-column, from the per-area
    qz output (positive up) of a single-column run."""
    qz = subsurface_column(handle, "/groundwater/qz", t, nz)
    return -float(qz[nz // 2])


def graded_depths(dz0: float, stretch: float, nz: int) -> np.ndarray:
    """Cell-centre depths below the surface for the dz dz_stretch^k mesh."""
    dz = dz0 * stretch ** np.arange(nz)
    tops = np.concatenate([[0.0], np.cumsum(dz)[:-1]])
    return tops + 0.5 * dz


def gate_g6(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g6-heat", args.work)
    tol = heat_tolerances()
    ok = True

    def report(sub: str, sub_ok: bool, msgs: list[str]) -> bool:
        for msg in msgs:
            print(f"  {msg}")
        print(f"g6({sub}):", "ok" if sub_ok else "FAIL")
        return sub_ok

    for config in ("g6a-bp.yaml", "g6b-ogata.yaml", "g6c-stallman.yaml",
                   "g6c-coupled.yaml"):
        if not validate_config(args.frehg, case_dir / config):
            ok &= report(config, False, [CAPABILITY_ABSENT_Q5])
    if not ok:
        print("g6 gate:", "FAIL")
        return 1

    # (a) B&P steady Peclet sweep. Both column ends are pinned CELLS, so
    # the isothermal planes sit at the pinned centres: L_eff = L - dz.
    rc_w, rc_b = 4.184e6, 2.8e6
    lam = 2.0
    alpha = lam / rc_b
    nz, dz, length, k_soil = 100, 0.1, 10.0, 1.0e-5
    l_eff = length - dz
    a_ok = True
    a_msgs: list[str] = []
    for pe in (-5.0, -1.0, 0.0, 1.0, 5.0):
        v_target = pe * alpha / l_eff
        q_target = v_target * rc_b / rc_w          # downward-positive [m/s]
        head_bottom = length - q_target * length / k_soil
        config = case_dir / f"g6a-pe{pe:+g}.yaml"
        with open(case_dir / "g6a-bp.yaml", encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        doc["simulation"]["id"] = f"g6a-pe{pe:+g}"
        for bc in doc["boundary_conditions"]:
            if bc["name"] == "head-bottom":
                bc["value"] = {"constant": float(head_bottom)}
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        t_end = int(doc["time"]["t_end"])
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            q_meas = measured_qz_down(handle, t_end, nz)
            temp = subsurface_column(handle, "/temperature/temperature", t_end, nz)
        (case_dir / "out" / "output.h5").unlink()
        f_ok, f_msgs = gate_heat.check_flux_window(
            q_meas, q_target, float(tol["g6a"]["flux_window_rel"]),
            float(tol["g6a"]["flux_floor_abs"]), f"g6a Pe={pe:+g}")
        a_msgs += f_msgs
        a_ok &= f_ok
        pe_meas = (q_meas * rc_w / rc_b) * l_eff / alpha
        depths = (np.arange(nz) + 0.5) * dz
        xi = (depths - 0.5 * dz) / l_eff
        t_ref = 20.0 + (10.0 - 20.0) * gate_heat.bp_profile(xi, pe_meas)
        err = float(np.max(np.abs(temp - t_ref)))
        p_ok, p_msgs = gate_heat.check_bp([pe_meas], [err], 10.0,
                                          float(tol["g6a"]["profile_tol_fraction"]))
        a_msgs += p_msgs
        a_ok &= p_ok
    ok &= report("a", a_ok, a_msgs)

    # (b) Ogata-Banks heat breakthrough on the graded vertical column.
    rc_eff = 2.0e6
    lam, k_soil = 2.2, 1.0e-5
    alpha = lam / rc_eff
    v_target = 1.5e-6
    q_target = v_target * rc_eff / rc_w
    nz, dz0, stretch = 145, 0.125, 1.02
    depths = graded_depths(dz0, stretch, nz)
    column_depth = float(dz0 * (stretch ** nz - 1.0) / (stretch - 1.0))
    b_ok = True
    b_msgs: list[str] = []
    results = {}
    for scheme in ("superbee", "upwind"):
        config = case_dir / f"g6b-{scheme}.yaml"
        with open(case_dir / "g6b-ogata.yaml", encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        doc["simulation"]["id"] = f"g6b-{scheme}"
        doc["temperature"]["scheme"]["advection"] = scheme
        for bc in doc["boundary_conditions"]:
            if bc["name"] == "head-bottom":
                bc["value"] = {"constant": float(column_depth -
                                                 q_target * column_depth / k_soil)}
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            times = output_times(handle, "/temperature/temperature")
            q_meas = measured_qz_down(handle, times[-1], nz)
            series = {t: subsurface_column(handle, "/temperature/temperature", t, nz)
                      for t in times}
        (case_dir / "out" / "output.h5").unlink()
        results[scheme] = (q_meas, times, series)
    q_meas = results["superbee"][0]
    f_ok, f_msgs = gate_heat.check_flux_window(
        q_meas, q_target, float(tol["g6b"]["flux_window_rel"]), 0.0, "g6b")
    b_msgs += f_msgs
    b_ok &= f_ok
    v_meas = q_meas * rc_w / rc_eff
    obs_targets = [1.0, 5.0, 10.0, 20.0, 50.0]
    obs_cells = [int(np.argmin(np.abs((depths - depths[0]) - x))) for x in obs_targets]
    for scheme in ("superbee", "upwind"):
        _, times, series = results[scheme]
        rel_errors = []
        # Gate from t >= 2e6 s: the pinned inlet CELL carries a half-cell
        # boundary-position ambiguity that dominates only while
        # sqrt(4 alpha t) spans a few cell widths (case README/comment).
        t_min = 2.0e6
        for cell in obs_cells:
            x = float(depths[cell] - depths[0])
            worst = 0.0
            for t in times:
                if t < t_min:
                    continue
                ref = gate_heat.ogata_banks(x, float(t), v_meas, alpha, 300.0, 330.0)
                worst = max(worst, abs(float(series[t][cell]) - ref) / 30.0)
            rel_errors.append(worst)
        points = [float(depths[c] - depths[0]) for c in obs_cells]
        c_ok, c_msgs = gate_heat.check_breakthrough(
            points, rel_errors, float(tol["g6b"]["breakthrough_tol_rel"]))
        if scheme == "superbee":
            b_msgs += [f"superbee (gated): {m}" for m in c_msgs]
            b_ok &= c_ok
        else:
            b_msgs += [f"upwind (recorded, not gated): {m}" for m in c_msgs]
    ok &= report("b", b_ok, b_msgs)

    # (c) Stallman diel damping: q_z in {0, +2e-6, -2e-6} m/s.
    rc_b = 2.96e6
    lam, k_soil = 2.0, 1.0e-4
    alpha = lam / rc_b
    nz, dz, length = 200, 0.01, 2.0
    period = 86400.0
    spinup = float(tol["g6c"]["spinup_cycles"]) * period
    gate_cells = [8, 18, 28]
    c_ok = True
    c_msgs: list[str] = []
    for q_target in (0.0, 2.0e-6, -2.0e-6):
        config = case_dir / f"g6c-q{q_target:+g}.yaml"
        with open(case_dir / "g6c-stallman.yaml", encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        doc["simulation"]["id"] = f"g6c-q{q_target:+g}"
        for bc in doc["boundary_conditions"]:
            if bc["name"] == "head-bottom":
                bc["value"] = {"constant": float(length - q_target * length / k_soil)}
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            times = np.array(output_times(handle, "/temperature/temperature"), dtype=float)
            q_meas = measured_qz_down(handle, int(times[-1]), nz)
            history = np.stack([
                subsurface_column(handle, "/temperature/temperature", int(t), nz)
                for t in times])
        (case_dir / "out" / "output.h5").unlink()
        f_ok, f_msgs = gate_heat.check_flux_window(
            q_meas, q_target, float(tol["g6c"]["flux_window_rel"]),
            float(tol["g6c"]["flux_floor_abs"]), f"g6c q={q_target:+g}")
        c_msgs += f_msgs
        c_ok &= f_ok
        sel = times >= spinup
        v_meas = q_meas * rc_w / rc_b
        amp_top, lag_top, _ = gate_heat.fit_sinusoid(times[sel], history[sel, 0], period)
        depths_z, amp_m, amp_r, lag_m, lag_r = [], [], [], [], []
        for cell in gate_cells:
            z_eff = (cell - 0) * dz
            amp, lag, _ = gate_heat.fit_sinusoid(times[sel], history[sel, cell], period)
            ratio_ref, lag_ref = gate_heat.stallman_amplitude_phase(
                z_eff, v_meas, alpha, period)
            depths_z.append(z_eff)
            amp_m.append(amp / amp_top)
            amp_r.append(ratio_ref)
            lag_m.append(lag - lag_top)
            lag_r.append(lag_ref)
        s_ok, s_msgs = gate_heat.check_stallman(
            depths_z, amp_m, amp_r, lag_m, lag_r,
            float(tol["g6c"]["amplitude_tol_rel"]), float(tol["g6c"]["phase_tol_s"]),
            f"g6c q={q_target:+g}")
        c_msgs += s_msgs
        c_ok &= s_ok
    ok &= report("c", c_ok, c_msgs)

    # (c, coupled) — the composed surface-fed variant (plan §4.2 note).
    config = case_dir / "g6c-coupled.yaml"
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    q_target = 2.0e-6
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        times = np.array(output_times(handle, "/temperature/temperature"), dtype=float)
        q_meas = measured_qz_down(handle, int(times[-1]), nz)
        history = np.stack([
            subsurface_column(handle, "/temperature/temperature", int(t), nz)
            for t in times])
    d_ok = True
    d_msgs: list[str] = []
    f_ok, f_msgs = gate_heat.check_flux_window(
        q_meas, q_target, float(tol["g6c_coupled"]["flux_window_rel"]), 0.0,
        "g6c coupled")
    d_msgs += f_msgs
    d_ok &= f_ok
    sel = times >= float(tol["g6c_coupled"]["spinup_cycles"]) * period
    v_meas = q_meas * rc_w / rc_b
    amp_top, lag_top, _ = gate_heat.fit_sinusoid(times[sel], history[sel, 0], period)
    depths_z, amp_m, amp_r, lag_m, lag_r = [], [], [], [], []
    for cell in gate_cells:
        z_eff = cell * dz
        amp, lag, _ = gate_heat.fit_sinusoid(times[sel], history[sel, cell], period)
        ratio_ref, lag_ref = gate_heat.stallman_amplitude_phase(z_eff, v_meas, alpha, period)
        depths_z.append(z_eff)
        amp_m.append(amp / amp_top)
        amp_r.append(ratio_ref)
        lag_m.append(lag - lag_top)
        lag_r.append(lag_ref)
    s_ok, s_msgs = gate_heat.check_stallman(
        depths_z, amp_m, amp_r, lag_m, lag_r,
        float(tol["g6c_coupled"]["amplitude_tol_rel"]),
        float(tol["g6c_coupled"]["phase_tol_s"]), "g6c coupled")
    d_msgs += s_msgs
    d_ok &= s_ok
    ok &= report("c coupled", d_ok, d_msgs)

    print("g6 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def temperature_audit_closure(output: Path) -> np.ndarray:
    """Per-step closure residuals of the temperature heat ledger,
    normalized by the basin heat content: d(surf_heat) minus the summed
    per-step deltas of the cumulative source columns."""
    with h5py.File(output, "r") as handle:
        table = handle["/monitor/temperature_audit"][:]
    surf = table[:, 1]
    cumulative = (table[:, 3] + table[:, 4] + table[:, 5] + table[:, 6] +
                  table[:, 7] + table[:, 8])
    residual = np.diff(surf) - np.diff(cumulative)
    scale = max(float(np.max(np.abs(surf))), 1.0e-12)
    return residual / scale


def gate_g7(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g7-heat", args.work)
    tol = heat_tolerances()
    ok = True

    def report(sub: str, sub_ok: bool, msgs: list[str]) -> bool:
        for msg in msgs:
            print(f"  {msg}")
        print(f"g7({sub}):", "ok" if sub_ok else "FAIL")
        return sub_ok

    for config in ("g7a-edinger.yaml", "g7a2-bulk.yaml", "g7b-channel.yaml"):
        if not validate_config(args.frehg, case_dir / config):
            ok &= report(config, False, [CAPABILITY_ABSENT_Q5])
    if not ok:
        print("g7 gate:", "FAIL")
        return 1

    # (a) stage 1: equilibrium relaxation + the per-step energy ledger.
    out_a = run_and_rename(args, case_dir, "g7a-edinger.yaml")
    with h5py.File(out_a, "r") as handle:
        tm, vm = monitor_series(handle, "t_probe")
    r_ok, r_msgs = gate_heat.check_relaxation(
        tm, vm, 30.0, 20.0, 30.0, 1.0, float(tol["g7a"]["relaxation_tol_rel"]))
    residuals = temperature_audit_closure(out_a)
    l_ok, l_msgs = gate_heat.check_energy_ledger(
        residuals, float(tol["g7a"]["ledger_tol_per_step"]))
    ok &= report("a", r_ok and l_ok, r_msgs + l_msgs)

    # (a) stage 2: full bulk-formula mode settles at the offline root.
    out_a2 = run_and_rename(args, case_dir, "g7a2-bulk.yaml")
    q_air = 0.6 * gate_heat.q_sat(25.0, 101.325)
    t_e = gate_heat.equilibrium_root(25.0, 101.325, 2.0, q_air, 150.0, 350.0)
    with h5py.File(out_a2, "r") as handle:
        tm, vm = monitor_series(handle, "t_probe")
    tail = vm[tm >= tm[-1] * 0.9]
    span = float(np.max(tail) - np.min(tail))
    s_ok = span <= float(tol["g7a2"]["steady_tol_c"])
    s_msgs = [f"g7a2 steadiness: trailing span {span:.2e} C "
              f"(allowed {tol['g7a2']['steady_tol_c']:g}) {'ok' if s_ok else 'FAIL'}"]
    e_ok, e_msgs = gate_heat.check_equilibrium_consistency(
        float(vm[-1]), t_e, float(tol["g7a2"]["equilibrium_tol_c"]))
    ok &= report("a2", s_ok and e_ok, s_msgs + e_msgs)

    # (b) channel thermal plume: steady profile + measured-inlet-convolved
    # transient breakthrough (the inlet cell is a 1000 s mixing volume, so
    # the ideal-step closed form is convolved with the measured inlet
    # history — the measured-flux principle applied to the inlet signal).
    config_b = case_dir / "g7b-channel.yaml"
    with open(config_b, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    doc["output"]["monitors"].append(
        {"name": "t_inlet", "i": 0, "j": 0, "variables": ["temperature_surface"]})
    with open(config_b, "w", encoding="utf-8") as handle:
        yaml.safe_dump(doc, handle, sort_keys=False)
    run_case(args.frehg, config_b, case_dir, args.mpiexec, args.ranks)
    nx, dx = 840, 100.0
    rho_c, depth_nom = 4.184e6, 1.0
    k_e, t_amb, d_l = 25.0, 15.0, 5.0
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        t_end = output_times(handle, "/temperature/temperature_surface")[-1]
        temp = handle[f"/temperature/temperature_surface/{t_end}"][:].reshape(nx)
        uu = handle[f"/surface/uu/{t_end}"][:].reshape(nx)
        eta = handle[f"/surface/eta/{t_end}"][:].reshape(nx)
        tm, vm = monitor_series(handle, "t_mid")
        ti, vi = monitor_series(handle, "t_inlet")
    depth = depth_nom + float(np.mean(eta[10:-10]))
    u_meas = float(np.mean(uu[10:-10]))
    v_ok = abs(u_meas - 0.5) <= float(tol["g7b"]["velocity_tol_rel"]) * 0.5
    b_msgs = [f"g7b velocity: measured u = {u_meas:.4f} m/s (target 0.5, "
              f"allowed {tol['g7b']['velocity_tol_rel']:.0%}) {'ok' if v_ok else 'FAIL'}"]
    decay = k_e / (rho_c * depth)
    x_c = (np.arange(nx) + 0.5) * dx
    excess0 = float(temp[0] - t_amb)
    steady_ref = t_amb + excess0 * gate_heat.channel_steady(
        x_c - x_c[0], u_meas, d_l, decay)
    c_ok, c_msgs = gate_heat.check_channel(
        x_c[:-2], temp[:-2], steady_ref[:-2], float(tol["g7b"]["steady_l2_tol"]),
        [0.0], [0.0], 1.0, excess0)
    # Transient: convolve the closed-form step response with the measured
    # inlet increments (superposition of the linear ADE).
    x_mid = float(x_c[420] - x_c[0])
    inlet = np.asarray(vi, dtype=float) - t_amb
    t_in = np.asarray(ti, dtype=float)
    ref_mid = []
    for t_now in tm:
        val = 0.0
        prev = 0.0
        for tj, sj in zip(t_in, inlet):
            if tj >= t_now:
                break
            step = sj - prev
            prev = sj
            if step != 0.0:
                val += step * gate_heat.channel_transient(
                    x_mid, float(t_now - tj), u_meas, d_l, decay)
        ref_mid.append(t_amb + val)
    err = float(np.max(np.abs(np.asarray(vm) - np.asarray(ref_mid)))) / abs(excess0)
    t_ok = err <= float(tol["g7b"]["transient_tol_rel"])
    b_msgs += c_msgs[:1]
    b_msgs.append(f"g7b transient (measured-inlet convolution): max "
                  f"{err:.2%} of the step (allowed "
                  f"{tol['g7b']['transient_tol_rel']:.0%}) {'ok' if t_ok else 'FAIL'}")
    ok &= report("b", v_ok and c_ok and t_ok, b_msgs)

    print("g7 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_g8(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g8-hrl", args.work)
    tol = heat_tolerances()
    if not validate_config(args.frehg, case_dir / "g8-hrl.yaml"):
        print(f"g8 gate: {CAPABILITY_ABSENT_Q5}")
        print("g8 gate: FAIL")
        return 1
    nx, nz, dx, dz = 40, 20, 0.05, 0.05
    width, h_eff, z_top_pin = 2.0, 0.95, 0.025
    lam, rc_w, beta_t, delta_t = 2.092, 4.184e6, 2.0e-4, 10.0
    alpha_e = lam / rc_w
    x_c = (np.arange(nx) + 0.5) * dx
    depths = (np.arange(nz) + 0.5) * dz
    histories = {}
    for label, k_soil in (("sub", 7.5e-3), ("super", 1.375e-2)):
        config = case_dir / f"g8-{label}.yaml"
        with open(case_dir / "g8-hrl.yaml", encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        doc["simulation"]["id"] = f"g8-{label}"
        for soil in doc["soil"]["types"]:
            soil["ksx"] = soil["ksy"] = soil["ksz"] = float(k_soil)
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        ra = gate_heat.rayleigh_number(k_soil, beta_t, delta_t, h_eff, alpha_e)
        print(f"  g8 {label}: K = {k_soil:g} -> Ra_eff = {ra:.2f} "
              f"(Ra_c = {gate_heat.RAYLEIGH_CRITICAL:.2f})")
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        amps, nus = [], []
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            for t in output_times(handle, "/temperature/temperature"):
                temp = subsurface_slice(handle, "/temperature/temperature", t, nx, nz)
                qz = subsurface_slice(handle, "/groundwater/qz", t, nx, nz)
                amps.append(gate_heat.hrl_mode_amplitude(
                    temp, x_c, depths, z_top_pin, h_eff, width))
                nus.append(gate_heat.hrl_nusselt(
                    temp, qz[:, nz // 2 - 1], nz // 2 - 1, dz, lam, rc_w,
                    delta_t, h_eff))
        (case_dir / "out" / "output.h5").unlink()
        histories[label] = (np.array(amps), np.array(nus))
    ok, msgs = gate_heat.check_hrl(histories["sub"][0], histories["super"][1],
                                   float(tol["g8"]["decay_factor"]),
                                   float(tol["g8"]["nusselt_min"]))
    msgs.append(f"recorded: subcritical trailing Nusselt "
                f"{float(histories['sub'][1][-1]):.4f} (conduction = 1); "
                f"supercritical mode amplitude end "
                f"{float(histories['super'][0][-1]):.3f} K (seed 0.5)")
    for msg in msgs:
        print(f"  {msg}")
    print("g8 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def heat_orient_ic(nx: int, ny: int, nz: int, dx: float, dy: float,
                   dz: float) -> np.ndarray:
    """The heat-battery base IC: 20 C + an off-axis warm blob, decaying
    with depth — injective enough that any orientation mix-up shows."""
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    depth = (np.arange(nz) + 0.5) * dz
    blob = 6.0 * np.exp(-(((x[None, :] - 1.1) ** 2 + (y[:, None] - 1.7) ** 2)
                          / 0.5 ** 2))
    profile = 1.0 - depth / (nz * dz)
    return 20.0 + blob[:, :, None] * profile[None, None, :]


def write_ic_raster(path: Path, field: np.ndarray) -> None:
    """(ny, nx, nz) -> the (j*nx + i)*nz + k flat file format."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("# heat orientation battery IC (generated by the harness)\n")
        for j in range(field.shape[0]):
            for i in range(field.shape[1]):
                for k in range(field.shape[2]):
                    handle.write(f"{field[j, i, k]:.10f}\n")


def gate_heat_orient(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "g6-heat", args.work)
    tol = heat_tolerances()
    base = case_dir / "orient-base.yaml"
    if not validate_config(args.frehg, base):
        print(f"heat-orient battery: {CAPABILITY_ABSENT_Q5}")
        print("heat-orient battery: FAIL")
        return 1
    nx = ny = 8
    nz = 6
    dx = dy = 0.5
    dz = 0.25
    length = nx * dx
    ic_base = heat_orient_ic(nx, ny, nz, dx, dy, dz)
    (case_dir / "input").mkdir(exist_ok=True)
    with open(base, encoding="utf-8") as handle:
        base_doc = yaml.safe_load(handle)
    t_fields = {}
    h_fields = {}
    # The §8.1 battery runs strict-mode (rank/orientation-invariant
    # preconditioning) so the tolerance bounds rounding, not solver noise.
    strict = STRICT_PETSC_OPTIONS_GW
    for name, (array_op, coord_op) in gate_heat.dihedral_table(length).items():
        ic = np.stack([array_op(ic_base[:, :, k]) for k in range(nz)], axis=2)
        write_ic_raster(case_dir / "input" / f"orient_ic_{name}.dat", ic)
        doc = yaml.safe_load(yaml.safe_dump(base_doc))
        doc["simulation"]["id"] = f"g6-orient-{name}"
        doc["initial_conditions"]["temperature"]["groundwater"] = {
            "file": f"input/orient_ic_{name}.dat"}
        for bc in doc["boundary_conditions"]:
            poly = [list(coord_op(float(p[0]), float(p[1]))) for p in
                    bc["region"]["polygon"]]
            bc["region"]["polygon"] = poly
        config = case_dir / f"orient-{name}.yaml"
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks, strict)
        t_out = int(base_doc["time"]["t_end"])
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            t_fields[name] = handle[f"/temperature/temperature/{t_out}"][:].reshape(
                ny, nx, nz)
            h_fields[name] = handle[f"/groundwater/hydraulic_head/{t_out}"][:].reshape(
                ny, nx, nz)
        (case_dir / "out" / "output.h5").unlink()
    ok_t, msgs_t = gate_heat.check_heat_orientations(
        t_fields, length, float(tol["orient"]["tol_abs"]), "temperature")
    ok_h, msgs_h = gate_heat.check_heat_orientations(
        h_fields, length, float(tol["orient"]["tol_abs"]), "head")
    for msg in msgs_t + msgs_h:
        print(f"  {msg}")
    ok = ok_t and ok_h
    print("heat-orient battery:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_heat_restart(args: argparse.Namespace) -> int:
    """Two-scalar restart determinism (v2 plan §10 risk register): the
    heat battery base case run straight vs checkpoint-restarted must
    agree bitwise-class (1e-12 rel) on temperature, head, and moisture."""
    case_dir = stage_case(args.repo, "g6-heat", args.work)
    base = case_dir / "orient-base.yaml"
    if not validate_config(args.frehg, base):
        print(f"heat-restart determinism: {CAPABILITY_ABSENT_Q5}")
        print("heat-restart determinism: FAIL")
        return 1
    nx = ny = 8
    nz = 6
    ic = heat_orient_ic(nx, ny, nz, 0.5, 0.5, 0.25)
    (case_dir / "input").mkdir(exist_ok=True)
    write_ic_raster(case_dir / "input" / "orient_ic_id.dat", ic)

    def use_ic(doc):
        doc["initial_conditions"]["temperature"]["groundwater"] = {
            "file": "input/orient_ic_id.dat"}
        doc["output"]["variables"]["groundwater"] = ["hydraulic_head", "water_content"]
        return doc

    def add_checkpoint(doc):
        use_ic(doc)
        doc["output"]["checkpoint"] = {"interval": 10000}
        doc["output"]["filename"] = "out/straight.h5"
        return doc

    def add_restart(doc):
        use_ic(doc)
        doc["restart"] = {"enabled": True, "file": "out/straight.h5", "time": 10000}
        doc["output"]["filename"] = "out/restarted.h5"
        return doc

    straight = case_dir / "straight.yaml"
    shutil.copy(base, straight)
    rewrite_config(straight, add_checkpoint)
    run_case(args.frehg, straight, case_dir, args.mpiexec, args.ranks)
    restarted = case_dir / "restarted.yaml"
    shutil.copy(base, restarted)
    rewrite_config(restarted, add_restart)
    run_case(args.frehg, restarted, case_dir, args.mpiexec, args.ranks)

    ok = True
    with h5py.File(case_dir / "out" / "straight.h5", "r") as a, \
            h5py.File(case_dir / "out" / "restarted.h5", "r") as b:
        for group in ("/temperature/temperature", "/groundwater/hydraulic_head",
                      "/groundwater/water_content"):
            t = 20000
            va = a[f"{group}/{t}"][:]
            vb = b[f"{group}/{t}"][:]
            scale = max(float(np.abs(va).max()), 1.0e-12)
            diff = float(np.abs(va - vb).max()) / scale
            print(f"  restart {group}@{t}: max rel diff = {diff:.3e}")
            if diff > 1.0e-12:
                ok = False
    print("heat-restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rank_invariance_heat(args: argparse.Namespace) -> int:
    """The §7.1 decomposition lane for the temperature module: the heat
    battery base case at 1/2/4 ranks, strict (rounding-only) and default
    (production-preconditioner drift) modes."""
    strict = args.mode == "strict"
    extra = STRICT_PETSC_OPTIONS_GW if strict else []
    limit = 1.0e-12 if strict else 1.0e-4
    nx = ny = 8
    nz = 6
    ic = heat_orient_ic(nx, ny, nz, 0.5, 0.5, 0.25)
    outputs = {}
    for ranks in (1, 2, 4):
        case_dir = stage_case(args.repo, "g6-heat", args.work / f"n{ranks}")
        if not validate_config(args.frehg, case_dir / "orient-base.yaml"):
            print(f"heat rank invariance: {CAPABILITY_ABSENT_Q5}")
            print(f"heat rank invariance ({args.mode}): FAIL")
            return 1
        (case_dir / "input").mkdir(exist_ok=True)
        write_ic_raster(case_dir / "input" / "orient_ic_id.dat", ic)

        def use_ic(doc):
            doc["initial_conditions"]["temperature"]["groundwater"] = {
                "file": "input/orient_ic_id.dat"}
            return doc

        rewrite_config(case_dir / "orient-base.yaml", use_ic)
        run_case(args.frehg, case_dir / "orient-base.yaml", case_dir, args.mpiexec,
                 ranks, extra)
        outputs[ranks] = case_dir / "out" / "output.h5"
    ok = True
    with h5py.File(outputs[1], "r") as base:
        for ranks in (2, 4):
            with h5py.File(outputs[ranks], "r") as other:
                worst = 0.0
                for group in ("/temperature/temperature", "/groundwater/hydraulic_head"):
                    for t in base[group]:
                        a = base[f"{group}/{t}"][:]
                        b = other[f"{group}/{t}"][:]
                        scale = max(float(np.abs(a).max()), 1.0e-12)
                        worst = max(worst, float(np.abs(a - b).max()) / scale)
                print(f"  n=1 vs n={ranks} [{args.mode}]: max rel diff = {worst:.3e} "
                      f"(allowed {limit:.0e})")
                if worst > limit:
                    ok = False
    print(f"heat rank invariance ({args.mode}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Q7 §8.1 backfill — the v1-module dihedral batteries (swe / gw / transport),
# on the g6-heat battery machinery (gate_heat.dihedral_table): the base
# configuration's polygons, rasters and IDs are pushed through the 8
# dihedral transforms, each variant is run in strict solver mode, and every
# listed field must map onto the base run's under the same transform to the
# strict tolerance. Base cases live in benchmarks/x-orient/. Waived legacy
# asymmetries are dormant in these bases BY DESIGN and documented in
# docs/theory/symmetry-exemptions.md — an over-tolerance deviation here is
# a bug, never exempt (plan §8.1).
# ---------------------------------------------------------------------------


def read_raster2d_file(path: Path) -> np.ndarray:
    """Flat-list 2D raster (j*nx + i, '#' comments, any number of
    whitespace-separated values per line — input-data.md) -> 1D values."""
    values = [float(tok) for line in path.read_text().splitlines()
              if line.strip() and not line.lstrip().startswith("#")
              for tok in line.split()]
    return np.asarray(values)


def write_raster2d_file(path: Path, field: np.ndarray) -> None:
    """(ny, nx) -> the flat-list j*nx + i format (input-data.md)."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("# x-orient battery raster (generated by the harness)\n")
        for j in range(field.shape[0]):
            for i in range(field.shape[1]):
                handle.write(f"{field[j, i]:.10f}\n")


def transform_polygons(doc: dict, coord_op) -> None:
    """Apply coord_op to every polygon in the configuration: all BC
    regions plus the optional rainfall/evaporation exclusion regions."""
    for bc in doc.get("boundary_conditions", []):
        bc["region"]["polygon"] = [list(coord_op(float(p[0]), float(p[1])))
                                   for p in bc["region"]["polygon"]]
    sw = doc.get("surface_water", {})
    for key in ("rainfall", "evaporation"):
        exclude = sw.get(key, {}).get("exclude")
        if exclude:
            exclude["polygon"] = [list(coord_op(float(p[0]), float(p[1])))
                                  for p in exclude["polygon"]]


def run_dihedral_battery(args: argparse.Namespace, label: str, base_name: str,
                         fields: dict[str, int], strict: list[str],
                         raster_keys: list[tuple] = (),
                         case: str = "x-orient",
                         tol_key: str | None = None,
                         tol_override: dict | None = None) -> int:
    """Shared §8.1 driver. `fields` maps HDF5 group paths to their nz;
    `raster_keys` lists (section, key) config paths whose 2D rasters are
    transformed alongside the polygons. Two tolerance classes per battery
    (tolerances/x-orient.yaml, or `tol_override` for Stage-3 batteries
    whose classes live in stage3-pd.yaml): `transpose` transforms preserve
    each domain edge's minus/plus class, so they gate scheme/indexing
    symmetry at machine precision; the six `mixed` transforms swap minus
    and plus edges, whose preserved legacy arithmetic differs by design —
    their bound is calibrated above the measured legacy floor and
    documented in docs/theory/symmetry-exemptions.md. Constant wind
    vectors are rotated with the domain (the coord_op's linear part), so
    a wind-driven base is a true dihedral identity."""
    case_dir = stage_case(args.repo, case, args.work / label)
    base = case_dir / base_name
    if tol_override is not None:
        tol = tol_override
    else:
        with open(HERE / "tolerances" / "x-orient.yaml",
                  encoding="utf-8") as handle:
            tol = yaml.safe_load(handle)["orient"][tol_key or label.split("-")[0]]
    with open(base, encoding="utf-8") as handle:
        base_doc = yaml.safe_load(handle)
    nx = int(base_doc["domain"]["nx"])
    ny = int(base_doc["domain"]["ny"])
    length = nx * float(base_doc["domain"]["dx"])
    assert nx == ny, "the dihedral table requires a square base"
    rasters = {}
    for path in raster_keys:
        node = base_doc
        for key in path:
            node = node[key]
        rasters[path] = read_raster2d_file(
            case_dir / node["file"]).reshape(ny, nx)
    results: dict[str, dict[str, np.ndarray]] = {group: {} for group in fields}
    for name, (array_op, coord_op) in gate_heat.dihedral_table(length).items():
        doc = yaml.safe_load(yaml.safe_dump(base_doc))
        doc["simulation"]["id"] = f"{label}-{name}"
        # Stage-3 T3 (W3.7): the axis-swapping transforms map periodic-x
        # onto periodic-y — the battery is what proves the implementation
        # axis-generic (plan §7 pd8).
        if name in ("transpose", "antitranspose", "rot90", "rot270"):
            per = doc.get("domain", {}).get("periodic")
            if per:
                doc["domain"]["periodic"] = {"x": bool(per.get("y", False)),
                                             "y": bool(per.get("x", False))}
        # Rotate a constant wind vector with the domain: the dihedral
        # coord_op is affine with an orthogonal linear part, so the image
        # of the vector is coord_op(u, v) - coord_op(0, 0).
        wind = doc.get("surface_water", {}).get("wind", {})
        if (wind.get("enabled")
                and isinstance(wind.get("u10"), dict)
                and "constant" in wind["u10"]):
            u = float(wind["u10"]["constant"])
            v = float(wind["v10"]["constant"])
            ox, oy = coord_op(0.0, 0.0)
            wx, wy = coord_op(u, v)
            wind["u10"] = {"constant": wx - ox}
            wind["v10"] = {"constant": wy - oy}
        for m, (path, field) in enumerate(rasters.items()):
            raster_path = f"input/orient_raster_{m}_{name}.dat"
            write_raster2d_file(case_dir / raster_path, array_op(field))
            node = doc
            for key in path[:-1]:
                node = node[key]
            node[path[-1]] = {"file": raster_path}
        transform_polygons(doc, coord_op)
        config = case_dir / f"{label}-{name}.yaml"
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks, strict)
        t_out = int(float(base_doc["time"]["t_end"]))
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            for group, nz in fields.items():
                results[group][name] = handle[f"{group}/{t_out}"][:].reshape(
                    ny, nx, nz)
        (case_dir / "out" / "output.h5").unlink()
    ok = True
    table = gate_heat.dihedral_table(length)
    for group, per_orient in results.items():
        base_field = per_orient["id"]
        short = group.rsplit("/", 1)[-1]
        for name, (array_op, _) in table.items():
            if name == "id":
                continue
            recorded = False
            if name == "transpose":
                bound = tol.get("per_field_transpose", {}).get(
                    group, tol["tol_abs_transpose"])
            else:
                bound = tol.get("per_field_mixed", {}).get(
                    group, tol["tol_abs_mixed"])
                recorded = bool(tol.get("mixed_recorded", False))
            pushed = np.stack([array_op(base_field[:, :, k])
                               for k in range(base_field.shape[2])], axis=2)
            dev = float(np.max(np.abs(pushed - per_orient[name])))
            good = recorded or dev <= bound
            print(f"  {short} orientation {name}: max |delta| = {dev:.3e} "
                  + (f"(RECORDED; legacy mixed class)" if recorded
                     else f"(allowed {bound:.1e})")
                  + ("" if good else " FAIL"))
            ok &= good
    print(f"{label} battery:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_swe_orient(args: argparse.Namespace) -> int:
    return run_dihedral_battery(
        args, "swe-orient", "swe-orient-base.yaml",
        {"/surface/eta": 1, "/surface/depth": 1},
        STRICT_PETSC_OPTIONS,
        raster_keys=[("domain", "bottom_elevation")])


def gate_gw_orient(args: argparse.Namespace) -> int:
    return run_dihedral_battery(
        args, "gw-orient", "gw-orient-base.yaml",
        {"/groundwater/hydraulic_head": 8, "/groundwater/water_content": 8},
        STRICT_PETSC_OPTIONS_GW)


def gate_transport_orient(args: argparse.Namespace) -> int:
    return run_dihedral_battery(
        args, "transport-orient", "transport-orient-base.yaml",
        {"/surface/eta": 1, "/transport/concentration_surface": 1},
        STRICT_PETSC_OPTIONS,
        raster_keys=[("initial_conditions", "surface", "eta"),
                     ("initial_conditions", "transport", "surface")])


def gate_surface_bc_sides(args: argparse.Namespace) -> int:
    """§8.2 side sweep for the surface footprint kinds: the conveyor case
    (single-cell discharge + riding scalar_value + eta outlet column,
    benchmarks/x-orient/conveyor-base.yaml) run in its four rotations, each
    judged by absolute criteria — closure, clamp, inflow plateau, mid-path
    front arrival, scalar bounds, outlet stage held. Three passes:
    uncoupled as committed, COUPLED (a soil column injected under the
    basin, seepage active — the coupled x side cells run by execution, not
    by same-code-path argument), and coupled with the eta outlet replaced
    by a transmissive `outflow` column (the coupled Outflow x side cells;
    the uncoupled ones belong to the analytic staircase gate)."""
    case_dir = stage_case(args.repo, "x-orient", args.work)
    with open(case_dir / "conveyor-base.yaml", encoding="utf-8") as handle:
        base_doc = yaml.safe_load(handle)
    length = int(base_doc["domain"]["nx"]) * float(base_doc["domain"]["dx"])
    d = float(base_doc["domain"]["dx"])
    table = gate_heat.dihedral_table(length)
    picks = {"west": "id", "south": "rot90", "east": "rot180", "north": "rot270"}
    outlet_sel = {"west": (slice(None), 11), "south": (11, slice(None)),
                  "east": (slice(None), 0), "north": (0, slice(None))}

    def cell_of(x: float, y: float) -> tuple[int, int]:
        return int(round(y / d - 0.5)), int(round(x / d - 0.5))

    def couple(doc):
        doc["modules"]["groundwater"] = True
        doc["coupling"] = {"mode": "sync"}
        doc["domain"]["nz"] = 4
        doc["domain"]["dz"] = 0.25
        doc["groundwater"] = {"specific_storage": 1.0e-5,
                              "timestep": {"dt_init": 0.5, "dt_min": 0.5,
                                           "dt_max": 0.5}}
        doc["soil"] = {"types": [{"name": "sand", "ksx": 1.0e-5,
                                  "ksy": 1.0e-5, "ksz": 1.0e-5,
                                  "theta_s": 0.35, "theta_r": 0.02,
                                  "vg_alpha": 3.0, "vg_n": 2.0}],
                       "map": {"constant": "sand"}}
        doc["initial_conditions"]["groundwater"] = {
            "water_table": {"constant": -0.7}}
        doc["initial_conditions"]["transport"]["groundwater"] = {"constant": 0.0}

    t_out = int(float(base_doc["time"]["t_end"]))
    # `outflow x coupled` is NOT a pass here: on this flat bed the
    # transmissive drop is degenerate (the analytic staircase gate owns the
    # uncoupled side sweep on its descending bed; the coupled exchange term
    # is side-agnostic and unit-covered by
    # CoupledModule tests + the swere-superslab validation record).
    passes = [("uncoupled", None), ("coupled", couple)]
    ok = True
    for pass_name, mutate in passes:
        coupled_pass = pass_name == "coupled"
        for side, transform in picks.items():
            coord_op = table[transform][1]
            doc = yaml.safe_load(yaml.safe_dump(base_doc))
            doc["simulation"]["id"] = f"conveyor-{pass_name}-{side}"
            if mutate is not None:
                mutate(doc)
            transform_polygons(doc, coord_op)
            doc["output"]["filename"] = f"out/conveyor-{pass_name}-{side}.h5"
            config = case_dir / f"conveyor-{pass_name}-{side}.yaml"
            with open(config, "w", encoding="utf-8") as handle:
                yaml.safe_dump(doc, handle, sort_keys=False)
            run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks,
                     STRICT_PETSC_OPTIONS)
            output = case_dir / f"out/conveyor-{pass_name}-{side}.h5"
            with h5py.File(output, "r") as handle:
                audit = handle["/monitor/mass_audit"]
                cols = {n: k for k, n in
                        enumerate(audit.attrs["columns"].decode().split(","))}
                row0, rowN = audit[0], audit[-1]
                seep = float(rowN[cols["seepage"]]) if "seepage" in cols else 0.0
                error = float(
                    (rowN[cols["volume"]] - row0[cols["volume"]])
                    + rowN[cols["boundary_outflow"]] - rowN[cols["rain"]]
                    + rowN[cols["evaporation"]] - rowN[cols["bc_inflow"]]
                    - rowN[cols["clamped"]]) - seep
                bc_in = float(rowN[cols["bc_inflow"]])
                clamped = float(rowN[cols["clamped"]])
                conc = handle[
                    f"/transport/concentration_surface/{t_out}"][:].reshape(12, 12)
                eta = handle[f"/surface/eta/{t_out}"][:].reshape(12, 12)
            inflow_cell = cell_of(*coord_op(1.0, 15.0))
            mid_cell = cell_of(*coord_op(13.0, 15.0))
            # The coupled pass carries two documented allowances, identical
            # on all four sides (side-invariance is the gated content): the
            # eta-held cells over live soil columns hold their stage to
            # ~3e-5 (seepage perturbs the held row within the solve
            # tolerance class), and the audit identity leaks ~1.3 % of the
            # injection through the eta-held/exchange bookkeeping — the
            # same unclosed eta-BC + hard-exchange audit term the b6 smoke
            # measures (recorded in dod-Q7; sync == subcycled there).
            closure_tol = (0.02 * max(bc_in, 1.0)) if coupled_pass \
                else 1.0e-6 * max(bc_in, 1.0)
            held_tol = 1.0e-4 if coupled_pass else 1.0e-9
            held = float(np.abs(eta[outlet_sel[side]]).max())
            checks = [
                ("closure", abs(error) <= closure_tol,
                 f"|{error:.3e}| m^3 vs {bc_in:.1f} injected"),
                ("clamp", clamped <= 1.0e-9, f"{clamped:.3e} m^3 minted"),
                ("inflow plateau", conc[inflow_cell] >= 0.999,
                 f"c{inflow_cell} = {conc[inflow_cell]:.6f}"),
                ("mid-path arrival", conc[mid_cell] >= 0.5,
                 f"c{mid_cell} = {conc[mid_cell]:.6f}"),
                ("bounds", -1.0e-9 <= float(conc.min()) and
                 float(conc.max()) <= 1.0 + 1.0e-9,
                 f"range [{conc.min():.3e}, {conc.max():.6f}]"),
                ("outlet stage held", held <= held_tol,
                 f"max |eta| = {held:.3e}"),
            ]
            for name, good, detail in checks:
                print(f"  {pass_name}/{side}: {name}: {detail}"
                      f"{'' if good else ' FAIL'}")
                ok &= good
    print("surface-bc-sides x-gate (discharge/scalar_value/eta x sides x "
          "coupling):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_wind_restart(args: argparse.Namespace) -> int:
    """§8.4 restart-mid-transient row for wind: the g10 seiche case
    checkpointed BEFORE the forcing step-down (wind cutoff at t = 30000 s),
    so the restarted leg re-samples the series across the step. Straight vs
    restarted must agree to the bitwise class (1e-12 rel)."""
    case_dir = stage_case(args.repo, "g9-wind", args.work)
    config = case_dir / "g10-seiche.yaml"
    if not validate_config(args.frehg, config):
        print("wind-restart: capability absent (Q6 wind schema)")
        return 1

    def add_checkpoint(doc):
        doc["output"]["checkpoint"] = {"interval": 26740}  # 7 x 3820, pre-cutoff
        doc["output"]["variables"]["surface"] = ["eta", "uu"]

    rewrite_config(config, add_checkpoint)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"
    restart_dir = stage_case(args.repo, "g9-wind", args.work / "restart")
    restart_config = restart_dir / "g10-seiche.yaml"

    def add_restart(doc):
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()),
                          "time": 26740}
        doc["output"]["filename"] = "out/restarted.h5"
        doc["output"]["variables"]["surface"] = ["eta", "uu"]

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)

    ok = True
    with h5py.File(full_output, "r") as full, \
            h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for var in ["eta", "uu"]:
            for t in ["30560", "34380", "38200"]:  # spanning the cutoff
                a = full[f"/surface/{var}/{t}"][:]
                b = part[f"/surface/{var}/{t}"][:]
                scale = max(float(np.abs(a).max()), 1.0e-12)
                diff = float(np.abs(a - b).max()) / scale
                print(f"  wind restart {var}@{t}: max rel diff = {diff:.3e}")
                if diff > 1.0e-12:
                    ok = False
    print("wind restart determinism (mid step-down):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_evap_restart(args: argparse.Namespace) -> int:
    """§8.4 restart-mid-transient row for evaporation: g4d (drain-to-dry)
    checkpointed mid-drying at t = 15000 s (dry-out at 20000 s), so the
    restarted leg carries the drying front across dry-out. Bitwise class."""
    case_dir = stage_case(args.repo, "g4-evap", args.work)
    config = case_dir / "g4d-drain.yaml"

    def add_checkpoint(doc):
        doc["output"]["checkpoint"] = {"interval": 15000}
        doc["output"]["filename"] = "out/straight.h5"

    rewrite_config(config, add_checkpoint)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "straight.h5"
    restart_dir = stage_case(args.repo, "g4-evap", args.work / "restart")
    restart_config = restart_dir / "g4d-drain.yaml"

    def add_restart(doc):
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()),
                          "time": 15000}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)

    ok = True
    with h5py.File(full_output, "r") as full, \
            h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for var in ["eta", "depth"]:
            for t in ["20000", "30000", "40000"]:  # spanning dry-out
                a = full[f"/surface/{var}/{t}"][:]
                b = part[f"/surface/{var}/{t}"][:]
                scale = max(float(np.abs(a).max()), 1.0e-12)
                diff = float(np.abs(a - b).max()) / scale
                print(f"  evap restart {var}@{t}: max rel diff = {diff:.3e}")
                if diff > 1.0e-12:
                    ok = False
    print("evaporation restart determinism (mid-drying):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_b6_subcycled(args: argparse.Namespace) -> int:
    """§8.3 gap closure for subcycling + transport (the b5/b6 gate split
    left the pair unexercised): b6-ss shortened to 1800 s under
    coupling.mode: subcycled, against its sync twin. Criteria: both
    transport_audit identities close (the pair's conservation content),
    salinity respects its bounds, and the surface volume-identity residual
    matches the sync twin's — the b6 regime (scaled-terrain, tidal eta BC,
    hard wetting) carries a pre-existing unclosed audit term (~30 % of the
    seepage flux on this 1800 s window, identical under sync, never gated
    by the b6 gates, recorded in dod-Q7), so the differential is the
    honest assertion: subcycling must not make the audit worse."""
    residuals = {}
    ledger_ok = True
    bounds_ok = True
    for mode in ("sync", "subcycled"):
        case_dir = stage_case(args.repo, "b6-kuan", args.work / mode)
        config = case_dir / "b6-kuan-ss.yaml"

        def shorten(doc, mode=mode):
            doc["coupling"]["mode"] = mode
            doc["time"]["t_end"] = 1800.0
            doc["time"]["output_interval"] = 900.0

        rewrite_config(config, shorten)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        with h5py.File(case_dir / "out-ss" / "output.h5", "r") as handle:
            audit = handle["/monitor/transport_audit"]
            cols = {n: k for k, n in
                    enumerate(audit.attrs["columns"].decode().split(","))}
            first, last = audit[0], audit[-1]

            def cum(name):
                return float(last[cols[name]] - first[cols[name]])

            d_surf = cum("surf_mass")
            d_subs = cum("subs_mass")
            surf_resid = d_surf - (cum("exchange") + cum("surf_source") +
                                   cum("surf_boundary") + cum("surf_adjust") +
                                   cum("surf_anchor"))
            subs_resid = d_subs - (-cum("exchange") + cum("subs_boundary") +
                                   cum("subs_adjust") + cum("subs_anchor"))
            scale = max(abs(d_surf) + abs(d_subs), 1.0)
            for name, resid in (("surface", surf_resid),
                                ("subsurface", subs_resid)):
                good = abs(resid) <= 1.0e-8 * scale
                print(f"  [{mode}] transport ledger [{name}]: residual "
                      f"{resid:.3e} (scale {scale:.3e})"
                      f"{'' if good else ' FAIL'}")
                ledger_ok &= good
            ma = handle["/monitor/mass_audit"]
            mcols = {n: k for k, n in
                     enumerate(ma.attrs["columns"].decode().split(","))}
            row0, rowN = ma[0], ma[-1]
            residuals[mode] = float(
                (rowN[mcols["volume"]] - row0[mcols["volume"]])
                + rowN[mcols["boundary_outflow"]] - rowN[mcols["rain"]]
                + rowN[mcols["evaporation"]] - rowN[mcols["bc_inflow"]]
                - rowN[mcols["seepage"]] - rowN[mcols["clamped"]])
            for group in ("/transport/concentration",
                          "/transport/concentration_surface"):
                field = handle[f"{group}/1800"][:]
                good = (float(field.min()) >= -1.0e-9
                        and float(field.max()) <= 35.0 + 1.0e-9)
                print(f"  [{mode}] bounds {group}: [{field.min():.3e}, "
                      f"{field.max():.4f}]{'' if good else ' FAIL'}")
                bounds_ok &= good
    drift = abs(residuals["subcycled"] - residuals["sync"])
    allowed = 0.1 * abs(residuals["sync"]) + 1.0e-6
    diff_ok = drift <= allowed
    print(f"  surface identity residual: sync {residuals['sync']:.4e}, "
          f"subcycled {residuals['subcycled']:.4e}, |drift| {drift:.3e} "
          f"(allowed {allowed:.3e}){'' if diff_ok else ' FAIL'}")
    ok = ledger_ok and bounds_ok and diff_ok
    print("b6 subcycled + transport smoke:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rank_invariance_tracer(args: argparse.Namespace) -> int:
    """The §7.1 decomposition lane for the (salinity) transport module —
    the twin of rank-invariance-heat, closing the 'salinity + MPI' §8.3
    cell (b6 is serial by design; g5 runs at 1 rank): the transport
    battery base at 1/2/4 ranks, strict and default modes, under BOTH
    advection schemes (its first execution caught the missing
    unconditional Fu/Fv halo exchange at 1.6e-2 under upwind / 1.6e-3
    under superbee — the flow itself never reads those halos, so no
    other lane could see it). Strict bound 1e-10: measured 9.0e-12 after
    the fix — the TVD ratio amplifies rounding near flat gradients, so
    the subsurface lanes' 1e-12 is not attainable on this path. A third
    lane runs superbee under transport.surface_update: consistent, whose
    flux volume reads the current Fu/Fv halos and whose limiter covers
    still cells (the 'surface_update + MPI' §8.3 cell)."""
    strict = args.mode == "strict"
    extra = STRICT_PETSC_OPTIONS if strict else []
    limit = 1.0e-10 if strict else 1.0e-4
    ok = True
    for scheme, update in (("superbee", "legacy"), ("upwind", "legacy"),
                           ("superbee", "consistent")):
        lane = scheme if update == "legacy" else f"{scheme}-{update}"
        outputs = {}
        for ranks in (1, 2, 4):
            case_dir = stage_case(args.repo, "x-orient",
                                  args.work / lane / f"n{ranks}")
            config = case_dir / "transport-orient-base.yaml"

            def set_scheme(doc, scheme=scheme, update=update):
                doc["transport"]["scheme"]["advection"] = scheme
                doc["transport"]["surface_update"] = update

            rewrite_config(config, set_scheme)
            run_case(args.frehg, config, case_dir, args.mpiexec, ranks, extra)
            outputs[ranks] = case_dir / "out" / "output.h5"
        with h5py.File(outputs[1], "r") as base:
            for ranks in (2, 4):
                with h5py.File(outputs[ranks], "r") as other:
                    worst = 0.0
                    for group in ("/surface/eta",
                                  "/transport/concentration_surface"):
                        for t in base[group]:
                            a = base[f"{group}/{t}"][:]
                            b = other[f"{group}/{t}"][:]
                            scale = max(float(np.abs(a).max()), 1.0e-12)
                            worst = max(worst,
                                        float(np.abs(a - b).max()) / scale)
                    print(f"  [{lane}] n=1 vs n={ranks} [{args.mode}]: "
                          f"max rel diff = {worst:.3e} (allowed {limit:.0e})")
                    ok &= worst <= limit
    print(f"tracer rank invariance ({args.mode}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rank_invariance_coupled_salt(args: argparse.Namespace) -> int:
    """The §7.1 decomposition lane for the consistent subsurface salinity
    update (transport.subsurface_update: consistent with
    groundwater.density_coupling.reallocation_gradient: density) — the
    'subsurface_update + MPI' §8.3 cell: x-orient/coupled-salt-base.yaml at
    1/2/4 ranks. The window volumes, the per-face post-allocation path, the
    vent volume and the transport-step mass basis all live in the halo
    layout, so a staging error shows at the interfaces immediately. The
    120 s horizon stays ahead of the saturation-front cell crossings that
    amplify reduction-order rounding in every coupled lane (P2 finding;
    measured here ~1e-3 by 3600 s for the legacy update too). Strict bound
    1e-7: measured 1.1e-8. Default bound 2e-2: measured 7.3e-3, set by the
    flow itself — the per-rank block-Jacobi/ICC preconditioner at rtol 1e-8
    moves the unsaturated water content by 4.6e-3..5.7e-3 between rank
    counts in the legacy update as well."""
    strict = args.mode == "strict"
    extra = (STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW) if strict else []
    limit = 1.0e-7 if strict else 2.0e-2
    outputs = {}
    for ranks in (1, 2, 4):
        case_dir = stage_case(args.repo, "x-orient", args.work / f"n{ranks}")
        config = case_dir / "coupled-salt-base.yaml"
        decomposition = {1: (1, 1), 2: (2, 1), 4: (2, 2)}[ranks]

        def set_ranks(doc, decomposition=decomposition):
            doc["domain"]["decomposition"] = {"mpi_nx": decomposition[0],
                                              "mpi_ny": decomposition[1]}

        rewrite_config(config, set_ranks)
        run_case(args.frehg, config, case_dir, args.mpiexec, ranks, extra)
        outputs[ranks] = case_dir / "out" / "output.h5"
    ok = True
    with h5py.File(outputs[1], "r") as base:
        salted = float(np.nanmax(base["/transport/concentration/120"][:]))
        moved = salted > 10.5
        print(f"  n=1: subsurface salinity reaches {salted:.3f} psu (needs > 10.5)"
              f"{'' if moved else ' FAIL'}")
        ok &= moved
        for ranks in (2, 4):
            with h5py.File(outputs[ranks], "r") as other:
                worst, where = 0.0, ""
                for group in ("/surface/eta", "/groundwater/water_content",
                              "/transport/concentration",
                              "/transport/concentration_surface"):
                    for t in base[group]:
                        a = base[f"{group}/{t}"][:]
                        b = other[f"{group}/{t}"][:]
                        scale = max(float(np.nanmax(np.abs(a))), 1.0e-12)
                        diff = float(np.nanmax(np.abs(a - b))) / scale
                        if diff > worst:
                            worst, where = diff, f"{group}/{t}"
                print(f"  n=1 vs n={ranks} [{args.mode}]: max rel diff = "
                      f"{worst:.3e} at {where} (allowed {limit:.0e})")
                ok &= worst <= limit
    print(f"coupled salinity rank invariance ({args.mode}):", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_wind_drydown(args: argparse.Namespace) -> int:
    """§8.3 pair closure wind + wetting/drying (the v1 lesson's named
    pair; g9c deliberately grazes without crossing): the g9c sloping
    basin under a 26 m/s wind, strong enough to dry the shallow tip.
    Criteria: cells actually dry, depths stay non-negative, the clamp
    mints nothing, and the volume budget closes."""
    case_dir = stage_case(args.repo, "g9-wind", args.work)
    config = case_dir / "g9c-slope.yaml"
    if not validate_config(args.frehg, config):
        print("wind-drydown: capability absent (Q6 wind schema)")
        return 1

    def strengthen(doc):
        doc["surface_water"]["wind"]["u10"] = {"constant": 26.0}
        doc["time"]["t_end"] = 20000.0
        doc["time"]["output_interval"] = 2000.0

    rewrite_config(config, strengthen)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        times = output_times(handle, "/surface/depth")
        dry_cells = 0
        min_depth = np.inf
        for t in times:
            depth = handle[f"/surface/depth/{t}"][:]
            dry_cells = max(dry_cells, int(np.sum(depth <= 1.0e-6)))
            min_depth = min(min_depth, float(depth.min()))
        ma = handle["/monitor/mass_audit"]
        cols = {n: k for k, n in
                enumerate(ma.attrs["columns"].decode().split(","))}
        row0, rowN = ma[0], ma[-1]
        resid = float((rowN[cols["volume"]] - row0[cols["volume"]])
                      + rowN[cols["boundary_outflow"]] - rowN[cols["rain"]]
                      + rowN[cols["evaporation"]] - rowN[cols["bc_inflow"]]
                      - rowN[cols["clamped"]])
        clamped = float(rowN[cols["clamped"]])
        vol = float(rowN[cols["volume"]])
    checks = [
        ("tip dries under wind", dry_cells >= 1, f"{dry_cells} dry cell(s)"),
        ("no negative depth", min_depth >= 0.0, f"min {min_depth:.3e} m"),
        # The below-bed clamp is the documented legacy-scheme mechanism at
        # hard drying fronts (output.md: "report it, don't hide it");
        # driving the tip dry with a 26 m/s wind engages it by design, so
        # the bound is relative (measured 5.3 m^3 of 10605 = 5e-4).
        ("clamp bounded", clamped <= 1.0e-3 * max(vol, 1.0),
         f"{clamped:.3e} m^3 minted ({clamped / max(vol, 1.0):.1e} of volume)"),
        ("closure", abs(resid) <= 1.0e-6 * max(vol, 1.0),
         f"|{resid:.3e}| m^3 of {vol:.1f}"),
    ]
    for name, good, detail in checks:
        print(f"  wind-drydown: {name}: {detail}{'' if good else ' FAIL'}")
        ok &= good
    print("wind + wetting/drying:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def temperature_ledger_residuals(handle) -> tuple[float, float, float]:
    """(surface residual, subsurface residual, scale) of the documented
    temperature_audit identities (output.md: the transport identities on
    surf_heat/subs_heat, with the atmospheric exchange in surf_atmos)."""
    audit = handle["/monitor/temperature_audit"]
    cols = {n: k for k, n in
            enumerate(audit.attrs["columns"].decode().split(","))}
    first, last = audit[0], audit[-1]

    def cum(name):
        return float(last[cols[name]] - first[cols[name]])

    surf = cum("surf_heat") - (cum("exchange") + cum("surf_source") +
                               cum("surf_boundary") + cum("surf_adjust") +
                               cum("surf_anchor") + cum("surf_atmos"))
    subs = cum("subs_heat") - (-cum("exchange") + cum("subs_boundary") +
                               cum("subs_adjust") + cum("subs_anchor"))
    scale = max(abs(cum("surf_heat")) + abs(cum("subs_heat")), 1.0)
    return surf, subs, scale


def gate_terrain_heat(args: argparse.Namespace) -> int:
    """§8.3 pair closure terrain-following + heat (every g6/g7/g8 case is
    a regular flat mesh): a sloped-bed saturated block on the scaled
    terrain mesh with lateral flow and a side thermal Dirichlet. The
    conservation criterion is the temperature ledger on the tilted
    metrics; bounds guard the tilted-face flux arithmetic."""
    work = args.work / "terrain-heat"
    (work / "input").mkdir(parents=True, exist_ok=True)
    nx, ny, nz = 8, 1, 6
    with open(work / "input" / "bed.dat", "w", encoding="utf-8") as handle:
        handle.write("# sloped bed for the terrain-heat pair gate\n")
        for i in range(nx):
            handle.write(f"{2.0 + 0.1 * i:.6f}\n")
    doc = {
        "simulation": {"id": "terrain-heat",
                       "title": "terrain-following thermal conservation"},
        "domain": {"nx": nx, "ny": ny, "nz": nz, "dx": 0.5, "dy": 0.5,
                   "dz": 0.25, "bottom_elevation": {"file": "input/bed.dat"},
                   "follow_terrain": True, "terrain_layers": "scaled"},
        "time": {"dt": 200.0, "t_end": 8000.0, "output_interval": 8000.0},
        "modules": {"groundwater": True, "temperature": True},
        "groundwater": {"specific_storage": 1.0e-6,
                        "timestep": {"dt_init": 200.0, "dt_min": 200.0,
                                     "dt_max": 200.0}},
        "soil": {"types": [{"name": "gravel", "ksx": 1.0e-3, "ksy": 1.0e-3,
                            "ksz": 1.0e-3, "theta_s": 0.35, "theta_r": 0.02,
                            "vg_alpha": 3.0, "vg_n": 2.0}],
                 "map": {"constant": "gravel"}},
        "temperature": {"scheme": {"advection": "superbee"},
                        "thermal_conductivity": 2.0,
                        "heat_capacity_water": 4.184e6,
                        "heat_capacity_solid": 2.0548e6,
                        "dispersivity": {"longitudinal": 0.0,
                                         "transverse": 0.0}},
        "initial_conditions": {
            "groundwater": {"water_table": {"constant": 3.0}},
            "temperature": {"groundwater": {"constant": 20.0}}},
        "boundary_conditions": [
            {"name": "head-in",
             "region": {"polygon": [[-0.1, -0.1], [0.3, -0.1], [0.3, 0.6],
                                    [-0.1, 0.6]]},
             "target": "groundwater_side", "kind": "head",
             "value": {"hydrostatic": {"eta": 3.05}}},
            {"name": "head-out",
             "region": {"polygon": [[3.7, -0.1], [4.1, -0.1], [4.1, 0.6],
                                    [3.7, 0.6]]},
             "target": "groundwater_side", "kind": "head",
             "value": {"hydrostatic": {"eta": 2.95}}},
            {"name": "warm-in",
             "region": {"polygon": [[-0.1, -0.1], [0.3, -0.1], [0.3, 0.6],
                                    [-0.1, 0.6]]},
             "target": "groundwater_side", "kind": "scalar_value",
             "scalar": "temperature", "value": {"constant": 25.0}}],
        "output": {"filename": "out/output.h5",
                   "variables": {"temperature": ["temperature"]}},
    }
    config = work / "terrain-heat.yaml"
    with open(config, "w", encoding="utf-8") as handle:
        yaml.safe_dump(doc, handle, sort_keys=False)
    run_case(args.frehg, config, work, args.mpiexec, args.ranks,
             STRICT_PETSC_OPTIONS_GW)
    ok = True
    with h5py.File(work / "out" / "output.h5", "r") as handle:
        surf, subs, scale = temperature_ledger_residuals(handle)
        field = handle["/temperature/temperature/8000"][:]
        t_lo, t_hi = float(field.min()), float(field.max())
    checks = [
        ("subsurface heat ledger", abs(subs) <= 1.0e-8 * scale,
         f"residual {subs:.3e} (scale {scale:.3e})"),
        ("bounds", 20.0 - 1.0e-9 <= t_lo and t_hi <= 25.0 + 1.0e-9,
         f"[{t_lo:.4f}, {t_hi:.4f}] C"),
        ("side heat arrived", t_hi >= 20.5,
         f"max {t_hi:.4f} C (Dirichlet 25)"),
    ]
    for name, good, detail in checks:
        print(f"  terrain-heat: {name}: {detail}{'' if good else ' FAIL'}")
        ok &= good
    print("terrain-following + heat:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_two_scalar_restart(args: argparse.Namespace) -> int:
    """§8.3 pair closures for two live scalars (transport + temperature
    co-resident — the V2-A17 two-scalar design's first co-execution) and
    for the density coupling with BOTH beta_s and beta_T active: the heat
    battery base with salinity on, checkpointed and restarted; straight vs
    restarted must agree bitwise-class on temperature, salinity, and head,
    and both scalar ledgers must close."""
    case_dir = stage_case(args.repo, "g6-heat", args.work)
    base = case_dir / "orient-base.yaml"
    if not validate_config(args.frehg, base):
        print("two-scalar-restart: capability absent (Q5 schema)")
        return 1
    nx = ny = 8
    nz = 6
    ic = heat_orient_ic(nx, ny, nz, 0.5, 0.5, 0.25)
    (case_dir / "input").mkdir(exist_ok=True)
    write_ic_raster(case_dir / "input" / "orient_ic_id.dat", ic)

    def two_scalar(doc):
        doc["initial_conditions"]["temperature"]["groundwater"] = {
            "file": "input/orient_ic_id.dat"}
        doc["modules"]["transport"] = True
        doc["initial_conditions"]["transport"] = {
            "groundwater": {"constant": 5.0}}
        doc["transport"] = {"scheme": {"advection": "superbee"},
                            "bounds": {"min": 0.0}}
        doc["output"]["variables"]["transport"] = ["concentration"]
        # both densities live: the committed base already carries beta_T
        doc["groundwater"]["density_coupling"]["beta_saline"] = 7.44e-4
        doc["boundary_conditions"].append(
            {"name": "salt-in",
             "region": {"polygon": [[-0.1, -0.1], [0.3, -0.1], [0.3, 4.1],
                                    [-0.1, 4.1]]},
             "target": "groundwater_side", "kind": "scalar_value",
             "scalar": "salinity", "value": {"constant": 10.0}})
        return doc

    def add_checkpoint(doc):
        two_scalar(doc)
        doc["output"]["checkpoint"] = {"interval": 10000}
        doc["output"]["filename"] = "out/straight.h5"

    def add_restart(doc):
        two_scalar(doc)
        doc["restart"] = {"enabled": True, "file": "out/straight.h5",
                          "time": 10000}
        doc["output"]["filename"] = "out/restarted.h5"

    straight = case_dir / "straight.yaml"
    shutil.copy(base, straight)
    rewrite_config(straight, add_checkpoint)
    run_case(args.frehg, straight, case_dir, args.mpiexec, args.ranks)
    restarted = case_dir / "restarted.yaml"
    shutil.copy(base, restarted)
    rewrite_config(restarted, add_restart)
    run_case(args.frehg, restarted, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "straight.h5", "r") as a, \
            h5py.File(case_dir / "out" / "restarted.h5", "r") as b:
        for group in ("/temperature/temperature", "/transport/concentration",
                      "/groundwater/hydraulic_head"):
            t = 20000
            va = a[f"{group}/{t}"][:]
            vb = b[f"{group}/{t}"][:]
            scale = max(float(np.abs(va).max()), 1.0e-12)
            diff = float(np.abs(va - vb).max()) / scale
            print(f"  two-scalar restart {group}@{t}: max rel diff = {diff:.3e}")
            ok &= diff <= 1.0e-12
        for name, audit, mass_col in (
                ("transport", "/monitor/transport_audit", "subs_mass"),
                ("temperature", "/monitor/temperature_audit", "subs_heat")):
            table = a[audit]
            cols = {n: k for k, n in
                    enumerate(table.attrs["columns"].decode().split(","))}
            first, last = table[0], table[-1]

            def cum(col):
                return float(last[cols[col]] - first[cols[col]])

            resid = cum(mass_col) - (-cum("exchange") + cum("subs_boundary")
                                     + cum("subs_adjust") + cum("subs_anchor"))
            scale = max(abs(cum(mass_col)), 1.0)
            good = abs(resid) <= 1.0e-8 * scale
            print(f"  {name} ledger: residual {resid:.3e} (scale {scale:.3e})"
                  f"{'' if good else ' FAIL'}")
            ok &= good
    print("two-scalar restart (transport + temperature, both betas):",
          "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_evap_heat(args: argparse.Namespace) -> int:
    """§8.3 pair closure for the two atmosphere consumers co-active (bulk
    evaporation + bulk surface heat exchange share the Q4/Q5 atm chain but
    no case runs both): g7a2-bulk with `evaporation.mode: bulk` added.
    Criteria: evaporation actually removes volume, the volume budget
    closes on the evaporating basin, the heat ledger closes, and the
    temperature stays finite and bounded by the forcing class."""
    case_dir = stage_case(args.repo, "g7-heat", args.work)
    config = case_dir / "g7a2-bulk.yaml"
    if not validate_config(args.frehg, config):
        print("evap-heat: capability absent (Q4/Q5 schema)")
        return 1

    def add_bulk_evap(doc):
        doc["surface_water"]["evaporation"] = {"mode": "bulk"}
        # the Q4 convention: a bulk-evaporation consumer requires the
        # prescribed skin temperature (heat exchange keeps using local T)
        doc["atmosphere"]["surface_temperature"] = {"constant": 22.0}
        doc["output"]["variables"]["surface"] = ["eta", "depth"]

    rewrite_config(config, add_bulk_evap)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        ma = handle["/monitor/mass_audit"]
        cols = {n: k for k, n in
                enumerate(ma.attrs["columns"].decode().split(","))}
        row0, rowN = ma[0], ma[-1]
        evap = float(rowN[cols["evaporation"]])
        resid = float((rowN[cols["volume"]] - row0[cols["volume"]])
                      + rowN[cols["boundary_outflow"]] - rowN[cols["rain"]]
                      + evap - rowN[cols["bc_inflow"]]
                      - rowN[cols["clamped"]])
        surf, _, scale = temperature_ledger_residuals(handle)
        times = output_times(handle, "/temperature/temperature_surface")
        temp = handle[f"/temperature/temperature_surface/{times[-1]}"][:]
    checks = [
        ("bulk evaporation active", evap > 0.0, f"{evap:.4e} m^3 evaporated"),
        ("closure", abs(resid) <= 1.0e-6 * max(evap, 1.0e-6),
         f"|{resid:.3e}| m^3 vs {evap:.3f} evaporated"),
        ("surface heat ledger", abs(surf) <= 1.0e-8 * scale,
         f"residual {surf:.3e} (scale {scale:.3e})"),
        ("temperature finite/bounded", np.isfinite(temp).all()
         and 0.0 <= float(temp.min()) and float(temp.max()) <= 60.0,
         f"[{temp.min():.3f}, {temp.max():.3f}] C"),
    ]
    for name, good, detail in checks:
        print(f"  evap-heat: {name}: {detail}{'' if good else ' FAIL'}")
        ok &= good
    print("bulk evaporation + bulk heat exchange (co-active atm):",
          "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_heat_subcycled(args: argparse.Namespace) -> int:
    """§8.3 pair closure subcycling + temperature (g6c-coupled is sync):
    the g6c coupled pond column shortened, sync vs subcycled twins; the
    temperature ledgers must close in both and the two modes' final
    temperature fields must agree to the coupling-window class."""
    fields = {}
    ok = True
    for mode in ("sync", "subcycled"):
        case_dir = stage_case(args.repo, "g6-heat", args.work / mode)
        config = case_dir / "g6c-coupled.yaml"
        if not validate_config(args.frehg, config):
            print("heat-subcycled: capability absent (Q5 schema)")
            return 1

        def shorten(doc, mode=mode):
            doc.setdefault("coupling", {})["mode"] = mode
            doc["time"]["t_end"] = 21600.0
            doc["time"]["output_interval"] = 21600.0

        rewrite_config(config, shorten)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            _, subs, scale = temperature_ledger_residuals(handle)
            good = abs(subs) <= 1.0e-8 * scale
            print(f"  [{mode}] temperature ledger: residual {subs:.3e} "
                  f"(scale {scale:.3e}){'' if good else ' FAIL'}")
            ok &= good
            fields[mode] = handle["/temperature/temperature/21600"][:]
    diff = float(np.abs(fields["sync"] - fields["subcycled"]).max())
    # The two coupling windows are different discretizations; the fields
    # must agree at the coupling-window truncation class, not bitwise.
    good = diff <= 0.5
    print(f"  sync vs subcycled max |dT| = {diff:.4f} K (allowed 0.5)"
          f"{'' if good else ' FAIL'}")
    ok &= good
    print("subcycling + temperature smoke:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Stage-3 default-bitwise gate (Stage-3 plan §3.3, label stage3_local): with
# no Stage-3 key, the listed per-PR cases must reproduce the T0 reference
# hashes dataset for dataset. BITWISE_PINS holds the v2.0-equivalent values
# of every flippable default (§2.6), injected so the gate survives a T8
# default flip; it is extended by the phase that adds each key.
# ---------------------------------------------------------------------------

BITWISE_CASES = [
    # (label, benchmark dir, config file)
    ("b1", "b1-sw", "b1-sw.yaml"),
    ("b2", "b2-gw", "b2-gw.yaml"),
    ("g4a", "g4-evap", "g4a-drawdown.yaml"),
    ("g4b", "g4-evap", "g4b-concentration.yaml"),
    ("g4d", "g4-evap", "g4d-drain.yaml"),
    ("g6c-coupled", "g6-heat", "g6c-coupled.yaml"),
    ("g7b", "g7-heat", "g7b-channel.yaml"),
    ("g10", "g9-wind", "g10-seiche.yaml"),
    ("swe-orient", "x-orient", "swe-orient-base.yaml"),
    ("gw-orient", "x-orient", "gw-orient-base.yaml"),
    ("transport-orient", "x-orient", "transport-orient-base.yaml"),
    ("conveyor", "x-orient", "conveyor-base.yaml"),
    ("coupled-salt", "x-orient", "coupled-salt-base.yaml"),
]
# dotted key -> v2.0-equivalent value (module-gated like the §3.2 overrides).
# Stage-3 T8 (W8.1, S3-A11): the table §2.6 promised was never populated by
# T1/T5 — the W8.1 evaluation found the gate riding the shipped default
# (b2/g6c-coupled went non-bitwise under the flipped binary). These pins
# make the gate test "the v2.0 configuration", not "the current default".
BITWISE_PINS: dict[str, object] = {
    "groundwater.reallocation": "legacy",
    "transport.surface_dryout": "discard",
}


def gate_default_bitwise(args: argparse.Namespace) -> int:
    sys.path.insert(0, str(HERE.parent.parent / "tools"))
    import h5_bitwise  # noqa: E402

    reference = HERE / "reference-hashes" / "stage3-T0.json"
    if not args.record and not reference.exists():
        print(f"default-bitwise: reference {reference} missing; run with --record")
        return 1
    if not args.record:
        doc = json.loads(reference.read_text())
        if doc["toolchain"] != h5_bitwise.toolchain_fingerprint():
            print("default-bitwise: toolchain differs from the recording "
                  "toolchain; SKIPPED (exit 77)")
            return 77
    allow = list(args.allow or [])
    ok = True
    prefix_doc = load_mass_audit_prefix_doc(args)
    for label, case_dir_name, config_name in BITWISE_CASES:
        case_dir = stage_case(args.repo, case_dir_name, args.work / label)
        config = case_dir / config_name
        with open(config, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        pinned = []
        for dotted, value in BITWISE_PINS.items():
            if not override_applies(doc, dotted):
                continue
            node = doc
            parts = dotted.split(".")
            for part in parts[:-1]:
                node = node.setdefault(part, {})
            node[parts[-1]] = value
            pinned.append(dotted)
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        if pinned:
            print(f"  pins [{label}]: {', '.join(pinned)}")
        run_case(args.frehg, config, case_dir, args.mpiexec, 1)
        output = case_dir / doc["output"]["filename"]
        # Stage-3 T9 (ws5): the mass_audit table is a declared exception
        # (it gains appended columns), but its PREFIX — the columns the T8
        # binary wrote — must stay bitwise.
        if not mass_audit_prefix_check(output, label, args, prefix_doc):
            ok = False
        if args.record:
            argv = ["hash", str(output), "--case", label, "--out", str(reference)]
            if h5_bitwise.main_from(argv) != 0:
                return 1
        else:
            argv = ["compare", str(output), str(reference), "--case", label]
            for glob in allow:
                argv += ["--allow", glob]
            code = h5_bitwise.main_from(argv)
            if code == 77:
                return 77
            ok = ok and (code == 0)
    save_mass_audit_prefix_doc(args, prefix_doc)
    if args.record:
        print(f"default-bitwise: reference recorded -> {reference}")
        return 0
    if args.record_massaudit:
        return 0
    print("default-bitwise:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# Stage-3 T9 (plan §12A ws5): the mass_audit PREFIX pin. The T8 binary's
# table hashes (and widths) per case live in stage3-T8-massaudit.json,
# recorded once with --record-massaudit on the frozen T8 binary; a T9 run
# must reproduce those columns bitwise (only appended columns may differ).
MASS_AUDIT_PREFIX = HERE / "reference-hashes" / "stage3-T8-massaudit.json"


def load_mass_audit_prefix_doc(args) -> dict:
    if getattr(args, "record_massaudit", False):
        return {"cases": {}}
    if MASS_AUDIT_PREFIX.exists():
        return json.loads(MASS_AUDIT_PREFIX.read_text())
    return {"cases": {}}


def save_mass_audit_prefix_doc(args, doc: dict) -> None:
    if getattr(args, "record_massaudit", False):
        sys.path.insert(0, str(HERE.parent.parent / "tools"))
        import h5_bitwise  # noqa: E402
        doc["toolchain"] = h5_bitwise.toolchain_fingerprint()
        MASS_AUDIT_PREFIX.write_text(json.dumps(doc, indent=1, sort_keys=True))
        print(f"mass_audit prefix reference recorded -> {MASS_AUDIT_PREFIX}")


def mass_audit_prefix_check(output: Path, label: str, args, doc: dict) -> bool:
    import hashlib
    with h5py.File(output, "r") as handle:
        if "/monitor/mass_audit" not in handle:
            return True
        table = handle["/monitor/mass_audit"][()]
    if getattr(args, "record_massaudit", False):
        doc["cases"][label] = {
            "ncols": int(table.shape[1]),
            "sha256": hashlib.sha256(np.ascontiguousarray(table).tobytes()).hexdigest()}
        return True
    ref = doc.get("cases", {}).get(label)
    if ref is None:
        print(f"  mass_audit prefix [{label}]: no T8 reference — FAIL")
        return False
    n = ref["ncols"]
    digest = hashlib.sha256(np.ascontiguousarray(table[:, :n]).tobytes()).hexdigest()
    same = digest == ref["sha256"]
    print(f"  mass_audit prefix [{label}]: first {n} of {table.shape[1]} columns "
          f"{'bitwise' if same else 'DIFFER'} vs the T8 binary")
    return same



# ---------------------------------------------------------------------------
# Stage-3 T1: the reallocation-mode gates rl2/rl3/rl4 (Stage-3 plan §5).
# Pure checks live in gate_realloc.py; tolerances in
# tolerances/stage3-rl.yaml; §6.3 negative battery in
# scripts/test_rl_gates.py. Each gate runs a `legacy` and an `off` lane by
# rewriting groundwater.reallocation in the staged configuration.
# ---------------------------------------------------------------------------

import gate_realloc  # noqa: E402

CAPABILITY_ABSENT_S3RL = ("capability absent: the schema rejects the T1 "
                          "target keys (groundwater.reallocation / the "
                          "realloc_* output variables)")


def rl_tolerances(gate: str) -> dict:
    with open(HERE / "tolerances" / "stage3-rl.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)[gate]


def rl_stage_lane(args: argparse.Namespace, case: str, config_name: str,
                  lane: str, extra_transform=None) -> Path | None:
    """Stage `case` into work/<lane>, set groundwater.reallocation = lane,
    and validate; None means the capability is absent (the red state)."""
    case_dir = stage_case(args.repo, case, args.work / lane)
    config = case_dir / config_name
    def set_lane(doc):
        doc.setdefault("groundwater", {})["reallocation"] = lane
        if extra_transform is not None:
            extra_transform(doc)
    rewrite_config(config, set_lane)
    if not validate_config(args.frehg, config):
        return None
    return config


def rl_field_sum(output: Path, var: str, t: int) -> float:
    """NaN-aware sum of a 3D groundwater diagnostic at output time t."""
    with h5py.File(output, "r") as handle:
        data = handle[f"/groundwater/{var}/{t}"][:]
    return float(np.nansum(data))


def rl_substeps(config: Path) -> int:
    """Groundwater substep count from the run record's timer tree."""
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    record_path = (config.parent / Path(doc["output"]["filename"])).parent / "run-record.yaml"
    with open(record_path, encoding="utf-8") as handle:
        record = yaml.safe_load(handle)
    counts = [entry["count"] for path, entry in record["timers"].items()
              if path.split("/")[-1] == "groundwater"]
    return int(max(counts))


def rl_soil_bounds(config: Path) -> tuple[float, float]:
    """(min theta_r, max theta_s) over the staged soil table."""
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    types = doc["soil"]["types"]
    return (min(float(t["theta_r"]) for t in types),
            max(float(t["theta_s"]) for t in types))


def gate_rl5_orient(args: argparse.Namespace) -> int:
    # S3-A4: the disabled-mode orientation battery runs on a VADOSE base
    # (benchmarks/rl-realloc/rl5-orient-base.yaml, no cell near saturation).
    # The saturated x-orient gw block is outside the mode's domain — there
    # the missing consistency restore amplifies rounding chaotically and NO
    # strict battery can pass (measured O(1); Reallocate.cpp header). The
    # harness sets reallocation: disabled through FREHG_STAGE3_OVERRIDES
    # (asserted into each run record by the r1 --expect hook).
    return run_dihedral_battery(
        args, "rl5-orient", "rl5-orient-base.yaml",
        {"/groundwater/hydraulic_head": 8, "/groundwater/water_content": 8},
        STRICT_PETSC_OPTIONS_GW, case="rl-realloc", tol_key="gw")


def gate_rl_restart(args: argparse.Namespace) -> int:
    """Stage-3 restart gate, rl lane (plan §3.4; W1.3): the per-cell
    cumulative diagnostics (realloc_created / realloc_dropped / vloss_cum)
    checkpoint and restart BITWISE — every dataset of the restarted rl2
    disabled run equals the uninterrupted run's. The diagnostics are
    requested through the case's own output list, so both sides carry the
    same checkpoint layout."""
    case_dir = stage_case(args.repo, "rl-realloc", args.work)
    config = case_dir / "rl2-tidal-box.yaml"

    def full_setup(doc):
        doc["groundwater"]["reallocation"] = "disabled"
        doc["output"]["checkpoint"] = {"interval": 7200}

    rewrite_config(config, full_setup)
    if not validate_config(args.frehg, config):
        print(f"rl restart gate: {CAPABILITY_ABSENT_S3RL}")
        print("rl restart determinism: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"

    restart_dir = stage_case(args.repo, "rl-realloc", args.work / "restart")
    restart_config = restart_dir / "rl2-tidal-box.yaml"

    def add_restart(doc):
        full_setup(doc)
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()),
                          "time": 7200}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)

    ok = True
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for var in ["hydraulic_head", "water_content", "realloc_created",
                    "realloc_dropped", "vloss"]:
            for t in output_times(part, f"/groundwater/{var}"):
                if t <= 7200:
                    continue
                a = full[f"/groundwater/{var}/{t}"][:]
                b = part[f"/groundwater/{var}/{t}"][:]
                same = np.array_equal(np.nan_to_num(a, nan=-9.0e99),
                                      np.nan_to_num(b, nan=-9.0e99))
                if not same or t == output_times(part, f"/groundwater/{var}")[-1]:
                    print(f"  restart {var}@{t}: {'bitwise' if same else 'DIFFERS'}")
                ok = ok and same
    print("rl restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rl2(args: argparse.Namespace) -> int:
    tol = rl_tolerances("rl2")
    lanes: dict[str, dict] = {}
    for lane in ("legacy", "disabled"):
        config = rl_stage_lane(args, "rl-realloc", "rl2-tidal-box.yaml", lane)
        if config is None:
            print(f"rl2 gate: {CAPABILITY_ABSENT_S3RL}")
            print("rl2 gate: FAIL")
            return 1
        run_case(args.frehg, config, config.parent, args.mpiexec, args.ranks)
        output = config.parent / "out" / "output.h5"
        with h5py.File(output, "r") as handle:
            audit = handle["/monitor/gw_mass_audit"][:]
            t_last = output_times(handle, "/groundwater/water_content")[-1]
            wc = handle[f"/groundwater/water_content/{t_last}"][:]
        lanes[lane] = {
            "audit": audit,
            "wc": wc,
            "created": rl_field_sum(output, "realloc_created", t_last),
            "dropped": rl_field_sum(output, "realloc_dropped", t_last),
            "vloss": rl_field_sum(output, "vloss", t_last),
            "substeps": rl_substeps(config),
            "config": config,
        }
    ok = True

    def report(label: str, result) -> None:
        nonlocal ok
        good, msg = result
        print(f"  {label}: {msg}{'' if good else ' FAIL'}")
        ok = ok and good

    theta_lo, theta_hi = rl_soil_bounds(lanes["disabled"]["config"])
    off = lanes["disabled"]
    legacy = lanes["legacy"]
    report("disabled ledger", gate_realloc.check_ledger(off["audit"], tol["ledger_rel"]))
    report("disabled realloc exactly zero",
           gate_realloc.check_realloc_exactly_zero(off["audit"]))
    report("disabled theta bounds",
           gate_realloc.check_theta_bounds(off["wc"], theta_lo, theta_hi,
                                           tol["theta_slack"]))
    # S3-A4 domain-of-applicability record: this box saturates (that is what
    # makes its legacy lane exercise reallocation), which is OUTSIDE the
    # disabled mode's vadose domain — its clamp cost and substep count are
    # RECORDED here as the owner's C-RL data, and BOUNDED on the vadose rl3
    # case where the mode is meant to run.
    d_off = gate_realloc.deltas(off["audit"])
    print(f"  disabled (recorded, S3-A4): vloss {d_off[gate_realloc.VLOSS]:.3e} m^3 "
          f"over T {gate_realloc.throughput(off['audit']):.3e} m^3; substeps "
          f"{off['substeps']} vs legacy {legacy['substeps']}")
    for lane in ("legacy", "disabled"):
        report(f"{lane} per-cell sums",
               gate_realloc.check_percell_sums(lanes[lane]["created"],
                                               lanes[lane]["dropped"],
                                               lanes[lane]["vloss"],
                                               lanes[lane]["audit"],
                                               tol["percell_rel"]))
        report(f"{lane} ledger", gate_realloc.check_ledger(lanes[lane]["audit"],
                                                           tol["ledger_rel"]))
    report("legacy churn (anti-vacuity)",
           gate_realloc.check_legacy_churn(legacy["audit"],
                                           tol["churn_floor_rel"]))
    print("rl2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rl3(args: argparse.Namespace) -> int:
    tol = rl_tolerances("rl3")
    results: dict[str, dict] = {}
    for lane in ("legacy", "disabled"):
        config = rl_stage_lane(args, "rl-realloc", "rl3-steady-column.yaml", lane)
        if config is None:
            print(f"rl3 gate: {CAPABILITY_ABSENT_S3RL}")
            print("rl3 gate: FAIL")
            return 1
        run_case(args.frehg, config, config.parent, args.mpiexec, args.ranks)
        output = config.parent / "out" / "output.h5"
        with h5py.File(output, "r") as handle:
            audit = handle["/monitor/gw_mass_audit"][:]
            times = output_times(handle, "/groundwater/water_content")
            window = times[int(0.9 * (len(times) - 1)):]
            first = handle[f"/groundwater/water_content/{window[0]}"][:]
            last = handle[f"/groundwater/water_content/{window[-1]}"][:]
            qz = handle[f"/groundwater/qz/{times[-1]}"][:]
        results[lane] = {"audit": audit, "first": first, "last": last, "qz": qz,
                         "substeps": rl_substeps(config)}
    ok = True

    def report(label: str, result) -> None:
        nonlocal ok
        good, msg = result
        print(f"  {label}: {msg}{'' if good else ' FAIL'}")
        ok = ok and good

    off = results["disabled"]
    report("disabled steadiness", gate_realloc.check_steady(off["first"], off["last"],
                                                       tol["steady_abs"]))
    report("disabled flux uniformity",
           gate_realloc.check_uniform_flux(off["qz"], -2.0e-6, tol["uniform_rel"]))
    report("disabled realloc exactly zero",
           gate_realloc.check_realloc_exactly_zero(off["audit"]))
    report("disabled ledger",
           gate_realloc.check_ledger(off["audit"], tol["ledger_rel"],
                                     tol["ledger_abs"]))
    # The vadose-domain cost bounds (S3-A4): vloss and substeps are gated
    # HERE, where the mode is applicable; rl2's saturated box records them.
    report("disabled vloss", gate_realloc.check_vloss(off["audit"], tol["vloss_rel"]))
    report("disabled substep cost",
           gate_realloc.check_substep_ratio(off["substeps"], results["legacy"]["substeps"],
                                            tol["substep_ratio"]))
    # Legacy lane: recorded, not gated (the restore's effect on a steady front).
    d = gate_realloc.deltas(results["legacy"]["audit"])
    churn = d[gate_realloc.REALLOC] + 2.0 * d[gate_realloc.REALLOC_DROPPED]
    _, legacy_uniform = gate_realloc.check_uniform_flux(results["legacy"]["qz"],
                                                        -2.0e-6, tol["uniform_rel"])
    print(f"  legacy (recorded): gross churn {churn:.3e} m^3; {legacy_uniform}")
    print("rl3 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_rl4(args: argparse.Namespace) -> int:
    tol = rl_tolerances("rl4")
    results: dict[str, dict] = {}
    for lane in ("legacy", "disabled"):
        def consistent(doc):
            doc.setdefault("transport", {})["subsurface_update"] = "consistent"
        config = rl_stage_lane(args, "x-orient", "coupled-salt-base.yaml", lane,
                               extra_transform=consistent)
        if config is None:
            print(f"rl4 gate: {CAPABILITY_ABSENT_S3RL}")
            print("rl4 gate: FAIL")
            return 1
        run_case(args.frehg, config, config.parent, args.mpiexec, args.ranks)
        output = config.parent / "out" / "output.h5"
        with h5py.File(output, "r") as handle:
            gw_audit = handle["/monitor/gw_mass_audit"][:]
            tr_audit = handle["/monitor/transport_audit"][:]
            s_max = 0.0
            for group in ("/transport/concentration", "/transport/concentration_surface"):
                if group in handle:
                    for t in handle[group]:
                        s_max = max(s_max, float(np.nanmax(handle[f"{group}/{t}"][:])))
        results[lane] = {"gw": gw_audit, "tr": tr_audit, "s_max": s_max}
    ok = True
    good, msg = gate_realloc.check_anchor(results["disabled"]["tr"], results["legacy"]["tr"],
                                          results["disabled"]["s_max"],
                                          results["disabled"]["gw"], tol["anchor_rel"],
                                          tol["anchor_improvement"],
                                          tol["anchor_legacy_floor_rel"])
    print(f"  disabled vs legacy subs_anchor: {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    print("rl4 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1



# ---------------------------------------------------------------------------
# Stage-3 T2: the open-boundary-closure gates ob1/ob3/ob5 (Stage-3 plan §6).
# Pure references and checks live in gate_openbc.py (the discrete steady
# chain); tolerances in tolerances/stage3-ob.yaml; §6.3 battery in
# scripts/test_ob_gates.py.
# ---------------------------------------------------------------------------

import gate_openbc  # noqa: E402

CAPABILITY_ABSENT_S3OB = ("capability absent: the schema rejects the T2 "
                          "target key (boundary_conditions[].inflow)")


def ob_tolerances(gate: str) -> dict:
    with open(HERE / "tolerances" / "stage3-ob.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)[gate]


def ob_stage_lane(args: argparse.Namespace, case_yaml: str, lane: str,
                  bc_name: str, u_ext: float | None = None) -> Path | None:
    """Stage ob-openbc into work/<lane> and set the named BC's inflow block:
    lane 'external' / 'zero_gradient', optionally with a constant u_ext."""
    case_dir = stage_case(args.repo, "ob-openbc", args.work / lane)
    config = case_dir / case_yaml

    def set_lane(doc):
        for bc in doc["boundary_conditions"]:
            if bc["name"] == bc_name:
                bc["inflow"] = {"closure": "external" if lane != "zero_gradient"
                                else "zero_gradient"}
                if u_ext is not None and lane != "zero_gradient":
                    bc["inflow"]["velocity"] = {"u": {"constant": u_ext}}

    rewrite_config(config, set_lane)
    if not validate_config(args.frehg, config):
        return None
    return config


def gate_ob1(args: argparse.Namespace) -> int:
    tol = ob_tolerances("ob1")
    ob1 = dict(nx=24, dx=1.0, dy=1.0, dt=0.25, bed=0.0, n_manning=0.01,
               eta_west=1.005)
    lanes = {"external": 0.0, "ext-moving": 0.1, "zero_gradient": None}
    results: dict[str, dict] = {}
    for lane, u_ext in lanes.items():
        config = ob_stage_lane(args, "ob1-channel.yaml", lane, "west-stage", u_ext)
        if config is None:
            print(f"ob1 gate: {CAPABILITY_ABSENT_S3OB}")
            print("ob1 gate: FAIL")
            return 1
        # Strict solver: the gated identity is the momentum chain, so the
        # free-surface solve must sit at machine precision (the default
        # bjacoby-icc rtol 1e-8 leaves a solver-class residual).
        run_case(args.frehg, config, config.parent, args.mpiexec, args.ranks,
                 STRICT_PETSC_OPTIONS)
        with h5py.File(config.parent / "out" / "output.h5", "r") as handle:
            times = output_times(handle, "/surface/uu")
            # Flat output index k holds uu(j, k+1) — face k+1 (the east face
            # of cell k+1); faces 1..nx-1 therefore sit at indices 0..nx-2.
            uu_last = handle[f"/surface/uu/{times[-1]}"][:].reshape(1, 24)
            eta_last = handle[f"/surface/eta/{times[-1]}"][:].reshape(1, 24)
            tail = times[int(0.9 * (len(times) - 1)):]
            u_in_tail = np.array([
                float(handle[f"/surface/uu/{t}"][:].reshape(1, 24)[0, 0])
                for t in tail])
        results[lane] = {"uu": uu_last[0], "eta": eta_last[0], "tail": u_in_tail}
    ok = True

    def report(label, result):
        nonlocal ok
        good, msg = result
        print(f"  {label}: {msg}{'' if good else ' FAIL'}")
        ok = ok and good

    chains = {}
    for lane, u_ext in (("external", 0.0), ("ext-moving", 0.1)):
        chain = gate_openbc.Chain(closure="external", u_ext=u_ext, **ob1)
        q, eta_chain, u_chain = gate_openbc.chain_solve(chain, 1.0)
        chains[lane] = (q, eta_chain, u_chain)
        # Model faces: uu field index i is the east face of cell i+1
        # (0-based flat output), so faces 1..nx-1 sit at indices 1..nx-1.
        model_u = results[lane]["uu"][0:ob1["nx"] - 1]
        report(f"{lane} velocity profile",
               gate_openbc.check_velocity_profile(model_u, u_chain[1:],
                                                  tol["velocity_rel"]))
        report(f"{lane} steadiness",
               gate_openbc.check_steadiness(results[lane]["tail"],
                                            tol["steady_abs"]))
        ratio_chain = ((eta_chain[1] - eta_chain[2])
                       / (u_chain[1]**2 / (2.0 * gate_openbc.GRAVITY)))
        report(f"{lane} entrance ratio",
               gate_openbc.check_entrance_ratio(
                   1.005, float(results[lane]["eta"][1]),
                   float(results[lane]["uu"][0]), ratio_chain,
                   tol["ratio_abs"]))
    u_a = float(results["external"]["uu"][0])
    u_c = float(results["zero_gradient"]["uu"][0])
    r_a = (1.005 - float(results["external"]["eta"][1])) / (
        u_a**2 / (2.0 * gate_openbc.GRAVITY))
    r_c = (1.005 - float(results["zero_gradient"]["eta"][1])) / (
        u_c**2 / (2.0 * gate_openbc.GRAVITY))
    report("zero_gradient contrast (the real negative lane)",
           gate_openbc.check_lane_contrast(u_a, u_c, tol["contrast_factor"],
                                           r_a, r_c, tol["ratio_cap"]))
    print("ob1 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_ob3(args: argparse.Namespace) -> int:
    tol = ob_tolerances("ob3")
    results: dict[str, dict] = {}
    for lane in ("external", "zero_gradient"):
        config = ob_stage_lane(args, "ob3-tidal-basin.yaml", lane, "tide-stage")
        if config is None:
            print(f"ob3 gate: {CAPABILITY_ABSENT_S3OB}")
            print("ob3 gate: FAIL")
            return 1
        # Strict solver: the ledger identity is gated at rounding; under the
        # default rtol 1e-8 the flux/volume mismatch is the solver residual
        # (first execution measured -6.6e-4 m^3 over X = 3.0e4 m^3).
        run_case(args.frehg, config, config.parent, args.mpiexec, args.ranks,
                 STRICT_PETSC_OPTIONS)
        with h5py.File(config.parent / "out" / "output.h5", "r") as handle:
            ma = handle["/monitor/mass_audit"][:]
            fr2 = 0.0
            for t in output_times(handle, "/surface/vv"):
                vv = handle[f"/surface/vv/{t}"][:]
                depth = handle[f"/surface/depth/{t}"][:]
                speed2 = float(np.nanmax(vv * vv))
                h_min = max(float(np.nanmin(depth)), 1.0e-3)
                fr2 = max(fr2, speed2 / (gate_openbc.GRAVITY * h_min))
        results[lane] = {"audit": ma, "fr2": fr2}
    ok = True
    for lane in ("external", "zero_gradient"):
        good, msg = gate_openbc.check_surface_ledger(results[lane]["audit"],
                                                     tol["ledger_rel"])
        print(f"  {lane} ledger: {msg}{'' if good else ' FAIL'}")
        ok = ok and good
    fr2 = max(results["external"]["fr2"], results["zero_gradient"]["fr2"])
    prisms = {lane: gate_openbc.prism_per_cycle(results[lane]["audit"], 3600.0)
              for lane in results}
    good, msg = gate_openbc.check_prism_agreement(
        prisms["external"], prisms["zero_gradient"], fr2,
        tol["prism_coef"], tol["prism_floor_rel"])
    print(f"  prism agreement: {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    print("ob3 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_ob5(args: argparse.Namespace) -> int:
    tol = ob_tolerances("ob5")
    results: dict[str, dict] = {}
    for lane in ("external", "zero_gradient"):
        config = ob_stage_lane(args, "ob5-gyre.yaml", lane, "sea-stage")
        if config is None:
            print(f"ob5 gate: {CAPABILITY_ABSENT_S3OB}")
            print("ob5 gate: FAIL")
            return 1
        run_case(args.frehg, config, config.parent, args.mpiexec, args.ranks)
        with h5py.File(config.parent / "out" / "output.h5", "r") as handle:
            times = output_times(handle, "/surface/uu")
            ke = []
            influx = 0.0
            prev_t = 0.0
            for t in times:
                uu = handle[f"/surface/uu/{t}"][:].reshape(30, 20)
                vv = handle[f"/surface/vv/{t}"][:].reshape(30, 20)
                depth = handle[f"/surface/depth/{t}"][:].reshape(30, 20)
                ke.append(gate_openbc.domain_ke(uu, vv, depth, 1.0, 1.0))
                # Faces into the stage row: vv at interior row 28 (0-based),
                # the north face of the last interior row.
                influx += gate_openbc.stage_ke_influx(vv[28, :], depth[28, :],
                                                      1.0, t - prev_t)
                prev_t = t
        results[lane] = {"ke": np.array(ke), "influx": influx}
    ok = True
    ext = results["external"]
    good, msg = gate_openbc.check_ke_decay(ext["ke"], tol["settle_frames"],
                                           tol["rise_rel"])
    print(f"  external KE decay: {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    good, msg = gate_openbc.check_ke_influx(ext["influx"], float(ext["ke"][0]),
                                            tol["influx_rel"])
    print(f"  external stage KE influx: {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    good, msg = gate_openbc.check_influx_contrast(
        results["zero_gradient"]["influx"], ext["influx"],
        tol["contrast_factor"], tol["influx_floor_j"])
    print(f"  zero_gradient contrast: {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    print("ob5 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1



# ---------------------------------------------------------------------------
# Stage-3 T3: the periodic-topology gates pd1-pd8 (Stage-3 plan §7).
# Pure references and checks live in gate_periodic.py; tolerances in
# tolerances/stage3-pd.yaml; §6.3 battery in scripts/test_pd_gates.py.
# ---------------------------------------------------------------------------

import gate_periodic  # noqa: E402

CAPABILITY_ABSENT_S3PD = ("capability absent: the schema rejects the T3 "
                          "target key (domain.periodic)")


def pd_tolerances(gate: str) -> dict:
    with open(HERE / "tolerances" / "stage3-pd.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)[gate]


def pd_transpose_doc(doc: dict) -> dict:
    """The pd1 transposed lane: swap nx/ny, dx/dy, periodic axes, wind
    components, wall polygons (x <-> y), and monitor indices."""
    doc = yaml.safe_load(yaml.safe_dump(doc))
    dom = doc["domain"]
    dom["nx"], dom["ny"] = dom["ny"], dom["nx"]
    dom["dx"], dom["dy"] = dom["dy"], dom["dx"]
    per = dom.get("periodic", {})
    dom["periodic"] = {"x": bool(per.get("y", False)), "y": bool(per.get("x", False))}
    wind = doc["surface_water"]["wind"]
    wind["u10"], wind["v10"] = wind["v10"], wind["u10"]
    for bc in doc.get("boundary_conditions", []):
        bc["region"]["polygon"] = [[pt[1], pt[0]] for pt in bc["region"]["polygon"]]
    for mon in doc.get("output", {}).get("monitors", []):
        mon["i"], mon["j"] = mon["j"], mon["i"]
    doc["simulation"]["id"] += "-T"
    return doc


def gate_pd_restart(args: argparse.Namespace) -> int:
    """Stage-3 restart gate, pd lane (plan §3.4, S3-A6): a periodic surface
    run (pd2-surface) checkpoints mid-run and restarts to the project's
    restart class — flow fields to 1e-12 ABSOLUTE (measured 1.9e-15; the
    wrapped halos and seam advection carry no restart-order state), the
    superbee tracer to its measured branch-amplification class (the walled
    control of the same case measures LARGER deviations on every field)."""
    tol = pd_tolerances("pd-restart")
    case_dir = stage_case(args.repo, "pd-periodic", args.work)
    config = case_dir / "pd2-surface.yaml"

    def full_setup(doc):
        doc["output"]["checkpoint"] = {"interval": 300}

    rewrite_config(config, full_setup)
    if not validate_config(args.frehg, config):
        print(f"pd restart gate: {CAPABILITY_ABSENT_S3PD}")
        print("pd restart determinism: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"

    restart_dir = stage_case(args.repo, "pd-periodic", args.work / "restart")
    restart_config = restart_dir / "pd2-surface.yaml"

    def add_restart(doc):
        full_setup(doc)
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()),
                          "time": 300}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for group in ("/surface/eta", "/surface/uu", "/surface/vv",
                      "/surface/depth", "/transport/concentration_surface"):
            tracer = "concentration" in group
            for t in output_times(part, group):
                if t <= 300:
                    continue
                a = full[f"{group}/{t}"][:]
                b = part[f"{group}/{t}"][:]
                good, msg = gate_periodic.check_restart_field(
                    a, b, f"{group.rsplit('/', 1)[-1]}@{t}",
                    tol_abs=None if tracer else tol["flow_abs"],
                    tol_rel=tol["tracer_rel"] if tracer else None)
                print(f"  restart {msg}{'' if good else ' FAIL'}")
                ok = ok and good
    print("pd restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_pd1(args: argparse.Namespace) -> int:
    tol = pd_tolerances("pd1")
    u_ref = gate_periodic.discrete_balance(10.0, 0.5, 0.05, 2.0)
    u_cont = gate_periodic.continuum_balance(10.0, 0.5, 0.05)
    print(f"  discrete balance {u_ref:.6f} m/s (continuum {u_cont:.6f}, recorded)")
    ok = True

    def report(label, result):
        nonlocal ok
        good, msg = result
        print(f"  {label}: {msg}{'' if good else ' FAIL'}")
        ok = ok and good

    for lane in ("strict", "default", "transposed"):
        case_dir = stage_case(args.repo, "pd-periodic", args.work / lane)
        config = case_dir / "pd1-alongshore.yaml"
        transposed = (lane == "transposed")
        if transposed:
            with open(config, encoding="utf-8") as handle:
                doc = yaml.safe_load(handle)
            with open(config, "w", encoding="utf-8") as handle:
                yaml.safe_dump(pd_transpose_doc(doc), handle, sort_keys=False)
        if not validate_config(args.frehg, config):
            print(f"pd1 gate: {CAPABILITY_ABSENT_S3PD}")
            print("pd1 gate: FAIL")
            return 1
        extra = [] if lane == "default" else STRICT_PETSC_OPTIONS
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks, extra)
        shape = (8, 4) if transposed else (4, 8)  # (ny, nx)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            times = output_times(handle, "/surface/eta")
            eta = handle[f"/surface/eta/{times[-1]}"][:].reshape(shape)
            vel_name = "vv" if transposed else "uu"
            vel = handle[f"/surface/{vel_name}/{times[-1]}"][:].reshape(shape)
            tail = times[-5:]
            mon = np.array([float(np.nanmax(np.abs(
                handle[f"/surface/{vel_name}/{t}"][:]))) for t in tail])
        # The along-shore faces: every face along the periodic axis (walls
        # make the cross-axis velocity ~0; the wall rows' along-axis faces
        # carry the same balance).
        if lane == "default":
            report("default-solver eta flatness (recorded class)",
                   gate_periodic.check_flat(eta, tol["default_eta_flat_abs"]))
            continue
        label = "strict" if lane == "strict" else "transposed (periodic y)"
        report(f"{label} eta flatness", gate_periodic.check_flat(eta, tol["eta_flat_abs"]))
        report(f"{label} u uniformity",
               gate_periodic.check_flat(vel, tol["u_flat_abs"], "u"))
        report(f"{label} balance",
               gate_periodic.check_balance(vel, u_ref, tol["balance_rel"]))
        report(f"{label} steadiness",
               gate_periodic.check_steady(mon, tol["steady_abs"]))
    print("pd1 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def pd_roll_inputs(case_dir: Path, config_name: str, k: int) -> Path:
    """Stage a rolled variant: every x-dependent raster rolled k columns and
    the config pointed at the rolled files."""
    config = case_dir / config_name
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    nx = int(doc["domain"]["nx"])
    ny = int(doc["domain"]["ny"])
    nz = int(doc["domain"]["nz"])

    def roll_raster(rel: str, is3d: bool) -> str:
        src = case_dir / rel
        values = [line for line in src.read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith("#")]
        flat = " ".join(values).split()
        if is3d:
            arr = np.array(flat).reshape(ny, nx, nz)
            arr = np.roll(arr, k, axis=1)
        else:
            arr = np.array(flat).reshape(ny, nx)
            arr = np.roll(arr, k, axis=1)
        out = case_dir / "input" / f"rolled_{k}_{Path(rel).name}"
        with open(out, "w", encoding="utf-8") as handle:
            handle.write(f"# rolled by {k} columns for the pd2 lane\n")
            for v in arr.reshape(-1):
                handle.write(f"{v}\n")
        return str(out.relative_to(case_dir))

    dom = doc["domain"]["bottom_elevation"]
    if isinstance(dom, dict) and "file" in dom:
        dom["file"] = roll_raster(dom["file"], False)
    soil = doc.get("soil", {}).get("map", {})
    if isinstance(soil, dict) and "file" in soil:
        soil["file"] = roll_raster(soil["file"], True)
    tr_ic = doc.get("initial_conditions", {}).get("transport", {})
    if isinstance(tr_ic.get("surface"), dict) and "file" in tr_ic["surface"]:
        tr_ic["surface"]["file"] = roll_raster(tr_ic["surface"]["file"], False)
    doc["simulation"]["id"] += f"-roll{k}"
    doc["output"]["filename"] = f"out-roll{k}/output.h5"
    rolled = case_dir / f"{Path(config_name).stem}-roll{k}.yaml"
    with open(rolled, "w", encoding="utf-8") as handle:
        yaml.safe_dump(doc, handle, sort_keys=False)
    return rolled


PD2_FIELDS = {
    "pd2-surface.yaml": [("/surface/eta", 1), ("/surface/uu", 1),
                         ("/surface/vv", 1), ("/surface/depth", 1),
                         ("/transport/concentration_surface", 1)],
    "pd2-gw.yaml": [("/groundwater/hydraulic_head", 8),
                    ("/groundwater/water_content", 8),
                    ("/groundwater/qy", 8), ("/groundwater/qz", 8)],
    "pd2-coupled.yaml": [("/surface/eta", 1), ("/surface/uu", 1),
                         ("/surface/vv", 1),
                         ("/groundwater/hydraulic_head", 6),
                         ("/groundwater/water_content", 6),
                         ("/transport/concentration", 6),
                         ("/transport/concentration_surface", 1)],
}


def pd2_strict_options(config_name: str) -> list[str]:
    if "coupled" in config_name:
        return STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW
    return (STRICT_PETSC_OPTIONS if "surface" in config_name
            else STRICT_PETSC_OPTIONS_GW)


def pd2_roll_bound(tol: dict, config_name: str, group: str) -> float:
    """S3-A6 per-field-class bounds: the branch-free gw lane holds the
    authored strict bound; the surface and coupled lanes' branches amplify
    the rolled solver-ordering rounding (tracers at the tracer class,
    everything else at the surface class — the coupled subsurface fields
    inherit the surface class through the seepage exchange)."""
    if "pd2-gw" in config_name:
        return tol["roll_abs"]
    if "concentration" in group:
        return tol["roll_tracer_abs"]
    return tol["roll_surface_abs"]


def gate_pd2(args: argparse.Namespace) -> int:
    tol = pd_tolerances("pd2")
    ok = True
    lanes = PD2_FIELDS
    if getattr(args, "pd2_lane", None):
        lanes = {name: fields for name, fields in PD2_FIELDS.items()
                 if args.pd2_lane in name}
    for config_name, fields in lanes.items():
        case_dir = stage_case(args.repo, "pd-periodic",
                              args.work / Path(config_name).stem)
        config = case_dir / config_name
        if not validate_config(args.frehg, config):
            print(f"pd2 gate: {CAPABILITY_ABSENT_S3PD}")
            print("pd2 gate: FAIL")
            return 1
        strict = pd2_strict_options(config_name)
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks, strict)
        base_out = case_dir / "out" / "output.h5"
        with h5py.File(base_out, "r") as handle:
            nx = 12
            doc_shape = None
            base = {}
            t_last = None
            for group, nz in fields:
                t_last = output_times(handle, group)[-1]
                data = handle[f"{group}/{t_last}"][:]
                ny = data.size // (nx * nz)
                shape = (ny, nx, nz) if nz > 1 else (ny, nx)
                base[group] = (data.reshape(ny, nx, nz) if nz > 1
                               else data.reshape(ny, nx))
            if "pd2-surface" in config_name:
                tracer_frames = [
                    handle[f"/transport/concentration_surface/{t}"][:].reshape(6, 12)
                    for t in output_times(handle, "/transport/concentration_surface")]
            if "coupled" in config_name:
                c_times = output_times(handle, "/transport/concentration")
                c_first = handle[f"/transport/concentration/{c_times[0]}"][:]
                c_last = handle[f"/transport/concentration/{c_times[-1]}"][:]
        if "pd2-surface" in config_name:
            good, msg = gate_periodic.check_seam_crossing(tracer_frames, 0,
                                                          tol["seam_floor"])
            print(f"  [{config_name}] {msg}{'' if good else ' FAIL'}")
            ok = ok and good
        if "coupled" in config_name:
            good, msg = gate_periodic.check_coupled_exchange(
                c_first, c_last, tol["coupled_exchange_floor"])
            print(f"  [{config_name}] {msg}{'' if good else ' FAIL'}")
            ok = ok and good
        for k in (3, 6):
            rolled_cfg = pd_roll_inputs(case_dir, config_name, k)
            run_case(args.frehg, rolled_cfg, case_dir, args.mpiexec, args.ranks,
                     strict)
            with h5py.File(case_dir / f"out-roll{k}" / "output.h5", "r") as handle:
                for group, nz in fields:
                    t_last = output_times(handle, group)[-1]
                    data = handle[f"{group}/{t_last}"][:]
                    ny = data.size // (12 * nz)
                    arr = (data.reshape(ny, 12, nz) if nz > 1
                           else data.reshape(ny, 12))
                    unrolled = np.roll(arr, -k, axis=1)
                    bound = pd2_roll_bound(tol, config_name, group)
                    good, msg = gate_periodic.check_roll(
                        base[group], unrolled, k, bound,
                        f"[{config_name}] {group.rsplit('/', 1)[-1]}")
                    print(f"  {msg}{'' if good else ' FAIL'}")
                    ok = ok and good
    print("pd2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_pd3(args: argparse.Namespace) -> int:
    tol = pd_tolerances("pd3")
    ok = True
    for config_name, fields in PD2_FIELDS.items():
        strict = pd2_strict_options(config_name)
        outputs = {}
        # The coupled lane adds a 2 x 2 decomposition (plan §7 pd3): both
        # wrapped-rank interfaces plus a y split through the tide rows.
        decomps = [(1, 1), (2, 1), (4, 1)]
        if "coupled" in config_name:
            decomps.append((2, 2))
        for px, py in decomps:
            n = px * py
            label = f"{Path(config_name).stem}-n{px}x{py}"
            case_dir = stage_case(args.repo, "pd-periodic", args.work / label)
            config = case_dir / config_name

            def set_decomp(doc):
                doc["domain"]["decomposition"] = {"mpi_nx": px, "mpi_ny": py}

            rewrite_config(config, set_decomp)
            if not validate_config(args.frehg, config):
                print(f"pd3 gate: {CAPABILITY_ABSENT_S3PD}")
                print("pd3 gate: FAIL")
                return 1
            run_case(args.frehg, config, case_dir, args.mpiexec, n, strict)
            outputs[(px, py)] = case_dir / "out" / "output.h5"
        with h5py.File(outputs[(1, 1)], "r") as base_h:
            for px, py in decomps[1:]:
                with h5py.File(outputs[(px, py)], "r") as other:
                    for group, _ in fields:
                        worst = 0.0
                        for t in output_times(base_h, group):
                            a = base_h[f"{group}/{t}"][:]
                            b = other[f"{group}/{t}"][:]
                            worst = max(worst,
                                        float(np.nanmax(np.abs(a - b))))
                        # S3-A6 per-field-class bounds (the pd2 classes):
                        # the gw lane holds strict 1e-12; the surface and
                        # coupled lanes' branches amplify the rank-dependent
                        # allreduce rounding (walled control measures ~10x
                        # worse than the periodic lane on every field).
                        if "pd2-gw" in config_name:
                            bound = tol["rank_abs"]
                        elif "concentration" in group:
                            bound = tol["rank_tracer_abs"]
                        else:
                            bound = tol["rank_surface_abs"]
                        good = worst <= bound
                        print(f"  [{config_name}] {group.rsplit('/', 1)[-1]} "
                              f"n=1 vs {px}x{py} (self-wrap vs wrapped "
                              f"ranks): max |delta| = {worst:.3e} "
                              f"(<= {bound:.0e}){'' if good else ' FAIL'}")
                        ok = ok and good
    print("pd3 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_pd4(args: argparse.Namespace) -> int:
    tol = pd_tolerances("pd4")
    case_dir = stage_case(args.repo, "pd-periodic", args.work)
    config = case_dir / "pd4-tidal-periodic.yaml"
    if not validate_config(args.frehg, config):
        print(f"pd4 gate: {CAPABILITY_ABSENT_S3PD}")
        print("pd4 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks,
             STRICT_PETSC_OPTIONS)
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/mass_audit"][:]
    good, msg = gate_periodic.check_surface_ledger(audit, tol["ledger_rel"])
    print(f"  periodic-x tidal ledger: {msg}{'' if good else ' FAIL'}")
    print("pd4 gate:", "PASS" if good else "FAIL")
    return 0 if good else 1


def gate_pd5(args: argparse.Namespace) -> int:
    tol = pd_tolerances("pd5")
    ok = True
    results = {}
    for lane in ("uniform", "columns"):
        case_dir = stage_case(args.repo, "pd-periodic", args.work / lane)
        config = case_dir / "pd5-confined.yaml"
        if lane == "columns":
            def vary_k(doc):
                # Periodic-by-construction x-column K variation: 3 types
                # cycling over the 6 columns.
                types = []
                for n, ks in enumerate((8.0e-6, 1.0e-5, 1.2e-5)):
                    types.append({"name": f"a{n}", "ksx": ks, "ksy": ks,
                                  "ksz": ks, "theta_s": 0.41, "theta_r": 0.05,
                                  "vg_alpha": 2.0, "vg_n": 1.8})
                doc["soil"]["types"] = types
                ids = []
                for j in range(10):
                    for i in range(6):
                        ids.extend([i % 3] * 5)
                raster = case_dir / "input" / "pd5_columns.dat"
                with open(raster, "w", encoding="utf-8") as handle:
                    handle.write("# pd5 periodic column ids\n")
                    for v in ids:
                        handle.write(f"{v}\n")
                doc["soil"]["map"] = {"file": "input/pd5_columns.dat"}
            rewrite_config(config, vary_k)
        if not validate_config(args.frehg, config):
            print(f"pd5 gate: {CAPABILITY_ABSENT_S3PD}")
            print("pd5 gate: FAIL")
            return 1
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks,
                 STRICT_PETSC_OPTIONS_GW)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            t = output_times(handle, "/groundwater/hydraulic_head")[-1]
            head = handle[f"/groundwater/hydraulic_head/{t}"][:].reshape(10, 6, 5)
            qx = handle[f"/groundwater/qx/{t}"][:]
            qy = handle[f"/groundwater/qy/{t}"][:]
        results[lane] = head
        good, msg = gate_periodic.check_linear_head(head, np.arange(10) + 0.5,
                                                    3.0, 2.8, tol["head_abs"])
        print(f"  [{lane}] {msg}{'' if good else ' FAIL'}")
        ok = ok and good
        qy_scale = float(np.nanmax(np.abs(qy)))
        good, msg = gate_periodic.check_zero_x_flux(qx, qy_scale, tol["qx_rel"])
        print(f"  [{lane}] {msg}{'' if good else ' FAIL'}")
        ok = ok and good
    worst = float(np.nanmax(np.abs(results["uniform"] - results["columns"])))
    good = worst <= tol["head_abs"]
    print(f"  uniform vs column-K head: max |delta| = {worst:.3e} "
          f"(<= {tol['head_abs']:.0e}){'' if good else ' FAIL'}")
    ok = ok and good
    print("pd5 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_pd6(args: argparse.Namespace) -> int:
    ok = True
    case_dir = stage_case(args.repo, "pd-periodic", args.work)
    # (a) x-edge-only BC under periodic x: setup fatal with the §13.B text.
    config = case_dir / "pd1-alongshore.yaml"

    def add_x_wall(doc):
        # Column 0's INTERIOR rows only (centres x = 5, y = 15/25): every
        # selected face is on the periodic x axis — the y-edge corner faces
        # would otherwise legitimately attach and defeat the probe.
        doc["boundary_conditions"].append({
            "name": "x-wall",
            "region": {"polygon": [[-1.0, 10.5], [5.5, 10.5], [5.5, 29.5],
                                    [-1.0, 29.5]]},
            "target": "surface", "kind": "velocity",
            "value": {"constant": 0.0}})

    rewrite_config(config, add_x_wall)
    if not validate_config(args.frehg, config):
        print(f"pd6 gate: {CAPABILITY_ABSENT_S3PD}")
        print("pd6 gate: FAIL")
        return 1
    result = subprocess.run(
        [str(args.mpiexec), "-n", "1", str(args.frehg), config.name],
        cwd=case_dir, capture_output=True, text=True, check=False)
    fatal = (result.returncode != 0 and
             "selects only faces on the periodic" in
             (result.stdout + result.stderr))
    print(f"  x-edge wall under periodic x: exit {result.returncode}, "
          f"message {'found' if fatal else 'MISSING'}"
          f"{'' if fatal else ' FAIL'}")
    ok = ok and fatal
    # (a2, Stage-3 T8 W8.7) the y-axis twin: pd1 flipped to periodic y,
    # with a probe BC on row 0's INTERIOR columns only (centres x = 15..65)
    # — its y-edge faces are periodic-skipped and, unlike the full-span
    # walls, it has no legitimate x-corner faces to attach (the same trap
    # the x probe's polygon avoids). Same fatal, axis-correct message.
    y_dir = stage_case(args.repo, "pd-periodic", args.work / "y-axis")
    config = y_dir / "pd1-alongshore.yaml"

    def flip_axis(doc):
        doc["domain"]["periodic"] = {"y": True}
        doc["boundary_conditions"].append({
            "name": "y-wall",
            "region": {"polygon": [[10.5, -1.0], [69.5, -1.0], [69.5, 5.5],
                                    [10.5, 5.5]]},
            "target": "surface", "kind": "velocity",
            "value": {"constant": 0.0}})
        doc["simulation"]["id"] += "-ywall"

    rewrite_config(config, flip_axis)
    result = subprocess.run(
        [str(args.mpiexec), "-n", "1", str(args.frehg), config.name],
        cwd=y_dir, capture_output=True, text=True, check=False)
    fatal = (result.returncode != 0 and
             "selects only faces on the periodic y" in
             (result.stdout + result.stderr))
    print(f"  y-edge wall under periodic y: exit {result.returncode}, "
          f"message {'found' if fatal else 'MISSING'}"
          f"{'' if fatal else ' FAIL'}")
    ok = ok and fatal
    # (b) seam warning on a deliberately non-periodic bed.
    warn_dir = stage_case(args.repo, "pd-periodic", args.work / "seam-warning")
    config = warn_dir / "pd2-surface.yaml"

    def break_seam(doc):
        doc["domain"]["bottom_elevation"] = {"file": "input/pd6_ramp.dat"}
        doc["time"]["t_end"] = 2.0
        doc["time"]["output_interval"] = 2.0
    ramp = ["# pd6 non-periodic ramp bed"]
    for j in range(6):
        for i in range(12):
            ramp.append(f"{-1.0 + 0.05 * i:.6f}")
    (warn_dir / "input" / "pd6_ramp.dat").write_text("\n".join(ramp) + "\n")
    rewrite_config(config, break_seam)
    result = subprocess.run(
        [str(args.mpiexec), "-n", "1", str(args.frehg), config.name],
        cwd=warn_dir, capture_output=True, text=True, check=False)
    warned = "the seam is a discontinuity" in (result.stdout + result.stderr)
    record = (warn_dir / "out" / "run-record.yaml")
    in_record = record.exists() and "seam is a discontinuity" in record.read_text()
    print(f"  seam warning on a ramp bed: log {'found' if warned else 'MISSING'}, "
          f"record {'found' if in_record else 'MISSING'}"
          f"{'' if warned and in_record else ' FAIL'}")
    ok = ok and warned and in_record
    print("pd6 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_pd8(args: argparse.Namespace) -> int:
    # The vadose rl5 base with periodic x through the 8-transform battery:
    # the W3.7 flag swap maps periodic-x variants onto periodic-y ones, so a
    # pass proves the implementation axis-generic (plan §7 / O3).
    case_dir = stage_case(args.repo, "rl-realloc", args.work / "stage")
    config = case_dir / "rl5-orient-base.yaml"

    def make_periodic(doc):
        doc["domain"]["periodic"] = {"x": True}
        # Drop the x-edge side heads (their faces are now interior); keep
        # the y-side flux probe for injectivity, add a y+ drain.
        doc["boundary_conditions"] = [
            bc for bc in doc["boundary_conditions"]
            if bc["name"] not in ("drain-west", "drain-east")]
        doc["boundary_conditions"].append({
            "name": "drain-north",
            "region": {"polygon": [[-0.1, 5.7], [6.1, 5.7], [6.1, 6.1],
                                    [-0.1, 6.1]]},
            "target": "groundwater_side", "kind": "head",
            "value": {"constant": -2.0}})

    rewrite_config(config, make_periodic)
    if not validate_config(args.frehg, config):
        print(f"pd8 gate: {CAPABILITY_ABSENT_S3PD}")
        print("pd8 gate: FAIL")
        return 1
    # Re-stage the modified base into the battery's expected location.
    battery_dir = args.work / "battery"
    if battery_dir.exists():
        shutil.rmtree(battery_dir)
    shutil.copytree(case_dir, battery_dir / "rl-realloc")
    args2 = argparse.Namespace(**vars(args))
    args2.work = battery_dir
    rc_gw = run_dihedral_battery(
        args2, "pd8-orient", "rl5-orient-base.yaml",
        {"/groundwater/hydraulic_head": 8, "/groundwater/water_content": 8},
        STRICT_PETSC_OPTIONS_GW, case="rl-realloc", tol_key="gw")
    # Surface lane (plan §7 pd8): the pd-periodic orient-base (periodic-x
    # torus bump, oblique wind rotated with the domain, tracer slab at the
    # seam, y walls). Scalar fields only — the dihedral push compares
    # cell-centered scalars; uu/vv are staggered vector components that
    # transform with sign flips and a component swap (the swe-orient
    # convention: a velocity asymmetry feeds eta within steps). S3-A6
    # field classes for `transpose` (the solver unknown-ordering change
    # seeds the same branch amplification the pd2 roll measures); the six
    # mixed transforms additionally swap the y walls' minus/plus legacy
    # edge arithmetic — the x-orient swe class (1e-2) applies.
    tol_surface = pd_tolerances("pd8")
    rc_sw = run_dihedral_battery(
        args, "pd8-surface", "pd8-orient-base.yaml",
        {"/surface/eta": 1, "/surface/depth": 1,
         "/transport/concentration_surface": 1},
        STRICT_PETSC_OPTIONS, case="pd-periodic",
        raster_keys=[("domain", "bottom_elevation"),
                     ("initial_conditions", "transport", "surface")],
        tol_override={
            "tol_abs_transpose": tol_surface["surface_transpose_abs"],
            "tol_abs_mixed": tol_surface["surface_mixed_abs"],
            "per_field_transpose": {
                "/transport/concentration_surface":
                    tol_surface["tracer_transpose_abs"]}})
    return rc_gw or rc_sw


# ---------------------------------------------------------------------------
# Stage-3 T4 diagnostics gates (plan §8): dg1 (interval-mean seepage rate
# vs the coupled ledger + the V3 anti-vacuity), dg2 (monitor decimation as
# a bitwise row subset, plus the checkpoint/t_end forced-row lane) and the
# dg restart lane (seepage_rate bitwise across a checkpoint).
# ---------------------------------------------------------------------------
import gate_diagnostics  # noqa: E402

CAPABILITY_ABSENT_S3DG = ("capability absent: the schema rejects the T4 "
                          "target keys (seepage_rate / output.monitor_interval)")


def dg_tolerances(gate: str) -> dict:
    with open(HERE / "tolerances" / "stage3-dg.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)[gate]


def gate_dg1(args: argparse.Namespace) -> int:
    tol = dg_tolerances("dg1")
    case_dir = stage_case(args.repo, "dg-diagnostics", args.work)
    config = case_dir / "dg1-seepage.yaml"
    if not validate_config(args.frehg, config):
        print(f"dg1 gate: {CAPABILITY_ABSENT_S3DG}")
        print("dg1 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    area = float(doc["domain"]["dx"]) * float(doc["domain"]["dy"])
    min_depth = float(doc["surface_water"]["min_depth"])
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        times = output_times(handle, "/surface/seepage_rate")
        rates = [handle[f"/surface/seepage_rate/{t}"][:] for t in times]
        snaps = [handle[f"/surface/seepage/{t}"][:] for t in times]
        depths = [handle[f"/surface/depth/{t}"][:] for t in times]
        audit = handle["/monitor/mass_audit"][:]
    # mass_audit columns (coupled): time, volume, rain, evap, outflow,
    # bc_inflow, seepage, clamped, stage_deposit — seepage is column 6.
    good, msgs = gate_diagnostics.check_rate_ledger(
        [float(t) for t in times], rates, audit[:, 0], audit[:, 6], area,
        tol["ledger_rel"])
    for msg in msgs:
        print(f"  {msg}")
    ok = ok and good
    good, msg = gate_diagnostics.check_rate_antivacuity(snaps, rates, depths,
                                                        min_depth)
    print(f"  {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    print("dg1 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def dg2_stage_lane(args: argparse.Namespace, name: str, transform) -> Path:
    case_dir = stage_case(args.repo, "dg-diagnostics", args.work / name)
    config = case_dir / "dg1-seepage.yaml"
    rewrite_config(config, transform)
    if not validate_config(args.frehg, config):
        return None
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    return case_dir / "out" / "output.h5"


def gate_dg2(args: argparse.Namespace) -> int:
    t_end = 3600.0

    def every_step(doc):
        doc["simulation"]["id"] += "-full"

    def decimated(doc):
        doc["simulation"]["id"] += "-dec60"
        doc["output"]["monitor_interval"] = 60

    def forced(doc):
        # 77 s divides neither the 1800 s checkpoint nor t_end = 3600 s:
        # the rows at 1800 and 3600 exist only through the forced-row
        # paths (checkpoint instants and t_end).
        doc["simulation"]["id"] += "-dec77"
        doc["output"]["monitor_interval"] = 77
        doc["output"]["checkpoint"] = {"interval": 1800}

    lanes = {}
    for name, transform in (("full", every_step), ("dec60", decimated),
                            ("dec77", forced)):
        out = dg2_stage_lane(args, name, transform)
        if out is None:
            print(f"dg2 gate: {CAPABILITY_ABSENT_S3DG}")
            print("dg2 gate: FAIL")
            return 1
        lanes[name] = out
    ok = True
    with h5py.File(lanes["full"], "r") as handle:
        tables = sorted(handle["/monitor"].keys())
        full = {name: handle[f"/monitor/{name}"][:] for name in tables}
    for lane in ("dec60", "dec77"):
        with h5py.File(lanes[lane], "r") as handle:
            for name in tables:
                dec = handle[f"/monitor/{name}"][:]
                good, msg = gate_diagnostics.check_decimated_subset(
                    full[name], dec, t_end, f"[{lane}] {name}")
                print(f"  {msg}{'' if good else ' FAIL'}")
                ok = ok and good
                if lane == "dec77":
                    expected = ([0.0] + [77.0 * k for k in
                                         range(1, int(t_end / 77.0) + 1)]
                                + [1800.0, t_end])
                    good, msg = gate_diagnostics.check_row_times(
                        dec[:, 0], expected, f"[{lane}] {name} schedule")
                    print(f"  {msg}{'' if good else ' FAIL'}")
                    ok = ok and good
    print("dg2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_dg_restart(args: argparse.Namespace) -> int:
    """Stage-3 restart gate, dg lane (plan §8): the interval-rate
    accumulators (acc_<name> + acc_t_prev) restart BITWISE — a seepage_rate
    frame written after the restart equals the continuous run's frame."""
    case_dir = stage_case(args.repo, "dg-diagnostics", args.work)
    config = case_dir / "dg1-seepage.yaml"

    def full_setup(doc):
        doc["output"]["checkpoint"] = {"interval": 1800}

    rewrite_config(config, full_setup)
    if not validate_config(args.frehg, config):
        print(f"dg restart gate: {CAPABILITY_ABSENT_S3DG}")
        print("dg restart determinism: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"

    restart_dir = stage_case(args.repo, "dg-diagnostics", args.work / "restart")
    restart_config = restart_dir / "dg1-seepage.yaml"

    def add_restart(doc):
        full_setup(doc)
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()),
                          "time": 1800}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(restart_config, add_restart)
    run_case(args.frehg, restart_config, restart_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for group in ("/surface/seepage_rate", "/surface/seepage",
                      "/surface/eta", "/surface/depth"):
            for t in output_times(part, group):
                if t <= 1800:
                    continue
                a = full[f"{group}/{t}"][:]
                b = part[f"{group}/{t}"][:]
                same = np.array_equal(a, b)
                if not same or t == output_times(part, group)[-1]:
                    print(f"  restart {group.rsplit('/', 1)[-1]}@{t}: "
                          f"{'bitwise' if same else 'DIFFERS'}")
                ok = ok and same
    print("dg restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Stage-3 T5 evapotranspiration gates (plan §9): et1 (new block vs the v2.0
# bulk path + the D1 coupled lane), et2 (coupled dry == uncoupled), et3 (the
# pond-to-soil switch + the v2.0 window-lag real negative), et4 (water
# balance), et7 (supply cap), et8 (salt lanes a/b/c/d), et9 (the tidal-flat
# smoke with the noET bitwise pin), et-x (orient battery + rank lane) and
# the et restart lanes.
# ---------------------------------------------------------------------------
import gate_et  # noqa: E402

CAPABILITY_ABSENT_S3ET = ("capability absent: the schema rejects the T5 "
                          "target keys (evapotranspiration / surface_dryout "
                          "/ atmosphere outputs)")


def et_tolerances(gate: str) -> dict:
    with open(HERE / "tolerances" / "stage3-et.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)[gate]


def h5_bitwise_groups(path_a: Path, path_b: Path, groups) -> tuple[bool, list]:
    msgs = []
    ok = True
    with h5py.File(path_a, "r") as a, h5py.File(path_b, "r") as b:
        for group in groups:
            if group not in a or group not in b:
                ok = False
                msgs.append(f"{group}: missing")
                continue
            for t in output_times(a, group):
                x = a[f"{group}/{t}"][:]
                y = b[f"{group}/{t}"][:]
                same = np.array_equal(np.nan_to_num(x, nan=-9.0e99),
                                      np.nan_to_num(y, nan=-9.0e99))
                if not same:
                    ok = False
                    msgs.append(f"{group}@{t}: DIFFERS (max |d| = "
                                f"{np.nanmax(np.abs(x - y)):.3e})")
    if ok:
        msgs.append(f"{len(list(groups))} groups bitwise")
    return ok, msgs


def gate_et1(args: argparse.Namespace) -> int:
    tol = et_tolerances("et1")
    ok = True
    # Lane a: the NEW uncoupled soil channel vs the v2.0 bulk path, bitwise.
    case_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "a-new")
    config = case_dir / "et1-column.yaml"
    if not validate_config(args.frehg, config):
        print(f"et1 gate: {CAPABILITY_ABSENT_S3ET}")
        print("et1 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    old_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "a-old")
    old_cfg = old_dir / "et1-column.yaml"

    def to_old(doc):
        del doc["evapotranspiration"]
        doc["groundwater"]["evaporation"] = {
            "mode": "bulk",
            "region": {"polygon": [[-1.0, -1.0], [2.0, -1.0], [2.0, 2.0],
                                    [-1.0, 2.0]]}}
        doc["simulation"]["id"] += "-v20"

    rewrite_config(old_cfg, to_old)
    run_case(args.frehg, old_cfg, old_dir, args.mpiexec, args.ranks)
    # S3-A8: the two kernels' fma contraction differs under the build's
    # default -ffp-contract=fast (bitwise through substep 15; measured
    # class in stage3-et.yaml).
    with h5py.File(case_dir / "out" / "output.h5", "r") as new_h, \
         h5py.File(old_dir / "out" / "output.h5", "r") as old_h:
        for group, bound in (("/groundwater/water_content", tol["a_theta_abs"]),
                             ("/groundwater/hydraulic_head", tol["a_head_abs"])):
            worst = 0.0
            for t in output_times(new_h, group):
                x = new_h[f"{group}/{t}"][:]
                y = old_h[f"{group}/{t}"][:]
                m = np.isfinite(x) & np.isfinite(y)
                worst = max(worst, float(np.max(np.abs(x[m] - y[m]))))
            good = worst <= bound
            print(f"  [a] new-vs-v2.0 {group.rsplit('/', 1)[-1]}: max |d| = "
                  f"{worst:.3e} (<= {bound:.0e}){'' if good else ' FAIL'}")
            ok = ok and good
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
    unreal = audit[-1][8]
    good = abs(unreal) <= 1.0e-30
    print(f"  [a] unrealised == 0 (cap never binds): {unreal:.3e}"
          f"{'' if good else ' FAIL'}")
    ok = ok and good

    # Lane b (S3-A8): the coupled D1 guard-mismatch lane.
    b_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "b-new")
    b_cfg = b_dir / "et1b-coupled.yaml"
    run_case(args.frehg, b_cfg, b_dir, args.mpiexec, args.ranks)
    bold_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "b-old")
    bold_cfg = bold_dir / "et1b-coupled.yaml"

    def b_to_old(doc):
        del doc["evapotranspiration"]
        doc["simulation"]["id"] += "-v20"
        doc.setdefault("boundary_conditions", []).append({
            "name": "evap-top",
            "region": {"polygon": [[-1.0, -1.0], [2.0, -1.0], [2.0, 2.0],
                                    [-1.0, 2.0]]},
            "target": "groundwater_top", "kind": "flux",
            "value": {"constant": 1.0e-6}})

    rewrite_config(bold_cfg, b_to_old)
    run_case(args.frehg, bold_cfg, bold_dir, args.mpiexec, args.ranks)
    with h5py.File(b_dir / "out" / "output.h5", "r") as new_h, \
         h5py.File(bold_dir / "out" / "output.h5", "r") as old_h:
        t1 = output_times(new_h, "/groundwater/hydraulic_head")[1]
        hn = new_h[f"/groundwater/hydraulic_head/{t1}"][:]
        ho = old_h[f"/groundwater/hydraulic_head/{t1}"][:]
        mask = np.isfinite(hn) & np.isfinite(ho)
        d1 = float(np.max(np.abs(hn[mask] - ho[mask])))
        audit_new = new_h["/monitor/gw_mass_audit"][:]
    good = d1 >= tol["d1_head_floor"]
    print(f"  [b] D1 head difference at t = {t1}: {d1:.3e} "
          f"(>= {tol['d1_head_floor']:.0e}){'' if good else ' FAIL'}")
    ok = ok and good
    # The new lane's ledger closes: dV = boundary_in - ss + realloc - vloss.
    dv = audit_new[-1][1] - audit_new[0][1]
    rhs = (audit_new[-1][2] - audit_new[-1][3] + audit_new[-1][4]
           - audit_new[-1][6])
    resid = abs(dv - rhs) / max(abs(audit_new[-1][7]), 1.0e-12)
    good = resid <= 1.0e-6
    print(f"  [b] new-path ledger residual / soil_evap: {resid:.3e}"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    print("et1 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et2(args: argparse.Namespace) -> int:
    tol = et_tolerances("et2")
    unc_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "uncoupled")
    unc_cfg = unc_dir / "et1-column.yaml"
    if not validate_config(args.frehg, unc_cfg):
        print(f"et2 gate: {CAPABILITY_ABSENT_S3ET}")
        print("et2 gate: FAIL")
        return 1
    run_case(args.frehg, unc_cfg, unc_dir, args.mpiexec, args.ranks,
             STRICT_PETSC_OPTIONS_GW)
    cpl_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "coupled")
    cpl_cfg = cpl_dir / "et2-coupled-dry.yaml"
    run_case(args.frehg, cpl_cfg, cpl_dir, args.mpiexec, args.ranks,
             STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW)
    ok = True
    with h5py.File(unc_dir / "out" / "output.h5", "r") as unc, \
         h5py.File(cpl_dir / "out" / "output.h5", "r") as cpl:
        worst = 0.0
        for group in ("/groundwater/hydraulic_head", "/groundwater/water_content"):
            for t in output_times(unc, group):
                a = unc[f"{group}/{t}"][:]
                b = cpl[f"{group}/{t}"][:]
                m = np.isfinite(a) & np.isfinite(b)
                scale = max(float(np.max(np.abs(a[m]))), 1.0e-12)
                worst = max(worst, float(np.max(np.abs(a[m] - b[m]))) / scale)
        good = worst <= tol["traj_rel"]
        print(f"  coupled-dry vs uncoupled trajectory: worst rel "
              f"{worst:.3e} (<= {tol['traj_rel']:.0e}){'' if good else ' FAIL'}")
        ok = ok and good
        depth_max = max(float(np.nanmax(cpl[f"/surface/depth/{t}"][:]))
                        for t in output_times(cpl, "/surface/depth"))
        good = depth_max <= 0.0
        print(f"  anti-vacuity: surface depth max = {depth_max:.3e} (== 0)"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
    print("et2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et3(args: argparse.Namespace) -> int:
    tol = et_tolerances("et3")
    area = 1.0
    h0 = 4.0e-4
    ok = True
    lanes = {}
    for lane in ("bulk", "prescribed", "v20"):
        case_dir = stage_case(args.repo, "et-evapotranspiration", args.work / lane)
        config = case_dir / "et3-pond-switch.yaml"
        if lane == "prescribed":
            def presc(doc):
                doc["evapotranspiration"]["potential"] = "prescribed"
                doc["evapotranspiration"]["prescribed"] = {"constant": 1.3e-6}
                doc["evapotranspiration"]["soil"]["limiter"] = "none"
                doc["simulation"]["id"] += "-presc"
            rewrite_config(config, presc)
        elif lane == "v20":
            def v20(doc):
                ew = gate_et.potential_bulk(
                    30.0, 101.3, gate_et.aero_resistance("liu2006", 8.0), 0.0,
                    gate_et.met_specific_humidity(25.0, 101.3, 0.3))
                del doc["evapotranspiration"]
                doc["surface_water"]["evaporation"] = {"mode": "bulk"}
                doc["boundary_conditions"] = [{
                    "name": "evap-top",
                    "region": {"polygon": [[-1.0, -1.0], [2.0, -1.0],
                                            [2.0, 2.0], [-1.0, 2.0]]},
                    "target": "groundwater_top", "kind": "flux",
                    "value": {"constant": float(ew)}}]
                doc["output"]["variables"].pop("atmosphere", None)
                doc["simulation"]["id"] += "-v20"
            rewrite_config(config, v20)
        if not validate_config(args.frehg, config):
            print(f"et3 gate [{lane}]: {CAPABILITY_ABSENT_S3ET}")
            print("et3 gate: FAIL")
            return 1
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks,
                 STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW)
        lanes[lane] = case_dir / "out" / "output.h5"

    qa = gate_et.met_specific_humidity(25.0, 101.3, 0.3)
    ew_bulk = gate_et.potential_bulk(30.0, 101.3,
                                     gate_et.aero_resistance("liu2006", 8.0),
                                     0.0, qa)
    for lane, ew in (("bulk", ew_bulk), ("prescribed", 1.3e-6)):
        with h5py.File(lanes[lane], "r") as handle:
            times = output_times(handle, "/surface/depth")
            audit = handle["/monitor/mass_audit"][:]
            depths = []
            cum_open = []
            cum_exch = []
            for t in times:
                depths.append(float(np.nanmax(handle[f"/surface/depth/{t}"][:])))
                row = audit[np.isclose(audit[:, 0], t, atol=1.0e-6)][0]
                cum_open.append(float(row[3]))   # evaporation column
                cum_exch.append(float(row[6]))   # seepage column
            good, msg = gate_et.check_pond_phase(
                [float(t) for t in times], depths, cum_open, h0, ew,
                cum_exch, area, tol["pond_rel"])
            print(f"  [{lane}] {msg}{'' if good else ' FAIL'}")
            ok = ok and good
            wet = handle[f"/atmosphere/wet_time/{times[-1]}"][:]
            soil = handle[f"/atmosphere/soil_evaporation_time/{times[-1]}"][:]
            good, msg = gate_et.check_no_double_count(wet, soil, float(times[-1]))
            print(f"  [{lane}] {msg}{'' if good else ' FAIL'}")
            ok = ok and good
            # (3) combined-rate identity from the AUDIT columns (an
            # impermeable saturated column funds the soil flux from Ss
            # storage, invisible in theta — the books are the observable),
            # asserted WHILE alpha_1 = 1 per the plan: the latest output up
            # to which the surface-extrapolated top moisture stays >= the
            # 0.375 cap mark (from the theta output; the bulk lane's top
            # cell desaturates past it near the end — measured rate decay
            # 1.6514e-6 -> 1.6413e-6 over the last two frames).
            gw_audit = handle["/monitor/gw_mass_audit"][:]
            t_ok = times[-1]
            if lane == "bulk":
                for n, t in enumerate(times):
                    th = handle[f"/groundwater/water_content/{t}"][:]
                    w_surf = 1.5 * float(th[0]) - 0.5 * float(th[1])
                    if w_surf < 0.3751:
                        t_ok = times[max(0, n - 1)]
                        break
            row_ok = int(np.nonzero(np.isclose(gw_audit[:, 0], float(t_ok),
                                               atol=1.0e-6))[0][0])
            cum_soil = float(gw_audit[row_ok][7])
            open_ok = cum_open[times.index(t_ok)]
            good, msg = gate_et.check_combined_rate(
                open_ok, cum_soil, ew, float(t_ok), area, tol["combined_rel"])
            print(f"  [{lane}] (alpha_1 = 1 through t = {t_ok}) {msg}"
                  f"{'' if good else ' FAIL'}")
            ok = ok and good
            # (4) continuity across the switch: the first full dry frame's
            # soil rate over E_w.
            t_d = h0 / ew
            full_dry = [t for t in times if float(t) - 70.0 >= t_d]
            rate = float(np.nanmax(
                handle[f"/atmosphere/soil_evaporation_rate/{full_dry[0]}"][:]))
            ratio = rate / ew
            good = abs(ratio - 1.0) <= tol["continuity_rel"]
            print(f"  [{lane}] post-switch soil rate / E_w = {ratio:.6f} "
                  f"(within {tol['continuity_rel']:.0e}){'' if good else ' FAIL'}")
            ok = ok and good
            # Anti-vacuity: >= 1 frame mid-transition.
            fracs = []
            prev = 0.0
            for t in times:
                w = float(np.nanmax(handle[f"/atmosphere/wet_time/{t}"][:]))
                fracs.append((w - prev) / 70.0 if t != times[0] else 0.0)
                prev = w
            partial = [f for f in fracs if 1.0e-6 < f < 1.0 - 1.0e-6]
            good = len(partial) >= 1
            print(f"  [{lane}] anti-vacuity: {len(partial)} partial-wet "
                  f"frame(s){'' if good else ' FAIL'}")
            ok = ok and good
    # The real negative: the v2.0 lane double counts (system volume drops
    # faster than E_w t by > v20_violation_min during/after the switch).
    with h5py.File(lanes["v20"], "r") as handle:
        times = output_times(handle, "/surface/depth")
        audit = handle["/monitor/mass_audit"][:]
        # v2.0's soil loss leaves through the exchange-accumulator strip,
        # i.e. as NEGATIVE audited seepage (the Darcy part is the ET
        # lane's measured 4e-6-relative class — negligible at 1e-3).
        combined = float(audit[-1][3]) - float(audit[-1][6])
        ratio = combined / (area * ew_bulk * float(times[-1]))
        good = abs(ratio - 1.0) > tol["v20_violation_min"]
        print(f"  [v20 negative] removal / (E_w t) = {ratio:.6f} "
              f"(must deviate > {tol['v20_violation_min']:.0e})"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
    print("et3 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et4(args: argparse.Namespace) -> int:
    tol = et_tolerances("et4")
    case_dir = stage_case(args.repo, "et-evapotranspiration", args.work)
    config = case_dir / "et4-balance.yaml"
    if not validate_config(args.frehg, config):
        print(f"et4 gate: {CAPABILITY_ABSENT_S3ET}")
        print("et4 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
    # Columns: t, vol, boundary_in, ss, realloc, dropped, vloss, soil, unreal.
    dv = audit[-1][1] - audit[0][1]
    applied = audit[-1][7]
    scale = max(abs(applied), 1.0e-12)
    # The closed column: boundary_in == -soil (the channel is the only BC).
    ident = abs(audit[-1][2] + applied) / scale
    good = ident <= tol["balance_rel"]
    print(f"  boundary_in == -soil_evaporation: rel {ident:.3e} "
          f"(<= {tol['balance_rel']:.0e}){'' if good else ' FAIL'}")
    ok = ok and good
    closure = abs(dv - (audit[-1][2] - audit[-1][3] + audit[-1][4]
                        - audit[-1][6])) / scale
    good = closure <= 1.0e-9
    print(f"  dS = boundary_in - ss + realloc - vloss: rel {closure:.3e} "
          f"(<= 1e-9){'' if good else ' FAIL'}")
    ok = ok and good
    good = abs(audit[-1][6]) == 0.0
    print(f"  vloss == 0: {audit[-1][6]:.3e}{'' if good else ' FAIL'}")
    ok = ok and good
    good = applied > 1.0e-4
    print(f"  anti-vacuity: soil channel removed {applied:.3e} m^3 "
          f"(> 1e-4){'' if good else ' FAIL'}")
    ok = ok and good
    print(f"  unrealised (recorded): {audit[-1][8]:.3e} m^3")
    print("et4 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et7(args: argparse.Namespace) -> int:
    tol = et_tolerances("et7")
    case_dir = stage_case(args.repo, "et-evapotranspiration", args.work)
    config = case_dir / "et7-cap.yaml"
    if not validate_config(args.frehg, config):
        print(f"et7 gate: {CAPABILITY_ABSENT_S3ET}")
        print("et7 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
        theta_min = min(float(np.nanmin(handle[f"/groundwater/water_content/{t}"][:]))
                        for t in output_times(handle, "/groundwater/water_content"))
    demand = 1.0e-5 * 1.0 * 3600.0
    applied = float(audit[-1][7])
    unreal = float(audit[-1][8])
    good, msg = gate_et.check_cap_trace(theta_min, 0.05, unreal, demand,
                                        applied, float(audit[-1][6]),
                                        tol["cap_abs"])
    print(f"  {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    print("et7 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et8(args: argparse.Namespace) -> int:
    tol = et_tolerances("et8")
    ok = True
    # Lane a: soil concentration without creation/destruction.
    a_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "a")
    a_cfg = a_dir / "et8a-soil-salt.yaml"
    if not validate_config(args.frehg, a_cfg):
        print(f"et8 gate: {CAPABILITY_ABSENT_S3ET}")
        print("et8 gate: FAIL")
        return 1
    run_case(args.frehg, a_cfg, a_dir, args.mpiexec, args.ranks)
    with h5py.File(a_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/transport_audit"][:]
        times = output_times(handle, "/transport/concentration")
        c0 = handle[f"/transport/concentration/{times[0]}"][:]
        th0 = handle[f"/groundwater/water_content/{times[0]}"][:]
        c1 = handle[f"/transport/concentration/{times[-1]}"][:]
        th1 = handle[f"/groundwater/water_content/{times[-1]}"][:]
    m0 = float(audit[0][2])
    good, msg = gate_et.check_conserved(float(audit[-1][2]), m0,
                                        tol["conserve_rel"], "[a] subs mass")
    print(f"  {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    good, msg = gate_et.check_cell_salt(c1, th1, c0, th0, tol["cell_rel"])
    print(f"  [a] {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    for name, col in (("subs_adjust", 9), ("subs_anchor", 10)):
        v = abs(float(audit[-1][col]))
        good = v <= tol["adjust_rel"] * m0
        print(f"  [a] {name} = {v:.3e} (<= {tol['adjust_rel']:.0e} M)"
              f"{'' if good else ' FAIL'}")
        ok = ok and good

    # Lane b: the pond store; lane c: the v2.0 delete loss (real negative,
    # staged WITHOUT the ET block — the schema forbids ET+salinity+delete
    # by construction, S3-A8).
    for lane in ("b", "c"):
        lane_dir = stage_case(args.repo, "et-evapotranspiration", args.work / lane)
        lane_cfg = lane_dir / "et8b-pond-store.yaml"
        if lane == "c":
            def to_delete(doc):
                del doc["evapotranspiration"]
                doc["surface_water"]["evaporation"] = {"mode": "bulk"}
                doc["transport"]["surface_dryout"] = "discard"
                doc["output"]["variables"]["transport"] = ["concentration_surface"]
                doc["output"]["variables"].pop("atmosphere", None)
                doc["simulation"]["id"] += "-delete"
            rewrite_config(lane_cfg, to_delete)
        run_case(args.frehg, lane_cfg, lane_dir, args.mpiexec, args.ranks)
        with h5py.File(lane_dir / "out" / "output.h5", "r") as handle:
            audit = handle["/monitor/transport_audit"][:]
            times = output_times(handle, "/transport/concentration_surface")
            s_h = []
            for t in times:
                s_v = float(np.nanmax(
                    handle[f"/transport/concentration_surface/{t}"][:]))
                h_v = float(np.nanmax(handle[f"/surface/depth/{t}"][:]))
                s_h.append((float(t), s_v, h_v))
        m_ref = 35.0 * 0.005
        if lane == "b":
            for t, s_v, h_v in s_h:
                if h_v <= 1.0e-5 or t == 0.0:
                    continue
                good, msg = gate_et.check_pond_salt(s_v, h_v, m_ref,
                                                    tol["pond_rel"] * 1.0e2,
                                                    f"[b] t={t:g}")
                if not good:
                    print(f"  {msg} FAIL")
                ok = ok and good
            # Store full at dry-out, mass conserved throughout.
            store_col = audit[:, 11]
            peak = float(np.max(store_col))
            good = abs(peak - m_ref) <= tol["conserve_rel"] * 1.0e2 * m_ref
            print(f"  [b] store peak = {peak:.6e} vs 35 h0 A = {m_ref:.6e}"
                  f"{'' if good else ' FAIL'}")
            ok = ok and good
            total = audit[:, 1] + store_col
            worst = float(np.max(np.abs(total - total[0]))) / m_ref
            good = worst <= tol["conserve_rel"] * 1.0e2
            print(f"  [b] total (mass + store) drift = {worst:.3e} rel"
                  f"{'' if good else ' FAIL'}")
            ok = ok and good
        else:
            # The v2.0 loss: the pond's salt lands in surf_adjust and the
            # total mass drops by ~35 h0 A — the gate asserts the LOSS.
            lost = float(audit[0][1] - audit[-1][1])
            good = lost > 0.9 * m_ref
            print(f"  [c negative] v2.0 delete loses {lost:.6e} of {m_ref:.6e}"
                  f"{'' if good else ' FAIL'}")
            ok = ok and good

    # Lane d: held-seepage salt.
    d_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "d")
    d_cfg = d_dir / "et8d-held-seepage.yaml"
    run_case(args.frehg, d_cfg, d_dir, args.mpiexec, args.ranks)
    with h5py.File(d_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/transport_audit"][:]
        times = output_times(handle, "/transport/salt_store")
        store_series = [float(np.nansum(handle[f"/transport/salt_store/{t}"][:]))
                        for t in times]
    # D9's claim is that NO salt is deleted: (1) the full audit identity
    # closes (every change is a booked term); (2) surf_adjust stays
    # EXACTLY zero — v2.0's deletion channel (the dry-cell discard) never
    # fires. subs_boundary (the bottom influx) and subs_anchor (the
    # theta-basis re-anchor of the upward-flow column, a measured legacy
    # property) are booked data, not losses.
    m0 = float(audit[0][1] + audit[0][2] + audit[0][11])
    lhs = float((audit[-1][1] + audit[-1][2] + audit[-1][11])) - m0
    rhs = float(audit[-1][4] + audit[-1][5] + audit[-1][6] + audit[-1][7]
                + audit[-1][8] + audit[-1][9] + audit[-1][10])
    resid = abs(lhs - rhs) / m0
    good = resid <= tol["conserve_rel"] * 1.0e2
    print(f"  [d] full audit identity residual = {resid:.3e} rel "
          f"(<= {tol['conserve_rel'] * 1.0e2:.0e}){'' if good else ' FAIL'}")
    ok = ok and good
    good = float(np.max(np.abs(audit[:, 6]))) == 0.0
    print(f"  [d] surf_adjust == 0 exactly (no deletion, the D9 observable): "
          f"max {float(np.max(np.abs(audit[:, 6]))):.3e}{'' if good else ' FAIL'}")
    ok = ok and good
    good = max(store_series) > 1.0e-4
    print(f"  [d] store holds salt while the seepage is held: peak "
          f"{max(store_series):.3e} (> 1e-4){'' if good else ' FAIL'}")
    ok = ok and good
    adj = abs(float(audit[-1][9]))
    good = adj <= tol["held_rel"] * m0
    print(f"  [d] subs_adjust = {adj:.3e} (<= {tol['held_rel']:.0e} M)"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    print("et8 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et9(args: argparse.Namespace) -> int:
    tol = et_tolerances("et9")
    ok = True
    outputs = {}
    records = {}
    for lane in ("et", "noet"):
        case_dir = stage_case(args.repo, "et-evapotranspiration", args.work / lane)
        config = case_dir / ("et9-tidal-flat.yaml" if lane == "et"
                             else "et9-noet.yaml")
        if lane == "et" and not validate_config(args.frehg, config):
            print(f"et9 gate: {CAPABILITY_ABSENT_S3ET}")
            print("et9 gate: FAIL")
            return 1
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        outputs[lane] = case_dir / "out" / "output.h5"
        with open(case_dir / "out" / "run-record.yaml", encoding="utf-8") as h:
            records[lane] = yaml.safe_load(h)
    bed = np.loadtxt(args.repo / "benchmarks" / "et-evapotranspiration" /
                     "input" / "et9_dem.dat").reshape(24, 4)
    with h5py.File(outputs["et"], "r") as handle:
        times = output_times(handle, "/surface/depth")
        t_end = float(times[-1])
        wet = handle[f"/atmosphere/wet_time/{times[-1]}"][:].reshape(24, 4)
        soil = handle[f"/atmosphere/soil_evaporation_time/{times[-1]}"][:].reshape(24, 4)
        # Zone integrals from the per-cell clocks and the rate outputs.
        sub = bed < -0.36
        supra = bed > 0.36
        inter = (bed >= -0.2) & (bed <= 0.2)
        good = float(soil[sub].max()) == 0.0
        print(f"  subtidal soil_evaporation_time max = {soil[sub].max():.3e} "
              f"(== 0){'' if good else ' FAIL'}")
        ok = ok and good
        good = float((t_end - wet[supra]).min()) >= 0.0 and \
               float(wet[supra].max()) == 0.0
        print(f"  supratidal wet_time max = {wet[supra].max():.3e} (== 0)"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
        good, msg = gate_et.check_no_double_count(wet, soil, t_end)
        print(f"  {msg}{'' if good else ' FAIL'}")
        ok = ok and good
        # Intertidal transitions: wet_time increments stop and restart.
        trans_ok = True
        wet_series = np.stack([
            handle[f"/atmosphere/wet_time/{t}"][:].reshape(24, 4)
            for t in times])
        for (j, i) in zip(*np.nonzero(inter)):
            incr = np.diff(wet_series[:, j, i])
            active = incr > 1.0e-9
            flips = int(np.sum(np.abs(np.diff(active.astype(int)))))
            if flips < 3:  # >= 1 wet->dry per cycle over 3 cycles
                trans_ok = False
        print(f"  intertidal wet/dry transitions per cell: "
              f"{'>= 1/cycle' if trans_ok else 'MISSING'}"
              f"{'' if trans_ok else ' FAIL'}")
        ok = ok and trans_ok
        # Anti-vacuity: both channels moved real volume.
        closure = records["et"].get("closure", {})
        et_cl = closure.get("evapotranspiration", {})
        vol_open = abs(float(et_cl.get("open_water_evaporation_m3", 0.0)))
        vol_soil = abs(float(et_cl.get("soil_evaporation_m3", 0.0)))
        good = vol_open > tol["zone_floor_m3"] and vol_soil > tol["zone_floor_m3"]
        print(f"  anti-vacuity: open {vol_open:.3e} / soil {vol_soil:.3e} m^3 "
              f"(> {tol['zone_floor_m3']:.0e}){'' if good else ' FAIL'}")
        ok = ok and good
        # Salt ledger: the FULL audit identity (the tide Dirichlet and the
        # side heads import/export salt legitimately through the booked
        # source/boundary columns; the et8d convention).
        audit = handle["/monitor/transport_audit"][:]
        m0 = float(audit[0][1] + audit[0][2] + audit[0][11])
        lhs = float(audit[-1][1] + audit[-1][2] + audit[-1][11]) - m0
        rhs = float(audit[-1][4] + audit[-1][5] + audit[-1][6] + audit[-1][7]
                    + audit[-1][8] + audit[-1][9] + audit[-1][10])
        drift = abs(lhs - rhs) / m0
        good = drift <= tol["ledger_salt_rel"]
        print(f"  salt audit identity residual = {drift:.3e} (<= "
              f"{tol['ledger_salt_rel']:.0e}){'' if good else ' FAIL'}")
        ok = ok and good
    # Substep-count ratio.
    def substeps(record):
        solver = record.get("solver", {}) or {}
        gws = solver.get("groundwater", {}) or {}
        return float(gws.get("solves", 0.0))
    s_et = substeps(records["et"])
    s_noet = substeps(records["noet"])
    if s_noet > 0:
        ratio = s_et / s_noet
        good = ratio <= tol["substep_ratio_max"]
        print(f"  substep ratio ET/noET = {ratio:.3f} (<= "
              f"{tol['substep_ratio_max']}){'' if good else ' FAIL'}")
        ok = ok and good
    # The noET lane is bitwise against its T5-start pin — defined only at
    # the pin's recording configuration (1 thread, the default aij LA
    # path). A p1 backend lane (--mat-type, OMP 4) changes the reduction
    # and solver algebra, so the comparison is skipped there, like
    # perf_baseline under sanitizers (S3-A11); the lane's tolerance-based
    # criteria above still gate.
    # Recorded on the T4 binary at the T5 start (g++-16 / arm64 Darwin);
    # off that toolchain h5_bitwise exits 77 and the skip is printed.
    pin = HERE / "reference-hashes" / "stage3-et9-noet.json"
    if args.mat_type is not None:
        print("  noET bitwise pin: SKIPPED under --mat-type "
              f"{args.mat_type} (recorded on the default path)")
    elif os.environ.get("FREHG_SANITIZED_BUILD"):
        # The h5_bitwise toolchain fingerprint is compiler identity only;
        # sanitizer instrumentation changes codegen (and so the chaotic
        # et9 trajectory) under the SAME compiler. run_sanitizers.sh
        # exports the flag (S3-A11).
        print("  noET bitwise pin: SKIPPED under FREHG_SANITIZED_BUILD "
              "(instrumented codegen; recorded on the default build)")
    elif pin.exists():
        result = subprocess.run(
            [sys.executable, str(args.repo / "tools" / "h5_bitwise.py"),
             "compare", str(outputs["noet"]), str(pin), "--case", "et9-noet",
             "--require-toolchain", "--allow", "monitor/mass_audit"],
            capture_output=True, text=True)
        good = result.returncode == 0
        tail = (result.stdout or result.stderr).strip().splitlines()
        print(f"  noET bitwise vs the T5-start pin: "
              f"{tail[-1] if tail else result.returncode}"
              f"{'' if good else ' FAIL'}")
        if result.returncode == 77:
            print("  (off-toolchain skip accepted)")
            good = True
        elif good:
            # Stage-3 T9 (ws5): the mass_audit table gains appended columns
            # (declared exception) — its T8 prefix must stay bitwise.
            good = mass_audit_prefix_check(outputs["noet"], "et9-noet", args,
                                           load_mass_audit_prefix_doc(args))
        ok = ok and good
    print("et9 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et_orient(args: argparse.Namespace) -> int:
    # The et-x dihedral battery on the coupled ET + store base (scalar
    # fields, the pd8/S3-A6 conventions). Classes adjudicated at first
    # green and recorded in stage3-et.yaml under a S3-A8 note if needed.
    return run_dihedral_battery(
        args, "et-orient", "orient-base.yaml",
        {"/surface/eta": 1, "/surface/depth": 1,
         "/transport/concentration_surface": 1, "/transport/salt_store": 1,
         "/groundwater/water_content": 6},
        STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW,
        case="et-evapotranspiration",
        raster_keys=[("domain", "bottom_elevation")],
        tol_override={"tol_abs_transpose": 1.0e-6,
                      "tol_abs_mixed": 1.0e-2,
                      "per_field_transpose": {
                          "/transport/concentration_surface": 1.0e-4,
                          "/transport/salt_store": 1.0e-4},
                      # S3-A8.9: T5's symmetry assertion is the TRANSPOSE
                      # class (measured exact: tracer 0.0, store 8e-20,
                      # theta 9e-15 — T5 adds no edge arithmetic). The six
                      # mixed lanes are RECORDED: the documented legacy
                      # live-east/north-edge-face asymmetry
                      # (symmetry-exemptions.md) lets the rotated pond
                      # drain its shoreline row (measured depth 5.9e-2,
                      # tracer 20 = the wet/dry status of shoreline
                      # cells), a pre-existing surface property whose own
                      # gate is the swe-orient battery.
                      "mixed_recorded": True})


def gate_et_rank(args: argparse.Namespace) -> int:
    # S3-A8: the strict rank identity runs on the trajectory-stable
    # orient base (et9's tide front is chaotic — the 1e-12 twin
    # reproduces its rank deltas exactly, so no per-field identity can
    # gate it; its ledgers/zones/complementarity stay in gate_et9).
    ok = True
    outputs = {}
    for n in (1, 2, 4):
        case_dir = stage_case(args.repo, "et-evapotranspiration",
                              args.work / f"n{n}")
        config = case_dir / "orient-base.yaml"

        def decomp(doc, n=n):
            doc["domain"]["decomposition"] = {"mpi_nx": n, "mpi_ny": 1}

        rewrite_config(config, decomp)
        if not validate_config(args.frehg, config):
            print(f"et rank gate: {CAPABILITY_ABSENT_S3ET}")
            print("et rank invariance: FAIL")
            return 1
        run_case(args.frehg, config, case_dir, args.mpiexec, n,
                 STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW)
        outputs[n] = case_dir / "out" / "output.h5"
    with h5py.File(outputs[1], "r") as base:
        for n in (2, 4):
            with h5py.File(outputs[n], "r") as other:
                for group, bound in (("/surface/eta", 1.0e-6),
                                     ("/surface/depth", 1.0e-6),
                                     ("/groundwater/water_content", 1.0e-6),
                                     ("/transport/concentration_surface", 1.0e-4),
                                     ("/transport/salt_store", 1.0e-4)):
                    worst = 0.0
                    for t in output_times(base, group):
                        a = base[f"{group}/{t}"][:]
                        b = other[f"{group}/{t}"][:]
                        m = np.isfinite(a) & np.isfinite(b)
                        worst = max(worst, float(np.max(np.abs(a[m] - b[m]))))
                    good = worst <= bound
                    print(f"  {group.rsplit('/', 1)[-1]} n=1 vs n={n}: "
                          f"max |delta| = {worst:.3e} (<= {bound:.0e})"
                          f"{'' if good else ' FAIL'}")
                    ok = ok and good
    print("et rank invariance:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et_restart(args: argparse.Namespace) -> int:
    """T5 restart lanes: et3-prescribed checkpointed inside the pond-to-
    soil transition, and et8b checkpointed while the store is full."""
    ok = True
    for name, base, ckpt, groups in (
            ("et3", "et3-pond-switch.yaml", 280,
             ["/surface/eta", "/surface/depth", "/atmosphere/wet_time",
              "/atmosphere/soil_evaporation_time",
              "/groundwater/water_content"]),
            ("et8b", "et8b-pond-store.yaml", 3600,
             ["/surface/eta", "/surface/depth",
              "/transport/concentration_surface", "/transport/salt_store"])):
        full_dir = stage_case(args.repo, "et-evapotranspiration",
                              args.work / f"{name}-full")
        config = full_dir / base

        def setup(doc, name=name, ckpt=ckpt):
            doc["output"]["checkpoint"] = {"interval": ckpt}
            if name == "et3":
                doc["evapotranspiration"]["potential"] = "prescribed"
                doc["evapotranspiration"]["prescribed"] = {"constant": 1.3e-6}
                doc["evapotranspiration"]["soil"]["limiter"] = "none"

        rewrite_config(config, setup)
        if not validate_config(args.frehg, config):
            print(f"et restart gate: {CAPABILITY_ABSENT_S3ET}")
            print("et restart determinism: FAIL")
            return 1
        run_case(args.frehg, config, full_dir, args.mpiexec, args.ranks)
        full_out = full_dir / "out" / "output.h5"
        re_dir = stage_case(args.repo, "et-evapotranspiration",
                            args.work / f"{name}-restart")
        re_cfg = re_dir / base

        def add_restart(doc, name=name, ckpt=ckpt, full_out=full_out):
            setup(doc, name=name, ckpt=ckpt)
            doc["restart"] = {"enabled": True,
                              "file": str(full_out.resolve()), "time": ckpt}
            doc["output"]["filename"] = "out/restarted.h5"

        rewrite_config(re_cfg, add_restart)
        run_case(args.frehg, re_cfg, re_dir, args.mpiexec, args.ranks)
        with h5py.File(full_out, "r") as full, \
             h5py.File(re_dir / "out" / "restarted.h5", "r") as part:
            for group in groups:
                for t in output_times(part, group):
                    if t <= ckpt:
                        continue
                    a = full[f"{group}/{t}"][:]
                    b = part[f"{group}/{t}"][:]
                    same = np.array_equal(np.nan_to_num(a, nan=-9.0e99),
                                          np.nan_to_num(b, nan=-9.0e99))
                    if not same or t == output_times(part, group)[-1]:
                        print(f"  [{name}] restart "
                              f"{group.rsplit('/', 1)[-1]}@{t}: "
                              f"{'bitwise' if same else 'DIFFERS'}")
                    ok = ok and same
    print("et restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Stage-3 T6 vegetation/transpiration gates (plan §10): et4-t (the ET water
# balance with a canopy), et5 (partition + root-zone location, with the
# vegetation-absent bitwise pin and the bare-contrast real negative), et6
# (root-fraction exactness incl. the root_uptake Field3 output), et7-t
# (Feddes stress + the T supply cap), et8-t (transpiration salt exclusion)
# and et9-veg (the vegetated tidal flat + O7) with its restart lane.
# ---------------------------------------------------------------------------

CAPABILITY_ABSENT_S3VG = ("capability absent: the schema rejects the T6 "
                          "target keys (vegetation / root_uptake / "
                          "transpiration outputs)")

ET9_MARSH = {
    "types": [{"name": "marsh", "lai": 2.5, "extinction": 0.5,
               "leaf_resistance": 100.0, "root_depth": 0.3,
               "root_profile": "linear",
               "feddes": {"h3": -4.0, "h4": -80.0}}],
    "map": {"file": "input/et9_veg.dat"},
}


def gate_et4t(args: argparse.Namespace) -> int:
    tol = et_tolerances("et4")
    case_dir = stage_case(args.repo, "et-evapotranspiration", args.work)
    config = case_dir / "et4-balance.yaml"

    def add_veg(doc):
        doc["vegetation"] = {
            "types": [{"name": "shrub", "lai": 2.0, "extinction": 0.5,
                       "leaf_resistance": 100.0, "root_depth": 0.4,
                       "root_profile": "linear",
                       "feddes": {"h3": -10.0, "h4": -150.0}}],
            "map": {"constant": "shrub"}}
        doc["output"]["variables"]["atmosphere"] = [
            "soil_evaporation_rate", "transpiration_rate",
            "transpiration_unrealised_rate"]
        doc["simulation"]["id"] += "-veg"

    rewrite_config(config, add_veg)
    if not validate_config(args.frehg, config):
        print(f"et4-t gate: {CAPABILITY_ABSENT_S3VG}")
        print("et4-t gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
    # Columns: t, vol, bin, ss, re, dr, vloss, soil, soil_un, T, T_un.
    dv = audit[-1][1] - audit[0][1]
    rhs = (audit[-1][2] - audit[-1][3] + audit[-1][4] - audit[-1][6]
           - audit[-1][9])
    scale = max(abs(audit[-1][7]) + abs(audit[-1][9]), 1.0e-12)
    resid = abs(dv - rhs) / scale
    good = resid <= 1.0e-9
    print(f"  dS = bin - ss + re - vloss - T: rel {resid:.3e} (<= 1e-9)"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    ident = abs(audit[-1][2] + audit[-1][7]) / scale
    good = ident <= tol["balance_rel"]
    print(f"  boundary_in == -soil_evaporation: rel {ident:.3e}"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    good = audit[-1][9] > 1.0e-4 and audit[-1][7] > 1.0e-5
    print(f"  anti-vacuity: T = {audit[-1][9]:.3e}, E = {audit[-1][7]:.3e} m^3"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    # The negative: dropping T from the identity must fail.
    bad = abs(dv - (rhs + audit[-1][9])) / scale
    good = bad > 1.0e-6
    print(f"  negative: identity without T misses by {bad:.3e} (> 1e-6)"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    print("et4-t gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et5(args: argparse.Namespace) -> int:
    ok = True
    # (i) the et2 column with a vegetation map of ALL ZEROS: bitwise
    # against the vegetation-absent run (the same case unmodified).
    base_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "base")
    base_cfg = base_dir / "et2-coupled-dry.yaml"
    if not validate_config(args.frehg, base_cfg):
        print(f"et5 gate: {CAPABILITY_ABSENT_S3VG}")
        print("et5 gate: FAIL")
        return 1
    run_case(args.frehg, base_cfg, base_dir, args.mpiexec, args.ranks,
             STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW)
    zero_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "zeros")
    zero_cfg = zero_dir / "et2-coupled-dry.yaml"

    def add_zero_map(doc):
        doc["vegetation"] = {
            "types": [{"name": "unused", "lai": 1.0, "extinction": 0.5,
                       "leaf_resistance": 100.0, "root_depth": 0.1,
                       "feddes": {"h3": -1.0, "h4": -10.0}}],
            "map": {"file": "input/et5_zeros.dat"}}
        doc["simulation"]["id"] += "-zeros"

    (zero_dir / "input" / "et5_zeros.dat").write_text(
        "# all-bare vegetation ids\n0\n")
    rewrite_config(zero_cfg, add_zero_map)
    run_case(args.frehg, zero_cfg, zero_dir, args.mpiexec, args.ranks,
             STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW)
    good, msgs = h5_bitwise_groups(
        base_dir / "out" / "output.h5", zero_dir / "out" / "output.h5",
        ["/groundwater/hydraulic_head", "/groundwater/water_content",
         "/surface/eta"])
    for m in msgs:
        print(f"  [i] all-zeros map vs vegetation-absent: {m}")
    ok = ok and good

    # (ii) the dense canopy: E_soil ~ 0, depletion spread over the roots.
    ii_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "canopy")
    ii_cfg = ii_dir / "et5-canopy.yaml"
    run_case(args.frehg, ii_cfg, ii_dir, args.mpiexec, args.ranks)
    with h5py.File(ii_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
        times = output_times(handle, "/groundwater/water_content")
        th0 = handle[f"/groundwater/water_content/{times[0]}"][:]
        th1 = handle[f"/groundwater/water_content/{times[-1]}"][:]
    e_soil = abs(float(audit[-1][7]))
    t_act = float(audit[-1][9])
    good = e_soil < 1.0e-12 * max(t_act, 1.0e-30) or e_soil < 1.0e-15
    print(f"  [ii] E_soil = {e_soil:.3e} vs T = {t_act:.3e} "
          f"(E < 1e-12 T){'' if good else ' FAIL'}")
    ok = ok and good
    depletion = (th0 - th1) * 0.02  # dz uniform
    total = float(np.nansum(depletion))
    top_share = float(depletion[0]) / total
    good = top_share <= 2.0 * (0.02 / 0.4)
    print(f"  [ii] top-cell share = {top_share:.4f} (<= {2.0 * 0.05:.2f})"
          f"{'' if good else ' FAIL'}")
    ok = ok and good

    # (iii) the bare contrast at (ii)'s T_pot: top-concentrated removal.
    t_pot = t_act / (86400.0 * 1.0)  # mean applied rate [m/s] over A = 1
    iii_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "bare")
    iii_cfg = iii_dir / "et5-canopy.yaml"

    def to_bare(doc, rate=t_pot):
        # S3-A9: the contrast runs 2 h, not the day — on the full horizon
        # the silt redistribution smears the top-face removal over the
        # whole column (measured: a uniform theta_r end state), erasing
        # the locality the lane contrasts. At 2 h the removal (~8 mm
        # water-equivalent) is top-concentrated before redistribution
        # catches up — the share the plan's criterion targets.
        del doc["vegetation"]
        doc["evapotranspiration"]["potential"] = "prescribed"
        doc["evapotranspiration"]["prescribed"] = {"constant": float(rate)}
        doc["evapotranspiration"]["open_water"] = False
        doc["evapotranspiration"]["soil"] = {"limiter": "none"}
        doc["time"]["t_end"] = 7200.0
        doc["time"]["output_interval"] = 7200.0
        doc["output"]["variables"]["groundwater"] = ["hydraulic_head",
                                                     "water_content"]
        doc["output"]["variables"].pop("atmosphere", None)
        doc["simulation"]["id"] += "-bare"

    rewrite_config(iii_cfg, to_bare)
    run_case(args.frehg, iii_cfg, iii_dir, args.mpiexec, args.ranks)
    with h5py.File(iii_dir / "out" / "output.h5", "r") as handle:
        times = output_times(handle, "/groundwater/water_content")
        b0 = handle[f"/groundwater/water_content/{times[0]}"][:]
        b1 = handle[f"/groundwater/water_content/{times[-1]}"][:]
    bdep = (b0 - b1) * 0.02
    btotal = float(np.nansum(bdep))
    bshare = float(bdep[0]) / btotal
    good = bshare >= 0.5
    print(f"  [iii negative] bare top-cell share = {bshare:.4f} (>= 0.5)"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    print("et5 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et6(args: argparse.Namespace) -> int:
    ok = True
    for profile in ("uniform", "linear"):
        case_dir = stage_case(args.repo, "et-evapotranspiration",
                              args.work / profile)
        config = case_dir / "et6-roots.yaml"
        if profile == "linear":
            def set_linear(doc):
                doc["vegetation"]["types"][0]["root_profile"] = "linear"
                doc["simulation"]["id"] += "-linear"
            rewrite_config(config, set_linear)
        if not validate_config(args.frehg, config):
            print(f"et6 gate: {CAPABILITY_ABSENT_S3VG}")
            print("et6 gate: FAIL")
            return 1
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            times = output_times(handle, "/groundwater/water_content")
            th0 = handle[f"/groundwater/water_content/{times[0]}"][:].reshape(1, 2, 14)
            th1 = handle[f"/groundwater/water_content/{times[-1]}"][:].reshape(1, 2, 14)
            uptake = handle[f"/groundwater/root_uptake/{times[-1]}"][:].reshape(1, 2, 14)
            audit = handle["/monitor/gw_mass_audit"][:]
        # dz-free exactness (S3-A9: the regular mesh trims partial top
        # cells, so no raster-side dz ladder is reliable): with Ss = 0
        # and K ~ 0 the corrector gives theta_new = theta_n - sink dtg
        # exactly, so per cell Delta theta_k == uptake_k * t (root_uptake
        # is the interval MEAN sink — the Field3 IntervalRates leg's
        # first consumer).
        t_pot = (1.0 - math.exp(-0.5 * 2.0)) * 2.0e-7
        t_total = 600.0
        worst = 0.0
        theta_scale = t_pot * t_total / 0.05  # the largest per-cell dtheta
        for i in range(2):
            dth = np.nan_to_num(th0[0, i] - th1[0, i])
            up = np.nan_to_num(uptake[0, i]) * t_total
            for k in range(14):
                worst = max(worst, abs(dth[k] - up[k]) / theta_scale)
        # The audit total: T_act = 2 columns x T_pot x t, nothing capped.
        # (No soil channel in et6: T is column 7 — vol..vloss, T, T_un.)
        t_act = float(audit[-1][7])
        worst = max(worst,
                    abs(t_act - 2.0 * t_pot * t_total) / (2.0 * t_pot * t_total))
        if profile == "uniform":
            # Fully-covered cells share ONE sink value T_pot / L
            # (dz-free). The zone-cutting cell sinks at its partial
            # fraction (the ladder cuts k = 4 at 0.594 of its thickness)
            # — exclude it with a strict 98% filter.
            vals = np.nan_to_num(uptake).ravel()
            vals = vals[vals > 0.98 * vals.max()]
            spread = float(vals.max() - vals.min()) / (t_pot / 0.33)
            worst = max(worst, spread)
            worst = max(worst, abs(float(vals.max()) - t_pot / 0.33) / (t_pot / 0.33))
        # S3-A9: K = 1e-15 (K = 0 is outside the vG solver) redistributes
        # at the 3e-10 class over the 10 substeps; bound 1e-9.
        good = worst <= 1.0e-9
        print(f"  [{profile}] exactness worst rel = {worst:.3e} (<= 1e-9)"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
    print("et6 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et10(args: argparse.Namespace) -> int:
    """Stage-3 T7 (plan §11): the Penman-Monteith potential at model
    level — the FAO-56 Example-19 grass column. (1) The audit T_act over
    the hour equals gate_et.pm_q_form() with the model's constants and
    the case's exact met inputs; (2) the rate reproduces FAO's published
    0.63 mm/h within the measured convention gap delta_conv (+ FAO's
    print rounding). Negative: the canopy albedo moved 0.23 -> 0.20 must
    break (1)."""
    tol = et_tolerances("et10")
    ok = True

    def run_lane(label, mutate=None):
        case_dir = stage_case(args.repo, "et-evapotranspiration",
                              args.work / label)
        config = case_dir / "et10-pm-grass.yaml"
        if mutate is not None:
            rewrite_config(config, mutate)
        if not validate_config(args.frehg, config):
            return None, None
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
        with open(config, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
            audit = handle["/monitor/gw_mass_audit"][:]
        # No soil channel: T_act is column 7 (vol..vloss, T, T_un);
        # area 1 m^2, 3600 s -> T_act [m^3] is the depth over the hour.
        return float(audit[-1][7]), doc

    def reference(doc):
        atm = doc["atmosphere"]
        veg = doc["vegetation"]["types"][0]
        p = atm["pressure"]["constant"]
        t_air = atm["air_temperature"]["constant"]
        qa = gate_et.met_specific_humidity(t_air, p,
                                           atm["relative_humidity"]["constant"])
        ra = gate_et.aero_resistance("fao56", atm["wind_speed"]["constant"])
        tau = math.exp(-veg["extinction"] * veg["lai"])
        rc = veg["leaf_resistance"] / (0.5 * veg["lai"])
        g_frac = doc["evapotranspiration"]["energy"]["ground_heat_fraction"]["canopy"]
        e = gate_et.pm_q_form(t_air, p, qa,
                              atm["shortwave_incoming"]["constant"],
                              atm["longwave_in"]["constant"],
                              veg["albedo"], g_frac, ra, rc)
        return (1.0 - tau) * e * 3600.0

    t_act, doc = run_lane("base")
    if t_act is None:
        print("et10 gate: capability absent: the schema rejects "
              "potential: penman_monteith / its energy keys")
        print("et10 gate: FAIL")
        return 1
    ref = reference(doc)
    rel = abs(t_act - ref) / ref
    good = rel <= tol["exact_rel"]
    print(f"  (1) T_act vs pm_q_form: rel {rel:.3e} "
          f"(<= {tol['exact_rel']:.0e}){'' if good else ' FAIL'}")
    ok = ok and good

    mm_h = t_act * 1000.0
    bound = tol["delta_conv_mm_h"] + tol["fao_rounding_mm_h"]
    diff = abs(mm_h - 0.63)
    good = diff <= bound
    print(f"  (2) vs FAO ET0: {mm_h:.5f} mm/h, |d| = {diff:.5f} "
          f"(<= {bound:.5f}){'' if good else ' FAIL'}")
    ok = ok and good

    # Real negative: albedo 0.23 -> 0.20 must break criterion (1)
    # against the UNCHANGED reference.
    def drop_albedo(doc):
        doc["vegetation"]["types"][0]["albedo"] = 0.20
        doc["simulation"]["id"] += "-neg"

    t_neg, _ = run_lane("neg", drop_albedo)
    rel_neg = abs(t_neg - ref) / ref
    good = rel_neg > tol["exact_rel"]
    print(f"  negative (albedo 0.20): rel {rel_neg:.3e} "
          f"(> {tol['exact_rel']:.0e}){'' if good else ' FAIL'}")
    ok = ok and good
    print("et10 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Stage-3 T9 paired water–salt store gates (plan §12A, S3-A12/O15): ws1
# (per-PR source bound on the compressed S3-A12 seepage face), ws1n (the
# owner's repro, nightly), ws2 (the water/salt ledgers, asserted inside the
# lanes), ws3 (rain/evaporation invariance and the pond closed forms), ws4
# (path coverage), ws-rank and ws-restart. ws5 (the pins) lives in
# default-bitwise and the et9 noET lane. Pure checks: gate_ws.py.
# ---------------------------------------------------------------------------
import gate_ws  # noqa: E402

CAPABILITY_ABSENT_S3WS = ("capability absent: the schema rejects the T9 target "
                          "keys (output surface: water_store / film_removal_rate)")


def ws_tolerances(gate: str) -> dict:
    with open(HERE / "tolerances" / "stage3-ws.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)[gate]


def ws_strict() -> list[str]:
    return STRICT_PETSC_OPTIONS + STRICT_PETSC_OPTIONS_GW


def ws_store_frames(handle) -> tuple[np.ndarray, np.ndarray]:
    times = output_times(handle, "/surface/water_store")
    vs = np.array([handle[f"/surface/water_store/{t}"][:].ravel() for t in times])
    ms = np.array([handle[f"/transport/salt_store/{t}"][:].ravel() for t in times])
    return np.nan_to_num(vs), np.nan_to_num(ms)


def ws_source_bound_lane(args, case_dir: Path, config: Path, tol: dict,
                         label: str, strict: bool = True) -> bool:
    """The ws1/ws1n assertions on one store lane plus its discard twin:
    (1) the source bound on every frame and monitor row; (2) no clip, CFL
    within the S3-A13 cap; (3) I2 per cell; anti-vacuity; ws2 ledgers."""
    ok = True
    extra = ws_strict() if strict else None
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks, extra)
    log = (case_dir / f"run_n{args.ranks}.log").read_text()
    out = case_dir / "out" / "output.h5"
    with open(case_dir / "out" / "run-record.yaml", encoding="utf-8") as h:
        record = yaml.safe_load(h)
    with h5py.File(out, "r") as handle:
        surf_max = [float(np.nanmax(handle[f"/transport/concentration_surface/{t}"][:]))
                    for t in output_times(handle, "/transport/concentration_surface")]
        subs_max = [float(np.nanmax(handle[f"/transport/concentration/{t}"][:]))
                    for t in output_times(handle, "/transport/concentration")]
        mon_max = []
        with open(config, encoding="utf-8") as h:
            monitors = yaml.safe_load(h)["output"].get("monitors", [])
        for mon in monitors:
            if "concentration_surface" in mon["variables"]:
                col = 1 + mon["variables"].index("concentration_surface")
                mon_max.append(float(np.nanmax(handle[f"/monitor/{mon['name']}"][:, col])))
        vs, ms = ws_store_frames(handle)
        ma = handle["/monitor/mass_audit"][:]
        ta = handle["/monitor/transport_audit"][:]
    for what, series in (("surface", surf_max), ("subsurface", subs_max),
                         ("monitor rows", mon_max)):
        good, msg = gate_ws.check_source_bound(series, tol["c_max"], tol["bound_rel"],
                                               f"{label} {what}")
        print(f"  {msg}{'' if good else ' FAIL'}")
        ok = ok and good
    cfl = [float(x) for x in re.findall(r"CFL = ([0-9.eE+-]+) exceeds 1", log)]
    clips = len(re.findall(r"clip", log, re.I))
    good = (max(cfl) if cfl else 0.0) <= tol["cfl_max"] and clips <= tol["clip_warnings"]
    print(f"  {label} CFL warnings {len(cfl)} (max {max(cfl) if cfl else 0.0:.4f} <= "
          f"{tol['cfl_max']}), clip warnings {clips}{'' if good else ' FAIL'}")
    ok = ok and good
    worst_excess = 0.0
    unpaired = 0.0
    for v, m in zip(vs, ms):
        good_f, msg = gate_ws.check_store_bound(m, v, tol["c_max"], tol["bound_rel"])
        worst_excess = max(worst_excess, float(np.max(m - tol["c_max"] * v)))
        unpaired = max(unpaired, float(np.max(np.abs(np.where(v == 0.0, m, 0.0)))))
        ok = ok and good_f
    good = worst_excess <= tol["c_max"] * float(vs.max(initial=0.0)) * 1.0e-12 and unpaired == 0.0
    print(f"  {label} I2: max (M_s - {tol['c_max']} V_s) = {worst_excess:.3e}, salt on "
          f"empty stores {unpaired:.3e} (== 0){'' if good else ' FAIL'}")
    ok = ok and good
    parked = float(vs.sum(axis=1).max(initial=0.0))
    cols = gate_ws.audit_columns(coupled=True, store=True)
    releases = int(np.sum(np.diff(ma[:, cols["store_release"]]) > 0.0))
    good = parked >= tol["parked_floor_m3"] and releases >= tol["release_rows_floor"]
    print(f"  {label} anti-vacuity: max Sigma V_s {parked:.3e} m^3 (>= "
          f"{tol['parked_floor_m3']:.1e}), release rows {releases} (>= "
          f"{tol['release_rows_floor']}){'' if good else ' FAIL'}")
    ok = ok and good
    surf = record.get("closure", {}).get("surface", {})
    gwc = record.get("closure", {}).get("groundwater", {})
    if "vent_thin_floor_m3" in tol:
        p1 = float(surf.get("film_removed_isolated_m3", 0.0))
        p4 = float(gwc.get("vent_thin_m3", 0.0))
        good = p4 >= tol["vent_thin_floor_m3"]
        print(f"  {label} path coverage: P1 isolated {p1:.3e} m^3 (recorded), P4 thin "
              f"vents {p4:.3e} m^3 (>= {tol['vent_thin_floor_m3']:.1e})"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
    tw = ws_tolerances("ws2")
    good, msg = gate_ws.check_water_identity(ma, cols, True, tw["water_rel"], tw["water_abs_m3"])
    print(f"  {label} ws2 {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    good, msg = gate_ws.check_salt_identity(ta, True, tw["salt_rel"])
    print(f"  {label} ws2 {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    return ok


def ws_discard_twin(args, case: str, config_name: str, workdir: Path,
                    label: str, strict: bool = True) -> bool:
    """ws2 on the discard twin: the discard identity at every row, and the
    real negative — the identity without film_removed fails."""
    case_dir = stage_case(args.repo, case, workdir)
    config = case_dir / config_name

    def to_discard(doc):
        doc["transport"]["surface_dryout"] = "discard"
        doc["output"]["variables"]["surface"] = [
            v for v in doc["output"]["variables"]["surface"] if v != "water_store"]
        doc["output"]["variables"]["transport"] = [
            v for v in doc["output"]["variables"]["transport"] if v != "salt_store"]
        doc["simulation"]["id"] += "-discard"

    rewrite_config(config, to_discard)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks,
             ws_strict() if strict else None)
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        ma = handle["/monitor/mass_audit"][:]
        ta = handle["/monitor/transport_audit"][:]
    tw = ws_tolerances("ws2")
    cols = gate_ws.audit_columns(coupled=True, store=False)
    good, msg = gate_ws.check_water_identity(ma, cols, False, tw["water_rel"], tw["water_abs_m3"])
    print(f"  {label} discard ws2 {msg}{'' if good else ' FAIL'}")
    ok = good
    bad = ma.copy()
    bad[:, cols["film_removed"]] = 0.0
    neg, msg = gate_ws.check_water_identity(bad, cols, False, tw["water_rel"], tw["water_abs_m3"])
    good = not neg
    print(f"  {label} discard negative (film_removed dropped): "
          f"{'fails as required' if good else 'STILL CLOSES — vacuous FAIL'}")
    ok = ok and good
    good, msg = gate_ws.check_salt_identity(ta, False, tw["salt_rel"])
    print(f"  {label} discard ws2 {msg}{'' if good else ' FAIL'}")
    return ok and good


def gate_ws1(args: argparse.Namespace) -> int:
    tol = ws_tolerances("ws1")
    case_dir = stage_case(args.repo, "ws-water-salt-store", args.work / "store")
    config = case_dir / "ws1-source-bound.yaml"
    if not validate_config(args.frehg, config):
        print(f"ws1 gate: {CAPABILITY_ABSENT_S3WS}")
        print("ws1 gate: FAIL")
        return 1
    ok = ws_source_bound_lane(args, case_dir, config, tol, "ws1")
    ok = ws_discard_twin(args, "ws-water-salt-store", "ws1-source-bound.yaml",
                         args.work / "discard", "ws1") and ok
    print("ws1 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_ws1n(args: argparse.Namespace) -> int:
    tol = ws_tolerances("ws1n")
    case_dir = stage_case(args.repo, "ws-water-salt-store/ws1n", args.work / "store")
    config = case_dir / "ws1n-h2-beach.yaml"
    if not validate_config(args.frehg, config):
        print(f"ws1n gate: {CAPABILITY_ABSENT_S3WS}")
        print("ws1n gate: FAIL")
        return 1
    # The owner's repro runs on the production (default) solver, as packaged.
    ok = ws_source_bound_lane(args, case_dir, config, tol, "ws1n", strict=False)
    print("ws1n gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def ws_run_lane(args, config_name: str, workdir: Path, mutate=None,
                strict: bool = False) -> tuple[Path, dict]:
    case_dir = stage_case(args.repo, "ws-water-salt-store", workdir)
    config = case_dir / config_name
    if mutate is not None:
        rewrite_config(config, mutate)
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks,
             ws_strict() if strict else None)
    with open(case_dir / "out" / "run-record.yaml", encoding="utf-8") as h:
        record = yaml.safe_load(h)
    return case_dir / "out" / "output.h5", record


def gate_ws3(args: argparse.Namespace) -> int:
    tol = ws_tolerances("ws3")
    ok = True
    probe = stage_case(args.repo, "ws-water-salt-store", args.work / "probe")
    if not validate_config(args.frehg, probe / "ws3bc-pond.yaml"):
        print(f"ws3 gate: {CAPABILITY_ABSENT_S3WS}")
        print("ws3 gate: FAIL")
        return 1
    sys.path.insert(0, str(HERE.parent.parent / "tools"))
    import h5_bitwise  # noqa: E402

    # (a) the never-drying rain/evaporation basin: store == discard bitwise
    # on every field (the audit tables: the shared prefix bitwise, the store
    # lane's appended columns aside).
    out_s, _ = ws_run_lane(args, "ws3a-basin.yaml", args.work / "a-store")

    def discard(doc):
        doc["transport"]["surface_dryout"] = "discard"
        doc["simulation"]["id"] += "-discard"

    out_d, _ = ws_run_lane(args, "ws3a-basin.yaml", args.work / "a-discard", discard)
    hs, hd = h5_bitwise.dataset_hashes(out_s), h5_bitwise.dataset_hashes(out_d)
    fields = sorted(k for k in hd if not k.startswith("monitor/"))
    same = all(hs.get(k) == hd[k] for k in fields)
    with h5py.File(out_s, "r") as a, h5py.File(out_d, "r") as b:
        prefix_ok = True
        for tab in ("monitor/mass_audit", "monitor/transport_audit"):
            ta, tb = a[tab][()], b[tab][()]
            prefix_ok = prefix_ok and np.array_equal(ta[:, :tb.shape[1]], tb)
        ts = output_times(a, "/transport/concentration_surface")
        c = [float(np.nanmean(a[f"/transport/concentration_surface/{t}"][:])) for t in ts]
        dmin = min(float(np.nanmin(a[f"/surface/depth/{t}"][:])) for t in ts)
    good = same and prefix_ok and (max(c) - min(c)) > tol["move_floor_psu"] and dmin > 0.0
    print(f"  (a) store vs discard: {len(fields)} field datasets "
          f"{'bitwise' if same else 'DIFFER'}, audit prefixes "
          f"{'bitwise' if prefix_ok else 'DIFFER'}; salinity moves "
          f"{max(c) - min(c):.3f} psu (> {tol['move_floor_psu']}), never dries "
          f"(min depth {dmin:.4f} m){'' if good else ' FAIL'}")
    ok = ok and good

    # (a-ET) the ET open-water basin plus rain: bitwise against its T8 hash.
    out_e, _ = ws_run_lane(args, "ws3aet-basin.yaml", args.work / "a-et")
    ref = HERE / "reference-hashes" / "stage3-ws3aet-T8.json"
    code = h5_bitwise.main_from(["compare", str(out_e), str(ref), "--case", "ws3aet",
                                 "--require-toolchain", "--allow", "monitor/mass_audit",
                                 "--allow", "monitor/transport_audit"])
    good = code in (0, 77)
    print(f"  (a-ET) vs the T8-binary hash: {'skipped off-toolchain' if code == 77 else ('PASS' if code == 0 else 'FAIL')}")
    ok = ok and good

    # (b)+(c) the ET pond: dry-out books, then the rewet closed form.
    out_b, _ = ws_run_lane(args, "ws3bc-pond.yaml", args.work / "bc")
    a_cell = 1.0
    v0, m0 = 0.005 * a_cell, 35.0 * 0.005 * a_cell
    with h5py.File(out_b, "r") as h:
        ts = output_times(h, "/surface/water_store")
        vs = np.array([float(h[f"/surface/water_store/{t}"][:].ravel()[0]) for t in ts])
        msx = np.array([float(h[f"/transport/salt_store/{t}"][:].ravel()[0]) for t in ts])
        dp = np.array([float(h[f"/surface/depth/{t}"][:].ravel()[0]) for t in ts])
        ma = h["/monitor/mass_audit"][:]
        ta = h["/monitor/transport_audit"][:]
        mon = h["/monitor/pond"][:]
    cols = gate_ws.audit_columns(coupled=False, store=True)
    k = max(i for i, t in enumerate(ts) if t < 7200)
    row = int(np.searchsorted(ma[:, 0], ts[k]))
    evap = float(ma[row, cols["evaporation"]] - ma[0, cols["evaporation"]])
    good = (dp[k] == 0.0 and abs(evap + vs[k] - v0) <= tol["exact_rel"] * v0
            and vs[k] <= 1.0e-5 * a_cell and abs(msx[k] - m0) <= tol["exact_rel"] * m0)
    print(f"  (b) dry at {ts[k]} s: Sigma E + V_s - V0 = {evap + vs[k] - v0:.2e}, V_s "
          f"{vs[k]:.3e} (<= min_depth A), M_s - M0 = {msx[k] - m0:.2e}"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    good, msg = gate_ws.check_monotone_while_dry(vs[:, None], (dp == 0.0)[:, None])
    print(f"  (b) {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    # (c) S3-A13: the ET channel keeps evaporating during the rain step, so
    # the rewet closed form carries that step's E dt: c = M0 / (V_s +
    # (h_r - E dt) A) (the plan's form assumed no concurrent evaporation).
    i = next(i for i in range(len(mon)) if mon[i, 0] > 7200.0 and mon[i, 1] > 0.0)
    h_r, e_dt = 1.0e-5 * 10.0, 1.5e-6 * 10.0
    expect = m0 / (vs[k] + (h_r - e_dt) * a_cell)
    tot = ta[:, 1] + ta[:, 11]
    good = (abs(mon[i, 2] - expect) <= tol["exact_rel"] * expect
            and float(np.max(np.abs(tot - tot[0]))) <= tol["exact_rel"] * tot[0])
    print(f"  (c) first wet step t = {mon[i, 0]:.0f} s: c = {mon[i, 2]:.10f} vs "
          f"M0/(V_s + (h_r - E dt) A) = {expect:.10f}; total salt drift "
          f"{float(np.max(np.abs(tot - tot[0]))) / tot[0]:.2e}{'' if good else ' FAIL'}")
    ok = ok and good

    # (d) legacy prescribed evaporation (no ET) to dryness.
    out_d4, _ = ws_run_lane(args, "ws3d-pond.yaml", args.work / "d")
    with h5py.File(out_d4, "r") as h:
        ts = output_times(h, "/surface/water_store")
        vs = np.array([float(h[f"/surface/water_store/{t}"][:].ravel()[0]) for t in ts])
        msx = np.array([float(h[f"/transport/salt_store/{t}"][:].ravel()[0]) for t in ts])
        dp = np.array([float(h[f"/surface/depth/{t}"][:].ravel()[0]) for t in ts])
        ma = h["/monitor/mass_audit"][:]
    evap = float(ma[-1, cols["evaporation"]] - ma[0, cols["evaporation"]])
    clamp = float(ma[-1, cols["clamped"]] - ma[0, cols["clamped"]])
    good = (dp[-1] == 0.0 and abs((evap - clamp) + vs[-1] - v0) <= tol["exact_rel"] * v0
            and vs[-1] <= 1.0e-5 * a_cell and abs(msx[-1] - m0) <= tol["exact_rel"] * m0)
    print(f"  (d) legacy evaporation: Sigma (E - clamped) + V_s - V0 = "
          f"{(evap - clamp) + vs[-1] - v0:.2e}, V_s {vs[-1]:.3e}, M_s - M0 = "
          f"{msx[-1] - m0:.2e}{'' if good else ' FAIL'}")
    ok = ok and good
    good, msg = gate_ws.check_monotone_while_dry(vs[:, None], (dp == 0.0)[:, None])
    print(f"  (d) {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    print("ws3 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_ws4(args: argparse.Namespace) -> int:
    tol = ws_tolerances("ws4")
    ok = True
    probe = stage_case(args.repo, "ws-water-salt-store", args.work / "probe")
    if not validate_config(args.frehg, probe / "ws4-p2-ebb.yaml"):
        print(f"ws4 gate: {CAPABILITY_ABSENT_S3WS}")
        print("ws4 gate: FAIL")
        return 1
    for name, path_key, c0, pair_tol, temp in (
            ("ws4-p2-ebb", "film_removed_small_m3", 35.0, tol["pair_slope_rel"], False),
            ("ws4-p2t-ebb", "film_removed_small_m3", 35.0, tol["pair_slope_rel"], True),
            ("ws4-p3-evap", "film_removed_clamp_m3", None, tol["pair_exact_rel"], False),
            ("ws4-p3t-evap", "film_removed_clamp_m3", None, tol["pair_exact_rel"], True)):
        out, record = ws_run_lane(args, f"{name}.yaml", args.work / name)
        surf = record["closure"]["surface"]
        with h5py.File(out, "r") as h:
            vs, ms = ws_store_frames(h)
            heat = h["/monitor/temperature_audit"][-1, 12] if temp else None
            t0 = h["/transport/concentration_surface/0"][:].ravel()
            ma = h["/monitor/mass_audit"][:]
        # The double-precision books (the run record prints 11 digits):
        # film_removed (every path, cumulative), store_volume (instant),
        # store_release (cumulative).
        cols = gate_ws.audit_columns(coupled=False, store=True)
        removed = float(ma[-1, cols["film_removed"]])
        released = float(ma[-1, cols["store_release"]])
        held = float(ma[-1, cols["store_volume"]])
        own = float(surf[path_key])
        good = abs(removed - held - released) <= tol["books_rel"] * max(removed, 1.0e-30) \
            and own >= tol["path_floor_m3"]
        print(f"  [{name}] own path {path_key} = {own:.3e} m^3 (>= "
              f"{tol['path_floor_m3']:.0e}); removed {removed:.6e} = held {held:.6e} + "
              f"released {released:.3e} (res {removed - held - released:.1e})"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
        mask = vs[-1] > 0.0
        if c0 is not None:
            good, msg = gate_ws.check_paired(ms[-1], vs[-1], c0, pair_tol)
        else:
            # P3 under evaporation: the salt crust (R4a) — the store holds the
            # whole M0; the paired part is c V_s at the concentrated value.
            m0 = float(np.nansum(t0)) * 0.005
            err = abs(float(ms[-1].sum()) - m0) / m0
            good = err <= pair_tol
            msg = f"crust + paired M_s = M0 (rel {err:.2e} <= {pair_tol:.0e})"
        print(f"  [{name}] {msg}{'' if good else ' FAIL'}")
        ok = ok and good
        if temp:
            # R9: heat parks with its water at 20 C (P3: exact; P2: within
            # the slope lane's uniformity class), no crust (R4 is salinity-only).
            ratio = float(heat) / (20.0 * held)
            bound = pair_tol if c0 is None else 1.0e-2
            good = abs(ratio - 1.0) <= bound
            print(f"  [{name}] R9 parked heat / (20 V_s) - 1 = {ratio - 1.0:.2e} "
                  f"(<= {bound:.0e}){'' if good else ' FAIL'}")
            ok = ok and good
    print("ws4 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_ws_restart(args: argparse.Namespace) -> int:
    """ws-x restart: ws1 checkpointed while >= 10 cells hold parked water;
    the continuation is bitwise on every field including water_store."""
    base = "ws1-source-bound.yaml"
    full_dir = stage_case(args.repo, "ws-water-salt-store", args.work / "full")
    config = full_dir / base
    if not validate_config(args.frehg, config):
        print(f"ws restart gate: {CAPABILITY_ABSENT_S3WS}")
        print("ws restart determinism: FAIL")
        return 1
    run_case(args.frehg, config, full_dir, args.mpiexec, args.ranks, ws_strict())
    full_out = full_dir / "out" / "output.h5"
    with h5py.File(full_out, "r") as h:
        cands = [t for t in output_times(h, "/surface/water_store")
                 if t % 1800 == 0 and 0 < t < 7200
                 and int(np.sum(h[f"/surface/water_store/{t}"][:] > 0.0)) >= 10]
    if not cands:
        print("  no checkpoint with >= 10 parked cells (vacuous) FAIL")
        print("ws restart determinism: FAIL")
        return 1
    ckpt = cands[0]
    re_dir = stage_case(args.repo, "ws-water-salt-store", args.work / "restart")
    re_cfg = re_dir / base

    def add_restart(doc):
        doc["restart"] = {"enabled": True, "file": str(full_out.resolve()), "time": ckpt}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(re_cfg, add_restart)
    run_case(args.frehg, re_cfg, re_dir, args.mpiexec, args.ranks, ws_strict())
    ok = True
    with h5py.File(full_out, "r") as full, \
            h5py.File(re_dir / "out" / "restarted.h5", "r") as part:
        nparked = int(np.sum(full[f"/surface/water_store/{ckpt}"][:] > 0.0))
        print(f"  checkpoint t = {ckpt} s with {nparked} parked cells")
        for group in ("/surface/eta", "/surface/depth", "/surface/water_store",
                      "/transport/concentration_surface", "/transport/salt_store",
                      "/transport/concentration"):
            worst_same = True
            for t in output_times(part, group):
                if t <= ckpt:
                    continue
                same = np.array_equal(np.nan_to_num(full[f"{group}/{t}"][:], nan=-9e99),
                                      np.nan_to_num(part[f"{group}/{t}"][:], nan=-9e99))
                worst_same = worst_same and same
            print(f"  restart {group.rsplit('/', 1)[-1]}: "
                  f"{'bitwise' if worst_same else 'DIFFERS'}")
            ok = ok and worst_same
    print("ws restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_ws_rank(args: argparse.Namespace) -> int:
    """ws-x rank (S3-A13). ws1 is trajectory-chaotic in BOTH modes (a 1e-12
    bed twin moves eta 6.8e-4 and flips a cell's surface salinity by 35 psu;
    the discard twin is as rank-sensitive), so the field-identity classes
    gate on the trajectory-stable paired-store base — the ws4 P2 slope —
    and ws1 gates the rank-INVARIANT guarantees at n = 2/4: the source
    bound, I2 and both ws2 ledgers. Per cell the parked volume inherits
    eta's rank class (V_s = A h_parked): A * rank_surface_abs (water) and
    c_max A rank_surface_abs (salt). The store TOTALS gate at
    store_sum_rel = 1e-9 relative: the plan's 1e-12 is below the strict
    solve's rounding at the slope's ~1e-5 m^3 scale (measured 5.4e-11)."""
    tol = pd_tolerances("pd3")
    tws = ws_tolerances("ws1")
    ok = True
    # (A) identity classes on the stable slope base, plus its chaos control.
    outs = {}
    for tag, n, twin in (("n1", 1, False), ("twin", 1, True), ("n2", 2, False),
                         ("n4", 4, False)):
        case_dir = stage_case(args.repo, "ws-water-salt-store", args.work / f"slope-{tag}")
        config = case_dir / "ws4-p2-ebb.yaml"

        def setup(doc, n=n):
            doc["time"]["output_interval"] = 60.0
            if n > 1:
                doc["domain"]["decomposition"] = {"mpi_nx": 1, "mpi_ny": n}

        rewrite_config(config, setup)
        if not validate_config(args.frehg, config):
            print(f"ws rank gate: {CAPABILITY_ABSENT_S3WS}")
            print("ws rank gate: FAIL")
            return 1
        if twin:
            bed = case_dir / "input" / "ws4_slope.dat"
            lines = bed.read_text().splitlines()
            k = next(i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")) + 4
            lines[k] = repr(float(lines[k]) + 1.0e-12)
            bed.write_text("\n".join(lines) + "\n")
        run_case(args.frehg, config, case_dir, args.mpiexec, n, STRICT_PETSC_OPTIONS)
        outs[tag] = case_dir / "out" / "output.h5"
    area = 1.0
    with h5py.File(outs["n1"], "r") as base:
        times = output_times(base, "/surface/eta")
        for tag in ("twin", "n2", "n4"):
            with h5py.File(outs[tag], "r") as other:
                for group, bound in (
                        ("/surface/eta", 1.0e-10 if tag == "twin" else tol["rank_surface_abs"]),
                        ("/transport/concentration_surface",
                         1.0e-10 if tag == "twin" else tol["rank_tracer_abs"]),
                        ("/surface/water_store",
                         area * (1.0e-10 if tag == "twin" else tol["rank_surface_abs"])),
                        ("/transport/salt_store",
                         tws["c_max"] * area * (1.0e-10 if tag == "twin" else tol["rank_surface_abs"]))):
                    worst = max(float(np.nanmax(np.abs(base[f"{group}/{t}"][:] -
                                                       other[f"{group}/{t}"][:])))
                                for t in times)
                    good = worst <= bound
                    what = "1e-12 bed twin (stability control)" if tag == "twin" else f"n={tag[1:]}"
                    print(f"  slope n=1 vs {what}: {group.rsplit('/', 1)[-1]} max |d| = "
                          f"{worst:.3e} (<= {bound:.0e}){'' if good else ' FAIL'}")
                    ok = ok and good
                for group in ("/surface/water_store", "/transport/salt_store"):
                    worst = 0.0
                    for t in times:
                        a = float(np.nansum(base[f"{group}/{t}"][:]))
                        b = float(np.nansum(other[f"{group}/{t}"][:]))
                        if a > 0.0:
                            worst = max(worst, abs(a - b) / a)
                    good = worst <= tws["store_sum_rel"]
                    print(f"  slope n=1 vs {tag}: Sigma {group.rsplit('/', 1)[-1]} rel "
                          f"{worst:.3e} (<= {tws['store_sum_rel']:.0e}){'' if good else ' FAIL'}")
                    ok = ok and good
                flips = 0
                for t in times:
                    a = base[f"/surface/water_store/{t}"][:] > 0.0
                    b = other[f"/surface/water_store/{t}"][:] > 0.0
                    flips += int(np.sum(a != b))
                good = flips == 0
                print(f"  slope n=1 vs {tag}: parking-timing flips {flips} (== 0)"
                      f"{'' if good else ' FAIL'}")
                ok = ok and good
    # (B) the rank-invariant guarantees on the chaotic ws1 at n = 2 and 4.
    for n in (2, 4):
        case_dir = stage_case(args.repo, "ws-water-salt-store", args.work / f"ws1-n{n}")
        config = case_dir / "ws1-source-bound.yaml"

        def decomp(doc, n=n):
            doc["domain"]["decomposition"] = {"mpi_nx": 1, "mpi_ny": n}

        rewrite_config(config, decomp)
        saved = args.ranks
        args.ranks = n
        try:
            good = ws_source_bound_lane(args, case_dir, config, tws, f"ws1 n={n}")
        finally:
            args.ranks = saved
        ok = ok and good
    print("ws rank gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Stage-3 T10 solve-consistent surface transport gates (plan §12B, O16/
# S3-A14): sx1 (uniform scalar through wetting/drying + the anchor), sx2
# (drying-face salt on the T5 et8d strip), sx3 (the ws1n anchor, nightly),
# sx4 (a tracer front), sx5 (temperature), sx-restart/sx-rank/sx-seam.
# Pure checks: gate_sx.py.
# ---------------------------------------------------------------------------
import gate_sx  # noqa: E402

CAPABILITY_ABSENT_S3SX = ("capability absent: the schema rejects the T10 "
                          "target token (transport.surface_update: "
                          "conservative)")


def sx_tolerances(gate: str) -> dict:
    with open(HERE / "tolerances" / "stage3-sx.yaml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)[gate]


def sx_frames(handle, group: str):
    """Concentration frames with their wet masks (depth > 0)."""
    times = output_times(handle, group)
    frames, wet = [], []
    for t in times:
        frames.append(np.nan_to_num(handle[f"{group}/{t}"][:].ravel()))
        wet.append(handle[f"/surface/depth/{t}"][:].ravel() > 0.0)
    return times, frames, wet


def gate_sx1(args: argparse.Namespace, case: str = "sx1-uniform-beach.yaml",
             label: str = "sx1") -> int:
    tol = sx_tolerances("sx1")
    case_dir = stage_case(args.repo, "sx-solve-consistent", args.work)
    config = case_dir / case
    if not validate_config(args.frehg, config):
        print(f"{label} gate: {CAPABILITY_ABSENT_S3SX}")
        print(f"{label} gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks, ws_strict())
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        times, frames, wet = sx_frames(handle, "/transport/concentration_surface")
        good, msg = gate_sx.check_uniform(frames, wet, tol["c0"],
                                          tol["uniform_abs"], "surface")
        print(f"  {msg}{'' if good else ' FAIL'}")
        ok = ok and good
        subs, swet = [], []
        for t in output_times(handle, "/transport/concentration"):
            c = handle[f"/transport/concentration/{t}"][:].ravel()
            m = np.isfinite(c)
            subs.append(np.nan_to_num(c))
            swet.append(m)
        good, msg = gate_sx.check_uniform(subs, swet, tol["c0"],
                                          tol["uniform_abs"], "subsurface")
        print(f"  {msg}{'' if good else ' FAIL'}")
        ok = ok and good
        ta = handle["/monitor/transport_audit"][:]
        m_scale = float(ta[0, 1] + ta[0, 2] + (ta[0, 11] if ta.shape[1] > 11 else 0.0))
        good, msg = gate_sx.check_anchor(float(ta[-1, 7]), m_scale,
                                         tol["anchor_rel"], "J2 anchor")
        print(f"  {msg}{'' if good else ' FAIL'}")
        ok = ok and good
        adj = abs(float(ta[-1, 6]))
        good = adj <= tol["adjust_rel"] * m_scale
        print(f"  C4 booked surf_adjust = {adj:.3e} (<= {tol['adjust_rel']:.0e} M)"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
        # Anti-vacuity: the beach wets and dries (parking + release events).
        ma = handle["/monitor/mass_audit"][:]
        cols = gate_ws.audit_columns(coupled=True, store=True)
        film = float(ma[-1, cols["film_removed"]])
        good = film > 0.0
        print(f"  anti-vacuity: film_removed = {film:.3e} m^3 (> 0)"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
        # Temperature instance (sx5 rides this body).
        if "/monitor/temperature_audit" in handle:
            tol5 = sx_tolerances("sx5")
            tt, tf, twet = sx_frames(handle, "/temperature/temperature_surface")
            good, msg = gate_sx.check_uniform(tf, twet, tol5["t0"],
                                              tol5["uniform_abs"], "temperature")
            print(f"  {msg}{'' if good else ' FAIL'}")
            ok = ok and good
            taT = handle["/monitor/temperature_audit"][:]
            h_scale = float(taT[0, 1] + taT[0, 2] + (taT[0, 12] if taT.shape[1] > 12 else 0.0))
            good, msg = gate_sx.check_anchor(float(taT[-1, 7]), h_scale,
                                             tol5["anchor_rel"], "heat anchor")
            print(f"  {msg}{'' if good else ' FAIL'}")
            ok = ok and good
    print(f"{label} gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_sx2(args: argparse.Namespace) -> int:
    tol = sx_tolerances("sx2")
    case_dir = stage_case(args.repo, "sx-solve-consistent", args.work)
    config = case_dir / "sx2-drying-strip.yaml"
    if not validate_config(args.frehg, config):
        print(f"sx2 gate: {CAPABILITY_ABSENT_S3SX}")
        print("sx2 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        ta = handle["/monitor/transport_audit"][:]
        # The ledger identity, not raw constancy: the bottom head BC feeds
        # 10-psu groundwater in (subs_boundary books it), so total salt
        # legitimately grows. total − Σ(booked cumulative terms) is the
        # unbooked residual J2/J4 guard.
        total = ta[:, 1] + ta[:, 2] + ta[:, 11]
        booked = (ta[:, 4] + ta[:, 5] + ta[:, 6] + ta[:, 7] + ta[:, 8] +
                  ta[:, 9] + ta[:, 10])
        resid = (total - total[0]) - (booked - booked[0])
        m_scale = float(total[0])
        worst = float(np.abs(resid).max())
        good = worst <= tol["total_rel"] * m_scale
        print(f"  ledger identity: worst unbooked residual {worst:.3e} of "
              f"{m_scale:.4e} (<= {tol['total_rel']:.0e}){'' if good else ' FAIL'}")
        ok = ok and good
        adj = abs(float(ta[-1, 6]))
        good = adj <= tol["adjust_rel"] * m_scale
        print(f"  J4 stranded salt: surf_adjust = {adj:.3e} (<= "
              f"{tol['adjust_rel']:.0e} M = {tol['adjust_rel'] * m_scale:.1e}; "
              f"red on T9: 1.09e-1){'' if good else ' FAIL'}")
        ok = ok and good
        # Anti-vacuity: the salt actually crossed the drying face — the
        # pit column banked real salt (its water parks below min_depth, so
        # the store, not the concentration field, carries the evidence).
        times = output_times(handle, "/transport/salt_store")
        peak = max(float(np.nan_to_num(handle[f"/transport/salt_store/{t}"][:]).ravel()[1])
                   for t in times)
        good = peak > tol["neighbor_floor_salt"]
        print(f"  anti-vacuity: pit-cell peak salt store {peak:.4f} psu m^3 "
              f"(> {tol['neighbor_floor_salt']}){'' if good else ' FAIL'}")
        ok = ok and good
    print("sx2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_sx3(args: argparse.Namespace) -> int:
    """ws1n (the owner's repro) under `conservative`, production solver."""
    tol = sx_tolerances("sx3")
    case_dir = stage_case(args.repo, "ws-water-salt-store/ws1n", args.work)
    config = case_dir / "ws1n-h2-beach.yaml"

    def to_conservative(doc):
        doc["transport"]["surface_update"] = "conservative"
        doc["simulation"]["id"] += "-conservative"

    rewrite_config(config, to_conservative)
    if not validate_config(args.frehg, config):
        print(f"sx3 gate: {CAPABILITY_ABSENT_S3SX}")
        print("sx3 gate: FAIL")
        return 1
    tws = ws_tolerances("ws1n")
    ok = ws_source_bound_lane(args, case_dir, config, tws, "sx3", strict=False)
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        ta = handle["/monitor/transport_audit"][:]
        m_scale = float(ta[0, 1] + ta[0, 2] + ta[0, 11])
        good, msg = gate_sx.check_anchor(float(ta[-1, 7]), m_scale,
                                         tol["anchor_rel"],
                                         "J2 anchor (red on T9: +28.4 = 2.4e-2 M)")
        print(f"  {msg}{'' if good else ' FAIL'}")
        ok = ok and good
    print("sx3 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_sx4(args: argparse.Namespace) -> int:
    tol = sx_tolerances("sx4")
    case_dir = stage_case(args.repo, "sx-solve-consistent", args.work)
    config = case_dir / "sx4-front-slope.yaml"
    if not validate_config(args.frehg, config):
        print(f"sx4 gate: {CAPABILITY_ABSENT_S3SX}")
        print("sx4 gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks,
             STRICT_PETSC_OPTIONS)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        times, frames, wet = sx_frames(handle, "/transport/concentration_surface")
        ta = handle["/monitor/transport_audit"][:]
        # Ledger identity (boundary/adjust/anchor are booked columns — the
        # slope's walls keep boundary at zero, but the overdraft pairing's
        # negative-residual path books into surf_adjust by design).
        total = ta[:, 1] + ta[:, 11]
        booked = ta[:, 4] + ta[:, 5] + ta[:, 6] + ta[:, 7]
        resid = (total - total[0]) - (booked - booked[0])
        m_scale = float(total[0])
        worst = float(np.abs(resid).max())
        good = worst <= tol["total_rel"] * m_scale
        print(f"  ledger identity: worst unbooked residual {worst:.3e} of "
              f"{m_scale:.4e} (<= {tol['total_rel']:.0e}){'' if good else ' FAIL'}")
        ok = ok and good
        adj = abs(float(ta[-1, 6]))
        good = adj <= tol["adjust_rel"] * m_scale
        print(f"  C4 booked surf_adjust = {adj:.3e} (<= "
              f"{tol['adjust_rel']:.0e} M){'' if good else ' FAIL'}")
        ok = ok and good
        good, msg = gate_sx.check_extrema(frames, wet, 5.0, 35.0,
                                          tol["extrema_tol"])
        print(f"  {msg}{'' if good else ' FAIL'}")
        ok = ok and good
        # Anti-vacuity: the front moves — a below-front cell rises.
        rise = max(float(f[6]) for f, w in zip(frames, wet) if w[6]) \
            if any(w[6] for w in wet) else 0.0
        good = rise > 5.0 + tol["front_floor_psu"]
        print(f"  anti-vacuity: cell 6 (initially 5 psu) peaks at {rise:.3f} "
              f"(> {5.0 + tol['front_floor_psu']}){'' if good else ' FAIL'}")
        ok = ok and good
    print("sx4 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_sx_restart(args: argparse.Namespace) -> int:
    """sx1 checkpointed mid-run; the continuation bitwise on every field."""
    base = "sx1-uniform-beach.yaml"
    full_dir = stage_case(args.repo, "sx-solve-consistent", args.work / "full")
    config = full_dir / base
    if not validate_config(args.frehg, config):
        print(f"sx restart gate: {CAPABILITY_ABSENT_S3SX}")
        print("sx restart determinism: FAIL")
        return 1
    run_case(args.frehg, config, full_dir, args.mpiexec, args.ranks, ws_strict())
    full_out = full_dir / "out" / "output.h5"
    ckpt = 1800
    re_dir = stage_case(args.repo, "sx-solve-consistent", args.work / "restart")
    re_cfg = re_dir / base

    def add_restart(doc):
        doc["restart"] = {"enabled": True, "file": str(full_out.resolve()),
                          "time": ckpt}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(re_cfg, add_restart)
    run_case(args.frehg, re_cfg, re_dir, args.mpiexec, args.ranks, ws_strict())
    ok = True
    with h5py.File(full_out, "r") as full, \
            h5py.File(re_dir / "out" / "restarted.h5", "r") as part:
        for group in ("/surface/eta", "/surface/depth", "/surface/water_store",
                      "/transport/concentration_surface", "/transport/salt_store",
                      "/transport/concentration"):
            same = all(np.array_equal(
                np.nan_to_num(full[f"{group}/{t}"][:], nan=-9e99),
                np.nan_to_num(part[f"{group}/{t}"][:], nan=-9e99))
                for t in output_times(part, group) if t > ckpt)
            print(f"  restart {group.rsplit('/', 1)[-1]}: "
                  f"{'bitwise' if same else 'DIFFERS'}")
            ok = ok and same
    print("sx restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_sx_rank(args: argparse.Namespace) -> int:
    """sx4 (trajectory-stable slope, S3-A13 doctrine) at mpi_ny = 1/2/4
    strict, plus a 1e-12 bed twin as the stability control."""
    tol = pd_tolerances("pd3")
    outs = {}
    for tag, n, twin in (("n1", 1, False), ("twin", 1, True), ("n2", 2, False),
                         ("n4", 4, False)):
        case_dir = stage_case(args.repo, "sx-solve-consistent",
                              args.work / f"slope-{tag}")
        config = case_dir / "sx4-front-slope.yaml"

        def setup(doc, n=n):
            if n > 1:
                doc["domain"]["decomposition"] = {"mpi_nx": 1, "mpi_ny": n}

        rewrite_config(config, setup)
        if not validate_config(args.frehg, config):
            print(f"sx rank gate: {CAPABILITY_ABSENT_S3SX}")
            print("sx rank gate: FAIL")
            return 1
        if twin:
            bed = case_dir / "input" / "ws4_slope.dat"
            lines = bed.read_text().splitlines()
            k = next(i for i, l in enumerate(lines)
                     if l.strip() and not l.startswith("#")) + 4
            lines[k] = repr(float(lines[k]) + 1.0e-12)
            bed.write_text("\n".join(lines) + "\n")
        run_case(args.frehg, config, case_dir, args.mpiexec, n, STRICT_PETSC_OPTIONS)
        outs[tag] = case_dir / "out" / "output.h5"
    ok = True
    with h5py.File(outs["n1"], "r") as base:
        times = output_times(base, "/surface/eta")
        for tag in ("twin", "n2", "n4"):
            with h5py.File(outs[tag], "r") as other:
                # Twin bounds: eta/water respond linearly to the 1e-12 bed
                # perturbation; the scalar's superbee weights amplify it
                # ~1e3 through the drying front (measured 1.28e-9 at first
                # green, S3-A15) — still trajectory-stable (a chaotic lane
                # is O(1)), so the scalar control binds at 1e-7.
                for group, bound in (
                        ("/surface/eta",
                         1.0e-10 if tag == "twin" else tol["rank_surface_abs"]),
                        ("/transport/concentration_surface",
                         1.0e-7 if tag == "twin" else tol["rank_tracer_abs"]),
                        ("/surface/water_store",
                         1.0e-10 if tag == "twin" else tol["rank_surface_abs"])):
                    worst = max(float(np.nanmax(np.abs(
                        base[f"{group}/{t}"][:] - other[f"{group}/{t}"][:])))
                        for t in times)
                    good = worst <= bound
                    what = ("1e-12 bed twin (stability control)"
                            if tag == "twin" else f"n={tag[1:]}")
                    print(f"  slope n=1 vs {what}: {group.rsplit('/', 1)[-1]} "
                          f"max |d| = {worst:.3e} (<= {bound:.0e})"
                          f"{'' if good else ' FAIL'}")
                    ok = ok and good
    print("sx rank gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et7t(args: argparse.Namespace) -> int:
    ok = True
    # (a) a column entirely below h4: T_act == 0, unrealised = demand,
    # the soil channel still works.
    a_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "a")
    a_cfg = a_dir / "et4-balance.yaml"

    def lane_a(doc):
        doc["initial_conditions"]["groundwater"] = {"head": {"constant": -5.0}}
        doc["vegetation"] = {
            "types": [{"name": "stressed", "lai": 2.0, "extinction": 0.5,
                       "leaf_resistance": 100.0, "root_depth": 0.4,
                       "feddes": {"h3": -1.0, "h4": -2.0}}],
            "map": {"constant": "stressed"}}
        doc["time"]["t_end"] = 7200.0
        doc["time"]["output_interval"] = 3600.0
        doc["output"]["variables"].pop("atmosphere", None)
        doc["simulation"]["id"] += "-stress"

    rewrite_config(a_cfg, lane_a)
    if not validate_config(args.frehg, a_cfg):
        print(f"et7-t gate: {CAPABILITY_ABSENT_S3VG}")
        print("et7-t gate: FAIL")
        return 1
    run_case(args.frehg, a_cfg, a_dir, args.mpiexec, args.ranks)
    with h5py.File(a_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
    t_act = float(audit[-1][9])
    t_un = float(audit[-1][10])
    e_soil = abs(float(audit[-1][7]))
    good = t_act == 0.0 and t_un > 0.0 and e_soil > 0.0
    print(f"  [a] below h4: T_act = {t_act:.3e} (== 0), unrealised = "
          f"{t_un:.3e} (> 0), E_soil = {e_soil:.3e} (> 0)"
          f"{'' if good else ' FAIL'}")
    ok = ok and good

    # (c) root zone at theta_r + 1e-4 under a huge demand: the floor holds
    # and the books balance exactly.
    c_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "c")
    c_cfg = c_dir / "et7-cap.yaml"

    def lane_c(doc):
        doc["evapotranspiration"]["open_water"] = False
        doc["evapotranspiration"].pop("soil", None)
        doc["vegetation"] = {
            "types": [{"name": "thirsty", "lai": 1.0e9, "extinction": 1.0,
                       "leaf_resistance": 0.0, "root_depth": 0.5,
                       "feddes": {"h3": -1000.0, "h4": -2000.0}}],
            "map": {"constant": "thirsty"}}
        doc["output"]["variables"].pop("atmosphere", None)
        doc["simulation"]["id"] += "-tcap"

    rewrite_config(c_cfg, lane_c)
    run_case(args.frehg, c_cfg, c_dir, args.mpiexec, args.ranks)
    with h5py.File(c_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/gw_mass_audit"][:]
        theta_min = min(float(np.nanmin(handle[f"/groundwater/water_content/{t}"][:]))
                        for t in output_times(handle, "/groundwater/water_content"))
    # gw_mass_audit here: ..., vloss(6), T(7), T_un(8) — no soil channel.
    t_act = float(audit[-1][7])
    t_un = float(audit[-1][8])
    demand = 1.0e-5 * 3600.0  # T_pot = E_pot at lai -> inf, A = 1
    good = theta_min >= 0.05 - 1.0e-14
    print(f"  [c] theta floor: min - theta_r = {theta_min - 0.05:.3e} "
          f"(>= -1e-14){'' if good else ' FAIL'}")
    ok = ok and good
    resid = abs((t_act + t_un) - demand) / demand
    good = resid <= 1.0e-12
    print(f"  [c] T + unrealised = demand: rel {resid:.3e} (<= 1e-12)"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    good = abs(float(audit[-1][6])) == 0.0
    print(f"  [c] vloss == 0: {audit[-1][6]:.3e}{'' if good else ' FAIL'}")
    ok = ok and good
    # (b) per-substep T_act <= T_pot holds by the per-cell min() —
    # asserted globally: cumulative T <= cumulative demand everywhere.
    mono = np.all(audit[:, 7] <= (audit[:, 7] + audit[:, 8]) + 1.0e-15)
    print(f"  [b] T_act <= T_pot per row: {'ok' if mono else 'FAIL'}")
    ok = ok and bool(mono)
    print("et7-t gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et8t(args: argparse.Namespace) -> int:
    tol = et_tolerances("et8")
    case_dir = stage_case(args.repo, "et-evapotranspiration", args.work)
    config = case_dir / "et8a-soil-salt.yaml"

    def add_veg(doc):
        doc["evapotranspiration"].pop("soil", None)  # T only
        doc["vegetation"] = {
            "types": [{"name": "halophyte", "lai": 2.0, "extinction": 0.5,
                       "leaf_resistance": 100.0, "root_depth": 0.3,
                       "root_profile": "linear",
                       "feddes": {"h3": -1000.0, "h4": -2000.0}}],
            "map": {"constant": "halophyte"}}
        doc["output"]["variables"]["groundwater"] = [
            "hydraulic_head", "water_content", "root_uptake"]
        doc["output"]["variables"].pop("atmosphere", None)
        doc["simulation"]["id"] += "-veg"

    rewrite_config(config, add_veg)
    if not validate_config(args.frehg, config):
        print(f"et8-t gate: {CAPABILITY_ABSENT_S3VG}")
        print("et8-t gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        audit = handle["/monitor/transport_audit"][:]
        times = output_times(handle, "/transport/concentration")
        c0 = handle[f"/transport/concentration/{times[0]}"][:]
        th0 = handle[f"/groundwater/water_content/{times[0]}"][:]
        c1 = handle[f"/transport/concentration/{times[-1]}"][:]
        th1 = handle[f"/groundwater/water_content/{times[-1]}"][:]
        uptake = handle[f"/groundwater/root_uptake/{times[-1]}"][:]
    m0 = float(audit[0][2])
    good, msg = gate_et.check_conserved(float(audit[-1][2]), m0,
                                        tol["conserve_rel"], "total salt")
    print(f"  {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    good, msg = gate_et.check_cell_salt(c1, th1, c0, th0, tol["cell_rel"])
    print(f"  {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    for name, col in (("subs_adjust", 9), ("subs_anchor", 10)):
        v = abs(float(audit[-1][col]))
        good = v <= tol["adjust_rel"] * m0
        print(f"  {name} = {v:.3e} (<= {tol['adjust_rel']:.0e} M)"
              f"{'' if good else ' FAIL'}")
        ok = ok and good
    # Salinity rise proportional to the uptake (correlation = 1): the
    # concentration factor c1/c0 - 1 ~ uptake_k dz/theta over root cells.
    m = np.isfinite(c1) & np.isfinite(c0) & (np.nan_to_num(uptake) > 0)
    if m.sum() >= 2:
        rise = (c1[m] / c0[m]) - 1.0
        frac = np.nan_to_num(uptake)[m]
        corr = float(np.corrcoef(rise, frac)[0, 1])
        # S3-A9: the multiplicative concentration factor is exp-nonlinear
        # in the uptake (c1/c0 = theta0/theta1), so the linear correlation
        # carries a quadratic linearization deviation (measured 6.9e-6
        # over the 2 h horizon); the EXACT statement is the per-cell
        # c*theta identity above. Bound 1e-4: 14x the measured class.
        good = corr >= 1.0 - 1.0e-4
        print(f"  salinity rise ~ rootFrac: corr = {corr:.10f} "
              f"(>= 1 - 1e-4){'' if good else ' FAIL'}")
        ok = ok and good
    print("et8-t gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et9veg(args: argparse.Namespace) -> int:
    case_dir = stage_case(args.repo, "et-evapotranspiration", args.work)
    config = case_dir / "et9-tidal-flat.yaml"

    def add_marsh(doc):
        doc["vegetation"] = dict(ET9_MARSH)
        doc["output"]["variables"]["atmosphere"] = [
            "evaporation_rate", "soil_evaporation_rate", "wet_time",
            "soil_evaporation_time", "transpiration_rate"]
        doc["simulation"]["id"] += "-veg"

    rewrite_config(config, add_marsh)
    if not validate_config(args.frehg, config):
        print(f"et9-veg gate: {CAPABILITY_ABSENT_S3VG}")
        print("et9-veg gate: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    ok = True
    bed = np.loadtxt(args.repo / "benchmarks" / "et-evapotranspiration" /
                     "input" / "et9_dem.dat").reshape(24, 4)
    veg = np.loadtxt(args.repo / "benchmarks" / "et-evapotranspiration" /
                     "input" / "et9_veg.dat").reshape(24, 4) > 0
    with h5py.File(case_dir / "out" / "output.h5", "r") as handle:
        times = output_times(handle, "/atmosphere/transpiration_rate")
        t_end = float(times[-1])
        wet = handle[f"/atmosphere/wet_time/{times[-1]}"][:].reshape(24, 4)
        soil = handle[f"/atmosphere/soil_evaporation_time/{times[-1]}"][:].reshape(24, 4)
        # Cumulative T per cell from the interval-mean frames.
        t_cum = np.zeros((24, 4))
        prev = 0.0
        for t in times:
            frame = handle[f"/atmosphere/transpiration_rate/{t}"][:].reshape(24, 4)
            t_cum += np.nan_to_num(frame) * (float(t) - prev)
            prev = float(t)
        # T during inundated instants (O7): the marsh band's lower rows
        # flood each cycle — find a frame where a marsh cell is wet AND
        # the T rate that interval is positive.
        o7 = False
        prev = 0.0
        for t in times:
            frame = handle[f"/atmosphere/transpiration_rate/{t}"][:].reshape(24, 4)
            depth = handle[f"/surface/depth/{t}"][:].reshape(24, 4)
            m = veg & (np.nan_to_num(depth) > 0.01)
            if m.any() and float(np.nan_to_num(frame)[m].max()) > 0.0:
                o7 = True
                break
        audit = handle["/monitor/gw_mass_audit"][:]
    good, msg = gate_et.check_no_double_count(wet, soil, t_end)
    print(f"  {msg}{'' if good else ' FAIL'}")
    ok = ok and good
    good = float(t_cum[veg].min()) > 0.0
    print(f"  T > 0 on every marsh cell: min cum = {t_cum[veg].min():.3e}"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    good = float(np.abs(t_cum[~veg]).max()) == 0.0
    print(f"  T == 0 outside the band: max = {np.abs(t_cum[~veg]).max():.3e}"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    print(f"  O7 (transpiration while inundated): "
          f"{'observed' if o7 else 'NOT OBSERVED'}{'' if o7 else ' FAIL'}")
    ok = ok and o7
    # gw_mass_audit here: ..., soil(7), soil_un(8), T(9), T_un(10).
    good = float(audit[-1][9]) > 1.0e-6
    print(f"  anti-vacuity: cumulative T = {audit[-1][9]:.3e} m^3 (> 1e-6)"
          f"{'' if good else ' FAIL'}")
    ok = ok and good
    print("et9-veg gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_et_veg_restart(args: argparse.Namespace) -> int:
    """T6 restart lane: et9-vegetated checkpointed mid-cycle; every field
    including the sink window accumulators restarts bitwise."""
    case_dir = stage_case(args.repo, "et-evapotranspiration", args.work)
    config = case_dir / "et9-tidal-flat.yaml"

    def full_setup(doc):
        doc["vegetation"] = dict(ET9_MARSH)
        doc["output"]["variables"]["atmosphere"] = [
            "wet_time", "soil_evaporation_time", "transpiration_rate"]
        doc["output"]["variables"]["groundwater"] = [
            "hydraulic_head", "water_content", "root_uptake"]
        doc["output"]["checkpoint"] = {"interval": 5400}
        doc["simulation"]["id"] += "-veg"

    rewrite_config(config, full_setup)
    if not validate_config(args.frehg, config):
        print(f"et-veg restart gate: {CAPABILITY_ABSENT_S3VG}")
        print("et-veg restart determinism: FAIL")
        return 1
    run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)
    full_output = case_dir / "out" / "output.h5"
    re_dir = stage_case(args.repo, "et-evapotranspiration", args.work / "restart")
    re_cfg = re_dir / "et9-tidal-flat.yaml"

    def add_restart(doc):
        full_setup(doc)
        doc["restart"] = {"enabled": True, "file": str(full_output.resolve()),
                          "time": 5400}
        doc["output"]["filename"] = "out/restarted.h5"

    rewrite_config(re_cfg, add_restart)
    run_case(args.frehg, re_cfg, re_dir, args.mpiexec, args.ranks)
    ok = True
    with h5py.File(full_output, "r") as full, \
         h5py.File(re_dir / "out" / "restarted.h5", "r") as part:
        for group in ("/surface/eta", "/groundwater/water_content",
                      "/groundwater/root_uptake", "/transport/salt_store",
                      "/atmosphere/wet_time", "/atmosphere/transpiration_rate"):
            for t in output_times(part, group):
                if t <= 5400:
                    continue
                a = full[f"{group}/{t}"][:]
                b = part[f"{group}/{t}"][:]
                same = np.array_equal(np.nan_to_num(a, nan=-9.0e99),
                                      np.nan_to_num(b, nan=-9.0e99))
                if not same or t == output_times(part, group)[-1]:
                    print(f"  restart {group.rsplit('/', 1)[-1]}@{t}: "
                          f"{'bitwise' if same else 'DIFFERS'}")
                ok = ok and same
    print("et-veg restart determinism:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gate", choices=["b1", "b2", "b3", "b4", "b5", "b6", "b1-restart",
                                         "b2-restart", "b5-restart", "b6-restart",
                                         "rank-invariance", "rank-invariance-b2",
                                         "rank-invariance-b5", "b6-gw-smoke",
                                         "smoke-b5", "smoke-b6", "smoke-g5",
                                         "smoke-g8", "g4", "g5",
                                         "g9", "g10", "wind-orient", "g6", "g7",
                                         "g8", "heat-orient", "heat-restart",
                                         "rank-invariance-heat",
                                         "outflow-staircase", "swe-orient",
                                         "gw-orient", "transport-orient",
                                         "surface-bc-sides", "wind-restart",
                                         "evap-restart", "b6-subcycled",
                                         "rank-invariance-tracer",
                                         "rank-invariance-coupled-salt",
                                         "wind-drydown", "terrain-heat",
                                         "two-scalar-restart", "evap-heat",
                                         "heat-subcycled", "default-bitwise",
                                         "rl2", "rl3", "rl4", "rl5-orient", "rl-restart",
                                         "ob1", "ob3", "ob5",
                                         "pd1", "pd2", "pd3", "pd4", "pd5",
                                         "pd6", "pd8", "pd-restart",
                                         "dg1", "dg2", "dg-restart",
                                         "et1", "et2", "et3", "et4", "et7",
                                         "et8", "et9", "et-orient", "et-rank",
                                         "et-restart",
                                         "et4-t", "et5", "et6", "et7-t",
                                         "et8-t", "et9-veg", "et-veg-restart", "et10",
                                         "ws1", "ws1n", "ws3", "ws4", "ws-restart",
                                         "ws-rank",
                                         "sx1", "sx2", "sx3", "sx4", "sx5",
                                         "sx-restart", "sx-rank"])
    parser.add_argument("--frehg", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--legacy", type=Path, required=True,
                        help="legacy benchmarks directory (goldens; never committed)")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--mode", choices=["strict", "default"], default="strict")
    parser.add_argument("--scenario", choices=["rain", "norain"], default="rain",
                        help="b5 gate scenario")
    parser.add_argument("--coupling", choices=["sync", "subcycled"], default="sync",
                        help="b5 gate synchronization mode")
    parser.add_argument("--variant", choices=["ss", "td"], default="ss",
                        help="b6 gate variant (no-tide / tidal)")
    parser.add_argument("--reuse-output", action="store_true",
                        help="b5 gate: re-apply the envelope checks to an existing "
                             "run in --work instead of re-running the hours-long "
                             "simulation")
    parser.add_argument("--t-end", type=float, default=None,
                        help="b5 gate: shorten the run to this horizon [s]; the "
                             "envelope checks clip to the covered window")
    parser.add_argument("--solver", choices=["bjacobi-icc", "amg", "gamg"], default=None,
                        help="g1 solver-invariance gate (v2 plan §2.3): run the gate "
                             "with this preconditioner selected via the solver YAML "
                             "block for both systems; pass/fail criteria unchanged")
    parser.add_argument("--transposed", action="store_true",
                        help="g5 gate: run the case as the y-z slice (the §9 Q4 "
                             "'transposed g5 slice' x-battery row); same criteria")
    parser.add_argument("--record", action="store_true",
                        help="default-bitwise gate: record the reference hashes "
                             "instead of comparing (Stage-3 plan W0.8)")
    parser.add_argument("--allow", action="append", default=[],
                        help="default-bitwise gate: dataset glob exempt from the "
                             "bitwise comparison (Stage-3 plan §2.6 declared "
                             "exceptions)")
    parser.add_argument("--record-massaudit", action="store_true",
                        help="default-bitwise (Stage-3 T9 ws5): record the T8 "
                             "mass_audit prefix reference instead of comparing; "
                             "run once with the frozen T8 binary")
    parser.add_argument("--pd2-lane", choices=["surface", "gw", "coupled"],
                        default=None,
                        help="pd2: run only this lane (Stage-3 T8 W8.4: the "
                             "p1 omp4 backend lanes run the coupled lane "
                             "alone); default = all three")
    parser.add_argument("--mat-type", choices=["aij", "aijkokkos"], default=None,
                        help="p1 backend-invariance gate (v2 plan §2B.3): run the "
                             "gate with this PETSc matrix/vector backend selected "
                             "via the solver YAML block; pass/fail criteria "
                             "unchanged")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)

    if args.solver is not None:
        if args.gate.startswith("rank-invariance") and args.mode == "strict":
            print("error: --solver conflicts with strict rank-invariance mode "
                  "(strict pins -pc_type jacobi on the PETSc command line)")
            return 2
        global SOLVER_OVERRIDE
        SOLVER_OVERRIDE = args.solver

    if args.mat_type is not None:
        global MAT_TYPE_OVERRIDE
        MAT_TYPE_OVERRIDE = args.mat_type

    if args.gate == "b1":
        return gate_b1(args)
    if args.gate == "b2":
        return gate_b2(args)
    if args.gate == "b3":
        return gate_b3(args)
    if args.gate == "b4":
        return gate_b4(args)
    if args.gate == "b5":
        return gate_b5(args)
    if args.gate == "b6":
        return gate_b6(args)
    if args.gate == "b6-restart":
        return gate_b6_restart(args)
    if args.gate == "b1-restart":
        return gate_b1_restart(args)
    if args.gate == "b2-restart":
        return gate_b2_restart(args)
    if args.gate == "b5-restart":
        return gate_b5_restart(args)
    if args.gate == "rank-invariance-b2":
        return gate_rank_invariance_b2(args)
    if args.gate == "rank-invariance-b5":
        return gate_rank_invariance_b5(args)
    if args.gate == "b6-gw-smoke":
        return gate_b6_gw_smoke(args)
    if args.gate == "smoke-b5":
        return gate_smoke_b5(args)
    if args.gate == "smoke-b6":
        return gate_smoke_b6(args)
    if args.gate == "smoke-g5":
        return gate_smoke_g5(args)
    if args.gate == "smoke-g8":
        return gate_smoke_g8(args)
    if args.gate == "g4":
        return gate_g4(args)
    if args.gate == "g5":
        return gate_g5(args)
    if args.gate == "g9":
        return gate_g9(args)
    if args.gate == "g10":
        return gate_g10(args)
    if args.gate == "wind-orient":
        return gate_wind_orient(args)
    if args.gate == "g6":
        return gate_g6(args)
    if args.gate == "g7":
        return gate_g7(args)
    if args.gate == "g8":
        return gate_g8(args)
    if args.gate == "heat-orient":
        return gate_heat_orient(args)
    if args.gate == "heat-restart":
        return gate_heat_restart(args)
    if args.gate == "rank-invariance-heat":
        return gate_rank_invariance_heat(args)
    if args.gate == "outflow-staircase":
        return gate_outflow_staircase(args)
    if args.gate == "swe-orient":
        return gate_swe_orient(args)
    if args.gate == "gw-orient":
        return gate_gw_orient(args)
    if args.gate == "transport-orient":
        return gate_transport_orient(args)
    if args.gate == "surface-bc-sides":
        return gate_surface_bc_sides(args)
    if args.gate == "wind-restart":
        return gate_wind_restart(args)
    if args.gate == "evap-restart":
        return gate_evap_restart(args)
    if args.gate == "b6-subcycled":
        return gate_b6_subcycled(args)
    if args.gate == "rank-invariance-tracer":
        return gate_rank_invariance_tracer(args)
    if args.gate == "rank-invariance-coupled-salt":
        return gate_rank_invariance_coupled_salt(args)
    if args.gate == "wind-drydown":
        return gate_wind_drydown(args)
    if args.gate == "terrain-heat":
        return gate_terrain_heat(args)
    if args.gate == "two-scalar-restart":
        return gate_two_scalar_restart(args)
    if args.gate == "evap-heat":
        return gate_evap_heat(args)
    if args.gate == "heat-subcycled":
        return gate_heat_subcycled(args)
    if args.gate == "default-bitwise":
        return gate_default_bitwise(args)
    if args.gate == "rl2":
        return gate_rl2(args)
    if args.gate == "rl3":
        return gate_rl3(args)
    if args.gate == "rl4":
        return gate_rl4(args)
    if args.gate == "rl5-orient":
        return gate_rl5_orient(args)
    if args.gate == "rl-restart":
        return gate_rl_restart(args)
    if args.gate == "ob1":
        return gate_ob1(args)
    if args.gate == "ob3":
        return gate_ob3(args)
    if args.gate == "ob5":
        return gate_ob5(args)
    if args.gate == "pd1":
        return gate_pd1(args)
    if args.gate == "pd2":
        return gate_pd2(args)
    if args.gate == "pd3":
        return gate_pd3(args)
    if args.gate == "pd4":
        return gate_pd4(args)
    if args.gate == "pd5":
        return gate_pd5(args)
    if args.gate == "pd6":
        return gate_pd6(args)
    if args.gate == "pd8":
        return gate_pd8(args)
    if args.gate == "pd-restart":
        return gate_pd_restart(args)
    if args.gate == "dg1":
        return gate_dg1(args)
    if args.gate == "dg2":
        return gate_dg2(args)
    if args.gate == "dg-restart":
        return gate_dg_restart(args)
    if args.gate == "et1":
        return gate_et1(args)
    if args.gate == "et2":
        return gate_et2(args)
    if args.gate == "et3":
        return gate_et3(args)
    if args.gate == "et4":
        return gate_et4(args)
    if args.gate == "et7":
        return gate_et7(args)
    if args.gate == "et8":
        return gate_et8(args)
    if args.gate == "et9":
        return gate_et9(args)
    if args.gate == "et-orient":
        return gate_et_orient(args)
    if args.gate == "et-rank":
        return gate_et_rank(args)
    if args.gate == "et-restart":
        return gate_et_restart(args)
    if args.gate == "et4-t":
        return gate_et4t(args)
    if args.gate == "et5":
        return gate_et5(args)
    if args.gate == "et6":
        return gate_et6(args)
    if args.gate == "et7-t":
        return gate_et7t(args)
    if args.gate == "et8-t":
        return gate_et8t(args)
    if args.gate == "et9-veg":
        return gate_et9veg(args)
    if args.gate == "et-veg-restart":
        return gate_et_veg_restart(args)
    if args.gate == "et10":
        return gate_et10(args)
    if args.gate == "ws1":
        return gate_ws1(args)
    if args.gate == "ws1n":
        return gate_ws1n(args)
    if args.gate == "ws3":
        return gate_ws3(args)
    if args.gate == "ws4":
        return gate_ws4(args)
    if args.gate == "ws-restart":
        return gate_ws_restart(args)
    if args.gate == "ws-rank":
        return gate_ws_rank(args)
    if args.gate == "sx1":
        return gate_sx1(args)
    if args.gate == "sx2":
        return gate_sx2(args)
    if args.gate == "sx3":
        return gate_sx3(args)
    if args.gate == "sx4":
        return gate_sx4(args)
    if args.gate == "sx5":
        return gate_sx1(args, case="sx5-uniform-heat.yaml", label="sx5")
    if args.gate == "sx-restart":
        return gate_sx_restart(args)
    if args.gate == "sx-rank":
        return gate_sx_rank(args)
    return gate_rank_invariance(args)


if __name__ == "__main__":
    sys.exit(main())
