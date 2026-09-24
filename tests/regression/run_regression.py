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
import os
import shutil
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


def apply_solver_override(case_dir: Path) -> None:
    """Write the SOLVER_OVERRIDE / MAT_TYPE_OVERRIDE into every staged YAML."""
    if SOLVER_OVERRIDE is None and MAT_TYPE_OVERRIDE is None:
        return
    for config in sorted(case_dir.glob("*.yaml")):
        with open(config, encoding="utf-8") as handle:
            doc = yaml.safe_load(handle)
        solver = doc.setdefault("solver", {})
        for system in ("surface", "groundwater"):
            entry = solver.setdefault(system, {})
            if SOLVER_OVERRIDE is not None:
                entry["preconditioner"] = SOLVER_OVERRIDE
            if MAT_TYPE_OVERRIDE is not None:
                entry["mat_type"] = MAT_TYPE_OVERRIDE
        with open(config, "w", encoding="utf-8") as handle:
            yaml.safe_dump(doc, handle, sort_keys=False)
        print(f"  solver override: {config.name} -> "
              f"preconditioner {SOLVER_OVERRIDE or '(default)'}"
              f" mat_type {MAT_TYPE_OVERRIDE or '(default)'}")


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
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout, end="")
    if result.returncode != 0:
        print(f"r1 run-record gate: FAIL for {config.name}")
        raise SystemExit(1)


def stage_case(repo: Path, case: str, workdir: Path) -> Path:
    """Copy a benchmark case directory into the scratch area."""
    source = repo / "benchmarks" / case
    target = workdir / case
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(source, target)
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


def gate_b6(args: argparse.Namespace) -> int:
    variant = args.variant
    case_dir = args.work / "b6-kuan"
    output = case_dir / f"out-{variant}" / "output.h5"
    if args.reuse_output and output.exists():
        print(f"b6 gate [{variant}]: reusing {output}")
    else:
        case_dir = stage_case(args.repo, "b6-kuan", args.work)
        config = case_dir / f"b6-kuan-{variant}.yaml"
        run_case(args.frehg, config, case_dir, args.mpiexec, args.ranks)

    with open(HERE / "tolerances" / "b6-kuan.yaml", encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)
    golden_dir = args.legacy / "b6-kuan" / f"out-{variant}-syncV4"
    exp_file = args.legacy / "b6-kuan" / "ex3_kuan_input" / (
        "Kuan_Exp_NoTide" if variant == "ss" else "Kuan_Exp_Tide")
    bath, dzf, nzf, actv = b6_geometry(args.legacy)
    ok = True

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
    if h_bad > 0.0:
        ok = False
    ds = np.abs(salt_m - salt_g)[finite]
    s_allow = np.maximum(g["salt_abs_floor"], g["salt_rel"] * np.abs(salt_g[finite]))
    s_bad = float(np.mean(ds > s_allow))
    print(f"b6-{variant} golden salinity: |ds| max {ds.max():.3f} psu, mean {ds.mean():.3f}, "
          f"{100.0 * s_bad:.2f} % of cells outside "
          f"(allowed {100.0 * g['salt_exceed_fraction']:.0f} %)")
    if s_bad > g["salt_exceed_fraction"]:
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
    if rel > spec["salt_mass_rel"]:
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
    with h5py.File(full_output, "r") as full, \
         h5py.File(restart_dir / "out" / "restarted.h5", "r") as part:
        for var in ["eta", "depth", "uu", "vv"]:
            for t in ["10800", "14400", "18000"]:
                a = full[f"/surface/{var}/{t}"][:]
                b = part[f"/surface/{var}/{t}"][:]
                scale = max(float(np.abs(a).max()), 1.0e-12)
                diff = float(np.abs(a - b).max()) / scale
                print(f"  restart {var}@{t}: max rel diff = {diff:.3e}")
                if diff > 1.0e-12:
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gate", choices=["b1", "b2", "b3", "b4", "b5", "b6", "b1-restart",
                                         "b2-restart", "b5-restart", "b6-restart",
                                         "rank-invariance", "rank-invariance-b2",
                                         "rank-invariance-b5", "b6-gw-smoke",
                                         "smoke-b5", "smoke-b6", "g4", "g5",
                                         "g9", "g10", "wind-orient", "g6", "g7",
                                         "g8", "heat-orient", "heat-restart",
                                         "rank-invariance-heat"])
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
    return gate_rank_invariance(args)


if __name__ == "__main__":
    sys.exit(main())
