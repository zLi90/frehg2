#!/usr/bin/env python3
"""Frehg2 scaling and solver-performance harness (v2 plan §7: the s-gates; §2.3: g2).

Modes, selected by --gate:

  strong        (s1) strong scaling over --ranks on a fixed grid; PASS iff
                parallel efficiency >= 70 % at 4 ranks (plan §10 P5 / A23,
                permanent per v2 plan §7.2 s1).
  weak          (s2/g3) weak scaling: the synthetic grid's nx is scaled by
                the rank count so per-rank work is constant (Kollet-2010
                protocol). PASS iff solver-time efficiency >= the hard
                floor (0.5; 0.8 is the target, warned) and mean solver
                iterations grow <= 15 % from the smallest to the largest
                rank count. Setup-time fraction is reported.
  threads       (s3) OpenMP thread scaling at 1 rank over --threads. The
                PETSc solve is MPI-only, so the assertion is honest about
                Amdahl: PASS iff *kernel* time (simulation minus solve)
                reaches >= 60 % efficiency at 4 threads AND the final
                mass-audit volume matches the 1-thread run to 1e-8
                relative. Solve-time flatness is warned, not gated
                (v2 plan risk register: timing soft until real runners).
  hybrid        (s4) rank x thread placements of 4 PEs: 4x1, 2x2, 1x4.
                PASS iff every placement's final mass-audit volume matches
                4x1 to 1e-8 relative; timings recorded.
  iters         (g2) iteration-flatness-vs-ranks gate (v2 plan §2.3): runs
                the synthetic case over --ranks with --solver, parses the
                per-system solver summaries, and PASSes iff
                mean_iters(largest) <= ratio_cap * mean_iters(smallest)
                for both systems, and (when tolerances record golden
                counts) mean <= recorded + max(band_abs, band_rel*recorded).
  perf-baseline (§7.1.2) per-module timer regression vs a tracked baseline
                (min over --repeats): FAIL if any module is > 25 % slower,
                warn > 10 %. --seed-baseline (re)writes the baseline file.
  none          measurement only (record, no assertions).

Common machinery: every run pins FI_PROVIDER=tcp (A19); OMP_NUM_THREADS is
1 unless the mode says otherwise; rank counts run back-to-back in the order
given (thermal consistency, A23); --repeats keeps the minimum simulation
time (the fanless-M3 measurement lottery, A23). --solver injects the v2
solver YAML block (preconditioner for both systems); --json writes the full
measurement artifact.

The gate value is the "simulation" timer (whole time loop, max across
ranks). Solver iteration statistics come from the end-of-run
"solver summary" lines (v2 plan §2.2 telemetry).
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import yaml

TIMER_LINE = re.compile(
    r"^\s{2}(?P<path>\S+)\s+(?P<count>\d+)\s+(?P<min>[\d.]+)\s+(?P<mean>[\d.]+)\s+(?P<max>[\d.]+)\s*$")

# End-of-run per-system telemetry line (v2 plan §2.2; printed by the driver
# on rank 0). Matched anywhere in the line so the logger prefix is free.
SOLVER_LINE = re.compile(
    r"solver summary (?P<system>fs|gw): solves=(?P<solves>\d+) "
    r"iters_mean=(?P<mean>[\d.]+) iters_max=(?P<max>\d+) "
    r"rebuilds=(?P<rebuilds>\d+) retries=(?P<retries>\d+) "
    r"setup_s=(?P<setup>[\d.eE+-]+) solve_s=(?P<solve>[\d.eE+-]+)")

# The amendment-A23 scaling-gate case: the b5 coupled physics (rain onto a
# draining box over a water table — the same predictor/corrector/
# reallocation/exchange kernels) at 16x the b5 cell count, so 4 ranks have
# real work. Fixed common step; closed box (no BC files needed). nx/ny are
# parameters so the weak-scaling mode can hold per-rank work constant.
SYNTHETIC_CONFIG = """\
simulation: {{id: p5-scaling-synthetic, title: "Scaling harness case (A23 physics)"}}
domain:
  nx: {nx}
  ny: {ny}
  nz: 25
  dx: 1.0
  dy: 1.0
  dz: 0.2
  bottom_elevation: {{constant: 0.0}}
  follow_terrain: false
time: {{dt: {dt}, t_end: {t_end}, output_interval: {output_interval}}}
modules: {{surface_water: true, groundwater: true, transport: false}}
surface_water:
  friction:
    law: manning
    coefficient: {{constant: 0.03}}
    thin_layer_depth: 0.01
  viscosity: {{x: 1.0e-6, y: 1.0e-6}}
  min_depth: 1.0e-6
  wetting_face_depth: 1.0e-3
  rainfall: {{constant: 5.5e-6}}
groundwater:
  timestep: {{dt_init: {dt}, dt_min: {dt}, dt_max: {dt}}}
  specific_storage: 1.0e-5
soil:
  types:
    - {{name: loam, ksx: 6.94e-5, ksy: 6.94e-5, ksz: 6.94e-5,
       theta_s: 0.4, theta_r: 0.08, vg_alpha: 1.0, vg_n: 2.0}}
  map: {{constant: loam}}
coupling: {{mode: sync}}
initial_conditions:
  surface: {{eta: {{constant: -10.0}}}}
  groundwater: {{water_table: {{constant: -2.5}}}}
output:
  filename: out/output.h5
  variables: {{surface: [depth]}}
  checkpoint: {{interval: 0}}
"""


def apply_solver(config: Path, solver: str | None, mat_type: str | None = None) -> None:
    """Inject the v2 solver block (both systems) into a staged YAML."""
    if solver is None and mat_type is None:
        return
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    block = doc.setdefault("solver", {})
    for system in ("surface", "groundwater"):
        entry = block.setdefault(system, {})
        if solver is not None:
            entry["preconditioner"] = solver
        if mat_type is not None:
            entry["mat_type"] = mat_type
    with open(config, "w", encoding="utf-8") as handle:
        yaml.safe_dump(doc, handle, sort_keys=False)


def stage_synthetic(work: Path, t_end: float, dt: float, output_interval: float,
                    nx: int = 404, ny: int = 220, solver: str | None = None,
                    mat_type: str | None = None,
                    tag: str = "synthetic") -> Path:
    target = work / tag
    if target.exists():
        shutil.rmtree(target)
    (target / "out").mkdir(parents=True)
    config = target / f"{tag}.yaml"
    config.write_text(SYNTHETIC_CONFIG.format(nx=nx, ny=ny, dt=dt, t_end=t_end,
                                              output_interval=output_interval),
                      encoding="utf-8")
    apply_solver(config, solver, mat_type)
    return config


def stage(repo: Path, work: Path, t_end: float, dt: float, output_interval: float,
          solver: str | None = None, mat_type: str | None = None) -> Path:
    source = repo / "benchmarks" / "b5-vcatchment"
    target = work / "b5-vcatchment"
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(source, target)
    config = target / "b5-vcatchment.yaml"
    with open(config, encoding="utf-8") as handle:
        doc = yaml.safe_load(handle)
    doc["time"]["dt"] = dt
    doc["time"]["t_end"] = t_end
    doc["time"]["output_interval"] = output_interval
    doc["groundwater"]["timestep"]["dt_init"] = dt
    doc["groundwater"]["timestep"]["dt_min"] = dt
    doc["groundwater"]["timestep"]["dt_max"] = dt
    doc["output"]["checkpoint"] = {"interval": 0}
    with open(config, "w", encoding="utf-8") as handle:
        yaml.safe_dump(doc, handle, sort_keys=False)
    apply_solver(config, solver, mat_type)
    return config


def final_volume(output: Path) -> float:
    """Last mass-audit volume (m^3): the s3/s4 cross-configuration check."""
    import h5py  # deferred: only the correctness modes need it
    with h5py.File(output, "r") as handle:
        table = handle["/monitor/mass_audit"][:]
    return float(table[-1, 1])


def run_one(frehg: Path, config: Path, mpiexec: Path, ranks: int, work: Path,
            threads: int = 1, extra_args: list[str] | None = None) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PROC_BIND="false",
               FI_PROVIDER="tcp")
    log = work / f"scaling_n{ranks}_t{threads}.log"
    cmd = [str(mpiexec), "-n", str(ranks), str(frehg), str(config)] + (extra_args or [])
    started = time.perf_counter()
    with open(log, "w", encoding="utf-8") as handle:
        result = subprocess.run(cmd, cwd=config.parent, stdout=handle,
                                stderr=subprocess.STDOUT, env=env, check=False)
    wall = time.perf_counter() - started
    if result.returncode != 0:
        print(f"error: frehg exited {result.returncode} at n={ranks} t={threads}; "
              f"log tail ({log}):")
        print("\n".join(log.read_text().splitlines()[-25:]))
        raise SystemExit(1)
    # Primary source: the run record (v2 plan §2A.4). The stdout scrape below
    # remains as a fallback for binaries predating the record.
    timers = {}
    solver = {}
    record_path = record_of(config)
    if record_path.exists():
        record = yaml.safe_load(record_path.read_text())
        for path, entry in (record.get("timers") or {}).items():
            timers[path] = {"min": float(entry["min_s"]), "mean": float(entry["mean_s"]),
                            "max": float(entry["max_s"]), "count": int(entry["count"])}
        for system, entry in (record.get("solver") or {}).items():
            if isinstance(entry, dict) and "iters_mean" in entry:
                solver[system] = entry
    if not timers:
        for line in log.read_text().splitlines():
            match = TIMER_LINE.match(line)
            if match:
                timers[match["path"]] = {k: float(match[k]) for k in ("min", "mean", "max")}
                timers[match["path"]]["count"] = int(match["count"])
    if not solver:
        for line in log.read_text().splitlines():
            smatch = SOLVER_LINE.search(line)
            if smatch:
                solver[smatch["system"]] = {
                    "solves": int(smatch["solves"]),
                    "iters_mean": float(smatch["mean"]),
                    "iters_max": int(smatch["max"]),
                    "rebuilds": int(smatch["rebuilds"]),
                    "retries": int(smatch["retries"]),
                    "setup_s": float(smatch["setup"]),
                    "solve_s": float(smatch["solve"]),
                }
    if "simulation" not in timers:
        print(f"error: no 'simulation' timer in {record_path} or {log}")
        raise SystemExit(1)
    return {"ranks": ranks, "threads": threads, "wall_s": wall, "timers": timers,
            "solver": solver, "output": str(config.parent / "out" / "output.h5"),
            "record": str(record_path)}


def record_of(config: Path) -> Path:
    """Path of the run record beside the configured HDF5 output."""
    doc = yaml.safe_load(config.read_text())
    output = Path(doc["output"]["filename"])
    return (config.parent / output).parent / "run-record.yaml"


def best_of(frehg: Path, config: Path, mpiexec: Path, ranks: int, work: Path,
            repeats: int, threads: int = 1, label: str = "") -> dict:
    best = None
    for attempt in range(repeats):
        print(f"running {label} at n={ranks} threads={threads} "
              f"(attempt {attempt + 1}/{repeats}) ...", flush=True)
        res = run_one(frehg, config, mpiexec, ranks, work, threads=threads)
        sim = res["timers"]["simulation"]["max"]
        print(f"  simulation {sim:.2f} s", flush=True)
        if best is None or sim < best["timers"]["simulation"]["max"]:
            best = res
    return best


def solve_seconds(res: dict) -> float:
    """Total time in the linear solves (all timer paths ending in /solve)."""
    return sum(t["max"] for p, t in res["timers"].items()
               if p == "solve" or p.endswith("/solve"))


def require_solver(res: dict, systems: tuple[str, ...] = ("fs", "gw")) -> bool:
    missing = [s for s in systems if s not in res["solver"]]
    if missing:
        print(f"error: no 'solver summary' telemetry for {missing} at "
              f"n={res['ranks']} — is the v2 solver telemetry implemented?")
        return False
    return True


def print_strong_table(results: list[dict], label: str) -> None:
    base = results[0]
    base_sim = base["timers"]["simulation"]["max"]
    base_ranks = base["ranks"]
    print(f"\nstrong scaling, {label} (baseline n={base_ranks}):")
    print(f"  {'ranks':>5} {'simulation[s]':>14} {'wall[s]':>10} "
          f"{'speedup':>8} {'efficiency':>11}")
    for res in results:
        sim = res["timers"]["simulation"]["max"]
        speedup = base_sim / sim
        eff = speedup * base_ranks / res["ranks"]
        res["speedup"] = speedup
        res["efficiency"] = eff
        print(f"  {res['ranks']:>5} {sim:>14.2f} {res['wall_s']:>10.2f} "
              f"{speedup:>8.2f} {eff:>10.1%}")


def print_module_table(results: list[dict], key: str = "ranks") -> None:
    print("\nper-module timers (max over ranks) [s]:")
    paths = sorted({p for res in results for p in res["timers"]})
    header = "  " + f"{'section':<40}" + "".join(
        f"{key[0] + '=' + str(r[key]):>12}" for r in results)
    print(header)
    for path in paths:
        row = f"  {path:<40}"
        for res in results:
            timer = res["timers"].get(path)
            row += f"{timer['max']:>12.2f}" if timer else f"{'-':>12}"
        print(row)


def gate_strong(results: list[dict]) -> int:
    four = next((res for res in results if res["ranks"] == 4), None)
    if four is None:
        return 0
    ok = four["efficiency"] >= 0.70
    print(f"s1 strong gate (4-rank efficiency >= 70 %): {four['efficiency']:.1%} — "
          f"{'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def gate_weak(results: list[dict]) -> int:
    base = results[0]
    last = results[-1]
    base_sim = base["timers"]["simulation"]["max"]
    ok = True
    print(f"\nweak scaling (per-rank work constant, baseline n={base['ranks']}):")
    print(f"  {'ranks':>5} {'simulation[s]':>14} {'efficiency':>11} "
          f"{'fs_it':>7} {'gw_it':>7} {'setup%':>7}")
    for res in results:
        if not require_solver(res):
            return 1
        sim = res["timers"]["simulation"]["max"]
        eff = base_sim / sim
        res["efficiency"] = eff
        setup = sum(res["solver"][s]["setup_s"] for s in ("fs", "gw"))
        solve = sum(res["solver"][s]["solve_s"] for s in ("fs", "gw"))
        frac = setup / (setup + solve) if setup + solve > 0 else 0.0
        print(f"  {res['ranks']:>5} {sim:>14.2f} {eff:>10.1%} "
              f"{res['solver']['fs']['iters_mean']:>7.1f} "
              f"{res['solver']['gw']['iters_mean']:>7.1f} {frac:>6.1%}")
    eff = last["efficiency"]
    if eff < 0.5:
        print(f"s2 weak gate: efficiency {eff:.1%} < 50 % hard floor — FAIL")
        ok = False
    elif eff < 0.8:
        print(f"s2 weak gate: efficiency {eff:.1%} below the 80 % target "
              f"(>= 50 % floor holds) — WARN")
    else:
        print(f"s2 weak gate: efficiency {eff:.1%} — ok")
    for system in ("fs", "gw"):
        growth = last["solver"][system]["iters_mean"] / \
            max(base["solver"][system]["iters_mean"], 1e-12) - 1.0
        passed = growth <= 0.15
        print(f"s2 weak gate: {system} iteration growth {growth:+.1%} "
              f"(allowed +15 %) — {'ok' if passed else 'FAIL'}")
        ok = ok and passed
    print("s2 weak gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_threads(results: list[dict]) -> int:
    base = results[0]
    if base["threads"] != 1:
        print("error: threads mode needs a 1-thread baseline first in --threads")
        return 2
    base_kernel = base["timers"]["simulation"]["max"] - solve_seconds(base)
    base_solve = solve_seconds(base)
    base_volume = final_volume(Path(base["output"]))
    ok = True
    print("\nOpenMP thread scaling at 1 rank (kernel = simulation - solve):")
    print(f"  {'threads':>7} {'simulation[s]':>14} {'kernel[s]':>10} "
          f"{'kernel eff':>11} {'solve[s]':>9}")
    for res in results:
        kernel = res["timers"]["simulation"]["max"] - solve_seconds(res)
        eff = base_kernel / (res["threads"] * kernel) if kernel > 0 else 0.0
        res["kernel_efficiency"] = eff
        print(f"  {res['threads']:>7} {res['timers']['simulation']['max']:>14.2f} "
              f"{kernel:>10.2f} {eff:>10.1%} {solve_seconds(res):>9.2f}")
        volume = final_volume(Path(res["output"]))
        rel = abs(volume - base_volume) / max(abs(base_volume), 1e-30)
        if rel > 1e-8:
            print(f"s3 threads gate: final volume at t={res['threads']} differs "
                  f"from 1-thread by {rel:.2e} (allowed 1e-8) — FAIL")
            ok = False
    four = next((res for res in results if res["threads"] == 4), None)
    if four is not None:
        passed = four["kernel_efficiency"] >= 0.60
        print(f"s3 threads gate (kernel efficiency at 4 threads >= 60 %): "
              f"{four['kernel_efficiency']:.1%} — {'ok' if passed else 'FAIL'}")
        ok = ok and passed
        drift = solve_seconds(four) / max(base_solve, 1e-12) - 1.0
        if abs(drift) > 0.25:
            # V2-A6: flat only holds on the host-aij lane (the ctest s3
            # entry); with mat_type aijkokkos the solve threads by design —
            # that lane is gated by p2 instead.
            print(f"s3 threads gate: solve time drifted {drift:+.1%} with threads "
                  f"(expected flat on the aij lane; aijkokkos is gated by p2) — WARN")
    print("s3 threads gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_p2_solver_threads(results: list[dict]) -> int:
    """p2 (v2 plan §2B.3): OpenMP thread scaling of the aijkokkos solve.

    Hard: gw iteration means constant within +-2 % across thread counts
    (thread parallelism must not perturb the algebra — the l1scaled-Jacobi
    smoother is thread-invariant by construction), and gw solve wall time at
    4 threads strictly below 1 thread (the solve must actually thread; the
    old MPI-only design fails exactly this). Recorded, not gated: solve
    thread efficiency, setup share, total-runtime efficiency.
    """
    base = results[0]
    if base["threads"] != 1:
        print("error: p2 needs a 1-thread baseline first in --threads")
        return 2
    if not require_solver(base, ("gw",)):
        return 2
    base_iters = base["solver"]["gw"]["iters_mean"]
    base_solve = float(base["solver"]["gw"]["solve_s"])
    base_volume = final_volume(Path(base["output"]))
    ok = True
    print("\np2 solver thread scaling at 1 rank (aijkokkos):")
    print(f"  {'threads':>7} {'gw iters':>9} {'gw solve[s]':>12} "
          f"{'solve eff':>10} {'simulation[s]':>14}")
    for res in results:
        iters = res["solver"]["gw"]["iters_mean"]
        solve = float(res["solver"]["gw"]["solve_s"])
        eff = base_solve / (res["threads"] * solve) if solve > 0 else 0.0
        print(f"  {res['threads']:>7} {iters:>9.2f} {solve:>12.2f} "
              f"{eff:>9.1%} {res['timers']['simulation']['max']:>14.2f}")
        drift = abs(iters - base_iters) / max(base_iters, 1e-12)
        if drift > 0.02:
            print(f"p2: gw iterations at t={res['threads']} drift {drift:.2%} "
                  f"from 1-thread (allowed 2 %) — FAIL")
            ok = False
        volume = final_volume(Path(res["output"]))
        rel = abs(volume - base_volume) / max(abs(base_volume), 1e-30)
        if rel > 1e-8:
            print(f"p2: final volume at t={res['threads']} differs from 1-thread "
                  f"by {rel:.2e} (allowed 1e-8) — FAIL")
            ok = False
    four = next((res for res in results if res["threads"] == 4), None)
    if four is not None:
        solve4 = float(four["solver"]["gw"]["solve_s"])
        passed = solve4 < base_solve
        print(f"p2 (gw solve at 4 threads < 1 thread): {solve4:.2f} vs "
              f"{base_solve:.2f} s — {'ok' if passed else 'FAIL'}")
        ok = ok and passed
    print("p2 solver-threads gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_p3_hybrid(results: list[dict], total_pes: int) -> int:
    """p3 (v2 plan §2B.3): hybrid MPI+OpenMP placements at fixed total PEs
    on the aijkokkos path. Hard: final volumes within 1e-8 across
    placements, and gw iteration means within 10 % of the all-MPI placement
    (a placement that changes the algebra is a decomposition bug, not a
    tuning result). Timings recorded; the winning ratio is printed."""
    base = results[0]
    if not require_solver(base, ("gw",)):
        return 2
    base_volume = final_volume(Path(base["output"]))
    base_iters = base["solver"]["gw"]["iters_mean"]
    ok = True
    best = None
    print(f"\np3 hybrid placements of {total_pes} PEs (aijkokkos):")
    print(f"  {'ranks':>5} {'threads':>7} {'gw iters':>9} {'simulation[s]':>14}")
    for res in results:
        sim = res["timers"]["simulation"]["max"]
        iters = res["solver"]["gw"]["iters_mean"]
        print(f"  {res['ranks']:>5} {res['threads']:>7} {iters:>9.2f} {sim:>14.2f}")
        if best is None or sim < best["timers"]["simulation"]["max"]:
            best = res
        volume = final_volume(Path(res["output"]))
        rel = abs(volume - base_volume) / max(abs(base_volume), 1e-30)
        if rel > 1e-8:
            print(f"p3: final volume at {res['ranks']}x{res['threads']} differs "
                  f"from {base['ranks']}x{base['threads']} by {rel:.2e} "
                  f"(allowed 1e-8) — FAIL")
            ok = False
        drift = abs(iters - base_iters) / max(base_iters, 1e-12)
        if drift > 0.10:
            print(f"p3: gw iterations at {res['ranks']}x{res['threads']} drift "
                  f"{drift:.1%} from the all-MPI placement (allowed 10 %) — FAIL")
            ok = False
    print(f"p3 fastest placement: {best['ranks']} ranks x {best['threads']} threads "
          f"({best['timers']['simulation']['max']:.2f} s)")
    print("p3 hybrid gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_hybrid(results: list[dict]) -> int:
    base = results[0]
    base_volume = final_volume(Path(base["output"]))
    ok = True
    print("\nhybrid placements of 4 PEs:")
    print(f"  {'ranks':>5} {'threads':>7} {'simulation[s]':>14}")
    for res in results:
        print(f"  {res['ranks']:>5} {res['threads']:>7} "
              f"{res['timers']['simulation']['max']:>14.2f}")
        volume = final_volume(Path(res["output"]))
        rel = abs(volume - base_volume) / max(abs(base_volume), 1e-30)
        if rel > 1e-8:
            print(f"s4 hybrid gate: final volume at {res['ranks']}x{res['threads']} "
                  f"differs from {base['ranks']}x{base['threads']} by {rel:.2e} "
                  f"(allowed 1e-8) — FAIL")
            ok = False
    print("s4 hybrid gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_iters(results: list[dict], tolerances: dict, solver: str) -> int:
    smallest, largest = results[0], results[-1]
    ratio_cap = float(tolerances.get("ratio_cap", 1.10))
    band_abs = float(tolerances.get("band_abs", 3))
    band_rel = float(tolerances.get("band_rel", 0.10))
    recorded = (tolerances.get("recorded") or {}).get(solver) or {}
    ok = True
    print(f"\ng2 iteration-flatness gate ({solver}, n={smallest['ranks']} -> "
          f"n={largest['ranks']}):")
    for res in results:
        if not require_solver(res):
            return 1
        print(f"  n={res['ranks']:>2}: fs mean {res['solver']['fs']['iters_mean']:.1f} "
              f"(max {res['solver']['fs']['iters_max']}, rebuilds "
              f"{res['solver']['fs']['rebuilds']}), gw mean "
              f"{res['solver']['gw']['iters_mean']:.1f} (max "
              f"{res['solver']['gw']['iters_max']}, rebuilds "
              f"{res['solver']['gw']['rebuilds']})")
    for system in ("fs", "gw"):
        small = smallest["solver"][system]["iters_mean"]
        large = largest["solver"][system]["iters_mean"]
        ratio = large / max(small, 1e-12)
        passed = ratio <= ratio_cap
        print(f"  {system}: iters({largest['ranks']})/iters({smallest['ranks']}) = "
              f"{ratio:.3f} (allowed {ratio_cap:.2f}) — {'ok' if passed else 'FAIL'}")
        ok = ok and passed
        cap = recorded.get(f"{system}_mean")
        if cap is not None:
            band = max(band_abs, band_rel * float(cap))
            passed = large <= float(cap) + band
            print(f"  {system}: mean {large:.1f} vs recorded {cap} + band {band:.1f} "
                  f"— {'ok' if passed else 'FAIL'}")
            ok = ok and passed
        else:
            print(f"  {system}: no recorded golden count for '{solver}' — measured "
                  f"mean {large:.1f} (seed tolerances/g2.yaml with it)")
    print("g2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_timer_coverage(results: list[dict]) -> int:
    """r2 (v2 plan §2A.3): the timer tree accounts for the time loop, the new
    halo/io sections are live in parallel runs, and the record is free."""
    ok = True
    for res in results:
        timers = res["timers"]
        sim = timers["simulation"]["max"]
        # Top-level cover of the time loop: sections under simulation/ with
        # no ancestor *section* between them and simulation. (Frame names may
        # embed slashes — "swe/free_surface" is one frame — so ancestry is
        # judged against the actual section set, not by slash counting.)
        paths = set(timers)
        children = {}
        for p, t in timers.items():
            if not p.startswith("simulation/"):
                continue
            inner = p[len("simulation/"):]
            has_ancestor = any(
                p != q and q.startswith("simulation/") and p.startswith(q + "/")
                for q in paths)
            if not has_ancestor:
                children[inner] = t
        covered = sum(t["max"] for t in children.values())
        coverage = covered / sim if sim > 0 else 0.0
        print(f"n={res['ranks']}: simulation {sim:.2f} s, top-level cover "
              f"{covered:.2f} s -> coverage {coverage:.1%}")
        for p, t in sorted(children.items()):
            print(f"    {p:<40} {t['max']:>9.2f} s")
        if res["ranks"] == 1:
            passed = coverage >= 0.90
            print(f"r2: coverage at 1 rank >= 90 %: {coverage:.1%} — "
                  f"{'ok' if passed else 'FAIL'}")
            ok = ok and passed
        halo = sum(t["max"] for p, t in timers.items()
                   if p.split("/")[-1] == "halo")
        io = sum(t["max"] for p, t in timers.items()
                 if "io/" in p and p.split("io/")[-1] in ("output", "checkpoint"))
        rr = sum(t["max"] for p, t in timers.items() if p.endswith("io/run_record"))
        print(f"    halo total {halo:.3f} s, io total {io:.3f} s, "
              f"run-record write {rr:.3f} s")
        if res["ranks"] > 1:
            passed = halo > 0.0 and io > 0.0
            print(f"r2: halo and io sections nonzero at {res['ranks']} ranks — "
                  f"{'ok' if passed else 'FAIL'}")
            ok = ok and passed
        overhead = rr / sim if sim > 0 else 0.0
        passed = overhead <= 0.01
        print(f"r2: run-record write {overhead:.2%} of the time loop (<= 1 %) — "
              f"{'ok' if passed else 'FAIL'}")
        ok = ok and passed
    print("r2 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def gate_perf_baseline(result: dict, baseline_path: Path, seed: bool) -> int:
    metrics = {p: t["max"] for p, t in result["timers"].items()}
    if seed:
        payload = {
            "machine": platform.platform(),
            "processor": platform.processor(),
            "note": "min-over-repeats 'max over ranks' timer seconds for the fixed "
                    "perf-baseline case (v2 plan §7.1.2); update only with a "
                    "justifying commit message",
            "timers_s": {p: round(v, 4) for p, v in metrics.items()},
        }
        baseline_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"perf baseline seeded at {baseline_path}")
        return 0
    if not baseline_path.exists():
        print(f"error: perf baseline {baseline_path} missing — run with "
              f"--seed-baseline once (v2 plan §7.1.2)")
        return 1
    baseline = json.loads(baseline_path.read_text())["timers_s"]
    ok = True
    print(f"\nper-module perf regression vs {baseline_path.name} "
          f"(fail > +25 %, warn > +10 %):")
    for path, seconds in sorted(metrics.items()):
        ref = baseline.get(path)
        if ref is None or ref < 0.05:
            continue  # new or too-small-to-time section: informational only
        drift = seconds / ref - 1.0
        flag = "ok"
        if drift > 0.25:
            flag = "FAIL"
            ok = False
        elif drift > 0.10:
            flag = "WARN"
        print(f"  {path:<40} {seconds:>9.2f} s vs {ref:>9.2f} s ({drift:+.1%}) {flag}")
    print("perf-baseline gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--frehg", type=Path, required=True)
    parser.add_argument("--mpiexec", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--case", choices=["synthetic", "b5"], default="synthetic")
    parser.add_argument("--gate",
                        choices=["strong", "weak", "threads", "hybrid", "iters",
                                 "p2-solver-threads", "p3-hybrid",
                                 "perf-baseline", "timer-coverage", "none"],
                        default="strong")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4],
                        help="thread counts for --gate threads (needs 1 first)")
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--t-end", type=float, default=None)
    parser.add_argument("--dt", type=float, default=2.0)
    parser.add_argument("--output-interval", type=float, default=None)
    parser.add_argument("--solver", choices=["bjacobi-icc", "amg", "gamg"], default=None,
                        help="inject the v2 solver block (both systems)")
    parser.add_argument("--mat-type", choices=["aij", "aijkokkos"], default=None,
                        help="inject solver.<system>.mat_type into the staged "
                             "YAML (p2/p3 force aijkokkos themselves)")
    parser.add_argument("--tolerances", type=Path, default=None,
                        help="g2: YAML with ratio_cap/band/recorded counts")
    parser.add_argument("--baseline", type=Path, default=None,
                        help="perf-baseline: tracked JSON baseline file")
    parser.add_argument("--seed-baseline", action="store_true")
    parser.add_argument("--json", type=Path, default=None,
                        help="write the full measurement artifact here")
    args = parser.parse_args()

    args.work.mkdir(parents=True, exist_ok=True)
    t_end = args.t_end if args.t_end is not None else \
        (60.0 if args.case == "synthetic" else 7200.0)
    output_interval = args.output_interval if args.output_interval is not None else \
        min(t_end, 3600.0)

    results: list[dict] = []
    if args.gate == "p2-solver-threads":
        # p2 (v2 plan §2B.3): the solve itself must thread on the aijkokkos
        # path, with thread-invariant algebra. Solver default: amg — the
        # scalable lane the gate exists to protect.
        threads = args.threads if args.threads[0] == 1 else [1] + args.threads
        config = stage_synthetic(args.work, t_end, args.dt, output_interval,
                                 solver=args.solver or "amg",
                                 mat_type="aijkokkos", tag="p2")
        label = "p2 synthetic 404x220x25 aijkokkos"
        for t in threads:
            results.append(best_of(args.frehg, config, args.mpiexec, 1, args.work,
                                   args.repeats, threads=t, label=label))
        print_module_table(results, key="threads")
        code = gate_p2_solver_threads(results)
    elif args.gate == "p3-hybrid":
        # p3 (v2 plan §2B.3): fixed total PEs, all rank x thread splits.
        total = max(args.ranks) if args.ranks else 8
        placements = []
        r = total
        while r >= 1:
            placements.append((r, total // r))
            r //= 2
        config = stage_synthetic(args.work, t_end, args.dt, output_interval,
                                 solver=args.solver or "amg",
                                 mat_type="aijkokkos", tag="p3")
        label = f"p3 synthetic 404x220x25 aijkokkos {total} PEs"
        for ranks, t in placements:
            results.append(best_of(args.frehg, config, args.mpiexec, ranks, args.work,
                                   args.repeats, threads=t, label=label))
        code = gate_p3_hybrid(results, total)
    elif args.gate == "threads":
        threads = args.threads if args.threads[0] == 1 else [1] + args.threads
        config = stage_synthetic(args.work, t_end, args.dt, output_interval,
                                 solver=args.solver, mat_type=args.mat_type)
        label = "synthetic 404x220x25"
        for t in threads:
            results.append(best_of(args.frehg, config, args.mpiexec, 1, args.work,
                                   args.repeats, threads=t, label=label))
        print_module_table(results, key="threads")
        code = gate_threads(results)
    elif args.gate == "hybrid":
        config = stage_synthetic(args.work, t_end, args.dt, output_interval,
                                 solver=args.solver, mat_type=args.mat_type)
        label = "synthetic 404x220x25"
        for ranks, t in ((4, 1), (2, 2), (1, 4)):
            results.append(best_of(args.frehg, config, args.mpiexec, ranks, args.work,
                                   args.repeats, threads=t, label=label))
        code = gate_hybrid(results)
    elif args.gate == "weak":
        # Per-rank unit: 202x110x25 (~555k cells); nx scales with ranks.
        for ranks in args.ranks:
            config = stage_synthetic(args.work, t_end, args.dt, output_interval,
                                     nx=202 * ranks, ny=110, solver=args.solver,
                                     mat_type=args.mat_type,
                                     tag=f"weak-n{ranks}")
            label = f"weak unit 202x110x25 x{ranks}"
            results.append(best_of(args.frehg, config, args.mpiexec, ranks, args.work,
                                   args.repeats, label=label))
        code = gate_weak(results)
    elif args.gate == "timer-coverage":
        # r2 (v2 plan §2A.3): short coupled run at 1 and 4 ranks; assertions
        # read the run record, not stdout.
        config = stage_synthetic(args.work, 20.0, args.dt, 20.0,
                                 solver=args.solver, mat_type=args.mat_type, tag="timer-coverage")
        for ranks in (1, 4):
            results.append(best_of(args.frehg, config, args.mpiexec, ranks, args.work,
                                   args.repeats, label="timer-coverage 404x220x25"))
        code = gate_timer_coverage(results)
    elif args.gate == "perf-baseline":
        if args.baseline is None:
            print("error: --gate perf-baseline needs --baseline")
            return 2
        # Fixed short case; min over repeats (>= 3 recommended).
        config = stage_synthetic(args.work, 20.0, args.dt, 20.0,
                                 solver=args.solver, mat_type=args.mat_type, tag="perf-baseline")
        result = best_of(args.frehg, config, args.mpiexec, 1, args.work,
                         max(args.repeats, 3), label="perf-baseline 404x220x25")
        results = [result]
        code = gate_perf_baseline(result, args.baseline, args.seed_baseline)
    else:  # strong / iters / none: fixed grid over --ranks
        if args.case == "synthetic":
            config = stage_synthetic(args.work, t_end, args.dt, output_interval,
                                     solver=args.solver, mat_type=args.mat_type)
            label = "synthetic 404x220x25 (A23 gate)"
        else:
            config = stage(args.repo, args.work, t_end, args.dt, output_interval,
                           solver=args.solver, mat_type=args.mat_type)
            label = "b5 rain/sync"
        for ranks in args.ranks:
            results.append(best_of(args.frehg, config, args.mpiexec, ranks, args.work,
                                   args.repeats, label=label))
        print_strong_table(results, label)
        print_module_table(results)
        if args.gate == "strong":
            code = gate_strong(results)
        elif args.gate == "iters":
            tolerances = {}
            if args.tolerances is not None and args.tolerances.exists():
                tolerances = yaml.safe_load(args.tolerances.read_text()) or {}
            code = gate_iters(results, tolerances, args.solver or "bjacobi-icc")
        else:
            code = 0

    if args.json is not None:
        payload = {
            "gate": args.gate,
            "solver": args.solver,
            "repeats": args.repeats,
            "selection": "minimum simulation time over repeats",
            "results": [
                {k: res[k] for k in
                 ("ranks", "threads", "wall_s", "solver")} |
                {"simulation_s": res["timers"]["simulation"]["max"],
                 "timers_max_s": {p: t["max"] for p, t in res["timers"].items()},
                 **({"efficiency": round(res["efficiency"], 4)}
                    if "efficiency" in res else {})}
                for res in results
            ],
        }
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(f"\nartifact written to {args.json}")

    # Legacy artifact kept for the A23 workflow.
    out = args.work / "scaling.yaml"
    with open(out, "w", encoding="utf-8") as handle:
        yaml.safe_dump({"gate": args.gate,
                        "results": [{"ranks": r["ranks"], "threads": r["threads"],
                                     "simulation_s": r["timers"]["simulation"]["max"]}
                                    for r in results]}, handle, sort_keys=False)
    return code


if __name__ == "__main__":
    sys.exit(main())
