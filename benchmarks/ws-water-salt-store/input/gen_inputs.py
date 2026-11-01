#!/usr/bin/env python3
"""gen_inputs.py — inputs for the Stage-3 T9 ws gates (plan §12A).

ws1 is the S3-A12 mechanism compressed into a per-PR case: a 2-D (nx = 1)
steep beach whose dry seepage face is fed by post-allocation vents during
the ebb, a heterogeneous K window over a uniform basement, a tidal stage
with the external closure, consistent updates and density coupling, no ET.
The geometry and fields follow the 2026Frehg2 H2 case the bug report
packaged (stage3/evidence/dryout-store), scaled down:

  bed      12.5 % beach from Z_TOP at y = 0 to Z_TOE at the toe, then a
           1 % terrace to the seaward end (the stage row);
  tide     eta = A sin(2 pi t / T) [m], 60-s samples;
  wt       linear Z_WT0 (y = 0) -> the MSL shoreline, kept <= bed - 0.05
           on dry ground, MSL under the sea (the H2 initial water table);
  salt     Ghyben-Herzberg interface z_i = -h / 0.02604 under the water
           table (35 psu below, fresh above), 35 under the sea; surface
           35 where the bed is below MSL;
  soil     the H2 log-K bins (24 types, sigma ln K = 1.0 around the base
           K_g = 7.78e-5 m/s) on an iid seeded field in the top NW layers,
           the base type below (ksx = 0: the 2-D section).

Run with --variant to write one of the named configurations; the default
writes the shipped ws1 inputs.
"""

from __future__ import annotations

import argparse
import math
import pathlib

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent

KG = 7.784475161383021e-05  # H2 base (basement) K [m/s]
N_BINS = 24

VARIANTS = {
    # name: ny, dy, nz, z_top, z_toe, z_end, amp, period, z_wt0, nw, seed
    "ws1": dict(ny=64, dy=0.5, nz=16, z_top=2.0, z_toe=-0.8, z_end=-0.92,
                amp=0.5, period=7200.0, z_wt0=0.45, nw=10, seed=20261031),
}


def bins():
    """The H2 log-K bin centres: 24 equal-probability bins of N(ln KG, 1)."""
    from statistics import NormalDist
    nd = NormalDist(math.log(KG), 1.0)
    edges = [nd.inv_cdf(q) for q in np.linspace(0.0, 1.0, N_BINS + 1)[1:-1]]
    centres = []
    lo = -math.inf
    for b in range(N_BINS):
        hi = edges[b] if b < N_BINS - 1 else math.inf
        p_lo = nd.cdf(lo) if lo != -math.inf else 0.0
        p_hi = nd.cdf(hi) if hi != math.inf else 1.0
        centres.append(nd.inv_cdf(0.5 * (p_lo + p_hi)))
        lo = hi
    return [math.exp(c) for c in centres], edges


def write(path: pathlib.Path, header: str, values, fmt: str, per_line: int = 10):
    with open(path, "w", encoding="utf-8") as out:
        out.write(f"# {header}\n")
        vals = list(values)
        for k in range(0, len(vals), per_line):
            out.write(" ".join(fmt % v for v in vals[k:k + per_line]) + "\n")


def generate(name: str, outdir: pathlib.Path) -> dict:
    p = VARIANTS[name]
    ny, dy, nz = p["ny"], p["dy"], p["nz"]
    y = (np.arange(ny) + 0.5) * dy
    beach_len = (p["z_top"] - p["z_toe"]) / 0.125
    bed = np.where(y <= beach_len, p["z_top"] - 0.125 * y,
                   p["z_toe"] - 0.01 * (y - beach_len))
    bed = np.maximum(bed, p["z_end"])
    shore = (p["z_top"] - 0.0) / 0.125  # MSL shoreline
    wt = p["z_wt0"] * np.clip(1.0 - y / shore, 0.0, None)
    wt = np.where(bed > 0.0, np.minimum(wt, bed - 0.05), 0.0)
    rng = np.random.default_rng(p["seed"])
    nw = p["nw"]
    ids = np.full((ny, nz), N_BINS, dtype=int)  # base type index
    ids[:, :nw] = rng.integers(0, N_BINS, size=(ny, nw))
    # Terrain-following scaled layers: depth fraction of each cell centre
    # (uniform here, the harness stretch is mild); the GH interface needs z.
    base = min(bed) - 2.0
    frac = (np.arange(nz) + 0.5) / nz
    salt = np.zeros((ny, nz))
    for j in range(ny):
        if bed[j] <= 0.0:
            salt[j, :] = 35.0
            continue
        z_i = -max(wt[j], 0.0) / 0.02604
        zc = bed[j] - frac * (bed[j] - base)
        salt[j, :] = np.where(zc < z_i, 35.0, 0.0)
    surf = np.where(bed < 0.0, 35.0, 0.0)
    outdir.mkdir(parents=True, exist_ok=True)
    write(outdir / f"{name}_dem.dat", f"{name} bed [m]: 12.5 % beach to {p['z_toe']}, 1 % terrace", bed, "%.6f")
    write(outdir / f"{name}_wt.dat", f"{name} initial water table [m] (H2 form)", wt, "%.6f")
    write(outdir / f"{name}_soil.dat", f"{name} soil ids (j*nx+i)*nz+k: iid log-K window in layers 1..{nw}, base below; seed {p['seed']}",
          ids.ravel(), "%d")
    write(outdir / f"{name}_salt_gw.dat", f"{name} subsurface salinity IC [psu], (j*nx+i)*nz+k: Ghyben-Herzberg", salt.ravel(), "%.4f")
    write(outdir / f"{name}_salt_surf.dat", f"{name} surface salinity IC [psu]: 35 where the bed is below MSL", surf, "%.1f")
    t = np.arange(0.0, 4.0 * p["period"] + 60.0, 60.0)
    with open(outdir / f"{name}_tide.dat", "w", encoding="utf-8") as out:
        out.write(f"# eta(t) = {p['amp']} sin(2 pi t / {p['period']}) [m], 60-s samples, 4 cycles\n")
        for tt in t:
            out.write(f"{tt:.0f} {p['amp'] * math.sin(2.0 * math.pi * tt / p['period']):.10f}\n")
    return {"ny": ny, "dy": dy, "nz": nz, "bed": bed, "seaward_y": ny * dy}


def soil_types_yaml():
    centres, _ = bins()
    rows = []
    for b, k in enumerate(centres):
        rows.append({"name": f"lnk{b:02d}", "ksx": 0.0, "ksy": k, "ksz": k,
                     "theta_s": 0.41, "theta_r": 0.15, "vg_alpha": 2.0, "vg_n": 1.41})
    rows.append({"name": "base", "ksx": 0.0, "ksy": KG, "ksz": KG,
                 "theta_s": 0.41, "theta_r": 0.15, "vg_alpha": 2.0, "vg_n": 1.41})
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", default="ws1")
    ap.add_argument("--out", type=pathlib.Path, default=HERE)
    a = ap.parse_args()
    info = generate(a.variant, a.out)
    print(f"{a.variant}: ny={info['ny']} nz={info['nz']} bed {info['bed'].max():.3f}..{info['bed'].min():.3f}")
