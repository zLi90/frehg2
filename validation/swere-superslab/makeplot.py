"""
Visualize and verify the swere-superslab coupled surface-subsurface benchmark
(SERGHEI-SWE-RE case5): a heterogeneous multi-soil vertical slab where 3 h of
rain infiltrates a layered profile and a low-conductivity block forces return
flow to the surface outlet over a 12 h event.

This script builds BOTH deliverables the case is judged on (one `makeplot.py`,
two PDFs), against the four inter-comparison codes shipped in `data/`
(ParFlow, CATHY, HGS, Cast3M):

  1. superslab_ponding_outflow.pdf   -- (a) surface ponding storage and
     (b) outlet flow rate vs time.
  2. superslab_saturation_profiles.pdf -- vertical saturation profiles at
     t = 3 h and 6 h, at x = 0, 8, 40 m (2 rows x 3 columns).

frehg2 signals
--------------
* Ponding storage [m3]  = `/monitor/mass_audit` `volume` column (surface water
  volume) vs time.
* Outlet flow rate [m3/h] = d/dt of the CUMULATIVE `boundary_outflow` column
  x 3600 (the exact BC flux the solver applied; the mass_audit is logged every
  dt so its derivative is clean -- same convention as the other SWE cases).
* Saturation [-] = `/groundwater/water_content` / theta_s (all three soils share
  theta_s = 0.1, per the README; max(water_content) ~= 0.1 confirms it).

HDF5 layout (frehg2, ny = 1 so the flat 3-D index (j*nx+i)*nz+k -> i*nz+k)
--------------------------------------------------------------------------
    /monitor/mass_audit   (nt, 8) cols time,volume,rain,evaporation,
                          boundary_outflow,bc_inflow,seepage,clamped
    /surface/{depth,eta,seepage}/<t>   (nx,) surface fields (unused here; the
                          monitor tables carry the integrated signals)
    /groundwater/water_content/<t>     (nx*nz,) index i*nz+k, reshape (nx,nz);
                          k = 0 is the TOP layer just below the surface.
    /groundwater/zcell/0  (nx*nz,) TRUE cell-centre z (terrain-following:
                          follow_terrain=true, so /grid/z_center is NOT usable);
                          depth below the local surface = zcell[i] - bottom[i].
    /grid/{x_center,bottom,dz}         (nx,) horizontal coord, DEM, layer height.

Reference data (`data/<code>-...csv`, comma-separated, no header)
----------------------------------------------------------------
    <code>-slab-ponding.csv   col0 = time [h], col1 = ponding storage [m3]
    <code>-slab-outflow.csv   col0 = time [h], col1 = outlet flow rate [m3/h]
    <code>-t{3,6}-x{0,8,40}.csv  col0 = saturation [-], col1 = depth below
                          surface [m] (positive down -> plotted at z = -depth)
(ParFlow only provides t6 at x = 0; the loader skips absent files.)

Reference: Maina et al. / SERGHEI-SWE-RE supplement, superslab (case5)
inter-comparison (ParFlow, CATHY, HGS, Cast3M).
"""

from pathlib import Path

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D


# -----------------------------------------------------------------------------
# Unified journal style (Elsevier): fonts 9-10 pt, editable-text PDF, palette.
# -----------------------------------------------------------------------------
cm = 1 / 2.54
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.linewidth": 0.8, "lines.linewidth": 1.4,
    "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
})
SINGLE, DOUBLE = 8.8 * cm, 17.0 * cm
C_SIM = "#1f5fa6"                                   # frehg2: solid blue line
REF_COLORS = ["#d1611f", "#2e8b57", "#7b3fa0", "#ff7f0e"]
MS, MEW = 5.0, 1.1


# -----------------------------------------------------------------------------
# Configuration (from swere-superslab.yaml / README)
# -----------------------------------------------------------------------------
output_file = Path("out/output.h5")
data_dir = Path("data")
theta_s = 0.1                       # saturated water content (all 3 soils)
x_targets = [0.0, 8.0, 39.0]        # profile locations [m] (outlet at x = 0)
t_targets_h = [3.0, 6.0]            # profile times [h]
out_pond = Path("superslab_ponding_outflow.pdf")
out_sat = Path("superslab_saturation_profiles.pdf")

# Inter-comparison codes: file stem, legend label, palette colour, marker.
REF_MODELS = [
    dict(key="PF",     label="ParFlow", color=REF_COLORS[0], marker="o"),
    dict(key="CATHY",  label="CATHY",   color=REF_COLORS[1], marker="^"),
    dict(key="HGS",    label="HGS",     color=REF_COLORS[2], marker="s"),
    dict(key="cast3m", label="Cast3M",  color=REF_COLORS[3], marker="D"),
]


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def load_csv(path):
    """Two-column comma CSV (no header) -> (col0, col1) arrays, or None if the
    file is absent (several reference series are only partially available)."""
    if not path.exists():
        return None
    d = np.atleast_2d(np.loadtxt(path, delimiter=","))
    return d[:, 0], d[:, 1]


def open_marker(ax, x, y, model, **kw):
    """Plot a reference series as open (hollow) colour-coded markers."""
    ax.plot(x, y, ls="none", marker=model["marker"], mfc="none",
            mec=model["color"], mew=MEW, ms=MS, **kw)


def ref_legend_handles():
    """Proxy handles: frehg2 solid line + one open marker per reference code."""
    handles = [Line2D([0], [0], color=C_SIM, lw=1.8, label="frehg2")]
    for m in REF_MODELS:
        handles.append(Line2D([0], [0], ls="none", marker=m["marker"], mfc="none",
                              mec=m["color"], mew=MEW, ms=MS, label=m["label"]))
    return handles


def containing_cell(coord, target):
    """Index of the finite-volume cell that CONTAINS `target` (half-open
    [edge_i, edge_{i+1})). This is used instead of nearest-centre because the
    superslab has material-block interfaces landing exactly on the sample
    locations (e.g. x = 40 m is a block boundary); a nearest-centre tie-break
    would silently sample the adjacent block. Cell centres are at (i+0.5)*dx,
    so the left edges are coord - dx/2."""
    dx = float(np.mean(np.diff(coord)))
    edges = coord - dx / 2.0
    return int(np.clip(np.searchsorted(edges, target, side="right") - 1, 0, coord.size - 1))


# -----------------------------------------------------------------------------
# Load frehg2 output
# -----------------------------------------------------------------------------
with h5py.File(output_file, "r") as f:
    x = f["/grid/x_center"][:]
    bottom = f["/grid/bottom"][:]
    dz = f["/grid/dz"][:]
    nx, nz = x.size, dz.size

    ma = f["/monitor/mass_audit"][:]           # time,volume,rain,evap,bout,...,seepage,clamped
    gma = f["/monitor/gw_mass_audit"][:]       # time,volume,...
    t_s = ma[:, 0]
    pond_vol = ma[:, 1]                         # surface water volume [m3]
    rain_cum = ma[:, 2]                         # cumulative rain [m3]
    bout_cum = ma[:, 4]                         # cumulative boundary outflow [m3]
    seep_cum = ma[:, 6]                         # cumulative surface<->gw seepage [m3]

    # Snapshot times + the terrain-following z of every cell centre.
    snap_t = sorted(int(t) for t in f["/groundwater/water_content"].keys())
    zcell = f["/groundwater/zcell/0"][:].reshape(nx, nz)      # i*nz+k

    # Vertical saturation profiles at the requested (time, x) pairs.
    profiles = {}
    for th in t_targets_h:
        t_snap = min(snap_t, key=lambda tt: abs(tt - th * 3600.0))
        wc = f[f"/groundwater/water_content/{t_snap}"][:].reshape(nx, nz)
        for xt in x_targets:
            i = containing_cell(x, xt)
            sat = wc[i] / theta_s                              # k = 0 top .. nz-1 bottom
            z_rel = zcell[i] - bottom[i]                       # depth below surface (<= 0)
            profiles[(th, xt)] = dict(sat=sat, z=z_rel, i=i, t=t_snap)

# Outlet flow rate [m3/h] = d/dt of cumulative outflow (dt is uniform & small).
q_flow = np.gradient(bout_cum, t_s) * 3600.0
t_h = t_s / 3600.0


# -----------------------------------------------------------------------------
# Console verification report + mass-balance sanity check
# -----------------------------------------------------------------------------
print("swere-superslab — coupled surface-subsurface superslab benchmark")
print(f"  grid                 : nx={nx}, nz={nz}, dz={dz[0]:.3f} m; DEM {bottom.min():.1f}-{bottom.max():.1f} m")
print(f"  run length           : {t_s[-1]/3600:.1f} h ({t_s.size} mass_audit rows, dt={t_s[1]-t_s[0]:.2f} s)")
print(f"  cumulative rain       : {rain_cum[-1]:.4f} m3")
print(f"  peak ponding storage  : {pond_vol.max():.4e} m3")
print(f"  peak outlet flow rate : {q_flow.max():.4e} m3/h")
print(f"  surface seepage (cum) : {seep_cum[-1]:.4e} m3")
d_gw = gma[-1, 1] - gma[0, 1]
print(f"  d(gw storage)         : {d_gw:+.4e} m3")

# Surface balance: rain_in should equal d(ponding) + outflow + seepage + evap.
residual = rain_cum[-1] - (pond_vol[-1] - pond_vol[0]) - bout_cum[-1] - seep_cum[-1] - ma[-1, 3]
degenerate = (pond_vol.max() < 1e-9) and (q_flow.max() < 1e-9) and (rain_cum[-1] > 1e-6)
if degenerate:
    print("\n  *** WARNING: ponding, outflow, and seepage are all ~0 while "
          f"{rain_cum[-1]:.2f} m3 of rain was injected and gw storage is static "
          f"(d={d_gw:+.2e} m3).")
    print(f"      The surface balance leaves {residual:.2f} m3 unaccounted — this run "
          "did NOT route rainfall through the surface/coupler (mass-balance issue).")
    print("      The frehg2 ponding/outflow curves below are therefore degenerate "
          "(flat 0); regenerate out/output.h5 from a valid coupled run.")
else:
    print(f"  surface mass residual : {residual:+.4e} m3 (rain - dV_surf - outflow - seepage - evap)")

# Saturation-profile sampling (containing cell for each nominal x location).
print("  saturation profiles   :")
for th in t_targets_h:
    cells = ", ".join(f"x={int(xt)}m->i={profiles[(th, xt)]['i']}"
                      f"(xc={x[profiles[(th, xt)]['i']]:.1f}m)" for xt in x_targets)
    print(f"    t={int(th)}h (snap {profiles[(th, x_targets[0])]['t']}s): {cells}")


# =============================================================================
# Figure 1: outlet flow rate vs time
# =============================================================================
fig1, axQ = plt.subplots(1, 1, figsize=(SINGLE, 7 * cm))

# Outlet flow rate -----------------------------------------------------------
axQ.plot(t_h, q_flow, color=C_SIM, lw=1.6, zorder=5)
for m in REF_MODELS:
    ref = load_csv(data_dir / f"{m['key']}-slab-outflow.csv")
    if ref is not None:
        open_marker(axQ, ref[0], ref[1], m, zorder=4)
axQ.set_xlabel("Time [h]")
axQ.set_ylabel(r"Outlet flow rate [m$^{3}$ h$^{-1}$]")
axQ.set_xlim(0, t_h[-1])
axQ.grid(True, lw=0.3, alpha=0.3)

fig1.legend(handles=ref_legend_handles(), loc="lower center", ncol=5,
            frameon=False, bbox_to_anchor=(0.5, -0.02))
fig1.tight_layout(rect=(0, 0.08, 1, 1))
fig1.savefig(out_pond, bbox_inches="tight")
print(f"\nwrote {out_pond.resolve()}")


# =============================================================================
# Figure 2: vertical saturation profiles (rows = time, cols = x location)
# =============================================================================
fig2, axes = plt.subplots(len(t_targets_h), len(x_targets),
                          figsize=(DOUBLE, 11 * cm), sharex=True, sharey=True)
tags = iter("(a) (b) (c) (d) (e) (f) (g) (h) (i)".split())

for r, th in enumerate(t_targets_h):
    for c, xt in enumerate(x_targets):
        ax = axes[r, c]
        pr = profiles[(th, xt)]
        # frehg2 profile (solid blue).
        ax.plot(pr["sat"], pr["z"], color=C_SIM, lw=1.6, zorder=5)
        # Reference codes (open markers); saturation vs z = -depth.
        for m in REF_MODELS:
            ref = load_csv(data_dir / f"{m['key']}-t{int(th)}-x{int(xt)}.csv")
            if ref is not None:
                open_marker(ax, ref[0], -ref[1], m, zorder=4)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(-5.0, 0.0)
        ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.grid(True, lw=0.3, alpha=0.3)
        ax.text(0.05, 0.05, f"{next(tags)} $t$={int(th)} h, $x$={int(xt)} m",
                transform=ax.transAxes, va="bottom", ha="left",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))
        if r == len(t_targets_h) - 1:
            ax.set_xlabel("Saturation [--]")
        if c == 0:
            ax.set_ylabel("Depth below surface [m]")

fig2.legend(handles=ref_legend_handles(), loc="lower center", ncol=5,
            frameon=False, bbox_to_anchor=(0.5, -0.02))
fig2.tight_layout(rect=(0, 0.05, 1, 1))
fig2.savefig(out_sat, bbox_inches="tight")
print(f"wrote {out_sat.resolve()}")

plt.show()
