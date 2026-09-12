"""
Combined verification figure for the 1-D Richards-equation benchmarks:
infiltration (left) and free drainage (right), each comparing frehg2 vertical
water-content profiles against a reference solution.

  (a) re-1d-infiltration : wetting front advancing downward; frehg2 vs the
      Warrick (1991) semi-analytical wetting-front positions.
  (b) re-1d-drainage     : column draining under gravity; frehg2 vs the
      laboratory/experimental profiles.

This single script builds the whole figure and is dropped identically into both
case directories, so either `re-1d-infiltration/` or `re-1d-drainage/` can
regenerate it. Both panels are read from their sibling directories via explicit
relative paths, so the two copies are interchangeable.

HDF5 layout (frehg2 groundwater output):
    /groundwater/water_content/<time_seconds>   flat NY*NX*NZ, index (j*NX+i)*NZ+k
    /groundwater/zcell/0                         static true cell-centre z, same layout
    /grid/{x_center,y_center,z_center}           grid sizes NX, NY, NZ
k = 0 is the top layer; a single vertical column is taken at (i, j) = (0, 0).
"""

from pathlib import Path

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
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
MS, MEW = 5.5, 1.0


# -----------------------------------------------------------------------------
# Per-panel configuration (read from the sibling case directories)
# -----------------------------------------------------------------------------
output_name = "output.h5"
i_profile, j_profile = 0, 0
out_pdf = Path("re1d_infiltration_drainage.pdf")

INFIL = dict(
    h5=Path("../re-1d-infiltration/out") / output_name,
    ind=[1, 2, 4],                 # 2nd/3rd/5th water-content snapshots
    ref_label="Warrick",
)
DRAIN = dict(
    h5=Path("../re-1d-drainage/out") / output_name,
    ind=[1, 4, 20, 99],            # 2nd/5th/21st/100th water-content snapshots
    ref_label="Experiment",
)

# Warrick (1991) infiltration reference: row 0 = water-content values (x), and
# each subsequent row = wetting-front depth (percent of column) for one time,
# converted to elevation z = -(1 - pct/100). Rows map to INFIL["ind"] in order.
WARRICK = np.array([
    [0.0825, 0.1650, 0.2475],
    [74.58, 75.20, 77.12],
    [61.32, 62.13, 64.54],
    [39.02, 40.04, 42.80],
])

# Free-drainage experimental profiles: columns [water content, elevation z (m)].
DRAIN_EXP = np.array([
    [0.16425893408572997, -0.39427695811338115],
    [0.17594208351217380, -1.16333166098120520],
    [0.17444522878559370, -1.91494584553626710],
    [0.18963217105749047, -2.70469608763175050],
    [0.20056584951282244, -3.42064903452543720],
    [0.19831494799621946, -0.40749017653666720],
    [0.20935061904675256, -1.15932787957338500],
    [0.21477473056186774, -1.88953739534432950],
    [0.23000590574990445, -2.69013779522488150],
    [0.24076542211528065, -3.44248409527955260],
    [0.23967180139554448, -0.43206954753778160],
    [0.25537581120438760, -1.13815941117059170],
    [0.26116779447669220, -1.89290813859370030],
    [0.27849156390251520, -2.71253455845181120],
    [0.27907726872312710, -0.42880342996067533],
    [0.29968045531271303, -1.13944635978069360],
    [0.30997544578109490, -1.89735367283997030],
    [0.31843818202447070, -2.69536168106330140],
    [0.31985985684911783, -3.43903197650730900],
])


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def time_label(seconds):
    """Compact hour/day label: hours below one day, days at or above."""
    if seconds < 86400.0:
        return f"{seconds / 3600.0:g} h"
    return f"{seconds / 86400.0:g} d"


def vertical_profile(flat_field, nx, nz, i, j):
    """One vertical column from a flattened groundwater field, index (j*NX+i)*NZ+k."""
    start = (j * nx + i) * nz
    return flat_field[start:start + nz]


def load_case(cfg):
    """Read the selected water-content profiles + z for one case, or None if the
    output file is missing (the panel then shows a 'not found' placeholder)."""
    fname = cfg["h5"]
    if not fname.exists():
        print(f"WARNING: {fname} not found — panel will be left empty.")
        return None
    with h5py.File(fname, "r") as fid:
        nx = fid["/grid/x_center"].shape[0]
        ny = fid["/grid/y_center"].shape[0]
        nz = fid["/grid/z_center"].shape[0]
        all_times = sorted(int(t) for t in fid["/groundwater/water_content"].keys())
        if max(cfg["ind"]) >= len(all_times):
            raise IndexError(f"Requested snapshot index {max(cfg['ind'])} but only "
                             f"{len(all_times)} snapshots in {fname}.")
        times = [all_times[n] for n in cfg["ind"]]
        z = vertical_profile(fid["/groundwater/zcell/0"][:], nx, nz, i_profile, j_profile)
        wc = [vertical_profile(fid[f"/groundwater/water_content/{t}"][:], nx, nz,
                               i_profile, j_profile) for t in times]
    print(f"  {fname}: {len(all_times)} snapshots, plotting "
          f"{[time_label(t) for t in times]}")
    return dict(times=times, z=z, wc=wc)


def time_colors(n):
    """n distinct colours along viridis (trimmed to avoid the pale extremes)."""
    return plt.cm.viridis(np.linspace(0.15, 0.85, n))


def draw_panel(ax, case, cfg, tag, ref_xy):
    """Plot frehg2 water-content profiles (colour = time, solid) plus the
    reference points (open circles). ref_xy is a list of (x, z) arrays, one per
    time for infiltration or a single array for drainage."""
    if case is None:
        ax.text(0.5, 0.5, "output not found", transform=ax.transAxes,
                ha="center", va="center", color="0.4")
        ax.set_xlabel(r"Water content $\theta$ [--]")
        ax.text(0.03, 0.96, tag, transform=ax.transAxes, va="top")
        return

    colors = time_colors(len(case["times"]))
    handles = []
    for k, (t, wc) in enumerate(zip(case["times"], case["wc"])):
        valid = np.isfinite(wc) & np.isfinite(case["z"])
        ax.plot(wc[valid], case["z"][valid], color=colors[k], lw=1.4, zorder=4)
        handles.append(Line2D([0], [0], color=colors[k], lw=1.4,
                              label=f"$t$ = {time_label(t)}"))

    # Reference points (open circles); per-time for infiltration, single set for
    # drainage. Colour-matched to the frehg2 time when a per-time set is given.
    if len(ref_xy) == len(case["times"]):
        for k, (rx, rz) in enumerate(ref_xy):
            ax.plot(rx, rz, ls="none", marker="o", mfc="white", mec="black",
                    mew=MEW, ms=MS, zorder=6)
    else:
        rx, rz = ref_xy[0]
        ax.plot(rx, rz, ls="none", marker="o", mfc="white", mec="black",
                mew=MEW, ms=MS, zorder=6)
    handles.append(Line2D([0], [0], ls="none", marker="o", mfc="white",
                          mec="black", mew=MEW, ms=MS, label=cfg["ref_label"]))

    ax.set_xlabel(r"Water content $\theta$ [--]")
    ax.set_ylabel("Elevation $z$ [m]")
    ax.grid(True, lw=0.3, alpha=0.3)
    ax.legend(handles=handles, loc="best", framealpha=0.92)
    ax.text(0.03, 0.96, tag, transform=ax.transAxes, va="top")


# -----------------------------------------------------------------------------
# Load both cases
# -----------------------------------------------------------------------------
print("re-1d infiltration + drainage — combined verification figure")
infil = load_case(INFIL)
drain = load_case(DRAIN)

# Reference overlays.
warrick_xy = [(WARRICK[0, :], -(1.0 - WARRICK[k + 1, :] / 100.0))
              for k in range(len(INFIL["ind"]))]      # per-time
drain_xy = [(DRAIN_EXP[:, 0], DRAIN_EXP[:, 1])]       # single set


# -----------------------------------------------------------------------------
# Figure: two panels (infiltration | drainage)
# -----------------------------------------------------------------------------
fig, (axL, axR) = plt.subplots(1, 2, figsize=(DOUBLE, 10 * cm))

draw_panel(axL, infil, INFIL, "(a)", warrick_xy)
draw_panel(axR, drain, DRAIN, "(b)", drain_xy)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nwrote {out_pdf.resolve()}")

plt.show()
