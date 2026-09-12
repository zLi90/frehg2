"""
Visualize and verify the swere-lateral-hillslope coupled surface-subsurface
benchmark (SERGHEI-SWE-RE case4): a 5-year lateral hillslope where saturated
subsurface drainage between two fixed-head boundaries, fed by light rainfall
recharge, sets the water table across a gently sloping 4 km strip.

It reproduces the two panels of the SERGHEI reference script `plot_wtdz.py`
in one figure (one `makeplot.py`, one PDF), against the Hydrus-2D reference:

  lateral_hillslope_watertable.pdf
    (a) water table at the domain centre (i = nx//2) vs time, and
    (b) water-table profile vs distance at the final snapshot.

frehg2 signal
-------------
* Water table [m] = elevation where the pressure head psi crosses zero, reported
  as a HEIGHT ABOVE THE AQUIFER BASE (0 at the slab bottom, 15 m at the surface)
  to match the SERGHEI reference datum. The vertical profile is scanned from the
  top (k = 0) downward for the first psi < 0 -> psi >= 0 crossing and linearly
  interpolated (identical to plot_wtdz.py). A column that is saturated to the
  surface reports the top cell centre; a fully unsaturated column reports NaN.

  NOTE: `/groundwater/hydraulic_head` is a MISLABELLED dataset -- it stores the
  pressure head psi (SERGHEI's `hd`), not the true head H = psi + z. This is
  confirmed here: the t = 0 field gives WT = 7.0 m at i = 0 and 0.9 m at the
  right end, exactly the prescribed left/right BC water tables (eta = -4.1 /
  -14.1 m -> 7.0 / 0.9 m above base).

HDF5 layout (frehg2, ny = 1 so the flat 3-D index (j*nx+i)*nz+k -> i*nz+k)
--------------------------------------------------------------------------
    /groundwater/hydraulic_head/<t>  (nx*nz,) pressure head psi, i*nz+k, reshape
                          (nx, nz); k = 0 is the TOP layer just below the surface.
    /groundwater/zcell/0  (nx*nz,) TRUE cell-centre z (terrain-following:
                          follow_terrain=true, so /grid/z_center is NOT usable).
    /grid/{x_center,bottom,dz}        (nx,) horizontal coord, DEM, layer height.

Reference data (comma-separated, in the case directory)
-------------------------------------------------------
    hydrus-time.csv  col0 = time [days], col1 = water table [m] at the centre.
                     (No header; a spurious negative-time row -- the digitised
                     t ~ 0 point -- is dropped.)
    HYDRUS2D.csv     col0 = distance x [m], col1 = water table [m] at 5 yr.
                     (Header `x, y`.)

Reference: Beegum et al. / SERGHEI-SWE-RE supplement, lateral hillslope (case4),
Hydrus-2D solution.
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
# Configuration (from swere-lateral-hillslope.yaml / README)
# -----------------------------------------------------------------------------
output_file = Path("out/output.h5")
head_dset = "/groundwater/hydraulic_head"          # stores pressure head psi
T_TARGET_DAYS = 157680000.0 / 86400.0              # 5 yr SERGHEI simLength (1825 d)
ic_frac = 0.5                                       # sample column = nx//2 (centre)

ref_time_file = Path("hydrus-time.csv")            # [days, WT] at the centre
ref_prof_file = Path("HYDRUS2D.csv")               # [x_m, WT] profile at 5 yr
REF_LABEL = "Hydrus-2D"
REF_C, REF_M = REF_COLORS[0], "^"

out_pdf = Path("lateral_hillslope_watertable.pdf")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def load_ref(path):
    """Two-column comma CSV -> (N, 2) array with any header / blank rows dropped
    (genfromtxt turns a text header into NaN, which the finite-mask removes), or
    None if the file is absent."""
    if not path.exists():
        return None
    d = np.atleast_2d(np.genfromtxt(path, delimiter=","))
    return d[np.all(np.isfinite(d), axis=1)]


def water_table(psi_col, z_col, base):
    """Water-table height above `base`, from the top-down first psi < 0 -> >= 0
    crossing (psi < 0 unsaturated above the table, psi >= 0 saturated below).
    Matches plot_wtdz.py's interpolation. Returns the top cell centre if the
    column is saturated to the surface, or NaN if it is fully unsaturated."""
    if psi_col[0] >= 0.0:                           # saturated to the surface
        return z_col[0] - base
    for k in range(psi_col.size - 1):
        h1, h2 = psi_col[k], psi_col[k + 1]
        if h1 < 0.0 and h2 >= 0.0:
            z1, z2 = z_col[k], z_col[k + 1]
            return (z1 - h1 * (z1 - z2) / (h1 - h2)) - base
    return np.nan                                   # no table within the slab


def plot_ref(ax, xv, yv):
    """Reference series as open (hollow) colour-coded markers."""
    ax.plot(xv, yv, ls="none", marker=REF_M, mfc="none", mec=REF_C,
            mew=MEW, ms=MS, zorder=4)


def legend_handles():
    """Proxy handles: frehg2 solid line + open reference marker."""
    return [
        Line2D([0], [0], color=C_SIM, lw=1.8, label="frehg2"),
        Line2D([0], [0], ls="none", marker=REF_M, mfc="none", mec=REF_C,
               mew=MEW, ms=MS, label=REF_LABEL),
    ]


# -----------------------------------------------------------------------------
# Load frehg2 output
# -----------------------------------------------------------------------------
with h5py.File(output_file, "r") as f:
    x = f["/grid/x_center"][:]
    bottom = f["/grid/bottom"][:]
    dz = f["/grid/dz"][:]
    nx, nz = x.size, dz.size
    dx = float(np.mean(np.diff(x)))
    H_slab = float(dz.sum())                         # 60 x 0.25 = 15.0 m
    base = bottom - H_slab                           # aquifer base per column

    zcell = f["/groundwater/zcell/0"][:].reshape(nx, nz)     # i*nz+k
    ic = int(round(ic_frac * (nx - 1)))              # centre column (i = 20)

    snaps = sorted(int(t) for t in f[head_dset].keys())
    t_days = np.array(snaps, dtype=float) / 86400.0

    # (a) water table at the centre, every snapshot.
    wt_centre = np.empty(len(snaps))
    for n, t in enumerate(snaps):
        psi = f[f"{head_dset}/{t}"][:].reshape(nx, nz)
        wt_centre[n] = water_table(psi[ic], zcell[ic], base[ic])

    # (b) water-table profile across x at the final snapshot.
    t_final = snaps[-1]
    psi_final = f[f"{head_dset}/{t_final}"][:].reshape(nx, nz)
    wt_prof = np.array([water_table(psi_final[i], zcell[i], base[i])
                        for i in range(nx)])

t_final_days = t_final / 86400.0

# Reference series.
ref_t = load_ref(ref_time_file)
if ref_t is not None:
    ref_t = ref_t[ref_t[:, 0] >= 0.0]               # drop digitised t ~ 0 outlier
    ref_t = ref_t[np.argsort(ref_t[:, 0])]
ref_p = load_ref(ref_prof_file)
if ref_p is not None:
    ref_p = ref_p[np.argsort(ref_p[:, 0])]


# -----------------------------------------------------------------------------
# Console verification report
# -----------------------------------------------------------------------------
print("swere-lateral-hillslope — coupled lateral hillslope subsurface flow")
print(f"  grid                 : nx={nx}, nz={nz}, dz={dz[0]:.3f} m; DEM "
      f"{bottom.min():.1f}-{bottom.max():.1f} m; slab H={H_slab:.2f} m")
print(f"  centre column        : i={ic} (x={x[ic]:.1f} m)")
print(f"  water table datum    : height above the aquifer base "
      f"(0 = slab bottom, {H_slab:.2f} m = surface)")
frac = t_final_days / T_TARGET_DAYS
print(f"  run coverage         : {len(snaps)} snapshots, t=0..{t_final_days:.1f} d "
      f"({t_final_days / 365.0:.2f} yr) of {T_TARGET_DAYS:.0f} d target "
      f"({100.0 * frac:.0f}%)")
if frac < 0.95:
    print("  *** NOTE: this is a PARTIAL run. Panel (b) compares the frehg2 field "
          f"at t={t_final_days:.0f} d against the Hydrus-2D 5-yr profile, so the")
    print("      two are at different times; the panel (a) overlap is the like-for-"
          "like check until the 5-yr run completes.")

print(f"  centre WT            : initial {wt_centre[0]:.3f} m -> "
      f"final {wt_centre[-1]:.3f} m (t={t_final_days:.0f} d)")
if ref_t is not None and ref_t.size:
    mask = ref_t[:, 0] <= t_final_days + 1e-9
    if mask.any():
        f_at_ref = np.interp(ref_t[mask, 0], t_days, wt_centre)
        mae = float(np.mean(np.abs(f_at_ref - ref_t[mask, 1])))
        print(f"  centre WT vs ref     : MAE={mae * 100:.1f} cm over "
              f"{int(mask.sum())} Hydrus points within the run window "
              f"(<= {t_final_days:.0f} d)")

print(f"  final profile WT     : min {np.nanmin(wt_prof):.2f} m "
      f"@i{int(np.nanargmin(wt_prof))} / max {np.nanmax(wt_prof):.2f} m "
      f"@i{int(np.nanargmax(wt_prof))}")
sat_cols = [i for i in (0, nx - 1) if wt_prof[i] >= H_slab - dz[0]]
if sat_cols:
    print(f"  *** boundary columns saturated to the surface: i={sat_cols} "
          f"(WT ~ {H_slab - dz[0]:.2f} m); the fixed-head BCs give the reference "
          "endpoints ~7.0 / ~1.2 m instead.")


# =============================================================================
# Figure: (a) water table at the centre vs time; (b) profile vs distance
# =============================================================================
fig, (axT, axB) = plt.subplots(2, 1, figsize=(SINGLE, 12 * cm))

# (a) Water table at the domain centre vs time ------------------------------
axT.plot(t_days / 365.0, wt_centre, color=C_SIM, lw=1.6, zorder=5)
if ref_t is not None and ref_t.size:
    plot_ref(axT, ref_t[:, 0] / 365.0, ref_t[:, 1])
axT.set_xlabel("Time [years]")
axT.set_ylabel("Water table [m]")
axT.set_xlim(0.0, T_TARGET_DAYS / 365.0)
axT.set_xticks(np.arange(0.0, T_TARGET_DAYS / 365.0 + 0.5, 1.0))
axT.grid(True, lw=0.3, alpha=0.3)
axT.text(0.03, 0.95, "(a)", transform=axT.transAxes, va="top", ha="left",
         path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

# (b) Water-table profile vs distance at the final snapshot -----------------
axB.plot(x, wt_prof, color=C_SIM, lw=1.6, zorder=5)
if ref_p is not None and ref_p.size:
    plot_ref(axB, ref_p[:, 0], ref_p[:, 1])
axB.set_xlabel("Distance $x$ [m]")
axB.set_ylabel("Water table [m]")
axB.set_xlim(0.0, nx * dx)
axB.grid(True, lw=0.3, alpha=0.3)
axB.text(0.03, 0.95, "(b)", transform=axB.transAxes, va="top", ha="left",
         path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])

fig.legend(handles=legend_handles(), loc="lower center", ncol=2,
           frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nwrote {out_pdf.resolve()}")

plt.show()
