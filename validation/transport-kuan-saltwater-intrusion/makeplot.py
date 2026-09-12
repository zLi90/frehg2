"""
Visualize the Kuan et al. (2012) beach-aquifer case: pore-water salinity, the
surface water level, and the beach topography. The figure has two stacked
panels sharing the cross-section axes:

    (a) ss  -- no-tide steady salt wedge                    (out-ss/)
    (b) td  -- tidal, plotted as the tidal (phase) MEAN     (out-td/)

Steady (ss) vs tidal (td)
-------------------------
The `ss` variant holds the sea at a fixed stage, so the field is (nearly) steady
and a single final snapshot is representative. The `td` variant forces the
seaward stage with the lab tide (period T = 62 s, mean 0.14 m, range
0.097-0.183 m); the salinity field then oscillates strongly with the tide.
Because the field is written every 600 s (~9.7 tidal cycles apart), a single
snapshot lands at an arbitrary tidal phase and is NOT representative. Instead
each 600 s snapshot samples a *different* tidal phase, so averaging the last
`n_avg` snapshots yields the tidal (phase) MEAN field -- the physically
meaningful steady picture of a tidally forced beach. That mean reveals the three
classic Kuan features the steady case lacks: a lower salt wedge, an upper saline
plume (USP) near the beach face, and a freshwater discharge tube between them.

Reference: Kuan, W.K., Jin, G., Xin, P., Robinson, C., Gibbes, B., Li, L.
(2012), "Tidal influence on seawater intrusion in unconfined coastal aquifers",
Water Resour. Res. 48, W02502, doi:10.1029/2011WR010678.

Geometry / HDF5 layout (frehg2 subsurface output, terrain-following)
--------------------------------------------------------------------
    /grid/bottom      (ny,)  beach-surface (bed) elevation = TOP of the aquifer
    /grid/y_center    (ny,)  landward (y=0) -> seaward cell-centre y [m]
    /groundwater/zcell/0     (ny*nz,) TRUE cell-centre z per cell (terrain
                             following: the nominal /grid/z_center is NOT used
                             because each column's 18 layers are stretched
                             between its bed and the ~flat aquifer base).
    /transport/concentration/<t>          (ny*nz,) index j*nz + k, reshape
                                           (ny, nz); k=0 just under the bed.
    /surface/{eta,depth}/<t>              (ny,) free-surface elevation & depth.
"""

from pathlib import Path

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


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
C_IFACE = "#d1201f"
C_BED_EDGE = "#6b4f2a"
C_WATER, C_WATER_EDGE = "#7fb6e6", "#12508f"


# -----------------------------------------------------------------------------
# Per-variant configuration (from kuan-ss.yaml / kuan-td.yaml / README)
# -----------------------------------------------------------------------------
s_sea = 35.0                       # seawater salinity [psu] (= density surrogate)
wet_tol = 1.0e-4                   # [m] surface depth above which a cell is "wet"
x_crop = 1.5                       # landward x-limit [m] (long/narrow domain)
out_pdf = Path("kuan_saltwedge_field.pdf")

CASES = {
    "ss": dict(
        h5=Path("out-ss/output.h5"),
        label="no-tide, steady",
        tidal=False,
        exp_csv=Path("plot-data-ss.csv"),
    ),
    "td": dict(
        h5=Path("out-td/output.h5"),
        label="tidal",
        tidal=True,
        tide_file=Path("input/tide.dat"),
        tide_period=62.0,           # tidal period [s] (from tide.dat)
        n_avg=8,                    # snapshots in the trailing tidal-mean window
        exp_csv=Path("plot-data-td.csv"),
    ),
}


# -----------------------------------------------------------------------------
# Tidal forcing statistics (mean sea level, high/low water) from tide.dat.
# -----------------------------------------------------------------------------
def tide_stats(path, t_lo, t_hi):
    """Mean / high / low stage over one representative tidal cycle. The forcing
    is periodic, so any full cycle gives the true envelope; fall back to the
    whole series if the requested window is unavailable."""
    d = np.loadtxt(path)
    t, e = d[:, 0], d[:, 1]
    sel = e[(t >= t_lo) & (t <= t_hi)]
    if sel.size < 2:
        sel = e
    return float(sel.mean()), float(sel.max()), float(sel.min())


# -----------------------------------------------------------------------------
# Load grid + the field to plot (single snapshot for ss, tidal-mean for td).
# -----------------------------------------------------------------------------
def load(cfg, n_avg):
    with h5py.File(cfg["h5"], "r") as f:
        y = f["/grid/y_center"][:]                     # (ny,) landward -> seaward
        bed = f["/grid/bottom"][:]                     # (ny,) beach surface = aquifer top
        ny = y.size
        nz = f["/grid/z_center"].shape[0]
        zcell = f["/groundwater/zcell/0"][:].reshape(ny, nz)   # (ny, nz) TRUE z

        times = sorted(int(t) for t in f["/transport/concentration"].keys())
        info = dict(y=y, bed=bed, ny=ny, nz=nz, zcell=zcell, times=times)

        if not cfg["tidal"]:
            # Steady: the final snapshot, plus the previous one for a drift check.
            t_end, t_prev = times[-1], times[-2]
            C = f[f"/transport/concentration/{t_end}"][:].reshape(ny, nz)
            C_prev = f[f"/transport/concentration/{t_prev}"][:].reshape(ny, nz)
            eta = f[f"/surface/eta/{t_end}"][:]
            depth = f[f"/surface/depth/{t_end}"][:]
            info.update(C=C, t_end=t_end, drift=np.max(np.abs(C - C_prev)),
                        drift_dt=t_end - t_prev, eta=eta, depth=depth)
        else:
            # Tidal: mean over the last n_avg snapshots (each a different phase).
            win = times[-n_avg:]
            stack = np.stack([f[f"/transport/concentration/{t}"][:].reshape(ny, nz)
                              for t in win])
            C = stack.mean(0)
            prev = times[-2 * n_avg:-n_avg]
            if len(prev) == n_avg:
                pstack = np.stack([f[f"/transport/concentration/{t}"][:].reshape(ny, nz)
                                   for t in prev])
                secular = np.max(np.abs(C - pstack.mean(0)))
            else:
                secular = np.nan
            info.update(C=C, win=win, t_end=times[-1], secular=secular)
    return info


# -----------------------------------------------------------------------------
# Build the terrain-following raster geometry: transpose to (nz, ny) and augment
# with a bed row (top) and base row (bottom) so the field fills the aquifer with
# no half-cell gap at the sloping surface.
# -----------------------------------------------------------------------------
def build_raster(info):
    Zc = info["zcell"].T                               # (nz, ny) cell-centre z
    Cc = info["C"].T                                   # (nz, ny) salinity
    base = Zc[-1, :] - 0.5 * (Zc[-2, :] - Zc[-1, :])   # bottom face (extrapolated)
    bed = info["bed"]

    Z_aug = np.vstack([bed[None, :], Zc, base[None, :]])
    C_aug = np.vstack([Cc[0:1, :], Cc, Cc[-1:, :]])
    Y_aug = np.broadcast_to(info["y"][None, :], Z_aug.shape)
    return dict(Zc=Zc, Cc=Cc, base=base, Z_aug=Z_aug, C_aug=C_aug, Y_aug=Y_aug)


def draw_interface(ax, Yc, Zc, Cc, level, lw=2.2, halo=3.4, zorder=5):
    cs = ax.contour(Yc, Zc, Cc, levels=[level], colors=C_IFACE,
                    linewidths=lw, zorder=zorder)
    if hasattr(cs, "set_path_effects"):
        cs.set_path_effects([pe.withStroke(linewidth=halo, foreground="white")])
    return cs


def load_experiment(path):
    """Experimental 18-psu interface points: CSV columns (x = seaward distance,
    y = elevation) -> return (dist, elev). May contain >1 branch."""
    if not path.exists():
        return None
    d = np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
    return d[:, 0], d[:, 1]


def contour_segments(cs):
    """Vertices of a ContourSet as a list of (Ni, 2) polyline arrays."""
    if hasattr(cs, "get_paths"):                       # Matplotlib >= 3.8
        return [p.vertices for p in cs.get_paths() if len(p.vertices) > 1]
    segs = []                                          # pragma: no cover (older mpl)
    for coll in cs.collections:
        segs += [seg for seg in coll.get_segments() if len(seg) > 1]
    return segs


def point_to_polyline(px, pz, segs):
    """Shortest Euclidean distance from point (px, pz) to any contour segment."""
    best = np.inf
    for v in segs:
        a, b = v[:-1], v[1:]                           # segment endpoints
        ab = b - a
        ap = np.column_stack([px - a[:, 0], pz - a[:, 1]])
        denom = (ab ** 2).sum(1)
        t = np.clip((ap * ab).sum(1) / np.where(denom > 0, denom, 1.0), 0.0, 1.0)
        proj = a + t[:, None] * ab
        best = min(best, np.hypot(px - proj[:, 0], pz - proj[:, 1]).min())
    return best


def interface_misfit(exp, cs):
    """Per-point distance from each experimental interface point to the simulated
    18-psu contour, and their mean (the interface MAE used by the b6 gate)."""
    segs = contour_segments(cs)
    if exp is None or not segs:
        return None, None
    xe, ze = exp
    dists = np.array([point_to_polyline(xe[i], ze[i], segs) for i in range(xe.size)])
    return dists, float(dists.mean())


# -----------------------------------------------------------------------------
# Draw one panel (ss or td) and return the pcolormesh handle for a shared bar.
# -----------------------------------------------------------------------------
def draw_panel(ax, key, tag, show_xlabel):
    cfg = CASES[key]
    tidal = cfg["tidal"]
    n_avg = cfg.get("n_avg", 8)
    info = load(cfg, n_avg)
    y, bed = info["y"], info["bed"]
    dy = float(np.round(np.diff(y).mean(), 6))
    C = info["C"]
    rast = build_raster(info)
    Yc = np.broadcast_to(y[None, :], rast["Zc"].shape)     # (nz, ny) for contouring
    half = 0.5 * s_sea                                      # 17.5 psu salt interface
    exp = load_experiment(cfg["exp_csv"])

    # 0.5-interface toe along the base (landward-most reach of 17.5 psu bottom row).
    base_sal = C[:, -1]
    toe_y = np.interp(half, base_sal, y) if base_sal.max() >= half else np.nan

    # ---- console verification report ----------------------------------------
    print(f"\nkuan-saltwater-intrusion ({cfg['label']}) — beach-aquifer salinity")
    if not tidal:
        print(f"  final time              : t = {info['t_end']} s (single snapshot)")
        print(f"  steadiness              : max|C(t) - C(t-{info['drift_dt']}s)| = "
              f"{info['drift']:.3f} psu")
    else:
        win = info["win"]
        T = cfg["tide_period"]
        span = win[-1] - win[0]
        print(f"  tidal-MEAN field        : average of {len(win)} snapshots "
              f"t = {win[0]}..{win[-1]} s (last {span} s ~= {span/T:.0f} cycles)")
        print(f"    secular drift         : max|mean(win) - mean(prev win)| = "
              f"{info['secular']:.3f} psu")
    print(f"  salinity range          : [{C.min():.2f}, {C.max():.2f}] psu (s_sea = {s_sea})")
    print(f"  {half:.0f}-psu interface toe on base: y = {toe_y:.3f} m "
          f"({y[-1]+dy/2 - toe_y:.3f} m from the sea)")

    # ---- field + interface --------------------------------------------------
    pcm = ax.pcolormesh(rast["Y_aug"], rast["Z_aug"], rast["C_aug"],
                        shading="gouraud", cmap="viridis", vmin=0.0, vmax=s_sea,
                        rasterized=True, zorder=2)
    draw_interface(ax, Yc, rast["Zc"], rast["Cc"], half)
    handles = [Line2D([0], [0], color=C_IFACE, lw=2.2, label=f"{half:.0f} psu interface")]

    # ---- surface water ------------------------------------------------------
    if not tidal:
        wet = info["depth"] > wet_tol
        eta = info["eta"]
        ax.fill_between(y, bed, eta, where=wet, interpolate=True, color=C_WATER,
                        alpha=0.85, zorder=3)
        ax.plot(y[wet], eta[wet], color=C_WATER_EDGE, lw=1.8, zorder=6)
        handles.append(Patch(facecolor=C_WATER, alpha=0.85, label="Sea"))
    else:
        t_end = info["t_end"]
        msl, hi, lo = tide_stats(cfg["tide_file"], t_end - cfg["tide_period"], t_end)
        reg_hi, reg_lo, reg_msl = bed < hi, bed < lo, bed < msl
        ax.fill_between(y, np.maximum(bed, lo), hi, where=reg_hi, color=C_WATER,
                        alpha=0.28, zorder=3)
        ax.fill_between(y, bed, lo, where=reg_lo, color=C_WATER, alpha=0.75, zorder=3)
        ax.plot(y[reg_msl], np.full(reg_msl.sum(), msl), color=C_WATER_EDGE,
                lw=1.6, zorder=6)
        handles += [Patch(facecolor=C_WATER, alpha=0.75, label="Sea (subtidal)"),
                    Patch(facecolor=C_WATER, alpha=0.28, label="Intertidal"),
                    Line2D([0], [0], color=C_WATER_EDGE, lw=1.6, label="Mean sea level")]

    # ---- topography + aquifer base ------------------------------------------
    ax.plot(y, bed, color=C_BED_EDGE, lw=2.2, zorder=7)
    ax.plot(y, rast["base"], color="0.4", lw=1.0, ls="--", zorder=3)

    # ---- experimental 18-psu interface (Kuan et al.) ------------------------
    if exp is not None:
        cs_main = ax.contour(Yc, rast["Zc"], rast["Cc"], levels=[half],
                             colors="none")            # for misfit geometry only
        dists, mae = interface_misfit(exp, cs_main)
        ax.plot(exp[0], exp[1], "o", mfc="white", mec="black", mew=1.0, ms=5.5,
                zorder=9)
        handles.append(Line2D([0], [0], ls="none", marker="o", mfc="white",
                              mec="black", mew=1.0, ms=5.5, label="Experiment"))
        print(f"  experiment overlay      : {exp[0].size} pts; MAE to simulated "
              f"{half:.0f}-psu contour = {mae*100:.2f} cm (max {dists.max()*100:.2f} cm)")

    # ---- axes ---------------------------------------------------------------
    ax.set_aspect("equal")
    ax.set_xlim(x_crop, y[-1] + dy / 2)
    ax.set_ylim(rast["base"].min() - 0.02, bed.max() + 0.03)
    if show_xlabel:
        ax.set_xlabel("Distance $y$ [m]")
    ax.set_ylabel("Elevation $z$ [m]")
    # Panel tag in the top-right white space (above the sloping beach); the
    # legend occupies the top-left fresh-water corner in both panels.
    ax.text(0.985, 0.94, tag, transform=ax.transAxes, va="top", ha="right",
            zorder=12, path_effects=[pe.withStroke(linewidth=2.0, foreground="white")])
    ax.legend(handles=handles, loc="upper left", framealpha=0.93, ncol=1)
    return pcm


# -----------------------------------------------------------------------------
# Figure: two stacked panels (ss on top, td below), shared salinity colorbar.
# -----------------------------------------------------------------------------
print("kuan-saltwater-intrusion — combined steady (ss) + tidal-mean (td) figure")
fig, (axT, axB) = plt.subplots(2, 1, figsize=(DOUBLE, 12 * cm),
                               sharex=True, sharey=True, constrained_layout=True)

draw_panel(axT, "ss", "(a)", show_xlabel=False)
pcm = draw_panel(axB, "td", "(b)", show_xlabel=True)

cbar = fig.colorbar(pcm, ax=(axT, axB), pad=0.02, fraction=0.05, aspect=28)
cbar.set_label("Salinity [psu]")

fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nwrote {out_pdf.resolve()}")

plt.show()
