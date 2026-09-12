"""
Visualize the henry-saltwater-intrusion result: the steady salt-wedge
concentration field and its 0.5*C_max isochlor at the end of the simulation.

Case
----
The classic Henry (1964) seawater-intrusion benchmark for variable-density
subsurface transport (groundwater + transport + density coupling, no surface
module). A 2 m (landward y- -> seaward y+) x 1 m (vertical z) confined, fully
saturated box (nx=1, ny=40, nz=20, dy=dz=0.05 m). Freshwater recharges the
inland face (flux 6.6e-5 m/s, s=0); hydrostatic seawater (s_sea=33.6, the
density surrogate for Henry's drho/rho=0.025) enters the seaward face. The flow
reaches a steady salt wedge: seawater along the base, fresh discharge over the
top.

Reference
---------
Henry, H. R. (1964), "Effects of dispersion on salt encroachment in coastal
aquifers", USGS Water-Supply Paper 1613-C, pp. C71-C84. Standard verification
metric: the 0.5-isochlor (50 % seawater) toe intrudes ~0.8-1.0 m inland of the
sea along the base in the a=2 domain.

The figure (single panel)
-------------------------
The steady (t = t_end) concentration field on the y-z cross-section
(pcolormesh), overlaid with the **0.5*C_max isochlor** (bold) that marks the
wedge front. Two reference 0.5-isochlors -- the Henry semi-analytical solution
(plot-data-henry.csv) and the SUTRA variable-density code (plot-data-sutra.csv)
-- are overlaid as markers for verification. Each reference's mean distance
(MAE) to the simulated isochlor is reported on the console / in the caption.

HDF5 layout (frehg2 subsurface output)
--------------------------------------
    /transport/concentration/<time_seconds>   flat length ny*nz, index j*nz + k
                                               -> reshape (ny, nz); k=0 top, k=19 bottom
    /grid/y_center   (ny,)  landward->seaward cell-centre y [m]
    /grid/z_center   (nz,)  cell-centre z [m] (negative down; -0.025 top .. -0.975 bottom)
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
C_IFACE = "#d1201f"
REF_COLORS = ["#d1611f", "#2e8b57", "#7b3fa0", "#ff7f0e"]
MS, MEW = 6.0, 1.0


# -----------------------------------------------------------------------------
# Configuration (from henry.yaml / README)
# -----------------------------------------------------------------------------
output_file = Path("out/output.h5")
s_sea = 33.6                 # seaward salinity (density surrogate), = bounds.max
sea_y = 2.0                  # seaward boundary position [m] (toe distance datum)
out_pdf = Path("henry_saltwedge_field.pdf")

# Reference 0.5-isochlors (50 % seawater) overlaid for model verification: the
# Henry semi-analytical solution and the SUTRA variable-density code. Each CSV
# is "x, y" = (seaward distance [m], elevation [m]) tracing the wedge front, so
# the columns map straight onto this figure's axes (no coordinate swap).
references = [
    dict(name="Henry semi-analytical", short="Henry",
         csv=Path("plot-data-henry.csv"),
         marker="o", mfc="white",       mec="black"),
    dict(name="SUTRA variable-density code", short="SUTRA",
         csv=Path("plot-data-sutra.csv"),
         marker="D", mfc=REF_COLORS[3], mec="black"),
]


# -----------------------------------------------------------------------------
# Reference-solution overlay helpers (loader + interface misfit)
# -----------------------------------------------------------------------------
def load_reference(path):
    """Digitized 0.5-isochlor points: CSV columns (x = seaward distance,
    y = elevation) -> (dist, elev) arrays. None if the file is absent."""
    if not path.exists():
        return None
    d = np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
    return d[:, 0], d[:, 1]


def contour_segments(cs):
    """Vertices of a ContourSet as a list of (Ni, 2) polyline arrays."""
    if hasattr(cs, "get_paths"):                        # Matplotlib >= 3.8
        return [p.vertices for p in cs.get_paths() if len(p.vertices) > 1]
    segs = []                                           # pragma: no cover (older mpl)
    for coll in cs.collections:
        segs += [seg for seg in coll.get_segments() if len(seg) > 1]
    return segs


def point_to_polyline(px, pz, segs):
    """Shortest Euclidean distance from point (px, pz) to any contour segment."""
    best = np.inf
    for v in segs:
        a, b = v[:-1], v[1:]                            # segment endpoints
        ab = b - a
        ap = np.column_stack([px - a[:, 0], pz - a[:, 1]])
        denom = (ab ** 2).sum(1)
        t = np.clip((ap * ab).sum(1) / np.where(denom > 0, denom, 1.0), 0.0, 1.0)
        proj = a + t[:, None] * ab
        best = min(best, np.hypot(px - proj[:, 0], pz - proj[:, 1]).min())
    return best


def interface_misfit(ref, cs):
    """Per-point distance from each reference 0.5-isochlor point to the simulated
    isochlor `cs`, and their mean (the interface MAE)."""
    segs = contour_segments(cs)
    if ref is None or not segs:
        return None, None
    xe, ze = ref
    dists = np.array([point_to_polyline(xe[i], ze[i], segs) for i in range(xe.size)])
    return dists, float(dists.mean())


# -----------------------------------------------------------------------------
# Load the final (steady) field + the previous snapshot for a steadiness check
# -----------------------------------------------------------------------------
with h5py.File(output_file, "r") as f:
    y = f["/grid/y_center"][:]                     # (ny,) landward -> seaward
    z = f["/grid/z_center"][:]                     # (nz,) top (-0.025) -> bottom (-0.975)
    ny, nz = y.size, z.size

    times = sorted(int(t) for t in f["/transport/concentration"].keys())
    t_end, t_prev = times[-1], times[-2]

    # flat index = j*nz + k  ->  reshape (ny, nz), then transpose to (nz, ny)
    # so rows follow z (k) and columns follow y (j) for meshgrid plotting.
    C = f[f"/transport/concentration/{t_end}"][:].reshape(ny, nz).T      # (nz, ny)
    C_prev = f[f"/transport/concentration/{t_prev}"][:].reshape(ny, nz).T


# -----------------------------------------------------------------------------
# Derived quantities + verification metrics
# -----------------------------------------------------------------------------
cmax = C.max()
half = 0.5 * cmax                                   # the requested 0.5*max level

Y, Z = np.meshgrid(y, z)                            # (nz, ny)

# 0.5-isochlor toe along the base: landward-most y where the bottom row reaches
# `half`. The base salinity increases monotonically seaward, so interpolate.
base = C[-1, :]                                     # bottom row (k = nz-1) vs y
toe_y = np.interp(half, base, y)                    # base increasing with y
toe_from_sea = sea_y - toe_y

# Steadiness: max change over the last output interval.
dC = np.max(np.abs(C - C_prev))

# Reference 0.5-isochlors (loaded now; misfit computed once the simulated
# isochlor contour exists in the figure section).
for ref in references:
    ref["pts"] = load_reference(ref["csv"])


# -----------------------------------------------------------------------------
# Console verification report
# -----------------------------------------------------------------------------
print("henry-saltwater-intrusion — steady salt-wedge field")
print(f"  final time            : t = {t_end} s")
print(f"  concentration range   : [{C.min():.2f}, {cmax:.2f}]  (s_sea = {s_sea})")
print(f"  steadiness            : max|C(t={t_end}) - C(t={t_prev})| = {dC:.3f}")
print(f"  0.5*C_max level        : {half:.2f}")
print(f"  0.5-isochlor toe (base): y = {toe_y:.3f} m  ->  {toe_from_sea:.3f} m from the sea")
print(f"                           (Henry reference band ~0.8-1.0 m; README ~0.73 m)")
for ref in references:
    n = 0 if ref["pts"] is None else ref["pts"][0].size
    print(f"  reference solution      : {ref['name']:<28} ({n} pts, {ref['csv'].name})")


# -----------------------------------------------------------------------------
# Figure: concentration field + 0.5*C_max isochlor + reference isochlors
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(DOUBLE, 7 * cm))

# (1) Concentration field (rasterized so the PDF stays small; text/lines vector).
pcm = ax.pcolormesh(Y, Z, C, shading="gouraud", cmap="viridis",
                    vmin=0.0, vmax=cmax, rasterized=True)
cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
cbar.set_label("Concentration $c$ [--]")

# (2) The 0.5*C_max isochlor (bold, the wedge front).
cs_half = ax.contour(Y, Z, C, levels=[half], colors=C_IFACE, linewidths=2.2)
if hasattr(cs_half, "set_path_effects"):
    cs_half.set_path_effects([pe.withStroke(linewidth=3.6, foreground="white")])
else:                                                        # pragma: no cover
    for coll in cs_half.collections:
        coll.set_path_effects([pe.withStroke(linewidth=3.6, foreground="white")])

# (3) Reference 0.5-isochlors (Henry semi-analytical + SUTRA) overlaid as
#     markers, with each reference's interface misfit (mean point-to-contour
#     distance) against the simulated 0.5*C_max isochlor `cs_half`.
sim_handle = Line2D([], [], color=C_IFACE, lw=2.2,
                    label=r"frehg2 (0.5$\,C_{\max}$)")
legend_handles = [sim_handle]
print("\n  0.5-isochlor verification (model vs reference solutions):")
for ref in references:
    if ref["pts"] is None:
        print(f"    {ref['name']:<28}: {ref['csv'].name} not found — skipped")
        continue
    xe, ze = ref["pts"]
    dists, mae = interface_misfit(ref["pts"], cs_half)
    ref["mae"] = mae
    h, = ax.plot(xe, ze, marker=ref["marker"], mfc=ref["mfc"], mec=ref["mec"],
                 mew=MEW, ms=MS, ls="none", zorder=9, label=ref["short"])
    legend_handles.append(h)
    print(f"    {ref['name']:<28}: {xe.size} pts, MAE = {mae*100:.2f} cm "
          f"(max {dists.max()*100:.2f} cm)")

ax.set_aspect("equal")                             # honest 2:1 cross-section
ax.set_xlim(0.0, sea_y)
ax.set_ylim(z.min() - 0.025, 0.0)
ax.set_xlabel("Distance $y$ [m]")
ax.set_ylabel("Elevation $z$ [m]")
ax.legend(handles=legend_handles, loc="upper left", framealpha=0.93)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nwrote {out_pdf.resolve()}")

plt.show()
