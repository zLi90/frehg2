"""
Visualize the swe-lake-at-rest (2D well-balancing / C-property) result.

Reads the frehg2 HDF5 output and draws a single vertical transect taken at
x = 0.5 * (domain width), i.e. the column i = nx // 2. Each of the two panels
shows the lake level (surface elevation eta) resting on the bed topography --
panel (a) at the start of the run, panel (b) at the end. A flat water surface
that is visually identical between the two panels is the "lake at rest"
signature: the still lake does not move. (The exact start/end drift is reported
on the console and quoted in the figure caption.)

HDF5 layout (frehg2 surface output):
    /surface/eta/<time_seconds>     flat length ny*nx, index = j*nx + i
    /surface/depth/<time_seconds>   same layout (used as the wet/dry mask)
    /grid/bottom                    static bed elevation, same layout
    /grid/x_center, /grid/y_center  cell-center coordinates [m]
"""

from pathlib import Path

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


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
C_SIM = "#1f5fa6"
C_BED, C_BED_EDGE = "#c8a06e", "#6b4f2a"
C_WATER, C_WATER_EDGE = "#7ba7d7", "#1f5fa6"


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
output_file = Path("out/output.h5")
x_fraction = 0.5      # transect location as a fraction of the x-extent
depth_tol = 1.0e-9    # wet/dry threshold [m]
reference_level = 3.0  # nominal lake surface eta [m] (datum line only)
out_pdf = Path("lake_at_rest_transect.pdf")


# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
with h5py.File(output_file, "r") as f:
    xc = f["/grid/x_center"][:]
    yc = f["/grid/y_center"][:]
    nx, ny = xc.size, yc.size

    # Static bed elevation, reshaped to (ny, nx) with index = j * nx + i.
    bed = f["/grid/bottom"][:].reshape(ny, nx)

    # Snapshot times (integer seconds) in numerical order; take first and last.
    times = sorted(int(t) for t in f["/surface/eta"].keys())
    t_start, t_end = times[0], times[-1]

    def surface_field(name, t):
        return f[f"/surface/{name}/{t}"][:].reshape(ny, nx)

    eta = {t: surface_field("eta", t) for t in (t_start, t_end)}
    depth = {t: surface_field("depth", t) for t in (t_start, t_end)}


# -----------------------------------------------------------------------------
# Extract the transect at x = x_fraction of the domain
# -----------------------------------------------------------------------------
i_slice = int(round(x_fraction * (nx - 1)))
x_slice = xc[i_slice]

bed_col = bed[:, i_slice]


def lake_level(t):
    """eta along the transect, masked to NaN where the cell is dry so the
    water-surface line breaks over emergent peaks instead of tracing the bed."""
    e = eta[t][:, i_slice].astype(float).copy()
    e[depth[t][:, i_slice] <= depth_tol] = np.nan
    return e


# -----------------------------------------------------------------------------
# Plot: two panels (start / end), shared axes
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 7 * cm), sharex=True, sharey=True)

y_bottom = min(bed_col.min(), 0.0) - 0.15
y_top = max(bed_col.max(), reference_level) + 0.35

panels = [(axes[0], t_start, "(a)"), (axes[1], t_end, "(b)")]

for ax, t, tag in panels:
    wet = depth[t][:, i_slice] > depth_tol
    eta_col = eta[t][:, i_slice]
    level = lake_level(t)

    # Topography: solid fill from the axis floor up to the bed.
    ax.fill_between(yc, y_bottom, bed_col, color=C_BED, zorder=1)
    ax.plot(yc, bed_col, color=C_BED_EDGE, lw=1.3, zorder=3, label="Bed")

    # Water body: fill between bed and eta only where the cell is wet.
    ax.fill_between(yc, bed_col, eta_col, where=wet, interpolate=True,
                    color=C_WATER, alpha=0.85, zorder=2)

    # Water surface line (broken over the emergent, dry peaks).
    ax.plot(yc, level, color=C_WATER_EDGE, lw=1.8, zorder=4,
            label=r"Lake level $\eta$")

    # Nominal still-water datum for reference.
    ax.axhline(reference_level, color="k", ls=":", lw=0.8, alpha=0.6, zorder=5)

    ax.text(0.03, 0.94, tag, transform=ax.transAxes, va="top", ha="left")
    ax.set_xlabel("Distance $y$ [m]")

axes[0].set_ylabel("Elevation $z$ [m]")
axes[0].set_ylim(y_bottom, y_top)
axes[0].set_xlim(yc.min(), yc.max())
axes[0].legend(loc="upper right", framealpha=0.9)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
print(f"wrote {out_pdf.resolve()}")

# End-vs-start drift, for the console verification report (goes in the caption).
drift_col = np.nanmax(np.abs(lake_level(t_end) - reference_level))
d = np.max(np.abs(eta[t_end] - eta[t_start]))
print(f"transect at x = {x_slice:.1f} m (i = {i_slice})")
print(f"max|eta - {reference_level:.1f}| on transect (t_end) = {drift_col:.3e} m")
print(f"max|eta(end) - eta(start)| over all cells = {d:.3e} m")

plt.show()
