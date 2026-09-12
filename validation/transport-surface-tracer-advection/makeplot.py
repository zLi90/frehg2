"""
Visualize and verify surface-tracer-advection: a passive scalar carried by a
steady, uniform open-channel flow, advected with the first-order UPWIND vs the
second-order TVD SUPERBEE scheme, both compared to the ideal analytic solution.

Case
----
A 1000 m x 1 m channel (nx = 200, dx = 5 m) on a mild 0.1 % slope. Manning
n = 0.03 gives normal depth h = 0.5 m and a uniform velocity u = q/h =
0.33202 / 0.5 = 0.664 m/s (unit discharge q = 0.332 m^2/s). The stage is matched
to the inflow so the flow is uniform from t = 0. At t = 0 the inflow
concentration steps 0 -> 1. With no physical diffusion
(surface_diffusivity ~ 1e-10) the exact answer is a sharp step front that
translates downstream at u:

    c_exact(x, t) = 1   for x <  u*t          (behind the front)
                    0   for x >= u*t          (ahead of the front)

Any spreading of the simulated front is therefore *purely numerical* diffusion
introduced by the advection scheme -- upwind smears it over many cells,
superbee keeps it sharp. Both must translate at the correct celerity u. (Front
positions and 10-90% smear widths are reported on the console / in the caption.)

Reference: for the dispersive analogue the profile is the Ogata & Banks (1961)
erfc solution (USGS Prof. Paper 411-A); here D ~ 0 so the ideal is the step.

The figure (single panel)
-------------------------
Concentration profiles c(x) at five snapshots (t = 200, 400, 600, 800, 1000 s).
  * colour  encodes time  -> the front marches downstream at u (advection),
  * dotted  = ideal analytic step (x = u*t),
  * dashed  = upwind  (1st order),
  * solid   = superbee (TVD 2nd order).

HDF5 layout (frehg2 output):
    /transport/concentration_surface/<time_seconds>   flat length nx (ny = 1)
    /surface/{depth,uu}/<time_seconds>                same layout
    /grid/x_center                                     cell-centre x [m]
"""

from pathlib import Path

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
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


# -----------------------------------------------------------------------------
# Configuration (from advect-*.yaml)
# -----------------------------------------------------------------------------
runs = {
    "upwind":   Path("out-upwind/output.h5"),
    "superbee": Path("out-superbee/output.h5"),
}
q_unit = 0.33202          # unit discharge [m^2/s]  (inflow value / 1 m width)
h_normal = 0.5            # normal depth [m]
u = q_unit / h_normal     # advection velocity [m/s]  = 0.66404

times_plot = [200, 400, 600, 800, 1000]   # snapshots [s]
out_pdf = Path("tracer_advection_comparison.pdf")

# Per-scheme line style (colour comes from the time colormap).
styles = {
    "upwind":   dict(ls="--", lw=1.4, name="Upwind (1st-order)"),
    "superbee": dict(ls="-",  lw=1.6, name="Superbee (TVD)"),
}
cmap = plt.cm.viridis
norm = Normalize(vmin=times_plot[0] - 60, vmax=times_plot[-1] + 20)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def analytic_step(x_edges_lo, x_edges_hi, t):
    """Ideal pure-advection profile drawn as a crisp step at x = u*t:
    four points [0 -> front (c=1)] then [front -> end (c=0)]."""
    xf = u * t
    xs = np.array([x_edges_lo, xf, xf, x_edges_hi])
    cs = np.array([1.0, 1.0, 0.0, 0.0])
    return xs, cs, xf


def rightmost_crossing(x, c, level):
    """x at which the (1 -> 0) profile last crosses `level`, linearly
    interpolated. Robust to small non-monotone wiggles in the tail."""
    s = (c[:-1] - level) * (c[1:] - level)
    idx = np.where(s <= 0.0)[0]
    if idx.size == 0:
        return np.nan
    i = idx[-1]
    c0, c1, x0, x1 = c[i], c[i + 1], x[i], x[i + 1]
    if c1 == c0:
        return x0
    return x0 + (level - c0) * (x1 - x0) / (c1 - c0)


# -----------------------------------------------------------------------------
# Load both runs
# -----------------------------------------------------------------------------
data = {}
for scheme, path in runs.items():
    with h5py.File(path, "r") as f:
        x = f["/grid/x_center"][:]
        conc = {t: f[f"/transport/concentration_surface/{t}"][:] for t in times_plot}
        uu_final = f["/surface/uu/1000"][:]
    data[scheme] = dict(x=x, conc=conc, uu_final=uu_final)

x = data["upwind"]["x"]
x_lo, x_hi = 0.0, x[-1] + 2.5     # domain edges for the analytic step


# -----------------------------------------------------------------------------
# Console verification report (final time t = 1000 s)
# -----------------------------------------------------------------------------
print("surface-tracer-advection — frehg2 vs ideal pure-advection step")
print(f"  advection velocity u = q/h = {q_unit}/{h_normal} = {u:.5f} m/s")
u_real = data["upwind"]["uu_final"][2:-2].mean()   # interior mean, avoid BC cells
print(f"  realized surface velocity (interior mean, t=1000) = {u_real:.5f} m/s "
      f"({100*(u_real-u)/u:+.2f}% vs analytic)")
print(f"  exact front at t = 1000 s : x_f = u*t = {u*1000:.1f} m\n")
print(f"  {'scheme':<9} {'front c=0.5':>12} {'10-90% smear':>14} {'L1 err':>10} {'c_max':>7}")
front_final = {}
for scheme in ("upwind", "superbee"):
    c = data[scheme]["conc"][1000]
    xf = rightmost_crossing(x, c, 0.5)
    x10 = rightmost_crossing(x, c, 0.1)     # downstream toe
    x90 = rightmost_crossing(x, c, 0.9)     # upstream shoulder
    smear = x10 - x90
    c_exact_cells = np.where(x < u * 1000.0, 1.0, 0.0)
    l1 = np.mean(np.abs(c - c_exact_cells))
    front_final[scheme] = xf
    print(f"  {scheme:<9} {xf:9.1f} m  {smear:9.1f} m  {l1:9.4f} {c.max():7.4f}")
print(f"\n  front error vs exact: upwind {front_final['upwind']-u*1000:+.1f} m, "
      f"superbee {front_final['superbee']-u*1000:+.1f} m (dx = 5 m)")


# -----------------------------------------------------------------------------
# Figure: single panel, all three (analytic / upwind / superbee) x 5 snapshots
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(DOUBLE, 6.5 * cm))

for t in times_plot:
    colour = cmap(norm(t))
    bold = (t == times_plot[-1])
    a_lw = 1.8 if bold else 1.0
    lw_scale = 1.3 if bold else 1.0
    alpha = 1.0 if bold else 0.7

    # Ideal analytic step (dotted).
    xs, cs, xf = analytic_step(x_lo, x_hi, t)
    ax.plot(xs, cs, ls=":", lw=a_lw, color=colour, alpha=alpha, zorder=3)

    # Numerical schemes.
    for scheme, st in styles.items():
        ax.plot(x, data[scheme]["conc"][t], ls=st["ls"],
                lw=st["lw"] * lw_scale, color=colour, alpha=alpha,
                zorder=5 if scheme == "superbee" else 4)

# Colourbar: time.
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, pad=0.015, ticks=times_plot)
cbar.set_label("Time $t$ [s]")
cbar.ax.set_yticklabels([f"{t}" for t in times_plot])

# Style / quantity legend (colour-neutral proxies).
proxies = [
    Line2D([0], [0], color="0.25", ls=":",  lw=1.6, label="Analytic (ideal step)"),
    Line2D([0], [0], color="0.25", ls="--", lw=1.6, label="Upwind (1st-order)"),
    Line2D([0], [0], color="0.25", ls="-",  lw=1.6, label="Superbee (TVD)"),
]
ax.legend(handles=proxies, loc="upper right", framealpha=0.93)

ax.set_xlim(0, 1000)
ax.set_ylim(-0.03, 1.08)
ax.set_xlabel("Distance along channel $x$ [m]")
ax.set_ylabel("Surface tracer concentration $c$ [--]")
ax.grid(True, lw=0.3, alpha=0.3)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nwrote {out_pdf.resolve()}")

plt.show()
