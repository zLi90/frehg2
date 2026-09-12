"""
Compare the swe-vcatchment frehg2 outlet hydrograph against the ParFlow, HGS
and tRIBS hydrographs from Maxwell et al. (2014).

Case
----
Classic tilted-V ("V-catchment") overland-flow benchmark: two 800 m x 1000 m
planar hillslopes (x-slope 5%) draining into a central 20 m channel (y-slope
2%) to a single south-edge outlet. Uniform rain 3.0e-6 m/s (= 1.8e-4 m/min)
falls for 0-90 min, then drainage. Manning n = 0.015 / 0.15, no infiltration.
Parameters match Maxwell et al. (2014) Table 3 (frehg2 uses a 10 m mesh).

Reference solution
------------------
Digitized outlet hydrographs for three of the intercomparison models, from
Figure 6 of:
  Maxwell, R.M., et al. (2014), "Surface-subsurface model intercomparison: A
  first set of benchmark results to diagnose integrated hydrology and
  feedbacks", Water Resources Research 50(2), 1531-1549,
  doi:10.1002/2013WR013725.
Each reference/plot-data-<model>.csv holds two columns: x = time [min],
y = outlet discharge [m^3/min].

Method
------
The frehg2 outlet discharge is the time-derivative of the cumulative
boundary-outflow volume recorded (every dt = 2 s) in /monitor/mass_audit --
i.e. the exact outflow flux the solver applied at the outlet BC. This is
preferred over a naive sum of depth*velocity at the outlet cells, which
overestimates the flux (cell-centre velocities are not the face-normal flux).
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
REF_COLORS = ["#d1611f", "#2e8b57", "#7b3fa0", "#ff7f0e"]
MARKERS = ["o", "s", "^", "D"]
MS, MEW = 6.0, 1.0


# -----------------------------------------------------------------------------
# Configuration (from swe-vcatchment.yaml)
# -----------------------------------------------------------------------------
output_file = Path("out/output.h5")
reference_dir = Path("reference")
rain_off = 5400.0                       # [s]  rain stops at 90 min
SEC_PER_MIN = 60.0
out_pdf = Path("vcatchment_hydrograph_comparison.pdf")

# Digitized Maxwell (2014) Figure 6 hydrographs: (label, file, marker, colour).
references = [
    ("ParFlow", reference_dir / "plot-data-Parflow.csv", MARKERS[0], REF_COLORS[0]),
    ("HGS",     reference_dir / "plot-data-HGS.csv",     MARKERS[1], REF_COLORS[1]),
    ("tRIBS",   reference_dir / "plot-data-tRIBS.csv",   MARKERS[2], REF_COLORS[2]),
]


# -----------------------------------------------------------------------------
# frehg2 outlet hydrograph from the cumulative boundary-outflow volume
# -----------------------------------------------------------------------------
with h5py.File(output_file, "r") as f:
    ma = f["/monitor/mass_audit"][:]
# columns: time, volume, rain, evaporation, boundary_outflow, bc_inflow, clamped
t_s = ma[:, 0]
out_cum = ma[:, 4]

t_min = t_s / SEC_PER_MIN
Q_min = np.gradient(out_cum, t_s) * SEC_PER_MIN     # outlet discharge [m^3/min]


# -----------------------------------------------------------------------------
# Load the reference hydrographs; report curve-wide agreement (not peak)
# -----------------------------------------------------------------------------
ref_data = []
print("swe-vcatchment — frehg2 vs Maxwell et al. (2014) ParFlow / HGS / tRIBS")
print(f"  frehg2 peak outflow : {Q_min.max():.1f} m^3/min "
      f"({Q_min.max()/60:.3f} m^3/s) at t = {t_min[np.argmax(Q_min)]:.1f} min")
for label, path, marker, colour in references:
    d = np.genfromtxt(path, delimiter=",", skip_header=1)
    rt, rq = d[:, 0], d[:, 1]                        # time [min], outflow [m^3/min]
    ref_data.append((label, rt, rq, marker, colour))
    # Mean absolute difference vs frehg2 sampled at the reference times.
    mad = np.mean(np.abs(np.interp(rt, t_min, Q_min) - rq))
    print(f"  {label:<8}: {len(rt):2d} pts, t = {rt.min():.0f}-{rt.max():.0f} min, "
          f"mean |frehg2 - {label}| = {mad:.1f} m^3/min")


# -----------------------------------------------------------------------------
# Figure: single-panel hydrograph comparison
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(SINGLE, 6.6 * cm))

# Rain-on context (0-90 min) as a light shaded band (no on-figure label).
ax.axvspan(0, rain_off / SEC_PER_MIN, color="#8fb4d8", alpha=0.12)

# frehg2 as a continuous line.
ax.plot(t_min, Q_min, color=C_SIM, lw=1.8, zorder=4, label="frehg2")

# Reference models as markers.
for label, rt, rq, marker, colour in ref_data:
    ax.scatter(rt, rq, s=34, marker=marker, facecolor=colour, edgecolor="white",
               linewidth=0.6, zorder=5, label=label)

ax.set_xlim(0, 200)
ax.set_ylim(0, 320)
ax.set_xlabel("Time [min]")
ax.set_ylabel(r"Outlet discharge [m$^3$ min$^{-1}$]")
ax.grid(True, lw=0.3, alpha=0.3)
ax.legend(loc="upper right", framealpha=0.92)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nwrote {out_pdf.resolve()}")

plt.show()
