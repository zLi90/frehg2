"""
Visualize and verify swe-macdonald-subcritical against the SWASHES analytic
solution.

Case
----
SWASHES / MacDonald steady *subcritical* flow in a long (L = 1000 m) channel of
varying bed, unit discharge q = 2 m^2/s, Manning n = 0.033, g = 9.81 m/s^2.
The bed z(x) is constructed by MacDonald's trick so that a prescribed
water-depth profile is an *exact* steady solution of the shallow-water
equations:

    h_exact(x) = (4/g)^(1/3) * ( 1 + 1/2 * exp( -16 (x/L - 1/2)^2 ) )

with a spatially constant unit discharge q(x) = 2 m^2/s everywhere. Endpoints
h_exact(0) = h_exact(L) = 0.748378 m (this is the imposed downstream depth in
the setup); the profile bulges to 1.1123 m at mid-channel.

Reference: Delestre et al. (2013), "SWASHES: a compilation of shallow-water
analytic solutions for hydraulic and environmental studies", Int. J. Numer.
Meth. Fluids 72(3), MacDonald-based steady solutions.

The figure (single panel) shows the converged longitudinal section: the frehg2
water-surface elevation vs the SWASHES analytic profile over the bed. The full
verification (depth L1/L2/Linf error norms, unit-discharge conservation, and the
subcritical Froude range) is reported on the console and quoted in the caption.

HDF5 layout (frehg2 surface output):
    /surface/{eta,depth,uu,vv}/<time_seconds>   flat ny*nx, index = j*nx + i
    /grid/bottom                                 static bed elevation, same layout
    /grid/x_center                               cell-center x [m]
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
C_EXACT = "#000000"
C_BED, C_BED_EDGE = "#c8a06e", "#6b4f2a"


# -----------------------------------------------------------------------------
# Configuration / physical parameters (from swe-macdonald-subcritical.yaml)
# -----------------------------------------------------------------------------
output_file = Path("out/output.h5")
g = 9.81              # gravity [m/s^2]           (surface_water.gravity)
q_target = 2.0        # unit discharge [m^2/s]    (inflow 8.0 m^3/s / 4 m width)
L = 1000.0            # channel length [m]        (nx * dx)
manning_n = 0.033     # for reference only
out_pdf = Path("macdonald_subcritical_verification.pdf")

# Acceptance tolerances for the console PASS/FAIL summary.
tol_depth_Linf = 1.0e-2   # [m]   max pointwise depth error
tol_q_rel = 0.02          # [-]   max relative unit-discharge deviation


# -----------------------------------------------------------------------------
# SWASHES / MacDonald closed-form subcritical depth
# -----------------------------------------------------------------------------
def macdonald_subcritical_depth(x):
    """Analytic steady depth h_exact(x) [m] for the SWASHES long-channel
    subcritical flume (q = 2 m^2/s, g = 9.81)."""
    return (4.0 / g) ** (1.0 / 3.0) * (1.0 + 0.5 * np.exp(-16.0 * (x / L - 0.5) ** 2))


# -----------------------------------------------------------------------------
# Load the frehg2 output (steady = last snapshot), averaged across the width
# -----------------------------------------------------------------------------
with h5py.File(output_file, "r") as f:
    x = f["/grid/x_center"][:]
    nx = x.size
    ny = f["/grid/y_center"].shape[0]

    def profile(name, t):
        """Longitudinal profile of a surface field, averaged over the ny rows."""
        return f[f"/surface/{name}/{t}"][:].reshape(ny, nx).mean(axis=0)

    bed = f["/grid/bottom"][:].reshape(ny, nx).mean(axis=0)

    times = sorted(int(t) for t in f["/surface/depth"].keys())
    t_end = times[-1]
    t_prev = times[-2] if len(times) > 1 else times[-1]

    eta = profile("eta", t_end)
    depth = profile("depth", t_end)
    u = profile("uu", t_end)
    v = profile("vv", t_end)
    depth_prev = profile("depth", t_prev)          # steadiness check


# -----------------------------------------------------------------------------
# Reference solution + derived quantities
# -----------------------------------------------------------------------------
h_exact = macdonald_subcritical_depth(x)
eta_exact = bed + h_exact

depth_err = depth - h_exact                        # [m]
q_sim = depth * u                                  # unit discharge [m^2/s]
froude = np.abs(u) / np.sqrt(g * np.maximum(depth, 1e-12))

# Error norms on depth.
L1 = np.mean(np.abs(depth_err))
L2 = np.sqrt(np.mean(depth_err ** 2))
Linf = np.max(np.abs(depth_err))
i_worst = int(np.argmax(np.abs(depth_err)))

q_dev = np.max(np.abs(q_sim - q_target))
steady = np.max(np.abs(depth - depth_prev))


# -----------------------------------------------------------------------------
# Console verification report
# -----------------------------------------------------------------------------
print(f"swe-macdonald-subcritical — verification vs SWASHES analytic solution")
print(f"  steady state       : |h(t={t_end}) - h(t={t_prev})|_inf = {steady:.2e} m")
print(f"  depth h(x)         : L1={L1:.3e}  L2={L2:.3e}  Linf={Linf:.3e} m"
      f"  (worst at x={x[i_worst]:.0f} m)")
print(f"  unit discharge q   : mean={q_sim.mean():.4f}  max|q-{q_target:.1f}|={q_dev:.3e} m^2/s"
      f"  ({100*q_dev/q_target:.2f}% of target)")
print(f"  Froude number      : min={froude.min():.3f}  max={froude.max():.3f}"
      f"  -> {'subcritical throughout' if froude.max() < 1 else 'NOT subcritical!'}")
print(f"  cross-flow |v|max   : {np.max(np.abs(v)):.2e} m/s (≈0 expected for 1D)")
ok = (Linf < tol_depth_Linf) and (q_dev / q_target < tol_q_rel) and (froude.max() < 1.0)
print(f"  RESULT             : {'PASS' if ok else 'CHECK'} "
      f"(depth Linf < {tol_depth_Linf} m, q within {100*tol_q_rel:.0f}%, Fr < 1)")


# -----------------------------------------------------------------------------
# Figure: single longitudinal section (bed + water surface: sim vs analytic)
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(DOUBLE, 6.5 * cm))

ax.fill_between(x, bed.min() - 0.3, bed, color=C_BED, zorder=1)
ax.plot(x, bed, color=C_BED_EDGE, lw=1.2, zorder=3, label="Bed $z(x)$")
ax.plot(x, eta_exact, color=C_EXACT, ls="--", lw=1.6, zorder=4,
        label=r"Analytic $\eta$ (SWASHES)")
ax.plot(x, eta, color=C_SIM, lw=1.4, zorder=5, label=r"frehg2 $\eta$")

ax.set_xlim(x.min(), x.max())
ax.set_xlabel("Distance $x$ [m]")
ax.set_ylabel(r"Water-surface elevation $\eta$ [m]")
ax.grid(True, lw=0.3, alpha=0.3)
ax.legend(loc="upper right", ncol=3, framealpha=0.9)

fig.tight_layout()
fig.savefig(out_pdf, bbox_inches="tight")
print(f"\nwrote {out_pdf.resolve()}")

plt.show()
