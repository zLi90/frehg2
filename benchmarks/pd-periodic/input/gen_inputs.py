#!/usr/bin/env python3
"""gen_inputs.py — pd-periodic raster generator (Stage-3 T3).

pd2-surface: a 12 x 6 channel (dx = dy = 5 m), periodic in x, with an
x-dependent bed bump, a tracer slab crossing the seam under wind, and
every x-dependent input EXACTLY periodic on the 60 m torus (the roll
identity the gate asserts is only exact for periodic inputs).
pd2-gw: a 12 x 6 x 8 confined block with a periodic K(x, y) modulation
expressed as 8 discrete soil types (ids quantize the modulation; the ids
are what must satisfy the roll identity, and a quantized periodic field
stays periodic). pd5 lanes reuse the pd2-gw map plus a uniform one.
Rolled variants (--roll K) shift every x-dependent raster K columns for
the pd2 roll-equivalence lanes; the gate generates them on the fly, this
script documents the convention. Rasters are flat lists j*nx + i
(2-D) / (j*nx + i)*nz + k (3-D), docs/user-guide/input-data.md.
"""
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
NX, NY, NZ = 12, 6, 8
LX = 60.0


def write(name, values, comment, per_line=12):
    lines = [f"# {comment}"]
    vals = [f"{v:.8f}" if isinstance(v, float) else str(v) for v in values]
    for n in range(0, len(vals), per_line):
        lines.append(" ".join(vals[n:n + per_line]))
    (HERE / name).write_text("\n".join(lines) + "\n")


# pd2-surface: bed bump b(x) = -1.0 + 0.2 exp(-(2 sin(pi (x-20)/60))^2)
# (periodic by construction on [0, 60)); tracer slab s = 10 in columns
# i = 10..11 (crosses the seam under +x wind).
bed = []
for j in range(NY):
    for i in range(NX):
        x = (i + 0.5) * 5.0
        bed.append(-1.0 + 0.2 * math.exp(-(2.0 * math.sin(math.pi * (x - 20.0) / LX))**2))
write("pd2_dem.dat", bed, "pd2 bed [m], j*nx+i: periodic bump on the 60 m torus")

tracer = [10.0 if i >= 10 else 0.0 for j in range(NY) for i in range(NX)]
write("pd2_tracer0.dat", tracer, "pd2 initial surface tracer: slab at the seam (i = 10..11)")

# pd2-gw: K modulation 1 + 0.5 sin(2 pi x / Lx) cos(pi y / Ly), quantized
# to 8 soil types (type k has K = 1e-5 * (0.5 + k/7)); EXACTLY periodic
# in x by construction of the ids.
ids = []
for j in range(NY):
    for i in range(NX):
        x = (i + 0.5) * 5.0
        y = (j + 0.5) * 5.0
        mod = 1.0 + 0.5 * math.sin(2.0 * math.pi * x / LX) * math.cos(math.pi * y / 30.0)
        k = min(7, max(0, int(round((mod - 0.5) / 1.0 * 7.0 / 1.0))))
        ids.extend([k] * NZ)
write("pd2_soil_id.dat", ids, "pd2 soil ids, (j*nx+i)*nz+k: quantized periodic K modulation", 16)

# pd2-coupled: a 12 x 8 (dx = dy = 2 m, Lx = 24 m) bed sloping down toward
# the north tide with an EXACTLY periodic x bump — the only x-dependent
# input of the coupled lane (ICs and soil are uniform), so the roll
# identity needs just this raster rolled. ALWAYS WET under the ±0.2 m
# tide (bed <= -0.95 m): the first authoring (bed up to -0.18 m) put a
# marginally-drying row under a vg_n = 1.23 soil and the trajectory was
# CHAOTIC (Lyapunov time ~40 s — a 1e-12 single-cell bed perturbation
# reproduced the roll deltas exactly: eta 3e-4, uu 1e-2, head 2.0 m), so
# no roll identity could gate it (S3-A6).
cbed = []
for j in range(8):
    for i in range(NX):
        x = (i + 0.5) * 2.0
        y = (j + 0.5) * 2.0
        cbed.append(-1.0 - 0.3 * (y / 16.0) + 0.05 * math.cos(2.0 * math.pi * x / 24.0))
write("pd2c_dem.dat", cbed, "pd2-coupled bed [m], j*nx+i: y slope + periodic x bump on the 24 m torus, always wet")

# pd8 orient-base (surface battery): an 8 x 8 (dx = dy = 5 m, 40 m torus)
# bed with an exactly-periodic x bump and a tracer slab crossing the seam
# under the oblique wind (u10 = 6, v10 = 3; the battery rotates the wind
# vector with the domain).
obed = []
for j in range(8):
    for i in range(8):
        x = (i + 0.5) * 5.0
        obed.append(-1.0 + 0.2 * math.exp(-(2.0 * math.sin(math.pi * (x - 20.0) / 40.0))**2))
write("pd8_dem.dat", obed, "pd8 orient bed [m], j*nx+i: periodic bump on the 40 m torus", 8)

otracer = [10.0 if i >= 6 else 0.0 for j in range(8) for i in range(8)]
write("pd8_tracer0.dat", otracer,
      "pd8 initial surface tracer: slab at the seam (i = 6..7)", 8)

print("pd rasters written")
