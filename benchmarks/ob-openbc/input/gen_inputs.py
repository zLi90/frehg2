#!/usr/bin/env python3
"""gen_inputs.py — ob-openbc raster generator (Stage-3 T2).

ob5: a 20 x 30 basin (dx = dy = 1 m) shoaling from 0.5 m depth at the
south wall to 3 m at the north stage row (bed elevation -0.5 at j = 0 to
-3.0 at j = 29, linear), with an initial along-shore jet uu = 0.3 m/s in
the 5 rows adjacent to the stage row and an x-asymmetric eta seed (+1 mm
on the west half of the jet band) so the zero_gradient lane's free-energy
feedback has a seed to amplify (plan ob5). Rasters are flat lists
j*nx + i (docs/user-guide/input-data.md).
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
NX, NY = 20, 30


def write(name: str, values, comment: str) -> None:
    lines = [f"# {comment}"]
    vals = [f"{v:.6f}" for v in values]
    for n in range(0, len(vals), 10):
        lines.append(" ".join(vals[n:n + 10]))
    (HERE / name).write_text("\n".join(lines) + "\n")


bed = [-0.5 - 2.5 * j / (NY - 1) for j in range(NY) for _ in range(NX)]
write("ob5_dem.dat", bed, "ob5 bed elevation [m], j*nx+i: -0.5 (south) -> -3.0 (north)")

uu0 = [0.3 if j >= NY - 5 else 0.0 for j in range(NY) for _ in range(NX)]
write("ob5_uu0.dat", uu0, "ob5 initial uu [m/s]: 0.3 in the 5 rows next to the stage row")

eta0 = [0.001 if (j >= NY - 5 and i < NX // 2) else 0.0
        for j in range(NY) for i in range(NX)]
write("ob5_eta0.dat", eta0, "ob5 initial eta [m]: +1 mm west-half seed in the jet band")

print("ob5 rasters written")
