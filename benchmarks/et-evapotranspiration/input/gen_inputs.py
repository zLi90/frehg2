#!/usr/bin/env python3
"""gen_inputs.py — et-evapotranspiration raster/series generator (Stage-3 T5).

et9: the tidal-flat bed b(y) = 0.40 - 0.02*y on a 4 x 24 (dx = dy = 2 m)
periodic-x strip (y_c = 1..47 m: 0.38 down to -0.54 m), uniform in x so
the roll/rank lanes stay exact. et8b: the rain pulse series (2 mm over
600 s starting at t = 7200 s). et-x: the 8 x 8 orient bed (corner pond
basin: a plane sloping toward the south-west corner, D4-symmetric inputs
NOT required — the battery transforms the raster). Rasters are flat
lists j*nx + i (docs/user-guide/input-data.md).
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent


def write(name, lines):
    (HERE / name).write_text("\n".join(lines) + "\n")


# et9 bed: 4 x 24, b = 0.40 - 0.02 * y_center (x-uniform).
rows = ["# et9 tidal-flat bed [m], j*nx+i: 0.40 - 0.02*y, x-uniform"]
for j in range(24):
    y = (j + 0.5) * 2.0
    b = 0.40 - 0.02 * y
    rows.append(" ".join(f"{b:.8f}" for _ in range(4)))
write("et9_dem.dat", rows)

# et8b rain pulse: 2 mm over 600 s from t = 7200 s (3.3333e-6 m/s).
rows = ["# et8b rain [m/s]: 2 mm pulse over 600 s at t = 7200 s"]
for t, v in ((0, 0.0), (7199, 0.0), (7200, 3.3333333333e-6),
             (7799, 3.3333333333e-6), (7800, 0.0), (10800, 0.0)):
    rows.append(f"{t} {v:.10e}")
write("et8b_rain.dat", rows)

# et-x orient bed: 8 x 8 (dx = dy = 2 m), plane sloping to the SW corner
# with the corner cells deepest (pond forms there from a flat eta IC).
# STEEP steps (0.1 m per index, S3-A8): with the eta IC at -0.10 every
# cell's wet/dry margin is >= 0.05 m - 5x the documented swe mixed
# asymmetry class (1e-2), so no mixed transform can flip a shoreline
# cell (the 0.02-step first authoring had 0.01 margins and the tracer
# flipped at the shoreline).
rows = ["# et-x orient bed [m], j*nx+i: steep plane to the SW corner"]
for j in range(8):
    vals = []
    for i in range(8):
        x = (i + 0.5) * 2.0
        y = (j + 0.5) * 2.0
        vals.append(f"{-0.75 + 0.05 * (x + y):.8f}")
    rows.append(" ".join(vals))
write("etx_dem.dat", rows)

print("et rasters written")
