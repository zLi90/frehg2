#!/usr/bin/env python3
"""gen_orient_ic.py — regenerate input/orient_ic.dat (the heat-battery
base IC: 20 C + an off-axis warm blob decaying with depth; the harness
regenerates transformed copies per dihedral variant). Committed-generator
rule, v2 plan §6.2; must match heat_orient_ic() in run_regression.py."""

import math
from pathlib import Path

NX = NY = 8
NZ = 6
DX = DY = 0.5
DZ = 0.25

lines = ["# heat orientation battery base IC (generator: gen_orient_ic.py)"]
for j in range(NY):
    y = (j + 0.5) * DY
    for i in range(NX):
        x = (i + 0.5) * DX
        blob = 6.0 * math.exp(-(((x - 1.1) ** 2 + (y - 1.7) ** 2) / 0.5 ** 2))
        for k in range(NZ):
            depth = (k + 0.5) * DZ
            profile = 1.0 - depth / (NZ * DZ)
            lines.append(f"{20.0 + blob * profile:.10f}")
out = Path(__file__).resolve().parent / "orient_ic.dat"
out.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"wrote {out}")
