#!/usr/bin/env python3
"""gen_hrl_ic.py — regenerate benchmarks/g8-hrl/input/ic_temperature.dat
(the g8 seeded initial temperature; v2 plan §6.2 committed-generator
rule).

Field: the conduction profile between the pinned-cell centres plus a
0.5 K seed of the exactly-critical admissible mode
    dT(x, z) = 0.5 cos(pi x / H) sin(pi zeta)
with zeta the fractional depth between the pinned planes (the m = 2
cosine mode of the 2H-wide box; V2-A17 — the earlier sine seed was
orthogonal to it). File order (j*nx + i)*nz + k, one value per line.
"""

from __future__ import annotations

import math
from pathlib import Path

NX, NY, NZ = 40, 1, 20
DX, DZ = 0.05, 0.05
H = 1.0
T_TOP, T_BOT = 20.0, 30.0
SEED = 0.5

Z_TOP_PIN = (0 + 0.5) * DZ            # depth of the pinned top-cell centre
Z_BOT_PIN = (NZ - 1 + 0.5) * DZ       # depth of the pinned bottom-cell centre
H_EFF = Z_BOT_PIN - Z_TOP_PIN         # 0.95 m


def temperature(i: int, k: int) -> float:
    x = (i + 0.5) * DX
    depth = (k + 0.5) * DZ
    zeta = (depth - Z_TOP_PIN) / H_EFF
    zeta = min(max(zeta, 0.0), 1.0)
    base = T_TOP + (T_BOT - T_TOP) * zeta
    return base + SEED * math.cos(math.pi * x / H) * math.sin(math.pi * zeta)


def main() -> None:
    out = Path(__file__).resolve().parent / "ic_temperature.dat"
    lines = ["# HRL initial temperature: conduction profile between the pinned-cell",
             "# centres + 0.5 K of the critical m = 2 mode cos(pi x/H) sin(pi zeta)",
             "# (j*nx+i)*nz+k order; generator: gen_hrl_ic.py (V2-A17)"]
    for j in range(NY):
        for i in range(NX):
            for k in range(NZ):
                lines.append(f"{temperature(i, k):.6f}")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {out} ({NX * NY * NZ} values, H_eff = {H_EFF:.2f} m)")


if __name__ == "__main__":
    main()
