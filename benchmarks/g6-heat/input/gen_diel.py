#!/usr/bin/env python3
"""gen_diel.py — regenerate input/diel_temperature.dat (the g6(c) diel
surface-temperature series: 20 + 5 sin(2 pi t / 86400) C at 600 s
sampling over 5 days; committed-generator rule, v2 plan §6.2)."""

import math
from pathlib import Path

lines = ["# g6(c) diel surface temperature: 20 + 5 sin(2 pi t / 86400) [C], 600 s sampling"]
for n in range(721):
    t = 600.0 * n
    lines.append(f"{t:.1f} {20.0 + 5.0 * math.sin(2.0 * math.pi * t / 86400.0):.6f}")
out = Path(__file__).resolve().parent / "diel_temperature.dat"
out.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"wrote {out}")
