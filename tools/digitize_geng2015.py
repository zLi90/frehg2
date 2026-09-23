#!/usr/bin/env python3
"""digitize_geng2015.py — extract the g5 reference curves from the Geng &
Boufadel (2015) PDF (v2 plan §6.4).

The gate figures turned out to be *vector line art*, not raster scans, so
"digitization" here is exact path extraction plus axis calibration — no
pixel picking:

  Fig. 3 (p. 431): the evaporation-rate curve is a stroked polyline of
    18,033 vertices in three chained subpaths (x strictly monotone) — the
    plotted data itself.
  Fig. 4 (p. 431): each moisture profile is a chain of 26 cubic Beziers,
    uniquely identified by its (dash pattern, gray level, stroke width)
    triple, which maps 1:1 onto the legend. Five of the six legend curves
    exist as vector paths (0/10/20/30/50 h); the 40 h curve is absent from
    the PDF's vector content (not needed: g5 gates 20 h and 50 h).
  Fig. 9b (p. 433): same construction; the base-case (K0 = 9e-4 m/s)
    salinity profiles are the two solid curves (black = 20 h, gray = 50 h).

Axis calibration is a linear map fitted to the printed tick-label centers;
the fit residual against every label is reported and asserted < 0.6 pt.
Self-checks (plan §6.3 reference sanity):
  - the Fig. 4 "0 h" curve must land at moisture ratio 1.000 +/- 0.01
    (the initial condition is exactly saturated);
  - the digitized E(0) must be within 5 % of the closed-form Table-1
    bulk-aerodynamic value 1.470e-7 m/s (V2-A13's derivation);
  - every profile's deep tail must approach its initial value.

Outputs (committed under benchmarks/g5-geng2015/reference/):
  fig3_evaporation_rate.csv   time_h, evaporation_m_per_s
  fig4_moisture_ratio.csv     elevation_m, moisture_{0,10,20,30,50}h
  fig9b_salinity.csv          elevation_m, salinity_{20,50}h_gL
  overlay_fig3.png, overlay_fig4.png, overlay_fig9b.png
                              digitized points re-plotted over the original
                              figure rendering (the §6.4 owner check)

Requires PyMuPDF and matplotlib; the paper PDF is NOT committed (see
--pdf). CI never runs this script — it is the reproducibility record.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import fitz  # PyMuPDF
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "benchmarks" / "g5-geng2015" / "reference"

# Closed-form E(0) from Table 1 via the plan §3.2 chain (V2-A13):
# Tetens e_sat(20C) = 2.3383 kPa -> q_sat = 0.014480 -> q_a = 0.2 q_sat;
# R_air = 94.909 s/m (Liu, U = 1 m/s); rho_a = 1.2041 kg/m^3 (ideal gas).
E0_CLOSED_FORM = 1.470e-7  # m/s


def label_centers(page: fitz.Page, region: fitz.Rect) -> list[tuple[float, float, float]]:
    """(value, x_center, y_center) of numeric text spans inside region."""
    out = []
    for block in page.get_text("dict")["blocks"]:
        if block["type"] != 0:
            continue
        for line in block["lines"]:
            for span in line["spans"]:
                x0, y0, x1, y1 = span["bbox"]
                if not region.contains(fitz.Point(0.5 * (x0 + x1), 0.5 * (y0 + y1))):
                    continue
                try:
                    value = float(span["text"].strip())
                except ValueError:
                    continue
                out.append((value, 0.5 * (x0 + x1), 0.5 * (y0 + y1)))
    return out


def linear_map(labels: list[tuple[float, float]], what: str) -> tuple[float, float]:
    """Least-squares data = a * page + b over (page_coord, data_value)."""
    page_c = np.array([p for p, _ in labels])
    data_v = np.array([v for _, v in labels])
    a, b = np.polyfit(page_c, data_v, 1)
    residual_pt = np.max(np.abs((data_v - b) / a - page_c))
    print(f"  {what}: {len(labels)} labels, max residual {residual_pt:.3f} pt")
    assert residual_pt < 0.6, f"{what}: axis labels are not collinear ({residual_pt:.2f} pt)"
    return float(a), float(b)


def bezier_points(items, samples: int = 25) -> np.ndarray:
    """Densely sample a drawing's cubic-Bezier ('c') and line ('l') items."""
    pts: list[tuple[float, float]] = []
    ts = np.linspace(0.0, 1.0, samples)
    for it in items:
        if it[0] == "c":
            p0, p1, p2, p3 = (np.array([p.x, p.y]) for p in it[1:5])
            for t in ts:
                q = ((1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 +
                     3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3)
                pts.append((q[0], q[1]))
        elif it[0] == "l":
            pts.append((it[1].x, it[1].y))
            pts.append((it[2].x, it[2].y))
    return np.array(pts)


def find_curve(page: fitz.Page, n_items: int, dash: str, gray: float, width: float,
               y_range: tuple[float, float]) -> np.ndarray:
    """Locate the unique stroked path with this legend signature."""
    hits = []
    for d in page.get_drawings():
        if d["type"] != "s" or len(d["items"]) != n_items:
            continue
        r = d["rect"]
        if not (y_range[0] <= r.y0 <= y_range[1]):
            continue
        if d.get("dashes") != dash:
            continue
        g = d["color"][0] if d.get("color") else -1.0
        if abs(g - gray) > 0.05 or abs((d.get("width") or 0.0) - width) > 0.05:
            continue
        hits.append(d)
    assert len(hits) == 1, f"signature ({dash!r}, {gray}, {width}) matched {len(hits)} paths"
    return bezier_points(hits[0]["items"])


def profile_to_grid(pts_page: np.ndarray, x_map, y_map, z_grid: np.ndarray,
                    clip: fitz.Rect) -> np.ndarray:
    """Curve -> value(z) on z_grid: clip to the axes box, sort by elevation,
    interpolate. Values outside the curve's z span are NaN."""
    ax, bx = x_map
    ay, by = y_map
    inside = ((pts_page[:, 0] >= clip.x0 - 0.5) & (pts_page[:, 0] <= clip.x1 + 0.5) &
              (pts_page[:, 1] >= clip.y0 - 0.5) & (pts_page[:, 1] <= clip.y1 + 0.5))
    pts = pts_page[inside]
    z = ay * pts[:, 1] + by
    v = ax * pts[:, 0] + bx
    order = np.argsort(z)
    z, v = z[order], v[order]
    # Collapse duplicate z (vertical runs) to their mean.
    zu, idx = np.unique(np.round(z, 6), return_inverse=True)
    vu = np.array([v[idx == k].mean() for k in range(zu.size)])
    out = np.interp(z_grid, zu, vu, left=np.nan, right=np.nan)
    out[(z_grid < zu[0]) | (z_grid > zu[-1])] = np.nan
    return out


def overlay(page: fitz.Page, clip: fitz.Rect, x_map, y_map, series, path: Path,
            xlabel: str, ylabel: str, logy: bool = False) -> None:
    """Digitized points re-plotted over the original figure rendering."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    zoom = 4
    pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), clip=clip)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    ax_, bx = x_map
    ay, by = y_map
    extent = [ax_ * clip.x0 + bx, ax_ * clip.x1 + bx,
              ay * clip.y1 + by, ay * clip.y0 + by]
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.imshow(img, extent=extent, aspect="auto", zorder=0)
    for label, xs, ys in series:
        ax.plot(xs, ys, "o", markersize=2.5, markerfacecolor="none",
                markeredgewidth=0.8, label=label, zorder=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if logy:
        pass  # figure axis is linear; keep the overlay linear to match
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"  wrote {path.name}")


def write_csv(path: Path, header: list[str], columns: list[np.ndarray],
              comment: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write(comment)
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in zip(*columns):
            writer.writerow(["" if (isinstance(v, float) and np.isnan(v)) else f"{v:.6g}"
                             for v in row])
    print(f"  wrote {path.name} ({columns[0].size} rows)")


COMMENT = ("# Digitized from Geng & Boufadel (2015), J. Hydrology 524:427-438,\n"
           "# by tools/digitize_geng2015.py (exact vector-path extraction; see\n"
           "# benchmarks/g5-geng2015/reference/DIGITIZATION.md). Code-to-code\n"
           "# reference (MARUN), not physical truth (v2 plan #3.4).\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path,
                        default=HERE.parent.parent / "references" /
                        "Geng_and_Boufadel_JH2015.pdf")
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    doc = fitz.open(args.pdf)

    ok = True

    # --- Fig. 3 (page index 4): evaporation rate vs time -------------------
    print("Fig. 3 — evaporation rate:")
    p431 = doc[4]
    xt = [(v, x) for v, x, y in label_centers(p431, fitz.Rect(80, 346, 270, 356))
          if v in (0, 10, 20, 30, 40, 50)]
    yt = [(v, y) for v, x, y in label_centers(p431, fitz.Rect(55, 255, 86, 350))]
    x_map3 = linear_map([(x, v) for v, x in xt], "fig3 x (hours)")
    y_map3 = linear_map([(y, v) for v, y in yt], "fig3 y (m/s)")
    # The three chained solid polylines (615 + 16803 + 612 'l' items).
    chain = []
    for d in p431.get_drawings():
        if d["type"] == "s" and len(d["items"]) in (615, 16803, 612):
            chain.append(bezier_points(d["items"], samples=1))
    assert len(chain) == 3, f"fig3: expected 3 subpaths, found {len(chain)}"
    pts = np.vstack(sorted(chain, key=lambda a: a[:, 0].min()))
    t_h = x_map3[0] * pts[:, 0] + x_map3[1]
    e_ms = y_map3[0] * pts[:, 1] + y_map3[1]
    order = np.argsort(t_h)
    t_h, e_ms = t_h[order], e_ms[order]
    t_grid = np.round(np.arange(0.1, 50.0 + 1e-9, 0.1), 3)
    e_grid = np.interp(t_grid, t_h, e_ms)
    e0 = float(e_grid[0])
    print(f"  E(0.1h) = {e0:.4e} m/s (closed form {E0_CLOSED_FORM:.3e}; "
          f"ratio {e0 / E0_CLOSED_FORM:.3f})")
    if not 0.95 <= e0 / E0_CLOSED_FORM <= 1.05:
        print("  SELF-CHECK FAIL: E(0) disagrees with the Table-1 closed form")
        ok = False
    e10 = float(np.interp(10.0, t_grid, e_grid))
    e50 = float(np.interp(50.0, t_grid, e_grid))
    print(f"  E(10h) = {e10:.4e}   E(50h) = {e50:.4e}   E50/E0 = {e50 / e0:.4e}")
    write_csv(args.out / "fig3_evaporation_rate.csv",
              ["time_h", "evaporation_m_per_s"], [t_grid, e_grid], COMMENT)
    clip3 = fitz.Rect(50, 255, 290, 370)
    overlay(p431, clip3, x_map3, y_map3,
            [("digitized ER", t_grid[::5], e_grid[::5])],
            args.out / "overlay_fig3.png", "time [h]", "evaporation flux [m/s]")

    # --- Fig. 4 (page index 4): moisture-ratio profiles --------------------
    print("Fig. 4 — moisture-ratio profiles:")
    xt4 = label_centers(p431, fitz.Rect(350, 78, 545, 88))
    yt4 = label_centers(p431, fitz.Rect(340, 85, 356, 188))
    x_map4 = linear_map([(x, v) for v, x, y in xt4], "fig4 x (moisture)")
    y_map4 = linear_map([(y, v) for v, x, y in yt4], "fig4 y (elevation)")
    clip4 = fitz.Rect(358.598, 90.863, 538.947, 182.818)  # the axes clip box
    # Legend signatures: (dash, gray, width). The 40 h curve has no vector
    # path in the PDF (checked exhaustively over all 26-item strokes).
    sig4 = {
        "0h": ("[] 0", 0.5758, 1.04),
        "10h": ("[] 0", 0.1367, 0.78),
        "20h": ("[ 2.3404 3.1205 0 3.1205 ] 0", 0.1367, 0.78),
        "30h": ("[ 0 3.7707 ] 0", 0.1367, 1.885),
        "50h": ("[ 2.0804 2.0804 ] 0", 0.7777, 1.04),
    }
    z_grid = np.round(np.arange(1.60, 2.00 + 1e-9, 0.005), 4)
    fig4 = {}
    for name, (dash, gray, width) in sig4.items():
        pts_c = find_curve(p431, 26, dash, gray, width, (85.0, 95.0))
        fig4[name] = profile_to_grid(pts_c, x_map4, y_map4, z_grid, clip4)
    # Self-check: the 0 h profile is the saturated initial condition.
    m0 = np.nanmean(fig4["0h"])
    print(f"  0h profile mean moisture = {m0:.4f} (must be 1.000 +/- 0.01)")
    if not 0.99 <= m0 <= 1.01:
        print("  SELF-CHECK FAIL: 0h calibration curve is off saturation")
        ok = False
    for name in ("10h", "20h", "30h", "50h"):
        deep = np.nanmean(fig4[name][z_grid <= 1.70])
        surf = fig4[name][~np.isnan(fig4[name])][-1]
        print(f"  {name:>4}: surface {surf:.3f}, deep tail {deep:.3f}")
        if abs(deep - 1.0) > 0.02:
            print(f"  SELF-CHECK FAIL: {name} deep tail should stay saturated")
            ok = False
    write_csv(args.out / "fig4_moisture_ratio.csv",
              ["elevation_m"] + [f"moisture_{n}" for n in sig4],
              [z_grid] + [fig4[n] for n in sig4],
              COMMENT + "# The 40 h legend curve has no vector path in the "
              "PDF; g5 gates 20 h and 50 h.\n")
    overlay(p431, fitz.Rect(330, 62, 545, 215), x_map4, y_map4,
            [(n, fig4[n][::4], z_grid[::4]) for n in sig4],
            args.out / "overlay_fig4.png", "moisture ratio", "elevation [m]")

    # --- Fig. 9b (page index 6): salinity profiles, base case --------------
    print("Fig. 9b — salinity profiles (K0 = 9e-4 m/s):")
    p433 = doc[6]
    # Fig. 9b's salinity axis sits on TOP of its panel (labels y ~ 372-379);
    # the elevation labels run down its left edge (y ~ 380-486).
    xt9 = [(v, x) for v, x, y in label_centers(p433, fitz.Rect(340, 370, 552, 380))
           if v in (20, 40, 60, 80, 100, 120)]
    yt9 = [(v, y) for v, x, y in label_centers(p433, fitz.Rect(330, 378, 348, 488))
           if 1.5 <= v <= 2.05]
    x_map9 = linear_map([(x, v) for v, x in xt9], "fig9b x (g/L)")
    y_map9 = linear_map([(y, v) for v, y in yt9], "fig9b y (elevation)")
    # Axes box in page coords from the calibrated maps.
    x0p = (20 - x_map9[1]) / x_map9[0]
    x1p = (120 - x_map9[1]) / x_map9[0]
    y0p = (2.0 - y_map9[1]) / y_map9[0]
    y1p = (1.6 - y_map9[1]) / y_map9[0]
    clip9 = fitz.Rect(min(x0p, x1p), min(y0p, y1p), max(x0p, x1p), max(y0p, y1p))
    print(f"  axes box {clip9}")
    # Base case = the solid curves in the lower (Fig. 9b) group.
    fig9 = {}
    for name, gray in (("20h", 0.1367), ("50h", 0.7777)):
        pts_c = find_curve(p433, 26, "[] 0", gray, 0.839, (clip9.y0 - 5, clip9.y1))
        fig9[name] = profile_to_grid(pts_c, x_map9, y_map9, z_grid, clip9)
    for name in ("20h", "50h"):
        deep = np.nanmean(fig9[name][z_grid <= 1.70])
        top = np.nanmax(fig9[name])
        print(f"  {name:>4}: deep tail {deep:.1f} g/L (initial 25), peak {top:.1f} g/L")
        if abs(deep - 25.0) > 2.0:
            print(f"  SELF-CHECK FAIL: {name} deep tail should hold the initial 25 g/L")
            ok = False
    if not (np.nanmax(fig9["50h"]) > 60.0):
        print("  SELF-CHECK FAIL: the 50h near-surface peak must exceed 60 g/L "
              "(the g5(ii) criterion reads it off this curve)")
        ok = False
    write_csv(args.out / "fig9b_salinity.csv",
              ["elevation_m", "salinity_20h_gL", "salinity_50h_gL"],
              [z_grid, fig9["20h"], fig9["50h"]],
              COMMENT + "# Base case K0 = 9e-4 m/s (the solid legend pair).\n")
    overlay(p433, fitz.Rect(clip9.x0 - 30, clip9.y0 - 12, clip9.x1 + 12, clip9.y1 + 40),
            x_map9, y_map9,
            [(n, fig9[n][::4], z_grid[::4]) for n in fig9],
            args.out / "overlay_fig9b.png", "salinity [g/L]", "elevation [m]")

    print("digitize_geng2015:", "OK" if ok else "SELF-CHECK FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
