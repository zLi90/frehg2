#!/usr/bin/env python3
"""Generate the committed x-orient battery inputs (v2 plan §8.1).

- swe_dem.dat: 10 x 10 flat-list bed (j*nx + i order, one value per line) —
  a plane sloping up to the east/north with an off-axis Gaussian mound, so
  the surface is injective under every dihedral transform (any orientation
  mix-up moves the pool and the drainage path).
- gw_soil_id.dat: 12 x 12 x 8 flat-list soil ids ((j*nx + i)*nz + k order) —
  two soils layered in z only (top half 0, bottom half 1), which is
  invariant under every in-plane transform by construction (the battery
  transforms polygons and the DEM; a z-only layering needs no transform).
- tr_eta0.dat / tr_blob.dat: 12 x 12 flat-list initial stage tilt and
  tracer blob for the transport battery — a walled sloshing basin whose
  flow advects an off-axis blob, gating the advection scheme's own
  orientation symmetry with no boundary throughflow (the BC side sweep
  lives in the outflow-staircase and surface-bc-sides gates instead).
- cs_dem.dat: 8 x 8 flat-list bed for the coupled-salinity decomposition
  lane (coupled-salt-base.yaml) — a plane rising 0.1 m per cell east and
  0.03 m per cell north from -0.3 m, so a 0.15 m stage floods the west
  half, leaves the east half dry, and drives lateral subsurface flow across
  both rank interfaces of a 2 x 2 decomposition.

The battery harness (tests/regression/run_regression.py, the *-orient
subcommands) transforms swe_dem.dat with the same dihedral table it applies
to the BC polygons; regenerating these files must reproduce them exactly
(they are committed so the case is runnable stand-alone).
"""

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def swe_dem() -> None:
    nx = ny = 10
    dx = dy = 8.0
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(x, y)  # (ny, nx)
    z = (0.004 * xx + 0.003 * yy
         + 0.15 * np.exp(-(((xx - 52.0) ** 2) + ((yy - 28.0) ** 2))
                         / (2.0 * 18.0 ** 2)))
    with open(HERE / "swe_dem.dat", "w", encoding="utf-8") as handle:
        handle.write("# x-orient SWE battery bed: tilted plane + off-axis "
                     "mound (gen_inputs.py)\n")
        for j in range(ny):
            for i in range(nx):
                handle.write(f"{z[j, i]:.10f}\n")


def gw_soil_id() -> None:
    nx = ny = 12
    nz = 8
    with open(HERE / "gw_soil_id.dat", "w", encoding="utf-8") as handle:
        handle.write("# x-orient gw battery soil ids: z-layered only "
                     "(k < 4 top soil 0, else soil 1) (gen_inputs.py)\n")
        for _j in range(ny):
            for _i in range(nx):
                for k in range(nz):
                    handle.write(f"{0 if k < nz // 2 else 1}\n")


def write2d(name: str, comment: str, field: np.ndarray) -> None:
    with open(HERE / name, "w", encoding="utf-8") as handle:
        handle.write(f"# {comment} (gen_inputs.py)\n")
        for j in range(field.shape[0]):
            for i in range(field.shape[1]):
                handle.write(f"{field[j, i]:.12f}\n")


def transport_ics() -> None:
    nx = ny = 12
    dx = dy = 2.0
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    xx, yy = np.meshgrid(x, y)  # (ny, nx)
    write2d("tr_eta0.dat", "x-orient transport battery initial stage tilt",
            0.008 * (xx - 12.0) / 12.0)
    write2d("tr_blob.dat", "x-orient transport battery tracer blob",
            np.exp(-(((xx - 7.0) ** 2) + ((yy - 15.0) ** 2)) / 8.0))


def coupled_salt_dem() -> None:
    nx = ny = 8
    ii, jj = np.meshgrid(np.arange(nx), np.arange(ny))  # (ny, nx)
    write2d("cs_dem.dat", "x-orient coupled-salinity lane bed: tilted plane",
            -0.3 + 0.1 * ii + 0.03 * jj)


if __name__ == "__main__":
    swe_dem()
    gw_soil_id()
    transport_ics()
    coupled_salt_dem()
    print("wrote swe_dem.dat, gw_soil_id.dat, tr_eta0.dat, tr_blob.dat, cs_dem.dat")
