#!/usr/bin/env python3
"""Generate gridded inputs for the swere-lateral-hillslope MPI strong-scaling case.

This case is a deliberate *extension* of validation/swere-lateral-hillslope
(SERGHEI case4) into a large, homogeneous 3D domain suitable for MPI strong
scaling. Three transforms relative to the parent case:

  1. x refined 100 m -> 10 m       : nx 40 -> 400  (physical width 4000 m kept)
  2. x-z plane duplicated in y      : ny 1 -> 800, dy = 10 m (homogeneous in y)
  3. run length 5 yr -> 1 hour      : set in the YAML (time.t_end), not here

Total grid = 400 x 800 x 60 = 19 200 000 cells.

Every field written here is an ANALYTIC continuation of the parent case, so the
refined + widened domain is the *same physical hillslope*, just discretised
finer and replicated in y (no y-gradient -> a clean, load-balanced scaling
workload):

  * Bed (dem.dat): the parent DEM is the line elev = 3.95 - 0.001*x_c
    (parent centres x_c = 50..3950 m give exactly 3.9..0.0). Evaluated at the
    refined cell centres.
  * Water table (water_table.dat): the parent initial head is hydrostatic about
    a water table that runs linearly from the two lateral head BCs, -4.1 m at
    the left (x_c=50) to -14.1 m at the right (x_c=3950). We hand frehg2 that
    water table directly (initial_conditions.groundwater.water_table) instead of
    a 300 MB 3D head file: RichardsSolver builds h = (WT + offset) - z_cell, the
    datum offset cancels, and the saturated-zone head matches the parent case
    exactly (the parent's -1.25 m unsaturated clamp is replaced by a clean
    hydrostatic profile above the table -- an initial condition the solver
    relaxes and which is, if anything, gentler for a benchmark start).

Both rasters are y-invariant: the same nx-long row is repeated ny times, in the
reader's (j*nx + i) file order (j = y-row, i = x-col; no ESRI north-up flip,
matching the parent case).

Rewrite the inputs in place:  python3 gen_inputs.py
"""

NX, NY, NZ = 400, 800, 60
DX, DY, DZ = 10.0, 10.0, 0.25

# Parent-case first-period rainfall rate [m/s] (rainfall 0.017 mm/h -> m/s).
# The 1-hour run stays entirely within this first period, so rain is constant.
RAIN_RATE = 4.722222222e-09


def x_center(i):
    """Cell-centre x of column i [m]."""
    return (i + 0.5) * DX


def bed(i):
    """Bed elevation at column i [m]  (parent line 3.95 - 0.001*x_c)."""
    return 3.95 - 0.001 * x_center(i)


def water_table(i):
    """Initial water-table elevation at column i [m].

    Linear through the parent lateral head BCs: -4.1 m at x_c=50 to
    -14.1 m at x_c=3950  ->  WT = -4.1 - (x_c - 50)/390.
    """
    return -4.1 - (x_center(i) - 50.0) / 390.0


def write_raster(path, values_per_column, cellsize):
    """Write an ESRI-header raster (nx*ny, y-invariant) in (j*nx+i) order.

    values_per_column: list of NX floats (the x-profile); repeated for every
    y-row because the domain is homogeneous in y.
    """
    row = " ".join("{:.4f}".format(v) for v in values_per_column)
    with open(path, "w") as f:
        f.write("ncols {}\n".format(NX))
        f.write("nrows {}\n".format(NY))
        f.write("xllcorner 0.0\n")
        f.write("yllcorner 0.0\n")
        f.write("cellsize {:.1f}\n".format(cellsize))
        f.write("nodata_value -9999\n")
        for _ in range(NY):
            f.write(row)
            f.write("\n")


def write_rain(path):
    with open(path, "w") as f:
        f.write("# swere-lateral-hillslope-scaling: constant rainfall [m/s]\n")
        f.write("# = SERGHEI case4 first-period rate; the 1-hour run stays within it\n")
        f.write("0.0 {:.9e}\n".format(RAIN_RATE))
        f.write("86400.0 {:.9e}\n".format(RAIN_RATE))
        f.write("90000.0 {:.9e}\n".format(RAIN_RATE))


def main():
    import os
    here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "input")
    os.makedirs(here, exist_ok=True)

    bed_profile = [bed(i) for i in range(NX)]
    wt_profile = [water_table(i) for i in range(NX)]

    write_raster(os.path.join(here, "dem.dat"), bed_profile, DX)
    write_raster(os.path.join(here, "water_table.dat"), wt_profile, DX)
    write_rain(os.path.join(here, "rain.dat"))

    print("wrote input/dem.dat, input/water_table.dat, input/rain.dat")
    print("  grid            : {} x {} x {} = {} cells".format(
        NX, NY, NZ, NX * NY * NZ))
    print("  bed  (x-profile): {:.4f} .. {:.4f} m".format(bed_profile[0], bed_profile[-1]))
    print("  WT   (x-profile): {:.4f} .. {:.4f} m".format(wt_profile[0], wt_profile[-1]))
    print("  raster values   : {} each (dem, water_table)".format(NX * NY))


if __name__ == "__main__":
    main()
