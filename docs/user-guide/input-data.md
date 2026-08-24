# Input data files

Two plain-text formats back the `file:` references in the configuration.
Paths resolve against the configuration file's directory.

## Gridded rasters (bathymetry, roughness, initial fields)

Used by `domain.bottom_elevation.file`, `surface_water.friction.coefficient.file`,
`initial_conditions.surface.eta.file`, etc. Two forms are auto-detected by the
first token:

**Flat list** — bare whitespace-separated values, `#` comments allowed,
ordered `j*nx + i` (i fastest). For a 200×10 grid that is 2000 values, the
first `nx` values being row `j = 0`:

```
# 200 x 10 bed, ordered j*nx + i
0.20 0.19 0.18 ... (200 values for j=0)
0.20 0.19 0.18 ... (200 values for j=1)
...
```

**Headered raster** — `key value` header lines (`ncols`, `nrows` required;
`xllcorner`, `yllcorner`, `cellsize`, `nodata_value` optional), then
`nrows·ncols` values. The **first data row is `j = 0`** (this is the
legacy-faithful convention, *not* ESRI north-up). `nodata_value` cells are
read as NaN.

```
ncols 200
nrows 10
cellsize 0.109725
0.20 0.19 ...
...
```

The value count must match `nx·ny` exactly, or the run aborts with the
observed-vs-expected counts.

**3D fields** (`soil.map.file`, 3D initial conditions) are flat lists in
`(j·nx + i)·nz + k` order — the same flattening the HDF5 output uses, `k = 0`
at the top layer.

## Time series (rainfall, evaporation, wind, time-varying BCs)

One `time value` pair per line, whitespace-separated, `#` for comments. Times
are seconds and must be strictly increasing; values are SI. Evaluation is
piecewise-linear and **clamped** to the first/last value outside the sampled
range:

```
# rainfall [m/s]: 5.5e-6 for the first 12000 s, then dry
0       5.5e-6
12000   5.5e-6
12001   0
18000   0
```

Because series are evaluated statelessly by time, restart reproduces them
exactly.
