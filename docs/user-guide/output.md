# Output reference

Everything lands in the single HDF5 file named by `output.filename`. The
layout (all field datasets are `float64`, snapshot times are keyed by integer
seconds):

```
/frehg2                     root group; attrs: version, git sha, config text, created
/grid/x_center     [NX]     cell-center x coordinates [m]
/grid/y_center     [NY]     cell-center y coordinates [m]
/grid/z_center     [NZ]     layer-center depths
/grid/bottom       [NY*NX]  bed elevation, index j*NX + i
/grid/ktop         [NY*NX]  top active layer per column
/grid/dz           [NZ]     layer thicknesses
/surface/<var>/<t> [NY*NX]     surface snapshot at time t, index j*NX + i
/groundwater/<var>/<t> [NY*NX*NZ]  subsurface snapshot, index (j*NX + i)*NZ + k
/groundwater/zcell/0 [NY*NX*NZ]    per-cell layer-center elevations (static)
/transport/<var>/<t>         scalar snapshots (concentration on either grid)
/monitor/<name>    (rows, 1+k)  extendable table: time + the monitored variables
/monitor/mass_audit (rows, 7+)  surface volume budget (see below)
/monitor/gw_mass_audit (rows, 7) subsurface volume budget (see below)
/monitor/transport_audit     scalar-mass budget (transport runs; see below)
/checkpoint/<t>/...          full prognostic state (only if checkpointing)
```

Surface fields (`eta`, `depth`, `uu`, `vv`, `seepage`) are 2D arrays of
`NY*NX` values flattened as `j*NX + i`. Groundwater fields
(`hydraulic_head`, `water_content`, `qx`, `qy`, `qz`) are 3D arrays
flattened as `(j*NX + i)*NZ + k`, `k = 0` at the top layer; cells outside
the active subsurface are NaN. Transport writes `concentration` (the
subsurface grid) and `concentration_surface` (the surface grid). Each
dataset carries `units`, `long_name`, and `time` attributes. `eta` is
absolute free-surface elevation; `depth` is `eta − bottom`; `seepage`
(coupled runs) is the surface-applied exchange rate [m/s], positive upward.

**Monitors.** `/monitor/<name>` is a growing table with one row per time
step. Its `columns` attribute is the CSV header, e.g. `time,depth,vv`, and
`i`/`j` attributes record the sampled cell.

**Mass audits.** `/monitor/mass_audit` (surface) is written whenever the
SWE module runs. Its columns (`columns` attribute:
`time,volume,rain,evaporation,boundary_outflow,bc_inflow[,seepage],clamped`
— `seepage` only in coupled runs) are the running water-volume budget,
reduced across all ranks. `clamped` measures the volume the legacy
below-bed clamp creates (a scheme defect reported as data; ~1.6 % of rain
in the b5 regime, effectively zero elsewhere). The table is the
mass-consistent source of boundary discharge — prefer it over
`uu·depth` when you need conserved outflow, because stored velocities carry a
legacy scaling that under-reads instantaneous flux.

`/monitor/gw_mass_audit` (subsurface) is written whenever the groundwater
module runs, with columns
`time,volume,boundary_in,ss_storage,realloc,realloc_dropped,vloss` — the
subsurface identity closes to rounding
(`Δvolume = boundary_in − ss_storage + realloc − vloss`). In coupled runs
the exchanged volume appears symmetrically: positive `seepage` in the
surface table is water the subsurface gave up.

`/monitor/transport_audit` (transport runs) is the scalar-mass budget:
every non-conservative piece of the legacy scheme is measured into its own
cumulative column (exchange, sources, boundary leaks, limiter/bounds clips,
and the ledger re-anchor terms); the closure identities hold to rounding in
every regime. The term-by-term definition is in the
[transport theory page](../theory/transport.md).

## Inspecting the file

The HDF5 command-line tools show structure and values:

```bash
h5ls -r out/output.h5                       # full tree
h5dump -d /surface/depth/1800 out/output.h5 # one snapshot
h5dump -a /monitor/outlet_q/columns out/output.h5
```

## Plotting with Python

Install `h5py`, `numpy`, `matplotlib` (`pip install h5py numpy matplotlib`).

A hydrograph from a monitor table:

```python
import h5py, matplotlib.pyplot as plt

def columns(dset):                      # the 'columns' attr is fixed-ASCII bytes
    raw = dset.attrs["columns"]
    return (raw.decode() if isinstance(raw, bytes) else raw).split(",")

with h5py.File("out/output.h5", "r") as f:
    d = f["/monitor/outlet_q"]          # columns: time, depth, uu
    cols = columns(d)
    t = d[:, 0]
    depth = d[:, cols.index("depth")]

plt.plot(t, depth)
plt.xlabel("time [s]"); plt.ylabel("depth [m]")
plt.title("outlet depth"); plt.savefig("hydrograph.png")
```

A 2D field snapshot as a map (reshape the flat `j*NX + i` array to `(NY, NX)`):

```python
import h5py, numpy as np, matplotlib.pyplot as plt

with h5py.File("out/output.h5", "r") as f:
    nx = f["/grid/x_center"].shape[0]
    ny = f["/grid/y_center"].shape[0]
    x = f["/grid/x_center"][:]; y = f["/grid/y_center"][:]
    depth = f["/surface/depth/1800"][:].reshape(ny, nx)  # j (rows), i (cols)

plt.pcolormesh(x, y, depth, shading="nearest")
plt.colorbar(label="depth [m]"); plt.xlabel("x [m]"); plt.ylabel("y [m]")
plt.gca().set_aspect("equal"); plt.savefig("depth_1800.png")
```

The mass balance over the run (should stay closed):

```python
import h5py
with h5py.File("out/output.h5", "r") as f:
    a = f["/monitor/mass_audit"][:]   # columns per its 'columns' attribute
t, vol, rain, evap, outflow, inflow = a.T[:6]
clamped = a.T[-1]                     # last column; 'seepage' precedes it in coupled runs
seep = a.T[6] if a.shape[1] == 8 else 0.0
residual = vol - (rain - evap - outflow + inflow + seep + clamped) - vol[0]
print("max mass-balance residual:", abs(residual).max())
```

For interactive exploration or 3D rendering, the file also opens directly in
**ParaView** and **VisIt** via their HDF5 readers, and in
[panoply](https://www.giss.nasa.gov/tools/panoply/) or any `h5py`/`xarray`
workflow.
