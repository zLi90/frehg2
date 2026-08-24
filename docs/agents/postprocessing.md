# Post-processing guide for AI assistants

**Purpose:** everything an LLM needs to write Python scripts that read
Frehg2 results and make plots. The output is one HDF5 file per run whose
layout is a frozen contract (it is pinned by the validation gates and
cannot change silently), so everything below is stable ground truth.

Required packages: `h5py`, `numpy`, `matplotlib`
(`pip install h5py numpy matplotlib`). The file also opens directly in
ParaView/VisIt via their HDF5 readers for interactive 3D work.

## 1. File layout

One file per run at the configured `output.filename` (conventionally
`out/output.h5`). All field datasets are `float64` 1-D arrays; snapshot
times are group keys of **integer seconds** (`str(int(t))`).

```
/frehg2                       root group; attrs: version, git sha, full config text, created
/grid/x_center      [NX]      cell-center x [m]
/grid/y_center      [NY]      cell-center y [m]
/grid/z_center      [NZ]      nominal layer-center depths
/grid/bottom        [NY*NX]   bed elevation [m], index j*NX + i
/grid/ktop          [NY*NX]   first active layer per column (nz = fully inactive)
/grid/dz            [NZ]      layer thicknesses [m]
/surface/<var>/<t>      [NY*NX]      surface snapshot, index j*NX + i
/groundwater/<var>/<t>  [NY*NX*NZ]   subsurface snapshot, index (j*NX + i)*NZ + k
/groundwater/zcell/0    [NY*NX*NZ]   per-cell layer-center elevations (static; use for
                                     terrain-following/partial-cell vertical coordinates)
/transport/<var>/<t>                 scalar snapshots (either grid; same flattening)
/monitor/<name>         (rows, 1+k)  per-time-step table; attrs: columns (CSV), i, j
/monitor/mass_audit                  surface volume budget      (see §5)
/monitor/gw_mass_audit               subsurface volume budget   (see §5)
/monitor/transport_audit             scalar-mass budget         (see §5)
/checkpoint/<t>/...                  restart state (ignore for plotting)
```

Variables and units:

| dataset | meaning | units |
|---|---|---|
| `surface/eta` | absolute free-surface elevation | m |
| `surface/depth` | water depth = eta − bottom (0 when dry) | m |
| `surface/uu`, `vv` | x / y velocity | m/s |
| `surface/seepage` | surface-applied exchange rate, **positive upward** (coupled runs only) | m/s |
| `groundwater/hydraulic_head` | pressure head | m |
| `groundwater/water_content` | volumetric moisture θ | – |
| `groundwater/qx,qy,qz` | Darcy fluxes | m/s |
| `transport/concentration` | subsurface scalar | user units (e.g. psu) |
| `transport/concentration_surface` | surface scalar | user units |

Every dataset carries `units`, `long_name`, and `time` attributes.

## 2. Universal reading patterns

```python
import h5py, numpy as np

f = h5py.File("out/output.h5", "r")
nx = f["/grid/x_center"].shape[0]
ny = f["/grid/y_center"].shape[0]
nz = f["/grid/dz"].shape[0]

# snapshot times: integer-second group keys, sort numerically
times = sorted(int(t) for t in f["/surface/depth"])          # or any written variable

# 2-D surface field -> (ny, nx) array [j, i]
depth = f[f"/surface/depth/{times[-1]}"][:].reshape(ny, nx)

# 3-D subsurface field -> (ny, nx, nz) array [j, i, k]; k = 0 at the surface, downward
theta = f[f"/groundwater/water_content/{times[-1]}"][:].reshape(ny, nx, nz)

# monitor / audit tables: the 'columns' attribute is the authoritative header
def table(dset):
    raw = dset.attrs["columns"]
    cols = (raw.decode() if isinstance(raw, bytes) else raw).split(",")
    return cols, dset[:]
```

Rules that prevent wrong plots:

- **Inactive subsurface cells are NaN** (above `ktop` in
  terrain-following runs) — mask, never zero-fill:
  `np.ma.masked_invalid(theta)`.
- Vertical coordinates for subsurface plots come from
  `/groundwater/zcell/0` (reshape like any 3-D field), **not** from
  `/grid/z_center`, whenever the mesh follows terrain or has partial
  cells.
- `depth`, not `eta`, tells you where water is (a dry domain has `eta`
  below the bed by construction).
- **Boundary discharge belongs to the mass-audit table**, not to
  `uu·depth` at the outlet: the stored velocities carry a legacy scaling
  that under-reads instantaneous flux. Differentiate the audit's
  cumulative `boundary_outflow` column for a conserved hydrograph.
- Monitor tables have one row **per time step** (high-resolution);
  field snapshots only exist at `output_interval`.

## 3. Canonical plot: hydrograph from a monitor

```python
with h5py.File("out/output.h5", "r") as f:
    cols, data = table(f["/monitor/outlet"])       # e.g. columns: time,depth,uu
    t = data[:, 0]
    depth = data[:, cols.index("depth")]
import matplotlib.pyplot as plt
plt.plot(t / 3600, depth)
plt.xlabel("time [h]"); plt.ylabel("depth [m]")
plt.savefig("hydrograph.png", dpi=150)
```

## 4. Canonical plot: 2-D map and vertical slice

```python
with h5py.File("out/output.h5", "r") as f:
    x = f["/grid/x_center"][:]; y = f["/grid/y_center"][:]
    depth = f["/surface/depth/3600"][:].reshape(ny, nx)
plt.pcolormesh(x, y, depth, shading="nearest")
plt.colorbar(label="depth [m]"); plt.gca().set_aspect("equal")
```

Vertical (y–z) slice at column i — e.g. a salinity wedge:

```python
with h5py.File("out/output.h5", "r") as f:
    salt = f["/transport/concentration/36000"][:].reshape(ny, nx, nz)
    zc   = f["/groundwater/zcell/0"][:].reshape(ny, nx, nz)
    y    = f["/grid/y_center"][:]
i = 0
S = np.ma.masked_invalid(salt[:, i, :])            # (ny, nz)
Z = zc[:, i, :]
Y = np.broadcast_to(y[:, None], Z.shape)
plt.pcolormesh(Y, Z, S, shading="nearest")          # z is elevation; no flip needed
plt.colorbar(label="salinity"); plt.xlabel("y [m]"); plt.ylabel("z [m]")
plt.contour(Y, Z, S, levels=[17.5], colors="k")     # e.g. the 50 % isohaline
```

## 5. The audit tables (budgets and hydrographs)

All three are `/monitor/*` tables (one row per step, `columns` attribute
authoritative). Volume/mass columns other than the instantaneous
`volume`/`*_mass` are **cumulative since t_start**.

**`mass_audit`** (surface; columns
`time,volume,rain,evaporation,boundary_outflow,bc_inflow[,seepage],clamped`
— `seepage` only in coupled runs): closure identity

```python
cols, a = table(f["/monitor/mass_audit"])
c = {name: a[:, k] for k, name in enumerate(cols)}
seep = c.get("seepage", 0.0)
residual = c["volume"] - c["volume"][0] - (c["rain"] - c["evaporation"]
           - c["boundary_outflow"] + c["bc_inflow"] + seep + c["clamped"])
print("max |residual| [m^3]:", np.abs(residual).max())
# healthy: |residual| small relative to the largest cumulative term.
# Surface-only runs close near rounding; hard-wetting coupled regimes
# close to ~1e-3 of the budget (the validation gates' allowance).
q_out = np.gradient(c["boundary_outflow"], c["time"])     # conserved outflow [m^3/s]
```

(`clamped` measures a documented legacy-scheme defect — volume created
by the below-bed clamp; report it, don't hide it.)

**`gw_mass_audit`** (subsurface; columns
`time,volume,boundary_in,ss_storage,realloc,realloc_dropped,vloss`):
`Δvolume = boundary_in − ss_storage + realloc − vloss` to rounding.

**`transport_audit`** (transport runs; columns `time,surf_mass,subs_mass,
exchange,surf_source,surf_boundary,surf_adjust,surf_anchor,subs_boundary,
subs_adjust,subs_anchor`): per-grid identities
`Δsurf_mass = exchange + surf_source + surf_boundary + surf_adjust + surf_anchor`
and `Δsubs_mass = −exchange + subs_boundary + subs_adjust + subs_anchor`.
The `*_anchor` terms measure the preserved legacy ledger lags (see
`docs/theory/transport.md`); they are data, not errors.

In coupled runs the exchange appears symmetrically: positive surface
`seepage` volume is water the subsurface gave up.

## 6. Provenance in the file

`/frehg2` root attributes carry the code version, git SHA, creation
time, and the **complete configuration text** of the run — scripts can
recover every parameter from the output file alone:

```python
cfg_text = f["/frehg2"].attrs["config"]
```

(Attribute names: `version`, `git_sha`, `config`, `created` — inspect
`dict(f['/frehg2'].attrs)` if in doubt.)
