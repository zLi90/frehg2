# Output reference

Everything lands in the single HDF5 file named by `output.filename`. The
layout (all field datasets are `float64`, snapshot times are keyed by integer
seconds):

```
/frehg2                     root group; attrs: version, git_sha (revision stamp), config text, created
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
/temperature/<var>/<t>       temperature snapshots (v2 Q5; either grid, degC)
/monitor/<name>    (rows, 1+k)  extendable table: time + the monitored variables
/monitor/mass_audit (rows, 7+)  surface volume budget (see below)
/monitor/gw_mass_audit (rows, 7) subsurface volume budget (see below)
/monitor/transport_audit     scalar-mass budget (transport runs; see below)
/monitor/temperature_audit   heat ledger (temperature runs; see below)
/checkpoint/<t>/...          full prognostic state (only if checkpointing)
```

Surface fields (`eta`, `depth`, `uu`, `vv`, `seepage`) are 2D arrays of
`NY*NX` values flattened as `j*NX + i`. Groundwater fields
(`hydraulic_head`, `water_content`, `qx`, `qy`, `qz`, and the Stage-3
cumulative per-cell diagnostics `realloc_created`, `realloc_dropped`,
`vloss` [m³], allocated only when requested) are 3D arrays
flattened as `(j*NX + i)*NZ + k`, `k = 0` at the top layer; cells outside
the active subsurface are NaN. Transport writes `concentration` (the
subsurface grid) and `concentration_surface` (the surface grid). Each
dataset carries `units`, `long_name`, and `time` attributes. `eta` is
absolute free-surface elevation; `depth` is `eta − bottom`; `seepage`
(coupled runs) is the surface-applied exchange rate [m/s], positive upward.

**Snapshot vs interval-mean rates (Stage-3 T4).** `seepage` is a PER-STEP
sample of the applied rate (the legacy qss observable): in a subcycled run
most surface steps apply no exchange, so the sample at an output instant is
usually 0 even while water is flowing — never integrate it. For analysis
use `seepage_rate`: the interval MEAN of the complete coupled surface gain
(applied seepage + bounce-back + vent deposits) in m/s, positive upward,
normalized over the true elapsed time between output writes, so
`sum over cells of rate · cell_area · Δt_out` equals the increment of the
audit's cumulative `seepage` column over the same interval exactly (the dg1
gate holds this to 1e-12). The frame at `t = t_start` is zeros (no interval
has elapsed). The accumulators restart bitwise (checkpoint fields
`acc_<name>` and the scalar `acc_t_prev`).

**The `/atmosphere` group (Stage-3 T5).** With the ET block,
`output.variables.atmosphere` selects `evaporation_rate`,
`soil_evaporation_rate` and `soil_evaporation_unrealised_rate`
(interval-MEAN rates in m/s, positive out — the same registry semantics
as `seepage_rate`: `sum(rate · area · Δt_out)` recovers the audited
volumes exactly) plus the cumulative complementarity clocks `wet_time`
and `soil_evaporation_time` [s] (per cell, `wet_time +
soil_evaporation_time ≤ t` always). Transport adds `salt_store`
[psu·m³] — the salt half of the per-cell paired store under
`transport.surface_dryout: store`; the transport audit's instantaneous
`surf_store` column is its global total, and
`Δ(surf_mass + surf_store + subs_mass)` closes on the unchanged ledger
RHS. Stage-3 T9 adds two surface variables: `water_store` [m³], the water
half of the paired store (requires `surface_dryout: store`), and
`film_removal_rate` [m/s], the interval-mean thin-layer removal per cell
(every mode: the cflLimiter rules P1 + P2; under the store also the clamp
film P3 and the end-of-step sweep). With `modules.temperature` under the
store, `temperature_audit` appends `surf_store`, the parked heat.

**Monitors.** `/monitor/<name>` is a growing table with one row per time
step. Its `columns` attribute is the CSV header, e.g. `time,depth,vv`, and
`i`/`j` attributes record the sampled cell.

**Row decimation (Stage-3 T4).** `output.monitor_interval: <seconds>`
decimates the ROW WRITES of every `/monitor` table (point monitors and the
audits alike). Audit columns still accumulate every step, so a decimated
run's row at time t is bitwise the row the every-step run writes at t;
rows are forced at checkpoint instants and at `t_end` (the dg2 gate).
Default 0 = one row per step, the v2.0 behavior.

**Mass audits.** `/monitor/mass_audit` (surface) is written whenever the
SWE module runs. Its columns (`columns` attribute:
`time,volume,rain,evaporation,boundary_outflow,bc_inflow[,seepage],clamped[,stage_deposit],film_removed,film_pending,stage_volume[,store_volume,store_release]`
— `seepage` and `stage_deposit` only in coupled runs; `stage_deposit`
(Stage-3 T2) is the cumulative coupled deposit on prescribed-stage cells,
removed through the reservoir under `inflow: {closure: external}` or
booked-in-place under `zero_gradient`, so the coupled surface identity
closes with stage conditions present) are the running water-volume budget,
reduced across all ranks. `clamped` measures the volume the legacy
below-bed clamp creates (a scheme defect reported as data; ~1.6 % of rain
in the b5 regime, effectively zero elsewhere) plus, since v2 Q4, the
signed volume of the rain/evaporation dry clamp (under the paired store
only the created part: the positive remainder parks instead).

Stage-3 T9 appends: `film_removed`, the cumulative counted thin-layer
removal (`theory/surface-water.md`, paths P1–P4); `film_pending`, the
instantaneous volume of sub-`min_depth` films still in η (reported dry by
the depth field, so outside `volume`; identically 0 under the store);
`stage_volume`, the instantaneous water on prescribed-stage cells (the
reservoir, outside the identity's domain); and under
`transport.surface_dryout: store`, `store_volume` (instantaneous ΣV_s) and
`store_release` (cumulative). The identity closes at **every** row:

    discard:  Δ(volume − stage_volume + film_pending) = rain − evaporation
              − boundary_outflow + bc_inflow [+ seepage] + clamped
              [− stage_deposit] − film_removed
    store:    Δ(volume − stage_volume + store_volume) = rain − evaporation
              − boundary_outflow + bc_inflow [+ seepage] + clamped
              [− stage_deposit]

(surface-only runs without stage conditions have `stage_volume` = 0, and
the T5 form `Δvolume = rain − evaporation − boundary_outflow + bc_inflow +
clamped` still closes when no film is ever removed — the g4(d) gate). The
run record's surface closure carries the per-path books
`film_removed_{isolated,small,clamp,sweep}_m3`, `film_pending_m3` and
under the store `store_volume_m3` / `store_release_m3`; the groundwater
closure adds `vent_thin_m3` (the P4 vents deposited under the store). The table is the
mass-consistent source of boundary discharge — prefer it over
`uu·depth` when you need conserved outflow, because stored velocities carry a
legacy scaling that under-reads instantaneous flux.

`/monitor/gw_mass_audit` (subsurface) is written whenever the groundwater
module runs, with columns
`time,volume,boundary_in,ss_storage,realloc,realloc_dropped,vloss[,evaporation]`
(`evaporation` only with `groundwater.evaporation` configured: the
cumulative actual evaporated volume through the bulk-soil zone, read off
the realized top-face flux — the g5(i) observable) — the
subsurface identity closes to rounding. Under
`groundwater.reallocation: disabled` (Stage-3 T1) the `realloc` and
`realloc_dropped` columns are exactly zero and the identity reduces to
`Δvolume = boundary_in − ss_storage − vloss`. It closes to rounding
(`Δvolume = boundary_in − ss_storage + realloc − vloss`; the evaporation
column is a sub-account of `boundary_in`, not an extra term). In coupled runs
the exchanged volume appears symmetrically: positive `seepage` in the
surface table is water the subsurface gave up.

`/monitor/transport_audit` (transport runs) is the scalar-mass budget:
every non-conservative piece of the legacy scheme is measured into its own
cumulative column (exchange, sources, boundary leaks, limiter/bounds clips,
and the ledger re-anchor terms); the closure identities hold to rounding in
every regime. The term-by-term definition is in the
[transport theory page](../theory/transport.md).

`/monitor/temperature_audit` (temperature runs, v2 Q5) is the same
identity on the temperature scalar — columns `surf_heat`, `subs_heat`,
`exchange`, `surf_source`, `surf_boundary`, `surf_adjust`, `surf_anchor`,
`surf_atmos`, `subs_boundary`, `subs_adjust`, `subs_anchor` — with the
atmospheric heat exchange in its own `surf_atmos` column and masses in
K·m³ (subsurface on the θ+κ retardation basis). The temperature module
also writes `temperature`/`temperature_surface` under `/temperature/`
(°C, selected by `output.variables.temperature`), and
`temperature_surface` is a valid point-monitor variable. See the
[temperature theory page](../theory/temperature.md).

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

## The run record (v2)

Every simulation writes `run-record.yaml` into the directory of the HDF5
output — a self-contained description of what was run and where the time
went (v2 plan §2A). The record is rewritten at every output flush and
finalized when the run completes, so a killed or crashed run still leaves a
truthful partial record (`provenance.finished: false` with the last
completed state).

Top-level sections:

| Key | Contents |
|---|---|
| `provenance` | frehg2 `version` + `git_sha` (the revision stamp, below), build type, hostname, MPI ranks and `[px, py]` decomposition, OpenMP threads, Kokkos backend, ISO-8601 start/end times, wall seconds, the input YAML path and its SHA-256, restart parentage (`restart_from`, `restart_time`) when restarted, `finished` |
| `configuration` | the fully **resolved** configuration — every default materialized; exactly what `frehg --resolve <input>` prints, byte-comparable against a re-resolve of the input |
| `modules` | module flags, coupling mode, density-coupling state |
| `boundary_conditions` | per BC: name, target, kind, value (series files carry their SHA-256), polygon vertex count and bounding box, global member-cell count |
| `timers` | the full hierarchical timer tree (count, min/mean/max seconds over ranks) — identical numbers to the end-of-run stdout table, plus the `init`, `io/output`, `io/checkpoint`, `io/run_record`, `halo` (message passing), and `monitors` sections |
| `closure` (additions) | Stage-3 T2 adds `stage_deposit_m3` (coupled runs) and, with any eta condition, `stage_ke_influx_J` + `stage_entrance_head_ratio` — the free-energy observables of the inflow closure (theory/surface-water.md); a zero_gradient run whose inflow imports kinetic energy without paying the entrance head warns at run end |
| `warnings` | non-fatal advisories, present only when there are any (Stage-3): configuration keys that are legal but inert under the current mode (also printed by `frehg --validate` as `WARNING:` lines) and runtime findings; a run with warnings is still a valid run |
| `solver` | per-system linear-solver telemetry: solves, iteration mean/max, rebuilds, retries, setup/solve seconds |
| `closure` | the final cumulative mass-audit budgets (surface and groundwater) |

`frehg --resolve config.yaml` prints the resolved configuration between
sentinel lines; `tools/check_run_record.py` validates a record (schema,
launch geometry, configuration round-trip) and is applied to every
regression gate run automatically (gate r1). `scripts/run_scaling.py` reads
its timings from the record rather than scraping stdout.

**Revision stamp.** `provenance.git_sha` names the source tree the binary
was built from; the HDF5 `/frehg2` `git_sha` attribute, the log header
(`frehg2 <version> (git <stamp>), ...`) and the usage text of `frehg` run
without arguments carry the same value. It is one of:

| Stamp | Meaning |
|---|---|
| `<sha12>` | the 12-digit abbreviated commit SHA (never a tag name); the tracked files matched that commit |
| `<sha12>-dirty` | tracked files differed from that commit, staged or not: the output cannot be reproduced from the commit alone. Untracked files do not count |
| `unknown` | the build found no git, or the source tree has no `.git` of its own (a release tarball, even one unpacked inside another checkout) |

The stamp is taken as each build starts, not when CMake configures, so an
incremental rebuild after a commit or an edit restamps the binary (v2 plan
V2-A20 finding 5); an unchanged stamp is not rewritten, so a no-op build
recompiles nothing. An edit made while a build runs shows up in the next
build's stamp. `check_run_record.py` (r1) checks that the field is present,
not that it is current; the unit test `unit.git_revision_stamp` gates its
currency.
