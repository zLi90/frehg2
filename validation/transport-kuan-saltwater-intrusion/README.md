# kuan-saltwater-intrusion — tidal saltwater intrusion in a beach aquifer

A fully **coupled surface + subsurface, density-driven** solute-transport case
(SWE + groundwater + transport + baroclinic coupling). It reproduces the Kuan et
al. laboratory sandbox experiment of **tidally forced seawater intrusion in a
sloping beach aquifer**: tidal flooding over the beach face builds an *upper
saline plume* that sits above the fresh submarine groundwater discharge, while a
classic salt wedge intrudes from the seaward boundary. This is the most complete
transport case in the suite — it exercises every path at once (surface tracer,
subsurface dispersion tensor, seepage exchange, and the baroclinic feedback into
the Darcy solve).

It is the validation-suite copy of the project's own **`b6-kuan`** gate
benchmark (`../../benchmarks/b6-kuan/`); the configs read that case's
`input/` directory directly (the data is identical — it is not duplicated
here). Two configurations share it:

| file | forcing | seaward stage |
|---|---|---|
| `kuan-ss.yaml` | no tide (steady) | held at mean sea level 0.15 m |
| `kuan-td.yaml` | tidal | stage series from the shared `tide.dat` (129 600 samples @ 2 s) |

## The problem

A 2-D vertical slice of a laboratory beach: `nx=1, ny=68, nz=18`
(`dy = 0.05 m`, `dz = 0.04 m`, `follow_terrain: true`, bed from
the shared `bathymetry.dat`). Homogeneous sand (`Ks = 5.26e-3 m/s`, `θs = 0.46`, van
Genuchten `α = 5.9`, `n = 2.68`). Seawater (35 psu) enters from the seaward (y+)
boundary as a hydrostatic head; freshwater (0 psu) enters inland (y-) as a small
constant flux (`2.77e-5 m/s`). Density coupling is **on**: salinity feeds the
`r_ρ = 1 + β_ρ·s` buoyancy ratio (β_ρ hardwired at 7.44e-4, giving Δρ/ρ ≈ 2.6 %
at 35 psu — real seawater) in the groundwater Darcy solve.

## Reference solution

- **Primary:** the tidal-beach sandbox salt-interface data of Kuan, W.K.,
  Jin, G., Xin, P., Robinson, C., Gibbes, B., Li, L. (2012), *Tidal
  influence on seawater intrusion in unconfined coastal aquifers*, Water
  Resources Research 48, W02502, doi:10.1029/2011WR010678 — the
  laboratory experiment that established upper-saline-plume formation
  under tidal forcing (its no-tide and tide cases are this case's `ss`
  and `td` variants). The legacy code's wet-cell salinity override is
  labeled "Kuan 2019" after the follow-on modeling work, but the
  experimental data are the 2012 paper's.
- **Secondary:** the frehg2 P4 gate goldens (`out-ss-syncV4`, `out-td-syncV4`).

Recorded P4 gate result (both variants PASS): golden head max-abs-diff
0.012 m (ss) / 0.013 m (td); interface MAE vs the Kuan experiment 0.033 m (ss) /
0.039 m (td); salinity stays in [0, 35] psu. See the gate README for the full
adjudication table.

## The tidal salinity condition (why `kuan-td` has a whole-tank Dirichlet)

Every wet surface cell of the `td` golden is exactly 35 psu — the legacy run held
the flooding film at ocean salinity (the commented "Kuan 2019" override,
`scalar.c:258-262`). Without it the fresh beach discharge dilutes the advancing
flood film and the upper saline plume — the experiment's defining feature — never
forms. `kuan-td.yaml` reproduces this as an explicit `sea-surface-salinity`
condition (a whole-tank surface `scalar_value`, i.e. a wet-cell Dirichlet).
`kuan-ss.yaml` keeps the tide-cell-only value matching its own steady golden.

## Notes

- `coupling.mode: sync` — the committed legacy input flags async, but the goldens
  are the sync runs; the config pins sync to match them.
- Salt enters from the **seaward y+** side, the only lateral boundary the
  transport module wires a scalar ghost to (`sideCodeYp`) — the same constraint
  documented in `../transport-ogata-banks-column/`.

## Cost & recommendation

`t_end = 36000 s`, `dt = 0.05 s` outer; the sync common step rides the 0.5 s
groundwater cap. Measured: a 600 s smoke run took ~5 s on the reference Apple M3, and the P4
gate recorded **~5–6 min (ss)** and **~13 min (td)** serial.

- **`kuan-ss`** — ~5–6 min: borderline, **runnable locally**.
- **`kuan-td`** — ~13 min: over the 10-min guideline, **prefer HPC** (or run
  locally overnight). The tidal forcing is what makes it the costlier variant.

Run with `FI_PROVIDER=tcp` if your MPICH/libfabric install emits fabric-probe warnings or hangs at exit.

## Run

```bash
cd validation/transport-kuan-saltwater-intrusion
OMP_NUM_THREADS=1 FI_PROVIDER=tcp ../../build/src/frehg kuan-ss.yaml
OMP_NUM_THREADS=1 FI_PROVIDER=tcp ../../build/src/frehg kuan-td.yaml
```

Outputs land in `out-ss/` and `out-td/`
(`transport/concentration`, `concentration_surface`).
