# Worked example: build a new case from scratch

This builds a small overland-flow case — a tilted 40 m × 5 m plane, uniform
rain for 10 minutes, draining out the east edge — from nothing.

**1. Make a case directory:**

```bash
mkdir -p mycase/input mycase/out
cd mycase
```

**2. Write the bed as a flat raster** `input/dem.dat`. Grid is 40×5
(`nx=40, ny=5`), sloping 1 % from west (high) to east (low), so bed
`= 0.4 − 0.01·x`. Any script that emits 200 values ordered `j*nx + i` works;
for a plane the value depends only on `i`:

```bash
python3 - <<'PY'
nx, ny, dx = 40, 5, 1.0
with open("input/dem.dat", "w") as fh:
    fh.write("# 40 x 5 bed, j*nx + i, 1% slope west->east\n")
    for j in range(ny):
        fh.write(" ".join(f"{0.4 - 0.01*((i+0.5)*dx):.4f}" for i in range(nx)) + "\n")
PY
```

**3. Write the rainfall series** `input/rain.dat` — 20 mm/h
(= 5.56e-6 m/s) for 600 s, then dry:

```
0      5.56e-6
600    5.56e-6
601    0
1200   0
```

**4. Write the configuration** `mycase.yaml`:

```yaml
simulation:
  id: mycase
  title: Tilted-plane overland flow

domain:
  nx: 40
  ny: 5
  nz: 1
  dx: 1.0
  dy: 1.0
  dz: 1.0
  bottom_elevation: {file: input/dem.dat}
  follow_terrain: false

time:
  dt: 0.2
  t_end: 1200
  output_interval: 60

modules:
  surface_water: true
  groundwater: false
  transport: false

surface_water:
  gravity: 9.81
  friction:
    law: manning
    coefficient: {constant: 0.03}
    thin_layer_depth: 0.1
  viscosity: {x: 1.0e-6, y: 1.0e-6}
  min_depth: 1.0e-6
  wetting_face_depth: 1.0e-6
  rainfall:
    series: {file: input/rain.dat}
  evaporation: {constant: 0.0}

initial_conditions:
  surface:
    eta: {constant: -1.0}        # dry start: below the bed everywhere

boundary_conditions:
  - name: east_outlet
    region: {polygon: [[39.5, -1.0], [41.0, -1.0], [41.0, 6.0], [39.5, 6.0]]}
    target: surface
    kind: outflow

output:
  filename: out/output.h5
  variables:
    surface: [eta, depth, uu, vv]
  monitors:
    - {name: outlet, i: 39, j: 2, variables: [depth, uu]}
```

**5. Validate, then run:**

```bash
build/src/frehg --validate mycase.yaml     # fix any reported problem first
OMP_NUM_THREADS=1 build/src/frehg mycase.yaml
```

**6. Plot the outlet hydrograph** with the monitor snippet from the
[output reference](output.md#plotting-with-python) (`out/output.h5`,
monitor `outlet`). You should see depth rise while it rains and recede
after 600 s.

To scale up, increase `nx`/`ny`, and — once the grid is large — run under MPI
with more threads. The physics keys stay the same.

**7. Make it coupled (optional).** Let the rain infiltrate a 1 m soil
column under the same plane: set `nz: 10` in `domain`, flip
`modules.groundwater: true`, and add the three subsurface sections —

```yaml
groundwater:
  timestep: {dt_init: 0.2, dt_min: 0.01, dt_max: 5.0}
  specific_storage: 1.0e-5

soil:
  types:
    - {name: loam, ksx: 2.89e-6, ksy: 2.89e-6, ksz: 2.89e-6,
       theta_s: 0.33, theta_r: 0.0, vg_alpha: 1.43, vg_n: 1.56}
  map: {constant: loam}

initial_conditions:
  surface:
    eta: {constant: -1.0}
  groundwater:
    moisture: {constant: 0.1}       # fairly dry soil
```

— then re-validate and re-run. The hydrograph now shows the
infiltration-delayed rise, and adding `seepage` to
`output.variables.surface` (plus `hydraulic_head` under
`variables.groundwater`) lets you watch the exchange itself. The complete
worked coupled benchmark is `benchmarks/b5-vcatchment/` — its README
documents every configuration decision against the published
intercomparison.
