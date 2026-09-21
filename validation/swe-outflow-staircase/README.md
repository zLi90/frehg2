# Transmissive outflow on a descending staircase (`outflow` × west edge)

A two-case minimal reproducer for the transmissive-outflow face-area defect
recorded as plan amendment **V2-A11**. It is the smallest configuration that
exhibits the defect, it has a closed-form answer, and it fails by ~500× against
the code as shipped — which is exactly what plan §6.1 asks of a gate authored
before the capability change that fixes it.

This is a *diagnostic* case, not a benchmark against published data. There is no
reference solution to digitize; the reference is Manning's equation.

## Setup

Both cases are the same domain:

| | |
|---|---|
| grid | 10 × 1 × 1, `dx = dy = 1 m` |
| bed | `bed_i = 0.1·i` — a 10 % staircase descending west, the same slope as `swere-superslab`'s DEM |
| friction | Manning, `n = 0.0036` |
| supply | steady `1e-4 m³/s` (`dy = 1 m`, so unit `q = 1e-4 m²/s`) |
| outlet | `kind: outflow` on the **west** face of cell `i = 0` |
| start | dry (`eta: {constant: -10.0}`) |
| run | `dt = 0.5 s`, `t_end = 4000 s` — long enough to reach steady state and drain |
| groundwater | off — the surface solver is the whole system under test |

They differ in one thing: **where the water enters.**

- **`upslope-fed.yaml`** — injected at the east-most cell (`i = 9`), so the
  water routes the length of the strip before reaching the outlet. This is the
  ordinary overland-flow arrangement.
- **`outlet-fed.yaml`** — injected into cell `i = 0` itself, the cell that
  carries the outflow BC. This stands in for the return flow that the
  superslab's groundwater delivers directly to its outlet cell.

The pair matters because the defect is **two-sided**: one wrong face area, with
the feeding direction setting the sign of the error.

## The analytic answer

At steady state the strip carries a uniform film at Manning normal depth and the
outlet passes exactly what is supplied:

```
h = (n·q / √S)^(3/5) = (0.0036 × 1e-4 / √0.1)^(3/5) = 2.715e-4 m
```

so ~0.27 mm of water everywhere, an outlet discharge of `1e-4 m³/s`, and — over
the 4000 s run — 0.400 m³ in and ~0.400 m³ out, with the below-bed clamp
contributing nothing.

## What actually happens

`src/swe/WetDry.cpp:144` and `:158` overwrite the west/south boundary face areas
with the *interior* face area, which the face kernel gauged over the higher of
two bed elevations. `applyOutflowCorrections`
([`src/swe/FreeSurface.cpp:311`](../../src/swe/FreeSurface.cpp)) squares it.
On a descending bed the outlet is therefore throttled by
`(deptx(j,1)/depth(j,0))²`. The east and north edges take their coefficient from
the cell's own `Sxp`/`Syp` and are correct.

| measured at t = 4000 s | stock | two lines deleted | analytic |
|---|---|---|---|
| **outlet-fed** — `depth[0]` | 0.1046 m (385× normal) | 2.039e-4 m (0.75×) | 2.715e-4 m |
| **outlet-fed** — outflow | identically 0 until t = 1001.5 s; 0.2907 m³ of 0.400 m³, leaving 0.1092 m³ trapped | 0.3997 m³ | 0.400 m³ |
| **upslope-fed** — outflow | 0.5530 m³ (**+38.3 %**) | 0.3953 m³ (−1.2 %, the film still on the strip) | 0.400 m³ |
| **upslope-fed** — clamped | 0.1573 m³ — **39.3 % of the input, minted** | 6.0e-16 m³ | 0 |
| **upslope-fed** — `depth[0]` | 6.78e-5 m (0.25×) | 5.56e-4 m (2.05×) | 2.715e-4 m |

Read the two cases together. Fed at the outlet, the borrowed area is *too
small*, the face is effectively shut, and the cell must fill to the 0.1 m sill
at `bed[1]` before it passes anything — the pool is one bed step deep, which is
the signature to look for in the field. Fed from upslope, the borrowed area is
*too large*, the outlet over-drains, η is driven below the bed, and the legacy
below-bed clamp manufactures the 39 % difference. Neither sign is visible
without the other case, so they are paired.

## Running

```bash
OMP_NUM_THREADS=1 ../../build/src/frehg upslope-fed.yaml
OMP_NUM_THREADS=1 ../../build/src/frehg outlet-fed.yaml
python3 check.py
```

Each case is a few seconds on one core. `check.py` asserts four things per case:
the outlet sits within a factor of 5 of normal depth rather than at the upslope
sill; outflow matches injection to 2 %; the below-bed clamp mints nothing; and a
directly-fed transmissive outlet passes water from the start rather than after a
long dead interval. It exits non-zero and names which assertions failed. The
depth band is deliberately an order-of-magnitude test — it has to separate
"normal depth scale" from "one bed step deep", which are 385× apart, and
cell-scale discretization on a 10-cell staircase moves the outlet depth by a
factor of ~2 by itself. The volume balance is the tight check.

**`check.py` fails today, by design** — 5 assertions, listed above. That is the
merged-failing gate §6.1 requires. Deleting the two assignments

```cpp
Asx(j, 0) = Asx(j, 1);   // src/swe/WetDry.cpp:144
Asy(0, i) = Asy(1, i);   // src/swe/WetDry.cpp:158
```

turns it green. Both directions were verified on 2026-09-21: 5 failures against
the stock build, a clean pass after deleting the two lines and rebuilding, and
the source restored afterwards. Note that both lines *overwrite a value that is
already correct* — the face kernel (`WetDry.cpp:49-57`) covers the halo column,
the bed ghost is a zero-gradient copy (`SurfaceSolver.cpp:214`), and the outflow
ghost sets `eta(j,0) = eta(j,1) − drop`, so `deptx(j,0)` is exactly the outlet
cell's own depth. The fix removes code; it adds none.

## Status

The fix is **diagnosed and verified but not landed.** Q0.3's deliverable was the
diagnosis; changing the west/south ghost rule is capability work that plan §6.1
puts behind a gate authored failing first, and §8.2 now names `Outflow × W` and
`Outflow × S` as the first deliverable of the x-gate work with this case as the
ready-made gate. The per-PR regression tier is 32/32 green with the two lines
removed — every gate that exercises `kind: outflow` (including b4-govindaraju,
the only one) uses an east-edge outlet, which is why nothing caught this.

The nightly tier has **not** been run against the fix. b5 uses `kind: eta`, so
its exposure is only the `Momentum.cpp:82` viscous term and transport diffusion
rather than the outflow coefficient — small, but to be measured, not assumed
(V2-A11 open item 2).

## Where the defect shows in real cases

Two validation cases sit on the defective path, both documented in their own
READMEs:

- [`swere-superslab`](../swere-superslab/README.md) — west-edge outlet;
  throttled ~1/812, the finding that came out of Q0.3.
- [`swe-vcatchment`](../swe-vcatchment/README.md) — south-edge outlet; a 0.192 m
  pool against a 0.2 m bed step, 1/252 throttle, 38.4 of 39.7 m³ passed.

Three other `kind: outflow` cases use east-edge outlets and are unaffected.
