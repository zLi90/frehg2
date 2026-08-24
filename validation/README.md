# Extended validation suites

The six `benchmarks/` cases are Frehg2's **blocking release gates** —
quantitative, CI-enforced, tied to the legacy goldens and the upgrade
plan's tolerances. The 19 cases in this directory are the **extended
validation** built on top of the released model: independently authored
cases probing the same physics against *additional* published
references, doubling as a library of ready-to-run example
configurations.

Two suites share this directory, distinguished by prefix:

| Prefix | Suite | Cases |
|---|---|---|
| `swe-*` | SERGHEI shallow-water ports ([Part I](#part-i--serghei-benchmark-ports-swe--re--swere-)) | 7 |
| `re-*` | SERGHEI Richards-equation ports (Part I) | 4 |
| `swere-*` | SERGHEI coupled SWE–RE ports (Part I) | 2 |
| `transport-*` | classic solute-transport benchmarks ([Part II](#part-ii--classic-solute-transport-benchmarks-transport-)) | 6 |

Every case directory is self-contained: a `README.md` (what the case is,
its reference solution, the measured result, and how to run it), one or
more frozen `*.yaml` configurations with a provenance header, and —
where needed — an `input/` directory of converted data files. Run any
case from inside its directory:

```bash
OMP_NUM_THREADS=1 ../../build/src/frehg <case>.yaml
```

`OMP_NUM_THREADS=1` is the fastest configuration for these small
benchmarks (a project finding from the validation gates); add
`FI_PROVIDER=tcp` if your MPICH/libfabric install emits fabric-probe
warnings or hangs at exit. `frehg --validate <case>.yaml` checks a
configuration against the schema without running it. Run outputs
(`out*/`) are not committed; each run regenerates them.

**Reference wall times** quoted below were measured serially
(`OMP_NUM_THREADS=1`) on an Apple M3 CPU; scale to your hardware. As a
rule of thumb, cases whose serial estimate exceeds ~10 minutes are
marked **HPC** — run them on a workstation batch queue or cluster rather
than interactively.

---

## Part I — SERGHEI benchmark ports (`swe-*`, `re-*`, `swere-*`)

These 13 cases are SERGHEI benchmark tests converted into Frehg2 v2
(YAML) validation setups. They answer one question: **how much of the
published SERGHEI verification suite can Frehg2 reproduce, and at what
cost?**

Source material is the published SERGHEI benchmark suites (their input
decks are **not** shipped in this repository; the per-case YAML headers
cite the exact source deck within each suite's layout):

| Suite | Paper | Source decks (SERGHEI supplements) |
|---|---|---|
| SERGHEI-SWE | Caviedes-Voullième et al., *GMD* **16**, 977–1008 (2023) | `serghei-swe-tests/` |
| SERGHEI-RE | Li et al., *GMD* **18**, 547–562 (2025) | `serghei-re-tests/` |
| SERGHEI-SWE-RE | Zheng et al., *GMD* (2026) | `serghei-swe-re-testcases/` |

### Feasibility analysis

Frehg2 is a structured-grid (Cartesian nx·ny·nz) finite-volume model
with a **semi-implicit** 2D shallow-water surface module, a **PCA**
variably-saturated Richards groundwater module (van Genuchten
hydraulics), and scalar transport. The gating question for each SERGHEI
case is whether its boundary conditions, forcing, and geometry fall
inside that envelope.

**Boundary-condition translation** (holds across all suites):

| SERGHEI | Frehg2 |
|---|---|
| SWE outer `REFLECTIVE` | closed edge (no BC entry) |
| SWE outer `TRANSMISSIVE`, `extbc` bctype 5 / 9 | `kind: outflow` |
| SWE `extbc` bctype 6 (imposed level) | `kind: eta` |
| SWE `extbc` bctype 7 (imposed discharge) | `kind: discharge` — **value is TOTAL volumetric [m³/s]**, so multiply SERGHEI's unit q [m²/s] by the inflow width (ny·dy) |
| SWE `extbc` bctype 10 (stage/wave hydrograph) | `kind: eta` + series |
| SWE outer `PERIODIC` | **unsupported** |
| GW `gwbc` bctype 2 (Dirichlet head) | `kind: head` (constant / hydrostatic / series) |
| GW `gwbc` bctype 3 (fixed flux) | `kind: flux` (constant) |
| GW `gwbc` bctype 6 (atmospheric flux series) | `kind: flux` (series) |
| GW `gwbc` bctype 9 (free drainage) | `kind: flux`, `value: {gravity: true}` |
| GW `gwbc` bctype 8 (surface coupling) | no explicit BC — Frehg2's coupler owns the exchange |

**What Frehg2 cannot represent:**

| Case(s) | Suite | Blocker |
|---|---|---|
| `column/small` | SWE | requires PERIODIC BC |
| `rousseau` | SWE | Darcy–Weisbach friction (Frehg2 has only Manning / Chézy) |
| `rainInfTest/hetInf`, Tatard/Thiès, any active-infiltration hydrology case | SWE | Frehg2's surface module has no infiltration source term |
| RE Case 6 (crop irrigation) | RE | root-water-uptake sink, soil-evaporation model, irrigation source term, partitioned atmospheric BC — none exist in Frehg2 |
| Lake Taihu, Lower-Triangle HighRes (~118 M cells), Malpasset | SWE / RE | HPC-scale and/or no shipped inputs |

**Already covered by the Frehg2 gate benchmarks** (not re-created here —
see `../benchmarks/`): SWE Govindaraju overland plane →
**b4-govindaraju**; RE Case 3 layered-soil infiltration (Kirkland 1992)
→ **b3-kirkland**; SWE-RE Case 3 tilted-V catchment →
**b5-vcatchment**.

**Feasible-with-caveats, deliberately excluded from this port:** all
SERGHEI **shock / transcritical / supercritical / tsunami-runup /
dam-break / Thacker** cases *run* in Frehg2, but its semi-implicit
solver will not reproduce shocks and supercritical transitions with the
fidelity of SERGHEI's explicit augmented-Roe scheme. Likewise the two
coupled inflow-bump cases (SWE-RE Case 1/2, 260²–520² cells with sync
coupling) are genuinely HPC-scale. These were left out to keep the
collection focused on cases that both map cleanly *and* validate
meaningfully at workstation cost; they can be added later if wanted.

### Cases included (Tier A — clean-feasible, non-duplicate)

| Directory | Suite | Modules | Physics / reference |
|---|---|---|---|
| `re-1d-infiltration` | RE 4.1 | GW | 1D infiltration into dry column — Warrick (1985) analytical |
| `re-1d-drainage` | RE 4.2 | GW | 1D gravity drainage of saturated column — Abeele (1984) |
| `re-2d-open-lateral` | RE 4.4 | GW | 2D infiltration, lateral outflow, fluctuating WT — Hydrus2D |
| `re-road-embankment` | RE 4.5 | GW | 2D infiltration into a road embankment — GFDM (Chávez-Negrete 2018) |
| `swe-bump-at-rest` | SWE 4.1.1 / A1 | SWE | lake-at-rest over a bump (C-property), emerged + immersed |
| `swe-lake-at-rest` | SWE | SWE | 2D still water over irregular bathymetry (well-balancing) |
| `swe-subcritical-bump` | SWE 4.1.2 | SWE | subcritical steady flow over a bump — SWASHES 3.1.3 |
| `swe-steady-plane` | SWE | SWE | steady normal-depth plane flow, slope2 + slope8 |
| `swe-macdonald-subcritical` | SWE 4.1.2 | SWE | MacDonald subcritical flume — SWASHES / MacDonald (1995) |
| `swe-vcatchment` | SWE | SWE | surface-only V-catchment overland flow |
| `swe-microtopo` | SWE | SWE | rainfall runoff over microtopography (3 variants) |
| `swere-lateral-hillslope` | SWE-RE 4 | SWE+GW | coupled lateral hillslope flow, 5-year run |
| `swere-superslab` | SWE-RE 5 | SWE+GW | coupled heterogeneous "superslab" infiltration/return flow |

### Compute cost

Calibrated per-cell-step cost (serial): **surface ≈ 1–4×10⁻⁷
s/cell-step** (the b3/b4 gates give ~1×10⁻⁷; small strips measure
2–4×10⁻⁷ because the per-step semi-implicit free-surface solve carries
fixed overhead that dominates at low cell counts), **groundwater ≈
2.0×10⁻⁶ s/cell-step**. Estimate = cell-count × steps × per-step cost.
Adaptive groundwater stepping is usually faster than `dt_max` implies,
so GW estimates are upper bounds.

Most Tier-A cases cost seconds to a few minutes. **Three are HPC-tier**
(included for completeness; use a batch queue):

| Case | Why HPC | Rough serial cost |
|---|---|---|
| `re-2d-open-lateral` | 2D Richards, `dt_max = 300 s` over a long horizon | ~42 min |
| `swere-superslab` | coupled sync run marches the small surface CFL step | ~15–30 min |
| `swere-lateral-hillslope` | coupled, 5-year horizon on the surface CFL step | ~2.3 h |

Measured reference wall times (serial spot-runs):

| Case | Wall time | Note |
|---|---|---|
| `re-1d-drainage` | 1.6 s | GW calibration anchor |
| `re-1d-infiltration` | 27 s | GW calibration anchor |
| `swe-lake-at-rest` | 2.2 s | surface anchor |
| `swe-steady-plane` (slope2/slope8) | ~1 s / <1 s | uniform normal depth reached |
| `swe-macdonald-subcritical` | ~7 s | outlet q within <1 % of target |
| `swe-bump-at-rest` (emerged/immersed) | 156 s / 29 s | C-property still-water tests |
| `swe-vcatchment` | 43 s | surface-only overland flow |
| `swe-microtopo` (s0.01/s0.02/s0.05) | 26 s / 6.8 min / 21 s | s0.02 is a 209×79 grid |
| `swe-subcritical-bump` | ~2.3 min | needs dt = 0.001 s (stability); most costly SWE case |

The one case without a measured runtime, `re-road-embankment` (2D
Richards), validates and is in the seconds-to-low-minutes range; its
README carries the estimate.

### Cross-cutting conversion findings (apply to any future SERGHEI port)

- **Discharge units.** Frehg2's `kind: discharge` value is a **total
  volumetric discharge [m³/s]**, whereas SERGHEI's `bctype 7` q is a
  **unit discharge [m²/s]**. Multiply by the inflow width (ny·dy).
  Getting this wrong over-feeds the domain by the width factor (here up
  to ~10×) and blows the solution up.
- **Steady-plane friction.** A "normal depth" test is meaningless
  without friction; SERGHEI's `steadyPlane/slope2` deck flags
  `friction: none` yet supplies the normal-flow q and depth for
  `roughness: 0.033`. We run it with n = 0.033 (see that case's README).
- **Frictionless started-from-rest transients** (e.g.
  `swe-subcritical-bump`) can require a much smaller `dt` than the
  SERGHEI deck for the semi-implicit free-surface solve to stay bounded.
- **Shallow (cm-scale) steady flows** land inside Frehg2's
  `thin_layer_depth` friction-regularization band (default 0.1 m),
  giving a ~20 % depth/discharge offset from the Manning analytic;
  deeper flows match to <1 %.

---

## Part II — Classic solute-transport benchmarks (`transport-*`)

These 6 cases are **classic literature benchmarks for scalar transport**
— the transport-module companion to the SERGHEI flow ports of Part I:
where those answer "can Frehg2 reproduce the SERGHEI flow benchmarks",
these answer **"can Frehg2 reproduce the canonical advection–dispersion
and variable-density transport benchmarks?"** They are
classic-literature cases, not SERGHEI ports — hence their own
`transport-` prefix. The sixth, `transport-kuan-saltwater-intrusion`, is
the validation-suite copy of the project's own `b6-kuan` gate benchmark
(it reads that case's `input/` data directly).

| # | Case | Modules exercised | Reference | Validates |
|---|---|---|---|---|
| 1 | [transport-surface-tracer-advection](transport-surface-tracer-advection/) | SWE + transport | pure advection `x_f = u·t` | surface advection schemes (upwind vs superbee) |
| 2 | [transport-ogata-banks-column](transport-ogata-banks-column/) | GW + transport | Ogata & Banks (1961) erfc | 1-D subsurface advection–dispersion |
| 3 | [transport-henry-saltwater-intrusion](transport-henry-saltwater-intrusion/) | GW + transport + density | Henry (1964) | variable-density salt wedge (semi-analytical) |
| 4 | [transport-goswami-clement-swi](transport-goswami-clement-swi/) | GW + transport + density | Goswami & Clement (2007) | variable-density salt wedge (laboratory) |
| 5 | [transport-elder-free-convection](transport-elder-free-convection/) | GW + transport + density | Elder (1967); Voss & Souza (1987) | density-driven free convection (no through-flow) |
| 6 | [transport-kuan-saltwater-intrusion](transport-kuan-saltwater-intrusion/) | SWE + GW + transport + density | Kuan et al. (2012) lab + P4 golden | fully coupled tidal beach transport |

Cases 1→6 form a deliberate complexity ladder: surface-only → subsurface
→ +density → +laboratory → +free convection → fully coupled.

### Measured results and reference cost

| Case | Key result | Wall time | Tier |
|---|---|---|---|
| tracer advection | front within 1 cell of `u·t`; superbee ~9× sharper than upwind | < 1 s ea. | workstation |
| ogata-banks | breakthrough `t50` −2 % vs `L/v`; front RMS 0.004 vs erfc | < 1 s | workstation |
| henry | steady Henry wedge; 0.5-isochlor toe 0.73 m (ref ~0.8–1.0) | ~6 s | workstation |
| goswami-clement | steady lab wedge; toe 0.285 m (ref ~0.24–0.26) | ~27 s | workstation |
| elder | symmetric two-lobe convection cell, 58 m penetration @ 10 yr | ~110 s | workstation\* |
| kuan `ss` | matches the P4 steady golden (head 0.012 m) | ~5–6 min | workstation |
| kuan `td` | matches the P4 tidal golden (head 0.013 m) | ~13 min | **HPC** |

\* Elder runs at workstation cost **only as a decaying-source
demonstration**; the canonical sustained-source Elder is HPC/research
(see its README and the limitation below). Approximate scaling from
these runs: surface transport ~1×10⁻⁷ s/cell-step, groundwater ~2×10⁻⁶
s/cell-step; density coupling adds ~2–4×.

### Three Frehg2 constraints that shape every transport case

Learned while building this suite — load-bearing for anyone extending it:

1. **Subsurface scalar injection is y+ only.** The transport module
   wires a prescribed side-ghost concentration solely on the **y+
   (north) lateral boundary** (`src/transport/SubsurfaceTransport.cpp`
   `sideCodeYp`, legacy `s_yp`) — there is *no* injection path on y−,
   x±, or the top/bottom faces. Any subsurface tracer/seawater inflow
   **must** enter from y+; flow runs north → south. This dictates the
   geometry of cases 2–5 (and is why Elder's top-face source cannot be
   reproduced).
2. **Density contrast is a compile-time constant.** The baroclinic law
   is `r_ρ = 1 + β_ρ·s` with `β_ρ` hardwired at **7.44e-4**. You cannot
   dial `Δρ/ρ` in the YAML — you set it through the seawater salinity
   surrogate `s = (Δρ/ρ)/β_ρ` (33.6 for Henry's 0.025, 35 for
   Goswami–Clement's 0.026, 269 for Elder's 0.2). The viscosity ratio is
   likewise fixed.
3. **Dispersion acts on volumetric face fluxes.** The dispersion term is
   applied to face volume-fluxes, not Darcy velocities (the face-area
   factor is not divided out), so a configured `dispersion.longitudinal`
   is **not** the physical dispersivity. Quantitative checks in this
   suite use quirk-independent metrics (advective `t50`, front celerity,
   wedge-toe *character*) rather than absolute isochlor spread.

### Coverage summary

- **Fully validated at workstation cost:** cases 1, 2, 3, 4, and
  `kuan-ss` — each matches its analytical / semi-analytical / laboratory
  / golden reference.
- **Qualitative demonstration only:** case 5 (Elder) — Frehg2 has no
  sustained top-face scalar Dirichlet, so only a decaying-source
  convection demonstration is possible; it reproduces the Elder
  two-plume *pattern* but not the steady field.
- **HPC recommended:** `kuan-td` (~13 min serial) and any quantitative
  Elder study.

---

## How these differ from `benchmarks/`

- `benchmarks/` cases are **gates**: the regression harness runs them
  and fails the build on any tolerance violation. They never change
  without a plan amendment.
- `validation/` cases are **evidence and examples**: each README records
  the measured result against its reference (and, where relevant, the
  model limitation it exposes). They are run manually, not by CI.

Per-case provenance (the exact source `.input` files for the SERGHEI
ports; the cited literature for the transport suite) is in each case's
YAML header. Both suites were authored against the v1.0.0 release
binaries.
