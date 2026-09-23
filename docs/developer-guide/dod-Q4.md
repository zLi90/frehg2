# Definition of Done — Q4 (Evaporation and salinity coupling)

Per v2 development plan §3 and §9, under the v1 plan's §11.4 DoD
conventions. Every item names its verification; commands were run green
on 2026-09-23 on macOS arm64 (**gcc-16.2** — Homebrew removed gcc-15
mid-phase, see "Environment note" below — MPICH 4.3, PETSc 3.25.1 +
hypre + Kokkos Kernels, Kokkos 5.1.1 OpenMP host backend;
`FI_PROVIDER=tcp`).

## Gates (blocking, §9 Q4 row)

- [x] **g4 — analytic drawdown + evaporative concentration** (per-PR,
      `regression.g4`): PASS, all four sub-criteria. (a) prescribed
      drawdown exact to 1.6e-14 m with closure 1.6e-14 m³; (b)
      concentration exact (2.6e-11 rel) with salt-mass drift 3.8e-13
      (was 6.6e-3 pre-fix — the velocity-wall flow-rate leak); (c) bulk
      mode end-to-end against the offline closed form E = 1.4691e-7 m/s
      (drawdown 3.8e-14 m, concentration 1.4e-12 rel); (d) dry-out
      positivity + audited-shortfall closure 9.5e-15 m³.
- [x] **g5 — Geng & Boufadel (2015) bare-soil salinization**
      (nightly-class, `regression.g5.geng2015`): PASS under the
      V2-A13/A14/A15 criteria. E(0) ratio 0.999 vs the Table-1 closed
      form; monotone decay 100 % of intervals; E(50h)/E(0) = 0.269
      (bound 0.5); extrapolated surface saturation 0.220 (band
      [0.0915, 0.5]); 50 h near-surface peak 116.8 g/L (> 60); salt-mass
      drift 6e-4 (bound 0.04); density contrast: 30 g/L edge at 1.800 m
      (β run, fingers) vs 1.988 m (β = 0), margin 0.188 m (≥ 0.05).
      Recorded, not gated (V2-A15): E(10h) factor 3.43 of digitized
      Fig. 3; moisture RMS 0.161/0.134 and salinity RMS 7.2/11.8 g/L at
      20/50 h vs the digitized (inconsistently-forced) Fig. 4/9b.
- [x] **Gate-first discipline (§6.1):** both gates authored, registered,
      and demonstrated failing in commit `4addaee` ("Q4 gates (red)")
      before the capability landed; §6.3 negative battery
      `unit.g45_gates.negative` (32 expectations) green throughout,
      including "the digitized reference passes its own criteria" and a
      corrupted-input failure per criterion.
- [x] **x-battery, transposed g5 slice** (`regression.g5.transposed`,
      nightly): the case rewritten as the y-z slice (nx/ny, dx/dy, and
      every polygon transposed) must pass the same criteria; finger
      positions are instability-set and deliberately not compared
      (§3.4).
- [x] **b1–b6 unchanged** under `prescribed`/`legacy_evap_allowance`
      defaults: full per-PR tier green post-capability (b1–b4 goldens,
      restarts, rank invariance, p1/p2/p3/p5, g1) plus the nightly b6
      ss/td pair — run under gcc-16 (doubling as the compiler-migration
      re-verification).
- [x] **atm/ module unit-covered:** `tests/unit/test_atm.cpp` — Tetens,
      q_sat (reproduces Table 1's own q_a = 2.9e-3), R_air, the closed-
      form E(0) = 1.4696e-7 m/s, α₁ (including the 0.0375 → 0.2
      equilibrium inversion and the 0.375 cap), condensation sign
      pass-through, MetForcing constants/series/relative-humidity
      folding.
- [x] **Mass audits extended:** the surface `clamped` column now absorbs
      the rain/evap dry clamp (signed), closing the surface identity
      through dry-out (g4(d)); `gw_mass_audit` gains the cumulative
      `evaporation` column (actual, off the realized top-face flux — the
      g5(i) observable). Documented in the output reference.
- [x] **p5** (pointer/forbidden discipline): `check_forbidden.sh` clean
      over the new code. **p1** (backend invariance) rides the per-PR
      tier (`regression.p1.*`, green).

## Deliverables (§3.2) — state

- [x] `src/atm/` bulk-aerodynamic module (header-only:
      `BulkAerodynamic.hpp` device formulas + `MetForcing.hpp` host
      sampling), one module, two consumers (swe, gw).
- [x] `atmosphere:` YAML block (schema + docs + resolved-config
      round-trip in one PR — the lockstep rule): air/surface
      temperature, pressure, exactly-one-of specific/relative humidity,
      wind speed, all SeriesOrConstant.
- [x] Open-water evaporation modes: `prescribed` (bitwise-legacy;
      optional `exclude` region — the rain-mask symmetry), `bulk` (wet
      cells only, at most the available depth, actuals audited;
      condensation deposits per V2-A13). The plan's third mode
      ("series") is the prescribed mode's series form + mask — no
      separate mode needed.
- [x] Soil evaporation potential→actual: `groundwater.evaporation`
      (mode bulk, region polygon; uncoupled runs only, schema-enforced),
      α₁(w_g) limiting with w_g evaluated at the **surface face**
      (1.5θ₀ − 0.5θ₁ — V2-A15: a cell mean half a cell down overstates
      surface moisture in a steep front), per-substep refresh, rides the
      legacy flux-top path unchanged.
- [x] Salinity coupling: `scalar_cauchy` (landed earlier, commit
      `03827ff`); the coupled `hi += 0.01` limiter allowance replaced by
      the same exact in-step factor behind
      `transport.legacy_evap_allowance` (default off; no golden pins the
      branch — b6 has no evaporating dry top — so the toggle is a
      documented escape hatch, not a live gate dependency).
- [x] Correctness fixes shipped with their gates: velocity conditions
      now correct face flow rates with the velocities
      (`WetDry.cpp::enforceVeloBc` — the g4(b) evaporation-induced edge
      salt leak, zero at rest); the evap/rain dry clamp audited.

## Owner touchpoint (§6.4)

- [x] Digitization delivered for visual check on 2026-09-23: three
      overlay plots (`benchmarks/g5-geng2015/reference/overlay_*.png`),
      CSVs + `DIGITIZATION.md` (exact vector extraction, axis residuals
      < 0.06 pt, self-checks green). **Owner sign-off: pending** — record
      it here when given. The V2-A15 reference-inconsistency finding
      (Fig. 4/9 forced ~4.28× harder than Table 1/Fig. 3; the wind-speed
      ratio) does not change the extraction, only what may be gated on.

## Amendments this phase

V2-A13 (g5(i) rederived from Fig. 3's resolvable content), V2-A14
(g5(iv) density inequality was inverted; corrected from Fig. 7/§3.5),
V2-A15 (reference figures mutually inconsistent; criteria re-anchored to
the self-consistent subset; internodal-conductivity finding; case-IC and
reallocation-mode corrections).

## Known limitations carried forward

1. Coupled runs cannot use `groundwater.evaporation` or `scalar_cauchy`
   (schema-rejected; §8.3: no gate exercises the combination). Lifting
   this needs a coupled evaporation gate first.
2. The E(t) decay-rate mismatch vs MARUN (factor ~3.4 at 10 h) is the
   upstream-K/starving-node internodal-conductivity difference; pinned
   by b2/b3 goldens, recorded not gated. Mesh-insensitive (2 mm top
   cells: −4 %).
3. The §8.1 general symmetry/orientation harness still does not exist;
   Q4 ships the transposed-g5 row only. Backfill remains x-gate work
   (with the §8.2 v1 matrix backfill and the V2-A11 outflow cells).

## Environment note

Homebrew auto-upgraded gcc 15.2 → 16.2 mid-phase (2026-09-23),
deleting `g++-15`. The build cache was migrated
(`CMAKE_CXX_COMPILER=g++-16`, stale `OpenMP_gomp_LIBRARY` cache entry
scrubbed); full rebuild zero-warning under `-Werror`, 176/176 unit
tests, and the full per-PR tier re-run green under gcc-16 — the
b1–b6 tolerances absorbed the codegen change without recalibration.
