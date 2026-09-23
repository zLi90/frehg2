# Q4 Report — Evaporation (open water + soil) and salinity coupling

**Phase:** v2 Q4 (v2 plan §3)
**Completed:** 2026-09-23
**Machine:** macOS arm64 (M3), gcc-15 → **gcc-16.2** mid-phase (Homebrew
upgrade; see the DoD's environment note), MPICH 4.3, PETSc 3.25.1,
Kokkos 5.1.1 OpenMP host backend.

## What the phase delivered

Physically-forced evaporation replaces the v1 prescribed-only stub. One
bulk-aerodynamic module (`src/atm/`: Tetens → q_sat → Liu R_air →
Mahfouf–Noilhan flux, all inputs configuration through the new
`atmosphere:` block) feeds two consumers: per-cell open-water evaporation
(`surface_water.evaporation: {mode: bulk}` — wet cells only, clamped by
available depth, actuals audited) and soil evaporation with
potential→actual limiting through the surface-layer soil relative
humidity α₁ (Geng & Boufadel Eq. 6), applied as the groundwater top-face
flux over a configured region. The salt side is the `scalar_cauchy`
zero-total-flux condition (this phase's first commit): water leaves, salt
stays, and the top-cell limiter admits the exact in-step concentration
factor. The coupled path's legacy `hi += 0.01` allowance is replaced by
the same exact factor behind `transport.legacy_evap_allowance`.

## Gate results

**g4** (per-PR, four analytic sub-criteria): drawdown exact to 1.6e-14 m;
concentration to 2.6e-11 relative with salt mass to 3.8e-13; the bulk
formula end-to-end against its offline closed form to 3.8e-14 m; dry-out
closure to 9.5e-15 m³ with the shortfall audited.

**g5** (nightly, code-to-code vs MARUN, criteria per V2-A13/A14/A15):

| criterion | measured | bound |
|---|---|---|
| E(0) vs Table-1 closed form | ratio 0.999 | ± 5 % |
| monotone decay after 3 h | 100 % of intervals | ≥ 99 % |
| E(50h)/E(0) | 0.269 | ≤ 0.5 |
| surface saturation (extrapolated, 50 h) | 0.220 | [0.0915, 0.5] |
| near-surface salinity peak (50 h, z > 1.9 m) | 116.8 g/L | > 60 |
| salt-mass drift | 6e-4 | ≤ 0.04 |
| density contrast (30 g/L edge, β vs β=0) | 0.188 m deeper | ≥ 0.05 m |

Recorded, not gated (V2-A15): E(10 h) at factor 3.43 of digitized Fig. 3;
moisture profile RMS 0.161/0.134 and salinity RMS 7.2/11.8 g/L at
20/50 h against the digitized Fig. 4/9b.

## The three findings that shaped the phase

1. **The reference's figures are mutually inconsistent (V2-A15).** In a
   closed domain, profile moisture deficit must equal cumulative
   evaporation; digitized Fig. 4 holds 5–8× more deficit than Fig. 3
   supplies, with interval increments in the constant ratio 4.2–4.3 —
   exactly R_air(1 m/s)/R_air(5 m/s) = 4.28 under the paper's own
   resistance law. The profile figures evidently come from a ~4.3×
   harder-forced run than Table 1 states. Fig. 3's scale is anchored by
   the Table-1 closed form at t → 0 (0.1 % agreement here), so the gate
   anchors on the self-consistent subset and records the rest. Found by
   the digitization's own mass audit — the §6.4 protocol paying for
   itself.
2. **The E(t) decay rate is an internodal-conductivity artifact in both
   codes (V2-A15).** Frehg2's legacy-pinned upstream-K faces resupply
   the drying skin at near-saturated conductivity (stage-1 evaporation
   persists); MARUN's surface node starves. Mesh-insensitive (2 mm top
   cells move E(10 h) by 4 %), and the face rule is pinned by the b2/b3
   goldens. The one principled free choice — α₁'s evaluation point —
   was moved to the surface-face extrapolation.
3. **v1's prescribed evaporation was authored-unexercised** (the
   V2-A11 pattern; b1 sets it to 0): under it, a "closed" basin refills
   through the legacy-transmissive default edges without bound, and a
   velocity-wall workaround leaked scalar because `enforceVeloBc`
   corrected velocities but not the flow rates transport snapshots.
   Both pinned by g4 and fixed; two more case-configuration hazards
   (uniform-moisture saturated IC; `reallocation_surplus: drop`
   discarding 24 mm under strong drying) are recorded in V2-A15 and the
   agent case-setup guide.

## Cost

g4 ≈ 2 min; g5 ≈ 29 min (two 50 h serial runs) + the transposed-slice
row ≈ 29 min — comfortably nightly-class. The digitization is a one-time
cost, reproducible from `tools/digitize_geng2015.py`.

## Carried forward

Coupled-run evaporation (schema-rejected until a coupled gate exists);
the §8.1 symmetry-harness backfill; the V2-A11 outflow x-gate cells; the
owner's §6.4 overlay sign-off recorded in the DoD when given.
