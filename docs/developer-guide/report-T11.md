# Report — Stage-3 phase T11 (predictor–corrector consistency)

Date: 2026-11-04. Companion: `dod-T11.md`. Plan §12C (added by S3-A16,
owner decision O17); amendment S3-A17. Owner commits (O10).

## 1. Why T11 existed

The groundwater predictor (the Richards matrix) and corrector (the
water-content update) disagreed on face conductivities, and the
saturated-contact restore θ := θ(h) turned the disagreement into water.
Found outside the gates by the 2026Frehg2 v30 wind × heterogeneity study
(`BUG_realloc_ks_inconsistency.md`): on the H2 beach (24-class lnK,
Ks 5.6e-6 … 1.1e-3 m/s), `reallocation: legacy` created 9.45 and dropped
7.44 m³ per metre of shore per tidal cycle — net +2.0 m³/m, 72 % of the
net exfiltration, steady from cycle 2, 99.9 % in cells with a Ks jump to
a neighbour. Verification found FOUR inherited predictor/corrector
disagreements and one derived force error:

- **M1** — the corrector evaluated x/y faces with the minus cell's Ks
  and z faces with the lower cell's Ks; the predictor used each flank's
  own Ks (ported faithfully from legacy `darcy_flux` vs
  `compute_K_face`, documented as a harmless quirk).
- **M2** — the ponded top face: predictor 0.5·(Ksz + K(h)), corrector
  Ksz (O17c: Ksz in both halves — "otherwise infiltration will be
  underestimated").
- **M3** — a pressurized dry column (h_top > rZp·dz/2): the predictor
  sealed the face the corrector drained — the artesian/seepage mode.
- **M4→P5** — the terrain-slope gravity term differed between halves
  AND carried the wrong density form: a salty layer at rest on a slope
  was pushed by K·(ρ_f/ρ₀−1)·sin α (the derived S1). One change — the
  z-gravity form r²v¹ in both halves — closes both.
- The original density suspicion (O17d) was corrected at W11.0:
  `updateWaterContent` already ρ-weights every term; the corrector IS
  the mass balance. What remained of (d) is P5.

O17: one face rule for both halves, behind
`groundwater.predictor_corrector: consistent` (default `legacy`,
bitwise), with `face_conductivity: arithmetic|harmonic` (harmonic =
thickness-weighted harmonic Ks × arithmetic kr) as the opt-in interior
mean.

## 2. What landed

| Piece | Where |
|---|---|
| Shared face rule — `faceK<HARM>` over both flanks' own Ks; harmonic guards a ≤ 0 flank to an exact 0 (the Ks = 0 wall seals by construction) | `src/gw/FaceConductivity.hpp` (W11.2) |
| Corrector `computeFluxesImpl<PC,HARM>`: all three axes take each flank's own Ks under PC (M1); `classifyCoupledTopImpl<PC>` adds mode 3; the mode-3 top branch books its clamp to the cumulative `pc_clamp_m3` (run record, MPI-reduced) | `Corrector.cpp`, `RichardsSolver.*`, `Simulation.cpp` |
| Predictor: ponded top kzF = Ksz under PC (M2); mode 3 joins the Dirichlet seepage leg (M3); terrain RHS in the z-gravity form rXp²·vXp (P5) | `Predictor.cpp` (W11.3/11.4) |
| Compile-time templates per the fma rule (S3-A8.1): PC=false trees are the literal legacy expressions; `(void)x;` first-capture guards before every `if constexpr` (the nvcc extended-lambda rule, T10 correction 2) | both files |
| Tokens: parse/resolve/describe; schema enums + the harmonic-requires-consistent rule; docs: groundwater theory (T11 section replaces the "harmless quirk" note), exchange-flux (mode 3), parameters (240 tokens), case-setup §9b, coverage (140 rows) | W11.5 |

The consistent corrector reads nothing from the predictor's loop; both
halves re-evaluate the SAME shared rule.

## 3. Gates (first green 2026-11-04; tolerances measured, `stage3-pc.yaml`)

| Gate | Checks | Red (frozen frehg-T10) | Green |
|---|---|---|---|
| pc0 (4 unit tests) | hetero saturated column: restore books nothing (both face rules); legacy churns (red control in-tree); homogeneous consistent ≡ legacy bitwise (K3, `EXPECT_EQ`); ponded top face | — | 267/267 |
| pc1 (per-PR) | rl2 lane under `consistent`: churn ≤ 0.01·T; gw ledger; ≥ 20 saturated Ks-contrast faces | churn/T **2.91** | **5.86e-3**, ledger 1.1e-13, 320 faces |
| pc2 (nightly, v30 H2 at T/4) | creation ≤ 0.01 m³/m/cycle; ledger ≤ 1e-2 (the coupled-subcycled class) | **6.97** (Ks-only patch: 0.045) | **0.0003** |
| pc2-h0 (nightly) | H0 ledger; creation recorded | 0.021 created | 0.0005; ledger 4.0e-3 |
| pc3 (per-PR) | every interior face at its exact series resistance (z, stretched z, x; harmonic); steady | arithmetic interface carries 3.4× the series flux | exact (ratio 1.0000, \|dh\| = 0.0); arithmetic lane records 25.5× (anti-vacuity) |
| pc4 (per-PR) | Ks = 0 wall: max \|q\| = 0.0 exactly (harmonic) | arithmetic leaks 5.36e-6 m³/s | **0.0** |
| pc5 (per-PR) | wetting front crosses a 100× interface (anti-stall for the harmonic kr) | — | θ below the interface 0.469 |
| pc6 (per-PR) | (a) M2 ponded: churn ratio consistent/legacy ≤ 0.05, K1 at saturation = 0; (b) M3 artesian: seepage = analytic K·ΔH/L, restore ≤ 1e-8 of exchanged | restore = 36 % of infiltration; face sealed, 0.33 m³/h minted | ratio 0.0433, sat-tail **exactly 0.0**; 5.0000e-6 vs 5e-6 (rel 4.1e-14), restore 4.3e-14 |
| pc7 (per-PR, O17d) | lock + heat: gw ledger, scalar identity, K1 front churn | flat θV drift exactly 0 (already consistent) | ledger ≤ 2.7e-15 m³, identity ≤ 4.0e-14·M, churn ≤ 1.9e-12·V |
| pc8 (per-PR, S1) | salty slab on the 10 % slope at rest ≤ 1e-6·Ks·A; fresh control | salty **1.01e-2**, fresh 1.1e-13 | salty **1.42e-13**, fresh 1.15e-13 |
| pc-restart (per-PR) | checkpoint bitwise incl. `realloc_dropped` | — | PASS |
| pc-rank (nightly) | n = 1/2/4 strict | — | head ≤ 5.3e-15 |

The negative battery (`unit.pc_gates.negative`) feeds every checker a
corruption modelled on its defect (the 2.91·T churn, the arithmetic
series flux, the leaking wall, the stalled front, a drifting ledger, the
S1 creep) — all reject.

## 4. Adjudications on the way (S3-A17)

1. **pc6a is ratio-gated, not zero-gated** — the unsaturated front's
   restore books the single-Picard chord error of C(h), a PCA property
   in both modes, not a face mismatch; it vanishes at saturation
   (measured exactly 0.0).
2. **The "M5" edge-leg suspicion is NOT a defect** — the probe's
   full-width side polygons triggered the documented v2 corner rule (six
   edge legs whose network solution reproduces the computed heads
   exactly); the predictor edge leg is exactly 2c·dtg (matrix dump).
   rl2's head-row churn is the chord class, inside pc1's cap.
3. **pc2's ledger class** is the coupled-subcycled residual (~5e-3 m³)
   every v30 run shows, not pc1's strict 1e-6·T.
4. pc0's stack keeps every layer saturated (ψ ∈ [2.0, 3.3]) and measures
   throughput from the mid-stack face flux (`boundaryIn` is NET).

## 5. W11.8 — every Stage-3 groundwater lane under the `consistent` override

15 lanes; 9 PASS unchanged (rl3, pd2, pd5, dg1, et8, et9-veg, ws4, sx2,
b6-gw-smoke). 6 fail BY DESIGN — floors, pins or bounds authored on the
legacy trajectory (the T10 W10.4 class); each lane keeps `legacy` and
stays green un-overridden:

| Lane | Failing sub-check | Why by design |
|---|---|---|
| rl2 | legacy-churn floor ≥ 0.1·T (churn 3.06e-3 vs 5.23e-2) | the floor IS M1's encoding; pc1 replaces it |
| rl4 | \|Δ anchor\| disabled ≤ legacy/5 (1.297 vs 2.317) | the ratio presumes the M1-fed legacy anchor |
| et9 | noET 1271-dataset bitwise pin | scheme-changing override vs a pinned trajectory (W10.4 precedent) |
| ws1 | anti-vacuity ΣV_s 6.17e-5 < 8.8e-5; P4 vents 1.9e-6 < 1.8e-4; identities 1.07e-9 vs 1e-10·X | less water minted at the swash face ⇒ less stored/vented (physical); the identity slip is rounding reordering — `pc_clamp_m3` = 2.2e-7 m³ is 3000× the residual and booked consistently |
| ws1n | ΣV_s 7.69e-5 < 1.0e-4 | same floor class |
| sx1 | max \|c−35\| 2.12e-6 vs 1e-6 psu | bound authored on the T10 exchange trajectory (1.9e-7 then); still 10⁵× under the red class |

## 6. W11.9 — should `consistent` become the default?

Comparison build (`build-s3-df`, the one-line default flip, source
reverted after), same toolchain. Measured:

- **`default-bitwise` un-pinned: 11/13 lanes stay bitwise.** Only
  `gw-orient` (heterogeneous gw: `hydraulic_head`, `water_content`,
  `gw_mass_audit`) and `coupled-salt` (gw fields, then surface
  depth/eta and both concentrations through the coupling) change. b2 and
  g6c-coupled stay bitwise — homogeneous soil on flat terrain is the K3
  identity even where coupled.
- **With the new `BITWISE_PINS` entry
  (`groundwater.predictor_corrector: legacy`): 13/13 bitwise** — the
  gate keeps testing the v2.0 configuration whatever the default (the
  S3-A11 pin rule).
- **b3** (berino/glendale, 41× contrast, legacy golden): **PASS with
  tolerances unchanged — and the fit IMPROVES 2×.** Flipped: h=0 max
  dist 0.085 m / rms 0.041 m, h=−400 max 0.091 / rms 0.057. Stock
  (`legacy`): h=0 max 0.160 / rms 0.073, h=−400 max 0.170 / rms 0.088
  (allowed 0.2 / 0.1). The consistent faces move the heterogeneous
  water table TOWARD the golden, not away from it.
- **b5** (vcatchment, the Ks = 0 outlet wall, rain/sync): **PASS with
  tolerances unchanged** — the flipped run covers 0–40.0 h (stopped at
  the owner's direction at t = 144 001 s of the 432 000 s horizon; the
  intermodel reference curves end at 39.8 h, so the envelope comparison
  is complete): discharge 0/59 reference times outside, peak +2.7 % of
  the envelope mean, integral +1.5 %; ponding 0/68 outside, peak
  −1.5 %, integral −7.1 % (allowed 10/15/10 %). Against the stock
  binary over the common 24 h: peak discharge +0.024 %, max discharge
  difference 1.3 % of peak, cumulative outflow +0.40 %, ponding ≤
  0.84 % — the consistent faces barely move this homogeneous-per-zone
  case (the wall itself ran `arithmetic`, where the predictor face
  already carried ½K_loam; the harmonic opt-in is what seals it, pc4).
  Budgets: surface residual −10.2 m³ vs stock −10.0 (allowance 11.1),
  subsurface 6.2e-8, clamp creation 184 vs 180 m³.
- **b6 ss/td** (density, terrain-following beach, legacy golden): the
  head fit is unchanged (ss max \|dh\| 0.0138 m, td 0.0130 m, 0 % of
  cells outside), but the salinity field leaves the golden: ss 9.15 %
  of cells outside (cap 6 %), tidally averaged salt mass +27.5 % off
  (cap 10 %); td 8.25 % and +25.3 %. Attribution: the case is
  homogeneous sand, so M1 is inert (K3), and `pc_clamp_m3` = 0.0 — the
  whole divergence is **P5**, the terrain-gravity term: `legacy` pushed
  the salty wedge seaward with the spurious K·(ρ/ρ₀−1)·sinα force on
  the terrain-following beach, and the golden (a legacy Newton run)
  carries the same force. Against the PHYSICAL experiment the td
  interface MAE **improves** (0.0365 m vs the golden run's 0.0430 m);
  the ss MAE worsens (0.0613 vs 0.0245 m, both far inside the 0.1067 m
  experiment cap). A flip would need the b6 salinity tolerances
  re-pinned (or re-goldened) — an owner call.
- **et9** noET pin: breaks, as it must (the same lane as W11.8).

**Report; the owner decides** (owner item 8c). Until then `legacy` stays
the default and every golden lane is bitwise-unchanged.

## 7. Defects found and fixed on the way (S3-A17)

1. The plan's §12C premise ("the corrector is a volume balance") was
   wrong — corrected at W11.0 against the code before any gate was
   authored; the D-M/D-V variant choice died with it.
2. pc4's 2-cell wall was vacuous ((0+0)/2 = 0 under arithmetic too) —
   1-cell wall.
3. pc7's raw Σs·θV assert would have re-asserted booked limiter terms
   (drift 3.1e-6 while the booked identity closed at 1.2e-14) —
   identity-form gate.
4. pc3's drafted reference forgot gravity (ΔH = Δψ + Δz); the measured
   flux was exactly right. The x-case's apparent 12/11 anomaly was the
   edge-leg convention — the gate asserts per-face resistances in the
   interior, self-contained.
5. pc6a never saturated as drafted (wt 0.5 under a 0.3 m pond on a 2 m
   column) — wt 1.6; pc6b's gate first read the wrong ledger
   (`boundary_in` nets bottom vs top) — the mass-audit seepage column.
6. pc0's draft head BCs desaturated the low-K layers and measured the
   chord error instead of M1 (and harmonic kr → 0 sealed the column) —
   saturated BCs.

## 8. State at close

Per-PR: pc0 ×4, pc1, pc3, pc4, pc5, pc6, pc7, pc8, pc-restart, the
negative battery and 8 validates green; nightly pc2, pc2-h0, pc-rank
green; 267/267 unit tests; `default_bitwise` + et9 pins green on the
stock binary. Docs and checkers green by exit code (doxygen caught five
undocumented T11 members before CI could). `legacy` is bitwise-unchanged
everywhere; `consistent` waits on the owner's W11.9 decision, and
`face_conductivity: harmonic` on the pc3–pc5 data (owner item 8d).
