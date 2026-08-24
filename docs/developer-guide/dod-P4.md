# Definition of Done — P4 (Scalar transport and density coupling)

Per upgrade plan §10 P4 and §11.4. Every item names the command that
verifies it; all commands were run green on 2026-08-21 on macOS arm64
(gcc-15.2, Apple clang 15, MPICH 4.3, PETSc 3.25.1, Kokkos 5.1.1, HDF5
1.14.5 parallel, deps at `/Users/zhili/Codes/local`; regression runs pin
`FI_PROVIDER=tcp` — amendment A19). Plan amendments A18–A20 (the coupled
live hydrostatic side boundary, the transport-port resolutions of legacy
undefined/rank-dependent behavior, and the b6 gate adjudications — the td
sea-surface salinity provenance, the nearest-crossing interface
extraction, and the data-adjudicated golden tolerances) are logged at the
end of `FREHG2_UPGRADE_PLAN.md` in this commit.

## Deliverables (plan §10 P4) — all implemented and tested

- [x] `frehg::transport` STATIC library (`src/transport/`: ScalarSolver,
      SurfaceTransport, SubsurfaceTransport, Dispersion, Limiters.hpp)
      wired into all three driver loops after the flow step (legacy
      solve.c:112-116) — `cmake --build build` (zero warnings, `-Werror`).
- [x] Upwind + TVD-superbee advection on both grids with the exact legacy
      face-value stencils and the serial golden's guards (rank-invariant:
      staged far neighbors and the 3D corner exchange replace the legacy
      rank-local degradations — amendment A19) —
      `build/tests/frehg_unit_tests --gtest_filter='Limiters.*:TransportModule.*'`.
- [x] Surface diffusion, the full anisotropic dispersion tensor with cross
      terms (on the legacy volumetric face fluxes — the calibrated scale
      the goldens embed), and the local min/max monotonicity clipping —
      provenance table in `docs/theory/transport.md`;
      `--gtest_filter='TransportModule.DispersionTensorDiagonalOnVerticalFlux'`.
- [x] Surface–subsurface scalar exchange through the seepage (upwind
      concentration + the two-point diffusive term over the top cell's
      thickness; the legacy denominator read unwritten memory — A19),
      conservative between the grids and audited —
      `--gtest_filter='TransportModule.CoupledExchangeMovesScalarConservatively'`.
- [x] Baroclinic activation of r_rho/r_visc in the Darcy fluxes (legacy
      update_rhovisc + baroclinic_face as `gw/Baroclinic.cpp`; the P2
      hooks were already multiplied everywhere) plus the coupled live
      hydrostatic side boundary (amendment A18) —
      `--gtest_filter='Baroclinic.*'`.
- [x] Scalar boundary conditions: wet-cell Dirichlet regions (tide/stage
      salinity, generalized per A19/A20), discharge-paired inflow
      concentrations, and groundwater_side ghost values —
      `--gtest_filter='TransportModule.EtaPairedScalarValuePrescribesStageSalinity:TransportModule.SideScalarValueSetsSubsurfaceGhosts'`.
- [x] Transport HDF5 output (`concentration`, `concentration_surface`) and
      the `/monitor/transport_audit` scalar-budget table (every
      non-conservative legacy piece measured: exchange, sources, boundary
      leaks, limiter clips, ledger re-anchors — docs/theory/transport.md);
      checkpoint/restart extends the P3 state with {s_surf, s_subs,
      s_fu_old, s_fv_old, s_dzz_top} —
      `ctest --test-dir build -R regression.b6_restart` (bitwise).

## Exit criteria (plan §10 P4)

- [x] **Limiter/dispersion unit tests** —
      `--gtest_filter='Limiters.*'` (superbee φ(r) at the plan §8.1 sample
      ratios and the tvd face values) and the module-level dispersion and
      clipping tests (14 new tests, 146 total).
- [x] **1D square-wave advection: superbee produces no new extrema beyond
      1e-12** — `--gtest_filter='TransportModule.SquareWaveAdvectionAddsNoNewExtrema'`
      (40 steps at CFL 0.3; also the upwind variant).
- [x] **Closed-domain scalar mass conservation ≤ 1e-8/step** —
      `--gtest_filter='TransportModule.StillWaterDiffusionConservesSurfaceScalarMass:TransportModule.WaterTableColumnConservesScalarMass'`
      (strict Σ s·V conservation in the quasi-steady regimes), with the
      audited closure identity holding to rounding in *every* regime
      including accelerating flow and the coupled exchange
      (`TransportModule.*ClosesTheScalarBudget`,
      `TransportModule.CoupledExchangeMovesScalarConservatively`; the
      legacy scheme's measured non-conservative terms are documented in
      docs/theory/transport.md).
- [x] **b6 gate passes for ss and td variants** —
      `ctest --test-dir build -R 'regression.b6\.'` (regression_nightly
      label; ~6/~13 min serial): golden head max |Δ| 0.0121/0.0129 m,
      salinity outsiders 3.51 %/4.82 % (allowed 6 % — amendment A20);
      **primary criteria at full plan-§9 strictness: interface MAE vs the
      Kuan experiment 0.0334 m (ss) and 0.0389 m (td) against caps
      0.107 m and 1.5× golden (golden's own MAEs: 0.0245/0.0430 — the td
      model beats the golden), tidally averaged salt mass +5.5 %/+8.8 %
      (allowed 10 %)**. The td configuration carries the golden's
      sea-surface salinity as an explicit condition (A20; b6 README).
- [x] **Salinity bounded in [0, s_boundary] everywhere at all output
      times** — asserted inside the b6 gate on both grids and both
      variants (range exactly [0, 35]; emergent — no `transport.bounds.max`
      is configured).

## Previous phases' criteria (plan §11.4)

- [x] **P0–P3 criteria still green** — `scripts/ci_build_and_test.sh`
      (zero-warning gcc build, all unit/mpi labels, b1–b4 gates, b5
      restart + rank-invariance lanes, the new b6 restart gate, forbidden
      scan, parameter docs, Doxygen exit 0). b1's worst achieved/allowed
      ratio and the b2/b3/b4 records are unchanged; coupled restart
      remains bitwise with the transport state included.
- [x] **ASan/UBSan lane clean** — `scripts/run_sanitizers.sh` (clang lane
      per report-P0 §4): unit + mpi labels and the b1/b2 regressions under
      both sanitizers, no findings.

## Additional gates (plan §11)

- [x] Fidelity provenance table for the transport —
      `docs/theory/transport.md` (legacy `file:line` → Frehg2 function for
      every preserved piece, the preserved quirks, the A19 resolutions,
      and the scalar budget identity); `docs/theory/groundwater.md` gains
      the baroclinic-activation section and the A18 hydrostatic rule.
- [x] Newly-found dead code recorded — four P4 rows in
      `docs/theory/removed-features.md` (dispersive in-branch doublings,
      scalar-mass ghost writes, the initial sm seed, qseepage_old).
- [x] Parameter docs in lockstep — `python3 scripts/check_parameter_docs.py`
      (scalar_value semantics, concentration_surface, the coupled
      hydrostatic note).
- [x] Forbidden scan clean — `scripts/check_forbidden.sh`.

## Known facts recorded for later phases

- The b6 golden-comparison tolerances are data-adjudicated (A20): the
  golden itself was produced with the dropped Newton scheme, the td
  variant with the legacy wet-cell salinity override, and the ss near-shore
  surface ramp is the fixed point of an uninitialized-memory exchange
  coefficient (A19) — the experiment remains the primary criterion, at
  full strictness.
- The transport audit's anchor terms measure the legacy ledger lags (flux
  volume one step behind the velocities; Vgflux pre-reallocation). Any
  future scheme change that closes them must re-gate b6.
- The libfabric sockets provider wedges MPI_Finalize on this platform
  (post-completion only); every harness run pins `FI_PROVIDER=tcp` (A19).
- b6 runtimes collapsed under the A10 sync common step (~6/~13 min serial
  vs the plan's nightly-class estimate); the nightly label is kept per
  plan §11.2, but per-PR promotion is available to P5 if wanted.
