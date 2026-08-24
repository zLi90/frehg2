# Removed features

Every capability of legacy Frehg 1.0 that Frehg2 deliberately does not
carry, with its rationale (upgrade plan §3.2). This registry is maintained
forever: anything found dead during implementation may be dropped **only**
by adding a row here in the same PR (plan §3.2, §11.1 rule 6). The
dropped-feature keyword tripwire in `scripts/check_forbidden.sh` fails CI if
any of these leak into `src/`.

| Feature | Legacy location | Rationale |
|---|---|---|
| Diffusive-wave approximation (`difuwave`) | branches in `momentum_source`, `shallowwater_mat_coeff`, `update_velocity` | Owner-directed removal. Unused by all benchmarks; removing it simplifies three hot functions. |
| Subgrid topography model (`use_subgrid`, ~20 `Data` members, 3 functions, lookup-table loader) | `shallowwater.c:1037-1349`, `initialize.c:1255` | Owner-directed removal. Its initializer call is already commented out in legacy (`initialize.c:68-69`) — dead code today. |
| Modified-Picard iterative GW scheme | SERGHEI `GwFunction.h` (gw_scheme=2); not present as such in legacy Frehg | Owner-directed removal. |
| Newton/FD-Jacobian GW path (`iter_solve=1`) | `groundwater.c:418-501`, `subroutines.c:198-279,451-541` | Owner decision (PCA only). Experimental: finite-difference Jacobian with inconsistent perturbation sizes (1e-3 vs 1e-8), silently falls back to PCA and mutates `param->iter_solve` mid-run on non-convergence. b6 (its only user) revalidates under PCA with relaxed tolerance + experimental data as the primary gate. |
| Modified van Genuchten (`use_mvg` branches) | `subroutines.c` guarded branches in C(h), K(h) | Unused by every benchmark (`use_mvg=0` everywhere). The **unguarded** `aev` cutoff in θ(h) is *kept* (live in b1/b2/b6). |
| Exponential retention fallback (`use_vg=0`, hard-coded `exp(0.1634·h)`) | `subroutines.c:296,399` | Undocumented site-specific curve; unused by benchmarks. |
| Aerodynamic evaporation model (`evap_model=1`) | `solve.c:302-309` | Hard-codes air temperature 20 °C, pressure, humidity, and a resistance fit — not a general model. Constant/time-series evaporation is kept. |
| `waterfall_velocity` weir treatment | `shallowwater.c:838` | Dead: its call site is commented out (`shallowwater.c:128`). |
| `waterfall_location` flag computation (`wtfx`/`wtfy`) | `shallowwater.c:687-716` | Found dead at P1 (plan §3.2 procedure): the flags' only consumer is the dropped `waterfall_velocity`, so computing them has no observable effect. The plan §10 P1 deliverable line is amended accordingly. The live `wtfh` face-wetting threshold (velocity limiters, `shallowwater.c:799-808`) is preserved — it is unrelated to the waterfall flags. |
| Free-surface `Sct == 0` fallback branch | `shallowwater.c:362-366` | Found unreachable at P1: `dept > 0` forces `Asz = dx·dy > 0`, so a wet cell's diagonal is always positive; only the dropped subgrid model could zero it. (The `dept == 0` row treatment itself was replaced at P3 by the zero-depth continuity closure — amendment A13, `docs/theory/surface-water.md`.) |
| `rain_sum` accumulator | `shallowwater.c:582`, `solve.c:118` | Found dead at P1: accumulated and reset every step, but every read is inside commented-out code. |
| `pseudo_seepage` async machinery | `groundwater.c:873`; async block commented at `solve.c:77-93` | Dead. Subcycled coupling is preserved without it. |
| `bctype_SW` edge codes | `configuration.c` | Parsed and never used; SWE boundaries are actually driven by tide/inflow regions. Replaced by the polygon BC system. |
| OpenMP remnants (`nthreads`, commented `omp` calls) | `FREHG.c:11`, config | Superseded by Kokkos. |
| Commented-out case-specific hacks ("Maina", "Henry", "Geng2015", "Kuan2019", "Toy delta", Warrick tabular curves, old assembly loops) | scattered | Undocumented experiment residue. Their *configurations* survive as benchmark YAMLs where relevant. |
| Hard clamp `s ∈ [0, 200]` on scalars | `scalar.c` | Replaced by configurable `transport.bounds` (default `[0, ∞)` plus monotonicity clipping); the magic 200 was a salinity-specific hack. Documented behavior change; b6 salinity never approaches the clamp. |
| All-mandatory flat `key = value` input parser | `configuration.c`, `utility.c` | Replaced by the strictly validated YAML v2 schema (plan §5.4); the parser also had a pointer-into-stack-buffer defect (`utility.c:41-59`). |
| Legacy `NT` parameter | `configuration.c:41` | Parsed and never used by the solver; its v1-draft descendant `max_step` is dropped by `tools/migrate_yaml_v1_to_v2.py` (rename table). |
| Post-allocation receive transfers (`allocate_recv`) and the send retries | `groundwater.c:1013-1014`, `:1026-1027`, `:983` (calls commented out), `:1428-1631` | Found dead at P2 (plan §3.2 procedure): the legacy author disabled the withdrawals with "limit recv to 1 adjacent cell to avoid instability" notes and never re-enabled the retries; withdrawing from a saturated neighbor destabilizes the front. The deficit direction of the consistency restore stays exactly as legacy left it (created volume is audited). The *send* direction is preserved and selectable (`groundwater.reallocation_surplus`, amendment A7). |
| Post-allocation `repeat` flag and `check_room` staleness | `groundwater.c:980`, `:890-900` | Found dead at P2: `repeat[0]` is set and never read; `check_room`'s pre-corrector room is superseded by the fresh room the sweep now maintains (amendment A7). |
| Adaptive-step truncation-error estimator and zone scan | `adaptive_time_step`, `groundwater.c:1662-1665` (`sat`/`unsat`), `:1670-1671` (`err`/`err_max`), with `hnm`/`dtn` state | Found dead at P2: the values are computed and never read (the controller acts on `dq_max` and the Courant cap only); `hnm` otherwise feeds only the dropped Newton path. |
| Terrain top-angle tables (`sintx`/`costx`/`sinty`/`costy`) | `map.c:287-290`, `:547-575` | Found dead at P2: their only consumer is the commented-out lateral-exchange block of `darcy_flux` (`subroutines.c:111-142`). |
| Newton scheme helpers (`compute_residual`, `compute_jacobian_fd`, `compute_dKdh`, `compute_dwcdh`, `h_incr`/`hp`/`resi`/`G*` Jacobian state) | `subroutines.c:198-279`, `:485-539`; `groundwater.c:418-501` | Dropped with the owner-directed Newton removal (plan §3.2); listed here because P2 confirmed nothing else reaches them. |
| `rss` / `r_async` monitor quantity | `solve.c:56-65`, `:149` | Found dead at P3 (plan §3.2 procedure): computed per step in coupled runs and written only to the legacy `r_async` monitor file; nothing in the physics reads it. It scaled with the dead `pseudo_seepage` machinery it was meant to diagnose. |
| Dispersive-flux in-branch boundary doublings | `scalar.c:786`, `:831` | Found dead at P4 (plan §3.2 procedure): both sit behind `actv == 1` guards on ghost cells whose `actv` is always 0 (`map.c:509-522`), so no input can reach them. The live edge doublings at `scalar.c:348-350` are preserved. |
| Scalar-mass ghost writes and exchanges (`sm_surf`/`sm_subs` edges) | `scalar.c:266-272`, `:333`, `:860-914`, MPI exchanges `:296`, `:495` | Found dead at P4: no live code reads the scalar-mass ghosts (the fluxes read concentrations); the coupled top-ghost write at `scalar.c:333` even indexes the 3D array with a 2D index. The concentration ghosts they travel with are preserved. |
| Surface scalar mass re-initialization remnant | `scalar.c:286-287` (commented), `initialize.c:630` | The initial `sm_surf = s Vs` seed is recomputed from scratch every step (`scalar.c:42`); the seed itself has no observable effect and is not ported. |
| `qseepage_old` lagged-seepage state | `groundwater.c:840` | Found dead at P4: written every coupled substep, read only by the dead `pseudo_seepage` interpolation (already dropped at P3 with the async machinery). |
