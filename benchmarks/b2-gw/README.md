# b2-gw — 1D vertical infiltration column vs Warrick (1971)

Groundwater-only gate case for phase P2 (plan §9): a 1 m column of Warrick's
Panoche clay loam discretized in 100 cells of 1 cm, initially at
theta = 0.033, ponded from above by a constant head of 0 (legacy
bctype_GW top = 1, htop = 0), no-flux everywhere else, marched on the
adaptive dtg from 1e-4 s to 2 s for 13 h.

## Gate

`tests/regression/run_regression.py b2` (ctest `regression.b2`):
element-wise against the legacy ASCII goldens (`legacy/benchmarks/b2-gw/out`,
converted on the fly) — head within max(1e-3 m, 1 %), theta within 0.005
absolute at every cell and output time — plus the wetting-front gate vs the
Warrick analytical profile: the depth of the theta = 0.165 crossing within
5 % of Warrick and within 1.2x the legacy golden's own error
(`tests/regression/tolerances/b2-gw.yaml` documents the measured margins).

`regression.b2_restart` runs the same case through a checkpoint at
t = 23400 s and requires the restarted run to match the uninterrupted one to
1e-12 relative. `regression.b2_rank_invariance.{strict,default}` replicate
the column to 4 x 4 (amendment A8: a 1 x 1 footprint cannot be decomposed;
lateral conductivity is zero under use_full3d = false, so per-column physics
is unchanged) and compare 1/2/4 ranks.

This case pins the legacy default `reallocation_surplus: drop` — the golden
embeds the legacy sweep's discard behavior (amendment A7).
