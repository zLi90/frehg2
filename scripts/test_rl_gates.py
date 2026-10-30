#!/usr/bin/env python3
"""test_rl_gates.py — the §6.3 self-check battery for the Stage-3 T1
reallocation gate criteria (rl2/rl3/rl4; Stage-3 plan §5).

1. Reference sanity (§6.3 rule 1) — the ledger identity is checked two
   independent ways on a synthetic consistent audit: the closed-form
   delta identity and a step-by-step reconstruction of the volume column
   from the cumulative increments.
2. Sensitivity (§6.3 rule 2) — every check fails on a corrupted input
   modeled on the physical failure it exists to catch.

Also runs the shared Stage-3 harness battery (plan §3.1) so the helpers
every rl gate relies on are provably able to fail.

Exits nonzero on any violated expectation.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tests" / "regression"))
sys.path.insert(0, str(ROOT / "scripts"))

import gate_realloc as rl  # noqa: E402
import test_s3_shared  # noqa: E402

failures: list[str] = []


def expect(condition: bool, label: str) -> None:
    print(f"  {label}: {'ok' if condition else 'FAIL'}")
    if not condition:
        failures.append(label)


def synthetic_audit(steps: int = 40, churn: float = 2.0e-3,
                    vloss_step: float = 1.0e-13, seed: int = 7) -> np.ndarray:
    """A consistent gw_mass_audit table: per-step increments drawn, the
    volume column RECONSTRUCTED from the identity so the ledger closes by
    construction (reference-sanity path 2)."""
    rng = np.random.default_rng(seed)
    audit = np.zeros((steps + 1, 7))
    audit[:, rl.TIME] = np.arange(steps + 1) * 30.0
    volume = 100.0
    audit[0, rl.VOLUME] = volume
    for n in range(1, steps + 1):
        d_bin = rng.normal(0.0, 1.0e-2)
        d_ss = rng.normal(0.0, 1.0e-5)
        d_re = rng.normal(0.0, churn / steps)
        d_dr = abs(rng.normal(0.0, churn / steps))
        d_vl = vloss_step
        audit[n, rl.BOUNDARY_IN] = audit[n - 1, rl.BOUNDARY_IN] + d_bin
        audit[n, rl.SS_STORAGE] = audit[n - 1, rl.SS_STORAGE] + d_ss
        audit[n, rl.REALLOC] = audit[n - 1, rl.REALLOC] + d_re
        audit[n, rl.REALLOC_DROPPED] = audit[n - 1, rl.REALLOC_DROPPED] + d_dr
        audit[n, rl.VLOSS] = audit[n - 1, rl.VLOSS] + d_vl
        volume += d_bin - d_ss + d_re - d_vl
        audit[n, rl.VOLUME] = volume
    return audit


print("ledger identity: two independent evaluations on a consistent table")
audit = synthetic_audit()
ok, msg = rl.check_ledger(audit, 1.0e-12)
expect(ok, f"closed-form identity closes ({msg})")
# Independent reconstruction: volume from cumulative columns directly.
recon = (audit[0, rl.VOLUME] + (audit[:, rl.BOUNDARY_IN] - audit[0, rl.BOUNDARY_IN])
         - (audit[:, rl.SS_STORAGE] - audit[0, rl.SS_STORAGE])
         + (audit[:, rl.REALLOC] - audit[0, rl.REALLOC])
         - (audit[:, rl.VLOSS] - audit[0, rl.VLOSS]))
expect(float(np.max(np.abs(recon - audit[:, rl.VOLUME]))) < 1.0e-12,
       "row-by-row reconstruction agrees with the volume column")

print("sensitivity: each check fails on its target corruption")
broken = audit.copy()
broken[-1, rl.VOLUME] += 1.0e-6
expect(not rl.check_ledger(broken, 1.0e-12)[0], "ledger: a 1e-6 m^3 hole fails")
tiny_t = synthetic_audit(churn=0.0, vloss_step=0.0)
tiny_t[:, rl.BOUNDARY_IN] *= 1.0e-3  # T ~ 2e-4 m^3; rebuild volume to match
tiny_t[:, rl.VOLUME] = (tiny_t[0, rl.VOLUME] +
                        (tiny_t[:, rl.BOUNDARY_IN] - tiny_t[0, rl.BOUNDARY_IN]) -
                        (tiny_t[:, rl.SS_STORAGE] - tiny_t[0, rl.SS_STORAGE]))
tiny_t[-1, rl.VOLUME] += 5.0e-13     # a rounding-class residual
expect(rl.check_ledger(tiny_t, 1.0e-12, 1.0e-12)[0],
       "ledger absolute floor covers rounding at low throughput (rl3 finding)")
expect(not rl.check_ledger(tiny_t, 1.0e-12, 1.0e-13)[0],
       "...and the same residual fails without the floor")

off_audit = synthetic_audit(churn=0.0, vloss_step=0.0)
expect(rl.check_realloc_exactly_zero(off_audit)[0], "exact-zero accepts true zeros")
tiny = off_audit.copy()
tiny[-1, rl.REALLOC] = 1.0e-300
expect(not rl.check_realloc_exactly_zero(tiny)[0],
       "exact-zero: even 1e-300 fails (exact, not approximate)")

expect(rl.check_vloss(off_audit, 1.0e-4)[0], "vloss bound accepts a quiet clamp")
leaky = off_audit.copy()
leaky[-1, rl.VLOSS] = 0.1 * rl.throughput(off_audit)
expect(not rl.check_vloss(leaky, 1.0e-4)[0], "vloss: a 10% clamp loss fails")

wc = np.full((8, 8, 12), 0.2)
wc[0, 0, 0] = np.nan  # a masked cell
expect(rl.check_theta_bounds(wc, 0.05, 0.45, 1.0e-12)[0],
       "theta bounds accept in-range moisture (NaN-masked)")
wet = wc.copy()
wet[3, 3, 3] = 0.45 + 1.0e-9
expect(not rl.check_theta_bounds(wet, 0.05, 0.45, 1.0e-12)[0],
       "theta bounds: 1e-9 above theta_s fails")

d = rl.deltas(audit)
created = d[rl.REALLOC] + d[rl.REALLOC_DROPPED]
dropped = d[rl.REALLOC_DROPPED]
expect(rl.check_percell_sums(created, dropped, d[rl.VLOSS], audit, 1.0e-12)[0],
       "per-cell sums accept the consistent decomposition")
expect(not rl.check_percell_sums(created + 1.0e-6, dropped, d[rl.VLOSS], audit,
                                 1.0e-12)[0],
       "per-cell sums: a 1e-6 m^3 booking mismatch fails")

expect(rl.check_substep_ratio(120, 100, 1.5)[0], "substep ratio accepts 1.2x")
expect(not rl.check_substep_ratio(200, 100, 1.5)[0], "substep ratio: 2.0x fails")

expect(rl.check_legacy_churn(audit, 1.0e-6)[0],
       "churn anti-vacuity accepts an active legacy lane")
expect(not rl.check_legacy_churn(off_audit, 1.0e-6)[0],
       "churn anti-vacuity: a lane with zero churn fails its own gate")

steady = np.full((1, 1, 40), 0.3)
expect(rl.check_steady(steady, steady, 1.0e-10)[0], "steadiness accepts a frozen window")
drift = steady + 1.0e-8
expect(not rl.check_steady(steady, drift, 1.0e-10)[0],
       "steadiness: a 1e-8 trailing drift fails")

qz = np.full(40, -2.0e-6)
expect(rl.check_uniform_flux(qz, -2.0e-6, 1.0e-8)[0], "uniformity accepts a uniform flux")
kink = qz.copy()
kink[20] = -2.0e-6 * (1.0 + 1.0e-6)
expect(not rl.check_uniform_flux(kink, -2.0e-6, 1.0e-8)[0],
       "uniformity: a 1e-6 relative kink fails")

gw_quiet = synthetic_audit(churn=0.0, vloss_step=1.0e-9)
tr_leg = np.zeros((5, 11))
tr_leg[:, rl.TR_SUBS_MASS] = 1.0e3
tr_leg[-1, rl.TR_SUBS_ANCHOR] = -4.0e-1   # an active legacy anchor (4e-4 of M)
tr_dis = np.zeros((5, 11))
tr_dis[:, rl.TR_SUBS_MASS] = 1.0e3
tr_dis[-1, rl.TR_SUBS_ANCHOR] = 2.0e-2    # 20x smaller, 2e-5 of M
expect(rl.check_anchor(tr_dis, tr_leg, 35.0, gw_quiet, 1.0e-4, 5.0, 6.0e-5)[0],
       "anchor bound accepts a 20x reduction within anchor_rel*M")
big = tr_dis.copy()
big[-1, rl.TR_SUBS_ANCHOR] = 2.0e-1       # only 2x smaller than legacy
expect(not rl.check_anchor(big, tr_leg, 35.0, gw_quiet, 1.0e-4, 5.0, 6.0e-5)[0],
       "anchor bound: an improvement below the 5x factor fails")
huge = tr_dis.copy()
huge[-1, rl.TR_SUBS_ANCHOR] = 2.0e-1
huge[:, rl.TR_SUBS_MASS] = 1.0e2          # 2e-3 of M: breaches anchor_rel too
expect(not rl.check_anchor(huge, tr_leg, 35.0, gw_quiet, 1.0e-4, 5.0, 6.0e-5)[0],
       "anchor bound: 2e-3 of the scalar mass fails anchor_rel")
quiet_leg = tr_leg.copy()
quiet_leg[-1, rl.TR_SUBS_ANCHOR] = 1.0e-9
expect(not rl.check_anchor(tr_dis, quiet_leg, 35.0, gw_quiet, 1.0e-4, 5.0, 6.0e-5)[0],
       "anchor anti-vacuity: a legacy lane with a negligible anchor fails")

print("shared Stage-3 harness battery (plan §3.1):")
failures.extend(test_s3_shared.run_shared_checks())

if failures:
    print(f"{len(failures)} rl gate-logic check(s) FAILED")
    sys.exit(1)
print("all rl gate-logic checks passed")
