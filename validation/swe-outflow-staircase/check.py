"""
Check the transmissive-outflow staircase reproducers against the analytic
answer (plan amendment V2-A11).

Both cases are the same 10 x 1 strip with a 10 % descending bed and a steady
1e-4 m3/s supply; they differ only in where the water enters. The steady state
is Manning normal depth everywhere and an outlet discharge equal to the
supply:

    h = (n q / sqrt(S))^(3/5) = (0.0036 * 1e-4 / sqrt(0.1))^(3/5) = 2.715e-4 m

Run both cases, then this script:

    OMP_NUM_THREADS=1 ../../build/src/frehg upslope-fed.yaml
    OMP_NUM_THREADS=1 ../../build/src/frehg outlet-fed.yaml
    python3 check.py

This script is a *negative test in its present form*: against the code as
shipped it FAILS, by ~500x on the outlet-fed depth. That is the point — it is
the gate the plan §8.2 BC kind x side matrix cell needs, and §6.1 requires it
to be merged failing before the fix lands. Deleting the two face-area copies
at src/swe/WetDry.cpp:144 and :158 turns it green; nothing else is needed.
"""

from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent

MANNING_N = 0.0036
UNIT_Q = 1.0e-04          # m2/s (dy = 1 m, so this is also the total m3/s)
SLOPE = 0.1
SUPPLY = 0.4              # m3 injected over the 4000 s run

NORMAL_DEPTH = (MANNING_N * UNIT_Q / np.sqrt(SLOPE)) ** 0.6

# The depth check is an order-of-magnitude band, not a friction gate. What it
# has to separate is "normal depth scale" from "one bed step deep", and those
# are ~385x apart. Cell-scale discretization on a 10-cell staircase moves the
# outlet depth by a factor of ~2 on its own (the fixed upslope-fed lane settles
# at 2.05x normal depth, the fixed outlet-fed lane at 0.75x), so the band is
# a factor of 5 either way. The volume balance below is the tight check.
DEPTH_BAND = (0.2, 5.0)
VOLUME_RTOL = 0.02
CLAMP_FRAC = 1.0e-3       # clamped volume must be negligible vs the supply

CASES = [
    ("upslope-fed", "out/upslope-fed.h5"),
    ("outlet-fed", "out/outlet-fed.h5"),
]


def check(name, path):
    """Return a list of failure strings for one case ([] = pass)."""
    fails = []
    with h5py.File(HERE / path, "r") as f:
        audit = f["/monitor/mass_audit"]
        # Read by name: the column set depends on which modules are on (a
        # groundwater run carries an extra `seepage`), so positions are not
        # portable even though this case always has gw off.
        cols = {n: k for k, n in enumerate(audit.attrs["columns"].decode().split(","))}
        ma = audit[:]
        snaps = sorted(int(t) for t in f["/surface/depth"].keys())
        depth = f[f"/surface/depth/{snaps[-1]}"][:].ravel()
        bottom = f["/grid/bottom"][:].ravel()

    time = ma[:, cols["time"]]
    outflow = ma[:, cols["boundary_outflow"]]
    t_end = float(time[-1])
    bout = float(outflow[-1])
    bc_in = float(ma[-1, cols["bc_inflow"]])
    clamped = float(ma[-1, cols["clamped"]])
    bed_step = float(bottom[1] - bottom[0])
    outlet_depth = float(depth[0])

    print(f"  {name}")
    print(f"    outlet depth   : {outlet_depth:.6e} m   (normal {NORMAL_DEPTH:.6e} m, "
          f"bed step {bed_step:.3f} m)")
    print(f"    outflow / in   : {bout:.6f} / {bc_in:.6f} m3")
    print(f"    clamped        : {clamped:.3e} m3")

    # 1. The outlet must sit at normal depth, not at the upslope sill.
    ratio = outlet_depth / NORMAL_DEPTH
    if not DEPTH_BAND[0] <= ratio <= DEPTH_BAND[1]:
        extra = ""
        if abs(outlet_depth - bed_step) < 0.25 * bed_step:
            extra = (f" -- this is the bed step {bed_step:.3f} m, i.e. the outlet is "
                     "held at the upslope sill (V2-A11)")
        fails.append(f"{name}: outlet depth {outlet_depth:.4e} m is {ratio:.1f}x "
                     f"normal depth {NORMAL_DEPTH:.4e} m{extra}")

    # 2. What went in must come out, once the film is subtracted.
    if not np.isclose(bout, bc_in, rtol=VOLUME_RTOL, atol=1e-3):
        sign = "over" if bout > bc_in else "under"
        fails.append(f"{name}: {sign}-drained -- {bout:.4f} m3 left against "
                     f"{bc_in:.4f} m3 injected ({(bout - bc_in) / bc_in:+.1%})")

    # 3. The below-bed clamp must not be making up the difference.
    if clamped > CLAMP_FRAC * SUPPLY:
        fails.append(f"{name}: below-bed clamp minted {clamped:.4e} m3 = "
                     f"{clamped / bc_in:.1%} of the injected volume")

    # 4. A transmissive outlet fed directly must pass water immediately; a
    #    long dead interval is the shut-face signature.
    first_flow = int(np.argmax(outflow > 1e-12))
    if outflow[first_flow] > 1e-12 and time[first_flow] > 0.05 * t_end:
        fails.append(f"{name}: no outflow at all until t = {time[first_flow]:.0f} s "
                     f"({time[first_flow] / t_end:.0%} of the run)")

    return fails


def main():
    print(f"swe-outflow-staircase -- Manning normal depth {NORMAL_DEPTH:.6e} m\n")
    missing = [p for _, p in CASES if not (HERE / p).exists()]
    if missing:
        raise SystemExit(f"missing output: {', '.join(missing)} -- run both yaml files first")

    fails = [msg for name, path in CASES for msg in check(name, path)]

    print()
    if fails:
        print(f"FAIL ({len(fails)}):")
        for msg in fails:
            print(f"  - {msg}")
        print("\nExpected against the code as shipped: this is the merged-failing "
              "gate for the plan §8.2 Outflow x {W,S} matrix cells (V2-A11).")
        raise SystemExit(1)
    print("PASS: both feeding directions settle at normal depth and conserve volume.")


if __name__ == "__main__":
    main()
