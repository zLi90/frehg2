"""gate_ws.py — pure reference and check functions for the Stage-3 T9 paired
water–salt store gates (plan §12A, ws0–ws5). No HDF5 I/O: run_regression.py
reads the files and passes arrays; scripts/test_ws_gates.py feeds corrupted
synthetic inputs to every checker (the §6.3 negative battery).

Column conventions (Monitor tables carry time in column 0):
  mass_audit       time, volume, rain, evaporation, boundary_outflow,
                   bc_inflow, [seepage], clamped, [stage_deposit],
                   film_removed, film_pending, stage_volume,
                   [store_volume, store_release]
  transport_audit  time, surf_mass, subs_mass, exchange, surf_source,
                   surf_boundary, surf_adjust, surf_anchor, subs_boundary,
                   subs_adjust, subs_anchor, [surf_store]
"""

from __future__ import annotations

import numpy as np


# --- closed forms -----------------------------------------------------------

def mix_release(mass: float, volume: float, store_mass: float,
                store_volume: float) -> float:
    """R5: the release mixes, c = (M + M_s) / (V + V_s)."""
    return (mass + store_mass) / (volume + store_volume)


def pond_rewet_concentration(c0: float, h0: float, area: float,
                             v_parked: float, h_rain: float) -> float:
    """et8b / ws3(c): the first wet step after a pond dried with a parked
    sub-min_depth remainder V_r rains h_rain:
    s = c0 h0 A / (V_r + h_rain A)  (S3-A12)."""
    return c0 * h0 * area / (v_parked + h_rain * area)


# --- identities -------------------------------------------------------------

def audit_columns(coupled: bool, store: bool) -> dict[str, int]:
    """Column index map of /monitor/mass_audit at T9 for a run's layout."""
    names = ["time", "volume", "rain", "evaporation", "boundary_outflow",
             "bc_inflow"]
    if coupled:
        names.append("seepage")
    names.append("clamped")
    if coupled:
        names.append("stage_deposit")
    names += ["film_removed", "film_pending", "stage_volume"]
    if store:
        names += ["store_volume", "store_release"]
    return {n: k for k, n in enumerate(names)}


def water_identity_residuals(rows: np.ndarray, cols: dict[str, int],
                             store: bool) -> np.ndarray:
    """Per-row residual of the mode's surface water identity (plan §12A):

    discard: Δ(volume − stage_volume + film_pending) = rain − evaporation − outflow
             + bc_inflow + seepage + clamped − stage_deposit − film_removed
    store:   Δ(volume − stage_volume + store_volume) = rain − evaporation − outflow
             + bc_inflow + seepage + clamped − stage_deposit

    (cumulative columns against row 0; film_pending is identically zero
    in store mode after the W9.2 sweep, and is carried either way).
    Returns lhs − rhs for every row."""
    def col(name):
        return rows[:, cols[name]] if name in cols else np.zeros(len(rows))

    def delta(name):
        c = col(name)
        return c - c[0]

    lhs = delta("volume") - delta("stage_volume") + delta("film_pending")
    rhs = (delta("rain") - delta("evaporation") - delta("boundary_outflow")
           + delta("bc_inflow") + delta("seepage") + delta("clamped")
           - delta("stage_deposit"))
    if store:
        lhs = lhs + delta("store_volume")
    else:
        rhs = rhs - delta("film_removed")
    return lhs - rhs


def water_identity_scale(rows: np.ndarray, cols: dict[str, int]) -> float:
    """The ws2 normaliser: cumulative |stage exchange| + rain + evaporation,
    the stage exchange being the BC face volumes plus the coupled exchange
    and the stage deposits (all final-row cumulative magnitudes)."""
    last, first = rows[-1], rows[0]
    total = 0.0
    for name in ("boundary_outflow", "bc_inflow", "seepage", "stage_deposit",
                 "rain", "evaporation"):
        if name in cols:
            total += abs(float(last[cols[name]] - first[cols[name]]))
    return total


def check_water_identity(rows: np.ndarray, cols: dict[str, int], store: bool,
                         tol_rel: float, abs_floor: float = 0.0) -> tuple[bool, str]:
    """|res| <= max(tol_rel * scale, abs_floor) at every row (S3-A13: the
    absolute floor is the coupled lanes' solve-residual class, ~3e-11 m^3,
    identical in the store and discard twins)."""
    res = water_identity_residuals(rows, cols, store)
    scale = water_identity_scale(rows, cols)
    worst = float(np.max(np.abs(res))) if len(res) else 0.0
    rel = worst / scale if scale > 0.0 else worst
    bound = max(tol_rel * scale, abs_floor)
    ok = worst <= bound
    mode = "store" if store else "discard"
    return ok, (f"water identity [{mode}] over {len(res)} rows: worst "
                f"|res| = {worst:.3e} m^3 = {rel:.3e} of scale {scale:.3e} "
                f"(<= max({tol_rel:.0e} X, {abs_floor:.0e} m^3))")


def salt_identity_residuals(audit: np.ndarray, has_store: bool) -> np.ndarray:
    """Δ(surf_mass + subs_mass [+ surf_store]) − Σ cumulative terms per row."""
    mass = audit[:, 1] + audit[:, 2]
    if has_store:
        mass = mass + audit[:, 11]
    lhs = mass - mass[0]
    terms = audit[:, 4:11].sum(axis=1)
    rhs = terms - terms[0]
    return lhs - rhs


def check_salt_identity(audit: np.ndarray, has_store: bool,
                        tol_rel: float) -> tuple[bool, str]:
    res = salt_identity_residuals(audit, has_store)
    mass0 = float(audit[0, 1] + audit[0, 2] + (audit[0, 11] if has_store else 0.0))
    scale = max(abs(mass0), float(np.max(np.abs(audit[:, 1] + audit[:, 2]))))
    worst = float(np.max(np.abs(res)))
    rel = worst / scale if scale > 0.0 else worst
    ok = rel <= tol_rel
    return ok, (f"salt identity over {len(res)} rows: worst {rel:.3e} of "
                f"M {scale:.3e} (<= {tol_rel:.0e})")


# --- bounds -----------------------------------------------------------------

def check_source_bound(maxima: list[float], c_max: float,
                       rel: float, label: str) -> tuple[bool, str]:
    """ws1 (1): every frame's maximum stays at or below c_max (1 + rel)."""
    worst = float(max(maxima)) if maxima else 0.0
    bound = c_max * (1.0 + rel)
    ok = worst <= bound
    return ok, (f"{label}: max over {len(maxima)} frames = {worst:.10g} "
                f"(<= {c_max} (1 + {rel:.0e}))")


def check_store_bound(store_mass: np.ndarray, store_volume: np.ndarray,
                      c_max: float, rel: float) -> tuple[bool, str]:
    """ws1 (3) / I2: M_s <= c_max V_s per cell (no salt without its water
    beyond the source bound; cells with V_s = 0 must hold no salt)."""
    excess = store_mass - c_max * store_volume * (1.0 + rel)
    # Absolute floor for rounding on empty cells: c_max * eps-scaled volume.
    floor = c_max * float(np.max(store_volume)) * 1.0e-15 if store_volume.size else 0.0
    worst = float(np.max(excess)) if excess.size else 0.0
    ok = worst <= floor
    return ok, (f"store bound M_s <= {c_max} V_s: worst excess "
                f"{worst:.3e} psu m^3 (<= {floor:.1e})")


def check_paired(store_mass: np.ndarray, store_volume: np.ndarray,
                 conc: float, rel: float) -> tuple[bool, str]:
    """ws4: on a uniform-concentration lane every parked cell carries
    exactly its water's salt, M_s = c V_s (1e-12)."""
    mask = store_volume > 0.0
    if not np.any(mask):
        return False, "pairing: no cell holds parked water (vacuous)"
    err = np.abs(store_mass[mask] - conc * store_volume[mask])
    worst = float(np.max(err / (conc * store_volume[mask])))
    ok = worst <= rel
    return ok, (f"pairing M_s = {conc} V_s on {int(mask.sum())} cells: "
                f"worst rel {worst:.3e} (<= {rel:.0e})")


def check_monotone_while_dry(store_volume_series: np.ndarray,
                             dry_series: np.ndarray) -> tuple[bool, str]:
    """ws3(b)/(d), R2/R7: V_s never decreases between two frames in which
    the cell stays dry (evaporation never debits the store)."""
    drops = 0
    worst = 0.0
    for k in range(1, len(store_volume_series)):
        both_dry = dry_series[k - 1] & dry_series[k]
        d = store_volume_series[k] - store_volume_series[k - 1]
        bad = both_dry & (d < 0.0)
        if np.any(bad):
            drops += int(np.sum(bad))
            worst = min(worst, float(np.min(d[bad])))
    ok = drops == 0
    return ok, (f"store volume monotone while dry: {drops} decreases "
                f"(worst {worst:.3e} m^3)")
