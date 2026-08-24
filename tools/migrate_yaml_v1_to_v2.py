#!/usr/bin/env python3
"""migrate_yaml_v1_to_v2.py - upgrade v1/experimental Frehg2 YAML drafts to
the frozen v2 schema (upgrade plan Section 5.4 and Section 6).

Usage:
    migrate_yaml_v1_to_v2.py INPUT.v1.yaml [-o OUTPUT.yaml]

Rename table (v1 key -> v2 key). This table is the authoritative migration
contract referenced by the plan ("per the rename table in its own header").

  simulation.id                       -> simulation.id
  simulation.title                    -> simulation.title
  simulation.mode                     -> (dropped: redundant with modules.*)
  grid.nx/ny/nz/dx/dy/dz              -> domain.nx/ny/nz/dx/dy/dz
  grid.dz_incre                       -> domain.dz_stretch
  grid.bot_z                          -> domain.bottom_elevation.constant
                                         (unless bathymetry.from_file)
  grid.bathymetry.from_file: true     -> domain.bottom_elevation.file:
                                         input/bathymetry.dat (convention)
  grid.follow_terrain                 -> domain.follow_terrain
  time.dt                             -> time.dt
  time.Tend                           -> time.t_end
  time.dt_out                         -> time.output_interval
  time.max_step                       -> (dropped: derived from legacy NT,
                                         which legacy Frehg parsed but never
                                         used as a step cap)
  modules.surface_water/groundwater   -> modules.surface_water/groundwater
  modules.solute                      -> modules.transport
  surface_water.gravity               -> surface_water.gravity
  surface_water.manning               -> surface_water.friction.coefficient.
                                         constant (+ friction.law: manning)
  surface_water.h_diffusion_ref       -> surface_water.friction.
                                         thin_layer_depth (legacy hD)
  surface_water.min_depth             -> surface_water.min_depth
  surface_water.waterfall_depth       -> surface_water.wetting_face_depth
                                         (legacy wtfh)
  surface_water.viscosity.x/y         -> surface_water.viscosity.x/y
  sources.surface.rainfall.from_file  -> surface_water.rainfall.series.file:
                                         input/rain.dat (convention)
  sources.surface.rainfall.rate       -> surface_water.rainfall.constant
                                         (when not from_file)
  groundwater.solver                  -> groundwater.scheme
  groundwater.full_3d                 -> groundwater.use_full3d
  groundwater.dt_min/dt_max           -> groundwater.timestep.dt_min/dt_max
                                         (+ dt_init = dt_min)
  groundwater.co_max                  -> groundwater.timestep.courant_max
  groundwater.specific_storage        -> groundwater.specific_storage
  groundwater.air_entry_value         -> soil.types[].aev
  groundwater.adaptive_dt             -> (dropped: dtg adaptivity is always
                                         on; plan Section 3.1)
  groundwater.use_corrector           -> (dropped: PCA always runs the
                                         corrector; plan Appendix A)
  groundwater.post_allocate           -> (dropped: post-allocation is always
                                         on; plan Appendix A)
  groundwater.use_vg                  -> (dropped: van Genuchten is the only
                                         retention model; plan Section 3.2)
  groundwater.bc_type_gw              -> (dropped: face codes are replaced by
                                         polygon boundary_conditions;
                                         plan Section 5.6)
  soil.uniform (+theta_s/theta_r/
       alpha/n/ksat)                  -> soil.types[0] {name: soil0,
                                         ksx=ksy=ksz=ksat, theta_s, theta_r,
                                         vg_alpha, vg_n, aev} and
                                         soil.map.constant: soil0
  initial_conditions.surface_water.eta-> initial_conditions.surface.eta.
                                         constant
  initial_conditions.groundwater.wc   -> initial_conditions.groundwater.
                                         moisture.constant
  io.dir                              -> output.filename: <dir>/output.h5
  output.format                       -> (dropped: HDF5 is the only format)
  output.variables[water_depth]       -> output.variables.surface[depth]
  output.variables[eta|uu|vv|seepage] -> output.variables.surface[...]
  output.variables[hydraulic_head|
       water_content|qx|qy|qz]        -> output.variables.groundwater[...]
  bc[].type: discharge|bc_discharge   -> boundary_conditions[].kind:
                                         discharge

Anything in the v1 file not covered above is an error: the tool refuses to
guess, so silent key loss is impossible.
"""

from __future__ import annotations

import argparse
import sys

import yaml


class MigrationError(Exception):
    """A v1 construct the rename table does not cover."""


def _take(mapping, key, default=None):
    """Pop a key from a dict-like YAML node; None if absent."""
    if mapping is None:
        return default
    return mapping.pop(key, default)


def _require_empty(mapping, where):
    if mapping:
        raise MigrationError(
            f"unrecognized v1 key(s) under '{where}': {sorted(mapping)}; "
            "the rename table in the tool header does not cover them")


def migrate(v1: dict) -> dict:
    """Translate a parsed v1 document into a v2 document (dict)."""
    v1 = dict(v1)  # shallow copy so pops do not mutate the caller's data
    v2: dict = {}
    notes: list[str] = []

    # -- simulation ---------------------------------------------------------
    sim = dict(_take(v1, "simulation") or {})
    v2["simulation"] = {"id": _take(sim, "id", "unnamed")}
    title = _take(sim, "title")
    if title is not None:
        v2["simulation"]["title"] = title
    if _take(sim, "mode") is not None:
        notes.append("simulation.mode dropped (redundant with modules.*)")
    _require_empty(sim, "simulation")

    # -- grid -> domain -----------------------------------------------------
    grid = dict(_take(v1, "grid") or {})
    domain = {}
    for key in ("nx", "ny", "nz", "dx", "dy", "dz"):
        value = _take(grid, key)
        if value is None:
            raise MigrationError(f"grid.{key} is required in a v1 file")
        domain[key] = value
    stretch = _take(grid, "dz_incre")
    if stretch is not None:
        domain["dz_stretch"] = stretch
    bathymetry = dict(_take(grid, "bathymetry") or {})
    bot_z = _take(grid, "bot_z")
    if bathymetry.pop("from_file", False):
        domain["bottom_elevation"] = {"file": "input/bathymetry.dat"}
        if bot_z is not None:
            notes.append("grid.bot_z dropped in favor of the bathymetry file")
    else:
        domain["bottom_elevation"] = {"constant": bot_z if bot_z is not None else 0.0}
    _require_empty(bathymetry, "grid.bathymetry")
    follow = _take(grid, "follow_terrain")
    if follow is not None:
        domain["follow_terrain"] = follow
    _require_empty(grid, "grid")
    v2["domain"] = domain

    # -- time ---------------------------------------------------------------
    time_v1 = dict(_take(v1, "time") or {})
    time_v2 = {"dt": _take(time_v1, "dt")}
    time_v2["t_end"] = _take(time_v1, "Tend")
    time_v2["output_interval"] = _take(time_v1, "dt_out")
    if _take(time_v1, "max_step") is not None:
        notes.append("time.max_step dropped (legacy NT was parsed but unused)")
    _require_empty(time_v1, "time")
    v2["time"] = time_v2

    # -- modules ------------------------------------------------------------
    modules_v1 = dict(_take(v1, "modules") or {})
    modules_v2 = {
        "surface_water": bool(_take(modules_v1, "surface_water", False)),
        "groundwater": bool(_take(modules_v1, "groundwater", False)),
        "transport": bool(_take(modules_v1, "solute", False)),
    }
    _require_empty(modules_v1, "modules")
    v2["modules"] = modules_v2

    # -- surface_water ------------------------------------------------------
    sw = dict(_take(v1, "surface_water") or {})
    if modules_v2["surface_water"] or sw:
        out = {}
        gravity = _take(sw, "gravity")
        if gravity is not None:
            out["gravity"] = gravity
        friction = {"law": "manning",
                    "coefficient": {"constant": _take(sw, "manning", 0.0)}}
        thin = _take(sw, "h_diffusion_ref")
        if thin is not None:
            friction["thin_layer_depth"] = thin
        out["friction"] = friction
        viscosity = _take(sw, "viscosity")
        if viscosity is not None:
            out["viscosity"] = viscosity
        out["min_depth"] = _take(sw, "min_depth")
        out["wetting_face_depth"] = _take(sw, "waterfall_depth")
        _require_empty(sw, "surface_water")
        v2["surface_water"] = out

    # -- sources -> surface_water.rainfall ---------------------------------
    sources = dict(_take(v1, "sources") or {})
    surface_sources = dict(sources.pop("surface", {}) or {})
    rainfall = dict(surface_sources.pop("rainfall", {}) or {})
    if rainfall:
        if rainfall.pop("from_file", False):
            rainfall.pop("rate", None)
            v2.setdefault("surface_water", {})["rainfall"] = {
                "series": {"file": "input/rain.dat"}}
        else:
            v2.setdefault("surface_water", {})["rainfall"] = {
                "constant": rainfall.pop("rate", 0.0)}
        _require_empty(rainfall, "sources.surface.rainfall")
    _require_empty(surface_sources, "sources.surface")
    _require_empty(sources, "sources")

    # -- groundwater --------------------------------------------------------
    gw = dict(_take(v1, "groundwater") or {})
    aev = _take(gw, "air_entry_value", 0.0)
    if modules_v2["groundwater"] or gw:
        out = {}
        scheme = _take(gw, "solver")
        if scheme is not None:
            out["scheme"] = scheme
        full3d = _take(gw, "full_3d")
        if full3d is not None:
            out["use_full3d"] = full3d
        dt_min = _take(gw, "dt_min")
        dt_max = _take(gw, "dt_max")
        timestep = {"dt_init": dt_min, "dt_min": dt_min, "dt_max": dt_max}
        co_max = _take(gw, "co_max")
        if co_max is not None:
            timestep["courant_max"] = float(co_max)
        out["timestep"] = timestep
        out["specific_storage"] = _take(gw, "specific_storage")
        for dropped, why in (
            ("adaptive_dt", "dtg adaptivity is always on"),
            ("use_corrector", "the PCA corrector always runs"),
            ("post_allocate", "post-allocation is always on"),
            ("use_vg", "van Genuchten is the only retention model"),
            ("use_mvg", "modified van Genuchten was removed (plan 3.2)"),
            ("bc_type_gw", "face codes replaced by polygon boundary_conditions"),
        ):
            if _take(gw, dropped) is not None:
                notes.append(f"groundwater.{dropped} dropped ({why})")
        _require_empty(gw, "groundwater")
        v2["groundwater"] = out

    # -- soil ----------------------------------------------------------------
    soil = dict(_take(v1, "soil") or {})
    if soil:
        if not soil.pop("uniform", False):
            raise MigrationError("only soil.uniform v1 files exist; refusing to guess")
        ksat = soil.pop("ksat")
        soil_type = {
            "name": "soil0",
            "ksx": ksat, "ksy": ksat, "ksz": ksat,
            "theta_s": soil.pop("theta_s"),
            "theta_r": soil.pop("theta_r"),
            "vg_alpha": soil.pop("alpha"),
            "vg_n": soil.pop("n"),
            "aev": aev,
        }
        _require_empty(soil, "soil")
        v2["soil"] = {"types": [soil_type], "map": {"constant": "soil0"}}

    # -- initial conditions --------------------------------------------------
    ic = dict(_take(v1, "initial_conditions") or {})
    ic_v2 = {}
    ic_sw = dict(ic.pop("surface_water", {}) or {})
    if ic_sw:
        ic_v2["surface"] = {"eta": {"constant": ic_sw.pop("eta")}}
        _require_empty(ic_sw, "initial_conditions.surface_water")
    ic_gw = dict(ic.pop("groundwater", {}) or {})
    if ic_gw:
        ic_v2["groundwater"] = {"moisture": {"constant": ic_gw.pop("wc")}}
        _require_empty(ic_gw, "initial_conditions.groundwater")
    _require_empty(ic, "initial_conditions")
    if ic_v2:
        v2["initial_conditions"] = ic_v2

    # -- boundary conditions -------------------------------------------------
    bcs = _take(v1, "bc") or _take(v1, "boundary_conditions") or []
    bcs_v2 = []
    for entry in bcs:
        entry = dict(entry)
        kind = entry.pop("type")
        if kind in ("discharge", "bc_discharge"):
            kind = "discharge"
        bcs_v2.append({
            "name": entry.pop("name"),
            "region": {"polygon": entry.pop("polygon")},
            "target": entry.pop("target", "surface"),
            "kind": kind,
            "value": {"constant": entry.pop("value")},
        })
        _require_empty(entry, "bc[]")
    if bcs_v2:
        v2["boundary_conditions"] = bcs_v2

    # -- io / output ---------------------------------------------------------
    io_v1 = dict(_take(v1, "io") or {})
    out_v1 = dict(_take(v1, "output") or {})
    out_dir = _take(io_v1, "dir", "out")
    _require_empty(io_v1, "io")
    if _take(out_v1, "format") is not None:
        notes.append("output.format dropped (HDF5 is the only format)")
    surface_vars = []
    gw_vars = []
    var_map = {"water_depth": ("surface", "depth"),
               "eta": ("surface", "eta"),
               "uu": ("surface", "uu"),
               "vv": ("surface", "vv"),
               "seepage": ("surface", "seepage"),
               "hydraulic_head": ("groundwater", "hydraulic_head"),
               "water_content": ("groundwater", "water_content"),
               "qx": ("groundwater", "qx"),
               "qy": ("groundwater", "qy"),
               "qz": ("groundwater", "qz")}
    for var in _take(out_v1, "variables", []) or []:
        if var not in var_map:
            raise MigrationError(f"output variable '{var}' has no v2 equivalent")
        group, renamed = var_map[var]
        (surface_vars if group == "surface" else gw_vars).append(renamed)
    _require_empty(out_v1, "output")
    output_v2 = {"filename": f"{out_dir}/output.h5"}
    variables = {}
    if surface_vars:
        variables["surface"] = surface_vars
    if gw_vars:
        variables["groundwater"] = gw_vars
    if variables:
        output_v2["variables"] = variables
    v2["output"] = output_v2

    _require_empty(v1, "(document root)")
    v2["_migration_notes"] = notes  # stripped before dumping; used by tests
    return v2


def render(v2: dict) -> str:
    """Serialize the migrated document with a provenance header."""
    notes = v2.pop("_migration_notes", [])
    header = ["# Migrated to the frozen v2 schema by tools/migrate_yaml_v1_to_v2.py."]
    for note in notes:
        header.append(f"# note: {note}")
    body = yaml.safe_dump(v2, sort_keys=False, default_flow_style=False)
    return "\n".join(header) + "\n" + body


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", help="v1 YAML file")
    parser.add_argument("-o", "--output", help="v2 YAML output (default: stdout)")
    args = parser.parse_args(argv)

    with open(args.input, "r", encoding="utf-8") as handle:
        v1 = yaml.safe_load(handle)
    try:
        v2 = migrate(v1)
    except MigrationError as error:
        print(f"migration failed: {error}", file=sys.stderr)
        return 1
    text = render(v2)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
