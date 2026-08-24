# Tutorial: adding a boundary-condition kind

Boundary conditions are split into **definition** (YAML → validated config),
**rasterization** (`bc/BoundarySet` selects member cells/faces once at
startup), and **application** (the owning physics module reads the member
list each step). Adding a kind touches every layer — and, per the plan's
binding rules (§11.1), the schema, the docs, and a test land **in the same
PR**.

Walkthrough, using a hypothetical surface kind `wall_overtopping` as the
example:

## 1. The enum

Add the kind to `BcKind` in `src/core/Config.hpp`:

```cpp
enum class BcKind { Eta, Discharge, Velocity, Outflow, Head, Flux,
                    ScalarValue, WallOvertopping };
```

## 2. The schema (`src/core/ConfigSchema.cpp`)

- Add the YAML spelling to the allowed-values list of
  `boundary_conditions[].kind` (search for the existing list
  `{"eta", "discharge", ...}`).
- Extend the cross-field checks: which `target`s the kind is valid for,
  whether `value` is required (see how `outflow` forbids a value and every
  other kind requires one), and any module dependency (e.g. `scalar_value`
  requires transport).
- Map the string to the enum in the config materialization.

Unknown-key rejection and nearest-key suggestions come for free from the
schema tree.

## 3. Rasterization semantics (`src/bc/BoundarySet`)

Most kinds reuse the existing member selection: footprint targets
(`surface`, `groundwater_top/bottom`) select every owned cell whose center
lies inside the polygon; `groundwater_side` selects domain-edge cells and
tags the edge face(s). Only touch `BoundarySet` if the new kind needs a
different member rule — and if so, extend `tests/unit/test_polygon.cpp`
with the new geometry cases, including a polygon spanning a rank boundary.

Per-BC MPI sub-communicators (for cross-section reductions like
distributing a total discharge) already exist on every rasterized BC.

## 4. Application in the owning module

Surface kinds are applied in `src/swe` (see how `eta` writes the stage and
its edge-face velocity correction, and how `discharge` distributes the
inflow over the member cross-section via the sub-communicator); subsurface
kinds in `src/gw` (ghost head / face flux in the predictor and the Darcy
flux); scalar kinds in `src/transport`. Follow the existing kind switch —
the modules read the `BoundarySet` member lists, never the YAML.

Keep the application in a Kokkos kernel free of host-pointer capture (the
CUDA compile lane rejects host lambdas capturing host memory).

## 5. Documentation and tests (same PR — enforced)

- `docs/user-guide/parameters.md`: the kind's row in the
  boundary-conditions table. `scripts/check_parameter_docs.py` (ctest
  `unit.parameter_docs`) fails the build if the schema and this file
  drift.
- A unit test exercising the new branch (`tests/unit/test_config.cpp` for
  schema acceptance/rejection; the owning module's test file for the
  physics), plus — if the kind changes a gated benchmark's configuration —
  the regression tolerance record.
- If the kind replaces or generalizes legacy behavior, add the provenance
  note in the relevant `docs/theory/*.md` table, and a row in
  `docs/theory/removed-features.md` if legacy behavior is dropped.

## Precedents to copy from

- `outflow` (amendment A2) — a kind with **no value**, edge-face-only
  semantics, added for the b4 gate.
- `scalar_value` (amendments A19/A20) — a kind whose meaning depends on
  what else covers the region (discharge-paired inflow concentration vs
  wet-cell Dirichlet), with a schema cross-check requiring transport.
- `flux: {gravity: true}` — a kind with a structured value form (free
  drainage), the pattern for non-scalar values.
