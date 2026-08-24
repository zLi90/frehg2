# Tutorial: adding an output variable

Field output is driven by the `output.variables.{surface,groundwater,transport}`
lists; the driver maps each configured name to a device view and writes it
at every output time through `io/Hdf5Output` (mask-gather → host mirror →
collective hyperslab write). Adding a variable touches four places — all in
the same PR (plan §11.1: schema, docs, and a test land together).

## 1. The schema (`src/core/ConfigSchema.cpp`)

Add the name to the allowed-values list of the grid it belongs to (search
for the existing `{"eta", "depth", "uu", "vv", "seepage"}` list and its
groundwater/transport siblings). If the variable only exists in some
configurations (like `seepage`, which is produced by the coupler), extend
the cross-field check that rejects it elsewhere — an impossible request
must fail validation, not produce silent zeros.

## 2. The driver mapping (`src/driver/Simulation.cpp`)

Two spots:

- the units/long-name helper at the top of the file — every dataset carries
  `units`, `long_name`, and `time` attributes (the HDF5 contract, plan §7);
- the output switch that resolves a configured name to a field view (search
  for `var == "seepage"`). Surface variables resolve to `Field2` views in
  `j·NX + i` layout; subsurface variables to `Field3` in
  `(j·NX + i)·NZ + k` with NaN above `ktop` — both handled by the shared
  write path, so the mapping is one `else if`.

If the quantity does not already exist as a persistent field (e.g. a
diagnostic that must be derived), compute it into a driver-owned scratch
view at output time — do not add state to a physics module for output's
sake.

## 3. Documentation

- `docs/user-guide/parameters.md`: the allowed-values row
  (`scripts/check_parameter_docs.py` fails the build on drift);
- `docs/user-guide/output.md`: the variable's meaning, units, and grid.

## 4. Tests

- Schema acceptance/rejection in `tests/unit/test_config.cpp` (accepted on
  the right grid, rejected where it cannot exist);
- a layout/attribute assertion in the HDF5 tests or the owning module's
  test if the variable carries new semantics.

## Things that are easy to get wrong

- **Restart is not output.** If the new quantity is *prognostic* (the next
  step reads it and it cannot be rebuilt from other checkpointed state), it
  must also join the checkpoint field list — see how the coupler's seepage
  accumulator (`seep_accum`) and the transport flow-rate snapshots
  (`s_fu_old`/`s_fv_old`) are handled — and the bitwise restart gates will
  tell you if you forgot.
- **Output times are integer seconds** (`/surface/<var>/<t>` with
  `<t> = str(int(round(t)))`); nothing about a new variable may bend that.
- **Monitors are separate.** Point monitors (`output.monitors`) have their
  own variable resolution in the driver; add the name there too if
  monitoring it per-step should work.
