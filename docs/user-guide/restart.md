# Checkpoint and restart

Frehg2 can checkpoint its full prognostic state into the output file and
resume from any checkpoint. Restart is **bitwise deterministic**: the
restarted run reproduces the uninterrupted run to the last bit at every
output time. This is enforced by regression gates for surface-only (b1),
groundwater-only (b2), coupled (b5), and coupled-with-transport (b6) runs.

## Writing checkpoints

```yaml
output:
  filename: out/output.h5
  checkpoint:
    interval: 3600        # [s]; 0 = off
```

Checkpoints land in the output file under `/checkpoint/<t>/...`. One is
always written at `t_end`, even with `interval: 0`. Under adaptive stepping
(groundwater-only and sync-coupled runs) the step generally crosses the
checkpoint interval mid-step: the group key `<t>` is the crossed
whole-second boundary, while the group's `t` attribute and the stored
adaptive-step scalars carry the **exact** state time — restart resumes from
the true state, not the label.

The checkpoint carries everything the time loop needs and nothing it can
recompute, including:

- the surface prognostic fields (`eta`, velocities, the previous-step stage
  `etan`, and — coupled — the seepage accumulator `seep_accum` and the
  adaptive surface step `dt_surface`);
- the subsurface head and water content plus the adaptive step `dtg`;
- the transport scalars on both grids and the transport-specific carried
  state (`s_fu_old`/`s_fv_old` — the pre-stage-correction flow-rate
  snapshots — and `s_dzz_top`, the carried top-cell dispersion
  coefficient), which cannot be rebuilt from the flow state.

Time series (rainfall, tides, hydrographs) are evaluated statelessly by
time, so they need no checkpoint state.

## Resuming

```yaml
restart:
  enabled: true
  file: out/output.h5     # the file holding the checkpoint
  time: 9000              # the checkpoint key to resume from [s]
```

Point `restart.file` at the earlier run's output file and `restart.time` at
one of its `/checkpoint/<t>` keys. Write the resumed run's output to a
**different** `output.filename` unless you intend to extend the same file.

The restarted run must use the same configuration otherwise (same grid,
modules, soils, boundary conditions); the schema does not attempt to verify
equivalence, and changing physics keys across a restart is undefined
behavior in the modeling sense — the run will proceed, but it is a new
experiment, not a continuation.

## Verifying determinism

The regression suite gates it, and you can check any pair of runs yourself:
every field at every common output time must match bitwise:

```bash
ctest --test-dir build -R 'b1_restart|b2_restart|b5_restart|b6_restart'
```
