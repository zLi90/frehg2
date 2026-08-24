# b4-govindaraju — overland-flow hydrograph (kinematic-wave reference)

Rainfall-runoff on a 21.95 m × 1.10 m plane (200 × 10 cells, 4 % slope),
Chezy friction C = 1.767, two identical rain pulses (50.8 / 101.6 mm/h).
The reference (`legacy/benchmarks/b4-govindaraju/ReferenceData/outflow.txt`)
is a kinematic-wave solution; the gate is metric-based (plan §9 b4) and run
by `tests/regression/run_regression.py b4`.

Legacy Frehg never ran this case — the committed legacy inputs are SERGHEI
format. Two decisions were adjudicated by the P1 gate:

1. **Outlet = `kind: outflow`** (plan amendment A2). The P0 provisional
   mapping (stage held at −10 m) keeps the outlet column permanently dry —
   the hydrograph measured there is identically zero — and a plain open
   boundary retains ~half the rainfall because the closed-edge fold removes
   the gravity forcing at the boundary face (measured: −48 % outflow
   volume). The transmissive outflow kind extrapolates the stage down the
   continued bed slope at the east boundary faces, restoring the interior
   momentum balance at the outlet; with it the outflow volume closes to
   −0.1 % of rainfall.

2. **The model hydrograph is the mass-consistent outflow rate**
   (d(cum outflow)/dt from `/monitor/mass_audit`), not `uu·depth` at the
   outlet column. The preserved legacy velocity update applies the drag
   factor twice (`shallowwater.c:745-746`), so stored velocities understate
   the conserved face fluxes by ~20 % here; the volume budget is what the
   scheme actually discharges (docs/theory/surface-water.md). Time-to-peak
   is the first crossing of 99 % of the peak: the two identical rain pulses
   produce twin peaks that differ at the 0.1 % level (also true of the
   reference), so a bare argmax is degenerate.

Gate results at P1 (allowed in parentheses): rel-L2 0.082 (0.15), peak
+0.3 % (±10 %), time-to-99 %-peak −7.9 % (±10 %), outflow volume −0.10 % of
rainfall (5 %) and −0.22 % of the reference volume (3 %).
