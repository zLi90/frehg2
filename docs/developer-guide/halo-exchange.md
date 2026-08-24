# Halo exchange

`core::HaloExchanger` implements the one communication pattern every module
shares. Its contract:

## Registration and exchange

- Modules register fields once at init: `add(name, field)`. 2D and 3D
  fields coexist in one exchanger.
- `exchangeAll()` sends **one coalesced message per neighbor** covering all
  registered fields (persistent device buffers sized once at
  registration).
- `exchange({names})` is the targeted variant for mid-step updates (e.g.
  the stage after the free-surface solve).
- Protocol: post `MPI_Irecv` ×4 → pack kernels ×4 → `Kokkos::fence` →
  `MPI_Isend` ×4 → `Waitall` → unpack ×4. Message tags pair direction
  `dir ^ 1`.

## What the exchanger does *not* do

- **Domain-edge halos are never written by the exchanger** — the BC system
  owns them. A rank at the domain edge leaves those ghost slots to the
  boundary-condition application.
- No halo in k (full columns are rank-local by the decomposition).

## Corner exchange

The base protocol carries no corner (diagonal) ghosts — all 5/7-point
stencils are corner-free. Two consumers need diagonals and use the opt-in
`exchangeWithCorners()` (amendment A1):

- the SWE `uy`/`vx` four-point velocity interpolation (2D);
- the subsurface dispersion cross terms, which read diagonal scalar
  neighbors (3D — the original "no 3D corners" premise was repealed in P4,
  amendment A19).

`exchangeWithCorners()` runs two sequential phases: west/east with
interior-row packs, then south/north with full-width rows that carry the
already-settled i-halo columns — corner values arrive with no extra
messages and no diagonal neighbors in the communicator.

## GPU-aware MPI

The constructor takes a `gpuAware` flag — **the only place in the code the
toggle exists** (CMake default `FREHG_GPU_AWARE_MPI`, runtime-overridable
via `runtime.gpu_aware_mpi`). Off means the pack/unpack buffers deep-copy
through persistent host mirrors; on passes device pointers to MPI. Physics
code cannot tell the difference.

## Guarantees and tests

Exchanged ghosts are **bitwise-exact** copies of the neighbor's interior
(no arithmetic in the path); `tests/mpi/test_halo.cpp` asserts this with
analytic f(i, j, k) fields at 1/2/4 ranks in both staging modes, corners
included. If you add a stencil that reads a ghost the exchanger does not
deliver (e.g. a wider stencil), the rank-invariance regression lanes are
the net that catches it — the legacy code's silent fallback-to-upwind at
rank interfaces is exactly the class of bug they exist for.
