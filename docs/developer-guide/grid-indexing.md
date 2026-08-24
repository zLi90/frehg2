# Grid, indexing, and global IDs

`core::Grid` owns the global extents (NX, NY, NZ), the spacings (dx, dy,
dz(k) with geometric stretch `dz_stretch`), the MPI decomposition, and the
compressed global IDs the linear systems and the output layer share.

## Decomposition

A **2D Cartesian decomposition over (i, j) only** — every rank owns full
z-columns, because the column is the coupling and moisture-reallocation
unit (the reallocation sweep is sequential per column; a z-decomposition
would put MPI inside the PCA corrector). Non-divisible extents are
supported: `nx_local = NX/Px + (rank_i < NX%Px)` — no divisibility
assumption anywhere (unit-tested with NX = 101 on 4 ranks, b5's extent).

Rank layout: `pi = rank % px`, `pj = rank / px`. The automatic layout uses
`MPI_Dims_create` with the larger rank count assigned to the larger extent;
`domain.decomposition.mpi_nx/mpi_ny` pins it explicitly.

## Field layout (binding convention)

Fields are `(nyl+2, nxl+2[, NZ])` `LayoutRight` Kokkos Views:

- index order `(j, i[, k])`, interior local indices `1..nyLocal` /
  `1..nxLocal`;
- halo width 1 in i and j, **no halo in k**;
- **k unit-stride**, `k = 0` at the land surface increasing downward —
  matching the column-sweep access pattern, the vertical-physics inner
  loops, and the pinned output order;
- flattening for files *and* HDF5: 2D `j·NX + i`, 3D `(j·NX + i)·NZ + k`
  (pinned by the b3 benchmark's published plotting script; asserted
  element-wise in tests).

## Active-cell masking

One source of truth: `ktop(j, i)` (the legacy `istop`). Cells
`k ∈ [ktop, NZ)` are active; `ktop == nz` means the column is fully
inactive in both the 2D and 3D systems. Kernels early-return on inactive
cells; the full box is allocated (accepted at benchmark scales). Dry
surface cells keep their matrix row with an identity closure, so the
sparsity pattern never changes.

## Compressed global IDs

`Grid::buildGlobalIds(ktop)` builds `gid2` (surface) and `gid3`
(subsurface): `PetscInt`, −1 for inactive cells, ordered **by rank-owned
blocks, then (j, i, k) with k innermost** — so PETSc row ownership matches
the decomposition and the §7 output flattening needs no transpose. Halo
gids are exchanged once at construction.

The intended assembly pattern (used by both physics modules): build COO
row/col arrays once from the gids (5-point surface stencil / 7-point
subsurface stencil; out-of-domain or inactive legs get column −1 with the
coefficient folded into the diagonal or RHS), then per solve fill the
values View in a Kokkos kernel and call
`LinearSystem::setValues` + `solve`. Duplicate COO indices are summed by
PETSc — symmetric scatter is expressed by writing each face coefficient to
both rows.

## Geometry

- Cell centers: `x = (i + 0.5)·dx`, `y = (j + 0.5)·dy` — polygons rasterize
  against these.
- Vertical spacing: `dz·dz_stretch^k`. With `domain.follow_terrain: true`
  the mesh follows the local bed; `terrain_layers: scaled` scales the
  column profile onto each column's depth (the legacy wedge rule),
  `terrain_layers: uniform` stacks the configured profile below each local
  bed (amendment A12; the b5 slab geometry).
- Raster inputs bind **first data row = j = 0** (file order, *not* ESRI
  north-up) — chosen because the legacy reference inputs were consumed that
  way; documented in `GridDataReader.hpp`.
