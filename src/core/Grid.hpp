/// \file Grid.hpp
/// \brief Structured grid, 2D block decomposition, and compressed global ids
///        (plan §5.2).
///
/// The decomposition is Cartesian over (i, j) only — every rank owns full
/// z-columns, because the column is the coupling and moisture-reallocation
/// unit. Non-divisible extents are supported with the block rule
/// \c n_local = N/P + (p < N%P).
///
/// Index conventions:
///  - fields are (j, i[, k]) with a one-cell halo in j and i; local interior
///    indices run 1..nyLocal / 1..nxLocal, halo at 0 and n+1;
///  - k runs 0..nz-1 from the land surface downward, unit stride;
///  - the cell (j, i) is active in the subsurface for k in [ktop(j,i), nz);
///    ktop == nz marks a fully inactive column (also inactive for the
///    surface solver, whose dry-but-active cells keep their matrix rows).
///
/// Global ids are compressed (inactive cells get -1) and ordered by
/// rank-owned blocks so PETSc row ownership matches the decomposition.

#ifndef FREHG_CORE_GRID_HPP
#define FREHG_CORE_GRID_HPP

#include "core/Config.hpp"
#include "core/Types.hpp"

#include <mpi.h>
#include <petscsystypes.h>

#include <array>
#include <vector>

namespace frehg {

/// Structured grid with 2D Cartesian block decomposition.
class Grid {
 public:
  /// Build the grid and decomposition.
  /// \param comm communicator over which the domain is decomposed.
  /// \param domain global extents, spacings, and decomposition request
  ///        (mpiNx/mpiNy of 0 means automatic via MPI_Dims_create).
  Grid(MPI_Comm comm, const DomainConfig& domain);

  /// \name Global extents and spacings
  ///@{
  int nx() const { return nx_; }              ///< global cells in i
  int ny() const { return ny_; }              ///< global cells in j
  int nz() const { return nz_; }              ///< global cells in k
  real_t dx() const { return dx_; }           ///< spacing in i [m]
  real_t dy() const { return dy_; }           ///< spacing in j [m]
  /// Layer thickness of level \p k [m]: dz * dz_stretch^k.
  real_t dzK(int k) const { return dz_[static_cast<std::size_t>(k)]; }
  /// Depth of the center of level \p k below the z = 0 reference [m]
  /// (negative values, decreasing with k).
  real_t zCenter(int k) const { return zCenter_[static_cast<std::size_t>(k)]; }
  /// Total subsurface thickness [m].
  real_t totalDepth() const { return totalDepth_; }
  /// x of the center of global column \p iGlobal [m].
  real_t xCenter(int iGlobal) const { return (static_cast<real_t>(iGlobal) + 0.5) * dx_; }
  /// y of the center of global row \p jGlobal [m].
  real_t yCenter(int jGlobal) const { return (static_cast<real_t>(jGlobal) + 0.5) * dy_; }
  ///@}

  /// \name Decomposition
  ///@{
  /// The communicator the domain is decomposed over.
  MPI_Comm comm() const { return comm_; }
  /// This process's rank in comm().
  int rank() const { return rank_; }
  /// Number of ranks in comm().
  int size() const { return size_; }
  int px() const { return px_; }              ///< ranks along i
  int py() const { return py_; }              ///< ranks along j
  int pi() const { return pi_; }              ///< this rank's block index along i
  int pj() const { return pj_; }              ///< this rank's block index along j
  int nxLocal() const { return nxLocal_; }    ///< owned cells in i
  int nyLocal() const { return nyLocal_; }    ///< owned cells in j
  int i0() const { return i0_; }              ///< global i of the first owned column
  int j0() const { return j0_; }              ///< global j of the first owned row
  /// Neighbor rank toward smaller i (MPI_PROC_NULL at the domain edge).
  int rankWest() const { return rankWest_; }
  /// Neighbor rank toward larger i (MPI_PROC_NULL at the domain edge).
  int rankEast() const { return rankEast_; }
  int rankSouth() const { return rankSouth_; }  ///< toward smaller j
  int rankNorth() const { return rankNorth_; }  ///< toward larger j
  ///@}

  /// Block size of dimension-\p N split over \p P parts at part \p p:
  /// N/P + (p < N%P). Exposed for unit testing (plan §8.1 test_grid).
  static int blockSize(int N, int P, int p);
  /// Global start offset of part \p p under the blockSize rule.
  static int blockStart(int N, int P, int p);

  /// \name Active-cell masking and compressed global ids
  ///@{
  /// Set the per-column top-active index and (re)build gid2/gid3.
  /// \param ktop host view sized (nyLocal, nxLocal) — interior cells only —
  ///        with values in [0, nz]; nz marks a fully inactive column.
  /// Collective on comm(). Without a call, all cells are active (ktop = 0).
  void buildGlobalIds(const HostField2<int>& ktop);

  /// Owned active surface cells on this rank / globally, and this rank's
  /// first global row index in the surface system.
  PetscInt activeCount2Local() const { return active2Local_; }
  /// \copydoc activeCount2Local
  PetscInt activeCount2Global() const { return active2Global_; }
  /// First global surface-system row owned by this rank.
  PetscInt offset2() const { return offset2_; }
  /// Owned active subsurface cells on this rank.
  PetscInt activeCount3Local() const { return active3Local_; }
  /// Global active subsurface cells.
  PetscInt activeCount3Global() const { return active3Global_; }
  /// First global subsurface-system row owned by this rank.
  PetscInt offset3() const { return offset3_; }

  /// Surface global ids, host copy, sized (nyLocal+2, nxLocal+2) including
  /// halo columns/rows exchanged from the owning neighbors; -1 = inactive or
  /// outside the domain.
  const HostField2<PetscInt>& gid2Host() const { return gid2Host_; }
  /// Subsurface global ids, host copy, sized (nyLocal+2, nxLocal+2, nz).
  const HostField3<PetscInt>& gid3Host() const { return gid3Host_; }
  /// Device copies of the global-id tables for assembly kernels.
  Field2<PetscInt> gid2() const { return gid2_; }
  /// \copydoc gid2
  Field3<PetscInt> gid3() const { return gid3_; }
  /// The ktop mask, host copy sized (nyLocal, nxLocal).
  const HostField2<int>& ktopHost() const { return ktopHost_; }
  ///@}

 private:
  void chooseDecomposition(const DecompositionConfig& request);
  void exchangeGidHalos();

  MPI_Comm comm_ = MPI_COMM_NULL;
  int rank_ = 0;
  int size_ = 1;

  int nx_ = 0, ny_ = 0, nz_ = 0;
  real_t dx_ = 0, dy_ = 0;
  std::vector<real_t> dz_;
  std::vector<real_t> zCenter_;
  real_t totalDepth_ = 0;

  int px_ = 1, py_ = 1, pi_ = 0, pj_ = 0;
  int nxLocal_ = 0, nyLocal_ = 0, i0_ = 0, j0_ = 0;
  int rankWest_ = MPI_PROC_NULL, rankEast_ = MPI_PROC_NULL;
  int rankSouth_ = MPI_PROC_NULL, rankNorth_ = MPI_PROC_NULL;

  HostField2<int> ktopHost_;
  HostField2<PetscInt> gid2Host_;
  HostField3<PetscInt> gid3Host_;
  Field2<PetscInt> gid2_;
  Field3<PetscInt> gid3_;
  PetscInt active2Local_ = 0, active2Global_ = 0, offset2_ = 0;
  PetscInt active3Local_ = 0, active3Global_ = 0, offset3_ = 0;
};

}  // namespace frehg

#endif  // FREHG_CORE_GRID_HPP
