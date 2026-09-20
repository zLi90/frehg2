/// \file HaloExchanger.hpp
/// \brief Persistent, coalescing halo exchange for registered fields
///        (plan §5.3).
///
/// Fields are registered once at module initialization; \c exchangeAll()
/// packs every registered field per neighbor into a single message, and
/// \c exchange({names}) restricts the transfer to a subset for mid-step
/// updates. Buffers are allocated once in device memory. The protocol is:
/// post MPI_Irecv (4), run pack kernels, fence, post MPI_Isend (4), wait,
/// unpack.
///
/// GPU-aware MPI is a single runtime toggle: when off, messages stage
/// through persistent host mirrors. This class is the only place the toggle
/// exists (plan §5.3). Domain-edge halos are never written by the exchanger;
/// the boundary-condition system owns them.
///
/// exchangeWithCorners() additionally fills the four halo corner cells
/// (plan §5.3 amendment): the SWE velocity interpolation uy/vx uses a
/// four-point stencil with one diagonal neighbor (A1), and the transport
/// dispersion cross terms average diagonal subsurface scalars (P4). It runs
/// two sequential phases — west/east with interior-row packs, then
/// south/north with full-width rows including the already-updated i-halo
/// columns — so corner values arrive without extra messages. Because the
/// wide rows carry the sender's i-halo columns, the boundary-condition
/// ghost fills must run before the exchange.

#ifndef FREHG_CORE_HALOEXCHANGER_HPP
#define FREHG_CORE_HALOEXCHANGER_HPP

#include "core/Grid.hpp"
#include "core/Types.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace frehg {

/// Coalescing halo exchanger over the four grid neighbors.
class HaloExchanger {
 public:
  /// \param grid decomposition to exchange over (kept by reference; the grid
  ///        must outlive the exchanger).
  /// \param gpuAwareMpi pass device buffers directly to MPI (true) or stage
  ///        through host mirrors (false).
  HaloExchanger(const Grid& grid, bool gpuAwareMpi);

  /// Register a 2D field. \p name must be unique; the field must be sized
  /// (nyLocal+2, nxLocal+2). Fields are captured by shallow View copy.
  void add(const std::string& name, const Field2<real_t>& field);

  /// Register a 3D field sized (nyLocal+2, nxLocal+2, nz).
  void add(const std::string& name, const Field3<real_t>& field);

  /// Exchange every registered field (one coalesced message per neighbor).
  void exchangeAll();

  /// Exchange a subset of registered fields by name; unknown names are
  /// fatal. Used for targeted mid-step updates such as exchange({"eta"}).
  void exchange(const std::vector<std::string>& names);

  /// Exchange a subset of registered 2D fields including the four halo
  /// corner cells (two sequential phases; see the class description).
  /// Requesting a 3D field here is fatal: no Frehg2 stencil needs 3D
  /// corners.
  void exchangeWithCorners(const std::vector<std::string>& names);

  /// \return number of registered fields.
  std::size_t fieldCount() const { return entries_.size(); }

  /// \return true when messages are passed to MPI from device memory.
  bool gpuAware() const { return gpuAware_; }

  // NOTE: everything down to the data members is implementation detail, kept
  // public only because nvcc forbids extended __host__ __device__ lambdas
  // (KOKKOS_LAMBDA) inside private member functions (CUDA C++ programming
  // guide, "Extended Lambda Restrictions"). Treat as private.
 public:
  /// One registered field: the unique name and the (shallow-copied) view it
  /// aliases — exactly one of \c f2 / \c f3 is active per \c is3d.
  struct Entry {
    std::string name;         ///< registration key (unique across add() calls)
    bool is3d = false;        ///< selects which view below is the live one
    Field2<real_t> f2;        ///< the 2D view (when !is3d)
    Field3<real_t> f3;        ///< the 3D view (when is3d)
  };

  // Per-direction packing geometry indices.
  static constexpr int kWest = 0;   ///< pack/unpack direction: west neighbor
  static constexpr int kEast = 1;   ///< pack/unpack direction: east neighbor
  static constexpr int kSouth = 2;  ///< pack/unpack direction: south neighbor
  static constexpr int kNorth = 3;  ///< pack/unpack direction: north neighbor

  /// Scalar count one field contributes to a message in one direction
  /// (rows x planes; \p wideRows includes the i-halo columns of the
  /// corner-filling south/north phase).
  std::size_t planeCount(const Entry& entry, int direction, bool wideRows) const;
  /// Size the per-neighbor device buffers (and host mirrors when staging)
  /// to the registered fields' worst-case coalesced message.
  void ensureCapacity();
  /// Pack one field's boundary cells for \p direction into the send buffer
  /// starting at scalar \p offset.
  void packEntry(const Entry& entry, int direction, std::size_t offset, bool wideRows);
  /// Unpack one field's received halo cells for \p direction from the recv
  /// buffer starting at scalar \p offset.
  void unpackEntry(const Entry& entry, int direction, std::size_t offset, bool wideRows);
  /// Run the Irecv/pack/fence/Isend/wait/unpack protocol for the selected
  /// entries over directions [\p dirBegin, \p dirEnd).
  void exchangeSelected(const std::vector<std::size_t>& selected, int dirBegin, int dirEnd,
                        bool wideRows);
  /// Map registration names to entry indices; an unknown name is fatal.
  std::vector<std::size_t> resolveNames(const std::vector<std::string>& names) const;

 private:
  const Grid& grid_;
  bool gpuAware_ = false;
  std::vector<Entry> entries_;

  using Buffer = Kokkos::View<real_t*, MemSpace>;
  using HostBuffer = Kokkos::View<real_t*, HostSpace>;
  std::array<Buffer, 4> sendBuf_;
  std::array<Buffer, 4> recvBuf_;
  std::array<HostBuffer, 4> sendHost_;
  std::array<HostBuffer, 4> recvHost_;
  std::array<int, 4> neighbor_ = {MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL};
};

}  // namespace frehg

#endif  // FREHG_CORE_HALOEXCHANGER_HPP
