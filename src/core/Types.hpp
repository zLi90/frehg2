/// \file Types.hpp
/// \brief Fundamental type aliases and the fatal-error exception for Frehg2.
///
/// Every Frehg2 component uses these aliases instead of naming Kokkos types
/// directly, so the floating-point type and the execution/memory spaces are
/// selected in exactly one place (upgrade plan §4). Managed/shared memory
/// spaces are prohibited by design: host access to device data goes through
/// \c Kokkos::create_mirror_view and \c Kokkos::deep_copy at I/O boundaries.

#ifndef FREHG_CORE_TYPES_HPP
#define FREHG_CORE_TYPES_HPP

#include <Kokkos_Core.hpp>

#include <stdexcept>
#include <string>

namespace frehg {

/// Floating-point type used for all physical fields. Double-only for the
/// v1.0 release; the alias is kept so precision experiments stay a one-line
/// change (plan §3.3).
using real_t = double;

/// Default device execution space (OpenMP/Serial on CPU builds, CUDA/HIP on
/// GPU builds). Physics code never names a backend explicitly.
using ExecSpace = Kokkos::DefaultExecutionSpace;

/// Memory space in which all simulation fields live.
using MemSpace = ExecSpace::memory_space;

/// Host mirror memory space used at I/O boundaries.
using HostSpace = Kokkos::HostSpace;

/// Two-dimensional cell field, indexed (j, i) with a one-cell halo in each
/// horizontal direction (plan §5.2).
template <class T>
using Field2 = Kokkos::View<T**, Kokkos::LayoutRight, MemSpace>;

/// Three-dimensional cell field, indexed (j, i, k): halo width 1 in j and i,
/// no halo in k, unit stride in k (plan §5.2).
template <class T>
using Field3 = Kokkos::View<T***, Kokkos::LayoutRight, MemSpace>;

/// Host-resident 2D counterpart, used for file I/O staging and tests.
template <class T>
using HostField2 = Kokkos::View<T**, Kokkos::LayoutRight, HostSpace>;

/// Host-resident 3D counterpart, used for file I/O staging and tests.
template <class T>
using HostField3 = Kokkos::View<T***, Kokkos::LayoutRight, HostSpace>;

/// Exception thrown by \c log::fatal after the error has been logged.
///
/// Library code never terminates the process directly; it reports through
/// the logger, which throws this type so that \c main (and unit tests)
/// control process termination.
class FatalError : public std::runtime_error {
 public:
  /// \param message the already-logged, human-readable error description.
  explicit FatalError(const std::string& message) : std::runtime_error(message) {}
};

}  // namespace frehg

#endif  // FREHG_CORE_TYPES_HPP
