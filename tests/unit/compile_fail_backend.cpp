/// \file compile_fail_backend.cpp
/// \brief p5 negative test (v2 plan §2B.3): this file MUST NOT compile.
///
/// It instantiates the backend-consistency invariant with the forbidden
/// pair — a non-host memory space against a PETSc without Kokkos support —
/// and the static_assert in PetscBackendConsistent must fire. The ctest
/// entry `unit.p5.backend_invariant_fires` builds this target expecting
/// failure (WILL_FAIL); if this file ever compiles, the invariant has
/// stopped protecting device builds and the gate goes red.

#include "core/LinearSystem.hpp"

namespace {

/// Stand-in for a device memory space: any type that is not
/// Kokkos::HostSpace exercises the same branch a CudaSpace would, without
/// needing a CUDA toolchain on the CPU CI lanes.
struct FakeDeviceSpace {};

// Forbidden combination: PetscHasKokkos = false with a device MemSpace.
[[maybe_unused]] constexpr bool mustNotCompile =
    frehg::PetscBackendConsistent<false, FakeDeviceSpace>::value;

}  // namespace

int main() { return 0; }
