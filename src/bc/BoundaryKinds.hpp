/// \file BoundaryKinds.hpp
/// \brief Enumerations shared by the boundary-condition system (plan §5.6).
///
/// The target/kind/value-form enumerations live in core/Config.hpp because
/// they are part of the configuration contract; this header re-exports them
/// under the bc component together with the face enumeration used by
/// rasterization, so physics modules can include a single bc header.

#ifndef FREHG_BC_BOUNDARYKINDS_HPP
#define FREHG_BC_BOUNDARYKINDS_HPP

#include "core/Config.hpp"

namespace frehg {

/// Cell faces a boundary condition can act on. Horizontal faces follow the
/// legacy bctype_GW ordering (x+, x-, y+, y-, bottom, top) decoded from
/// groundwater.c:224-288 and initialize.c:1129-1131.
enum class BcFace {
  XPlus,   ///< i + 1/2 face on the domain's east edge
  XMinus,  ///< i - 1/2 face on the west edge
  YPlus,   ///< j + 1/2 face on the north edge
  YMinus,  ///< j - 1/2 face on the south edge
  Bottom,  ///< k = nz - 1 lower face
  Top      ///< k = ktop upper face
};

/// One rasterized member cell of a boundary region on this rank.
/// Indices are local interior indices including the halo offset
/// (1..nyLocal / 1..nxLocal).
struct BcCell {
  int j = 0;                 ///< local row index (halo layout)
  int i = 0;                 ///< local column index (halo layout)
  BcFace face = BcFace::Top; ///< which cell face the condition acts on
};

}  // namespace frehg

#endif  // FREHG_BC_BOUNDARYKINDS_HPP
