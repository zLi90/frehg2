/// \file TerrainMetric.hpp
/// \brief Subsurface mesh geometry: bathymetry, active-column masking, and
///        the terrain-following metrics (plan §10 P2; legacy map.c:196-616).
///
/// The subsurface box is anchored at the global maximum bed elevation (the
/// legacy bot1d rule, map.c:229-234): its top sits at max(bath) and it spans
/// the configured nz layers of thickness dz * dz_stretch^k downward. Two
/// vertical meshes are supported, selected by domain.follow_terrain:
///
///  - **Regular** (legacy map.c:357-393): fixed layer interfaces; cells whose
///    bottom lies at or above the local bed are inactive. The layer the bed
///    crosses keeps a shortened dz when at least a quarter layer remains,
///    otherwise it is merged into the layer below — exactly the legacy
///    partial-cell rules. ktop(j, i) is the first active layer (nz for a
///    fully inactive column).
///  - **Terrain-following** (legacy map.c:316-353): per-face slope angles
///    (sin/cos, map.c:530-575) tilt the lateral fluxes and face areas. Two
///    vertical extent rules (domain.terrain_layers, amendment A12): *scaled*
///    columns carry all nz layers between the local bed and the box bottom,
///    scaled per column (the legacy rule — b6's flat-bottom tank); *uniform*
///    columns carry the configured dz * dz_stretch^k profile below their own
///    bed (the terrain-parallel slab of the b5 reference, SERGHEI
///    GwInit.h:109-117).
///
/// Elevations are lifted into the same non-negative frame the surface module
/// uses (-min(bath), legacy initialize.c:142-147) so P3's coupling shares one
/// datum; zcell output subtracts the offset again (legacy map.c:417).
///
/// The driver builds this before Grid::buildGlobalIds (the ktop mask comes
/// from here) and hands it to the groundwater module, which treats it as
/// read-only geometry.

#ifndef FREHG_GW_TERRAINMETRIC_HPP
#define FREHG_GW_TERRAINMETRIC_HPP

#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Types.hpp"

namespace frehg::gw {

/// Read-only subsurface mesh geometry shared by the groundwater module.
class TerrainMetric {
 public:
  /// Read the bathymetry, build the mask and per-cell geometry, and exchange
  /// interface ghosts. Collective on grid.comm(). Registers the fields
  /// "gw_bath", "gw_dz3d", and "gw_bot3d" with \p halo (exchanged once here;
  /// the geometry is static).
  TerrainMetric(const Grid& grid, const FrehgConfig& config, HaloExchanger& halo);

  /// Per-column first active layer, sized (nyLocal, nxLocal) — the input for
  /// Grid::buildGlobalIds.
  const HostField2<int>& ktop() const { return ktop_; }

  /// The elevation lift applied to all elevations (see the file comment).
  real_t elevationOffset() const { return offset_; }

  /// z of the subsurface box top in the offset frame (= global max bath).
  real_t boxTop() const { return boxTop_; }

  /// True when the terrain-following mesh is active.
  bool followTerrain() const { return followTerrain_; }

  /// \name Device geometry fields (halo layout; ghosts filled)
  ///@{
  const Field2<real_t>& bath() const { return bath_; }    ///< bed elevation, offset frame [m]
  const Field3<real_t>& dz3d() const { return dz3d_; }    ///< cell thickness [m]
  const Field3<real_t>& bot3d() const { return bot3d_; }  ///< cell bottom elevation, offset frame [m]
  /// Face area of the x+ face of cell (j, i, k) [m^2]; slot (j, 0, k) holds
  /// the west boundary/interface face.
  const Field3<real_t>& areaX() const { return ax_; }
  /// \copydoc areaX
  const Field3<real_t>& areaY() const { return ay_; }
  const Field3<real_t>& sinX() const { return sinx_; }  ///< x+ face slope sine (terrain)
  const Field3<real_t>& cosX() const { return cosx_; }  ///< x+ face slope cosine
  const Field3<real_t>& sinY() const { return siny_; }  ///< y+ face slope sine
  const Field3<real_t>& cosY() const { return cosy_; }  ///< y+ face slope cosine
  ///@}

  /// Horizontal (z) face area dx * dy [m^2] (legacy Az; uniform).
  real_t areaZ() const { return az_; }

  /// Cell-center elevations relative to the un-shifted datum, interior cells,
  /// sized (nyLocal, nxLocal, nz) — the /groundwater/zcell/0 output
  /// (legacy map.c:417).
  const HostField3<real_t>& zCellHost() const { return zcell_; }

 private:
  void readBathymetry(const Grid& grid, const FrehgConfig& config, HaloExchanger& halo);
  void buildRegularMesh(const Grid& grid);
  void buildTerrainMesh(const Grid& grid);
  void buildMetrics(const Grid& grid, HaloExchanger& halo);

  bool followTerrain_ = false;
  bool uniformLayers_ = false;  ///< domain.terrain_layers: uniform (A12)
  real_t offset_ = 0.0;
  real_t boxTop_ = 0.0;
  real_t az_ = 0.0;
  real_t stretch_ = 1.0;  ///< dz_stretch (terrain-following layer weights)

  HostField2<int> ktop_;
  HostField3<real_t> zcell_;
  Field2<real_t> bath_;
  Field3<real_t> dz3d_, bot3d_;
  Field3<real_t> ax_, ay_;
  Field3<real_t> sinx_, cosx_, siny_, cosy_;
};

}  // namespace frehg::gw

#endif  // FREHG_GW_TERRAINMETRIC_HPP
