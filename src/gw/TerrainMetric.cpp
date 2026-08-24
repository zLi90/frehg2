/// \file TerrainMetric.cpp
/// \brief Subsurface mesh construction (legacy map.c:196-616 provenance).

#include "gw/TerrainMetric.hpp"

#include "core/Logger.hpp"
#include "io/GridDataReader.hpp"

#include <cmath>
#include <vector>

namespace frehg::gw {

namespace {

/// Copy adjacent interior values into domain-edge ghosts of a 2D host field.
template <class HostView2>
void fillEdgeGhosts2(const Grid& grid, HostView2& host) {
  const std::size_t nyl = static_cast<std::size_t>(grid.nyLocal());
  const std::size_t nxl = static_cast<std::size_t>(grid.nxLocal());
  for (std::size_t j = 1; j <= nyl; ++j) {
    if (grid.rankWest() == MPI_PROC_NULL) {
      host(j, 0) = host(j, 1);
    }
    if (grid.rankEast() == MPI_PROC_NULL) {
      host(j, nxl + 1) = host(j, nxl);
    }
  }
  for (std::size_t i = 0; i <= nxl + 1; ++i) {
    if (grid.rankSouth() == MPI_PROC_NULL) {
      host(0, i) = host(1, i);
    }
    if (grid.rankNorth() == MPI_PROC_NULL) {
      host(nyl + 1, i) = host(nyl, i);
    }
  }
}

/// Copy adjacent interior columns into domain-edge ghosts of a 3D host field
/// (legacy build_subsurf_map copies bot3d/dz3d into the boundary ghosts,
/// map.c:509-527).
template <class HostView3>
void fillEdgeGhosts3(const Grid& grid, HostView3& host) {
  const std::size_t nyl = static_cast<std::size_t>(grid.nyLocal());
  const std::size_t nxl = static_cast<std::size_t>(grid.nxLocal());
  const std::size_t nz = host.extent(2);
  for (std::size_t k = 0; k < nz; ++k) {
    for (std::size_t j = 1; j <= nyl; ++j) {
      if (grid.rankWest() == MPI_PROC_NULL) {
        host(j, 0, k) = host(j, 1, k);
      }
      if (grid.rankEast() == MPI_PROC_NULL) {
        host(j, nxl + 1, k) = host(j, nxl, k);
      }
    }
    for (std::size_t i = 0; i <= nxl + 1; ++i) {
      if (grid.rankSouth() == MPI_PROC_NULL) {
        host(0, i, k) = host(1, i, k);
      }
      if (grid.rankNorth() == MPI_PROC_NULL) {
        host(nyl + 1, i, k) = host(nyl, i, k);
      }
    }
  }
}

}  // namespace

TerrainMetric::TerrainMetric(const Grid& grid, const FrehgConfig& config, HaloExchanger& halo)
    : followTerrain_(config.domain.followTerrain),
      uniformLayers_(config.domain.terrainLayers == DomainConfig::TerrainLayers::Uniform),
      az_(grid.dx() * grid.dy()), stretch_(config.domain.dzStretch) {
  const std::size_t ny2 = static_cast<std::size_t>(grid.nyLocal()) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(grid.nxLocal()) + 2;
  const std::size_t nz = static_cast<std::size_t>(grid.nz());

  bath_ = Field2<real_t>("gw_bath", ny2, nx2);
  dz3d_ = Field3<real_t>("gw_dz3d", ny2, nx2, nz);
  bot3d_ = Field3<real_t>("gw_bot3d", ny2, nx2, nz);
  ax_ = Field3<real_t>("gw_area_x", ny2, nx2, nz);
  ay_ = Field3<real_t>("gw_area_y", ny2, nx2, nz);
  sinx_ = Field3<real_t>("gw_sin_x", ny2, nx2, nz);
  cosx_ = Field3<real_t>("gw_cos_x", ny2, nx2, nz);
  siny_ = Field3<real_t>("gw_sin_y", ny2, nx2, nz);
  cosy_ = Field3<real_t>("gw_cos_y", ny2, nx2, nz);
  ktop_ = HostField2<int>("gw_ktop", static_cast<std::size_t>(grid.nyLocal()),
                          static_cast<std::size_t>(grid.nxLocal()));
  zcell_ = HostField3<real_t>("gw_zcell", static_cast<std::size_t>(grid.nyLocal()),
                              static_cast<std::size_t>(grid.nxLocal()), nz);

  readBathymetry(grid, config, halo);

  if (followTerrain_) {
    buildTerrainMesh(grid);
  } else {
    buildRegularMesh(grid);
  }
  buildMetrics(grid, halo);
}

void TerrainMetric::readBathymetry(const Grid& grid, const FrehgConfig& config,
                                   HaloExchanger& halo) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  auto host = Kokkos::create_mirror_view(bath_);
  Kokkos::deep_copy(host, 0.0);

  const FileOrConstant& source = config.domain.bottomElevation;
  if (source.fromFile) {
    const std::vector<real_t> global =
        io::readRaster2D(config.resolvePath(source.file), grid.nx(), grid.ny());
    for (int j = 1; j <= nyl; ++j) {
      for (int i = 1; i <= nxl; ++i) {
        const std::size_t g = static_cast<std::size_t>(grid.j0() + j - 1) *
                                  static_cast<std::size_t>(grid.nx()) +
                              static_cast<std::size_t>(grid.i0() + i - 1);
        host(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) = global[g];
      }
    }
  } else {
    for (int j = 1; j <= nyl; ++j) {
      for (int i = 1; i <= nxl; ++i) {
        host(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) = source.constant;
      }
    }
  }

  // The same non-negative lift the surface module applies
  // (initialize.c:142-147); one datum shared by both modules.
  real_t localMin = 0.0;
  real_t localMax = 0.0;
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      const real_t b = host(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
      localMin = (b < localMin) ? b : localMin;
      localMax = (b > localMax) ? b : localMax;
    }
  }
  real_t globalMin = 0.0;
  MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, grid.comm());
  offset_ = (globalMin < 0.0) ? -globalMin : 0.0;
  localMax += offset_;
  MPI_Allreduce(&localMax, &boxTop_, 1, MPI_DOUBLE, MPI_MAX, grid.comm());
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      host(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) += offset_;
    }
  }
  Kokkos::deep_copy(bath_, host);

  halo.add("gw_bath", bath_);
  halo.exchange({"gw_bath"});
  Kokkos::deep_copy(host, bath_);
  fillEdgeGhosts2(grid, host);
  Kokkos::deep_copy(bath_, host);
}

void TerrainMetric::buildRegularMesh(const Grid& grid) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  const int nz = grid.nz();
  auto bathHost = Kokkos::create_mirror_view(bath_);
  Kokkos::deep_copy(bathHost, bath_);
  auto dzHost = Kokkos::create_mirror_view(dz3d_);
  auto botHost = Kokkos::create_mirror_view(bot3d_);
  Kokkos::deep_copy(dzHost, 0.0);
  Kokkos::deep_copy(botHost, 0.0);

  // Fixed layer interfaces anchored at the box top (legacy bot1d,
  // map.c:229-234).
  std::vector<real_t> bot1d(static_cast<std::size_t>(nz));
  real_t level = boxTop_;
  for (int k = 0; k < nz; ++k) {
    level -= grid.dzK(k);
    bot1d[static_cast<std::size_t>(k)] = level;
  }
  if (grid.rank() == 0 && level > 0.0) {
    // Legacy warns when the box bottom sits above the lowest bed
    // (map.c:261-262); columns whose bed is below the box end up empty.
    log::warn(log::msg() << "subsurface box bottom (" << level - offset_
                         << " m) is above the minimum bed elevation; low columns are truncated");
  }

  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      const std::size_t jj = static_cast<std::size_t>(j);
      const std::size_t ii = static_cast<std::size_t>(i);
      const real_t bath = bathHost(jj, ii);
      // Nominal geometry, then the legacy bed-crossing adjustments
      // (map.c:357-384): the crossing layer keeps a shortened dz when at
      // least a quarter layer remains, otherwise it merges into the layer
      // below.
      std::vector<int> active(static_cast<std::size_t>(nz), 0);
      for (int k = 0; k < nz; ++k) {
        const std::size_t kk = static_cast<std::size_t>(k);
        botHost(jj, ii, kk) = bot1d[kk];
        dzHost(jj, ii, kk) = grid.dzK(k);
      }
      for (int k = 0; k < nz; ++k) {
        const std::size_t kk = static_cast<std::size_t>(k);
        if (botHost(jj, ii, kk) >= bath - 1.0e-5) {
          active[kk] = 0;
          if (botHost(jj, ii, kk) - bath < grid.dzK(k)) {
            botHost(jj, ii, kk) = bath;
          }
        } else {
          active[kk] = 1;
          if (bath - botHost(jj, ii, kk) <= grid.dzK(k)) {
            if (bath - botHost(jj, ii, kk) >= 0.25 * grid.dzK(k)) {
              dzHost(jj, ii, kk) = bath - botHost(jj, ii, kk);
            } else if (k + 1 < nz) {
              active[kk] = 0;
              dzHost(jj, ii, kk + 1) += bath - botHost(jj, ii, kk);
            }
          }
        }
      }
      int top = nz;
      for (int k = 0; k < nz; ++k) {
        if (active[static_cast<std::size_t>(k)] == 1) {
          top = k;
          break;
        }
      }
      ktop_(jj - 1, ii - 1) = top;
    }
  }
  Kokkos::deep_copy(dz3d_, dzHost);
  Kokkos::deep_copy(bot3d_, botHost);
}

void TerrainMetric::buildTerrainMesh(const Grid& grid) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  const int nz = grid.nz();
  auto bathHost = Kokkos::create_mirror_view(bath_);
  Kokkos::deep_copy(bathHost, bath_);
  auto dzHost = Kokkos::create_mirror_view(dz3d_);
  auto botHost = Kokkos::create_mirror_view(bot3d_);
  Kokkos::deep_copy(dzHost, 0.0);
  Kokkos::deep_copy(botHost, 0.0);

  const real_t boxBottom = boxTop_ - grid.totalDepth();
  // Geometric layer weights (legacy incre_sum, map.c:337-340); with
  // dz_stretch = 1 every layer gets (bath - boxBottom)/nz.
  real_t increSum = 0.0;
  real_t weight = 1.0;
  const real_t stretch = stretch_;
  for (int k = 0; k < nz; ++k) {
    increSum += weight;
    weight *= stretch;
  }

  bool shallow = false;
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      const std::size_t jj = static_cast<std::size_t>(j);
      const std::size_t ii = static_cast<std::size_t>(i);
      const real_t bath = bathHost(jj, ii);
      if (bath <= boxBottom) {
        shallow = true;
      }
      // Layer profile: scaled columns span from the local bed to the common
      // box bottom (legacy map.c:316-353); uniform columns carry the
      // configured dz * dz_stretch^k profile below their own bed — the
      // terrain-parallel slab of the b5 reference (SERGHEI GwInit.h:109-117;
      // amendment A12).
      real_t dz = uniformLayers_ ? grid.dzK(0) : (bath - boxBottom) / increSum;
      real_t bot = bath;
      for (int k = 0; k < nz; ++k) {
        const std::size_t kk = static_cast<std::size_t>(k);
        if (uniformLayers_) {
          dz = grid.dzK(k);
        }
        bot -= dz;
        dzHost(jj, ii, kk) = dz;
        botHost(jj, ii, kk) = bot;
        if (!uniformLayers_) {
          dz *= stretch;
        }
      }
      ktop_(jj - 1, ii - 1) = 0;
    }
  }
  if (shallow && !uniformLayers_) {
    log::fatal(
        "domain.follow_terrain requires every bed elevation to sit above the "
        "subsurface box bottom (max(bath) - sum of dz layers); lower the box "
        "or disable follow_terrain");
  }
  Kokkos::deep_copy(dz3d_, dzHost);
  Kokkos::deep_copy(bot3d_, botHost);
}

void TerrainMetric::buildMetrics(const Grid& grid, HaloExchanger& halo) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  const int nz = grid.nz();

  // Interface ghosts carry the true neighbor geometry (legacy copied its own
  // edge cells into interface ghosts, map.c:509-527, which makes face areas
  // rank-dependent under follow_terrain — a §2.1-class defect fixed on port;
  // single-rank results are identical).
  halo.add("gw_dz3d", dz3d_);
  halo.add("gw_bot3d", bot3d_);
  halo.exchange({"gw_dz3d", "gw_bot3d"});
  {
    auto dzHost = Kokkos::create_mirror_view(dz3d_);
    auto botHost = Kokkos::create_mirror_view(bot3d_);
    Kokkos::deep_copy(dzHost, dz3d_);
    Kokkos::deep_copy(botHost, bot3d_);
    fillEdgeGhosts3(grid, dzHost);
    fillEdgeGhosts3(grid, botHost);
    Kokkos::deep_copy(dz3d_, dzHost);
    Kokkos::deep_copy(bot3d_, botHost);
  }

  auto dzHost = Kokkos::create_mirror_view(dz3d_);
  auto botHost = Kokkos::create_mirror_view(bot3d_);
  Kokkos::deep_copy(dzHost, dz3d_);
  Kokkos::deep_copy(botHost, bot3d_);
  auto axHost = Kokkos::create_mirror_view(ax_);
  auto ayHost = Kokkos::create_mirror_view(ay_);
  auto sinxHost = Kokkos::create_mirror_view(sinx_);
  auto cosxHost = Kokkos::create_mirror_view(cosx_);
  auto sinyHost = Kokkos::create_mirror_view(siny_);
  auto cosyHost = Kokkos::create_mirror_view(cosy_);
  Kokkos::deep_copy(sinxHost, 0.0);
  Kokkos::deep_copy(cosxHost, 1.0);
  Kokkos::deep_copy(sinyHost, 0.0);
  Kokkos::deep_copy(cosyHost, 1.0);

  const real_t dx = grid.dx();
  const real_t dy = grid.dy();

  if (followTerrain_) {
    // Per-face slope angles from cell-center elevation differences (legacy
    // map.c:531-545); boundary faces keep sin = 0, cos = 1 (legacy default).
    for (int j = 0; j <= nyl + 1; ++j) {
      for (int i = 0; i <= nxl; ++i) {
        const int gi = grid.i0() + i - 1;
        if (gi < 0 || gi >= grid.nx() - 1) {
          continue;
        }
        const std::size_t jj = static_cast<std::size_t>(j);
        const std::size_t ii = static_cast<std::size_t>(i);
        for (int k = 0; k < nz; ++k) {
          const std::size_t kk = static_cast<std::size_t>(k);
          const real_t zc = botHost(jj, ii, kk) + 0.5 * dzHost(jj, ii, kk);
          const real_t zp = botHost(jj, ii + 1, kk) + 0.5 * dzHost(jj, ii + 1, kk);
          const real_t hdiff = std::fabs(zp - zc);
          const real_t dist = std::sqrt(hdiff * hdiff + dx * dx);
          sinxHost(jj, ii, kk) = hdiff / dist;
          cosxHost(jj, ii, kk) = dx / dist;
        }
      }
    }
    for (int j = 0; j <= nyl; ++j) {
      const int gj = grid.j0() + j - 1;
      if (gj < 0 || gj >= grid.ny() - 1) {
        continue;
      }
      for (int i = 0; i <= nxl + 1; ++i) {
        const std::size_t jj = static_cast<std::size_t>(j);
        const std::size_t ii = static_cast<std::size_t>(i);
        for (int k = 0; k < nz; ++k) {
          const std::size_t kk = static_cast<std::size_t>(k);
          const real_t zc = botHost(jj, ii, kk) + 0.5 * dzHost(jj, ii, kk);
          const real_t zp = botHost(jj + 1, ii, kk) + 0.5 * dzHost(jj + 1, ii, kk);
          const real_t hdiff = std::fabs(zp - zc);
          const real_t dist = std::sqrt(hdiff * hdiff + dy * dy);
          sinyHost(jj, ii, kk) = hdiff / dist;
          cosyHost(jj, ii, kk) = dy / dist;
        }
      }
    }
  }

  // Face areas: terrain-following averages the two flanking thicknesses and
  // tilts by cos (legacy map.c:577-580); the regular mesh keeps the legacy
  // one-sided rule Ax = dy * dz3d(cell) (map.c:407-408). Ghost-slot faces
  // (i = 0 / j = 0) hold the west/south boundary or interface face.
  for (int j = 0; j <= nyl + 1; ++j) {
    for (int i = 0; i <= nxl + 1; ++i) {
      const std::size_t jj = static_cast<std::size_t>(j);
      const std::size_t ii = static_cast<std::size_t>(i);
      for (int k = 0; k < nz; ++k) {
        const std::size_t kk = static_cast<std::size_t>(k);
        if (followTerrain_) {
          const real_t dzE =
              (i <= nxl) ? dzHost(jj, ii + 1, kk) : dzHost(jj, ii, kk);
          const real_t dzN =
              (j <= nyl) ? dzHost(jj + 1, ii, kk) : dzHost(jj, ii, kk);
          axHost(jj, ii, kk) = 0.5 * (dzHost(jj, ii, kk) + dzE) * cosxHost(jj, ii, kk) * dy;
          ayHost(jj, ii, kk) = 0.5 * (dzHost(jj, ii, kk) + dzN) * cosyHost(jj, ii, kk) * dx;
        } else {
          axHost(jj, ii, kk) = dy * dzHost(jj, ii, kk);
          ayHost(jj, ii, kk) = dx * dzHost(jj, ii, kk);
        }
      }
    }
  }

  Kokkos::deep_copy(ax_, axHost);
  Kokkos::deep_copy(ay_, ayHost);
  Kokkos::deep_copy(sinx_, sinxHost);
  Kokkos::deep_copy(cosx_, cosxHost);
  Kokkos::deep_copy(siny_, sinyHost);
  Kokkos::deep_copy(cosy_, cosyHost);

  // zcell output values (legacy map.c:417): cell centers back in the
  // un-shifted datum.
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      for (int k = 0; k < nz; ++k) {
        const std::size_t jj = static_cast<std::size_t>(j);
        const std::size_t ii = static_cast<std::size_t>(i);
        const std::size_t kk = static_cast<std::size_t>(k);
        zcell_(jj - 1, ii - 1, kk) =
            botHost(jj, ii, kk) + 0.5 * dzHost(jj, ii, kk) - offset_;
      }
    }
  }
}

}  // namespace frehg::gw
