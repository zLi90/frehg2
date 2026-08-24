/// \file SurfaceSolver.cpp
/// \brief SurfaceSolver construction, initial state, and step orchestration.

#include "swe/SurfaceSolver.hpp"

#include "core/Logger.hpp"
#include "core/Timer.hpp"
#include "io/GridDataReader.hpp"
#include "swe/SweFormulas.hpp"

#include <string>
#include <vector>

namespace frehg::swe {

namespace {

/// Fill an interior field from a global flat vector (j*nx + i ordering).
void assignInteriorFromGlobal(const Grid& grid, const Field2<real_t>& field,
                              const std::vector<real_t>& global) {
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  for (int j = 1; j <= grid.nyLocal(); ++j) {
    for (int i = 1; i <= grid.nxLocal(); ++i) {
      const std::size_t g = static_cast<std::size_t>(grid.j0() + j - 1) *
                                static_cast<std::size_t>(grid.nx()) +
                            static_cast<std::size_t>(grid.i0() + i - 1);
      host(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) = global[g];
    }
  }
  Kokkos::deep_copy(field, host);
}

/// Materialize a FileOrConstant into interior cells.
void assignFileOrConstant(const Grid& grid, const Field2<real_t>& field,
                          const FileOrConstant& source, const FrehgConfig& config) {
  if (source.fromFile) {
    const std::vector<real_t> global =
        io::readRaster2D(config.resolvePath(source.file), grid.nx(), grid.ny());
    assignInteriorFromGlobal(grid, field, global);
    return;
  }
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  const real_t value = source.constant;
  Field2<real_t> f = field;
  Kokkos::parallel_for(
      "swe_fill_constant",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) { f(j, i) = value; });
}

}  // namespace

SurfaceSolver::SurfaceSolver(const Grid& grid, const FrehgConfig& config,
                             const BoundarySet& boundaries, HaloExchanger& halo)
    : grid_(grid), halo_(halo) {
  const SurfaceWaterConfig& sw = config.surfaceWater;
  dt_ = config.time.dt;
  gravity_ = sw.gravity;
  viscX_ = sw.viscosityX;
  viscY_ = sw.viscosityY;
  minDepth_ = sw.minDepth;
  wettingFaceDepth_ = sw.wettingFaceDepth;
  thinLayerDepth_ = sw.friction.thinLayerDepth;
  chezy_ = (sw.friction.law == FrictionConfig::Law::Chezy);
  windCfg_ = sw.wind;

  rainIsSeries_ = sw.rainfall.fromSeries;
  rainConstant_ = sw.rainfall.constant;
  if (rainIsSeries_) {
    rainSeries_ = TimeSeries::fromFile(config.resolvePath(sw.rainfall.file));
  }
  evapIsSeries_ = sw.evaporation.fromSeries;
  evapConstant_ = sw.evaporation.constant;
  if (evapIsSeries_) {
    evapSeries_ = TimeSeries::fromFile(config.resolvePath(sw.evaporation.file));
  }
  if (windCfg_.enabled) {
    if (windCfg_.speed.fromSeries) {
      windSpeedSeries_ = TimeSeries::fromFile(config.resolvePath(windCfg_.speed.file));
    }
    if (windCfg_.direction.fromSeries) {
      windDirectionSeries_ = TimeSeries::fromFile(config.resolvePath(windCfg_.direction.file));
    }
  }

  const auto ny2 = static_cast<std::size_t>(grid_.nyLocal()) + 2;
  const auto nx2 = static_cast<std::size_t>(grid_.nxLocal()) + 2;
  const auto alloc = [&](const char* name) { return Field2<real_t>(name, ny2, nx2); };
  eta_ = alloc("swe_eta");
  etan_ = alloc("swe_etan");
  dept_ = alloc("swe_dept");
  deptx_ = alloc("swe_deptx");
  depty_ = alloc("swe_depty");
  bottom_ = alloc("swe_bottom");
  uu_ = alloc("swe_uu");
  vv_ = alloc("swe_vv");
  uy_ = alloc("swe_uy");
  vx_ = alloc("swe_vx");
  Ex_ = alloc("swe_Ex");
  Ey_ = alloc("swe_Ey");
  Dx_ = alloc("swe_Dx");
  Dy_ = alloc("swe_Dy");
  CDx_ = alloc("swe_CDx");
  CDy_ = alloc("swe_CDy");
  Vs_ = alloc("swe_Vs");
  Vsx_ = alloc("swe_Vsx");
  Vsy_ = alloc("swe_Vsy");
  Asx_ = alloc("swe_Asx");
  Asy_ = alloc("swe_Asy");
  Asz_ = alloc("swe_Asz");
  Aszx_ = alloc("swe_Aszx");
  Aszy_ = alloc("swe_Aszy");
  Fu_ = alloc("swe_Fu");
  Fv_ = alloc("swe_Fv");
  cflx_ = alloc("swe_cflx");
  cfly_ = alloc("swe_cfly");
  cflActive_ = alloc("swe_cfl_active");
  Sxp_ = alloc("swe_Sxp");
  Sxm_ = alloc("swe_Sxm");
  Syp_ = alloc("swe_Syp");
  Sym_ = alloc("swe_Sym");
  Sct_ = alloc("swe_Sct");
  Srhs_ = alloc("swe_Srhs");
  etaBcValue_ = alloc("swe_eta_bc_value");
  isEtaBc_ = alloc("swe_is_eta_bc");
  rainMask_ = alloc("swe_rain_mask");
  frictionCoef_ = alloc("swe_friction_coef");

  halo_.add("swe_eta", eta_);
  halo_.add("swe_dept", dept_);
  halo_.add("swe_uu", uu_);
  halo_.add("swe_vv", vv_);
  halo_.add("swe_Ex", Ex_);
  halo_.add("swe_Ey", Ey_);
  halo_.add("swe_Dx", Dx_);
  halo_.add("swe_Dy", Dy_);
  halo_.add("swe_cfl_active", cflActive_);
  halo_.add("swe_Fu", Fu_);
  halo_.add("swe_Fv", Fv_);
  halo_.add("swe_bottom", bottom_);
  halo_.add("swe_is_eta_bc", isEtaBc_);
  halo_.add("swe_eta_bc_value", etaBcValue_);

  readBathymetry(config);
  buildBoundaryLists(boundaries);
  buildRainMask(config);
  assignFileOrConstant(grid_, frictionCoef_, sw.friction.coefficient, config);

  system_ = std::make_unique<LinearSystem>(grid_.comm(), "fs_", grid_.activeCount2Local(),
                                           grid_.activeCount2Global());
  buildCooPattern();

  applyInitialConditions(config);
}

void SurfaceSolver::readBathymetry(const FrehgConfig& config) {
  assignFileOrConstant(grid_, bottom_, config.domain.bottomElevation, config);
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  Field2<real_t> bottom = bottom_;

  // Legacy lifts all elevations so the bed is non-negative
  // (read_bathymetry, initialize.c:142-147). The shift changes how dry
  // cells wet (the regularized dry row is not shift-invariant), so it is
  // part of the preserved algorithm, not an I/O convention.
  real_t localMin = 0.0;
  Kokkos::parallel_reduce(
      "swe_bottom_min",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& value) {
        if (bottom(j, i) < value) {
          value = bottom(j, i);
        }
      },
      Kokkos::Min<real_t>(localMin));
  real_t globalMin = 0.0;
  MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, grid_.comm());
  offset_ = (globalMin < 0.0) ? -globalMin : 0.0;
  const real_t offset = offset_;
  Kokkos::parallel_for(
      "swe_bottom_offset",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) { bottom(j, i) += offset; });

  // Interface halos come from the neighbors; domain-edge halos copy the
  // adjacent interior cell (legacy boundary_bath, initialize.c:158-176).
  halo_.exchange({"swe_bottom"});
  const bool westEdge = (grid_.rankWest() == MPI_PROC_NULL);
  const bool eastEdge = (grid_.rankEast() == MPI_PROC_NULL);
  const bool southEdge = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool northEdge = (grid_.rankNorth() == MPI_PROC_NULL);
  Kokkos::parallel_for(
      "swe_bottom_edge_ghosts", Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nyl + 1),
      KOKKOS_LAMBDA(const int j) {
        if (westEdge) {
          bottom(j, 0) = bottom(j, 1);
        }
        if (eastEdge) {
          bottom(j, nxl + 1) = bottom(j, nxl);
        }
      });
  Kokkos::parallel_for(
      "swe_bottom_edge_ghosts_j",
      Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nxl + 1),
      KOKKOS_LAMBDA(const int i) {
        if (southEdge) {
          bottom(0, i) = bottom(1, i);
        }
        if (northEdge) {
          bottom(nyl + 1, i) = bottom(nyl, i);
        }
      });
}

void SurfaceSolver::buildRainMask(const FrehgConfig& config) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  hasRainExclusion_ = !config.surfaceWater.rainfallExcludePolygon.empty();
  auto host = Kokkos::create_mirror_view(rainMask_);
  Kokkos::deep_copy(host, 0.0);
  if (hasRainExclusion_) {
    const Polygon poly(config.surfaceWater.rainfallExcludePolygon);
    for (int j = 1; j <= nyl; ++j) {
      for (int i = 1; i <= nxl; ++i) {
        const real_t xc = grid_.xCenter(grid_.i0() + i - 1);
        const real_t yc = grid_.yCenter(grid_.j0() + j - 1);
        host(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) =
            poly.contains(xc, yc) ? 0.0 : 1.0;
      }
    }
  } else {
    for (int j = 1; j <= nyl; ++j) {
      for (int i = 1; i <= nxl; ++i) {
        host(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) = 1.0;
      }
    }
  }
  Kokkos::deep_copy(rainMask_, host);
}

void SurfaceSolver::buildBoundaryLists(const BoundarySet& boundaries) {
  for (const BoundaryCondition& bc : boundaries.all()) {
    if (bc.target() != BcTarget::Surface) {
      continue;
    }
    if (bc.kind() == BcKind::ScalarValue) {
      continue;  // consumed by the transport module (P4)
    }
    DeviceBcList list;
    list.bc = &bc;
    const std::size_t n = bc.cells().size();
    list.j = Kokkos::View<int*, MemSpace>("swe_bc_j", n);
    list.i = Kokkos::View<int*, MemSpace>("swe_bc_i", n);
    list.face = Kokkos::View<int*, MemSpace>("swe_bc_face", n);
    auto hj = Kokkos::create_mirror_view(list.j);
    auto hi = Kokkos::create_mirror_view(list.i);
    auto hf = Kokkos::create_mirror_view(list.face);
    for (std::size_t m = 0; m < n; ++m) {
      hj(m) = bc.cells()[m].j;
      hi(m) = bc.cells()[m].i;
      hf(m) = static_cast<int>(bc.cells()[m].face);
    }
    Kokkos::deep_copy(list.j, hj);
    Kokkos::deep_copy(list.i, hi);
    Kokkos::deep_copy(list.face, hf);

    switch (bc.kind()) {
      case BcKind::Eta:
        etaBcs_.push_back(list);
        break;
      case BcKind::Discharge:
        dischargeBcs_.push_back(list);
        break;
      case BcKind::Velocity:
        velocityBcs_.push_back(list);
        break;
      case BcKind::Outflow: {
        // The transmissive stage extrapolation continues the local bed slope
        // across the boundary face: precompute the bed drop from the inward
        // neighbor to the member cell (clamped at zero — an adverse slope
        // releases nothing).
        list.bedDrop = Kokkos::View<real_t*, MemSpace>("swe_bc_bed_drop", n);
        auto hd = Kokkos::create_mirror_view(list.bedDrop);
        auto bottomHost = Kokkos::create_mirror_view(bottom_);
        Kokkos::deep_copy(bottomHost, bottom_);
        for (std::size_t m = 0; m < n; ++m) {
          const std::size_t j = static_cast<std::size_t>(hj(m));
          const std::size_t i = static_cast<std::size_t>(hi(m));
          real_t inward = bottomHost(j, i);
          switch (static_cast<BcFace>(hf(m))) {
            case BcFace::XPlus:
              inward = bottomHost(j, i - 1);
              break;
            case BcFace::XMinus:
              inward = bottomHost(j, i + 1);
              break;
            case BcFace::YPlus:
              inward = bottomHost(j - 1, i);
              break;
            case BcFace::YMinus:
              inward = bottomHost(j + 1, i);
              break;
            default:
              break;
          }
          const real_t drop = inward - bottomHost(j, i);
          hd(m) = (drop > 0.0) ? drop : 0.0;
        }
        Kokkos::deep_copy(list.bedDrop, hd);
        outflowBcs_.push_back(list);
        break;
      }
      default:
        log::fatal(log::msg() << "surface boundary condition '" << bc.name()
                              << "' has a kind the surface module does not accept");
    }
  }

  // Mark eta-condition cells once; the marker is static (regions are static)
  // and its interface halos are refreshed by one exchange here.
  Field2<real_t> marker = isEtaBc_;
  for (const DeviceBcList& list : etaBcs_) {
    auto lj = list.j;
    auto li = list.i;
    Kokkos::parallel_for(
        "swe_mark_eta_bc",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m) { marker(lj(m), li(m)) = 1.0; });
  }
  halo_.exchange({"swe_is_eta_bc"});
}

void SurfaceSolver::buildCooPattern() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const std::size_t nCells = static_cast<std::size_t>(nyl) * static_cast<std::size_t>(nxl);
  const std::size_t nCoo = 5 * nCells;
  cooRows_ = Kokkos::View<PetscInt*, MemSpace>("swe_coo_rows", nCoo);
  cooCols_ = Kokkos::View<PetscInt*, MemSpace>("swe_coo_cols", nCoo);
  cooValues_ = Kokkos::View<real_t*, MemSpace>("swe_coo_values", nCoo);
  rhsVec_ = Kokkos::View<real_t*, MemSpace>("swe_rhs", static_cast<std::size_t>(nyl) *
                                                           static_cast<std::size_t>(nxl));
  solVec_ = Kokkos::View<real_t*, MemSpace>("swe_sol", rhsVec_.extent(0));

  auto hostRows = Kokkos::create_mirror_view(cooRows_);
  auto hostCols = Kokkos::create_mirror_view(cooCols_);
  const HostField2<PetscInt>& gid = grid_.gid2Host();
  const PetscInt offset = grid_.offset2();
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      const std::size_t jj = static_cast<std::size_t>(j);
      const std::size_t ii = static_cast<std::size_t>(i);
      const PetscInt row = gid(jj, ii);
      const std::size_t base = 5 * static_cast<std::size_t>(row - offset);
      // Value ordering fixed as [diag, xm, xp, ym, yp]; the fill kernel in
      // FreeSurface.cpp writes the same slots.
      for (std::size_t leg = 0; leg < 5; ++leg) {
        hostRows(base + leg) = row;
      }
      hostCols(base + 0) = row;
      hostCols(base + 1) = gid(jj, ii - 1);
      hostCols(base + 2) = gid(jj, ii + 1);
      hostCols(base + 3) = gid(jj - 1, ii);
      hostCols(base + 4) = gid(jj + 1, ii);
    }
  }
  Kokkos::deep_copy(cooRows_, hostRows);
  Kokkos::deep_copy(cooCols_, hostCols);
  system_->setPattern(cooRows_, cooCols_);
}

void SurfaceSolver::applyInitialConditions(const FrehgConfig& config) {
  assignFileOrConstant(grid_, eta_, config.initialConditions.surface.eta, config);
  if (config.initialConditions.surface.hasUu) {
    assignFileOrConstant(grid_, uu_, config.initialConditions.surface.uu, config);
  }
  if (config.initialConditions.surface.hasVv) {
    assignFileOrConstant(grid_, vv_, config.initialConditions.surface.vv, config);
  }
  // The initial surface enters the internal offset frame and eta below the
  // bed is lifted onto it (legacy ic_surface, initialize.c:480-489,544-545).
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t offset = offset_;
  Field2<real_t> eta = eta_;
  Field2<real_t> bottom = bottom_;
  Kokkos::parallel_for(
      "swe_ic_clamp",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        eta(j, i) += offset;
        if (eta(j, i) < bottom(j, i)) {
          eta(j, i) = bottom(j, i);
        }
      });
  halo_.exchange({"swe_uu", "swe_vv"});
  beginStep(config.time.tStart);
  refreshDerivedState();
}

const Field2<real_t>& SurfaceSolver::etaAbsolute() {
  if (etaOut_.extent(0) == 0) {
    etaOut_ = Field2<real_t>("swe_eta_out", eta_.extent(0), eta_.extent(1));
  }
  const real_t offset = offset_;
  Field2<real_t> eta = eta_;
  Field2<real_t> out = etaOut_;
  Kokkos::parallel_for(
      "swe_eta_absolute",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>(
          {0, 0}, {eta_.extent(0), eta_.extent(1)}),
      KOKKOS_LAMBDA(const std::size_t j, const std::size_t i) {
        out(j, i) = eta(j, i) - offset;
      });
  return etaOut_;
}

void SurfaceSolver::refreshDerivedState() {
  enforceSurfBc(false);
  halo_.exchange({"swe_eta"});
  updateDepth();
  updateGeometry();
  updateDragCoef();
  // Flow rates from the restored (post-limiter) velocities: the west/south
  // stage-boundary edge-slot reconstruction below reads them. cflx/cfly
  // keep their restored values — they belong to the completed step.
  {
    const int nyl = grid_.nyLocal();
    const int nxl = grid_.nxLocal();
    Field2<real_t> uu = uu_, vv = vv_, Asx = Asx_, Asy = Asy_, Fu = Fu_, Fv = Fv_;
    Kokkos::parallel_for(
        "swe_refresh_flow_rates",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i) {
          Fu(j, i) = uu(j, i) * Asx(j, i);
          Fv(j, i) = vv(j, i) * Asy(j, i);
        });
  }
  enforceVeloBc(VeloBcApply::Refresh);
  halo_.exchangeWithCorners({"swe_uu", "swe_vv"});
  fillVelocityGhostCorners();
  interpolateVelocity();
  Kokkos::deep_copy(cflActive_, 0.0);
}

void SurfaceSolver::beginStep(real_t t) {
  Kokkos::deep_copy(etan_, eta_);
  audit_ = SurfaceStepAudit{};

  rain_ = rainIsSeries_ ? rainSeries_.value(t) : rainConstant_;
  evap_ = evapIsSeries_ ? evapSeries_.value(t) : evapConstant_;
  if (windCfg_.enabled) {
    windSpeed_ = windCfg_.speed.fromSeries ? windSpeedSeries_.value(t) : windCfg_.speed.constant;
    windDirection_ = windCfg_.direction.fromSeries ? windDirectionSeries_.value(t)
                                                   : windCfg_.direction.constant;
  }

  for (DeviceBcList& list : etaBcs_) {
    // Prescribed stages enter the internal offset frame (legacy
    // get_current_bc adds the offset, solve.c:202-204).
    list.current = list.bc->value(t) + offset_;
    Field2<real_t> value = etaBcValue_;
    auto lj = list.j;
    auto li = list.i;
    const real_t current = list.current;
    Kokkos::parallel_for(
        "swe_eta_bc_value",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m) { value(lj(m), li(m)) = current; });
  }
  for (DeviceBcList& list : dischargeBcs_) {
    list.current = list.bc->value(t);
  }
  for (DeviceBcList& list : velocityBcs_) {
    list.current = list.bc->value(t);
  }
}

void SurfaceSolver::solveFreeSurface() {
  Timer::Scoped timer("swe/free_surface");
  enforceSurfBc(true);
  halo_.exchange({"swe_eta"});
  momentumSource();
  if (etaBcs_.empty()) {
    halo_.exchange({"swe_Ex", "swe_Ey", "swe_Dx", "swe_Dy"});
  } else {
    halo_.exchange({"swe_Ex", "swe_Ey", "swe_Dx", "swe_Dy", "swe_eta_bc_value"});
  }
  assembleRhs();
  assembleCoefficients();
  applyOutflowCorrections();
  fillAndSolve();
  accumulateBoundaryFluxes();
  enforceSurfBc(true);
  cflLimiter();
  evapRain();
}

void SurfaceSolver::updateVelocity() {
  Timer::Scoped timer("swe/velocity");
  halo_.exchange({"swe_eta", "swe_cfl_active"});
  updateDepth();
  updateGeometry();
  updateDragCoef();
  updateVelocityField();
  applyVelocityLimiters();
  enforceVeloBc(VeloBcApply::Step);
  halo_.exchangeWithCorners({"swe_uu", "swe_vv"});
  fillVelocityGhostCorners();
  interpolateVelocity();
}

real_t SurfaceSolver::ownedVolume() const {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t cellArea = grid_.dx() * grid_.dy();
  Field2<real_t> dept = dept_;
  real_t volume = 0.0;
  Kokkos::parallel_reduce(
      "swe_owned_volume",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sum) { sum += dept(j, i) * cellArea; },
      volume);
  return volume;
}

real_t SurfaceSolver::maxCfl() const {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  Field2<real_t> cflx = cflx_;
  Field2<real_t> cfly = cfly_;
  real_t result = 0.0;
  Kokkos::parallel_reduce(
      "swe_max_cfl",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& value) {
        if (cflx(j, i) > value) {
          value = cflx(j, i);
        }
        if (cfly(j, i) > value) {
          value = cfly(j, i);
        }
      },
      Kokkos::Max<real_t>(result));
  return result;
}

}  // namespace frehg::swe
