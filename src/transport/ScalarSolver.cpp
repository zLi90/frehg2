/// \file ScalarSolver.cpp
/// \brief Transport-module setup, boundary staging, initial conditions,
///        the per-step orchestration, and the end-of-step snapshots.

#include "transport/ScalarSolver.hpp"

#include "core/Logger.hpp"
#include "core/Timer.hpp"
#include "io/GridDataReader.hpp"

#include <cmath>

namespace frehg::transport {

namespace {

/// Materialize a FileOrConstant into a 2D field's interior cells.
void assign2(const Grid& grid, const Field2<real_t>& field, const FileOrConstant& source,
             const FrehgConfig& config) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  if (source.fromFile) {
    const std::vector<real_t> global =
        io::readRaster2D(config.resolvePath(source.file), grid.nx(), grid.ny());
    auto host = Kokkos::create_mirror_view(field);
    Kokkos::deep_copy(host, field);
    for (int j = 0; j < nyl; ++j) {
      for (int i = 0; i < nxl; ++i) {
        host(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1) =
            global[static_cast<std::size_t>(grid.j0() + j) * static_cast<std::size_t>(grid.nx()) +
                   static_cast<std::size_t>(grid.i0() + i)];
      }
    }
    Kokkos::deep_copy(field, host);
    return;
  }
  const real_t value = source.constant;
  Field2<real_t> f = field;
  Kokkos::parallel_for(
      "transport_fill_constant_2d",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) { f(j, i) = value; });
}

/// Materialize a FileOrConstant into a 3D field's interior cells
/// ((j*nx + i)*nz + k file order, the §7 flattening).
void assign3(const Grid& grid, const Field3<real_t>& field, const FileOrConstant& source,
             const FrehgConfig& config) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  const int nz = grid.nz();
  if (source.fromFile) {
    const std::vector<real_t> global =
        io::readField3D(config.resolvePath(source.file), grid.nx(), grid.ny(), nz);
    auto host = Kokkos::create_mirror_view(field);
    Kokkos::deep_copy(host, field);
    for (int j = 0; j < nyl; ++j) {
      for (int i = 0; i < nxl; ++i) {
        for (int k = 0; k < nz; ++k) {
          const std::size_t g =
              (static_cast<std::size_t>(grid.j0() + j) * static_cast<std::size_t>(grid.nx()) +
               static_cast<std::size_t>(grid.i0() + i)) *
                  static_cast<std::size_t>(nz) +
              static_cast<std::size_t>(k);
          host(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1,
               static_cast<std::size_t>(k)) = global[g];
        }
      }
    }
    Kokkos::deep_copy(field, host);
    return;
  }
  const real_t value = source.constant;
  Field3<real_t> f = field;
  Kokkos::parallel_for(
      "transport_fill_constant_3d",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) { f(j, i, k) = value; });
}

/// Per-column top layer from the (static) active mask, including halo
/// columns — the stencil helpers evaluate at interface cells too. A free
/// function rather than constructor code: nvcc requires an extended
/// __host__ __device__ lambda's enclosing function to have a takeable
/// address, which a constructor never has.
void computeTopLayer(const Grid& grid, const Field2<int>& kTop) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  const int nzc = grid.nz();
  Field2<int> f = kTop;
  Field3<PetscInt> gid = grid.gid3();
  Kokkos::parallel_for(
      "transport_ktop",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 0}, {nyl + 2, nxl + 2}),
      KOKKOS_LAMBDA(const int j, const int i) {
        int top = -1;
        for (int k = 0; k < nzc; ++k) {
          if (gid(j, i, k) >= 0) {
            top = k;
            break;
          }
        }
        f(j, i) = top;
      });
}

}  // namespace

ScalarSpec ScalarSpec::salinity(const FrehgConfig& config) {
  const TransportConfig& tr = config.transport;
  ScalarSpec spec;
  spec.prefix = "s";
  spec.field = BcScalar::Salinity;
  spec.isTemperature = false;
  spec.superbee = (tr.scheme.advection == TransportSchemeConfig::Advection::Superbee);
  spec.difuX = tr.surfaceDiffusivityX;
  spec.difuY = tr.surfaceDiffusivityY;
  spec.dispLon = tr.dispersionLongitudinal;
  spec.dispLat = tr.dispersionTransverse;
  spec.dispMol = tr.dispersionMolecular;
  spec.hasBoundMin = true;
  spec.boundMin = tr.boundMin;
  spec.hasBoundMax = tr.hasBoundMax;
  spec.boundMax = tr.boundMax;
  spec.legacyEvapAllowance = tr.legacyEvapAllowance;
  return spec;
}

ScalarSpec ScalarSpec::temperature(const FrehgConfig& config) {
  const TemperatureConfig& tp = config.temperature;
  ScalarSpec spec;
  spec.prefix = "t";
  spec.field = BcScalar::Temperature;
  spec.isTemperature = true;
  spec.superbee = (tp.scheme.advection == TransportSchemeConfig::Advection::Superbee);
  spec.difuX = tp.surfaceDiffusivityX;
  spec.difuY = tp.surfaceDiffusivityY;
  spec.dispLon = tp.dispersivityLongitudinal;
  spec.dispLat = tp.dispersivityTransverse;
  // The molecular slot carries the effective thermal diffusivity
  // alpha_e = lambda_eff / (rho c)_w (used as-is in the tensor).
  spec.dispMol = tp.thermalConductivity / tp.heatCapacityWater;
  spec.hasBoundMin = tp.hasBoundMin;
  spec.boundMin = tp.boundMin;
  spec.hasBoundMax = tp.hasBoundMax;
  spec.boundMax = tp.boundMax;
  spec.legacyEvapAllowance = false;
  spec.kappaFactor = tp.heatCapacitySolid / tp.heatCapacityWater;
  spec.heatCapacityWater = tp.heatCapacityWater;
  spec.exchange = tp.surfaceExchange.mode;
  spec.equilibriumT = tp.surfaceExchange.equilibriumTemperature;
  spec.equilibriumK = tp.surfaceExchange.equilibriumCoefficient;
  return spec;
}

ScalarSolver::ScalarSolver(const Grid& grid, const FrehgConfig& config,
                           const BoundarySet& boundaries, HaloExchanger& halo,
                           const SurfaceWiring& surface, const SubsurfaceWiring& subsurface,
                           const CouplingWiring& coupling, const ScalarSpec& spec)
    : grid_(grid), halo_(halo), surf_(surface), subs_(subsurface), cpl_(coupling),
      spec_(spec) {
  superbee_ = spec_.superbee;
  difuX_ = spec_.difuX;
  difuY_ = spec_.difuY;
  dispLon_ = spec_.dispLon;
  dispLat_ = spec_.dispLat;
  dispMol_ = spec_.dispMol;
  hasBoundMin_ = spec_.hasBoundMin;
  hasBoundMax_ = spec_.hasBoundMax;
  // The open-bound sentinels replace the legacy hard limits (plan §3.2):
  // the limiter guards compare against them exactly as legacy compared
  // against s_lim_hi/s_lim_lo, and the final clamps only apply with a
  // configured bound (temperature defaults both bounds open).
  boundMin_ = hasBoundMin_ ? spec_.boundMin : -1.0e30;
  boundMax_ = hasBoundMax_ ? spec_.boundMax : 1.0e30;
  legacyEvapAllowance_ = spec_.legacyEvapAllowance;
  if (spec_.exchange == SurfaceExchangeConfig::Mode::Bulk) {
    met_ = atm::MetForcing(config.atmosphere, config);
  }

  const std::size_t ny2 = static_cast<std::size_t>(grid_.nyLocal()) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(grid_.nxLocal()) + 2;
  const std::size_t nz = static_cast<std::size_t>(grid_.nz());

  // Field and halo names carry the spec prefix: "s_*" is the v1 salinity
  // layout byte-for-byte (checkpoint compatibility); the temperature
  // instance registers a parallel "t_*" set.
  const std::string p = spec_.prefix + "_";
  if (surf_.active) {
    sSurf_ = Field2<real_t>(p + "surf", ny2, nx2);
    smSurf_ = Field2<real_t>(p + "sm_surf", ny2, nx2);
    sSurfKp_ = Field2<real_t>(p + "surf_kp", ny2, nx2);
    sseepage_ = Field2<real_t>(p + "seepage", ny2, nx2);
    vsn_ = Field2<real_t>(p + "vsn", ny2, nx2);
    vflux_ = Field2<real_t>(p + "vflux", ny2, nx2);
    fuOld_ = Field2<real_t>(p + "fu_old", ny2, nx2);
    fvOld_ = Field2<real_t>(p + "fv_old", ny2, nx2);
    sMinS_ = Field2<real_t>(p + "min_surf", ny2, nx2);
    sMaxS_ = Field2<real_t>(p + "max_surf", ny2, nx2);
    sFarXm2_ = Field2<real_t>(p + "far_xm_2d", ny2, nx2);
    sFarXp2_ = Field2<real_t>(p + "far_xp_2d", ny2, nx2);
    sFarYm2_ = Field2<real_t>(p + "far_ym_2d", ny2, nx2);
    sFarYp2_ = Field2<real_t>(p + "far_yp_2d", ny2, nx2);
    halo_.add(p + "surf", sSurf_);
    halo_.add(p + "fu_old", fuOld_);
    halo_.add(p + "fv_old", fvOld_);
    halo_.add(p + "far_xm_2d", sFarXm2_);
    halo_.add(p + "far_xp_2d", sFarXp2_);
    halo_.add(p + "far_ym_2d", sFarYm2_);
    halo_.add(p + "far_yp_2d", sFarYp2_);
  }
  if (subs_.active) {
    dzzTop_ = Field2<real_t>(p + "dzz_top", ny2, nx2);
    sSubs_ = Field3<real_t>(p + "subs", ny2, nx2, nz);
    smSubs_ = Field3<real_t>(p + "sm_subs", ny2, nx2, nz);
    sMin3_ = Field3<real_t>(p + "min_subs", ny2, nx2, nz);
    sMax3_ = Field3<real_t>(p + "max_subs", ny2, nx2, nz);
    dxx_ = Field3<real_t>(p + "dxx", ny2, nx2, nz);
    dyy_ = Field3<real_t>(p + "dyy", ny2, nx2, nz);
    dzz_ = Field3<real_t>(p + "dzz", ny2, nx2, nz);
    dxy_ = Field3<real_t>(p + "dxy", ny2, nx2, nz);
    dxz_ = Field3<real_t>(p + "dxz", ny2, nx2, nz);
    dyz_ = Field3<real_t>(p + "dyz", ny2, nx2, nz);
    sFarXm3_ = Field3<real_t>(p + "far_xm_3d", ny2, nx2, nz);
    sFarXp3_ = Field3<real_t>(p + "far_xp_3d", ny2, nx2, nz);
    sFarYm3_ = Field3<real_t>(p + "far_ym_3d", ny2, nx2, nz);
    sFarYp3_ = Field3<real_t>(p + "far_yp_3d", ny2, nx2, nz);
    kzLower_ = Field3<real_t>(p + "kz_lower", ny2, nx2, nz);
    kTop_ = Field2<int>(p + "ktop", ny2, nx2);
    cauchyTop_ = Field2<int>(p + "cauchy_top", ny2, nx2);
    halo_.add(p + "subs", sSubs_);
    halo_.add(p + "dxx", dxx_);
    halo_.add(p + "dyy", dyy_);
    halo_.add(p + "dzz", dzz_);
    halo_.add(p + "dxy", dxy_);
    halo_.add(p + "dxz", dxz_);
    halo_.add(p + "dyz", dyz_);
    halo_.add(p + "far_xm_3d", sFarXm3_);
    halo_.add(p + "far_xp_3d", sFarXp3_);
    halo_.add(p + "far_ym_3d", sFarYm3_);
    halo_.add(p + "far_yp_3d", sFarYp3_);
    halo_.add(p + "kz_lower", kzLower_);
    // The flow module's face conductivities and fluxes feed the transport
    // stencils at interface cells (the legacy exchanged ghosts); registering
    // the module's own views under transport names lets the transport step
    // refresh exactly the halos it reads. (With both scalars active the
    // kx/ky views are registered twice under both prefixes and exchanged
    // once per instance step — redundant but tiny, and it keeps the
    // instances independent.)
    halo_.add(p + "gw_kx", subs_.kx);
    halo_.add(p + "gw_ky", subs_.ky);
  }

  buildBoundaryLists(boundaries);
  applyInitialConditions(config);

  if (subs_.active) {
    computeTopLayer(grid_, kTop_);
  }

  // Initial ghost/derived state: cell pins and side ghosts at t_start, the
  // exported top-cell scalar, and the end-of-step snapshots (legacy
  // initialize.c:587-588 and 630 seed Vsn and the scalar mass from the
  // initial state; Fu/Fv start at zero, matching the zero-initialized
  // legacy flow rates).
  if (subs_.active) {
    applyCellPins(config.time.tStart);
    audit_ = TransportAudit{};  // pin deltas at t_start are initial state
    enforceSubsurfaceBc(config.time.tStart);
    halo_.exchangeWithCorners({spec_.prefix + "_subs"});
  }
  if (surf_.active) {
    halo_.exchange({spec_.prefix + "_surf"});
  }
  snapshotEndOfStep();
}

void ScalarSolver::buildBoundaryLists(const BoundarySet& boundaries) {
  const auto stage = [](const std::vector<BcCell>& cells) {
    ScalarBcList list;
    const std::size_t n = cells.size();
    list.j = Kokkos::View<int*, MemSpace>("s_bc_j", n);
    list.i = Kokkos::View<int*, MemSpace>("s_bc_i", n);
    list.face = Kokkos::View<int*, MemSpace>("s_bc_face", n);
    auto hj = Kokkos::create_mirror_view(list.j);
    auto hi = Kokkos::create_mirror_view(list.i);
    auto hf = Kokkos::create_mirror_view(list.face);
    for (std::size_t m = 0; m < n; ++m) {
      hj(m) = cells[m].j;
      hi(m) = cells[m].i;
      hf(m) = static_cast<int>(cells[m].face);
    }
    Kokkos::deep_copy(list.j, hj);
    Kokkos::deep_copy(list.i, hi);
    Kokkos::deep_copy(list.face, hf);
    return list;
  };

  for (const BoundaryCondition& bc : boundaries.all()) {
    if (bc.kind() == BcKind::ScalarCauchy && spec_.isTemperature) {
      // The zero-total-scalar-flux relation is a SALT condition (water
      // leaves, salt stays); heat leaves with the water — the temperature
      // instance's donor-value top face is that physics, so the condition
      // does not bind here (schema restricts it to salinity anyway).
      continue;
    }
    if (bc.kind() == BcKind::ScalarCauchy) {
      // The v2 §3.2 zero-total-scalar-flux top condition (Geng & Boufadel
      // Eq. (7)): mark the member columns; the subsurface step reads the
      // mask at the top face (no scalar crosses) and in the limiter update
      // (the exact concentration/dilution allowance). Schema restricts the
      // kind to target groundwater_top in uncoupled groundwater runs.
      if (!subs_.active) {
        log::fatal(log::msg() << "scalar condition '" << bc.name()
                              << "': scalar_cauchy requires the groundwater module");
      }
      auto host = Kokkos::create_mirror_view(cauchyTop_);
      Kokkos::deep_copy(host, cauchyTop_);
      for (const BcCell& cell : bc.cells()) {
        host(static_cast<std::size_t>(cell.j), static_cast<std::size_t>(cell.i)) = 1;
      }
      Kokkos::deep_copy(cauchyTop_, host);
      continue;
    }
    if (bc.kind() != BcKind::ScalarValue) {
      continue;
    }
    if (bc.scalarField() != spec_.field) {
      // The other registered scalar's condition (v2 Q5 selector).
      continue;
    }
    if (bc.target() == BcTarget::GroundwaterTop ||
        bc.target() == BcTarget::GroundwaterBottom) {
      // Cell-pinning Dirichlet rows (v2 Q5, V2-A17): the member columns'
      // top/bottom cells are re-imposed after each update — the b6 tide
      // rule applied vertically. Temperature-only in v2.0 (schema).
      if (!subs_.active) {
        log::fatal(log::msg() << "scalar condition '" << bc.name()
                              << "': groundwater targets require the groundwater module");
      }
      ScalarBcList list = stage(bc.cells());
      list.bc = &bc;
      if (bc.target() == BcTarget::GroundwaterTop) {
        topPin_.push_back(list);
      } else {
        botPin_.push_back(list);
      }
      continue;
    }
    if (bc.target() == BcTarget::Surface) {
      if (!surf_.active) {
        log::fatal(log::msg() << "scalar condition '" << bc.name()
                              << "': surface scalar_value requires the surface module");
      }
      // Partition the member cells by the flow condition they ride on:
      // discharge-condition cells take the legacy inflow mass source
      // (scalar.c:180-195); every other member is a wet-cell Dirichlet —
      // the legacy tide rule (scalar.c:275-282) generalized to any region
      // (the b6 tidal variant prescribes the flooding seawater over the
      // whole tank; its golden ran the legacy wet-cell salinity override —
      // the "Kuan 2019" block, scalar.c:258-262 — see the b6 README).
      std::vector<BcCell> etaCells;
      std::vector<BcCell> perDischarge;
      const BoundaryCondition* discharge = nullptr;
      for (const BcCell& cell : bc.cells()) {
        const BoundaryCondition* owner = nullptr;
        for (const BoundaryCondition& other : boundaries.all()) {
          if (other.kind() != BcKind::Discharge || other.target() != BcTarget::Surface) {
            continue;
          }
          for (const BcCell& oc : other.cells()) {
            if (oc.j == cell.j && oc.i == cell.i) {
              owner = &other;
              break;
            }
          }
          if (owner != nullptr) {
            break;
          }
        }
        if (owner == nullptr) {
          etaCells.push_back(cell);
          continue;
        }
        if (discharge != nullptr && discharge != owner) {
          log::fatal(log::msg() << "scalar condition '" << bc.name()
                                << "': members span multiple discharge conditions; split the "
                                << "scalar condition per inflow");
        }
        discharge = owner;
        perDischarge.push_back(cell);
      }
      if (!etaCells.empty() || bc.cells().empty()) {
        ScalarBcList list = stage(etaCells);
        list.bc = &bc;
        surfaceDirichlet_.push_back(list);
      }
      if (!perDischarge.empty()) {
        ScalarBcList list = stage(perDischarge);
        list.bc = &bc;
        list.discharge = discharge;
        list.dischargeCells = discharge->globalCellCount();
        surfaceInflow_.push_back(list);
      }
    } else if (bc.target() == BcTarget::GroundwaterSide) {
      if (!subs_.active) {
        log::fatal(log::msg() << "scalar condition '" << bc.name()
                              << "': groundwater_side scalar_value requires the groundwater "
                              << "module");
      }
      ScalarBcList list = stage(bc.cells());
      list.bc = &bc;
      sideGhost_.push_back(list);
    } else {
      log::fatal(log::msg() << "scalar condition '" << bc.name()
                            << "': scalar_value supports targets surface and groundwater_side");
    }
  }
}

void ScalarSolver::applyInitialConditions(const FrehgConfig& config) {
  const TransportInitialConfig& ic = spec_.isTemperature
                                         ? config.initialConditions.temperature
                                         : config.initialConditions.transport;
  if (surf_.active) {
    assign2(grid_, sSurf_, ic.surface, config);
  }
  if (subs_.active) {
    assign3(grid_, sSubs_, ic.groundwater, config);
    // Inactive cells hold no scalar (legacy scalar_groundwater:483).
    const int nyl = grid_.nyLocal();
    const int nxl = grid_.nxLocal();
    const int nz = grid_.nz();
    Field3<real_t> s = sSubs_;
    Field3<PetscInt> gid = grid_.gid3();
    Kokkos::parallel_for(
        "transport_ic_mask",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
            {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
        KOKKOS_LAMBDA(const int j, const int i, const int k) {
          if (gid(j, i, k) < 0) {
            s(j, i, k) = 0.0;
          }
        });
  }
}

void ScalarSolver::applyCellPins(real_t t) {
  if (!subs_.active || (topPin_.empty() && botPin_.empty())) {
    return;
  }
  const int nz = grid_.nz();
  Field3<real_t> s = sSubs_;
  Field3<real_t> wc = subs_.wc, wcs = subs_.wcs, dz3d = subs_.dz3d;
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> kTop = kTop_;
  const real_t az = subs_.az;
  const real_t kappaFactor = spec_.kappaFactor;
  real_t pinned = 0.0;
  const auto applyList = [&](const ScalarBcList& list, const bool top) {
    const real_t value = list.bc->value(t);
    auto lj = list.j;
    auto li = list.i;
    real_t delta = 0.0;
    Kokkos::parallel_reduce(
        "transport_cell_pin", Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, lj.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m, real_t& sum) {
          const int j = lj(m);
          const int i = li(m);
          int k = -1;
          if (top) {
            k = kTop(j, i);
          } else {
            for (int kk = nz - 1; kk >= 0; --kk) {
              if (gid(j, i, kk) >= 0) {
                k = kk;
                break;
              }
            }
          }
          if (k < 0) {
            return;
          }
          const real_t basis =
              (wc(j, i, k) + kappaFactor * (1.0 - wcs(j, i, k))) * az * dz3d(j, i, k);
          sum += (value - s(j, i, k)) * basis;
          s(j, i, k) = value;
        },
        delta);
    pinned += delta;
  };
  for (const ScalarBcList& list : topPin_) {
    applyList(list, true);
  }
  for (const ScalarBcList& list : botPin_) {
    applyList(list, false);
  }
  // The pinned-cell reset is boundary-driven mass (the vertical analogue of
  // the surface tide reset, which books into surfSource).
  audit_.subsBoundary += pinned;
}

void ScalarSolver::snapshotEndOfStep() {
  if (!surf_.active) {
    return;
  }
  // Vsn and the volume_by_flux flow rates of the *next* step are the final
  // depth and flow-rate state of this one (see the file comment in
  // ScalarSolver.hpp). Full-extent copies: the halo and physical-edge slots
  // carry exactly what the flow module left there (exchanged interface
  // values; zero-initialized physical edges — the legacy ghost flow rates).
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t cellArea = grid_.dx() * grid_.dy();
  Field2<real_t> vsn = vsn_, dept = surf_.dept;
  Kokkos::parallel_for(
      "transport_snapshot_vsn",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 0}, {nyl + 2, nxl + 2}),
      KOKKOS_LAMBDA(const int j, const int i) { vsn(j, i) = dept(j, i) * cellArea; });
  Kokkos::deep_copy(fuOld_, surf_.fu);
  Kokkos::deep_copy(fvOld_, surf_.fv);
}

void ScalarSolver::step(real_t t, real_t dt, real_t dtgLast, real_t rain, real_t evap) {
  Timer::Scoped timer("transport");
  audit_ = TransportAudit{};
  if (surf_.active) {
    stepSurface(t, dt, rain, evap);
  }
  if (subs_.active) {
    stepSubsurface(t, dt, dtgLast);
  }
  snapshotEndOfStep();
}

void ScalarSolver::refreshDerivedState(real_t t) {
  const std::string p = spec_.prefix + "_";
  if (surf_.active) {
    halo_.exchange({p + "surf", p + "fu_old", p + "fv_old"});
  }
  if (subs_.active) {
    // The cell pins, side ghosts, and the exported top-cell scalar are
    // functions of the restored scalar and the boundary values at the
    // completed step's end time — exactly what the uninterrupted run's
    // BC pass wrote (the restored pinned cells already hold these values;
    // re-imposing is idempotent).
    applyCellPins(t);
    audit_ = TransportAudit{};
    enforceSubsurfaceBc(t);
    halo_.exchangeWithCorners({p + "subs"});
  }
  if (surf_.active) {
    // vsn reconstructs from the restored depth; fuOld/fvOld were restored
    // (their end-of-step values are pre-correction flow rates a restart
    // cannot recompute from the checkpointed velocities).
    const int nyl = grid_.nyLocal();
    const int nxl = grid_.nxLocal();
    const real_t cellArea = grid_.dx() * grid_.dy();
    Field2<real_t> vsn = vsn_, dept = surf_.dept;
    Kokkos::parallel_for(
        "transport_refresh_vsn",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {0, 0}, {nyl + 2, nxl + 2}),
        KOKKOS_LAMBDA(const int j, const int i) { vsn(j, i) = dept(j, i) * cellArea; });
  }
}

real_t ScalarSolver::ownedSurfaceMass() const {
  if (!surf_.active) {
    return 0.0;
  }
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t cellArea = grid_.dx() * grid_.dy();
  Field2<real_t> s = sSurf_, dept = surf_.dept;
  real_t mass = 0.0;
  Kokkos::parallel_reduce(
      "transport_surf_mass",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sum) {
        sum += s(j, i) * dept(j, i) * cellArea;
      },
      mass);
  return mass;
}

real_t ScalarSolver::ownedSubsurfaceMass() const {
  if (!subs_.active) {
    return 0.0;
  }
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t az = subs_.az;
  const real_t kappaFactor = spec_.kappaFactor;
  Field3<real_t> s = sSubs_, wc = subs_.wc, wcs = subs_.wcs, dz3d = subs_.dz3d;
  Field3<PetscInt> gid = grid_.gid3();
  real_t mass = 0.0;
  Kokkos::parallel_reduce(
      "transport_subs_mass",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sum) {
        if (gid(j, i, k) >= 0) {
          // Salinity: theta V (kappa = 0, the legacy saltmass_subs basis).
          // Temperature: (theta + kappa) V — the heat-content basis of the
          // thermal retardation (V2-A17).
          sum += s(j, i, k) * (wc(j, i, k) + kappaFactor * (1.0 - wcs(j, i, k))) * az *
                 dz3d(j, i, k);
        }
      },
      mass);
  return mass;
}

}  // namespace frehg::transport
