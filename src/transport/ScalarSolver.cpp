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

}  // namespace

ScalarSolver::ScalarSolver(const Grid& grid, const FrehgConfig& config,
                           const BoundarySet& boundaries, HaloExchanger& halo,
                           const SurfaceWiring& surface, const SubsurfaceWiring& subsurface,
                           const CouplingWiring& coupling)
    : grid_(grid), halo_(halo), surf_(surface), subs_(subsurface), cpl_(coupling) {
  const TransportConfig& tr = config.transport;
  superbee_ = (tr.scheme.advection == TransportSchemeConfig::Advection::Superbee);
  difuX_ = tr.surfaceDiffusivityX;
  difuY_ = tr.surfaceDiffusivityY;
  dispLon_ = tr.dispersionLongitudinal;
  dispLat_ = tr.dispersionTransverse;
  dispMol_ = tr.dispersionMolecular;
  boundMin_ = tr.boundMin;
  hasBoundMax_ = tr.hasBoundMax;
  // The open-bound sentinel replaces the legacy hard 200 (plan §3.2): the
  // limiter guards compare against it exactly as legacy compared against
  // s_lim_hi, and the final clamp only applies with a configured maximum.
  boundMax_ = hasBoundMax_ ? tr.boundMax : 1.0e30;

  const std::size_t ny2 = static_cast<std::size_t>(grid_.nyLocal()) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(grid_.nxLocal()) + 2;
  const std::size_t nz = static_cast<std::size_t>(grid_.nz());

  if (surf_.active) {
    sSurf_ = Field2<real_t>("s_surf", ny2, nx2);
    smSurf_ = Field2<real_t>("s_sm_surf", ny2, nx2);
    sSurfKp_ = Field2<real_t>("s_surf_kp", ny2, nx2);
    sseepage_ = Field2<real_t>("s_seepage", ny2, nx2);
    vsn_ = Field2<real_t>("s_vsn", ny2, nx2);
    vflux_ = Field2<real_t>("s_vflux", ny2, nx2);
    fuOld_ = Field2<real_t>("s_fu_old", ny2, nx2);
    fvOld_ = Field2<real_t>("s_fv_old", ny2, nx2);
    sMinS_ = Field2<real_t>("s_min_surf", ny2, nx2);
    sMaxS_ = Field2<real_t>("s_max_surf", ny2, nx2);
    sFarXm2_ = Field2<real_t>("s_far_xm_2d", ny2, nx2);
    sFarXp2_ = Field2<real_t>("s_far_xp_2d", ny2, nx2);
    sFarYm2_ = Field2<real_t>("s_far_ym_2d", ny2, nx2);
    sFarYp2_ = Field2<real_t>("s_far_yp_2d", ny2, nx2);
    halo_.add("s_surf", sSurf_);
    halo_.add("s_fu_old", fuOld_);
    halo_.add("s_fv_old", fvOld_);
    halo_.add("s_far_xm_2d", sFarXm2_);
    halo_.add("s_far_xp_2d", sFarXp2_);
    halo_.add("s_far_ym_2d", sFarYm2_);
    halo_.add("s_far_yp_2d", sFarYp2_);
  }
  if (subs_.active) {
    dzzTop_ = Field2<real_t>("s_dzz_top", ny2, nx2);
    sSubs_ = Field3<real_t>("s_subs", ny2, nx2, nz);
    smSubs_ = Field3<real_t>("s_sm_subs", ny2, nx2, nz);
    sMin3_ = Field3<real_t>("s_min_subs", ny2, nx2, nz);
    sMax3_ = Field3<real_t>("s_max_subs", ny2, nx2, nz);
    dxx_ = Field3<real_t>("s_dxx", ny2, nx2, nz);
    dyy_ = Field3<real_t>("s_dyy", ny2, nx2, nz);
    dzz_ = Field3<real_t>("s_dzz", ny2, nx2, nz);
    dxy_ = Field3<real_t>("s_dxy", ny2, nx2, nz);
    dxz_ = Field3<real_t>("s_dxz", ny2, nx2, nz);
    dyz_ = Field3<real_t>("s_dyz", ny2, nx2, nz);
    sFarXm3_ = Field3<real_t>("s_far_xm_3d", ny2, nx2, nz);
    sFarXp3_ = Field3<real_t>("s_far_xp_3d", ny2, nx2, nz);
    sFarYm3_ = Field3<real_t>("s_far_ym_3d", ny2, nx2, nz);
    sFarYp3_ = Field3<real_t>("s_far_yp_3d", ny2, nx2, nz);
    kzLower_ = Field3<real_t>("s_kz_lower", ny2, nx2, nz);
    kTop_ = Field2<int>("s_ktop", ny2, nx2);
    halo_.add("s_subs", sSubs_);
    halo_.add("s_dxx", dxx_);
    halo_.add("s_dyy", dyy_);
    halo_.add("s_dzz", dzz_);
    halo_.add("s_dxy", dxy_);
    halo_.add("s_dxz", dxz_);
    halo_.add("s_dyz", dyz_);
    halo_.add("s_far_xm_3d", sFarXm3_);
    halo_.add("s_far_xp_3d", sFarXp3_);
    halo_.add("s_far_ym_3d", sFarYm3_);
    halo_.add("s_far_yp_3d", sFarYp3_);
    halo_.add("s_kz_lower", kzLower_);
    // The flow module's face conductivities and fluxes feed the transport
    // stencils at interface cells (the legacy exchanged ghosts); registering
    // the module's own views under transport names lets the transport step
    // refresh exactly the halos it reads.
    halo_.add("s_gw_kx", subs_.kx);
    halo_.add("s_gw_ky", subs_.ky);
  }

  buildBoundaryLists(boundaries);
  applyInitialConditions(config);

  if (subs_.active) {
    // Per-column top layer from the (static) active mask, including halo
    // columns — the stencil helpers evaluate at interface cells too.
    const int nyl = grid_.nyLocal();
    const int nxl = grid_.nxLocal();
    const int nzc = grid_.nz();
    Field2<int> kTop = kTop_;
    Field3<PetscInt> gid = grid_.gid3();
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
          kTop(j, i) = top;
        });
  }

  // Initial ghost/derived state: side ghosts at t_start, the exported
  // top-cell scalar, and the end-of-step snapshots (legacy
  // initialize.c:587-588 and 630 seed Vsn and the scalar mass from the
  // initial state; Fu/Fv start at zero, matching the zero-initialized
  // legacy flow rates).
  if (subs_.active) {
    enforceSubsurfaceBc(config.time.tStart);
    halo_.exchangeWithCorners({"s_subs"});
  }
  if (surf_.active) {
    halo_.exchange({"s_surf"});
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
    if (bc.kind() != BcKind::ScalarValue) {
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
  const TransportInitialConfig& ic = config.initialConditions.transport;
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
  if (surf_.active) {
    halo_.exchange({"s_surf", "s_fu_old", "s_fv_old"});
  }
  if (subs_.active) {
    // The side ghosts and the exported top-cell scalar are functions of the
    // restored scalar and the boundary values at the completed step's end
    // time — exactly what the uninterrupted run's BC pass wrote.
    enforceSubsurfaceBc(t);
    halo_.exchangeWithCorners({"s_subs"});
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
  Field3<real_t> s = sSubs_, wc = subs_.wc, dz3d = subs_.dz3d;
  Field3<PetscInt> gid = grid_.gid3();
  real_t mass = 0.0;
  Kokkos::parallel_reduce(
      "transport_subs_mass",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sum) {
        if (gid(j, i, k) >= 0) {
          sum += s(j, i, k) * wc(j, i, k) * az * dz3d(j, i, k);
        }
      },
      mass);
  return mass;
}

}  // namespace frehg::transport
