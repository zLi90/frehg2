/// \file RichardsSolver.cpp
/// \brief RichardsSolver construction, initial state, boundary staging, and
///        step orchestration (legacy solve_groundwater, groundwater.c:57-199).

#include "gw/RichardsSolver.hpp"

#include "core/Logger.hpp"
#include "core/Timer.hpp"
#include "io/GridDataReader.hpp"
#include "gw/VanGenuchten.hpp"

#include <cmath>
#include <string>
#include <vector>

namespace frehg::gw {

namespace {

/// Copy adjacent interior values into domain-edge ghosts of a 3D device
/// field (soil parameters are static; ghosts mirror the interior so face
/// formulas at domain edges collapse to the interior cell's closure — the
/// legacy arrays left these ghosts uninitialized, a §2.1-class hazard).
void mirrorEdgeGhosts(const Grid& grid, const Field3<real_t>& field) {
  const int nyl = grid.nyLocal();
  const int nxl = grid.nxLocal();
  const int nz = static_cast<int>(field.extent(2));
  const bool west = (grid.rankWest() == MPI_PROC_NULL);
  const bool east = (grid.rankEast() == MPI_PROC_NULL);
  const bool south = (grid.rankSouth() == MPI_PROC_NULL);
  const bool north = (grid.rankNorth() == MPI_PROC_NULL);
  Field3<real_t> f = field;
  Kokkos::parallel_for(
      "gw_mirror_edge_ghosts",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 0}, {nyl + 2, nz}),
      KOKKOS_LAMBDA(const int j, const int k) {
        if (west) {
          f(j, 0, k) = f(j, 1, k);
        }
        if (east) {
          f(j, nxl + 1, k) = f(j, nxl, k);
        }
      });
  Kokkos::parallel_for(
      "gw_mirror_edge_ghosts_j",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {0, 0}, {nxl + 2, nz}),
      KOKKOS_LAMBDA(const int i, const int k) {
        if (south) {
          f(0, i, k) = f(1, i, k);
        }
        if (north) {
          f(nyl + 1, i, k) = f(nyl, i, k);
        }
      });
}

}  // namespace

RichardsSolver::RichardsSolver(const Grid& grid, const FrehgConfig& config,
                               const TerrainMetric& mesh, const BoundarySet& boundaries,
                               HaloExchanger& halo)
    : grid_(grid), mesh_(mesh), halo_(halo) {
  const GroundwaterConfig& gw = config.groundwater;
  ss_ = gw.specificStorage;
  useFull3d_ = gw.useFull3d;
  surplusRedistribute_ =
      (gw.reallocationSurplus == GroundwaterConfig::ReallocationSurplus::Redistribute);
  dtMin_ = gw.timestep.dtMin;
  dtMax_ = gw.timestep.dtMax;
  dqGrow_ = gw.timestep.dqGrow;
  dqShrink_ = gw.timestep.dqShrink;
  courantMax_ = gw.timestep.courantMax;
  dtgNext_ = gw.timestep.dtInit;

  const std::size_t ny2 = static_cast<std::size_t>(grid_.nyLocal()) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(grid_.nxLocal()) + 2;
  const std::size_t nz = static_cast<std::size_t>(grid_.nz());
  const auto alloc = [&](const char* name) { return Field3<real_t>(name, ny2, nx2, nz); };
  const auto allocF = [&](const char* name) { return Field3<real_t>(name, ny2, nx2, nz + 1); };
  h_ = alloc("gw_h");
  hn_ = alloc("gw_hn");
  wc_ = alloc("gw_wc");
  wcn_ = alloc("gw_wcn");
  ch_ = alloc("gw_ch");
  ksx_ = alloc("gw_ksx");
  ksy_ = alloc("gw_ksy");
  ksz_ = alloc("gw_ksz");
  vga_ = alloc("gw_vga");
  vgn_ = alloc("gw_vgn");
  wcs_ = alloc("gw_wcs");
  wcr_ = alloc("gw_wcr");
  aev_ = alloc("gw_aev");
  kx_ = alloc("gw_kx");
  ky_ = alloc("gw_ky");
  kzF_ = allocF("gw_kz_face");
  qx_ = alloc("gw_qx");
  qy_ = alloc("gw_qy");
  qzF_ = allocF("gw_qz_face");
  rRho_ = alloc("gw_r_rho");
  rRhoXp_ = alloc("gw_r_rho_xp");
  rRhoYp_ = alloc("gw_r_rho_yp");
  rRhoZp_ = allocF("gw_r_rho_zp");
  rVisc_ = alloc("gw_r_visc");
  rViscXp_ = alloc("gw_r_visc_xp");
  rViscYp_ = alloc("gw_r_visc_yp");
  rViscZp_ = allocF("gw_r_visc_zp");
  vloss_ = alloc("gw_vloss");
  room_ = alloc("gw_room");
  sendUp_ = alloc("gw_send_up");
  sendDown_ = alloc("gw_send_down");
  qOut_ = alloc("gw_q_out");

  // Density/viscosity ratio hooks: present in every formula, held at 1
  // until P4 activates density coupling (plan §10 P2).
  Kokkos::deep_copy(rRho_, 1.0);
  Kokkos::deep_copy(rRhoXp_, 1.0);
  Kokkos::deep_copy(rRhoYp_, 1.0);
  Kokkos::deep_copy(rRhoZp_, 1.0);
  Kokkos::deep_copy(rVisc_, 1.0);
  Kokkos::deep_copy(rViscXp_, 1.0);
  Kokkos::deep_copy(rViscYp_, 1.0);
  Kokkos::deep_copy(rViscZp_, 1.0);
  Kokkos::deep_copy(vloss_, 0.0);

  topCode_ = Field2<int>("gw_top_code", ny2, nx2);
  botCode_ = Field2<int>("gw_bot_code", ny2, nx2);
  cplMode_ = Field2<int>("gw_cpl_mode", ny2, nx2);
  Kokkos::deep_copy(cplMode_, 0);
  topValue_ = Field2<real_t>("gw_top_value", ny2, nx2);
  botValue_ = Field2<real_t>("gw_bot_value", ny2, nx2);
  sideCodeXm_ = Field2<int>("gw_side_code_xm", ny2, nx2);
  sideCodeXp_ = Field2<int>("gw_side_code_xp", ny2, nx2);
  sideCodeYm_ = Field2<int>("gw_side_code_ym", ny2, nx2);
  sideCodeYp_ = Field2<int>("gw_side_code_yp", ny2, nx2);
  sideValueXm_ = Field2<real_t>("gw_side_value_xm", ny2, nx2);
  sideValueXp_ = Field2<real_t>("gw_side_value_xp", ny2, nx2);
  sideValueYm_ = Field2<real_t>("gw_side_value_ym", ny2, nx2);
  sideValueYp_ = Field2<real_t>("gw_side_value_yp", ny2, nx2);

  halo_.add("gw_h", h_);
  halo_.add("gw_wc", wc_);
  halo_.add("gw_ksx", ksx_);
  halo_.add("gw_ksy", ksy_);
  halo_.add("gw_ksz", ksz_);
  halo_.add("gw_vga", vga_);
  halo_.add("gw_vgn", vgn_);
  halo_.add("gw_wcs", wcs_);
  halo_.add("gw_wcr", wcr_);
  halo_.add("gw_aev", aev_);

  stageSoil(config);
  buildBoundaryLists(boundaries, config);
  updateBoundaryValues(config.time.tStart);
  applyInitialConditions(config);

  system_ = std::make_unique<LinearSystem>(
      grid_.comm(), "gw_", grid_.activeCount3Local(), grid_.activeCount3Global(),
      SolverSettings{1.0e-8, 1.0e-14, 1000});
  buildCooPattern();
}

void RichardsSolver::stageSoil(const FrehgConfig& config) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const SoilConfig& soil = config.soil;

  // Resolve the per-cell type index: uniform by name, or the 3D id field
  // ((j*nx + i)*nz + k order; ids index soil.types).
  std::vector<int> typeIndex(
      static_cast<std::size_t>(nyl) * static_cast<std::size_t>(nxl) * static_cast<std::size_t>(nz),
      0);
  if (soil.map.fromFile) {
    const std::vector<real_t> raw = io::readField3D(config.resolvePath(soil.map.file),
                                                    grid_.nx(), grid_.ny(), nz);
    for (int j = 0; j < nyl; ++j) {
      for (int i = 0; i < nxl; ++i) {
        for (int k = 0; k < nz; ++k) {
          const std::size_t g =
              (static_cast<std::size_t>(grid_.j0() + j) * static_cast<std::size_t>(grid_.nx()) +
               static_cast<std::size_t>(grid_.i0() + i)) *
                  static_cast<std::size_t>(nz) +
              static_cast<std::size_t>(k);
          const real_t value = raw[g];
          const int id = static_cast<int>(std::llround(value));
          if (std::fabs(value - static_cast<real_t>(id)) > 1.0e-12 || id < 0 ||
              id >= static_cast<int>(soil.types.size())) {
            log::fatal(log::msg()
                       << "soil.map.file: id " << value << " is not a valid index into "
                       << soil.types.size() << " soil type(s)");
          }
          typeIndex[(static_cast<std::size_t>(j) * static_cast<std::size_t>(nxl) +
                     static_cast<std::size_t>(i)) *
                        static_cast<std::size_t>(nz) +
                    static_cast<std::size_t>(k)] = id;
        }
      }
    }
  } else {
    int id = -1;
    for (std::size_t n = 0; n < soil.types.size(); ++n) {
      if (soil.types[n].name == soil.map.constantName) {
        id = static_cast<int>(n);
        break;
      }
    }
    if (id < 0) {
      log::fatal(log::msg() << "soil.map.constant: type '" << soil.map.constantName
                            << "' is not defined in soil.types");
    }
    for (int& t : typeIndex) {
      t = id;
    }
  }

  const auto fill = [&](const Field3<real_t>& field, auto member) {
    auto host = Kokkos::create_mirror_view(field);
    Kokkos::deep_copy(host, 0.0);
    for (int j = 0; j < nyl; ++j) {
      for (int i = 0; i < nxl; ++i) {
        for (int k = 0; k < nz; ++k) {
          const SoilType& type =
              soil.types[static_cast<std::size_t>(typeIndex[(static_cast<std::size_t>(j) *
                                                                 static_cast<std::size_t>(nxl) +
                                                             static_cast<std::size_t>(i)) *
                                                                static_cast<std::size_t>(nz) +
                                                            static_cast<std::size_t>(k)])];
          host(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1,
               static_cast<std::size_t>(k)) = member(type);
        }
      }
    }
    Kokkos::deep_copy(field, host);
  };
  fill(ksx_, [](const SoilType& t) { return t.ksx; });
  fill(ksy_, [](const SoilType& t) { return t.ksy; });
  fill(ksz_, [](const SoilType& t) { return t.ksz; });
  fill(vga_, [](const SoilType& t) { return t.vgAlpha; });
  fill(vgn_, [](const SoilType& t) { return t.vgN; });
  fill(wcs_, [](const SoilType& t) { return t.thetaS; });
  fill(wcr_, [](const SoilType& t) { return t.thetaR; });
  fill(aev_, [](const SoilType& t) { return t.aev; });

  // Interface halos carry the true neighbor soil; domain-edge ghosts mirror
  // the interior (the soil map is static, so once is enough).
  halo_.exchange({"gw_ksx", "gw_ksy", "gw_ksz", "gw_vga", "gw_vgn", "gw_wcs", "gw_wcr", "gw_aev"});
  for (const Field3<real_t>& f : {ksx_, ksy_, ksz_, vga_, vgn_, wcs_, wcr_, aev_}) {
    mirrorEdgeGhosts(grid_, f);
  }
}

void RichardsSolver::applyInitialConditions(const FrehgConfig& config) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const GroundwaterInitialConfig& ic = config.initialConditions.groundwater;

  // Materialize the configured field: head/moisture are 3D flat lists in
  // (j*nx+i)*nz+k order; a water table is a per-column 2D raster (legacy
  // init_wt_abs reads one elevation per column, initialize.c:1097-1115).
  const bool waterTable = (ic.form == GroundwaterInitialConfig::Form::WaterTable);
  std::vector<real_t> global;
  if (ic.value.fromFile) {
    global = waterTable
                 ? io::readRaster2D(config.resolvePath(ic.value.file), grid_.nx(), grid_.ny())
                 : io::readField3D(config.resolvePath(ic.value.file), grid_.nx(), grid_.ny(), nz);
  }
  const auto icValue = [&](int j, int i, int k) -> real_t {
    if (!ic.value.fromFile) {
      return ic.value.constant;
    }
    const std::size_t column =
        static_cast<std::size_t>(grid_.j0() + j - 1) * static_cast<std::size_t>(grid_.nx()) +
        static_cast<std::size_t>(grid_.i0() + i - 1);
    if (waterTable) {
      return global[column];
    }
    return global[column * static_cast<std::size_t>(nz) + static_cast<std::size_t>(k)];
  };

  auto hHost = Kokkos::create_mirror_view(h_);
  auto wcHost = Kokkos::create_mirror_view(wc_);
  auto bathHost = Kokkos::create_mirror_view(mesh_.bath());
  auto botHost = Kokkos::create_mirror_view(mesh_.bot3d());
  auto dzHost = Kokkos::create_mirror_view(mesh_.dz3d());
  const auto soilHost = [&](const Field3<real_t>& f) {
    auto host = Kokkos::create_mirror_view(f);
    Kokkos::deep_copy(host, f);
    return host;
  };
  auto vgaHost = soilHost(vga_);
  auto vgnHost = soilHost(vgn_);
  auto wcsHost = soilHost(wcs_);
  auto wcrHost = soilHost(wcr_);
  auto aevHost = soilHost(aev_);
  Kokkos::deep_copy(hHost, 0.0);
  Kokkos::deep_copy(wcHost, 0.0);
  Kokkos::deep_copy(bathHost, mesh_.bath());
  Kokkos::deep_copy(botHost, mesh_.bot3d());
  Kokkos::deep_copy(dzHost, mesh_.dz3d());

  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      for (int k = 0; k < nz; ++k) {
        const std::size_t jj = static_cast<std::size_t>(j);
        const std::size_t ii = static_cast<std::size_t>(i);
        const std::size_t kk = static_cast<std::size_t>(k);
        const VgSoil soil{vgaHost(jj, ii, kk), vgnHost(jj, ii, kk), wcsHost(jj, ii, kk),
                          wcrHost(jj, ii, kk), aevHost(jj, ii, kk)};
        const real_t zc = botHost(jj, ii, kk) + 0.5 * dzHost(jj, ii, kk);
        const real_t value = icValue(j, i, k);
        real_t h = 0.0;
        real_t wc = 0.0;
        switch (ic.form) {
          case GroundwaterInitialConfig::Form::Moisture:
            // Legacy ic_subsurface:1032-1059: saturated cells get the
            // hydrostatic head from the bed, unsaturated ones invert the
            // retention curve.
            wc = value;
            h = (wc >= soil.thetaS) ? (bathHost(jj, ii) - zc) : headFromWaterContent(soil, wc);
            break;
          case GroundwaterInitialConfig::Form::Head:
            // Legacy ic_subsurface:1062-1069.
            h = value;
            wc = waterContentFromHead(soil, h);
            break;
          case GroundwaterInitialConfig::Form::WaterTable:
            // Legacy ic_subsurface:1097-1115 (init_wt_abs): pressure head
            // hydrostatic about the table elevation (configured in the
            // un-shifted datum).
            h = (value + mesh_.elevationOffset()) - zc;
            wc = waterContentFromHead(soil, h);
            break;
        }
        hHost(jj, ii, kk) = h;
        wcHost(jj, ii, kk) = wc;
      }
    }
  }
  Kokkos::deep_copy(h_, hHost);
  Kokkos::deep_copy(wc_, wcHost);

  halo_.exchange({"gw_h", "gw_wc"});
  enforceHeadBc();
  enforceMoistureBc();
  Kokkos::deep_copy(hn_, h_);
  Kokkos::deep_copy(wcn_, wc_);
}

void RichardsSolver::buildBoundaryLists(const BoundarySet& boundaries,
                                        const FrehgConfig& config) {
  (void)config;
  // Marker defaults: no-flux everywhere (legacy bctype_GW code 0).
  Kokkos::deep_copy(topCode_, static_cast<int>(GwBcCode::NoFlux));
  Kokkos::deep_copy(botCode_, static_cast<int>(GwBcCode::NoFlux));
  Kokkos::deep_copy(sideCodeXm_, static_cast<int>(GwBcCode::NoFlux));
  Kokkos::deep_copy(sideCodeXp_, static_cast<int>(GwBcCode::NoFlux));
  Kokkos::deep_copy(sideCodeYm_, static_cast<int>(GwBcCode::NoFlux));
  Kokkos::deep_copy(sideCodeYp_, static_cast<int>(GwBcCode::NoFlux));
  Kokkos::deep_copy(topValue_, 0.0);
  Kokkos::deep_copy(botValue_, 0.0);
  Kokkos::deep_copy(sideValueXm_, 0.0);
  Kokkos::deep_copy(sideValueXp_, 0.0);
  Kokkos::deep_copy(sideValueYm_, 0.0);
  Kokkos::deep_copy(sideValueYp_, 0.0);

  for (const BoundaryCondition& bc : boundaries.all()) {
    if (bc.target() != BcTarget::GroundwaterTop && bc.target() != BcTarget::GroundwaterBottom &&
        bc.target() != BcTarget::GroundwaterSide) {
      continue;
    }
    if (bc.kind() == BcKind::ScalarValue) {
      continue;  // consumed by the transport module (plan §10 P4)
    }
    GwBcCode code = GwBcCode::NoFlux;
    if (bc.kind() == BcKind::Head) {
      code = (bc.valueForm() == BcValueConfig::Form::Hydrostatic) ? GwBcCode::HeadHydrostatic
                                                                  : GwBcCode::Head;
    } else if (bc.kind() == BcKind::Flux) {
      code = (bc.valueForm() == BcValueConfig::Form::Gravity) ? GwBcCode::Gravity : GwBcCode::Flux;
    } else {
      log::fatal(log::msg() << "groundwater boundary condition '" << bc.name()
                            << "' must be of kind head or flux");
    }

    GwDeviceBcList list;
    list.bc = &bc;
    list.code = static_cast<int>(code);
    const std::size_t n = bc.cells().size();
    list.j = Kokkos::View<int*, MemSpace>("gw_bc_j", n);
    list.i = Kokkos::View<int*, MemSpace>("gw_bc_i", n);
    list.face = Kokkos::View<int*, MemSpace>("gw_bc_face", n);
    auto hj = Kokkos::create_mirror_view(list.j);
    auto hi = Kokkos::create_mirror_view(list.i);
    auto hf = Kokkos::create_mirror_view(list.face);
    for (std::size_t m = 0; m < n; ++m) {
      hj(m) = bc.cells()[m].j;
      hi(m) = bc.cells()[m].i;
      BcFace face = bc.cells()[m].face;
      if (bc.target() == BcTarget::GroundwaterTop) {
        face = BcFace::Top;
      } else if (bc.target() == BcTarget::GroundwaterBottom) {
        face = BcFace::Bottom;
      }
      hf(m) = static_cast<int>(face);
    }
    Kokkos::deep_copy(list.j, hj);
    Kokkos::deep_copy(list.i, hi);
    Kokkos::deep_copy(list.face, hf);
    bcLists_.push_back(list);

    // Static code markers (values refresh per step).
    Field2<int> topCode = topCode_;
    Field2<int> botCode = botCode_;
    Field2<int> xm = sideCodeXm_, xp = sideCodeXp_, ym = sideCodeYm_, yp = sideCodeYp_;
    auto lj = list.j;
    auto li = list.i;
    auto lf = list.face;
    const int codeInt = list.code;
    Kokkos::parallel_for(
        "gw_mark_bc", Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, n),
        KOKKOS_LAMBDA(const std::size_t m) {
          switch (static_cast<BcFace>(lf(m))) {
            case BcFace::Top:
              topCode(lj(m), li(m)) = codeInt;
              break;
            case BcFace::Bottom:
              botCode(lj(m), li(m)) = codeInt;
              break;
            case BcFace::XMinus:
              xm(lj(m), li(m)) = codeInt;
              break;
            case BcFace::XPlus:
              xp(lj(m), li(m)) = codeInt;
              break;
            case BcFace::YMinus:
              ym(lj(m), li(m)) = codeInt;
              break;
            case BcFace::YPlus:
              yp(lj(m), li(m)) = codeInt;
              break;
          }
        });
  }
}

void RichardsSolver::updateBoundaryValues(real_t t) {
  for (const GwDeviceBcList& list : bcLists_) {
    real_t value = list.bc->value(t);
    if (static_cast<GwBcCode>(list.code) == GwBcCode::HeadHydrostatic) {
      // The reference stage is configured in the un-shifted datum; ghost
      // heads are formed in the offset frame (enforce_head_bc:769-788).
      value += mesh_.elevationOffset();
    }
    Field2<real_t> topValue = topValue_;
    Field2<real_t> botValue = botValue_;
    Field2<real_t> xm = sideValueXm_, xp = sideValueXp_, ym = sideValueYm_, yp = sideValueYp_;
    auto lj = list.j;
    auto li = list.i;
    auto lf = list.face;
    const real_t v = value;
    Kokkos::parallel_for(
        "gw_update_bc_value", Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, list.j.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m) {
          switch (static_cast<BcFace>(lf(m))) {
            case BcFace::Top:
              topValue(lj(m), li(m)) = v;
              break;
            case BcFace::Bottom:
              botValue(lj(m), li(m)) = v;
              break;
            case BcFace::XMinus:
              xm(lj(m), li(m)) = v;
              break;
            case BcFace::XPlus:
              xp(lj(m), li(m)) = v;
              break;
            case BcFace::YMinus:
              ym(lj(m), li(m)) = v;
              break;
            case BcFace::YPlus:
              yp(lj(m), li(m)) = v;
              break;
          }
        });
  }
}

void RichardsSolver::buildCooPattern() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const std::size_t nActive = static_cast<std::size_t>(grid_.activeCount3Local());
  const std::size_t nCoo = 7 * nActive;
  cooRows_ = Kokkos::View<PetscInt*, MemSpace>("gw_coo_rows", nCoo);
  cooCols_ = Kokkos::View<PetscInt*, MemSpace>("gw_coo_cols", nCoo);
  cooValues_ = Kokkos::View<real_t*, MemSpace>("gw_coo_values", nCoo);
  rhsVec_ = Kokkos::View<real_t*, MemSpace>("gw_rhs", nActive);
  solVec_ = Kokkos::View<real_t*, MemSpace>("gw_sol", nActive);

  auto hostRows = Kokkos::create_mirror_view(cooRows_);
  auto hostCols = Kokkos::create_mirror_view(cooCols_);
  const HostField3<PetscInt>& gid = grid_.gid3Host();
  const PetscInt offset = grid_.offset3();
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      for (int k = 0; k < nz; ++k) {
        const std::size_t jj = static_cast<std::size_t>(j);
        const std::size_t ii = static_cast<std::size_t>(i);
        const std::size_t kk = static_cast<std::size_t>(k);
        const PetscInt row = gid(jj, ii, kk);
        if (row < 0) {
          continue;
        }
        const std::size_t base = 7 * static_cast<std::size_t>(row - offset);
        // Value ordering fixed as [diag, xm, xp, ym, yp, zm, zp]; the fill
        // kernel in Predictor.cpp writes the same slots. Legs whose
        // neighbor is inactive or outside the domain get column -1 (ignored
        // by PETSc; their coefficients live in the diagonal/RHS rules).
        for (std::size_t leg = 0; leg < 7; ++leg) {
          hostRows(base + leg) = row;
        }
        hostCols(base + 0) = row;
        hostCols(base + 1) = gid(jj, ii - 1, kk);
        hostCols(base + 2) = gid(jj, ii + 1, kk);
        hostCols(base + 3) = gid(jj - 1, ii, kk);
        hostCols(base + 4) = gid(jj + 1, ii, kk);
        hostCols(base + 5) = (k > 0) ? gid(jj, ii, kk - 1) : -1;
        hostCols(base + 6) = (k + 1 < nz) ? gid(jj, ii, kk + 1) : -1;
      }
    }
  }
  Kokkos::deep_copy(cooRows_, hostRows);
  Kokkos::deep_copy(cooCols_, hostCols);
  system_->setPattern(cooRows_, cooCols_);
}

void RichardsSolver::step(real_t t, real_t dtg) {
  Timer::Scoped timer("groundwater");
  audit_ = GwStepAudit{};
  dtgCurrent_ = dtg;
  updateBoundaryValues(t);

  // Step-start state refresh (legacy groundwater.c:68-82): interior h/wc
  // changed after the last exchange (reallocation, clamp), so halos and
  // domain-edge ghosts are refreshed before h^n / θ^n are frozen. The
  // baroclinic ratios update first from the scalar of the last transport
  // step (legacy baroclinic_face runs once per gw step, groundwater.c:87;
  // hoisted ahead of the ghost refresh so the hydrostatic side ghosts see
  // this step's ratios — see enforceHeadBc).
  halo_.exchange({"gw_h", "gw_wc"});
  updateBaroclinicFaces();
  enforceHeadBc();
  enforceMoistureBc();
  Kokkos::deep_copy(hn_, h_);
  Kokkos::deep_copy(wcn_, wc_);

  // Predictor.
  if (cpl_.active) {
    classifyCoupledTop(dtg);
  }
  computeFaceConductivity();
  assembleSystem(dtg);
  fillAndSolve();
  enforceHeadBc();
  halo_.exchange({"gw_h"});

  // Corrector (legacy groundwater.c:134-157; the coupled top bookkeeping is
  // groundwater_flux:835-868 and runs between the fluxes and the θ update,
  // exactly as legacy orders it).
  computeFaceConductivity();
  computeFluxes(dtg);
  if (cpl_.active) {
    applyCoupledTopBookkeeping(dtg);
  }
  updateWaterContent(dtg);
  halo_.exchange({"gw_wc"});

  // Post-allocation and the final clamp (legacy groundwater.c:159-189).
  reallocateWaterContent();
  finalizeWaterContent();
  enforceMoistureBc();
  halo_.exchange({"gw_wc"});

  accumulateBoundaryFlux(dtg);
  adaptTimeStep(dtg);
}

void RichardsSolver::refreshDerivedState() {
  halo_.exchange({"gw_h", "gw_wc"});
  // With an attached scalar the transport refresh already restored the
  // scalar ghosts, so the ratios (and through them the hydrostatic side
  // ghosts) reproduce the running run's values bitwise.
  updateBaroclinicFaces();
  enforceHeadBc();
  enforceMoistureBc();
  Kokkos::deep_copy(hn_, h_);
  Kokkos::deep_copy(wcn_, wc_);
}

real_t RichardsSolver::ownedVolume() const {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  Field3<real_t> wc = wc_;
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<PetscInt> gid = grid_.gid3();
  const real_t az = mesh_.areaZ();
  real_t volume = 0.0;
  Kokkos::parallel_reduce(
      "gw_owned_volume",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sum) {
        if (gid(j, i, k) >= 0) {
          sum += wc(j, i, k) * az * dz3d(j, i, k);
        }
      },
      volume);
  return volume;
}

const Field3<real_t>& RichardsSolver::fluxXPerArea() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  Field3<real_t> out = qOut_, qx = qx_, ax = mesh_.areaX();
  Kokkos::parallel_for(
      "gw_flux_x_out",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        out(j, i, k) = (ax(j, i, k) > 0.0) ? qx(j, i, k) / ax(j, i, k) : 0.0;
      });
  return qOut_;
}

const Field3<real_t>& RichardsSolver::fluxYPerArea() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  Field3<real_t> out = qOut_, qy = qy_, ay = mesh_.areaY();
  Kokkos::parallel_for(
      "gw_flux_y_out",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        out(j, i, k) = (ay(j, i, k) > 0.0) ? qy(j, i, k) / ay(j, i, k) : 0.0;
      });
  return qOut_;
}

const Field3<real_t>& RichardsSolver::fluxZPerArea() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  Field3<real_t> out = qOut_, qzF = qzF_;
  const real_t az = mesh_.areaZ();
  Kokkos::parallel_for(
      "gw_flux_z_out",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        // The face below each cell, as the legacy qz output reports it.
        out(j, i, k) = qzF(j, i, k + 1) / az;
      });
  return qOut_;
}

void RichardsSolver::accumulateBoundaryFlux(real_t dtg) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  Field3<real_t> qx = qx_, qy = qy_, qzF = qzF_;
  Field3<PetscInt> gid = grid_.gid3();

  // Net inflow through domain-boundary faces [m^3]. Positive q is flow in
  // the -axis direction (x/y) and upward (z), the legacy darcy_flux
  // conventions; the signs below turn each into "into the domain".
  real_t inflow = 0.0;
  Kokkos::parallel_reduce(
      "gw_boundary_inflow",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sum) {
        if (gid(j, i, k) < 0) {
          return;
        }
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        const bool top = (k == 0) || (gid(j, i, k - 1) < 0);
        if (top) {
          sum -= qzF(j, i, k) * dtg;
        }
        if (k == nz - 1) {
          sum += qzF(j, i, k + 1) * dtg;
        }
        if (gi == 0) {
          sum -= qx(j, i - 1, k) * dtg;
        }
        if (gi == nxGlobal - 1) {
          sum += qx(j, i, k) * dtg;
        }
        if (gj == 0) {
          sum -= qy(j - 1, i, k) * dtg;
        }
        if (gj == nyGlobal - 1) {
          sum += qy(j, i, k) * dtg;
        }
      },
      inflow);
  audit_.boundaryIn = inflow;
}

}  // namespace frehg::gw
