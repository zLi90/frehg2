/// \file Evapotranspiration.cpp
/// \brief The Stage-3 T5 ET soil channel (plan §9 W5.4): per-substep
///        potential x limiter x dry-fraction flux, supply-capped, with
///        the time-basis complementarity clocks.
///
/// The channel fixes the v2.0 bulk-soil-evaporation defects on the ET
/// path (the legacy groundwater.evaporation path is byte-identical when
/// the block is absent): D1 (the predictor/corrector guard mismatch —
/// qVap is read guard-free on both sides), D2 (the audit books the
/// channel's own flux, not the whole top face), D4 (the supply cap), and
/// the window-lag double count (the dry fraction reads the SURFACE
/// channel's wet clock at substep resolution, so no instant is covered
/// by both channels).

#include "core/Logger.hpp"
#include "atm/PenmanMonteith.hpp"
#include "bc/Polygon.hpp"
#include "gw/RichardsSolver.hpp"
#include "io/GridDataReader.hpp"

namespace frehg::gw {

void RichardsSolver::setupEvapotranspiration(const FrehgConfig& config) {
  const EvapotranspirationConfig& et = config.evapotranspiration;
  if (!et.enabled || !et.soil) {
    return;
  }
  etActive_ = true;
  etCfg_ = et;
  etAeroLaw_ = config.atmosphere.aeroResistance;
  etAero_ = {config.atmosphere.aeroRefHeight, config.atmosphere.aeroZ0m,
             config.atmosphere.aeroZ0h, config.atmosphere.aeroDisplacement};
  if (et.potential == EvapotranspirationConfig::Potential::Bulk ||
      et.potential == EvapotranspirationConfig::Potential::PenmanMonteith) {
    met_ = atm::MetForcing(config.atmosphere, config);
  } else if (et.prescribed.fromSeries) {
    etPrescribedSeries_ = TimeSeries::fromFile(config.resolvePath(et.prescribed.file));
  }

  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const std::size_t ny2 = static_cast<std::size_t>(nyl) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(nxl) + 2;
  etColumn_ = Field2<int>("gw_et_column", ny2, nx2);
  etKtop_ = Field2<int>("gw_et_ktop", ny2, nx2);
  etThetaFc_ = Field2<real_t>("gw_et_theta_fc", ny2, nx2);
  etAlbedo_ = Field2<real_t>("gw_et_albedo", ny2, nx2);
  qVap_ = Field2<real_t>("gw_et_qvap", ny2, nx2);
  wetClockLast_ = Field2<real_t>("gw_et_wet_clock_last", ny2, nx2);
  soilTime_ = Field2<real_t>("gw_et_soil_time", ny2, nx2);
  etAppliedDepth_ = Field2<real_t>("gw_et_applied_depth", ny2, nx2);
  etUnrealisedDepth_ = Field2<real_t>("gw_et_unrealised_depth", ny2, nx2);
  Kokkos::deep_copy(qVap_, 0.0);
  Kokkos::deep_copy(wetClockLast_, 0.0);
  Kokkos::deep_copy(soilTime_, 0.0);
  Kokkos::deep_copy(etAppliedDepth_, 0.0);
  Kokkos::deep_copy(etUnrealisedDepth_, 0.0);

  // Per-column membership (region or all), ktop, and the field capacity
  // of the top cell's soil type (the configured theta_fc, else
  // theta(h = -3.3 m) from van Genuchten; W5.1). Fatal on overlap with a
  // configured groundwater_top BC — the same rule (and reason) as the Q4
  // zone (RichardsSolver.cpp): the two would prescribe the same face.
  const bool all = et.soilPolygon.empty();
  const std::vector<std::array<real_t, 2>> dummy = {
      {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
  const Polygon poly(all ? dummy : et.soilPolygon);
  auto hostCol = Kokkos::create_mirror_view(etColumn_);
  auto hostKtop = Kokkos::create_mirror_view(etKtop_);
  auto hostFc = Kokkos::create_mirror_view(etThetaFc_);
  auto hostAlb = Kokkos::create_mirror_view(etAlbedo_);
  auto hostCode = Kokkos::create_mirror_view(topCode_);
  Kokkos::deep_copy(hostCol, 0);
  Kokkos::deep_copy(hostKtop, 0);
  Kokkos::deep_copy(hostFc, 0.0);
  Kokkos::deep_copy(hostAlb, 0.0);
  Kokkos::deep_copy(hostCode, topCode_);
  const HostField2<int>& ktop = mesh_.ktop();
  long members = 0;
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      const real_t xc = grid_.xCenter(grid_.i0() + i - 1);
      const real_t yc = grid_.yCenter(grid_.j0() + j - 1);
      if (!all && !poly.contains(xc, yc)) {
        continue;
      }
      const auto ju = static_cast<std::size_t>(j);
      const auto iu = static_cast<std::size_t>(i);
      if (hostCode(ju, iu) != static_cast<int>(GwBcCode::NoFlux)) {
        log::fatal(log::msg()
                   << "evapotranspiration.soil: the region overlaps a configured "
                      "groundwater_top boundary condition at cell ("
                   << grid_.j0() + j - 1 << ", " << grid_.i0() + i - 1
                   << "); the two prescribe the same face");
      }
      hostCol(ju, iu) = 1;
      const int k0 = ktop(ju - 1, iu - 1);
      hostKtop(ju, iu) = k0;
      hostFc(ju, iu) = thetaFcHost_(ju, iu, static_cast<std::size_t>(k0));
      hostAlb(ju, iu) = albedoHost_(ju, iu, static_cast<std::size_t>(k0));
      ++members;
    }
  }
  Kokkos::deep_copy(etColumn_, hostCol);
  Kokkos::deep_copy(etKtop_, hostKtop);
  Kokkos::deep_copy(etThetaFc_, hostFc);
  Kokkos::deep_copy(etAlbedo_, hostAlb);

  long global = 0;
  MPI_Allreduce(&members, &global, 1, MPI_LONG, MPI_SUM, grid_.comm());
  if (global == 0) {
    log::fatal(log::msg() << "evapotranspiration.soil: the region selects no "
                             "column anywhere in the domain; check the polygon "
                             "against cell centers");
  }
}

void RichardsSolver::updateEvapotranspiration(real_t t, real_t dtg) {
  if (!etActive_) {
    return;
  }
  // The potential chain is host-evaluated once per substep (one met
  // sample serves every column; alpha_barton's alpha_1 enters per column
  // inside the humidity difference below).
  real_t prefactor = 0.0;  // rho_a / ((r_a) rho_w) [m/s per unit humidity]
  real_t qsat = 0.0;
  real_t qa = 0.0;
  real_t ePrescribed = 0.0;
  // Stage-3 T7 (W7.3): the PM soil channel's met scalars (r_s = 0; the
  // per-column soil albedo enters in-kernel; G = gFractionSoil * R_n).
  real_t pmTair = 0.0;
  real_t pmP = 0.0;
  real_t pmQa = 0.0;
  real_t pmRs = 0.0;
  real_t pmLw = 0.0;
  real_t pmRa = 1.0;
  const bool bulk = etCfg_.potential == EvapotranspirationConfig::Potential::Bulk;
  const bool pm =
      etCfg_.potential == EvapotranspirationConfig::Potential::PenmanMonteith;
  if (bulk) {
    const atm::MetSample met = met_.sample(t);
    const real_t rA = atm::aeroResistance(etAeroLaw_, met.windSpeed, etAero_);
    qsat = atm::saturatedSpecificHumidity(met.surfaceTemperatureC, met.pressureKpa);
    qa = met.airSpecificHumidity;
    prefactor = atm::airDensity(met.surfaceTemperatureC, met.pressureKpa) /
                (rA * atm::kWaterDensity);
  } else if (pm) {
    const atm::MetSample met = met_.sample(t);
    pmTair = met.airTemperatureC;
    pmP = met.pressureKpa;
    pmQa = met.airSpecificHumidity;
    pmRs = met.shortwaveIncoming;
    pmLw = met.longwaveIn;
    pmRa = atm::aeroResistance(etAeroLaw_, met.windSpeed, etAero_);
  } else {
    ePrescribed = etCfg_.prescribed.fromSeries ? etPrescribedSeries_.value(t)
                                               : etCfg_.prescribed.constant;
  }
  const real_t pmG = etCfg_.energy.gFractionSoil;

  const bool coupled = cpl_.active;
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t az = mesh_.areaZ();
  const bool limitAlpha =
      etCfg_.limiter == EvapotranspirationConfig::Limiter::AlphaBarton;
  const bool limitBeta =
      etCfg_.limiter == EvapotranspirationConfig::Limiter::BetaLeePielke;

  Field2<int> col = etColumn_, ktop = etKtop_;
  Field2<int> cplMode = coupled ? cplMode_ : Field2<int>();
  const bool veg = vegActive_;
  Field2<real_t> vegTauK = veg ? vegTau_ : Field2<real_t>();
  Field2<real_t> thetaFc = etThetaFc_, qVap = qVap_;
  Field2<real_t> albedoF = etAlbedo_;
  Field2<real_t> wetClockLast = wetClockLast_, soilTime = soilTime_;
  Field2<real_t> appliedDepth = etAppliedDepth_, unrealDepth = etUnrealisedDepth_;
  Field2<real_t> wetClockStart = coupled ? cpl_.wetClockStepStart : Field2<real_t>();
  Field2<real_t> wetFrac = coupled ? cpl_.wetFrac : Field2<real_t>();
  Field3<real_t> wcn = wcn_, wc = wc_, wcr = wcr_, wcs = wcs_;
  Field3<real_t> dz3d = mesh_.dz3d();
  // Substep window [a, b] inside the surface step starting at stepStart:
  // tSub (the coupler's clock) has already advanced to b = t.
  const real_t stepStart = coupled ? cpl_.stepStart : 0.0;

  real_t demand = 0.0;
  real_t applied = 0.0;
  Kokkos::parallel_reduce(
      "gw_et_update",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumD, real_t& sumA) {
        if (col(j, i) == 0) {
          return;
        }
        const int k0 = ktop(j, i);
        // (1) The potential soil rate E_s under the configured limiter.
        real_t eS;
        if (bulk) {
          real_t alpha1 = 1.0;
          if (limitAlpha) {
            // alpha_barton: alpha_1(w_g) INSIDE the humidity difference,
            // with w_g the v2.0 surface extrapolation of the top two
            // theta values (updateBulkEvaporation's lines verbatim; et1
            // asserts the chain bitwise against the old path).
            const real_t w0 = wcn(j, i, k0);
            const real_t w1 = (k0 + 1 < nz) ? wcn(j, i, k0 + 1) : w0;
            real_t wSurf = 1.5 * w0 - 0.5 * w1;
            if (wSurf < wcr(j, i, k0)) {
              wSurf = wcr(j, i, k0);
            }
            if (wSurf > wcs(j, i, k0)) {
              wSurf = wcs(j, i, k0);
            }
            alpha1 = atm::soilRelativeHumidity(wSurf);
          }
          eS = prefactor * (alpha1 * qsat - qa);
        } else if (pm) {
          // T7: PM with r_s = 0 at the column's top soil albedo (the
          // limiter below is beta or none; alpha_barton schema-rejected).
          eS = atm::penmanMonteithRate(pmTair, pmP, pmQa, pmRs, pmLw,
                                       albedoF(j, i), pmG, pmRa, 0.0);
        } else {
          eS = ePrescribed;
        }
        if (limitBeta) {
          // beta_lee_pielke multiplies the potential by beta(w_g, fc) on
          // the same surface-extrapolated moisture.
          const real_t w0 = wcn(j, i, k0);
          const real_t w1 = (k0 + 1 < nz) ? wcn(j, i, k0 + 1) : w0;
          real_t wSurf = 1.5 * w0 - 0.5 * w1;
          if (wSurf < wcr(j, i, k0)) {
            wSurf = wcr(j, i, k0);
          }
          if (wSurf > wcs(j, i, k0)) {
            wSurf = wcs(j, i, k0);
          }
          eS *= atm::limiterBeta(wSurf, thetaFc(j, i));
        }
        if (veg && vegTauK(j, i) < 1.0) {
          // Stage-3 T6 (W6.3): the ground rate under a canopy is scaled
          // by the Beer-Lambert transmission (exactly 1.0 when bare).
          eS *= vegTauK(j, i);
        }
        // (2) The dry fraction of this substep from the surface channel's
        // wet clock (uncoupled: always 1). W_b reconstructs the clock at
        // the substep end b = t by linear growth within the surface step;
        // the consumed-clock marker makes each wet second count once.
        real_t fDry = 1.0;
        if (coupled) {
          const real_t wB =
              wetClockStart(j, i) + wetFrac(j, i) * (t - stepStart);
          real_t fWet = (wB - wetClockLast(j, i)) / dtg;
          if (fWet < 0.0) {
            fWet = 0.0;
          }
          if (fWet > 1.0) {
            fWet = 1.0;
          }
          wetClockLast(j, i) = wB;
          fDry = 1.0 - fWet;
        }
        // (3) The applied flux: supply-capped on the top cell's start-of-
        // substep moisture; a coupled wet column (mode != 0) is the
        // surface channel's to evaporate — the soil channel stands down.
        const real_t demandRate = fDry * eS;
        real_t q = 0.0;
        if (!coupled || cplMode(j, i) == 0) {
          q = atm::supplyCap(demandRate, wcn(j, i, k0), wcr(j, i, k0),
                             wcs(j, i, k0), dz3d(j, i, k0), dtg);
        }
        qVap(j, i) = q;
        // (4) The soil-time clock: the dry fraction is covered iff the
        // channel actually ran.
        if (q != 0.0) {
          soilTime(j, i) += fDry * dtg;
        }
        // (5) The books (positive out) and the rate-output contributors.
        appliedDepth(j, i) += q * dtg;
        unrealDepth(j, i) += (demandRate - q) * dtg;
        sumD += demandRate * az * dtg;
        sumA += q * az * dtg;
      },
      demand, applied);
  etDemand_ = demand;
  etApplied_ = applied;
  audit_.etDemand += demand;
  audit_.etApplied += applied;
}

void RichardsSolver::resetEtDepths() {
  if (!etActive_) {
    return;
  }
  Kokkos::deep_copy(etAppliedDepth_, 0.0);
  Kokkos::deep_copy(etUnrealisedDepth_, 0.0);
}



void RichardsSolver::setupVegetation(const FrehgConfig& config) {
  const VegetationConfig& veg = config.vegetation;
  if (!veg.enabled) {
    return;
  }
  vegActive_ = true;
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const std::size_t ny2 = static_cast<std::size_t>(nyl) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(nxl) + 2;
  vegTau_ = Field2<real_t>("gw_veg_tau", ny2, nx2);
  vegTPotFactor_ = Field2<real_t>("gw_veg_tpot_factor", ny2, nx2);
  vegRc_ = Field2<real_t>("gw_veg_rc", ny2, nx2);
  vegAlbedo_ = Field2<real_t>("gw_veg_albedo", ny2, nx2);
  vegCol_ = Field2<int>("gw_veg_col", ny2, nx2);
  vegH1_ = Field2<real_t>("gw_veg_h1", ny2, nx2);
  vegH2_ = Field2<real_t>("gw_veg_h2", ny2, nx2);
  vegH3_ = Field2<real_t>("gw_veg_h3", ny2, nx2);
  vegH4_ = Field2<real_t>("gw_veg_h4", ny2, nx2);
  vegHasH12_ = Field2<int>("gw_veg_has_h12", ny2, nx2);
  rootFrac_ = Field3<real_t>("gw_veg_root_frac", ny2, nx2,
                             static_cast<std::size_t>(nz));
  tPot_ = Field2<real_t>("gw_veg_tpot", ny2, nx2);
  sink_ = Field3<real_t>("gw_veg_sink", ny2, nx2, static_cast<std::size_t>(nz));
  sinkStepWin_ =
      Field3<real_t>("gw_veg_sink_win", ny2, nx2, static_cast<std::size_t>(nz));
  trAppliedDepth_ = Field2<real_t>("gw_veg_tr_applied", ny2, nx2);
  trUnrealisedDepth_ = Field2<real_t>("gw_veg_tr_unrealised", ny2, nx2);
  Kokkos::deep_copy(vegTau_, 1.0);
  Kokkos::deep_copy(vegTPotFactor_, 0.0);
  Kokkos::deep_copy(vegRc_, 0.0);
  Kokkos::deep_copy(vegAlbedo_, 0.0);
  Kokkos::deep_copy(vegCol_, 0);
  Kokkos::deep_copy(rootFrac_, 0.0);
  Kokkos::deep_copy(tPot_, 0.0);
  Kokkos::deep_copy(sink_, 0.0);
  Kokkos::deep_copy(sinkStepWin_, 0.0);
  Kokkos::deep_copy(trAppliedDepth_, 0.0);
  Kokkos::deep_copy(trUnrealisedDepth_, 0.0);

  // Resolve the id map (0 = bare, k >= 1 = types[k-1]) and the optional
  // LAI override (host, static).
  std::vector<int> ids(static_cast<std::size_t>(grid_.nx()) *
                           static_cast<std::size_t>(grid_.ny()),
                       0);
  if (veg.mapFromFile) {
    ids = io::readIntRaster2D(config.resolvePath(veg.mapFile), grid_.nx(),
                              grid_.ny());
    for (const int id : ids) {
      if (id < 0 || id > static_cast<int>(veg.types.size())) {
        log::fatal(log::msg() << "vegetation.map.file: id " << id
                              << " is not 0 (bare) or a valid index into "
                              << veg.types.size() << " vegetation type(s)");
      }
    }
  } else {
    int id = -1;
    for (std::size_t n = 0; n < veg.types.size(); ++n) {
      if (veg.types[n].name == veg.mapConstantName) {
        id = static_cast<int>(n) + 1;
        break;
      }
    }
    std::fill(ids.begin(), ids.end(), id);
  }
  std::vector<real_t> laiOverride;
  if (veg.hasLaiFile) {
    laiOverride =
        io::readRaster2D(config.resolvePath(veg.laiFile), grid_.nx(), grid_.ny());
  }

  auto hTau = Kokkos::create_mirror_view(vegTau_);
  auto hFac = Kokkos::create_mirror_view(vegTPotFactor_);
  auto hRc = Kokkos::create_mirror_view(vegRc_);
  auto hAlb = Kokkos::create_mirror_view(vegAlbedo_);
  auto hCol = Kokkos::create_mirror_view(vegCol_);
  auto hH1 = Kokkos::create_mirror_view(vegH1_);
  auto hH2 = Kokkos::create_mirror_view(vegH2_);
  auto hH3 = Kokkos::create_mirror_view(vegH3_);
  auto hH4 = Kokkos::create_mirror_view(vegH4_);
  auto hHas = Kokkos::create_mirror_view(vegHasH12_);
  auto hFrac = Kokkos::create_mirror_view(rootFrac_);
  Kokkos::deep_copy(hTau, 1.0);
  Kokkos::deep_copy(hFac, 0.0);
  Kokkos::deep_copy(hRc, 0.0);
  Kokkos::deep_copy(hAlb, 0.0);
  Kokkos::deep_copy(hCol, 0);
  Kokkos::deep_copy(hFrac, 0.0);

  // Per-cell root fractions from the closed-form density integrals on
  // the terrain-following/masked column (depth measured from the ktop
  // cell's top face; shallow columns renormalise with a warning count).
  auto hDz = Kokkos::create_mirror_view(mesh_.dz3d());
  Kokkos::deep_copy(hDz, mesh_.dz3d());
  const HostField2<int>& ktopHost = mesh_.ktop();
  auto hGid = grid_.gid3Host();
  long renormalised = 0;
  long vegCells = 0;
  for (int j = 1; j <= nyl; ++j) {
    for (int i = 1; i <= nxl; ++i) {
      const auto ju = static_cast<std::size_t>(j);
      const auto iu = static_cast<std::size_t>(i);
      const std::size_t g =
          static_cast<std::size_t>(grid_.j0() + j - 1) *
              static_cast<std::size_t>(grid_.nx()) +
          static_cast<std::size_t>(grid_.i0() + i - 1);
      const int id = ids[g];
      if (id == 0) {
        continue;
      }
      const VegetationType& type = veg.types[static_cast<std::size_t>(id - 1)];
      real_t lai = type.lai;
      if (veg.hasLaiFile) {
        lai = laiOverride[g];
      }
      const real_t tau = atm::beerLambert(type.extinction, lai);
      hTau(ju, iu) = tau;
      hFac(ju, iu) = 1.0 - tau;
      hRc(ju, iu) = lai > 0.0 ? atm::canopyResistance(type.leafResistance, lai)
                              : 0.0;
      hAlb(ju, iu) = type.hasAlbedo ? type.albedo : 0.0;
      hCol(ju, iu) = 1;
      hH1(ju, iu) = type.h1;
      hH2(ju, iu) = type.h2;
      hH3(ju, iu) = type.h3;
      hH4(ju, iu) = type.h4;
      hHas(ju, iu) = type.hasH12 ? 1 : 0;
      ++vegCells;
      const bool linear = type.profile == VegetationType::RootProfile::Linear;
      const int k0 = ktopHost(ju - 1, iu - 1);
      real_t depth = 0.0;
      real_t total = 0.0;
      for (int k = k0; k < nz; ++k) {
        if (hGid(ju, iu, static_cast<std::size_t>(k)) < 0) {
          break;
        }
        const real_t dz = hDz(ju, iu, static_cast<std::size_t>(k));
        const real_t frac = atm::rootFractionIntegral(linear, type.rootDepth,
                                                      depth, depth + dz);
        hFrac(ju, iu, static_cast<std::size_t>(k)) = frac;
        total += frac;
        depth += dz;
        if (depth >= type.rootDepth) {
          break;
        }
      }
      if (total <= 0.0) {
        log::fatal(log::msg()
                   << "vegetation: the root zone of type '" << type.name
                   << "' at cell (" << grid_.j0() + j - 1 << ", "
                   << grid_.i0() + i - 1 << ") covers no active soil cell");
      }
      if (total < 1.0 - 1.0e-12) {
        ++renormalised;
      }
      for (int k = k0; k < nz; ++k) {
        hFrac(ju, iu, static_cast<std::size_t>(k)) /= total;
      }
    }
  }
  Kokkos::deep_copy(vegTau_, hTau);
  Kokkos::deep_copy(vegTPotFactor_, hFac);
  Kokkos::deep_copy(vegRc_, hRc);
  Kokkos::deep_copy(vegAlbedo_, hAlb);
  Kokkos::deep_copy(vegCol_, hCol);
  Kokkos::deep_copy(vegH1_, hH1);
  Kokkos::deep_copy(vegH2_, hH2);
  Kokkos::deep_copy(vegH3_, hH3);
  Kokkos::deep_copy(vegH4_, hH4);
  Kokkos::deep_copy(vegHasH12_, hHas);
  Kokkos::deep_copy(rootFrac_, hFrac);

  long sums[2] = {vegCells, renormalised};
  long global[2] = {0, 0};
  MPI_Allreduce(sums, global, 2, MPI_LONG, MPI_SUM, grid_.comm());
  if (global[0] == 0) {
    // An all-bare map is LEGAL (the et5(i) vegetation-absent pin): warn
    // and leave the module inert — every guard reads vegActive_, so the
    // run is bitwise the vegetation-absent one.
    if (grid_.rank() == 0) {
      log::warn("vegetation: the map selects no vegetated cell anywhere in "
                "the domain (all ids are 0); the block is inert");
    }
    vegActive_ = false;
    return;
  }
  if (global[1] > 0 && grid_.rank() == 0) {
    log::warn(log::msg() << "vegetation: " << global[1]
                         << " column(s) are shallower than the root depth; "
                            "their root fractions were renormalised");
  }
  // ET consumers read the potential chain: a vegetated run without the
  // soil channel still needs the met/prescribed machinery staged.
  if (!etActive_) {
    etCfg_ = config.evapotranspiration;
    if (etCfg_.potential == EvapotranspirationConfig::Potential::Bulk ||
        etCfg_.potential == EvapotranspirationConfig::Potential::PenmanMonteith) {
      met_ = atm::MetForcing(config.atmosphere, config);
      etAeroLaw_ = config.atmosphere.aeroResistance;
      etAero_ = {config.atmosphere.aeroRefHeight, config.atmosphere.aeroZ0m,
                 config.atmosphere.aeroZ0h, config.atmosphere.aeroDisplacement};
    } else if (etCfg_.prescribed.fromSeries) {
      etPrescribedSeries_ =
          TimeSeries::fromFile(config.resolvePath(etCfg_.prescribed.file));
    }
  }
}

void RichardsSolver::resetTrDepths() {
  if (!vegActive_) {
    return;
  }
  Kokkos::deep_copy(trAppliedDepth_, 0.0);
  Kokkos::deep_copy(trUnrealisedDepth_, 0.0);
  Kokkos::deep_copy(sinkStepWin_, 0.0);
}




void RichardsSolver::updateTranspiration(real_t t, real_t dtg) {
  if (!vegActive_) {
    return;
  }
  // T_pot = (1 - tau) E_x(r_s = r_c) per column (W6.3): bulk evaluates
  // the canopy at T_s with the column's r_c inside the chain; PM (T7)
  // evaluates the q-form at T_air with r_s = r_c and the column's canopy
  // albedo; prescribed uses E_pot and ignores r_s. One met sample per
  // substep.
  const bool bulk = etCfg_.potential == EvapotranspirationConfig::Potential::Bulk;
  const bool pm =
      etCfg_.potential == EvapotranspirationConfig::Potential::PenmanMonteith;
  real_t qsat = 0.0;
  real_t qa = 0.0;
  real_t rhoOverRw = 0.0;
  real_t rA = 1.0;
  real_t ePrescribed = 0.0;
  real_t pmTair = 0.0;
  real_t pmP = 0.0;
  real_t pmRs = 0.0;
  real_t pmLw = 0.0;
  if (bulk) {
    const atm::MetSample met = met_.sample(t);
    rA = atm::aeroResistance(etAeroLaw_, met.windSpeed, etAero_);
    qsat = atm::saturatedSpecificHumidity(met.surfaceTemperatureC, met.pressureKpa);
    qa = met.airSpecificHumidity;
    rhoOverRw = atm::airDensity(met.surfaceTemperatureC, met.pressureKpa) /
                atm::kWaterDensity;
  } else if (pm) {
    const atm::MetSample met = met_.sample(t);
    rA = atm::aeroResistance(etAeroLaw_, met.windSpeed, etAero_);
    pmTair = met.airTemperatureC;
    pmP = met.pressureKpa;
    qa = met.airSpecificHumidity;
    pmRs = met.shortwaveIncoming;
    pmLw = met.longwaveIn;
  } else {
    ePrescribed = etCfg_.prescribed.fromSeries ? etPrescribedSeries_.value(t)
                                               : etCfg_.prescribed.constant;
  }
  const real_t pmG = etCfg_.energy.gFractionCanopy;

  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t az = mesh_.areaZ();
  Field2<int> col = vegCol_, hasH12 = vegHasH12_;
  Field2<real_t> fac = vegTPotFactor_, rc = vegRc_;
  Field2<real_t> vegAlb = vegAlbedo_;
  Field2<real_t> h1 = vegH1_, h2 = vegH2_, h3 = vegH3_, h4 = vegH4_;
  Field2<real_t> tPot = tPot_;
  Field2<real_t> trApplied = trAppliedDepth_, trUnreal = trUnrealisedDepth_;
  Field3<real_t> frac = rootFrac_, sink = sink_, sinkWin = sinkStepWin_;
  Field3<real_t> hn = hn_, wcn = wcn_, wcr = wcr_;
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<PetscInt> gid = grid_.gid3();

  real_t demand = 0.0;
  real_t applied = 0.0;
  Kokkos::parallel_reduce(
      "gw_veg_transpiration",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumD, real_t& sumA) {
        if (col(j, i) == 0) {
          return;
        }
        real_t tp;
        if (bulk) {
          tp = fac(j, i) * rhoOverRw / (rA + rc(j, i)) * (qsat - qa);
        } else if (pm) {
          tp = fac(j, i) * atm::penmanMonteithRate(pmTair, pmP, qa, pmRs, pmLw,
                                                   vegAlb(j, i), pmG, rA,
                                                   rc(j, i));
        } else {
          tp = fac(j, i) * ePrescribed;
        }
        if (tp < 0.0) {
          tp = 0.0;  // a condensing canopy does not transpire
        }
        tPot(j, i) = tp;
        // The volumetric sink per root cell (O7: no inundation test):
        //   S_k = min(alpha_F(h^n_k) rootFrac_k T_pot / dz_k,
        //             (theta^n_k - theta_r) / dtg)   [1/s].
        for (int k = 0; k < nz; ++k) {
          if (gid(j, i, k) < 0 || frac(j, i, k) == 0.0) {
            sink(j, i, k) = 0.0;
            continue;
          }
          const real_t alpha =
              atm::feddes(hn(j, i, k), h1(j, i), h2(j, i), h3(j, i), h4(j, i),
                          hasH12(j, i) != 0);
          const real_t want = alpha * frac(j, i, k) * tp / dz3d(j, i, k);
          const real_t supply =
              Kokkos::fmax(0.0, wcn(j, i, k) - wcr(j, i, k)) / dtg;
          const real_t sValue = Kokkos::fmin(want, supply);
          sink(j, i, k) = sValue;
          sinkWin(j, i, k) += sValue * dtg;
          const real_t vol = az * dz3d(j, i, k);
          sumD += frac(j, i, k) * tp * az * dtg;
          sumA += sValue * vol * dtg;
          trApplied(j, i) += sValue * dz3d(j, i, k) * dtg;
          trUnreal(j, i) +=
              (frac(j, i, k) * tp - sValue * dz3d(j, i, k)) * dtg;
        }
      },
      demand, applied);
  trDemand_ = demand;
  trApplied_ = applied;
  audit_.trDemand += demand;
  audit_.trApplied += applied;
}


}  // namespace frehg::gw
