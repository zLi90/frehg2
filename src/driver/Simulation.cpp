/// \file Simulation.cpp
/// \brief Driver implementation: setup, time loops, outputs, restart.

#include "driver/Simulation.hpp"

#include "core/Logger.hpp"
#include "core/Timer.hpp"

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <type_traits>

namespace frehg::driver {

namespace {

/// Attributes for the §7 surface datasets.
io::VarMeta surfaceVarMeta(const std::string& var) {
  if (var == "eta") {
    return {"m", "free-surface elevation"};
  }
  if (var == "depth") {
    return {"m", "water depth"};
  }
  if (var == "uu") {
    return {"m s-1", "x-face velocity"};
  }
  if (var == "vv") {
    return {"m s-1", "y-face velocity"};
  }
  if (var == "seepage") {
    return {"m s-1", "surface-applied seepage rate (positive upward; a "
                     "PER-STEP SAMPLE — integrate seepage_rate instead)"};
  }
  if (var == "water_store") {
    return {"m3", "paired water store V_s: thin-layer water parked with its "
                  "salt, inert until release (Stage-3 T9)"};
  }
  return {"", var};
}

/// Attributes for the §7 groundwater datasets.
io::VarMeta groundwaterVarMeta(const std::string& var) {
  if (var == "hydraulic_head") {
    return {"m", "pressure head"};
  }
  if (var == "water_content") {
    return {"1", "volumetric water content"};
  }
  if (var == "qx") {
    return {"m s-1", "x-face Darcy flux (positive toward -x)"};
  }
  if (var == "qy") {
    return {"m s-1", "y-face Darcy flux (positive toward -y)"};
  }
  if (var == "qz") {
    return {"m s-1", "lower-face Darcy flux (positive upward)"};
  }
  if (var == "realloc_created") {
    return {"m3", "cumulative post-allocation creation per cell (Stage-3 W1.3)"};
  }
  if (var == "realloc_dropped") {
    return {"m3", "cumulative post-allocation discard per cell (Stage-3 W1.3)"};
  }
  if (var == "vloss") {
    return {"m3", "cumulative finalize-clamp volume per cell (signed; + = removed)"};
  }
  return {"", var};
}

/// Attributes for the §7 transport datasets.
io::VarMeta transportVarMeta(const std::string& var) {
  if (var == "concentration") {
    return {"psu", "subsurface scalar concentration"};
  }
  if (var == "concentration_surface") {
    return {"psu", "surface scalar concentration"};
  }
  return {"", var};
}

/// Attributes for the §7 temperature datasets (v2 Q5).
io::VarMeta temperatureVarMeta(const std::string& var) {
  if (var == "temperature") {
    return {"degC", "subsurface temperature"};
  }
  if (var == "temperature_surface") {
    return {"degC", "surface temperature"};
  }
  return {"", var};
}

bool resolveGpuAware(const RuntimeConfig& runtime) {
  switch (runtime.gpuAwareMpi) {
    case RuntimeConfig::GpuAwareMpi::On:
      return true;
    case RuntimeConfig::GpuAwareMpi::Off:
      return false;
    case RuntimeConfig::GpuAwareMpi::Auto:
      break;
  }
#if defined(FREHG_GPU_AWARE_MPI_DEFAULT)
  return true;
#else
  return false;
#endif
}

}  // namespace

Simulation::Simulation(MPI_Comm comm, const FrehgConfig& config,
                       const std::string& inputPath)
    : config_(config), grid_(comm, config.domain) {
  if (!config_.modules.surfaceWater && !config_.modules.groundwater) {
    log::fatal("no module enabled: set modules.surface_water or modules.groundwater");
  }
  // GPU release status (v2 plan §2B.4): device execution is 'experimental'
  // until the owner-run p6 acceptance bundle passes on real hardware —
  // compile/link-verified and CPU-physics-verified, but CI never executed
  // this backend. The notice is part of the release contract, not a warning
  // about a known defect.
  if constexpr (!std::is_same_v<MemSpace, Kokkos::HostSpace>) {
    log::info(log::msg() << "device build (" << Kokkos::DefaultExecutionSpace::name()
                         << "): GPU execution is EXPERIMENTAL pending the p6 "
                            "acceptance bundle (docs/developer-guide/"
                            "gpu-acceptance.md, v2 plan §2B.4)");
  }
  const bool coupled = config_.modules.surfaceWater && config_.modules.groundwater;

  if (!config_.solver.petscOptionsFile.empty()) {
    const std::string path = config_.resolvePath(config_.solver.petscOptionsFile);
    const PetscErrorCode err =
        PetscOptionsInsertFile(comm, nullptr, path.c_str(), PETSC_TRUE);
    if (err != PETSC_SUCCESS) {
      log::fatal(log::msg() << "failed to read PETSc options file '" << path << "'");
    }
  }

  halo_ = std::make_unique<HaloExchanger>(grid_, resolveGpuAware(config_.runtime));

  if (config_.modules.groundwater) {
    // The subsurface mesh derives the active-column mask from the
    // bathymetry; surface-only runs mask nothing (every column active from
    // the surface).
    mesh_ = std::make_unique<gw::TerrainMetric>(grid_, config_, *halo_);
    grid_.buildGlobalIds(mesh_->ktop());
  } else {
    HostField2<int> ktop("driver_ktop", static_cast<std::size_t>(grid_.nyLocal()),
                         static_cast<std::size_t>(grid_.nxLocal()));
    Kokkos::deep_copy(ktop, 0);
    grid_.buildGlobalIds(ktop);
  }

  boundaries_ = std::make_unique<BoundarySet>(grid_, config_.boundaryConditions,
                                              config_.configDir);
  if (config_.modules.surfaceWater) {
    if (coupled && grid_.activeCount2Global() <
                       static_cast<long>(grid_.nx()) * static_cast<long>(grid_.ny())) {
      // The surface free-surface system assembles a row for every (j, i)
      // cell; fully masked columns (local bed below the subsurface box
      // bottom) have no surface row. No P3 benchmark reaches this (b5 and
      // b6 keep every column active), so it fails loudly instead of
      // running an untested path (plan §11.1 rule 3).
      log::fatal(
          "coupled runs require every column active in the surface system "
          "(fully masked columns are not part of the P3-validated scope)");
    }
    surface_ = std::make_unique<swe::SurfaceSolver>(grid_, config_, *boundaries_, *halo_);
  }
  if (config_.modules.groundwater) {
    gw_ = std::make_unique<gw::RichardsSolver>(grid_, config_, *mesh_, *boundaries_, *halo_);
    dtg_ = config_.groundwater.timestep.dtInit;
  }
  if (coupled) {
    coupler_ = std::make_unique<coupling::Coupler>(grid_, config_, *surface_, *gw_);
  }
  if (surface_ && gw_ && gw_->vegActive() && surface_->etActive()) {
    // Stage-3 T6 (W6.3): the canopy transmission scales the open-water
    // channel's rate per cell (exactly 1.0 on bare cells — the T5
    // outputs stay bitwise without vegetation).
    surface_->attachCanopyFactor(gw_->vegCanopyTau());
  }
  if (config_.modules.transport || config_.modules.temperature) {
    buildTransport();
  }

  output_ = std::make_unique<io::Hdf5Output>(grid_, config_.output.filename, config_.rawText);
  // Run provenance record (v2 plan §2A): static sections captured now,
  // rewritten at every output flush and finalized after run().
  runRecord_ = std::make_unique<io::RunRecord>(
      grid_.comm(), config_,
      inputPath.empty() ? std::string("<in-memory>") : inputPath, *boundaries_,
      grid_.px(), grid_.py());
  // Non-fatal validation advisories (Stage-3 W0.9): visible in the log and
  // persisted in the record's `warnings:` list (seeded by the constructor).
  if (grid_.rank() == 0) {
    for (const std::string& warning : config_.validationWarnings) {
      log::warn(log::msg() << "configuration: " << warning);
    }
  }
  // Stage-3 T3 (W3.5): seam-discontinuity diagnostic on each periodic
  // axis. A periodic wrap continues the field across the seam; a jump
  // larger than the largest interior neighbor jump means the input is a
  // discontinuity, not a periodic continuation (the H2-class hazard of
  // plan §13.B). Checked on the bed elevation of whichever module carries
  // it (the halo exchange has already wrapped the ghost columns/rows).
  if (grid_.periodicX() || grid_.periodicY()) {
    const auto diagnose = [&](const Field2<real_t>& field, const char* name) {
      auto host = Kokkos::create_mirror_view(field);
      Kokkos::deep_copy(host, field);
      const std::size_t nyl = static_cast<std::size_t>(grid_.nyLocal());
      const std::size_t nxl = static_cast<std::size_t>(grid_.nxLocal());
      // [0] seam-x, [1] seam-y, [2] interior max jump.
      real_t local[3] = {0.0, 0.0, 0.0};
      for (std::size_t j = 1; j <= nyl; ++j) {
        for (std::size_t i = 1; i + 1 <= nxl; ++i) {
          local[2] = std::max(local[2], std::fabs(host(j, i + 1) - host(j, i)));
        }
      }
      for (std::size_t i = 1; i <= nxl; ++i) {
        for (std::size_t j = 1; j + 1 <= nyl; ++j) {
          local[2] = std::max(local[2], std::fabs(host(j + 1, i) - host(j, i)));
        }
      }
      if (grid_.periodicX() && grid_.i0() == 0) {
        for (std::size_t j = 1; j <= nyl; ++j) {
          // The wrapped halo column 0 carries global column nx-1.
          local[0] = std::max(local[0], std::fabs(host(j, 1) - host(j, 0)));
        }
      }
      if (grid_.periodicY() && grid_.j0() == 0) {
        for (std::size_t i = 1; i <= nxl; ++i) {
          local[1] = std::max(local[1], std::fabs(host(1, i) - host(0, i)));
        }
      }
      real_t global[3] = {0.0, 0.0, 0.0};
      MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_MAX, grid_.comm());
      if (grid_.rank() != 0) {
        return;
      }
      for (int axis = 0; axis < 2; ++axis) {
        if (global[axis] > global[2] && global[axis] > 0.0) {
          std::ostringstream msg;
          msg << "domain.periodic." << (axis == 0 ? "x" : "y") << ": " << name
              << " jumps by " << global[axis] << " across the "
              << (axis == 0 ? "x" : "y") << " seam (largest interior jump "
              << global[2] << "); the seam is a discontinuity, not a "
              << "periodic continuation";
          log::warn(msg.str());
          runRecord_->addWarning(msg.str());
        }
      }
    };
    if (surface_) {
      diagnose(surface_->bottom(), "bottom_elevation");
    } else if (mesh_) {
      diagnose(mesh_->bath(), "bottom_elevation");
    }
  }
  {
    HostField2<real_t> bottomInterior("driver_bottom",
                                      static_cast<std::size_t>(grid_.nyLocal()),
                                      static_cast<std::size_t>(grid_.nxLocal()));
    const Field2<real_t>& bottomField =
        config_.modules.surfaceWater ? surface_->bottom() : mesh_->bath();
    const real_t shift = config_.modules.surfaceWater ? surface_->elevationOffset()
                                                      : mesh_->elevationOffset();
    auto bottomHost = Kokkos::create_mirror_view(bottomField);
    Kokkos::deep_copy(bottomHost, bottomField);
    for (int j = 0; j < grid_.nyLocal(); ++j) {
      for (int i = 0; i < grid_.nxLocal(); ++i) {
        bottomInterior(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) =
            bottomHost(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1) - shift;
      }
    }
    if (config_.modules.groundwater) {
      output_->writeGridMeta(bottomInterior, mesh_->zCellHost());
    } else {
      output_->writeGridMeta(bottomInterior);
    }
  }
  checkpoint_ = std::make_unique<io::Checkpoint>(*output_);

  // Stage-3 T4 (W4.1/W4.2): interval-mean rate outputs. Only requested
  // variables allocate; today's sole producer is the coupler's per-step
  // gain (seepage_rate).
  rates_ = std::make_unique<io::IntervalRates>(grid_, config_.time.tStart);
  {
    const auto& sv = config_.output.surfaceVariables;
    if (coupler_ && std::find(sv.begin(), sv.end(), "seepage_rate") != sv.end()) {
      rates_->add("surface", "seepage_rate",
                  {"m s-1", "interval-mean coupled surface gain (applied "
                            "seepage + bounce-back + vent deposits, positive "
                            "upward)"});
    }
    // Stage-3 T9 (W9.1): the per-cell counted thin-layer removal (P1 + P2
    // always; + P3's positive part and the sweep under the paired store) —
    // the per-cell, per-interval form of the film_removed column.
    if (surface_ && std::find(sv.begin(), sv.end(), "film_removal_rate") != sv.end()) {
      rates_->add("surface", "film_removal_rate",
                  {"m s-1", "interval-mean thin-layer removal (cflLimiter "
                            "isolated + small rules; plus the clamp film and "
                            "the end-of-step sweep under the paired store)"});
    }
    // Stage-3 T5 (W5.10): the /atmosphere interval-mean ET rates.
    const auto& av = config_.output.atmosphereVariables;
    const auto wants = [&av](const char* name) {
      return std::find(av.begin(), av.end(), name) != av.end();
    };
    if (surface_ && surface_->etActive() && wants("evaporation_rate")) {
      rates_->add("atmosphere", "evaporation_rate",
                  {"m s-1", "interval-mean open-water evaporation (positive "
                            "out; negative = condensation)"});
    }
    if (gw_ && gw_->vegActive()) {
      if (wants("transpiration_rate")) {
        rates_->add("atmosphere", "transpiration_rate",
                    {"m s-1", "interval-mean applied transpiration "
                              "(positive out)"});
      }
      if (wants("transpiration_unrealised_rate")) {
        rates_->add("atmosphere", "transpiration_unrealised_rate",
                    {"m s-1", "interval-mean unrealised transpiration "
                              "demand (demand minus applied)"});
      }
      const auto& gv = config_.output.groundwaterVariables;
      if (std::find(gv.begin(), gv.end(), "root_uptake") != gv.end()) {
        // Stage-3 T6 (S3-A7's Field3 leg, first consumer): the
        // interval-mean volumetric uptake sink.
        rates_->add3("groundwater", "root_uptake",
                     {"s-1", "interval-mean root-uptake sink"});
      }
    }
    if (gw_ && gw_->etActive()) {
      if (wants("soil_evaporation_rate")) {
        rates_->add("atmosphere", "soil_evaporation_rate",
                    {"m s-1", "interval-mean applied soil evaporation "
                              "(positive out)"});
      }
      if (wants("soil_evaporation_unrealised_rate")) {
        rates_->add("atmosphere", "soil_evaporation_unrealised_rate",
                    {"m s-1", "interval-mean unrealised soil-evaporation "
                              "demand (demand minus applied)"});
      }
    }
  }

  for (const MonitorConfig& mc : config_.output.monitors) {
    monitors_.push_back(std::make_unique<io::Monitor>(*output_, mc));
  }
  // The volume budgets of plan §9 are recorded as built-in monitor tables:
  // /monitor/mass_audit for the surface module (columns volume, rain,
  // evaporation, boundary_outflow, bc_inflow) and /monitor/gw_mass_audit for
  // the subsurface (columns volume, boundary_in, ss_storage, realloc, realloc_dropped,
  // vloss); cumulative volumes in m^3, reduced over all ranks.
  if (config_.modules.surfaceWater) {
    MonitorConfig auditConfig;
    auditConfig.name = "mass_audit";
    auditConfig.i = 0;
    auditConfig.j = 0;
    auditConfig.variables = {"volume", "rain", "evaporation", "boundary_outflow", "bc_inflow"};
    if (coupled) {
      // Cumulative net volume the surface received from the subsurface
      // (applied seepage + the bounce-back and vent deposits).
      auditConfig.variables.push_back("seepage");
    }
    // Cumulative volume created by the legacy below-bed clamp — a defect
    // measure (docs/theory/surface-water.md), audited so the closure
    // identity holds to its bound with the defect visible, not hidden in
    // the residual. Appended last so earlier columns keep their indices.
    auditConfig.variables.push_back("clamped");
    if (coupled) {
      // Stage-3 T2 (W2.4): cumulative coupled deposit on prescribed-stage
      // cells (+ = external to the surface: removed through the reservoir
      // under closure external, or deleted unbooked-no-more by the next
      // stage enforcement under zero_gradient). Appended last.
      auditConfig.variables.push_back("stage_deposit");
    }
    // Stage-3 T9 (W9.1, plan §12A): the cumulative counted thin-layer
    // removal (P1 + P2 in every mode; + P3's positive part + the sweep
    // under the paired store) and the INSTANTANEOUS sub-min_depth film
    // volume still in eta (dry to updateDepth, so outside `volume`). With
    // them the surface identity closes at every row:
    //   discard  D(volume + film_pending) = flows + clamped - stage_deposit
    //                                       - film_removed
    //   store    D(volume + store_volume) = flows + clamped - stage_deposit
    // (film_pending is identically 0 after the W9.2 sweep). Appended last
    // (plan §2.6 exception 1).
    auditConfig.variables.push_back("film_removed");
    auditConfig.variables.push_back("film_pending");
    // T9 (S3-A13): the instantaneous water on prescribed-stage cells — the
    // reservoir storage the identity's left-hand side subtracts.
    auditConfig.variables.push_back("stage_volume");
    if (surface_->pairedStore()) {
      // W9.2: the instantaneous parked volume and the cumulative release.
      auditConfig.variables.push_back("store_volume");
      auditConfig.variables.push_back("store_release");
    }
    massAudit_ = std::make_unique<io::Monitor>(*output_, auditConfig);
  }
  if (config_.modules.groundwater) {
    MonitorConfig auditConfig;
    auditConfig.name = "gw_mass_audit";
    auditConfig.i = 0;
    auditConfig.j = 0;
    auditConfig.variables = {"volume", "boundary_in", "ss_storage", "realloc",
                             "realloc_dropped", "vloss"};
    if (config_.groundwater.evaporation.enabled) {
      // v2 Q4: cumulative actual evaporated volume through the bulk-soil-
      // evaporation zone (appended so earlier columns keep their indices;
      // the g5 gate reads it as column 7).
      auditConfig.variables.push_back("evaporation");
    }
    if (config_.evapotranspiration.enabled && config_.evapotranspiration.soil) {
      // Stage-3 T5 (W5.10): the ET soil channel's cumulative applied and
      // unrealised volumes [m^3]; soil_evaporation is a sub-account of
      // boundary_in (the closure identity is unchanged).
      auditConfig.variables.push_back("soil_evaporation");
      auditConfig.variables.push_back("soil_evaporation_unrealised");
    }
    if (config_.vegetation.enabled) {
      // Stage-3 T6 (W6.6): the transpiration sink's cumulative applied
      // and unrealised volumes [m^3]. The closure identity becomes
      // dV = boundary_in - ss_storage + realloc - vloss - transpiration
      // (the sink leaves through no boundary face).
      auditConfig.variables.push_back("transpiration");
      auditConfig.variables.push_back("transpiration_unrealised");
    }
    gwMassAudit_ = std::make_unique<io::Monitor>(*output_, auditConfig);
  }
  if (transport_) {
    // Scalar-mass budget table (plan §10 P4): instantaneous masses plus the
    // cumulative exchange/boundary/adjustment terms, reduced over ranks —
    // the transport analogue of /monitor/mass_audit.
    MonitorConfig auditConfig;
    auditConfig.name = "transport_audit";
    auditConfig.i = 0;
    auditConfig.j = 0;
    auditConfig.variables = {"surf_mass",     "subs_mass",   "exchange",
                             "surf_source",   "surf_boundary", "surf_adjust",
                             "surf_anchor",   "subs_boundary", "subs_adjust",
                             "subs_anchor"};
    if (transport_->hasSurfaceStore()) {
      // Stage-3 T5 (W5.7): the INSTANTANEOUS store total [scalar m^3],
      // appended last. Identity: Delta(surf_mass + surf_store) closes on
      // the same RHS as before (the store transfers are internal).
      auditConfig.variables.push_back("surf_store");
    }
    transportAudit_ = std::make_unique<io::Monitor>(*output_, auditConfig);
  }
  if (temperature_) {
    // The heat ledger (v2 Q5): the transport-audit identity on the
    // temperature scalar, with the atmospheric exchange in its own
    // surf_atmos column (the g7(a) closure criterion). Masses are
    // K m^3 (surface) and K m^3 on the theta + kappa basis (subsurface).
    MonitorConfig auditConfig;
    auditConfig.name = "temperature_audit";
    auditConfig.i = 0;
    auditConfig.j = 0;
    auditConfig.variables = {"surf_heat",     "subs_heat",     "exchange",
                             "surf_source",   "surf_boundary", "surf_adjust",
                             "surf_anchor",   "surf_atmos",    "subs_boundary",
                             "subs_adjust",   "subs_anchor"};
    if (temperature_->hasSurfaceStore()) {
      // Stage-3 T9 (R9): the instantaneous parked heat, appended last —
      // the heat identity's left-hand side gains it as the salt one does.
      auditConfig.variables.push_back("surf_store");
    }
    temperatureAudit_ = std::make_unique<io::Monitor>(*output_, auditConfig);
  }
}

void Simulation::buildTransport() {
  transport::SurfaceWiring surfWiring;
  if (surface_) {
    surfWiring.active = true;
    surfWiring.minDepth = surface_->minDepth();
    surfWiring.dept = surface_->depth();
    surfWiring.etan = surface_->etaStart();
    surfWiring.bottom = surface_->bottom();
    surfWiring.uu = surface_->uu();
    surfWiring.vv = surface_->vv();
    surfWiring.fu = surface_->flowRateX();
    surfWiring.fv = surface_->flowRateY();
    surfWiring.asx = surface_->faceAreaX();
    surfWiring.asy = surface_->faceAreaY();
    surfWiring.rainMask = surface_->rainApplyMask();
    if (surface_->etActive()) {
      // Stage-3 T5 (W5.7 D8): the per-cell rain/evaporation depths.
      surfWiring.etActive = true;
      surfWiring.rainAdded = surface_->rainAdded();
      surfWiring.evapRemoved = surface_->evapRemoved();
    }
    if (surface_->solveVolumesActive()) {
      // Stage-3 T10 (W10.2): the solve's continuity face volumes.
      surfWiring.solveVolumes = true;
      surfWiring.qxSolve = surface_->solveVolumeX();
      surfWiring.qySolve = surface_->solveVolumeY();
      surfWiring.isEtaBc = surface_->etaBcMask();
    }
    if (surface_->pairedStore()) {
      // Stage-3 T9 (W9.3): the paired store's per-step surface books.
      surfWiring.pairedStore = true;
      surfWiring.parkedVol = surface_->parkedVolume();
      surfWiring.releasedVol = surface_->releasedVolume();
      surfWiring.evapStep = surface_->evaporatedVolume();
    }
  }
  transport::SubsurfaceWiring subsWiring;
  // transport.subsurface_update: consistent reads the groundwater module's
  // window volumes (all substeps, post-allocation transfers on every face
  // crossed); switch the accumulation on before wiring the views.
  const bool transportWindow =
      gw_ && config_.modules.transport && config_.transport.consistentSubsurfaceUpdate;
  if (transportWindow) {
    gw_->enableTransportWindow();
  }
  if (gw_) {
    subsWiring.active = true;
    subsWiring.wc = gw_->waterContent();
    subsWiring.wcn = gw_->waterContentStart();
    subsWiring.qx = gw_->fluxXVolumetric();
    subsWiring.qy = gw_->fluxYVolumetric();
    subsWiring.qzF = gw_->fluxZFaceVolumetric();
    subsWiring.kx = gw_->faceConductivityX();
    subsWiring.ky = gw_->faceConductivityY();
    subsWiring.kzF = gw_->faceConductivityZFace();
    subsWiring.wcs = gw_->soilThetaS();
    subsWiring.ksz = gw_->soilKsz();
    subsWiring.dz3d = mesh_->dz3d();
    subsWiring.ax = mesh_->areaX();
    subsWiring.ay = mesh_->areaY();
    subsWiring.cosx = mesh_->cosX();
    subsWiring.cosy = mesh_->cosY();
    subsWiring.az = mesh_->areaZ();
    subsWiring.topCode = gw_->topBcCode();
    subsWiring.topValue = gw_->topBcValue();
    subsWiring.sideCodeYp = gw_->sideBcCodeYp();
    subsWiring.sideCodeYm = gw_->sideBcCodeYm();
    subsWiring.sideCodeXm = gw_->sideBcCodeXm();
    subsWiring.sideCodeXp = gw_->sideBcCodeXp();
    if (transportWindow) {
      subsWiring.qxWindow = gw_->windowVolumeX();
      subsWiring.qyWindow = gw_->windowVolumeY();
      subsWiring.qzWindow = gw_->windowVolumeZ();
    }
    if (gw_->etActive()) {
      // Stage-3 T5 (W5.8): the top-cell vapour concentration factor.
      subsWiring.etActive = true;
      subsWiring.etColumn = gw_->etColumn();
      if (transportWindow) {
        subsWiring.vapWindow = gw_->windowVapourVolume();
      }
    }
    if (gw_->vegActive() && transportWindow) {
      // Stage-3 T6 (W6.5): the per-cell root-uptake factor.
      subsWiring.vegActive = true;
      subsWiring.sinkWindow = gw_->windowSinkVolume();
    }
  }
  transport::CouplingWiring cplWiring;
  if (coupler_) {
    cplWiring.active = true;
    cplWiring.qss = coupler_->seepageRate();
    cplWiring.gain = coupler_->gain();
    cplWiring.bounce = coupler_->bounce();
    if (transportWindow) {
      cplWiring.ventVolume = gw_->windowVentVolume();
    }
    // Stage-3 T5 (W5.9): the held-seepage salt path — ET with the store.
    if (config_.evapotranspiration.enabled && config_.modules.transport &&
        config_.transport.surfaceDryout == TransportConfig::SurfaceDryout::Store) {
      coupler_->enableHeldSeepageSalt();
      cplWiring.heldSeepage = true;
      cplWiring.exfNew = coupler_->exfNew();
      cplWiring.seepAccum = coupler_->seepAccum();
    }
  }
  if (config_.modules.transport) {
    transport_ = std::make_unique<transport::ScalarSolver>(
        grid_, config_, *boundaries_, *halo_, surfWiring, subsWiring, cplWiring,
        transport::ScalarSpec::salinity(config_));
  }
  if (config_.modules.temperature) {
    temperature_ = std::make_unique<transport::ScalarSolver>(
        grid_, config_, *boundaries_, *halo_, surfWiring, subsWiring, cplWiring,
        transport::ScalarSpec::temperature(config_));
  }
  if (gw_) {
    // Baroclinic activation (plan §10 P4; thermal term v2 Q5): the
    // density/viscosity ratios read the registered scalars at every
    // subsurface step.
    gw_->attachScalar(
        transport_ ? transport_->subsurfaceScalar() : Field3<real_t>(),
        (transport_ && surface_) ? transport_->surfaceScalar() : Field2<real_t>(),
        config_.groundwater.densityCoupling.enabled);
    if (temperature_) {
      gw_->attachTemperature(
          temperature_->subsurfaceScalar(),
          surface_ ? temperature_->surfaceScalar() : Field2<real_t>());
    }
    gw_->setDensityCoefficients(config_.groundwater.densityCoupling);
  }
}

void Simulation::stepTransport(real_t t, real_t dt, real_t dtgLast) {
  if (!transport_ && !temperature_) {
    return;
  }
  const real_t rain = surface_ ? surface_->currentRain() : 0.0;
  const real_t evap = surface_ ? surface_->currentEvaporation() : 0.0;
  // Salinity first, temperature second (the legacy scalar block order
  // generalized; the baroclinic ratios read both AFTER the transport
  // block, inside the next flow step, so the order only affects which
  // scalar sees the other's pre-step ghosts — neither does).
  if (transport_) {
    transport_->step(t, dt, dtgLast, rain, evap);
  }
  if (temperature_) {
    temperature_->step(t, dt, dtgLast, rain, evap);
  }
  if (gw_) {
    // The window volumes are consumed (no-op unless the consistent
    // subsurface update switched the accumulation on).
    gw_->resetTransportWindow();
  }
}

io::Checkpoint::Fields2 Simulation::checkpointFields2() const {
  io::Checkpoint::Fields2 fields;
  if (surface_) {
    fields = {{"eta", surface_->eta()},
              {"etan", surface_->etaStart()},
              {"uu", surface_->uu()},
              {"vv", surface_->vv()},
              {"cflx", surface_->cflX()},
              {"cfly", surface_->cflY()}};
  }
  if (coupler_) {
    // The dry-cell hold rule makes the seepage accumulator prognostic; qss
    // only feeds the "seepage" output but restarting it keeps restarted
    // outputs identical (plan §10 P3 restart determinism).
    fields.push_back({"seep_accum", coupler_->seepAccum()});
    fields.push_back({"qss", coupler_->seepageRate()});
  }
  if (rates_ && !rates_->empty()) {
    // Stage-3 T4 (W4.1): the interval-rate accumulators restart bitwise
    // (acc_<name>; the interval start rides as the scalar acc_t_prev).
    // Runs without rate outputs keep the v2.0 checkpoint layout.
    for (const auto& named : rates_->checkpointFields()) {
      fields.push_back(named);
    }
  }
  if (surface_ && surface_->etActive()) {
    // Stage-3 T5 (W5.10): the complementarity clocks are prognostic.
    fields.push_back({"et_wet_clock", surface_->wetClock()});
  }
  if (surface_ && surface_->pairedStore()) {
    // Stage-3 T9 (W9.2): the parked water is prognostic (released later).
    // Runs without the paired store keep the T8 checkpoint layout.
    fields.push_back({"swe_water_store", surface_->waterStore()});
  }
  if (gw_ && gw_->etActive()) {
    fields.push_back({"et_wet_clock_last", gw_->etWetClockLast()});
    fields.push_back({"et_soil_time", gw_->etSoilTime()});
  }
  if (transport_) {
    if (surface_) {
      // The scalar plus the transport state a restart cannot reconstruct:
      // the pre-correction flow-rate snapshots of volume_by_flux (see
      // ScalarSolver.hpp) and — with a subsurface — the carried top-cell
      // dispersion coefficient of the surface exchange term.
      fields.push_back({"s_surf", transport_->surfaceScalar()});
      fields.push_back({"s_fu_old", transport_->flowRateSnapshotX()});
      fields.push_back({"s_fv_old", transport_->flowRateSnapshotY()});
      if (transport_->hasSurfaceStore()) {
        // Stage-3 T5 (W5.7): the store is prognostic (O4).
        fields.push_back({"s_store", transport_->surfaceStore()});
      }
    }
    if (gw_) {
      fields.push_back({"s_dzz_top", transport_->dispersionTopSnapshot()});
    }
  }
  if (temperature_) {
    // The temperature instance's parallel restart state (v2 Q5): the
    // "t_*" layout mirrors "s_*", so a salinity-only checkpoint stays
    // byte-identical to v1 and a two-scalar checkpoint is the union.
    if (surface_) {
      fields.push_back({"t_surf", temperature_->surfaceScalar()});
      fields.push_back({"t_fu_old", temperature_->flowRateSnapshotX()});
      fields.push_back({"t_fv_old", temperature_->flowRateSnapshotY()});
      if (temperature_->hasSurfaceStore()) {
        fields.push_back({"t_store", temperature_->surfaceStore()});  // T9 (R9)
      }
    }
    if (gw_) {
      fields.push_back({"t_dzz_top", temperature_->dispersionTopSnapshot()});
    }
  }
  return fields;
}

io::Checkpoint::Fields3 Simulation::checkpointFields3() const {
  io::Checkpoint::Fields3 fields;
  if (gw_) {
    fields = {{"h", gw_->head()}, {"wc", gw_->waterContent()}};
    if (gw_->hasReallocDiagnostics()) {
      // Stage-3 W1.3: the cumulative diagnostics restart bitwise. The flags
      // follow the requested outputs, so both sides of a restart carry the
      // same checkpoint layout; non-diagnostic runs are byte-identical to
      // the v2.0 layout.
      fields.push_back({"realloc_created", gw_->reallocCreated()});
      fields.push_back({"realloc_dropped", gw_->reallocDroppedCell()});
    }
    const auto& gwVars = config_.output.groundwaterVariables;
    if (gw_->hasReallocDiagnostics() ||
        std::find(gwVars.begin(), gwVars.end(), "vloss") != gwVars.end()) {
      fields.push_back({"vloss_cum", gw_->vlossField()});
    }
    if (transport_) {
      fields.push_back({"s_subs", transport_->subsurfaceScalar()});
    }
    if (transport_ || temperature_) {
      // S3-A9 (a pre-existing gap the T6 vegetated restart lane caught):
      // a subcycled window can take ZERO substeps right after a restart,
      // leaving the transport to read the LAST SUBSTEP'S per-step fluxes
      // and dtg (the documented stale-flux legacy quirk) — which must
      // therefore restart with the state. Runs without transport keep
      // the prior layout.
      fields.push_back({"qx_last", gw_->fluxXVolumetric()});
      fields.push_back({"qy_last", gw_->fluxYVolumetric()});
      fields.push_back({"qzf_last", gw_->fluxZFaceVolumetric()});
      // ... and the face conductivities the transport's limiter gates
      // read (same zero-substep staleness class).
      fields.push_back({"kx_last", gw_->faceConductivityX()});
      fields.push_back({"ky_last", gw_->faceConductivityY()});
      fields.push_back({"kzf_last", gw_->faceConductivityZFace()});
    }
    if (rates_ && !rates_->empty()) {
      // Stage-3 T6: the Field3 interval-rate accumulators (acc_root_uptake).
      for (const auto& named : rates_->checkpointFields3()) {
        fields.push_back(named);
      }
    }
    if (temperature_) {
      fields.push_back({"t_subs", temperature_->subsurfaceScalar()});
    }
  }
  return fields;
}

real_t Simulation::restoreFromCheckpoint() {
  const std::string path = config_.resolvePath(config_.restart.file);
  io::Hdf5Output source(grid_, path);
  io::Checkpoint reader(source);
  const io::Checkpoint::Header header =
      reader.read(config_.restart.time, checkpointFields2(), checkpointFields3());
  if (surface_) {
    if (coupler_) {
      // The refresh reconstructs the stage-boundary edge-slot velocities
      // from the completed step's flow rates, eta^n, and dt; sync-coupled
      // runs marched on the adaptive step, so the configured time.dt the
      // solver was built with is not that dt.
      const auto dtSurf = header.scalars.find("dt_surface");
      if (dtSurf == header.scalars.end()) {
        log::fatal("checkpoint is missing the coupled-step scalar 'dt_surface'");
      }
      surface_->setTimeStep(dtSurf->second);
    }
    surface_->refreshDerivedState();
  }
  if (transport_) {
    // Before the groundwater refresh: the baroclinic ratios and the coupled
    // hydrostatic ghosts read the restored scalar's ghosts.
    transport_->refreshDerivedState(header.t);
  }
  if (temperature_) {
    temperature_->refreshDerivedState(header.t);
  }
  if (gw_) {
    gw_->refreshDerivedState();
    const auto dtg = header.scalars.find("dtg");
    if (dtg == header.scalars.end()) {
      log::fatal("checkpoint is missing the adaptive-step scalar 'dtg'");
    }
    dtg_ = dtg->second;
  }
  if (coupler_) {
    const auto lag = header.scalars.find("tgw_lag");
    if (lag == header.scalars.end()) {
      log::fatal("checkpoint is missing the coupled-clock scalar 'tgw_lag'");
    }
    coupler_->restore(header.t, dtg_, lag->second);
  }
  if (gw_ && (transport_ || temperature_)) {
    const auto dtgLast = header.scalars.find("dtg_last");
    if (dtgLast == header.scalars.end()) {
      log::fatal("checkpoint is missing the transport scalar 'dtg_last' "
                 "(the source run predates the S3-A9 restart fix)");
    }
    gw_->restoreLastDtg(dtgLast->second);
  }
  if (rates_ && !rates_->empty()) {
    const auto tPrev = header.scalars.find("acc_t_prev");
    if (tPrev == header.scalars.end()) {
      log::fatal("checkpoint is missing the interval-rate scalar 'acc_t_prev' "
                 "(the source run did not write the requested rate outputs)");
    }
    rates_->restore(tPrev->second);
  }
  log::info(log::msg() << "restart: resumed from " << path << " at t = " << header.t
                       << " s (step " << header.step << ")");
  return header.t;
}

void Simulation::addFilmRateContribution() {
  if (rates_ && rates_->has("film_removal_rate")) {
    rates_->addInto("film_removal_rate", surface_->filmRemovedDepth());
  }
}

void Simulation::addEtRateContributions() {
  // Stage-3 T5 (W5.10): feed the registered atmosphere rates once per
  // outer step — the surface channel's removal and the soil channel's
  // depth windows (accumulated over this step's substeps, then reset).
  if (!rates_ || rates_->empty()) {
    if (gw_ && gw_->etActive()) {
      gw_->resetEtDepths();
    }
    return;
  }
  const auto& av = config_.output.atmosphereVariables;
  const auto wants = [&av](const char* name) {
    return std::find(av.begin(), av.end(), name) != av.end();
  };
  if (surface_ && surface_->etActive() && wants("evaporation_rate")) {
    rates_->addInto("evaporation_rate", surface_->evapRemoved());
  }
  if (gw_ && gw_->etActive()) {
    if (wants("soil_evaporation_rate")) {
      rates_->addInto("soil_evaporation_rate", gw_->etAppliedDepth());
    }
    if (wants("soil_evaporation_unrealised_rate")) {
      rates_->addInto("soil_evaporation_unrealised_rate", gw_->etUnrealisedDepth());
    }
    gw_->resetEtDepths();
  }
  if (gw_ && gw_->vegActive()) {
    if (wants("transpiration_rate")) {
      rates_->addInto("transpiration_rate", gw_->trAppliedDepth());
    }
    if (wants("transpiration_unrealised_rate")) {
      rates_->addInto("transpiration_unrealised_rate", gw_->trUnrealisedDepth());
    }
    const auto& gv = config_.output.groundwaterVariables;
    if (std::find(gv.begin(), gv.end(), "root_uptake") != gv.end()) {
      rates_->addInto3("root_uptake", gw_->sinkStepWindow());
    }
    gw_->resetTrDepths();
  }
}

void Simulation::writeOutputs(real_t label, real_t tActual) {
  Timer::Scoped timer("io/output");
  const real_t t = label;
  for (const std::string& var : config_.output.surfaceVariables) {
    if (!surface_) {
      break;
    }
    const Field2<real_t>* field = nullptr;
    if (var == "eta") {
      // eta leaves the internal offset frame on output (legacy write_output
      // subtracts the offset for surf files, utility.c:220-222).
      field = &surface_->etaAbsolute();
    } else if (var == "depth") {
      field = &surface_->depth();
    } else if (var == "uu") {
      field = &surface_->uu();
    } else if (var == "vv") {
      field = &surface_->vv();
    } else if (var == "seepage" && coupler_) {
      field = &coupler_->seepageRate();
    } else if (var == "water_store" && surface_->pairedStore()) {
      field = &surface_->waterStore();
    }
    // "seepage_rate" and "film_removal_rate" are written by the
    // IntervalRates registry below.
    if (field != nullptr) {
      output_->writeField2("surface", var, t, *field, surfaceVarMeta(var));
    }
  }
  for (const std::string& var : config_.output.groundwaterVariables) {
    if (!gw_) {
      break;
    }
    const Field3<real_t>* field = nullptr;
    if (var == "hydraulic_head") {
      field = &gw_->head();
    } else if (var == "water_content") {
      field = &gw_->waterContent();
    } else if (var == "qx") {
      field = &gw_->fluxXPerArea();
    } else if (var == "qy") {
      field = &gw_->fluxYPerArea();
    } else if (var == "qz") {
      field = &gw_->fluxZPerArea();
    } else if (var == "realloc_created" && gw_->hasReallocDiagnostics()) {
      field = &gw_->reallocCreated();
    } else if (var == "realloc_dropped" && gw_->hasReallocDiagnostics()) {
      field = &gw_->reallocDroppedCell();
    } else if (var == "vloss") {
      field = &gw_->vlossField();
    }
    if (field != nullptr) {
      output_->writeField3("groundwater", var, t, *field, groundwaterVarMeta(var));
    }
  }
  for (const std::string& var : config_.output.transportVariables) {
    if (!transport_) {
      break;
    }
    if (var == "concentration" && gw_) {
      output_->writeField3("transport", var, t, transport_->subsurfaceScalar(),
                           transportVarMeta(var));
    } else if (var == "concentration_surface" && surface_) {
      output_->writeField2("transport", var, t, transport_->surfaceScalar(),
                           transportVarMeta(var));
    } else if (var == "salt_store" && surface_ && transport_->hasSurfaceStore()) {
      output_->writeField2("transport", var, t, transport_->surfaceStore(),
                           {"psu m3", "per-cell surface salt store (O4)"});
    }
  }
  for (const std::string& var : config_.output.temperatureVariables) {
    if (!temperature_) {
      break;
    }
    if (var == "temperature" && gw_) {
      output_->writeField3("temperature", var, t, temperature_->subsurfaceScalar(),
                           temperatureVarMeta(var));
    } else if (var == "temperature_surface" && surface_) {
      output_->writeField2("temperature", var, t, temperature_->surfaceScalar(),
                           temperatureVarMeta(var));
    }
  }
  for (const std::string& var : config_.output.atmosphereVariables) {
    // Stage-3 T5 (W5.10): the cumulative time clocks are plain field
    // outputs (the *_rate variables are the registry's below).
    const Field2<real_t>* field = nullptr;
    if (var == "wet_time" && surface_ && surface_->etActive()) {
      field = &surface_->wetClock();
    } else if (var == "soil_evaporation_time" && gw_ && gw_->etActive()) {
      field = &gw_->etSoilTime();
    }
    if (field != nullptr) {
      output_->writeField2("atmosphere", var, t, *field,
                           {"s", var == "wet_time"
                                     ? "cumulative open-water-covered time"
                                     : "cumulative soil-evaporation time"});
    }
  }
  if (rates_ && !rates_->empty()) {
    // Datasets keyed by the output label; the mean normalizes over the
    // TRUE elapsed state time so the rate integrates back to the audited
    // volume even when an adaptive step crosses the label (plan §8 W4.1).
    rates_->writeAndReset(*output_, label, tActual < 0.0 ? label : tActual);
  }
}

void Simulation::recordMonitors(real_t t) {
  for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
    if (!monitor->ownsPoint()) {
      monitor->record(t, {});
      continue;
    }
    const MonitorConfig& mc = monitor->config();
    const std::size_t j = static_cast<std::size_t>(mc.j - grid_.j0()) + 1;
    const std::size_t i = static_cast<std::size_t>(mc.i - grid_.i0()) + 1;
    std::vector<real_t> values;
    values.reserve(mc.variables.size());
    for (const std::string& var : mc.variables) {
      const Field2<real_t>* field = nullptr;
      real_t shift = 0.0;
      if (surface_) {
        if (var == "eta") {
          field = &surface_->eta();
          shift = surface_->elevationOffset();
        } else if (var == "depth") {
          field = &surface_->depth();
        } else if (var == "uu") {
          field = &surface_->uu();
        } else if (var == "vv") {
          field = &surface_->vv();
        }
      }
      if (temperature_ && var == "temperature_surface") {
        field = &temperature_->surfaceScalar();
      }
      if (transport_ && var == "concentration_surface") {
        field = &transport_->surfaceScalar();
      }
      if (field == nullptr) {
        log::fatal(log::msg() << "monitor '" << mc.name << "': variable '" << var
                              << "' is not available in this run's module set");
      }
      auto cell = Kokkos::subview(*field, j, i);
      auto host = Kokkos::create_mirror_view(cell);
      Kokkos::deep_copy(host, cell);
      values.push_back(host() - shift);
    }
    monitor->record(t, values);
  }
}

void Simulation::recordMassAudit(real_t t, bool writeRow) {
  const swe::SurfaceStepAudit& audit = surface_->audit();
  const real_t seepage = coupler_ ? coupler_->audit().surfaceGain : 0.0;
  // Stage-3 T9 (W9.1): the counted removal — P1 + P2 always; P3's
  // positive part and the sweep are removals only under the paired store
  // (under discard P3 stays inside clamped, bitwise).
  const bool paired = surface_->pairedStore();
  const real_t film = audit.filmIsolated + audit.filmSmall +
                      (paired ? audit.filmClamp + audit.filmSweep : 0.0);
  const real_t local[13] = {surface_->ownedVolume(), audit.rainVolume, audit.evapVolume,
                            audit.boundaryOutflow, audit.bcInflow, seepage,
                            audit.clampVolume, audit.stageDeposit, film,
                            surface_->pendingFilmVolume(),
                            paired ? surface_->ownedStoreVolume() : 0.0,
                            audit.storeRelease,
                            surface_->stageVolume()};
  real_t global[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  MPI_Reduce(local, global, 13, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
  // Per-path run-record books (W9.1), reduced on every rank's audit.
  {
    const real_t lp[4] = {audit.filmIsolated, audit.filmSmall, audit.filmClamp,
                          audit.filmSweep};
    real_t gp[4] = {0, 0, 0, 0};
    MPI_Reduce(lp, gp, 4, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
    if (grid_.rank() == 0) {
      cumFilmIsolated_ += gp[0];
      cumFilmSmall_ += gp[1];
      cumFilmClamp_ += gp[2];
      cumFilmSweep_ += gp[3];
    }
  }
  if (grid_.rank() == 0) {
    cumRain_ += global[1];
    cumEvap_ += global[2];
    cumOutflow_ += global[3];
    cumBcInflow_ += global[4];
    cumSeepage_ += global[5];
    cumClamped_ += global[6];
    cumStageDeposit_ += global[7];
    cumFilmRemoved_ += global[8];
    cumStoreRelease_ += global[11];
    if (writeRow) {
      std::vector<real_t> row = {global[0], cumRain_, cumEvap_, cumOutflow_, cumBcInflow_};
      if (coupler_) {
        row.push_back(cumSeepage_);
      }
      row.push_back(cumClamped_);
      if (coupler_) {
        row.push_back(cumStageDeposit_);
      }
      row.push_back(cumFilmRemoved_);
      row.push_back(global[9]);
      row.push_back(global[12]);
      if (paired) {
        row.push_back(global[10]);
        row.push_back(cumStoreRelease_);
      }
      massAudit_->record(t, row);
    }
  } else if (writeRow) {
    massAudit_->record(t, {});
  }
}

void Simulation::recordGwMassAudit(real_t t, bool writeRow) {
  // Coupled runs report the window aggregate (all subcycled substeps of the
  // surface step); groundwater-only runs report the single step.
  const gw::GwStepAudit& audit = coupler_ ? coupler_->audit().gw : gw_->audit();
  const real_t local[12] = {gw_->ownedVolume(), audit.boundaryIn, audit.ssStorage,
                            audit.reallocAdjust, audit.reallocDropped, audit.vloss,
                            audit.evap,          audit.etApplied,
                            audit.etDemand - audit.etApplied,
                            audit.trApplied,
                            audit.trDemand - audit.trApplied,
                            audit.cplVentThin};
  real_t global[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  MPI_Reduce(local, global, 12, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
  if (grid_.rank() == 0) {
    cumGwBoundary_ += global[1];
    cumGwStorage_ += global[2];
    cumGwRealloc_ += global[3];
    cumGwDropped_ += global[4];
    cumGwVloss_ += global[5];
    cumGwEvap_ += global[6];
    cumGwEtApplied_ += global[7];
    cumGwEtUnrealised_ += global[8];
    cumGwTrApplied_ += global[9];
    cumGwTrUnrealised_ += global[10];
    cumGwVentThin_ += global[11];
    if (writeRow) {
      std::vector<real_t> row = {global[0], cumGwBoundary_, cumGwStorage_, cumGwRealloc_,
                                 cumGwDropped_, cumGwVloss_};
      if (config_.groundwater.evaporation.enabled) {
        row.push_back(cumGwEvap_);
      }
      if (config_.evapotranspiration.enabled && config_.evapotranspiration.soil) {
        row.push_back(cumGwEtApplied_);
        row.push_back(cumGwEtUnrealised_);
      }
      if (config_.vegetation.enabled) {
        row.push_back(cumGwTrApplied_);
        row.push_back(cumGwTrUnrealised_);
      }
      gwMassAudit_->record(t, row);
    }
  } else if (writeRow) {
    gwMassAudit_->record(t, {});
  }
}

void Simulation::recordTransportAudit(real_t t, bool writeRow) {
  const transport::TransportAudit& audit = transport_->audit();
  const real_t local[11] = {transport_->ownedSurfaceMass(), transport_->ownedSubsurfaceMass(),
                            audit.exchange,     audit.surfSource, audit.surfBoundary,
                            audit.surfAdjust,   audit.surfAnchor, audit.subsBoundary,
                            audit.subsAdjust,   audit.subsAnchor,
                            transport_->surfaceStoreTotal()};
  real_t global[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  MPI_Reduce(local, global, 11, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
  if (grid_.rank() == 0) {
    cumTrExchange_ += global[2];
    cumTrSurfSource_ += global[3];
    cumTrSurfBoundary_ += global[4];
    cumTrSurfAdjust_ += global[5];
    cumTrSurfAnchor_ += global[6];
    cumTrSubsBoundary_ += global[7];
    cumTrSubsAdjust_ += global[8];
    cumTrSubsAnchor_ += global[9];
    if (writeRow) {
      std::vector<real_t> row = {global[0], cumTrExchange_, cumTrSurfSource_,
                                 cumTrSurfBoundary_, cumTrSurfAdjust_, cumTrSurfAnchor_,
                                 cumTrSubsBoundary_, cumTrSubsAdjust_, cumTrSubsAnchor_};
      row.insert(row.begin() + 1, global[1]);
      if (transport_->hasSurfaceStore()) {
        row.push_back(global[10]);
      }
      transportAudit_->record(t, row);
    }
  } else if (writeRow) {
    transportAudit_->record(t, {});
  }
}

void Simulation::recordTemperatureAudit(real_t t, bool writeRow) {
  const transport::TransportAudit& audit = temperature_->audit();
  const real_t local[12] = {temperature_->ownedSurfaceMass(),
                            temperature_->ownedSubsurfaceMass(),
                            audit.exchange,     audit.surfSource, audit.surfBoundary,
                            audit.surfAdjust,   audit.surfAnchor, audit.surfAtmos,
                            audit.subsBoundary, audit.subsAdjust, audit.subsAnchor,
                            temperature_->surfaceStoreTotal()};
  real_t global[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  MPI_Reduce(local, global, 12, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
  if (grid_.rank() == 0) {
    cumTpExchange_ += global[2];
    cumTpSurfSource_ += global[3];
    cumTpSurfBoundary_ += global[4];
    cumTpSurfAdjust_ += global[5];
    cumTpSurfAnchor_ += global[6];
    cumTpSurfAtmos_ += global[7];
    cumTpSubsBoundary_ += global[8];
    cumTpSubsAdjust_ += global[9];
    cumTpSubsAnchor_ += global[10];
    if (writeRow) {
      std::vector<real_t> row = {global[0], global[1], cumTpExchange_,
                                 cumTpSurfSource_, cumTpSurfBoundary_,
                                 cumTpSurfAdjust_, cumTpSurfAnchor_, cumTpSurfAtmos_,
                                 cumTpSubsBoundary_, cumTpSubsAdjust_,
                                 cumTpSubsAnchor_};
      if (temperature_->hasSurfaceStore()) {
        row.push_back(global[11]);  // T9 (R9): the parked heat [K m^3]
      }
      temperatureAudit_->record(t, row);
    }
  } else if (writeRow) {
    temperatureAudit_->record(t, {});
  }
}

bool Simulation::monitorRowDue(real_t t, bool forced) {
  const real_t mi = config_.output.monitorInterval;
  if (mi <= 0.0) {
    return true;
  }
  if (t >= nextMonitorRow_ - 1.0e-9) {
    // Advance past every boundary this step crossed (the output-schedule
    // rule; W4.3 restart recomputes nextMonitorRow_ from t0 the same way).
    nextMonitorRow_ = (std::floor(t / mi + 1.0e-9) + 1.0) * mi;
    return true;
  }
  return forced;
}

void Simulation::writeCheckpoint(real_t t, long step, real_t labelTime) {
  Timer::Scoped timer("io/checkpoint");
  io::Checkpoint::Scalars scalars;
  if (gw_) {
    scalars["dtg"] = coupler_ ? coupler_->currentDtg() : dtg_;
  }
  if (coupler_) {
    scalars["tgw_lag"] = coupler_->gwLag(t);
    scalars["dt_surface"] = coupler_->lastSurfaceDt();
  }
  if (rates_ && !rates_->empty()) {
    scalars["acc_t_prev"] = rates_->tPrev();
  }
  if (gw_ && (transport_ || temperature_)) {
    scalars["dtg_last"] = gw_->lastDtg();  // S3-A9
  }
  checkpoint_->write(t, step, scalars, checkpointFields2(), checkpointFields3(), labelTime);
}

void Simulation::run() {
  Timer::Scoped total("simulation");
  real_t t0 = config_.time.tStart;
  {
    // Everything before the time loop proper: restart restore or the
    // initial output/monitor writes (v2 plan §2A: the 'init' section; the
    // module/grid construction time is outside "simulation" and appears in
    // provenance.wall_seconds minus the loop).
    Timer::Scoped init("init");
    if (config_.restart.enabled) {
      t0 = restoreFromCheckpoint();
    } else {
      writeOutputs(t0);
      recordMonitors(t0);
      if (surface_) {
        recordMassAudit(t0);
      }
      if (gw_) {
        recordGwMassAudit(t0);
      }
      if (transport_) {
        recordTransportAudit(t0);
      }
      if (temperature_) {
        recordTemperatureAudit(t0);
      }
    }
  }
  // W4.3: the monitor row schedule starts at t0 exactly like the output
  // schedule (restart recomputes it from the restored t0).
  if (config_.output.monitorInterval > 0.0) {
    nextMonitorRow_ = (std::floor(t0 / config_.output.monitorInterval + 1.0e-9) + 1.0) *
                      config_.output.monitorInterval;
  }
  if (coupler_) {
    runCoupledLoop(t0);
  } else if (surface_) {
    runSurfaceLoop(t0);
  } else {
    runGroundwaterLoop(t0);
  }

  // End-of-run per-system solver telemetry (v2 plan §2.2.6): one line per
  // system on rank 0, parsed by the g2/s-gate harness — and the same values
  // feed the run record's solver section (v2 plan §2A), so the two cannot
  // disagree.
  std::ostringstream solverYaml;
  const auto logSolverSummary = [this, &solverYaml](const char* name,
                                                    const SolverTelemetry& t,
                                                    const std::string& preconditioner,
                                                    const std::string& matType,
                                                    const std::string& amgCoarsen,
                                                    const std::string& amgRelax) {
    // Iteration counts are KSP-collective and identical on every rank; the
    // setup/solve times are reduced to their max so the summary reports the
    // critical path.
    double local[2] = {t.setupSeconds, t.solveSeconds};
    double reduced[2] = {0.0, 0.0};
    MPI_Allreduce(local, reduced, 2, MPI_DOUBLE, MPI_MAX, grid_.comm());
    if (t.solves == 0) {
      return;
    }
    const double mean =
        static_cast<double>(t.totalIterations) / static_cast<double>(t.solves);
    log::info(log::msg() << "solver summary " << name << ": solves=" << t.solves
                         << " iters_mean=" << mean << " iters_max=" << t.maxIterations
                         << " rebuilds=" << t.rebuilds << " retries=" << t.retries
                         << " setup_s=" << reduced[0] << " solve_s=" << reduced[1]);
    // preconditioner/mat_type record the *resolved* solver backend (a
    // -<prefix>mat_type override is read back from PETSc), so a run is
    // self-describing about which linear algebra it actually got (v2 plan
    // §2B.5; the p2/p3 gates key on mat_type).
    solverYaml << name << ":\n"
               << "  preconditioner: " << preconditioner << "\n"
               << "  mat_type: " << matType << "\n";
    if (!amgCoarsen.empty()) {
      solverYaml << "  amg_coarsen_type: " << amgCoarsen << "\n"
                 << "  amg_relax_type: " << amgRelax << "\n";
    }
    solverYaml << "  solves: " << t.solves << "\n"
               << "  iters_mean: " << mean << "\n"
               << "  iters_max: " << t.maxIterations << "\n"
               << "  rebuilds: " << t.rebuilds << "\n"
               << "  retries: " << t.retries << "\n"
               << "  setup_s: " << reduced[0] << "\n"
               << "  solve_s: " << reduced[1] << "\n";
  };
  if (surface_) {
    logSolverSummary("fs", surface_->solverTelemetry(),
                     config_.solver.surface.preconditioner, surface_->solverMatType(),
                     surface_->solverAmgCoarsenType(), surface_->solverAmgRelaxType());
  }
  if (gw_) {
    logSolverSummary("gw", gw_->solverTelemetry(),
                     config_.solver.groundwater.preconditioner, gw_->solverMatType(),
                     gw_->solverAmgCoarsenType(), gw_->solverAmgRelaxType());
  }
  runRecord_->setSolver(solverYaml.str());
}

void Simulation::finalizeRunRecord() { flushRunRecord(true); }

void Simulation::flushRunRecord(bool finished) {
  // Closure: the final cumulative budgets already reduced to rank 0 by the
  // audit recorders (the same numbers the mass-audit monitor tables carry).
  std::ostringstream closure;
  closure.setf(std::ios::scientific);
  closure.precision(10);
  if (surface_) {
    closure << "surface:\n"
            << "  rain_m3: " << cumRain_ << "\n"
            << "  evaporation_m3: " << cumEvap_ << "\n"
            << "  boundary_outflow_m3: " << cumOutflow_ << "\n"
            << "  bc_inflow_m3: " << cumBcInflow_ << "\n"
            << "  seepage_m3: " << cumSeepage_ << "\n"
            << "  clamped_m3: " << cumClamped_ << "\n";
    if (coupler_) {
      closure << "  stage_deposit_m3: " << cumStageDeposit_ << "\n";
    }
    // Stage-3 T9 (W9.1): the thin-layer removals per path (P1 isolated,
    // P2 small, P3 clamp — the clamp's positive part, informational under
    // discard where it stays inside clamped — and the W9.2 sweep), the
    // counted total the film_removed column carries, and under the paired
    // store the parked volume still held and the cumulative release.
    closure << "  film_removed_m3: " << cumFilmRemoved_ << "\n"
            << "  film_removed_isolated_m3: " << cumFilmIsolated_ << "\n"
            << "  film_removed_small_m3: " << cumFilmSmall_ << "\n"
            << "  film_removed_clamp_m3: " << cumFilmClamp_ << "\n"
            << "  film_removed_sweep_m3: " << cumFilmSweep_ << "\n";
    {
      const real_t local[2] = {surface_->ownedStoreVolume(), surface_->pendingFilmVolume()};
      real_t global[2] = {0, 0};
      MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, grid_.comm());
      closure << "  film_pending_m3: " << global[1] << "\n";
      if (surface_->pairedStore()) {
        closure << "  store_volume_m3: " << global[0] << "\n"
                << "  store_release_m3: " << cumStoreRelease_ << "\n";
      }
    }
    // Stage-3 T2 (W2.6): the stage-boundary free-energy observables. The
    // KE influx is the report_V1 engine metric (zero at a quiet boundary);
    // the ratio is the entrance head actually paid per v^2/2g (~1 under
    // closure external, ~0 under zero_gradient with fast inflow). The
    // zero_gradient warning on these is evaluated here, run end.
    {
      real_t local[3] = {surface_->stageKeInflux(), surface_->stageHeadRatioNum(),
                         surface_->stageHeadRatioDen()};
      real_t global[3] = {0, 0, 0};
      MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_SUM, grid_.comm());
      if (grid_.rank() == 0) {
        closure << "  stage_ke_influx_J: " << global[0] << "\n";
        if (global[2] > 0.0) {
          closure << "  stage_entrance_head_ratio: " << global[1] / global[2] << "\n";
        }
      }
      if (finished && grid_.rank() == 0 && !surface_->hasExternalClosure() &&
          global[2] > 0.0 && global[1] / global[2] < 0.5 && global[0] > 1.0) {
        const std::string warning =
            "eta boundaries: inflow imported " + std::to_string(global[0]) +
            " J of kinetic energy while paying " +
            std::to_string(global[1] / global[2]) +
            " of the v^2/2g entrance head; if the water beyond the boundary "
            "is still, use inflow: {closure: external} (Stage-3 T2)";
        log::warn(warning);
        runRecord_->addWarning(warning);
      }
    }
  }
  if (gw_) {
    closure << "groundwater:\n"
            << "  reallocation: "
            << (config_.groundwater.reallocation ==
                        GroundwaterConfig::Reallocation::Off
                    ? "disabled"
                    : "legacy")
            << "\n"
            << "  boundary_in_m3: " << cumGwBoundary_ << "\n"
            << "  ss_storage_m3: " << cumGwStorage_ << "\n"
            << "  realloc_m3: " << cumGwRealloc_ << "\n"
            << "  realloc_dropped_m3: " << cumGwDropped_ << "\n"
            << "  vloss_m3: " << cumGwVloss_ << "\n";
    if (surface_ && surface_->pairedStore() && coupler_) {
      // Stage-3 T9 (W9.4): vents the legacy rule would have discarded,
      // deposited for the sweep to park with their salt (P4).
      closure << "  vent_thin_m3: " << cumGwVentThin_ << "\n";
    }
    if (config_.groundwater.predictorCorrector ==
        GroundwaterConfig::PredictorCorrector::Consistent) {
      // Stage-3 T11 (P4): the mode-3 guard clamp's cumulative volume
      // (expected 0; a nonzero value is the C-SEEP chatter window).
      real_t clampLocal = gw_->pcClampVolume();
      real_t clamp = 0.0;
      MPI_Allreduce(&clampLocal, &clamp, 1, MPI_DOUBLE, MPI_SUM, grid_.comm());
      closure << "  pc_clamp_m3: " << clamp << "\n";
    }
  }
  if (config_.evapotranspiration.enabled) {
    // Stage-3 T5 (W5.10): the ET closure lines. open_water is the
    // mass_audit evaporation column under ET; the soil pair are the gw
    // audit's T5 columns; salt_store_final is the surface store's
    // run-end content (0 without the store).
    closure << "evapotranspiration:\n";
    if (surface_ && surface_->etActive()) {
      closure << "  open_water_evaporation_m3: " << cumEvap_ << "\n";
    }
    if (gw_ && gw_->etActive()) {
      closure << "  soil_evaporation_m3: " << cumGwEtApplied_ << "\n"
              << "  soil_evaporation_unrealised_m3: " << cumGwEtUnrealised_ << "\n";
    }
    if (gw_ && gw_->vegActive()) {
      closure << "  transpiration_m3: " << cumGwTrApplied_ << "\n"
              << "  transpiration_unrealised_m3: " << cumGwTrUnrealised_ << "\n";
    }
    if (transport_ && surface_) {
      real_t localStore = transport_->surfaceStoreTotal();
      real_t globalStore = 0.0;
      MPI_Allreduce(&localStore, &globalStore, 1, MPI_DOUBLE, MPI_SUM, grid_.comm());
      closure << "  salt_store_final: " << globalStore << "\n";
    }
  }
  runRecord_->setClosure(closure.str());
  runRecord_->flush(finished);
}

void Simulation::runCoupledLoop(real_t t0) {
  const bool sync = (config_.coupling.mode == CouplingConfig::Mode::Sync);
  const real_t tEnd = config_.time.tEnd;
  const real_t interval = config_.time.outputInterval;
  const real_t checkpointInterval = config_.output.checkpointInterval;

  // Sync-coupled runs march on the common adaptive step (amendment A10), so
  // the loop uses the crossed-boundary output/checkpoint semantics of the
  // adaptive groundwater loop (amendment A6); with the fixed subcycled dt
  // the crossed boundaries coincide with the aligned steps.
  real_t nextOutput = (std::floor(t0 / interval + 1.0e-9) + 1.0) * interval;
  real_t nextCheckpoint = (checkpointInterval > 0.0)
                              ? (std::floor(t0 / checkpointInterval + 1.0e-9) + 1.0) *
                                    checkpointInterval
                              : 0.0;
  real_t lastCheckpointTime = -1.0;

  log::info(log::msg() << "time loop: coupled (" << (sync ? "sync" : "subcycled")
                       << "), dt = " << coupler_->nextDt() << " s to t = " << tEnd << " s");

  real_t t = t0;
  long step = 0;
  while (t < tEnd - 1.0e-9) {
    const real_t dt = coupler_->nextDt();
    t += dt;
    ++step;
    coupler_->step(t, dt);
    if (rates_ && rates_->has("seepage_rate")) {
      // W4.2: the step's complete coupled surface gain (applied seepage +
      // bounce-back + vent deposits) — the same numbers recordMassAudit
      // reduces into the cumulative seepage column (the dg1 identity).
      rates_->addInto("seepage_rate", coupler_->gain());
    }
    addEtRateContributions();
    addFilmRateContribution();
    stepTransport(t, dt, gw_->lastDtg());

    const bool checkpointDue =
        checkpointInterval > 0.0 && t >= nextCheckpoint - 1.0e-9;
    const bool rowDue = monitorRowDue(t, checkpointDue || t >= tEnd - 1.0e-9);
    {
      Timer::Scoped monitorsTimer("monitors");
      if (rowDue) {
        recordMonitors(t);
      }
      recordMassAudit(t, rowDue);
      recordGwMassAudit(t, rowDue);
      if (transport_) {
        recordTransportAudit(t, rowDue);
      }
      if (temperature_) {
        recordTemperatureAudit(t, rowDue);
      }
    }

    real_t cflLocal = surface_->maxCfl();
    real_t cfl = 0.0;
    MPI_Allreduce(&cflLocal, &cfl, 1, MPI_DOUBLE, MPI_MAX, grid_.comm());
    if (cfl > 1.0 && grid_.rank() == 0) {
      log::warn(log::msg() << "CFL = " << cfl << " exceeds 1 at t = " << t << " s");
    }

    if (t >= nextOutput - 1.0e-9) {
      const real_t label = nextOutput;
      // Advance past every boundary this step crossed (see the
      // groundwater loop's note).
      nextOutput = (std::floor(t / interval + 1.0e-9) + 1.0) * interval;
      writeOutputs(label, t);
      for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
        monitor->flush();
      }
      massAudit_->flush();
      gwMassAudit_->flush();
      if (transportAudit_) {
        transportAudit_->flush();
      }
      if (temperatureAudit_) {
        temperatureAudit_->flush();
      }
      output_->flush();
      flushRunRecord(false);
      log::info(log::msg() << "output written at t = " << t << " s, dt = " << dt
                           << " s (fs " << surface_->lastSolve().iterations << " it, gw "
                           << gw_->lastSolve().iterations << " it)");
    }
    if (checkpointDue) {
      const real_t label = nextCheckpoint;
      nextCheckpoint =
          (std::floor(t / checkpointInterval + 1.0e-9) + 1.0) * checkpointInterval;
      lastCheckpointTime = t;
      writeCheckpoint(t, step, label);
    }
  }

  if (checkpointInterval > 0.0 && lastCheckpointTime != t) {
    writeCheckpoint(t, step, std::floor(t + 0.5));
  }
  for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
    monitor->flush();
  }
  massAudit_->flush();
  gwMassAudit_->flush();
  if (transportAudit_) {
    transportAudit_->flush();
  }
  if (temperatureAudit_) {
    temperatureAudit_->flush();
  }
  output_->flush();
  flushRunRecord(false);
}

void Simulation::runSurfaceLoop(real_t t0) {
  const real_t dt = config_.time.dt;
  const real_t tEnd = config_.time.tEnd;
  const real_t interval = config_.time.outputInterval;
  const real_t checkpointInterval = config_.output.checkpointInterval;
  long lastOutputIndex = static_cast<long>(std::llround(t0 / interval));
  long lastCheckpointIndex =
      (checkpointInterval > 0.0) ? static_cast<long>(std::llround(t0 / checkpointInterval)) : 0;
  real_t lastCheckpointTime = -1.0;

  const long nSteps = static_cast<long>(std::ceil((tEnd - t0) / dt - 1.0e-9));
  log::info(log::msg() << "time loop: " << nSteps << " steps of dt = " << dt << " s to t = "
                       << tEnd << " s");

  for (long step = 1; step <= nSteps; ++step) {
    const real_t t = t0 + static_cast<real_t>(step) * dt;
    surface_->beginStep(t);
    surface_->solveFreeSurface();
    // Stage-3 T9 (W9.2): the paired store's end-of-step sweep and release
    // after the last eta change of the step (no-op unless the store is on).
    surface_->sweepAndRelease();
    surface_->updateVelocity();
    addEtRateContributions();
    addFilmRateContribution();
    stepTransport(t, dt, 0.0);

    bool checkpointDue = false;
    if (checkpointInterval > 0.0) {
      const long ci = static_cast<long>(std::llround(t / checkpointInterval));
      checkpointDue = ci > lastCheckpointIndex &&
                      std::fabs(t - static_cast<real_t>(ci) * checkpointInterval) <
                          0.5 * dt;
    }
    const bool rowDue = monitorRowDue(t, checkpointDue || step == nSteps);
    {
      Timer::Scoped monitorsTimer("monitors");
      if (rowDue) {
        recordMonitors(t);
      }
      recordMassAudit(t, rowDue);
      if (transport_) {
        recordTransportAudit(t, rowDue);
      }
      if (temperature_) {
        recordTemperatureAudit(t, rowDue);
      }
    }

    real_t cflLocal = surface_->maxCfl();
    real_t cfl = 0.0;
    MPI_Allreduce(&cflLocal, &cfl, 1, MPI_DOUBLE, MPI_MAX, grid_.comm());
    if (cfl > 1.0 && grid_.rank() == 0) {
      log::warn(log::msg() << "CFL = " << cfl << " exceeds 1 at t = " << t << " s");
    }

    const long outputIndex = static_cast<long>(std::llround(t / interval));
    if (outputIndex > lastOutputIndex &&
        std::fabs(t - static_cast<real_t>(outputIndex) * interval) < 0.5 * dt) {
      lastOutputIndex = outputIndex;
      writeOutputs(static_cast<real_t>(outputIndex) * interval);
      for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
        monitor->flush();
      }
      massAudit_->flush();
      if (transportAudit_) {
        transportAudit_->flush();
      }
      if (temperatureAudit_) {
        temperatureAudit_->flush();
      }
      output_->flush();
      flushRunRecord(false);
      log::info(log::msg() << "output written at t = " << t << " s (solver "
                           << surface_->lastSolve().iterations << " it, residual "
                           << surface_->lastSolve().residualNorm << ")");
    }
    if (checkpointInterval > 0.0) {
      const long checkpointIndex = static_cast<long>(std::llround(t / checkpointInterval));
      if (checkpointIndex > lastCheckpointIndex &&
          std::fabs(t - static_cast<real_t>(checkpointIndex) * checkpointInterval) < 0.5 * dt) {
        lastCheckpointIndex = checkpointIndex;
        lastCheckpointTime = static_cast<real_t>(checkpointIndex) * checkpointInterval;
        writeCheckpoint(lastCheckpointTime, step);
      }
    }
  }

  // A final checkpoint is always written when checkpointing is on (plan §6),
  // unless the interval already produced one at t_end.
  if (checkpointInterval > 0.0 && lastCheckpointTime != tEnd) {
    writeCheckpoint(tEnd, nSteps);
  }
  for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
    monitor->flush();
  }
  massAudit_->flush();
  if (transportAudit_) {
    transportAudit_->flush();
  }
  if (temperatureAudit_) {
    temperatureAudit_->flush();
  }
  output_->flush();
  flushRunRecord(false);
}

void Simulation::runGroundwaterLoop(real_t t0) {
  const real_t tEnd = config_.time.tEnd;
  const real_t interval = config_.time.outputInterval;
  const real_t checkpointInterval = config_.output.checkpointInterval;

  // Groundwater-only runs march on the adaptive dtg (legacy feeds the
  // adapted dtg back into the outer step, solve_groundwater:193); outputs
  // are labeled with the crossed multiple of output_interval, holding the
  // state within one dtg of it — exactly the legacy trigger semantics
  // (solve.c:121-126) expressed against boundaries instead of rounding.
  real_t nextOutput =
      (std::floor(t0 / interval + 1.0e-9) + 1.0) * interval;
  real_t nextCheckpoint = (checkpointInterval > 0.0)
                              ? (std::floor(t0 / checkpointInterval + 1.0e-9) + 1.0) *
                                    checkpointInterval
                              : 0.0;
  real_t lastCheckpointTime = -1.0;

  log::info(log::msg() << "time loop: adaptive dtg from " << dtg_ << " s to t = " << tEnd
                       << " s");

  real_t t = t0;
  long step = 0;
  while (t < tEnd - 1.0e-9) {
    const real_t dtg = dtg_;
    t += dtg;
    ++step;
    gw_->step(t, dtg);
    dtg_ = gw_->nextDt();
    addEtRateContributions();
    stepTransport(t, dtg, gw_->lastDtg());

    const bool checkpointDue =
        checkpointInterval > 0.0 && t >= nextCheckpoint - 1.0e-9;
    const bool rowDue = monitorRowDue(t, checkpointDue || t >= tEnd - 1.0e-9);
    {
      Timer::Scoped monitorsTimer("monitors");
      if (rowDue) {
        recordMonitors(t);
      }
      recordGwMassAudit(t, rowDue);
      if (transport_) {
        recordTransportAudit(t, rowDue);
      }
      if (temperature_) {
        recordTemperatureAudit(t, rowDue);
      }
    }

    if (t >= nextOutput - 1.0e-9) {
      const real_t label = nextOutput;
      // Advance past every boundary this step crossed (a step larger than
      // the interval crosses several; bumping by one interval per step
      // would let the labels drift behind the state — found at P3).
      nextOutput = (std::floor(t / interval + 1.0e-9) + 1.0) * interval;
      writeOutputs(label);
      for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
        monitor->flush();
      }
      gwMassAudit_->flush();
      if (transportAudit_) {
        transportAudit_->flush();
      }
      if (temperatureAudit_) {
        temperatureAudit_->flush();
      }
      output_->flush();
      flushRunRecord(false);
      log::info(log::msg() << "output written at t = " << t << " s, dtg = " << dtg_
                           << " s (solver " << gw_->lastSolve().iterations << " it, residual "
                           << gw_->lastSolve().residualNorm << ")");
    }
    if (checkpointDue) {
      const real_t label = nextCheckpoint;
      nextCheckpoint =
          (std::floor(t / checkpointInterval + 1.0e-9) + 1.0) * checkpointInterval;
      lastCheckpointTime = t;
      writeCheckpoint(t, step, label);
    }
  }

  if (checkpointInterval > 0.0 && lastCheckpointTime != t) {
    writeCheckpoint(t, step, std::floor(t + 0.5));
  }
  for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
    monitor->flush();
  }
  gwMassAudit_->flush();
  if (transportAudit_) {
    transportAudit_->flush();
  }
  if (temperatureAudit_) {
    temperatureAudit_->flush();
  }
  output_->flush();
  flushRunRecord(false);
}

}  // namespace frehg::driver
