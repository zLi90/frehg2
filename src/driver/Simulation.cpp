/// \file Simulation.cpp
/// \brief Driver implementation: setup, time loops, outputs, restart.

#include "driver/Simulation.hpp"

#include "core/Logger.hpp"
#include "core/Timer.hpp"

#include <petscsys.h>

#include <algorithm>
#include <cmath>
#include <string>

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
    return {"m s-1", "surface-applied seepage rate (positive upward)"};
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

Simulation::Simulation(MPI_Comm comm, const FrehgConfig& config)
    : config_(config), grid_(comm, config.domain) {
  if (!config_.modules.surfaceWater && !config_.modules.groundwater) {
    log::fatal("no module enabled: set modules.surface_water or modules.groundwater");
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
  if (config_.modules.transport) {
    buildTransport();
  }

  output_ = std::make_unique<io::Hdf5Output>(grid_, config_.output.filename, config_.rawText);
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
    massAudit_ = std::make_unique<io::Monitor>(*output_, auditConfig);
  }
  if (config_.modules.groundwater) {
    MonitorConfig auditConfig;
    auditConfig.name = "gw_mass_audit";
    auditConfig.i = 0;
    auditConfig.j = 0;
    auditConfig.variables = {"volume", "boundary_in", "ss_storage", "realloc",
                             "realloc_dropped", "vloss"};
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
    transportAudit_ = std::make_unique<io::Monitor>(*output_, auditConfig);
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
  }
  transport::SubsurfaceWiring subsWiring;
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
  }
  transport::CouplingWiring cplWiring;
  if (coupler_) {
    cplWiring.active = true;
    cplWiring.qss = coupler_->seepageRate();
  }
  transport_ = std::make_unique<transport::ScalarSolver>(grid_, config_, *boundaries_, *halo_,
                                                         surfWiring, subsWiring, cplWiring);
  if (gw_) {
    // Baroclinic activation (plan §10 P4): the density/viscosity ratios read
    // the transport scalar at every subsurface step.
    gw_->attachScalar(transport_->subsurfaceScalar(),
                      surface_ ? transport_->surfaceScalar() : Field2<real_t>(),
                      config_.groundwater.densityCoupling.enabled);
  }
}

void Simulation::stepTransport(real_t t, real_t dt, real_t dtgLast) {
  if (!transport_) {
    return;
  }
  const real_t rain = surface_ ? surface_->currentRain() : 0.0;
  const real_t evap = surface_ ? surface_->currentEvaporation() : 0.0;
  transport_->step(t, dt, dtgLast, rain, evap);
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
  if (transport_) {
    if (surface_) {
      // The scalar plus the transport state a restart cannot reconstruct:
      // the pre-correction flow-rate snapshots of volume_by_flux (see
      // ScalarSolver.hpp) and — with a subsurface — the carried top-cell
      // dispersion coefficient of the surface exchange term.
      fields.push_back({"s_surf", transport_->surfaceScalar()});
      fields.push_back({"s_fu_old", transport_->flowRateSnapshotX()});
      fields.push_back({"s_fv_old", transport_->flowRateSnapshotY()});
    }
    if (gw_) {
      fields.push_back({"s_dzz_top", transport_->dispersionTopSnapshot()});
    }
  }
  return fields;
}

io::Checkpoint::Fields3 Simulation::checkpointFields3() const {
  io::Checkpoint::Fields3 fields;
  if (gw_) {
    fields = {{"h", gw_->head()}, {"wc", gw_->waterContent()}};
    if (transport_) {
      fields.push_back({"s_subs", transport_->subsurfaceScalar()});
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
  log::info(log::msg() << "restart: resumed from " << path << " at t = " << header.t
                       << " s (step " << header.step << ")");
  return header.t;
}

void Simulation::writeOutputs(real_t t) {
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
    }
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
    }
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

void Simulation::recordMassAudit(real_t t) {
  const swe::SurfaceStepAudit& audit = surface_->audit();
  const real_t seepage = coupler_ ? coupler_->audit().surfaceGain : 0.0;
  const real_t local[7] = {surface_->ownedVolume(), audit.rainVolume, audit.evapVolume,
                           audit.boundaryOutflow, audit.bcInflow, seepage,
                           audit.clampVolume};
  real_t global[7] = {0, 0, 0, 0, 0, 0, 0};
  MPI_Reduce(local, global, 7, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
  if (grid_.rank() == 0) {
    cumRain_ += global[1];
    cumEvap_ += global[2];
    cumOutflow_ += global[3];
    cumBcInflow_ += global[4];
    cumSeepage_ += global[5];
    cumClamped_ += global[6];
    std::vector<real_t> row = {global[0], cumRain_, cumEvap_, cumOutflow_, cumBcInflow_};
    if (coupler_) {
      row.push_back(cumSeepage_);
    }
    row.push_back(cumClamped_);
    massAudit_->record(t, row);
  } else {
    massAudit_->record(t, {});
  }
}

void Simulation::recordGwMassAudit(real_t t) {
  // Coupled runs report the window aggregate (all subcycled substeps of the
  // surface step); groundwater-only runs report the single step.
  const gw::GwStepAudit& audit = coupler_ ? coupler_->audit().gw : gw_->audit();
  const real_t local[6] = {gw_->ownedVolume(), audit.boundaryIn, audit.ssStorage,
                           audit.reallocAdjust, audit.reallocDropped, audit.vloss};
  real_t global[6] = {0, 0, 0, 0, 0, 0};
  MPI_Reduce(local, global, 6, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
  if (grid_.rank() == 0) {
    cumGwBoundary_ += global[1];
    cumGwStorage_ += global[2];
    cumGwRealloc_ += global[3];
    cumGwDropped_ += global[4];
    cumGwVloss_ += global[5];
    gwMassAudit_->record(t, {global[0], cumGwBoundary_, cumGwStorage_, cumGwRealloc_,
                             cumGwDropped_, cumGwVloss_});
  } else {
    gwMassAudit_->record(t, {});
  }
}

void Simulation::recordTransportAudit(real_t t) {
  const transport::TransportAudit& audit = transport_->audit();
  const real_t local[10] = {transport_->ownedSurfaceMass(), transport_->ownedSubsurfaceMass(),
                            audit.exchange,     audit.surfSource, audit.surfBoundary,
                            audit.surfAdjust,   audit.surfAnchor, audit.subsBoundary,
                            audit.subsAdjust,   audit.subsAnchor};
  real_t global[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  MPI_Reduce(local, global, 10, MPI_DOUBLE, MPI_SUM, 0, grid_.comm());
  if (grid_.rank() == 0) {
    cumTrExchange_ += global[2];
    cumTrSurfSource_ += global[3];
    cumTrSurfBoundary_ += global[4];
    cumTrSurfAdjust_ += global[5];
    cumTrSurfAnchor_ += global[6];
    cumTrSubsBoundary_ += global[7];
    cumTrSubsAdjust_ += global[8];
    cumTrSubsAnchor_ += global[9];
    transportAudit_->record(t, {global[0], global[1], cumTrExchange_, cumTrSurfSource_,
                                cumTrSurfBoundary_, cumTrSurfAdjust_, cumTrSurfAnchor_,
                                cumTrSubsBoundary_, cumTrSubsAdjust_, cumTrSubsAnchor_});
  } else {
    transportAudit_->record(t, {});
  }
}

void Simulation::writeCheckpoint(real_t t, long step, real_t labelTime) {
  io::Checkpoint::Scalars scalars;
  if (gw_) {
    scalars["dtg"] = coupler_ ? coupler_->currentDtg() : dtg_;
  }
  if (coupler_) {
    scalars["tgw_lag"] = coupler_->gwLag(t);
    scalars["dt_surface"] = coupler_->lastSurfaceDt();
  }
  checkpoint_->write(t, step, scalars, checkpointFields2(), checkpointFields3(), labelTime);
}

void Simulation::run() {
  Timer::Scoped total("simulation");
  real_t t0 = config_.time.tStart;
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
  }
  if (coupler_) {
    runCoupledLoop(t0);
  } else if (surface_) {
    runSurfaceLoop(t0);
  } else {
    runGroundwaterLoop(t0);
  }
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
    stepTransport(t, dt, gw_->lastDtg());

    recordMonitors(t);
    recordMassAudit(t);
    recordGwMassAudit(t);
    if (transport_) {
      recordTransportAudit(t);
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
      writeOutputs(label);
      for (const std::unique_ptr<io::Monitor>& monitor : monitors_) {
        monitor->flush();
      }
      massAudit_->flush();
      gwMassAudit_->flush();
      if (transportAudit_) {
        transportAudit_->flush();
      }
      output_->flush();
      log::info(log::msg() << "output written at t = " << t << " s, dt = " << dt
                           << " s (fs " << surface_->lastSolve().iterations << " it, gw "
                           << gw_->lastSolve().iterations << " it)");
    }
    if (checkpointInterval > 0.0 && t >= nextCheckpoint - 1.0e-9) {
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
  output_->flush();
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
    surface_->updateVelocity();
    stepTransport(t, dt, 0.0);

    recordMonitors(t);
    recordMassAudit(t);
    if (transport_) {
      recordTransportAudit(t);
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
      output_->flush();
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
  output_->flush();
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
    stepTransport(t, dtg, gw_->lastDtg());

    recordMonitors(t);
    recordGwMassAudit(t);
    if (transport_) {
      recordTransportAudit(t);
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
      output_->flush();
      log::info(log::msg() << "output written at t = " << t << " s, dtg = " << dtg_
                           << " s (solver " << gw_->lastSolve().iterations << " it, residual "
                           << gw_->lastSolve().residualNorm << ")");
    }
    if (checkpointInterval > 0.0 && t >= nextCheckpoint - 1.0e-9) {
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
  output_->flush();
}

}  // namespace frehg::driver
