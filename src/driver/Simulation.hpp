/// \file Simulation.hpp
/// \brief The time loop and module orchestration (plan §4 driver layer).
///
/// P1 + P2 + P3 + P4 scope: the full legacy solve() loop (solve.c:25-166) —
/// surface-water-only, groundwater-only, and coupled runs, each with the
/// optional scalar-transport block (solve.c:112-116) after the flow step.
/// Groundwater-only runs march on the adaptive subsurface step dtg exactly
/// as legacy does (solve.c:37 initializes dtg,
/// and solve_groundwater:193 feeds the adapted dtg back into the outer step
/// for runs without a surface module); sync-coupled runs march both modules
/// on the same adaptive step (amendment A10), subcycled runs keep the fixed
/// surface dt and let the coupler subcycle the subsurface. The driver owns
/// the grid, subsurface mesh, boundary set, coupler, output file, monitors,
/// and checkpointing, and mediates every module's I/O.

#ifndef FREHG_DRIVER_SIMULATION_HPP
#define FREHG_DRIVER_SIMULATION_HPP

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Types.hpp"
#include "coupling/Coupler.hpp"
#include "gw/RichardsSolver.hpp"
#include "gw/TerrainMetric.hpp"
#include "io/Checkpoint.hpp"
#include "io/Hdf5Output.hpp"
#include "io/IntervalRates.hpp"
#include "io/Monitor.hpp"
#include "io/RunRecord.hpp"
#include "swe/SurfaceSolver.hpp"
#include "transport/ScalarSolver.hpp"

#include <mpi.h>

#include <memory>
#include <vector>

namespace frehg::driver {

/// One configured model run.
class Simulation {
 public:
  /// Build the grid, boundary conditions, modules, and output file.
  /// Collective on \p comm. \p inputPath is the configuration file the
  /// run was launched with (recorded in the run record, v2 plan §2A).
  Simulation(MPI_Comm comm, const FrehgConfig& config,
             const std::string& inputPath = std::string());

  /// Execute the time loop from t_start (or the restart time) to t_end,
  /// writing outputs, monitors, and checkpoints per the configuration.
  void run();

  /// Write the final run record (provenance.finished = true). Collective;
  /// call after run() returns so the record carries the completed
  /// "simulation" timer (v2 plan §2A).
  void finalizeRunRecord();

 private:
  void runSurfaceLoop(real_t t0);
  void runGroundwaterLoop(real_t t0);
  void runCoupledLoop(real_t t0);
  void buildTransport();
  void stepTransport(real_t t, real_t dt, real_t dtgLast);
  /// Write field outputs under the dataset key \p label; \p tActual (< 0:
  /// same as the label) is the true state time the interval-rate outputs
  /// normalize over (Stage-3 T4 W4.1 — adaptive loops can cross a label).
  void writeOutputs(real_t label, real_t tActual = -1.0);
  void recordMonitors(real_t t);
  /// The audit recorders ACCUMULATE every step; \p writeRow decimates the
  /// row write only (Stage-3 T4 W4.3, output.monitor_interval).
  void recordMassAudit(real_t t, bool writeRow = true);
  void recordGwMassAudit(real_t t, bool writeRow = true);
  void recordTransportAudit(real_t t, bool writeRow = true);
  void recordTemperatureAudit(real_t t, bool writeRow = true);
  /// True when a monitor/audit row is due at \p t: every step with
  /// monitor_interval 0/absent, else at the first step reaching the next
  /// interval multiple, or \p forced (a checkpoint instant or t_end).
  bool monitorRowDue(real_t t, bool forced);
  /// Stage-3 T5 (W5.10): per-outer-step ET contributions to the
  /// atmosphere interval rates (and the soil depth-window reset).
  void addEtRateContributions();
  /// Stage-3 T9 (W9.1): add the step's counted thin-layer removal to the
  /// film_removal_rate accumulator (when requested).
  void addFilmRateContribution();
  void writeCheckpoint(real_t t, long step, real_t labelTime = -1.0);
  /// Refresh the run record's dynamic sections (solver telemetry, closure)
  /// and rewrite it. Collective (v2 plan §2A).
  void flushRunRecord(bool finished);
  real_t restoreFromCheckpoint();
  io::Checkpoint::Fields2 checkpointFields2() const;
  io::Checkpoint::Fields3 checkpointFields3() const;

  FrehgConfig config_;
  Grid grid_;
  std::unique_ptr<HaloExchanger> halo_;
  std::unique_ptr<gw::TerrainMetric> mesh_;
  std::unique_ptr<BoundarySet> boundaries_;
  std::unique_ptr<swe::SurfaceSolver> surface_;
  std::unique_ptr<gw::RichardsSolver> gw_;
  std::unique_ptr<coupling::Coupler> coupler_;
  std::unique_ptr<transport::ScalarSolver> transport_;
  std::unique_ptr<transport::ScalarSolver> temperature_;  ///< second scalar (Q5)
  std::unique_ptr<io::Hdf5Output> output_;
  std::unique_ptr<io::Checkpoint> checkpoint_;
  std::unique_ptr<io::RunRecord> runRecord_;
  std::vector<std::unique_ptr<io::Monitor>> monitors_;
  std::unique_ptr<io::IntervalRates> rates_;  ///< interval-mean rate outputs (T4)
  std::unique_ptr<io::Monitor> massAudit_;
  std::unique_ptr<io::Monitor> gwMassAudit_;
  std::unique_ptr<io::Monitor> transportAudit_;
  std::unique_ptr<io::Monitor> temperatureAudit_;

  /// Current adaptive subsurface step [s] (groundwater-only runs).
  real_t dtg_ = 0.0;
  /// Next monitor/audit row-write boundary [s] (monitor_interval > 0).
  real_t nextMonitorRow_ = 0.0;

  // Running global volume budgets, maintained on the audit owner rank.
  real_t cumRain_ = 0.0;
  real_t cumEvap_ = 0.0;
  real_t cumOutflow_ = 0.0;
  real_t cumBcInflow_ = 0.0;
  real_t cumSeepage_ = 0.0;
  real_t cumClamped_ = 0.0;  ///< running below-bed clamp creation [m^3]
  /// Running coupled deposit on prescribed-stage cells [m^3] (Stage-3 T2).
  real_t cumStageDeposit_ = 0.0;
  /// Stage-3 T9 (W9.1/W9.2): the running counted thin-layer removal, its
  /// per-path books (run-record lines; P3 is informational under discard)
  /// and the store's running release [m^3].
  real_t cumFilmRemoved_ = 0.0;
  real_t cumFilmIsolated_ = 0.0;
  real_t cumFilmSmall_ = 0.0;
  real_t cumFilmClamp_ = 0.0;
  real_t cumFilmSweep_ = 0.0;
  real_t cumStoreRelease_ = 0.0;
  real_t cumGwVentThin_ = 0.0;  ///< W9.4 thin vents deposited (P4) [m^3]
  real_t cumGwBoundary_ = 0.0;
  real_t cumGwStorage_ = 0.0;
  real_t cumGwRealloc_ = 0.0;
  real_t cumGwDropped_ = 0.0;
  real_t cumGwVloss_ = 0.0;
  real_t cumGwEvap_ = 0.0;
  real_t cumGwEtApplied_ = 0.0;      ///< ET soil channel applied [m^3] (T5)
  real_t cumGwEtUnrealised_ = 0.0;   ///< ET soil channel unrealised [m^3] (T5)
  real_t cumGwTrApplied_ = 0.0;      ///< transpiration applied [m^3] (T6)
  real_t cumGwTrUnrealised_ = 0.0;   ///< transpiration unrealised [m^3] (T6)
  real_t cumTrExchange_ = 0.0;
  real_t cumTrSurfSource_ = 0.0;
  real_t cumTrSurfBoundary_ = 0.0;
  real_t cumTrSurfAdjust_ = 0.0;
  real_t cumTrSurfAnchor_ = 0.0;
  real_t cumTrSubsBoundary_ = 0.0;
  real_t cumTrSubsAdjust_ = 0.0;
  real_t cumTrSubsAnchor_ = 0.0;
  // Temperature-instance cumulative budgets (the heat ledger, v2 Q5).
  real_t cumTpExchange_ = 0.0;
  real_t cumTpSurfSource_ = 0.0;
  real_t cumTpSurfBoundary_ = 0.0;
  real_t cumTpSurfAdjust_ = 0.0;
  real_t cumTpSurfAnchor_ = 0.0;
  real_t cumTpSurfAtmos_ = 0.0;
  real_t cumTpSubsBoundary_ = 0.0;
  real_t cumTpSubsAdjust_ = 0.0;
  real_t cumTpSubsAnchor_ = 0.0;
};

}  // namespace frehg::driver

#endif  // FREHG_DRIVER_SIMULATION_HPP
