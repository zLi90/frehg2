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
#include "io/Monitor.hpp"
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
  /// Collective on \p comm.
  Simulation(MPI_Comm comm, const FrehgConfig& config);

  /// Execute the time loop from t_start (or the restart time) to t_end,
  /// writing outputs, monitors, and checkpoints per the configuration.
  void run();

 private:
  void runSurfaceLoop(real_t t0);
  void runGroundwaterLoop(real_t t0);
  void runCoupledLoop(real_t t0);
  void buildTransport();
  void stepTransport(real_t t, real_t dt, real_t dtgLast);
  void writeOutputs(real_t t);
  void recordMonitors(real_t t);
  void recordMassAudit(real_t t);
  void recordGwMassAudit(real_t t);
  void recordTransportAudit(real_t t);
  void writeCheckpoint(real_t t, long step, real_t labelTime = -1.0);
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
  std::unique_ptr<io::Hdf5Output> output_;
  std::unique_ptr<io::Checkpoint> checkpoint_;
  std::vector<std::unique_ptr<io::Monitor>> monitors_;
  std::unique_ptr<io::Monitor> massAudit_;
  std::unique_ptr<io::Monitor> gwMassAudit_;
  std::unique_ptr<io::Monitor> transportAudit_;

  /// Current adaptive subsurface step [s] (groundwater-only runs).
  real_t dtg_ = 0.0;

  // Running global volume budgets, maintained on the audit owner rank.
  real_t cumRain_ = 0.0;
  real_t cumEvap_ = 0.0;
  real_t cumOutflow_ = 0.0;
  real_t cumBcInflow_ = 0.0;
  real_t cumSeepage_ = 0.0;
  real_t cumClamped_ = 0.0;  ///< running below-bed clamp creation [m^3]
  real_t cumGwBoundary_ = 0.0;
  real_t cumGwStorage_ = 0.0;
  real_t cumGwRealloc_ = 0.0;
  real_t cumGwDropped_ = 0.0;
  real_t cumGwVloss_ = 0.0;
  real_t cumTrExchange_ = 0.0;
  real_t cumTrSurfSource_ = 0.0;
  real_t cumTrSurfBoundary_ = 0.0;
  real_t cumTrSurfAdjust_ = 0.0;
  real_t cumTrSurfAnchor_ = 0.0;
  real_t cumTrSubsBoundary_ = 0.0;
  real_t cumTrSubsAdjust_ = 0.0;
  real_t cumTrSubsAnchor_ = 0.0;
};

}  // namespace frehg::driver

#endif  // FREHG_DRIVER_SIMULATION_HPP
