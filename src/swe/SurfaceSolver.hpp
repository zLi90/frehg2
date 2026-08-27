/// \file SurfaceSolver.hpp
/// \brief Semi-implicit theta-scheme shallow-water module (plan §10 P1).
///
/// The solver preserves legacy Frehg's SWE algorithm exactly (plan §3.1;
/// provenance per function in the .cpp files): CFL-damped upwind advection,
/// central eddy viscosity, point-implicit Manning drag with the thin-layer
/// exponent switch (plus the plan §5.8 Chezy extension), quadratic wind
/// stress with thin-layer attenuation, the implicit 5-point free-surface
/// system, min-depth wetting/drying with the one-cell wetting limiter and
/// higher-of-two-bottoms face depths, and rain/evaporation applied to eta.
///
/// A time step runs in two phases, mirroring the legacy call graph
/// (solve.c:51 / solve.c:97) so the groundwater module can slot between them
/// in P3:
///   beginStep(t)        - save eta^n, evaluate forcings at t
///   solveFreeSurface()  - legacy solve_shallowwater (shallowwater.c:53)
///   updateVelocity()    - legacy shallowwater_velocity (shallowwater.c:95)
///
/// Differences from legacy, all recorded in docs/theory/surface-water.md:
/// the free-surface system is one global PETSc solve instead of per-rank
/// blocks with lagged Dirichlet interfaces (rank invariance, plan §8.2);
/// eta-condition columns are eliminated to keep the matrix SPD (identical
/// solution); the uy/vx interpolation uses the intended four-point stencil
/// (the legacy index arithmetic wrapped rows and read out of bounds at
/// domain edges); rainfall exclusion is an explicit configured region
/// instead of the hardcoded last-row skip.

#ifndef FREHG_SWE_SURFACESOLVER_HPP
#define FREHG_SWE_SURFACESOLVER_HPP

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/LinearSystem.hpp"
#include "core/TimeSeries.hpp"
#include "core/Types.hpp"

#include <memory>
#include <string>
#include <vector>

namespace frehg::swe {

/// One boundary condition's member cells staged for device kernels.
struct DeviceBcList {
  const BoundaryCondition* bc = nullptr;      ///< source condition (host)
  Kokkos::View<int*, MemSpace> j;             ///< local row indices (halo layout)
  Kokkos::View<int*, MemSpace> i;             ///< local column indices
  Kokkos::View<int*, MemSpace> face;          ///< BcFace as int (velocity/outflow kinds)
  Kokkos::View<real_t*, MemSpace> bedDrop;    ///< outflow kind: bed drop across the edge cell [m]
  real_t current = 0.0;                       ///< value at the current step time
};

/// This rank's share of the per-step volume budget [m^3]; the driver reduces
/// across ranks and appends the running totals to the mass-audit monitor.
struct SurfaceStepAudit {
  real_t rainVolume = 0.0;       ///< rain volume added to eta
  real_t evapVolume = 0.0;       ///< prescribed evaporative volume
  real_t boundaryOutflow = 0.0;  ///< outflow through open edges and eta-condition cells
  real_t bcInflow = 0.0;         ///< inflow injected by discharge/velocity conditions
  /// Volume created by the legacy below-bed clamp (shallowwater.c:507-511
  /// lifts eta onto the bed with no compensating flux). Negligible in
  /// b1-b4's regimes; b5's hard-drawn films make it a measurable legacy
  /// mass defect (~1 % of rain), so the audit measures it instead of
  /// leaving it in the closure residual.
  real_t clampVolume = 0.0;
};

/// The shallow-water module.
class SurfaceSolver {
 public:
  /// Allocate fields, read the bathymetry, rasterize surface conditions to
  /// device lists, set the initial state, and fix the free-surface COO
  /// pattern. Collective on grid.comm().
  /// \param grid decomposition (kept by reference, must outlive the solver);
  ///        buildGlobalIds must already have been called.
  /// \param config the loaded configuration.
  /// \param boundaries rasterized boundary conditions (kept by reference).
  /// \param halo exchanger the solver registers its fields with (kept by
  ///        reference; the solver owns the registered names' lifecycle).
  SurfaceSolver(const Grid& grid, const FrehgConfig& config, const BoundarySet& boundaries,
                HaloExchanger& halo);

  /// Save eta^n and evaluate rain/evaporation/wind/boundary series at \p t.
  /// Legacy evaluates forcings at the end-of-step time (solve.c:43-49).
  void beginStep(real_t t);

  /// Set the step length used by every subsequent phase. Surface-only and
  /// subcycled-coupling runs keep the configured time.dt; sync-coupled runs
  /// march on the common adaptive step (legacy solve.c:193 feeds the adapted
  /// dtg back into dt — amendment A10) and call this before beginStep.
  void setTimeStep(real_t dt) { dt_ = dt; }

  /// Recompute cell and face depths from the current eta (legacy calls
  /// update_depth at the end of solve_shallowwater, shallowwater.c:87, so
  /// the groundwater phase sees the post-solve depth). Interior values are
  /// authoritative; halo depths refresh in updateVelocity after the eta
  /// exchange.
  void refreshDepth() { updateDepth(); }

  /// The implicit free-surface phase: momentum source, matrix assembly and
  /// solve, wetting limiter, rain/evaporation (legacy solve_shallowwater).
  void solveFreeSurface();

  /// The velocity phase: depth and face geometry, drag update, velocity
  /// update with wet/dry limiters, ghost fills, and the uy/vx interpolation
  /// (legacy shallowwater_velocity).
  void updateVelocity();

  /// \name State access (device views; interiors are authoritative)
  /// The coupler writes eta and depth through these views exactly where
  /// legacy did (subsurface_source, the saturation bounce-back, and the
  /// reallocation vent all mutate the shared surface state mid-step).
  ///@{
  const Field2<real_t>& eta() const { return eta_; }     ///< free surface [m], offset frame
  const Field2<real_t>& depth() const { return dept_; }  ///< cell-center depth [m]
  const Field2<real_t>& uu() const { return uu_; }       ///< x+ face velocity [m/s]
  const Field2<real_t>& vv() const { return vv_; }       ///< y+ face velocity [m/s]
  const Field2<real_t>& bottom() const { return bottom_; }  ///< bed elevation [m], offset frame
  ///@}

  /// The internal elevation shift: legacy lifts all elevations by
  /// -min(bottom) so they are non-negative (read_bathymetry,
  /// initialize.c:142-147). Preserved as the shared elevation frame of both
  /// modules (the coupler's depth/head exchange assumes it); outputs
  /// subtract it. The legacy dry-cell row made the shift algorithmically
  /// significant — the P3 continuity closure (amendment A13) is
  /// translation-invariant, and b1 gates green either way.
  real_t elevationOffset() const { return offset_; }

  /// eta - elevationOffset() staged into a scratch field for output.
  const Field2<real_t>& etaAbsolute();

  /// \name Transport-module state access (plan §10 P4; wired by the driver)
  ///@{
  const Field2<real_t>& flowRateX() const { return Fu_; }  ///< uu Asx [m^3/s]
  const Field2<real_t>& flowRateY() const { return Fv_; }  ///< vv Asy [m^3/s]
  const Field2<real_t>& faceAreaX() const { return Asx_; } ///< x+ face area [m^2]
  const Field2<real_t>& faceAreaY() const { return Asy_; } ///< y+ face area [m^2]
  /// 1.0 on eta-condition cells (the transport Dirichlet pairing).
  const Field2<real_t>& etaBcMask() const { return isEtaBc_; }
  /// 1.0 where rain applies (the configured exclusion region is 0).
  const Field2<real_t>& rainApplyMask() const { return rainMask_; }
  real_t minDepth() const { return minDepth_; }      ///< legacy min_dept
  real_t currentRain() const { return rain_; }       ///< rate of the step [m/s]
  real_t currentEvaporation() const { return evap_; }  ///< rate of the step [m/s]
  ///@}

  /// \name Restart state (plan §7): prognostic fields beyond eta/uu/vv
  ///@{
  const Field2<real_t>& cflX() const { return cflx_; }  ///< |u| dt / dx of the last step
  const Field2<real_t>& cflY() const { return cfly_; }  ///< |v| dt / dy of the last step
  /// eta at the start of the last step. Checkpointed so a restart refresh
  /// can reconstruct the stage-boundary velocity correction's west/south
  /// *edge-slot* writes bitwise (they feed the uy/vx interpolation but
  /// live in halo slots the §7 interior-only checkpoint cannot carry;
  /// found at P3 — the reconstruction gap flipped the outlet limit-cycle
  /// phase and broke coupled restart determinism).
  const Field2<real_t>& etaStart() const { return etan_; }
  ///@}

  /// Rebuild every derived field (depth, geometry, drag, uy/vx) from the
  /// prognostic state after a restart read. Collective on grid.comm().
  void refreshDerivedState();

  /// Convergence data of the last free-surface solve.
  const SolveStats& lastSolve() const { return lastSolve_; }

  /// This rank's volume-budget contributions of the current step.
  const SurfaceStepAudit& audit() const { return audit_; }

  /// This rank's owned surface volume Sum(dept dx dy) [m^3]
  /// (legacy get_mass, solve.c:363).
  real_t ownedVolume() const;

  /// This rank's maximum face CFL number of the last step (legacy monitors
  /// it; the driver reduces and warns).
  real_t maxCfl() const;

  // NOTE: the step-phase helpers below are implementation detail, kept
  // public only because nvcc forbids extended __host__ __device__ lambdas
  // (KOKKOS_LAMBDA) inside private member functions (CUDA C++ programming
  // guide, "Extended Lambda Restrictions"). Treat as private.
 public:
  // Initialization helpers (SurfaceSolver.cpp).
  void readBathymetry(const FrehgConfig& config);
  void buildRainMask(const FrehgConfig& config);
  void applyInitialConditions(const FrehgConfig& config);
  void buildBoundaryLists(const BoundarySet& boundaries);
  void buildCooPattern();

  // Free-surface phase pieces (FreeSurface.cpp).
  /// Eta clamp, prescribed stages, and ghost fills. A restart refresh
  /// skips the stage prescription: the checkpointed eta already satisfies
  /// it up to the sources applied after the last solve (rain on the
  /// stage cells), and re-prescribing changes eta^n of the first restarted
  /// step — a rounding-class seed the coupled trajectory chaos amplifies
  /// (found at P3; b1's rain exclusion on its tide row masked it).
  void enforceSurfBc(bool prescribeStage);
  void assembleRhs();
  void assembleCoefficients();
  void applyOutflowCorrections();
  void fillAndSolve();
  void accumulateBoundaryFluxes();

  // Momentum pieces (Momentum.cpp).
  void momentumSource();
  void updateDragCoef();
  void updateVelocityField();
  void interpolateVelocity();

  // Wet/dry and geometry pieces (WetDry.cpp).
  void updateDepth();
  void updateGeometry();
  void cflLimiter();
  void applyVelocityLimiters();
  /// How enforceVeloBc applies the stage-boundary velocity correction.
  enum class VeloBcApply {
    /// During steps: all four global-edge branches.
    Step,
    /// During a restart refresh: only the west/south *edge-slot* writes
    /// are reconstructed (their inputs — the interior Fu/Fv and the
    /// restored eta^n — are bitwise recoverable, and the slots themselves
    /// are halo cells the checkpoint cannot carry). The east/north
    /// branches write interior slots the checkpoint restored already;
    /// re-running them over post-correction state would perturb them
    /// (found at P3; b1's north-edge tide row pins the skip).
    Refresh,
  };
  /// Velocity ghost fills, prescribed faces, and the mass-consistent
  /// stage-boundary face velocities per \p apply.
  void enforceVeloBc(VeloBcApply apply);
  void fillVelocityGhostCorners();

  // Sources (SurfaceSources.cpp).
  void evapRain();

 private:
  const Grid& grid_;
  HaloExchanger& halo_;

  // Configuration extracts.
  real_t dt_ = 0.0;
  real_t gravity_ = 9.81;
  real_t viscX_ = 0.0, viscY_ = 0.0;
  real_t minDepth_ = 0.0;
  real_t wettingFaceDepth_ = 0.0;
  real_t thinLayerDepth_ = 0.0;
  real_t offset_ = 0.0;
  bool chezy_ = false;
  WindConfig windCfg_;
  bool rainIsSeries_ = false, evapIsSeries_ = false;
  real_t rainConstant_ = 0.0, evapConstant_ = 0.0;
  TimeSeries rainSeries_, evapSeries_;
  TimeSeries windSpeedSeries_, windDirectionSeries_;
  bool hasRainExclusion_ = false;

  // Current-step forcing values.
  real_t rain_ = 0.0, evap_ = 0.0;
  real_t windSpeed_ = 0.0, windDirection_ = 0.0;

  // Fields (nyLocal+2, nxLocal+2), (j, i) with one-cell halos.
  Field2<real_t> eta_, etan_, dept_, deptx_, depty_;
  Field2<real_t> bottom_;
  Field2<real_t> uu_, vv_, uy_, vx_;
  Field2<real_t> Ex_, Ey_, Dx_, Dy_, CDx_, CDy_;
  Field2<real_t> Vs_, Vsx_, Vsy_, Asx_, Asy_, Asz_, Aszx_, Aszy_;
  Field2<real_t> Fu_, Fv_;
  Field2<real_t> cflx_, cfly_, cflActive_;
  Field2<real_t> Sxp_, Sxm_, Syp_, Sym_, Sct_, Srhs_;
  Field2<real_t> etaBcValue_;
  Field2<real_t> isEtaBc_;      ///< 1.0 on eta-condition cells (halo-exchanged once)
  Field2<real_t> rainMask_;     ///< 1.0 where rain applies (exclusion region = 0)
  Field2<real_t> frictionCoef_; ///< Manning n or Chezy C per cell
  Field2<real_t> etaOut_;       ///< scratch for offset-corrected eta output

  // Boundary conditions staged for kernels.
  std::vector<DeviceBcList> etaBcs_;
  std::vector<DeviceBcList> dischargeBcs_;
  std::vector<DeviceBcList> velocityBcs_;
  std::vector<DeviceBcList> outflowBcs_;

  // Free-surface linear system.
  std::unique_ptr<LinearSystem> system_;
  Kokkos::View<PetscInt*, MemSpace> cooRows_, cooCols_;
  Kokkos::View<real_t*, MemSpace> cooValues_;
  Kokkos::View<real_t*, MemSpace> rhsVec_, solVec_;
  SolveStats lastSolve_;

  SurfaceStepAudit audit_;
};

}  // namespace frehg::swe

#endif  // FREHG_SWE_SURFACESOLVER_HPP
