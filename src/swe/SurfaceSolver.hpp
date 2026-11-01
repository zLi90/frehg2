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

#include "atm/MetForcing.hpp"
#include "atm/EtFormulas.hpp"
#include "swe/WindForcing.hpp"
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
  /// Stage-3 T2 (W2.4): coupled deposits (applied seepage, bounce-back,
  /// vent) that landed on prescribed-stage cells after the solve
  /// [m^3, + = removed from / external to the surface]. Under closure
  /// external the deposit is removed before the velocity update (the
  /// stage cell is an external reservoir); under zero_gradient it is
  /// booked only and the dynamics stay bitwise v2.0 — either way the
  /// coupled ledger closes with stage conditions present (the b6
  /// "pre-existing unclosed audit term" of dod-Q7).
  real_t stageDeposit = 0.0;
  /// \name Stage-3 T9 (W9.1, plan §12A): positive thin-layer removals of
  /// the step [m^3], per path. P1 = cflLimiter's isolated-cell rule, P2 =
  /// its small-depth rule (both booked in every mode — at T8 no audit
  /// column held them, the S3-A12 ledger gap), P3 = the positive part of
  /// the evaprain clamp (informational under discard, where it stays
  /// inside clampVolume bitwise), sweep = the W9.2 end-of-step park.
  /// Under the paired store (W9.2) every path parks instead of deleting.
  ///@{
  real_t filmIsolated = 0.0;  ///< P1 volume removed (or parked)
  real_t filmSmall = 0.0;     ///< P2 volume removed (or parked)
  real_t filmClamp = 0.0;     ///< P3 positive part (parked in store mode)
  real_t filmSweep = 0.0;     ///< W9.2 sweep park (store mode only)
  real_t storeRelease = 0.0;  ///< W9.2 water returned from the store
  ///@}
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

  /// Stage-3 T2: the momentum source kernel, compiled once per closure
  /// regime (Ext = false folds to the v2.0 expression trees — Momentum.cpp).
  template <bool Ext>
  void momentumSourceImpl();
  /// Stage-3 T2: the velocity interpolation kernel, same Ext dispatch.
  template <bool Ext>
  void interpolateVelocityImpl();

  /// Stage-3 T2 (W2.4): book this step's coupled deposits \p deposit [m]
  /// (the coupler's gain field) on prescribed-stage cells into
  /// audit().stageDeposit, and remove them from eta on EXTERNAL-closure
  /// stage cells (the reservoir absorbs them before updateVelocity can
  /// read the kick). The coupler calls it between the seepage application
  /// and the velocity update. No eta condition: a no-op.
  void routeStageDeposits(const Field2<real_t>& deposit);

  /// \name Stage-3 T9 (W9.1/W9.2, plan §12A): thin-layer removals and the
  /// paired water store
  ///@{
  /// This step's counted thin-layer removal per cell [m of depth]: P1 + P2
  /// in every mode, plus P3's positive part and the sweep under the paired
  /// store — the film_removed column's per-cell source.
  const Field2<real_t>& filmRemovedDepth() const { return filmRemoved_; }
  /// Instantaneous volume of films 0 < eta - bottom <= min_depth present in
  /// eta [m^3, rank-local]: updateDepth reports them dry, so `volume` omits
  /// them until P2 removes them at the next solve. Identically zero after
  /// the W9.2 sweep (the R8 predicate made exact).
  real_t pendingFilmVolume() const;
  /// Instantaneous water on prescribed-stage cells, Sigma dept A [m^3,
  /// rank-local]. boundary_outflow books the flux INTO those cells (they
  /// are external reservoirs), so the surface identity's domain is
  /// volume - stage_volume; at T8 the stage storage rode the tide unbooked
  /// and the identity closed only at mean-stage instants (S3-A12/S3-A13).
  real_t stageVolume() const;
  /// True when the paired water–salt store is active (surface water +
  /// transport + surface_dryout: store, plan §12A).
  bool pairedStore() const { return pairedStore_; }
  /// W9.2 end-of-step sweep and release, store mode only: park every film
  /// 0 < eta - bottom <= min_depth, then return the store to cells holding
  /// water (R5) or to dry cells whose store reaches wetting_face_depth
  /// (R6). Called after the last eta change of the surface step (after the
  /// coupler's seepage and stage routing; directly after the free-surface
  /// phase when uncoupled). A no-op when the store is inactive.
  void sweepAndRelease();
  const Field2<real_t>& waterStore() const { return waterStore_; }      ///< V_s per cell [m^3]
  Field2<real_t>& waterStoreMutable() { return waterStore_; }           ///< restart load
  const Field2<real_t>& parkedVolume() const { return parkedVol_; }     ///< parked this step [m^3]
  const Field2<real_t>& releasedVolume() const { return releasedVol_; } ///< released this step [m^3]
  /// Water evaporation removed from each cell this step [m^3] (store mode;
  /// the R4(a) crust predicate for every evaporation mode).
  const Field2<real_t>& evaporatedVolume() const { return evapStep_; }
  real_t ownedStoreVolume() const;  ///< Sigma V_s [m^3, rank-local]
  ///@}

  /// \name Stage-3 T10 (W10.1, plan §12B): the solve's continuity face
  /// volumes. Recorded right after fillAndSolve (post-solve eta, the
  /// start-of-step face geometry the solve assembled with), only under
  /// transport.surface_update: conservative — the default path takes no
  /// extra exchange and no extra kernel.
  ///@{
  bool solveVolumesActive() const { return recordSolveVolumes_; }
  /// Eastward face volume of the step [m^3]; slot (j, i) is the face
  /// between cells i and i+1, slot (j, 0) the west interface/edge face.
  const Field2<real_t>& solveVolumeX() const { return Qx_; }
  /// Northward face volume of the step [m^3]; slot layout as Fv.
  const Field2<real_t>& solveVolumeY() const { return Qy_; }
  ///@}

  /// \name Stage-3 T2 (W2.6): entrance-head diagnostics (rank-local)
  ///@{
  real_t stageKeInflux() const { return stageKeInflux_; }          ///< cumulative KE imported through stage faces [J]
  real_t stageHeadRatioNum() const { return stageHeadRatioNum_; }  ///< entrance-head ratio numerator [m^4/s]
  real_t stageHeadRatioDen() const { return stageHeadRatioDen_; }  ///< entrance-head ratio denominator [m^4/s]
  bool hasExternalClosure() const { return hasExternal_; }         ///< any BC uses closure external
  ///@}

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

  /// \name Stage-3 T5 (plan §9 W5.3): the ET open-water channel state
  /// (valid when the ET block enables the surface channel).
  ///@{
  bool etActive() const { return etOpenWater_; }            ///< ET surface channel on
  /// Stage-3 T6 (W6.3): attach the canopy transmission field (the gw
  /// module's vegCanopyTau; scales the open-water rate per cell).
  void attachCanopyFactor(const Field2<real_t>& tau) {
    hasCanopy_ = true;
    canopyTau_ = tau;
  }
  const Field2<real_t>& evapRemoved() const { return evapRemoved_; }  ///< this step's removal [m]
  const Field2<real_t>& rainAdded() const { return rainAdded_; }      ///< this step's masked rain [m]
  const Field2<real_t>& wetFrac() const { return wetFrac_; }          ///< wet fraction of the step [-]
  const Field2<real_t>& wetClock() const { return wetClock_; }        ///< cumulative wet time [s]
  const Field2<real_t>& wetClockStepStart() const { return wetClockStepStart_; }  ///< clock at step start [s]
  ///@}
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
  /// Cumulative fs_ solver telemetry (v2 plan §2.2.6; the g2/s-gates parse
  /// the driver's end-of-run summary line built from this).
  const SolverTelemetry& solverTelemetry() const { return system_->telemetry(); }

  /// \return the resolved PETSc matrix backend ("aij" | "aijkokkos") for
  /// the run record's solver section (v2 plan §2B.5).
  const std::string& solverMatType() const { return system_->matType(); }

  /// \return the resolved BoomerAMG coarsening/smoother (empty unless the
  /// preconditioner is "amg") for the run record (v2 plan §2B.2 B2).
  const std::string& solverAmgCoarsenType() const { return system_->amgCoarsenType(); }
  /// \return the resolved BoomerAMG relaxation (smoother) choice; same
  /// contract and run-record destination as solverAmgCoarsenType().
  const std::string& solverAmgRelaxType() const { return system_->amgRelaxType(); }

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
  /// Read/stage the bathymetry (cell bottom elevations, with the
  /// load-bearing elevation offset) from the config inputs.
  void readBathymetry(const FrehgConfig& config);
  /// Stage the rain coverage mask (whole domain or rain polygon) on device.
  void buildRainMask(const FrehgConfig& config);
  /// Stage the evaporation application mask (v2 Q4: the prescribed-mode
  /// exclude region; all-ones without one — bitwise the unmasked legacy).
  void buildEvapMask(const FrehgConfig& config);
  /// Set the initial eta/velocity state from the config and derive the
  /// consistent depth and face geometry.
  void applyInitialConditions(const FrehgConfig& config);
  /// Stage the surface boundary-condition member-cell lists on device.
  void buildBoundaryLists(const BoundarySet& boundaries);
  /// Build the fixed 5-point COO sparsity pattern handed to
  /// MatSetValuesCOO on every assembly.
  void buildCooPattern();

  // Free-surface phase pieces (FreeSurface.cpp).
  /// Eta clamp, prescribed stages, and ghost fills. A restart refresh
  /// skips the stage prescription: the checkpointed eta already satisfies
  /// it up to the sources applied after the last solve (rain on the
  /// stage cells), and re-prescribing changes eta^n of the first restarted
  /// step — a rounding-class seed the coupled trajectory chaos amplifies
  /// (found at P3; b1's rain exclusion on its tide row masked it).
  void enforceSurfBc(bool prescribeStage);
  /// Free-surface RHS from the current eta and explicit momentum fluxes.
  void assembleRhs();
  /// Free-surface matrix coefficients from the face flow areas and the
  /// semi-implicit gravity/drag factors.
  void assembleCoefficients();
  /// Apply the outflow-condition row corrections before the solve.
  void applyOutflowCorrections();
  /// Fill the PETSc system and run the eta solve.
  void fillAndSolve();
  /// Accumulate the boundary face volumes of the step (inflow, outflow,
  /// prescribed stages) into the mass audit.
  void accumulateBoundaryFluxes();

  // Momentum pieces (Momentum.cpp).
  /// Explicit momentum operators Fu/Fv (advection, diffusion, sources)
  /// for the semi-implicit step.
  void momentumSource();
  /// Update the drag coefficients from the current depths (Manning form,
  /// with the drag limit-cycle clamp found at P3).
  void updateDragCoef();
  /// Back-substitute the face velocities from the solved eta gradient.
  void updateVelocityField();
  /// Cell-centered uy/vx interpolation (the four-point stencil with one
  /// diagonal neighbor; needs the corner halos).
  void interpolateVelocity();

  // Wet/dry and geometry pieces (WetDry.cpp).
  /// Depth from eta and bathymetry with the wet/dry classification
  /// against the thin-layer threshold.
  void updateDepth();
  /// Face wet areas and geometry factors from the updated depths.
  void updateGeometry();
  /// Track the per-face CFL numbers of the step (feeds maxCfl()).
  void cflLimiter();
  /// Stage-3 T9 (W9.1): book cflLimiter's removals per path from the eta
  /// snapshot taken before it (P1 = cells it flagged isolated, P2 = the
  /// rest). A separate read-only pass: the pinned limiter kernels stay
  /// byte-identical (the S3-A8.1/A9.7 fma rule). Under the paired store
  /// the booked depth also parks into waterStore_.
  void bookLimiterRemovals();
  /// Clamp non-physical face velocities (the legacy limiter set; clamped
  /// volume is audited).
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
  /// Stage-3 T10 (W10.1): fill Qx_/Qy_ from the assembled continuity —
  /// Q = dt As E - S (eta_nb - eta) on the post-solve eta — with the
  /// edge folds, velocity-BC faces and outflow drops applied (§12B C1).
  void recordSolveFaceVolumes();
  /// Velocity ghost fills, prescribed faces, and the mass-consistent
  /// stage-boundary face velocities per \p apply.
  void enforceVeloBc(VeloBcApply apply);
  /// Fill the four velocity halo corners the uy/vx interpolation stencil
  /// reads (pairs with HaloExchanger::exchangeWithCorners()).
  void fillVelocityGhostCorners();

  // Sources (SurfaceSources.cpp).
  /// Apply the rain and evaporation depth sources of the step (their
  /// volumes feed the mass audit).
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
  WindForcing windForcing_;  ///< per-step wind evaluation (v2 Q6)
  bool rainIsSeries_ = false, evapIsSeries_ = false;
  real_t rainConstant_ = 0.0, evapConstant_ = 0.0;
  TimeSeries rainSeries_, evapSeries_;
  bool hasRainExclusion_ = false;
  /// Open-water evaporation mode and its met forcing (v2 Q4, plan §3.2).
  SurfaceWaterConfig::EvapMode evapMode_ = SurfaceWaterConfig::EvapMode::Prescribed;
  atm::MetForcing met_;

  // Current-step forcing values.
  real_t rain_ = 0.0, evap_ = 0.0;
  // Stage-3 T5 (W5.3): the ET open-water channel. etOpenWater_ switches
  // the evapRain() branch; the fields feed the coupler's dry-fraction
  // clock, the transport dilution (D8) and the atmosphere outputs.
  bool etOpenWater_ = false;
  EvapotranspirationConfig etCfg_;
  TimeSeries etPrescribedSeries_;    ///< prescribed-potential series
  /// Stage-3 T7 (W7.3): the PM open-water channel's r_a law (the ET
  /// consumers' configured law; the bulk channel keeps evaporationRate's
  /// pinned liu2006 chain).
  atm::AeroParams etAero_;
  AtmosphereConfig::AeroResistance etAeroLaw_ =
      AtmosphereConfig::AeroResistance::Liu2006;
  bool hasCanopy_ = false;           ///< a canopy factor is attached (T6)
  Field2<real_t> canopyTau_;         ///< Beer-Lambert transmission [-]
  Field2<real_t> evapRemoved_;       ///< per-cell removal this step [m]
  Field2<real_t> rainAdded_;         ///< per-cell masked rain this step [m]
  Field2<real_t> wetFrac_;           ///< wet fraction of this step [-]
  Field2<real_t> wetClock_;          ///< cumulative wet time [s] (checkpointed)
  Field2<real_t> wetClockStepStart_; ///< wetClock before this step's increment
  /// \name Stage-3 T9 (W9.1/W9.2): thin-layer removal books and the
  /// paired water store (plan §12A). etaPreLimit_ and filmRemoved_ exist in
  /// every mode (W9.1 books P1/P2 always); the store fields only when
  /// pairedStore_ (empty views otherwise — every store path is inert).
  ///@{
  Field2<real_t> etaPreLimit_;   ///< eta before cflLimiter (P1/P2 snapshot)
  Field2<real_t> filmRemoved_;   ///< counted removal this step [m]
  bool pairedStore_ = false;     ///< surface + transport + surface_dryout: store
  Field2<real_t> waterStore_;    ///< V_s [m^3] (checkpointed: swe_water_store)
  Field2<real_t> parkedVol_;     ///< volume parked this step [m^3]
  Field2<real_t> releasedVol_;   ///< volume released this step [m^3]
  Field2<real_t> etaPreEvap_;    ///< eta before the evaporation pass
  /// Stage-3 T10 (W10.1): the solve's face volumes (empty unless
  /// recordSolveVolumes_).
  bool recordSolveVolumes_ = false;
  Field2<real_t> Qx_, Qy_;
  Field2<real_t> evapStep_;      ///< water evaporation removed this step [m^3]
  ///@}
  real_t windSpeed_ = 0.0;   ///< |U10| of the step [m/s]
  real_t windOmega_ = 0.0;   ///< wind direction of the step, radians from +x
  real_t windCd_ = 0.0;      ///< Cd of the step under the configured law

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
  /// \name Stage-3 T2: the external inflow closure (plan §6/§13.A)
  /// Allocated/registered unconditionally (the halo set is fixed at
  /// construction); the kernels read them only when hasExternal_, so runs
  /// without an external closure take the v2.0 code paths bitwise.
  ///@{
  bool hasExternal_ = false;    ///< any eta BC declares closure external
  /// 1.0 on stage cells whose BC is external (static, corner-exchanged).
  Field2<real_t> isExternalStage_;
  /// 1.0 on reservoir x+/y+ faces: both cells external-stage, or one
  /// external-stage and the other outside the global domain (static,
  /// corner-exchanged; the momentum stencils read them, W2.3).
  Field2<real_t> resFaceX_, resFaceY_;
  /// The external reservoir's velocity at the reservoir faces, refreshed
  /// per step from each BC's inflow series (corner-exchanged).
  Field2<real_t> uResX_, vResY_;
  /// Entrance-head diagnostics (W2.6): cumulative stage-face KE influx [J]
  /// and the paid-head ratio accumulators (rank-local; the run record
  /// reduces them).
  real_t stageKeInflux_ = 0.0;
  real_t stageHeadRatioNum_ = 0.0;
  real_t stageHeadRatioDen_ = 0.0;
  ///@}
  Field2<real_t> rainMask_;     ///< 1.0 where rain applies (exclusion region = 0)
  Field2<real_t> evapMask_;     ///< 1.0 where evaporation applies (v2 Q4 mask)
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
  /// Allreduced wet-mask checksum of the previous solve (v2 plan §2.2.4):
  /// a change forces an AMG-hierarchy rebuild. Only computed when the
  /// selected preconditioner reuses its hierarchy.
  long long wetMaskHash_ = -1;

  SurfaceStepAudit audit_;
};

}  // namespace frehg::swe

#endif  // FREHG_SWE_SURFACESOLVER_HPP
