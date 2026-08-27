/// \file ScalarSolver.hpp
/// \brief The scalar-transport module (plan §10 P4): explicit FV transport
///        of one scalar on the surface and subsurface grids.
///
/// One transport step reproduces the legacy per-step scalar block
/// (solve.c:112-116) exactly: the surface scalar advances first
/// (scalar_shallowwater, scalar.c:25-298), then the subsurface scalar
/// (scalar_groundwater, scalar.c:303-498), both over the surface step dt
/// with the flow state the flow modules left behind — the surface flow
/// rates and volumes of the completed step, and the subsurface Darcy fluxes
/// of the window's last substep.
///
/// The module owns the scalar fields and all transport scratch; it reads
/// the flow modules' state through driver-wired views (plan §4: physics
/// libraries never include each other) and hands its scalar views back to
/// the groundwater module for the baroclinic activation
/// (RichardsSolver::attachScalar).
///
/// Legacy state carried across the step boundary and its Frehg2 form:
///  - Vsn, the surface volume at the previous velocity phase (the legacy
///    geometry update, shallowwater.c:1050): snapshotted at the end of
///    each transport step from the final depth (vsn_);
///  - Fu/Fv as read by volume_by_flux (shallowwater.c:1003): legacy calls
///    it before the velocity update, so the divergence volume uses the
///    *previous* step's flow rates — snapshotted at the end of each
///    transport step (fuOld_/fvOld_, checkpointed: the end-of-step Fu is
///    the pre-stage-correction value a restart cannot reconstruct from the
///    checkpointed post-correction velocities);
///  - s_surfkP, the top-cell subsurface scalar the surface exchange reads
///    (enforce_scalar_bc, scalar.c:903): refreshed in the subsurface BC
///    pass (reconstructed from the restored scalar on restart);
///  - Dzz at each column's top cell, read by the surface exchange
///    (scalar.c:137): the tensor is rebuilt inside every subsurface step
///    from the last substep's fluxes, so the top-cell value is snapshotted
///    for the next surface step (dzzTop_, checkpointed — the fluxes it came
///    from are not restart state).

#ifndef FREHG_TRANSPORT_SCALARSOLVER_HPP
#define FREHG_TRANSPORT_SCALARSOLVER_HPP

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/Types.hpp"

#include <vector>

namespace frehg::transport {

/// Surface-module state the transport reads (filled by the driver from
/// SurfaceSolver accessors; all views alias the module's storage).
struct SurfaceWiring {
  bool active = false;        ///< true when the surface module runs
  real_t minDepth = 0.0;      ///< surface_water.min_depth (legacy min_dept)
  Field2<real_t> dept;        ///< cell depth of the completed step
  Field2<real_t> etan;        ///< eta at the start of the completed step
  Field2<real_t> bottom;      ///< bed elevation, offset frame
  Field2<real_t> uu;          ///< x+ face velocity (limiter gates)
  Field2<real_t> vv;          ///< y+ face velocity
  Field2<real_t> fu;          ///< x+ face flow rate [m^3/s] (advection)
  Field2<real_t> fv;          ///< y+ face flow rate [m^3/s]
  Field2<real_t> asx;         ///< x+ face area [m^2] (diffusion, limiter gates)
  Field2<real_t> asy;         ///< y+ face area [m^2]
  Field2<real_t> rainMask;    ///< 1.0 where rain applies (dilution volume)
};

/// Subsurface-module state the transport reads (RichardsSolver /
/// TerrainMetric accessors; face-slot layouts as documented there).
struct SubsurfaceWiring {
  bool active = false;        ///< true when the groundwater module runs
  Field3<real_t> wc;          ///< θ after the completed step
  Field3<real_t> wcn;         ///< θ at the start of the last substep (Vgn)
  Field3<real_t> qx;          ///< volumetric x face flux [m^3/s]
  Field3<real_t> qy;          ///< volumetric y face flux [m^3/s]
  Field3<real_t> qzF;         ///< z face fluxes, plane k above cell k
  Field3<real_t> kx;          ///< x face conductivity (limiter gates)
  Field3<real_t> ky;          ///< y face conductivity
  Field3<real_t> kzF;         ///< z face conductivities
  Field3<real_t> wcs;         ///< θs (dispersion porosity factor)
  Field3<real_t> ksz;         ///< saturated Kz (the legacy Ksz > 0 gates)
  Field3<real_t> dz3d;        ///< cell thickness
  Field3<real_t> ax;          ///< x+ face area [m^2]
  Field3<real_t> ay;          ///< y+ face area [m^2]
  Field3<real_t> cosx;        ///< x face slope cosine (dispersive distances)
  Field3<real_t> cosy;        ///< y face slope cosine
  real_t az = 0.0;            ///< horizontal face area dx dy
  Field2<int> topCode;        ///< configured top BC code (GwBcCode)
  Field2<real_t> topValue;    ///< configured top flux (legacy qtop)
  Field2<int> sideCodeYp;     ///< y+ side code (limiter ghost gates)
  Field2<int> sideCodeYm;     ///< y- side code
};

/// Coupler state the transport reads in coupled runs.
struct CouplingWiring {
  bool active = false;  ///< true when the coupler runs
  Field2<real_t> qss;   ///< surface-applied seepage rate of the step [m/s]
};

/// This rank's scalar-mass budget contributions of one transport step
/// [concentration * m^3]. The driver reduces them into the
/// /monitor/transport_audit table. Every non-conservative piece of the
/// legacy scheme is measured, so the closure identities
///   Δ(Σ s dept A) = exchange + surfSource + surfBoundary + surfAdjust
///                   + surfAnchor
///   Δ(Σ s θ V)    = -exchange + subsBoundary + subsAdjust + subsAnchor
/// hold to rounding (the P3 mass-audit convention: defects reported as
/// data, not hidden in the residual). The anchor terms measure the legacy
/// ledger inconsistencies: the surface flux volume lags the velocity update
/// by one step (volume_by_flux runs before update_velocity,
/// shallowwater.c:112-113), and the subsurface flux volume predates the
/// reallocation/clamp θ adjustments — both vanish in quasi-steady flow,
/// which is where the plan's conservation criterion gates.
struct TransportAudit {
  real_t exchange = 0.0;      ///< seepage scalar mass into the surface
  real_t surfSource = 0.0;    ///< inflow mass + Dirichlet (tide) reset delta
  /// Surface advective mass through global-edge faces (positive in). The
  /// closed-edge fold keeps the *water* in, but the boundary-face velocity
  /// is a live computed value and the scalar advects across it against the
  /// zero-gradient ghost — a legacy leak path, measured here.
  real_t surfBoundary = 0.0;
  real_t surfAdjust = 0.0;    ///< surface limiter/bounds/dry-zero mass delta
  real_t surfAnchor = 0.0;    ///< surface ledger re-anchor (actual - flux volume)
  /// Subsurface boundary-face scalar mass in, including the coupled top
  /// interface's one-sided dispersive gain (scalar.c:362 adds it to the
  /// subsurface with no surface counterpart — a measured legacy defect).
  real_t subsBoundary = 0.0;
  real_t subsAdjust = 0.0;    ///< subsurface limiter/bounds mass delta
  real_t subsAnchor = 0.0;    ///< subsurface ledger re-anchor (θ_new V - Vgflux)
};

/// The scalar-transport module (single scalar, plan Appendix A note).
class ScalarSolver {
 public:
  /// Allocate fields, stage the scalar boundary conditions, and apply the
  /// initial conditions (ghosts enforced, halos exchanged). Collective on
  /// grid.comm().
  ScalarSolver(const Grid& grid, const FrehgConfig& config, const BoundarySet& boundaries,
               HaloExchanger& halo, const SurfaceWiring& surface,
               const SubsurfaceWiring& subsurface, const CouplingWiring& coupling);

  /// Advance the scalar by the surface step \p dt, evaluating series-valued
  /// scalar conditions at time \p t. \p dtgLast is the subsurface window's
  /// last substep (the legacy Vgflux step); \p rain and \p evap the surface
  /// module's forcing rates of the completed step [m/s].
  void step(real_t t, real_t dt, real_t dtgLast, real_t rain, real_t evap);

  /// Rebuild ghosts, derived state, and the end-of-step snapshots from the
  /// restored {s_surf, s_subs, fuOld, fvOld, dzzTop} after a checkpoint
  /// read at time \p t. Collective on grid.comm(). Must run before the
  /// groundwater module's own refresh (its baroclinic ratios read the
  /// scalar ghosts).
  void refreshDerivedState(real_t t);

  /// \name State access (device views; interiors are authoritative)
  ///@{
  /// Surface scalar concentration (output variable concentration_surface).
  const Field2<real_t>& surfaceScalar() const { return sSurf_; }
  /// Subsurface scalar concentration (output variable concentration).
  const Field3<real_t>& subsurfaceScalar() const { return sSubs_; }
  /// The surface-side scalar exchange rate of the last step
  /// [concentration m/s] (legacy sseepage observable).
  const Field2<real_t>& scalarSeepage() const { return sseepage_; }
  ///@}

  /// \name Restart state beyond the scalar fields (see the file comment)
  ///@{
  /// Pre-correction flow-rate snapshot (checkpoint field "s_fu_old").
  const Field2<real_t>& flowRateSnapshotX() const { return fuOld_; }
  /// Pre-correction flow-rate snapshot (checkpoint field "s_fv_old").
  const Field2<real_t>& flowRateSnapshotY() const { return fvOld_; }
  /// Carried top-cell Dzz (checkpoint field "s_dzz_top").
  const Field2<real_t>& dispersionTopSnapshot() const { return dzzTop_; }
  ///@}

  /// This rank's scalar-budget contributions of the last step.
  const TransportAudit& audit() const { return audit_; }

  /// This rank's owned surface scalar mass Sum(s dept dx dy).
  real_t ownedSurfaceMass() const;
  /// This rank's owned subsurface scalar mass Sum(s θ V) (the legacy
  /// saltmass_subs monitor observable).
  real_t ownedSubsurfaceMass() const;

 private:
  /// One staged scalar boundary list (member cells on device).
  struct ScalarBcList {
    const BoundaryCondition* bc = nullptr;  ///< the scalar_value condition
    const BoundaryCondition* discharge = nullptr;  ///< paired inflow (surface)
    Kokkos::View<int*, MemSpace> j;
    Kokkos::View<int*, MemSpace> i;
    Kokkos::View<int*, MemSpace> face;  ///< BcFace as int (side conditions)
    long dischargeCells = 0;  ///< paired condition's global member count
  };

  // NOTE: the step-phase helpers below are implementation detail, kept
  // public only because nvcc forbids extended __host__ __device__ lambdas
  // (KOKKOS_LAMBDA) inside private member functions (CUDA C++ programming
  // guide, "Extended Lambda Restrictions"). Treat as private.
 public:
  // ScalarSolver.cpp
  void buildBoundaryLists(const BoundarySet& boundaries);
  void applyInitialConditions(const FrehgConfig& config);
  void snapshotEndOfStep();

  // SurfaceTransport.cpp (scalar_shallowwater, scalar.c:25-298)
  void stepSurface(real_t t, real_t dt, real_t rain, real_t evap);

  // SubsurfaceTransport.cpp (scalar_groundwater, scalar.c:303-498)
  void stepSubsurface(real_t t, real_t dt, real_t dtgLast);
  void enforceSubsurfaceBc(real_t t);

  // Dispersion.cpp (dispersion_tensor, scalar.c:958-1003)
  void updateDispersionTensor();

  // Far-neighbor staging for the superbee stencils at rank interfaces
  // (the legacy guards fall back to upwind at *rank-local* edges, which is
  // rank-count-dependent; Frehg2 stages shifted copies so decomposed runs
  // reproduce the serial golden stencil — docs/theory/transport.md).
  void stageSurfaceFarNeighbors();
  void stageSubsurfaceFarNeighbors();

 private:
  const Grid& grid_;
  HaloExchanger& halo_;

  SurfaceWiring surf_;
  SubsurfaceWiring subs_;
  CouplingWiring cpl_;

  // Configuration extracts.
  bool superbee_ = false;
  real_t difuX_ = 0.0, difuY_ = 0.0;      ///< surface diffusivities
  real_t dispLon_ = 0.0, dispLat_ = 0.0;  ///< dispersivities [m]
  real_t dispMol_ = 0.0;                  ///< molecular diffusivity [m^2/s]
  real_t boundMin_ = 0.0;                 ///< transport.bounds.min
  real_t boundMax_ = 0.0;                 ///< transport.bounds.max (or the
                                          ///< open-bound sentinel)
  bool hasBoundMax_ = false;

  // Surface fields (halo layout).
  Field2<real_t> sSurf_, smSurf_, sSurfKp_, sseepage_;
  Field2<real_t> vsn_, vflux_, fuOld_, fvOld_;
  Field2<real_t> sMinS_, sMaxS_;
  Field2<real_t> sFarXm2_, sFarXp2_, sFarYm2_, sFarYp2_;
  Field2<real_t> dzzTop_;

  // Subsurface fields.
  Field3<real_t> sSubs_, smSubs_;
  Field3<real_t> sMin3_, sMax3_;
  Field3<real_t> dxx_, dyy_, dzz_, dxy_, dxz_, dyz_;
  Field3<real_t> sFarXm3_, sFarXp3_, sFarYm3_, sFarYp3_;
  /// Lower-face conductivity restaged into cell planes (kzF plane k+1) so
  /// the halo exchange can carry it to interface cells (the exchanger packs
  /// nz planes; kzF has nz+1).
  Field3<real_t> kzLower_;
  Field2<int> kTop_;  ///< per-column first active layer (device scratch)

  // Boundary lists.
  std::vector<ScalarBcList> surfaceDirichlet_;  ///< eta-paired (tide) cells
  std::vector<ScalarBcList> surfaceInflow_;     ///< discharge-paired cells
  std::vector<ScalarBcList> sideGhost_;         ///< groundwater_side ghosts

  TransportAudit audit_;
};

}  // namespace frehg::transport

#endif  // FREHG_TRANSPORT_SCALARSOLVER_HPP
