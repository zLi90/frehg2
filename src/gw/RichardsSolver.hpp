/// \file RichardsSolver.hpp
/// \brief The PCA mixed-form Richards module (plan §10 P2).
///
/// One step reproduces the legacy solve_groundwater sequence exactly
/// (groundwater.c:57-199, PCA path — the Newton branch is dropped, plan
/// §3.2):
///
///  1. save h^n and θ^n, evaluate C(h^n), refresh halos and boundary ghosts;
///  2. **predictor** — face conductivities (compute_K_face,
///     groundwater.c:202-335), the 7-point linear head system with storage
///     C(h) + Ss·θ/θs (groundwater_mat_coeff/:503, groundwater_rhs/:572),
///     one global PETSc CG solve (prefix "gw_"), boundary-ghost head
///     enforcement (enforce_head_bc/:751);
///  3. **corrector** — face conductivities again, Darcy fluxes
///     (darcy_flux, subroutines.c:27-195), θ update by flux divergence
///     (update_water_content/:903);
///  4. **post-allocation** — the θ/h consistency restore and the
///     over-saturation redistribution (reallocate_water_content /
///     check_head_gradient / allocate_send, :959-1424, always on per plan
///     Appendix A; surplus handling selected by
///     groundwater.reallocation_surplus — amendment A7);
///  5. final θ clamp to [θr, θs] with volume-loss accounting
///     (groundwater.c:165-186) and moisture-ghost enforcement;
///  6. the adaptive dtg controller (adaptive_time_step/:1651).
///
/// The density/viscosity face-ratio hooks (r_rho, r_visc) are present in
/// every formula (plan §10 P2 "present-but-inactive"); P4 activates them
/// from the transport scalar (attachScalar + Baroclinic.cpp) when
/// groundwater.density_coupling is enabled, and they stay at 1 otherwise.
///
/// Fidelity notes and the fixed legacy defects are tabulated in
/// docs/theory/groundwater.md.

#ifndef FREHG_GW_RICHARDSSOLVER_HPP
#define FREHG_GW_RICHARDSSOLVER_HPP

#include "bc/BoundarySet.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/HaloExchanger.hpp"
#include "core/LinearSystem.hpp"
#include "core/TimeSeries.hpp"
#include "core/Types.hpp"
#include "gw/TerrainMetric.hpp"

#include <memory>
#include <string>
#include <vector>

namespace frehg::gw {

/// \name Baroclinic coupling constants (plan §3.1 item 4/5)
/// The legacy density and viscosity ratios of a saline pore fluid
/// (update_rhovisc, scalar.c:928-930): r_rho = 1 + BETA_RHO s and
/// r_visc = 1 / (1 + BETA_VISC s) with s the salinity [psu]. Named
/// constants with the legacy values (docs/theory/groundwater.md).
///@{
inline constexpr real_t kBaroclinicBetaRho = 0.000744;   ///< r_rho slope
inline constexpr real_t kBaroclinicBetaVisc = 0.0022;    ///< r_visc slope
///@}

/// Boundary-condition codes staged into per-column / per-edge marker fields.
/// The numeric values decode the legacy bctype_GW codes where one exists.
enum class GwBcCode : int {
  NoFlux = 0,          ///< legacy code 0 (default everywhere)
  Head = 1,            ///< legacy code 1 (top htop / bottom hbot)
  Flux = 2,            ///< legacy code 2 (top qtop / bottom qbot / side qyp, qym)
  Gravity = 3,         ///< legacy code 3, bottom only (free drainage)
  HeadHydrostatic = 4  ///< sides only: ghost head = eta_ref - z (enforce_head_bc:769-788)
};

/// One boundary condition's member cells staged for per-step value updates.
struct GwDeviceBcList {
  const BoundaryCondition* bc = nullptr;   ///< source condition (host)
  Kokkos::View<int*, MemSpace> j;          ///< local row indices (halo layout)
  Kokkos::View<int*, MemSpace> i;          ///< local column indices
  Kokkos::View<int*, MemSpace> face;       ///< BcFace as int (side conditions)
  int code = 0;                            ///< GwBcCode as int
};

/// This rank's share of the per-step subsurface volume budget [m^3]. The
/// driver reduces across ranks into the /monitor/gw_mass_audit table; the
/// closed-domain conservation test (plan §10 P2 exit) asserts on it.
struct GwStepAudit {
  real_t boundaryIn = 0.0;     ///< net inflow through domain-boundary faces
  real_t ssStorage = 0.0;      ///< water moved into compressible (Ss) storage
  real_t reallocAdjust = 0.0;  ///< net θ change by the post-allocation step
  /// Volume the post-allocation could not place and discarded (>= 0): the
  /// adjacent-cell surplus in drop mode, lateral split fractions, and
  /// send-walk overflow. Included in reallocAdjust; separated here so gates
  /// and tests can see pure losses.
  real_t reallocDropped = 0.0;
  real_t vloss = 0.0;          ///< volume removed (+) or created (-) by the final clamp
  /// \name Coupled-run terms (zero in groundwater-only runs; plan §10 P3)
  ///@{
  /// Net upward volume through coupled top faces this step [m^3], the
  /// Darcy-flux part (legacy qseepage accumulation, groundwater.c:857).
  real_t cplExchanged = 0.0;
  /// Volume returned to the surface when infiltration meets a full top cell
  /// under a dry surface (legacy groundwater.c:844-851) [m^3, >= 0].
  real_t cplBounce = 0.0;
  /// Volume the post-allocation send walk vented onto the surface (legacy
  /// allocate_send, groundwater.c:1185-1207) [m^3, >= 0].
  real_t cplVent = 0.0;
  /// Evaporative volume removed from the seepage accumulator (legacy
  /// groundwater.c:858-863) [m^3, >= 0]: water that left the subsurface
  /// through the top face but evaporated instead of ponding.
  real_t cplEvap = 0.0;
  ///@}
};

/// Driver-owned coupled state wired in by the coupler (plan §4: cross-module
/// data flows through driver-owned state; §10 P3 Exchange.cpp). All Field2
/// views alias surface-module or coupler storage in the shared halo layout.
/// The groundwater module reads eta/dept for the wet-cell Dirichlet top and
/// mutates eta/dept/avail/gain exactly where legacy did (the bounce-back and
/// the reallocation vent write the surface state mid-substep).
struct GwCoupling {
  bool active = false;         ///< true once the coupler attached the views
  real_t minDepth = 0.0;       ///< surface min_depth [m] (legacy min_dept)
  Field2<real_t> eta;          ///< surface elevation, offset frame (RW)
  Field2<real_t> dept;         ///< surface depth (RW; wet/dry switch input)
  Field2<real_t> seepAccum;    ///< exchanged-not-yet-applied depth [m] (RW)
  Field2<real_t> avail;        ///< remaining infiltrable depth this window [m]
  Field2<real_t> gain;         ///< surface deposits this surface step [m]
};

/// The groundwater module.
class RichardsSolver {
 public:
  /// Allocate fields, stage the soil map, apply initial conditions, rasterize
  /// groundwater boundary conditions, and fix the 7-point COO pattern.
  /// Collective on grid.comm().
  /// \param grid decomposition (kept by reference); buildGlobalIds must
  ///        already have been called with the mesh's ktop.
  /// \param config the loaded configuration.
  /// \param mesh subsurface geometry (kept by reference).
  /// \param boundaries rasterized boundary conditions (kept by reference).
  /// \param halo exchanger the solver registers its fields with.
  RichardsSolver(const Grid& grid, const FrehgConfig& config, const TerrainMetric& mesh,
                 const BoundarySet& boundaries, HaloExchanger& halo);

  /// Wire the coupled surface state in (called once by the coupler before
  /// the first step; plan §10 P3). With an active coupling the top boundary
  /// is coupler-owned: wet columns take the Dirichlet surface depth with the
  /// saturated face conductivity, dry columns behave as a seepage face
  /// (exfiltration allowed, infiltration prohibited — the legacy
  /// bctype_GW[5] = 2 semantics with qtop from any configured flux BC).
  void attachCoupling(const GwCoupling& coupling) { cpl_ = coupling; }

  /// Wire the transport module's scalar state in (called once by the driver
  /// before the first step when transport is active; plan §10 P4). With
  /// \p baroclinic the density/viscosity ratios r_rho = 1 + kBaroclinicBetaRho
  /// s and r_visc = 1/(1 + kBaroclinicBetaVisc s) and their face means are
  /// evaluated from the scalar at the start of every step (legacy
  /// update_rhovisc, scalar.c:921-955, and baroclinic_face,
  /// groundwater.c:338-415); without it the ratios stay at 1 exactly as in
  /// P2/P3. \p sSurf is the surface scalar (the ghost value above coupled
  /// top faces, enforce_scalar_bc, scalar.c:901); pass an empty view in
  /// groundwater-only runs.
  void attachScalar(const Field3<real_t>& sSubs, const Field2<real_t>& sSurf,
                    bool baroclinic) {
    sSubs_ = sSubs;
    sSurfGhost_ = sSurf;
    baroclinic_ = baroclinic;
  }

  /// Advance the subsurface state by \p dtg, evaluating series-valued
  /// boundary conditions at time \p t (the end-of-step time, matching the
  /// legacy call order, solve.c:43-49).
  void step(real_t t, real_t dtg);

  /// The adaptive controller's step size for the next step (legacy dtg after
  /// adaptive_time_step; already clamped and min-reduced across ranks).
  real_t nextDt() const { return dtgNext_; }

  /// \name State access (device views; interiors are authoritative)
  ///@{
  const Field3<real_t>& head() const { return h_; }           ///< pressure head [m]
  const Field3<real_t>& waterContent() const { return wc_; }  ///< volumetric moisture
  ///@}

  /// Per-unit-area x Darcy flux staged for output [m/s]. The internal flux
  /// fields carry the legacy face-area factor [m^3/s]; the accessors stage
  /// q/A into a scratch field (invalidated by the next accessor call).
  const Field3<real_t>& fluxXPerArea();
  /// Per-unit-area y Darcy flux staged for output [m/s].
  const Field3<real_t>& fluxYPerArea();
  /// Per-unit-area z Darcy flux at the face below each cell, positive
  /// upward [m/s] — the convention of the legacy qz output.
  const Field3<real_t>& fluxZPerArea();

  /// Rebuild ghost/derived state from {h, wc} after a restart read.
  /// Collective on grid.comm(). With an attached scalar the transport
  /// module's own refresh must run first (the baroclinic face ratios and
  /// coupled ghost heads read the restored scalar's ghosts).
  void refreshDerivedState();

  /// \name Transport-module state access (plan §10 P4; wired by the driver)
  /// The scalar transport reads the corrector's face fluxes, the pre-step
  /// moisture, and the face conductivities of the completed step.
  ///@{
  /// θ at the start of the last executed step (legacy wcn / Vgn basis).
  const Field3<real_t>& waterContentStart() const { return wcn_; }
  /// Volumetric x face flux [m^3/s] at the x+ face of each cell (slot i = 0
  /// is the west boundary/interface face); positive toward -x.
  const Field3<real_t>& fluxXVolumetric() const { return qx_; }
  /// Volumetric y face flux [m^3/s]; positive toward -y.
  const Field3<real_t>& fluxYVolumetric() const { return qy_; }
  /// Volumetric z face flux [m^3/s], plane k the face above cell k (plane
  /// nz the bottom boundary face); positive upward.
  const Field3<real_t>& fluxZFaceVolumetric() const { return qzF_; }
  /// Face conductivities of the completed step (the legacy Kx/Ky/Kz gates
  /// of the scalar limiter, scalar.c:380-434); layouts match the fluxes.
  const Field3<real_t>& faceConductivityX() const { return kx_; }
  /// \copydoc faceConductivityX
  const Field3<real_t>& faceConductivityY() const { return ky_; }
  /// \copydoc faceConductivityX
  const Field3<real_t>& faceConductivityZFace() const { return kzF_; }
  /// Saturated moisture θs per cell (the dispersion tensor's porosity
  /// factor, scalar.c:965).
  const Field3<real_t>& soilThetaS() const { return wcs_; }
  /// Saturated vertical conductivity per cell (the legacy param->Ksz > 0
  /// gates of the surface scalar limiter, scalar.c:175, :217).
  const Field3<real_t>& soilKsz() const { return ksz_; }
  /// Configured top boundary code / flux value per column (the legacy qtop
  /// quirks of the scalar limiter, scalar.c:452-458 and the top-face
  /// advection rules).
  const Field2<int>& topBcCode() const { return topCode_; }
  /// \copydoc topBcCode
  const Field2<real_t>& topBcValue() const { return topValue_; }
  /// Side boundary codes at the y edges (the sea-side limiter ghost rule,
  /// scalar.c:397-400).
  const Field2<int>& sideBcCodeYp() const { return sideCodeYp_; }
  /// \copydoc sideBcCodeYp
  const Field2<int>& sideBcCodeYm() const { return sideCodeYm_; }
  /// The dtg of the last executed step (the legacy Vgflux step,
  /// groundwater.c:1639).
  real_t lastDtg() const { return dtgCurrent_; }
  ///@}

  /// \name Baroclinic state (diagnostics/tests; held at 1 without an
  /// attached scalar — plan §10 P2 "present-but-inactive", activated in P4)
  ///@{
  const Field3<real_t>& densityRatio() const { return rRho_; }    ///< r_rho per cell
  const Field3<real_t>& viscosityRatio() const { return rVisc_; } ///< r_visc per cell
  const Field3<real_t>& densityRatioFaceY() const { return rRhoYp_; }  ///< y+ face mean
  const Field3<real_t>& densityRatioFaceZ() const { return rRhoZp_; }  ///< z face mean
  ///@}

  /// Convergence data of the last predictor solve.
  const SolveStats& lastSolve() const { return lastSolve_; }

  /// This rank's volume-budget contributions of the last step.
  const GwStepAudit& audit() const { return audit_; }

  /// This rank's owned subsurface water volume Sum(θ dx dy dz) [m^3]
  /// (legacy get_mass, solve.c:372-376).
  real_t ownedVolume() const;

  // NOTE: the step-phase helpers below are implementation detail, kept
  // public only because nvcc forbids extended __host__ __device__ lambdas
  // (KOKKOS_LAMBDA) inside private member functions (CUDA C++ programming
  // guide, "Extended Lambda Restrictions"). Treat as private.
 public:
  // Initialization (RichardsSolver.cpp).
  void stageSoil(const FrehgConfig& config);
  void applyInitialConditions(const FrehgConfig& config);
  void buildBoundaryLists(const BoundarySet& boundaries, const FrehgConfig& config);
  void buildCooPattern();
  void updateBoundaryValues(real_t t);

  // Predictor (Predictor.cpp).
  void computeFaceConductivity();
  void assembleSystem(real_t dtg);
  void fillAndSolve();
  void enforceHeadBc();

  // Baroclinic activation (Baroclinic.cpp): r_rho/r_visc from the attached
  // scalar and their face means (legacy update_rhovisc + baroclinic_face).
  void updateBaroclinicFaces();

  // Corrector (Corrector.cpp).
  void classifyCoupledTop(real_t dtg);
  void computeFluxes(real_t dtg);
  void applyCoupledTopBookkeeping(real_t dtg);
  void updateWaterContent(real_t dtg);
  void enforceMoistureBc();
  void finalizeWaterContent();
  void accumulateBoundaryFlux(real_t dtg);

  // Post-allocation (Reallocate.cpp).
  void reallocateWaterContent();

  // Adaptive stepping (AdaptiveStep.cpp).
  void adaptTimeStep(real_t dtg);

 private:
  const Grid& grid_;
  const TerrainMetric& mesh_;
  HaloExchanger& halo_;

  /// Coupled surface state (inactive views in groundwater-only runs).
  GwCoupling cpl_;

  /// Transport-module scalar views (empty until attachScalar; plan §10 P4).
  Field3<real_t> sSubs_;      ///< subsurface scalar with ghosts
  Field2<real_t> sSurfGhost_; ///< surface scalar (coupled top ghost value)
  bool baroclinic_ = false;   ///< groundwater.density_coupling.enabled

  // Configuration extracts.
  real_t ss_ = 0.0;              ///< specific storage [1/m]
  bool useFull3d_ = true;
  bool surplusRedistribute_ = false;  ///< groundwater.reallocation_surplus (amendment A7)
  real_t dtMin_ = 0.0, dtMax_ = 0.0;
  real_t dqGrow_ = 0.01, dqShrink_ = 0.02, courantMax_ = 2.0;
  real_t dtgNext_ = 0.0;

  // State fields (nyLocal+2, nxLocal+2, nz), (j, i, k) with one-cell halos
  // in j and i. Face fields: Kx/qx at the x+ face of the cell (slot i = 0 is
  // the west boundary/interface face); KzF/qzF sized nz+1 with plane k being
  // the face *above* cell k and plane nz the bottom boundary face.
  Field3<real_t> h_, hn_, wc_, wcn_, ch_;
  Field3<real_t> ksx_, ksy_, ksz_, vga_, vgn_, wcs_, wcr_, aev_;
  Field3<real_t> kx_, ky_, kzF_;
  Field3<real_t> qx_, qy_, qzF_;
  Field3<real_t> rRho_, rRhoXp_, rRhoYp_, rRhoZp_;
  Field3<real_t> rVisc_, rViscXp_, rViscYp_, rViscZp_;
  Field3<real_t> vloss_;
  Field3<real_t> room_, sendUp_, sendDown_;  ///< post-allocation scratch
  Field3<real_t> qOut_;  ///< scratch for the per-area flux accessors
  real_t dtgCurrent_ = 0.0;  ///< the step's dtg (post-allocation flux bookkeeping)

  // Boundary-condition marker fields (2D, halo layout). Codes are GwBcCode
  // as int; values are the current prescribed value per column/edge cell.
  Field2<int> topCode_, botCode_;
  Field2<real_t> topValue_, botValue_;
  /// Coupled top-exchange mode per column, refreshed each substep
  /// (amendment A13; SERGHEI GwBC.h:277-288): 0 = dry seepage face,
  /// 1 = wet Dirichlet (capacity-limited), 2 = wet supply-limited flux.
  Field2<int> cplMode_;
  Field2<int> sideCodeXm_, sideCodeXp_, sideCodeYm_, sideCodeYp_;
  Field2<real_t> sideValueXm_, sideValueXp_, sideValueYm_, sideValueYp_;
  std::vector<GwDeviceBcList> bcLists_;

  // Linear system.
  std::unique_ptr<LinearSystem> system_;
  Kokkos::View<PetscInt*, MemSpace> cooRows_, cooCols_;
  Kokkos::View<real_t*, MemSpace> cooValues_;
  Kokkos::View<real_t*, MemSpace> rhsVec_, solVec_;
  SolveStats lastSolve_;

  GwStepAudit audit_;
};

}  // namespace frehg::gw

#endif  // FREHG_GW_RICHARDSSOLVER_HPP
