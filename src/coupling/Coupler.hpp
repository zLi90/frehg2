/// \file Coupler.hpp
/// \brief The surface-subsurface coupler (plan §10 P3): synchronization
///        modes, the seepage exchange, and the coupled volume audit.
///
/// The coupler owns the exchange state (plan §4: cross-module data flows
/// through driver-owned coupled state) and reproduces the legacy coupled
/// step (solve.c:25-166):
///
///   beginStep / solveFreeSurface          (surface eta phase)
///   refreshDepth                          (update_depth, shallowwater.c:87)
///   groundwater phase                     (sync: one lockstep gw step;
///                                          subcycled: solve.c:70-94 window)
///   applySeepage                          (subsurface_source,
///                                          shallowwater.c:639-683)
///   updateVelocity                        (surface velocity phase)
///
/// **Sync mode** marches surface and subsurface in lockstep on the common
/// adaptive step: legacy feeds the adapted dtg back into the surface dt
/// (solve_groundwater:193), starting from time.dt (solve.c:37) and clamped
/// to [dt_min, dt_max]. Frehg2 adapts between whole steps so one consistent
/// dt covers both phases — legacy swapped dt mid-step between the eta solve
/// and the velocity update, a §2.1-class inconsistency (amendment A10).
///
/// **Subcycled mode** keeps the fixed surface dt and advances the subsurface
/// with its adaptive dtg while the next substep still fits in the surface
/// step (solve.c:72-75); the subsurface clock may lag the surface clock by
/// up to one dtg, carried across steps (and through restarts).
///
/// **Exchange bookkeeping** (amendment A9): the per-column accumulator
/// carries the exchanged volume per unit area (sum of q dtg / Az over
/// substeps) and the surface applies the whole accumulator each step —
/// identical to legacy's qss dt in sync mode, and conservative in subcycled
/// mode where legacy's rate-based application over-draws when dtg < dt. The
/// infiltration limit debits a per-window available-depth budget seeded from
/// the post-solve depth and credited by every coupled deposit; the legacy
/// porosity factor in the limit (vseep = |q| dtg wcs, subroutines.c:166-168)
/// is dropped per the plan §5.7 unit resolution. Legacy's hold rule is kept:
/// accumulated seepage onto a dry cell is withheld until it exceeds
/// min_depth (shallowwater.c:668-675).

#ifndef FREHG_COUPLING_COUPLER_HPP
#define FREHG_COUPLING_COUPLER_HPP

#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/Types.hpp"
#include "gw/RichardsSolver.hpp"
#include "swe/SurfaceSolver.hpp"

namespace frehg::coupling {

/// This rank's share of the coupled volume budget of one surface step
/// [m^3]. The closure identity asserted by the coupled mass test
/// (plan §5.7, §8.1) is, with S the surface module's own audit terms:
///   d(V_surf + V_gw + accumVolume) = S + gwBoundaryExternal + gw.cplBounce
///                                    - discarded - gw.cplEvap
///                                    - gw.ssStorage + gw.reallocAdjust
///                                    - gw.vloss
/// where gwBoundaryExternal = gw.boundaryIn + gw.cplExchanged + gw.cplVent
/// (the top-face exchange is internal to the coupled system).
struct CouplingAudit {
  gw::GwStepAudit gw;        ///< subsurface audit summed over the window
  real_t applied = 0.0;      ///< net volume applySeepage added to eta (+-)
  real_t discarded = 0.0;    ///< infiltration remainder cut by the eta >= bottom clamp
  real_t surfaceGain = 0.0;  ///< all coupled surface deposits (applied + bounce + vent)
  real_t accumVolume = 0.0;  ///< end-of-step in-transit volume (held seepage)
  long substeps = 0;         ///< subsurface steps taken this surface step
};

/// Orchestrates one coupled time step and owns the exchange fields.
class Coupler {
 public:
  /// Allocate the exchange fields and wire the coupled views into the
  /// groundwater module. Both solvers must outlive the coupler.
  Coupler(const Grid& grid, const FrehgConfig& config, swe::SurfaceSolver& surface,
          gw::RichardsSolver& gw);

  /// Advance the coupled system to time \p t with surface step \p dt
  /// (= nextDt() in sync mode, the fixed time.dt in subcycled mode).
  /// Collective on grid.comm().
  void step(real_t t, real_t dt);

  /// The next surface step: the adapted common step in sync mode (legacy
  /// solve_groundwater:193), the fixed configured dt in subcycled mode.
  real_t nextDt() const { return sync_ ? dtg_ : dt_; }

  /// The current adaptive subsurface step [s] (checkpoint scalar "dtg").
  real_t currentDtg() const { return dtg_; }

  /// Surface dt of the most recently completed step (config time.dt until a
  /// step runs). Checkpointed so a restart refresh can reconstruct the
  /// stage-boundary edge-slot velocities with the dt that produced them.
  real_t lastSurfaceDt() const { return lastSurfaceDt_; }

  /// Subsurface clock lag t - t_gw [s] (checkpoint scalar "tgw_lag";
  /// nonzero only in subcycled mode).
  real_t gwLag(real_t t) const { return sync_ ? 0.0 : (t - tSub_); }

  /// Restore the adaptive state after a checkpoint read. \p t is the
  /// restart time, \p dtg the stored step, \p lag the stored clock lag.
  void restore(real_t t, real_t dtg, real_t lag);

  /// This rank's coupled budget of the last step.
  const CouplingAudit& audit() const { return audit_; }

  /// Surface-applied seepage rate of the last step [m/s] (the legacy qss
  /// observable; output variable "seepage" and checkpoint field "qss").
  const Field2<real_t>& seepageRate() const { return qss_; }

  /// Exchanged-not-yet-applied depth [m] (checkpoint field "seep_accum";
  /// nonzero when the dry-cell hold rule is active).
  const Field2<real_t>& seepAccum() const { return seepAccum_; }

 private:
  // Exchange.cpp kernels.
  void beginWindow();
  void applySeepage(real_t dt);
  void reduceAudit();

  const Grid& grid_;
  swe::SurfaceSolver& surface_;
  gw::RichardsSolver& gw_;

  bool sync_ = true;
  real_t dt_ = 0.0;        ///< configured surface step
  real_t lastSurfaceDt_ = 0.0;  ///< surface dt of the last completed step
  real_t dtg_ = 0.0;       ///< current adaptive subsurface step
  real_t tSub_ = 0.0;      ///< subsurface clock (subcycled mode)
  real_t minDepth_ = 0.0;  ///< surface min_depth (legacy min_dept)

  Field2<real_t> seepAccum_;  ///< exchanged-not-applied depth [m]
  Field2<real_t> avail_;      ///< window infiltration budget [m]
  Field2<real_t> gain_;       ///< coupled surface deposits this step [m]
  Field2<real_t> qss_;        ///< applied seepage rate [m/s]

  CouplingAudit audit_;
};

}  // namespace frehg::coupling

#endif  // FREHG_COUPLING_COUPLER_HPP
