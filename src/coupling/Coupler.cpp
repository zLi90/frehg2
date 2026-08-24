/// \file Coupler.cpp
/// \brief Coupled step orchestration (legacy solve.c:25-166; the coupled
///        sequence and synchronization modes are documented in Coupler.hpp).

#include "coupling/Coupler.hpp"

#include "core/Logger.hpp"
#include "core/Timer.hpp"

namespace frehg::coupling {

Coupler::Coupler(const Grid& grid, const FrehgConfig& config, swe::SurfaceSolver& surface,
                 gw::RichardsSolver& gw)
    : grid_(grid), surface_(surface), gw_(gw) {
  sync_ = (config.coupling.mode == CouplingConfig::Mode::Sync);
  dt_ = config.time.dt;
  // Sync starts the common adaptive step from the configured dt (legacy
  // solve.c:37 sets dtg = dt; dt_init is a subcycled-mode knob there too);
  // subcycled starts the subsurface at dt_init.
  dtg_ = sync_ ? dt_ : config.groundwater.timestep.dtInit;
  lastSurfaceDt_ = dt_;
  tSub_ = config.time.tStart;
  minDepth_ = config.surfaceWater.minDepth;

  const std::size_t ny2 = static_cast<std::size_t>(grid_.nyLocal()) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(grid_.nxLocal()) + 2;
  seepAccum_ = Field2<real_t>("cpl_seep_accum", ny2, nx2);
  avail_ = Field2<real_t>("cpl_avail", ny2, nx2);
  gain_ = Field2<real_t>("cpl_gain", ny2, nx2);
  qss_ = Field2<real_t>("cpl_qss", ny2, nx2);
  Kokkos::deep_copy(seepAccum_, 0.0);
  Kokkos::deep_copy(avail_, 0.0);
  Kokkos::deep_copy(gain_, 0.0);
  Kokkos::deep_copy(qss_, 0.0);

  gw::GwCoupling coupling;
  coupling.active = true;
  coupling.minDepth = minDepth_;
  coupling.eta = surface_.eta();
  coupling.dept = surface_.depth();
  coupling.seepAccum = seepAccum_;
  coupling.avail = avail_;
  coupling.gain = gain_;
  gw_.attachCoupling(coupling);
}

void Coupler::restore(real_t t, real_t dtg, real_t lag) {
  dtg_ = dtg;
  tSub_ = t - lag;
}

void Coupler::step(real_t t, real_t dt) {
  Timer::Scoped timer("coupled_step");
  lastSurfaceDt_ = dt;
  surface_.setTimeStep(dt);
  surface_.beginStep(t);
  surface_.solveFreeSurface();
  // The groundwater phase reads the post-solve depth (legacy update_depth at
  // the end of solve_shallowwater, shallowwater.c:87).
  surface_.refreshDepth();

  beginWindow();
  audit_ = CouplingAudit{};
  const auto takeGwStep = [&](real_t tt, real_t dtg) {
    gw_.step(tt, dtg);
    dtg_ = gw_.nextDt();
    const gw::GwStepAudit& a = gw_.audit();
    audit_.gw.boundaryIn += a.boundaryIn;
    audit_.gw.ssStorage += a.ssStorage;
    audit_.gw.reallocAdjust += a.reallocAdjust;
    audit_.gw.reallocDropped += a.reallocDropped;
    audit_.gw.vloss += a.vloss;
    audit_.gw.cplExchanged += a.cplExchanged;
    audit_.gw.cplBounce += a.cplBounce;
    audit_.gw.cplVent += a.cplVent;
    audit_.gw.cplEvap += a.cplEvap;
    ++audit_.substeps;
  };

  if (sync_) {
    // Lockstep: one subsurface step of the common dt (legacy
    // sync_coupling = 1 path, solve.c:68-69).
    tSub_ = t;
    takeGwStep(t, dt);
  } else {
    // Subcycled window (solve.c:72-75): advance while a whole dtg fits.
    // Every rank iterates identically — dtg is min-reduced inside the
    // subsurface step — so the loop is collective-safe.
    while (tSub_ + dtg_ <= t + 1.0e-9) {
      tSub_ += dtg_;
      takeGwStep(tSub_, dtg_);
    }
  }

  applySeepage(dt);
  surface_.updateVelocity();
  reduceAudit();
}

}  // namespace frehg::coupling
