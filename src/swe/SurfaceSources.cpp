/// \file SurfaceSources.cpp
/// \brief Rainfall and evaporation applied directly to the free surface.
///
/// Provenance: evaprain (shallowwater.c:577-636). Legacy hardcoded a skip of
/// the last global row when adding rain (shallowwater.c:596, a b1-specific
/// outlet-row rule); Frehg2 expresses the same behavior as the configured
/// surface_water.rainfall.exclude region (schema amendment, plan §6), so b1
/// stays golden-faithful while cases without the quirk rain uniformly.
/// Coupled runs rain on every non-excluded cell too: the legacy coupled
/// branch (shallowwater.c:602-612) rained only on already-wet cells,
/// silently discarding rain over dry land — a conservation defect no golden
/// pins (b6 is rainless; b5 was never a legacy-Frehg run) that would make
/// the b5 rain gate unreachable (amendment A11).

#include "swe/SurfaceSolver.hpp"

namespace frehg::swe {

void SurfaceSolver::evapRain() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t dt = dt_;
  const real_t rain = rain_;
  const real_t evap = evap_;
  const real_t minDepth = minDepth_;
  const real_t cellArea = grid_.dx() * grid_.dy();
  Field2<real_t> eta = eta_, bottom = bottom_, mask = rainMask_, emask = evapMask_;

  // Rain lands on every non-excluded cell, wet or dry
  // (shallowwater.c:583-601).
  real_t rainAdded = 0.0;
  const bool et = etOpenWater_;
  Field2<real_t> etRemoved = et ? evapRemoved_ : Field2<real_t>();
  Field2<real_t> etRain = et ? rainAdded_ : Field2<real_t>();
  Field2<real_t> wetFrac = et ? wetFrac_ : Field2<real_t>();
  Field2<real_t> wetClock = et ? wetClock_ : Field2<real_t>();
  Field2<real_t> wetClockStart = et ? wetClockStepStart_ : Field2<real_t>();
  if (et) {
    // Stage-3 T5 (W5.4): the coupler's dry-fraction clock reads the
    // step-start value and this step's wetFrac to reconstruct the wet
    // time at any instant inside the step.
    Kokkos::deep_copy(wetClockStart, wetClock);
  }
  Kokkos::parallel_reduce(
      "swe_rain",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sum) {
        eta(j, i) += rain * dt * mask(j, i);
        if (et) {
          etRain(j, i) = rain * dt * mask(j, i);
        }
        sum += rain * dt * mask(j, i) * cellArea;
      },
      rainAdded);
  audit_.rainVolume += rainAdded;

  real_t evapRemoved = 0.0;
  const bool canopy = et && hasCanopy_;
  Field2<real_t> tau = canopy ? canopyTau_ : Field2<real_t>();
  if (et) {
    // Stage-3 T5 (W5.3): the ET open-water channel. Per cell, after the
    // rain: a wet cell loses min(E_w dt, depth) (E_w < 0 condenses onto
    // water unconditionally); a dry cell loses nothing — the soil channel
    // owns it (W5.4). The wet fraction min(1, depth/(E_w dt)) is the time
    // share the surface channel covered; the wet clock integrates it, and
    // the soil channel's dry fraction is its complement on each substep —
    // complementarity on a time basis by construction (plan §9).
    Kokkos::parallel_reduce(
        "swe_et_open_water",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i, real_t& sum) {
          const real_t depth = eta(j, i) - bottom(j, i);
          // Stage-3 T6 (W6.3): the canopy transmission scales the RATE
          // (the wet-coverage time is physical: a shaded wet cell is
          // still wet, so the clock integrates against the shaded rate).
          const real_t eEff = canopy ? tau(j, i) * evap : evap;
          real_t removed = 0.0;
          real_t frac = 0.0;
          if (depth > 0.0) {
            if (eEff > 0.0) {
              removed = Kokkos::fmin(eEff * dt, depth);
              frac = Kokkos::fmin(1.0, depth / (eEff * dt));
            } else {
              removed = eEff * dt;  // condensation deposits in full
              frac = 1.0;
            }
            eta(j, i) -= removed;
          }
          etRemoved(j, i) = removed;
          wetFrac(j, i) = frac;
          wetClock(j, i) += frac * dt;
          sum += removed * cellArea;
        },
        evapRemoved);
  } else if (evapMode_ == SurfaceWaterConfig::EvapMode::Bulk) {
    // v2 Q4 bulk mode (plan §3.2): wet cells only, at most the available
    // depth per cell — no volume creation, so the audited value IS the
    // actual removal. A negative rate (condensation, V2-A13) deposits.
    Kokkos::parallel_reduce(
        "swe_evap_bulk",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i, real_t& sum) {
          const real_t depth = eta(j, i) - bottom(j, i);
          if (depth > 0.0) {
            const real_t removed = Kokkos::fmin(evap * dt, depth);
            eta(j, i) -= removed;
            sum += removed * cellArea;
          }
        },
        evapRemoved);
  } else {
    // Prescribed mode: subtracted unconditionally over the non-excluded
    // cells; the clamp below absorbs over-drying (shallowwater.c:614-626;
    // no exclude region -> the mask is 1.0 everywhere and the arithmetic
    // is bitwise the legacy form).
    Kokkos::parallel_reduce(
        "swe_evap",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i, real_t& sum) {
          eta(j, i) -= evap * dt * emask(j, i);
          sum += evap * dt * emask(j, i) * cellArea;
        },
        evapRemoved);
  }
  audit_.evapVolume += evapRemoved;

  // Negative and sub-threshold depths are dried (shallowwater.c:628-635).
  // v2 Q4: the clamp's signed volume change is measured into the clamp
  // audit (positive = created) instead of silently entering the closure
  // residual — the g4(d) audited-shortfall requirement. The eta arithmetic
  // is unchanged.
  real_t clamped = 0.0;
  Kokkos::parallel_reduce(
      "swe_evaprain_clamp",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& acc) {
        const real_t diff = eta(j, i) - bottom(j, i);
        if (diff < minDepth) {
          acc += (bottom(j, i) - eta(j, i)) * cellArea;
          eta(j, i) = bottom(j, i);
        }
      },
      clamped);
  audit_.clampVolume += clamped;
}

}  // namespace frehg::swe
