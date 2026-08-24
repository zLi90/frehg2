/// \file Momentum.cpp
/// \brief Momentum source terms, drag coefficient, velocity update, and the
///        uy/vx interpolation.
///
/// Provenance: momentum_source (shallowwater.c:151-256), wind_source
/// (shallowwater.c:258-292), update_drag_coef (shallowwater.c:977-994),
/// update_velocity (shallowwater.c:719-835), interp_velocity
/// (shallowwater.c:962-974). The interpolation uses the intended four-point
/// stencil; the legacy flat-index arithmetic wrapped to the far end of the
/// row at i = 0 and read out of bounds at j = 0 (both index-arithmetic
/// defects of the class listed in plan §2.1, fixed on port and recorded in
/// docs/theory/surface-water.md).

#include "swe/SurfaceSolver.hpp"
#include "swe/SweFormulas.hpp"

namespace frehg::swe {

void SurfaceSolver::momentumSource() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t dt = dt_;
  const real_t viscx = viscX_;
  const real_t viscy = viscY_;
  Field2<real_t> uu = uu_, vv = vv_, uy = uy_, vx = vx_;
  Field2<real_t> cflx = cflx_, cfly = cfly_;
  Field2<real_t> Vsx = Vsx_, Vsy = Vsy_, Asx = Asx_, Asy = Asy_, Aszx = Aszx_, Aszy = Aszy_;
  Field2<real_t> CDx = CDx_, CDy = CDy_;
  Field2<real_t> Ex = Ex_, Ey = Ey_, Dx = Dx_, Dy = Dy_;
  Field2<real_t> deptx = deptx_, depty = depty_;

  const bool wind = windCfg_.enabled;
  const real_t windDrag = windCfg_.cd;
  const real_t windAttenuation = windCfg_.attenuationDepth;
  const real_t hD = thinLayerDepth_;
  // Wind direction measured from north plus the grid-to-north rotation,
  // converted to radians from +x (legacy wind_source, shallowwater.c:263-265
  // with pi truncated to the legacy literal 3.1415926).
  const real_t legacyPi = 3.1415926;
  const real_t omega = (windDirection_ + windCfg_.northAngle) * legacyPi / 180.0;
  const real_t windSpeed = windSpeed_;

  Kokkos::parallel_for(
      "swe_momentum_source",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const real_t u = uu(j, i);
        const real_t v = vv(j, i);
        const real_t vxc = vx(j, i);
        const real_t uyc = uy(j, i);

        // Upwind advection (shallowwater.c:158-165).
        real_t advX =
            (0.5 / dx) * ((u + Kokkos::fabs(u)) * (u - uu(j, i - 1)) +
                          (u - Kokkos::fabs(u)) * (uu(j, i + 1) - u)) +
            (0.5 / dy) * ((vxc + Kokkos::fabs(vxc)) * (u - uu(j - 1, i)) +
                          (vxc - Kokkos::fabs(vxc)) * (uu(j + 1, i) - u));
        real_t advY =
            (0.5 / dx) * ((uyc + Kokkos::fabs(uyc)) * (v - vv(j, i - 1)) +
                          (uyc - Kokkos::fabs(uyc)) * (vv(j, i + 1) - v)) +
            (0.5 / dy) * ((v + Kokkos::fabs(v)) * (v - vv(j - 1, i)) +
                          (v - Kokkos::fabs(v)) * (vv(j + 1, i) - v));
        advX = cflDampedAdvection(advX, u, cflx(j, i));
        advY = cflDampedAdvection(advY, v, cfly(j, i));

        // Central eddy viscosity with the legacy face-area asymmetry
        // (shallowwater.c:170-190): the x-part of difX uses Asx(j,i) on both
        // faces, the y-part of difY uses Asy(j,i) on both faces.
        real_t difX = 0.0;
        real_t difY = 0.0;
        if (Vsx(j, i) > 0.0) {
          difX = (viscx / Vsx(j, i) / dx) * (Asx(j, i) * (uu(j, i + 1) - u) -
                                             Asx(j, i) * (u - uu(j, i - 1))) +
                 (viscy / Vsx(j, i) / dy) * (Asy(j, i) * (uu(j + 1, i) - u) -
                                             Asy(j - 1, i) * (u - uu(j - 1, i)));
        }
        if (Vsy(j, i) > 0.0) {
          difY = (viscx / Vsy(j, i) / dx) * (Asx(j, i) * (vv(j, i + 1) - v) -
                                             Asx(j, i - 1) * (v - vv(j, i - 1))) +
                 (viscy / Vsy(j, i) / dy) * (Asy(j, i) * (vv(j + 1, i) - v) -
                                             Asy(j, i) * (v - vv(j - 1, i)));
        }

        // Point-implicit drag (shallowwater.c:192-202).
        const real_t velx = Kokkos::sqrt(u * u + vxc * vxc);
        const real_t vely = Kokkos::sqrt(uyc * uyc + v * v);
        const real_t facdx = (Vsx(j, i) > 0.0) ? Aszx(j, i) / Vsx(j, i) : 0.0;
        const real_t facdy = (Vsy(j, i) > 0.0) ? Aszy(j, i) / Vsy(j, i) : 0.0;
        Dx(j, i) = pointImplicitFactor(dt, CDx(j, i), velx, facdx);
        Dy(j, i) = pointImplicitFactor(dt, CDy(j, i), vely, facdy);

        // Momentum source; wind enters before the drag factor is applied
        // (shallowwater.c:204-210).
        real_t ex = u + dt * (difX - advX);
        real_t ey = v + dt * (difY - advY);
        if (wind) {
          const real_t tau = windStress(windDrag, windSpeed, u, v, omega);
          const real_t tauX = windStressAttenuated(tau, deptx(j, i), hD, windAttenuation);
          const real_t tauY = windStressAttenuated(tau, depty(j, i), hD, windAttenuation);
          if (deptx(j, i) > 0.0) {
            ex += dt * tauX * Kokkos::cos(omega) / (deptx(j, i) * kWaterDensity);
          }
          if (depty(j, i) > 0.0) {
            ey += dt * tauY * Kokkos::sin(omega) / (depty(j, i) * kWaterDensity);
          }
        }
        Ex(j, i) = ex * Dx(j, i);
        Ey(j, i) = ey * Dy(j, i);
      });
}

void SurfaceSolver::updateDragCoef() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t cellArea = grid_.dx() * grid_.dy();
  const real_t gravity = gravity_;
  const real_t hD = thinLayerDepth_;
  const real_t minDepth = minDepth_;
  const bool chezy = chezy_;
  Field2<real_t> Vs = Vs_, CDx = CDx_, CDy = CDy_, coefField = frictionCoef_;
  Field2<real_t> deptx = deptx_, depty = depty_;

  // Legacy update_drag_coef only writes where Vs > 0 (shallowwater.c:984),
  // leaving dry cells with a stale (initially zero) coefficient: a face
  // between a dry cell and a wet neighbor then carries water with *no*
  // drag, and the point-implicit update accumulates gravity unopposed
  // (b5's dry outlet-adjacent faces reached km/s; no earlier benchmark
  // keeps a persistently dry cell on a live face). Dry cells evaluate the
  // same law at their deepest adjacent face depth — the water actually
  // flowing across them (a bed-level outfall reservoir sees the upstream
  // depth; a hairline film face sees the film) — floored at min_depth.
  // Wet cells are bit-identical (wet means depth > min_depth).
  // Amendment A13.
  Kokkos::parallel_for(
      "swe_drag_coef",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        real_t depth = Vs(j, i) / cellArea;
        if (Vs(j, i) <= 0.0) {
          depth = Kokkos::fmax(
              Kokkos::fmax(deptx(j, i), deptx(j, i - 1)),
              Kokkos::fmax(depty(j, i), depty(j - 1, i)));
          depth = Kokkos::fmax(depth, minDepth);
        }
        const real_t cd = chezy ? chezyDrag(gravity, coefField(j, i))
                                : manningDrag(gravity, coefField(j, i), depth, hD);
        CDx(j, i) = cd;
        CDy(j, i) = cd;
      });
}

void SurfaceSolver::updateVelocityField() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t coef = gravity_ * dt_;
  Field2<real_t> uu = uu_, vv = vv_, Ex = Ex_, Ey = Ey_, Dx = Dx_, Dy = Dy_;
  Field2<real_t> Vsx = Vsx_, Vsy = Vsy_, Asx = Asx_, Asy = Asy_, eta = eta_;

  // The stored velocity applies the drag factor to the full expression even
  // though Ex already carries one factor — the legacy inconsistency is part
  // of the preserved scheme (update_velocity, shallowwater.c:745-746).
  Kokkos::parallel_for(
      "swe_update_velocity",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const real_t effhx = (Vsx(j, i) > 0.0) ? Asx(j, i) / Vsx(j, i) : 0.0;
        const real_t effhy = (Vsy(j, i) > 0.0) ? Asy(j, i) / Vsy(j, i) : 0.0;
        uu(j, i) = (Ex(j, i) - coef * effhx * (eta(j, i + 1) - eta(j, i))) * Dx(j, i);
        vv(j, i) = (Ey(j, i) - coef * effhy * (eta(j + 1, i) - eta(j, i))) * Dy(j, i);
      });
}

void SurfaceSolver::interpolateVelocity() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  Field2<real_t> uu = uu_, vv = vv_, uy = uy_, vx = vx_;

  Kokkos::parallel_for(
      "swe_interp_velocity",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        uy(j, i) = 0.25 * (uu(j, i) + uu(j, i - 1) + uu(j + 1, i) + uu(j + 1, i - 1));
        vx(j, i) = 0.25 * (vv(j, i) + vv(j - 1, i) + vv(j, i + 1) + vv(j - 1, i + 1));
      });
}

}  // namespace frehg::swe
