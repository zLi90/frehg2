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

/// Stage-3 T2 (plan §6 W2.3): the advection stencils with the external
/// reservoir state. Ext is a compile-time switch so the Ext = false
/// instantiation folds every selector away and keeps the v2.0 expression
/// trees verbatim (the default-bitwise gate verifies the result
/// empirically). Under Ext = true:
///  - a stencil value read ACROSS a reservoir face is the reservoir's
///    velocity (uRes/vRes at that face);
///  - SELF-advection across a reservoir face takes the conservative
///    d(u^2/2) form, so a steady entrance from a reservoir at u_ext costs
///    exactly (u^2 - u_ext^2)/2g of head (the 2026-10 diagnosis's 0.5
///    weight, generalized to a moving reservoir);
///  - outflow and every non-reservoir face keep the legacy upwind form.
template <bool Ext>
struct ExtStencil {
  Field2<real_t> uu;    ///< face-centered x velocity
  Field2<real_t> vv;    ///< face-centered y velocity
  Field2<real_t> uRes;  ///< reservoir u at x reservoir faces
  Field2<real_t> vRes;  ///< reservoir v at y reservoir faces
  Field2<real_t> resX;  ///< x reservoir-face mask (nonzero = reservoir)
  Field2<real_t> resY;  ///< y reservoir-face mask (nonzero = reservoir)

  /// The x-velocity the advection stencil reads at face (j, i).
  KOKKOS_INLINE_FUNCTION real_t uA(int j, int i) const {
    if constexpr (Ext) {
      return resX(j, i) != 0.0 ? uRes(j, i) : uu(j, i);
    } else {
      return uu(j, i);
    }
  }
  /// The y-velocity the advection stencil reads at face (j, i).
  KOKKOS_INLINE_FUNCTION real_t vA(int j, int i) const {
    if constexpr (Ext) {
      return resY(j, i) != 0.0 ? vRes(j, i) : vv(j, i);
    } else {
      return vv(j, i);
    }
  }
};

void SurfaceSolver::momentumSource() {
  if (hasExternal_) {
    momentumSourceImpl<true>();
  } else {
    momentumSourceImpl<false>();
  }
}

template <bool Ext>
void SurfaceSolver::momentumSourceImpl() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t dt = dt_;
  const real_t viscx = viscX_;
  const real_t viscy = viscY_;
  Field2<real_t> uu = uu_, vv = vv_, uy = uy_, vx = vx_;
  const ExtStencil<Ext> st{uu_, vv_, uResX_, vResY_, resFaceX_, resFaceY_};
  Field2<real_t> uRes = uResX_, vRes = vResY_, resX = resFaceX_, resY = resFaceY_;
  Field2<real_t> cflx = cflx_, cfly = cfly_;
  Field2<real_t> Vsx = Vsx_, Vsy = Vsy_, Asx = Asx_, Asy = Asy_, Aszx = Aszx_, Aszy = Aszy_;
  Field2<real_t> CDx = CDx_, CDy = CDy_;
  Field2<real_t> Ex = Ex_, Ey = Ey_, Dx = Dx_, Dy = Dy_;
  Field2<real_t> deptx = deptx_, depty = depty_;

  const bool wind = windCfg_.enabled;
  // Per-step wind sample (v2 Q6, WindForcing): speed, direction in
  // radians from +x (the legacy compass arithmetic preserved on the
  // speed/direction path; atan2 on the component path), and Cd under the
  // configured law (the legacy constant Cw by default).
  const real_t windDrag = windCd_;
  const real_t windAttenuation = windCfg_.attenuationDepth;
  const real_t hD = thinLayerDepth_;
  const real_t omega = windOmega_;
  const real_t windSpeed = windSpeed_;

  Kokkos::parallel_for(
      "swe_momentum_source",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        // nvcc: an extended __host__ __device__ lambda may not FIRST capture
        // a variable inside an if-constexpr branch (CUDA C++ guide, "Extended
        // Lambda Restrictions"), so the Ext branch's captures are named
        // here first. Discarded expressions: no code, no arithmetic.
        (void)st;
        (void)uRes;
        (void)vRes;
        (void)resX;
        (void)resY;
        (void)dx;
        (void)dy;
        const real_t u = uu(j, i);
        const real_t v = vv(j, i);
        const real_t vxc = vx(j, i);
        const real_t uyc = uy(j, i);

        // Upwind advection (shallowwater.c:158-165). Under Ext the
        // self-advection across a reservoir face takes the conservative
        // d(u^2/2) form and the cross/neighbor reads take the reservoir
        // state (the ExtStencil header comment); Ext = false folds to the
        // legacy expressions verbatim.
        real_t advSelfX;
        real_t advSelfY;
        if constexpr (Ext) {
          if (u > 0.0 && resX(j, i - 1) != 0.0) {
            const real_t uE = uRes(j, i - 1);
            advSelfX = (u * u - uE * uE) / (2.0 * dx);
          } else if (u < 0.0 && resX(j, i + 1) != 0.0) {
            const real_t uE = uRes(j, i + 1);
            advSelfX = (uE * uE - u * u) / (2.0 * dx);
          } else {
            advSelfX = (0.5 / dx) * ((u + Kokkos::fabs(u)) * (u - st.uA(j, i - 1)) +
                                     (u - Kokkos::fabs(u)) * (st.uA(j, i + 1) - u));
          }
          if (v > 0.0 && resY(j - 1, i) != 0.0) {
            const real_t vE = vRes(j - 1, i);
            advSelfY = (v * v - vE * vE) / (2.0 * dy);
          } else if (v < 0.0 && resY(j + 1, i) != 0.0) {
            const real_t vE = vRes(j + 1, i);
            advSelfY = (vE * vE - v * v) / (2.0 * dy);
          } else {
            advSelfY = (0.5 / dy) * ((v + Kokkos::fabs(v)) * (v - st.vA(j - 1, i)) +
                                     (v - Kokkos::fabs(v)) * (st.vA(j + 1, i) - v));
          }
        } else {
          advSelfX = (0.5 / dx) * ((u + Kokkos::fabs(u)) * (u - uu(j, i - 1)) +
                                   (u - Kokkos::fabs(u)) * (uu(j, i + 1) - u));
          advSelfY = (0.5 / dy) * ((v + Kokkos::fabs(v)) * (v - vv(j - 1, i)) +
                                   (v - Kokkos::fabs(v)) * (vv(j + 1, i) - v));
        }
        real_t advX =
            advSelfX +
            (0.5 / dy) * ((vxc + Kokkos::fabs(vxc)) * (u - st.uA(j - 1, i)) +
                          (vxc - Kokkos::fabs(vxc)) * (st.uA(j + 1, i) - u));
        real_t advY =
            (0.5 / dx) * ((uyc + Kokkos::fabs(uyc)) * (v - st.vA(j, i - 1)) +
                          (uyc - Kokkos::fabs(uyc)) * (st.vA(j, i + 1) - v)) +
            advSelfY;
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
  if (hasExternal_) {
    interpolateVelocityImpl<true>();
  } else {
    interpolateVelocityImpl<false>();
  }
}

template <bool Ext>
void SurfaceSolver::interpolateVelocityImpl() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  Field2<real_t> uu = uu_, vv = vv_, uy = uy_, vx = vx_;
  const ExtStencil<Ext> st{uu_, vv_, uResX_, vResY_, resFaceX_, resFaceY_};

  // uy/vx feed only momentumSource; under Ext the reservoir faces enter
  // the four-point interpolation at the reservoir state.
  Kokkos::parallel_for(
      "swe_interp_velocity",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        uy(j, i) = 0.25 * (st.uA(j, i) + st.uA(j, i - 1) + st.uA(j + 1, i) +
                           st.uA(j + 1, i - 1));
        vx(j, i) = 0.25 * (st.vA(j, i) + st.vA(j - 1, i) + st.vA(j, i + 1) +
                           st.vA(j - 1, i + 1));
      });
}

}  // namespace frehg::swe
