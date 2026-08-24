/// \file Corrector.cpp
/// \brief The PCA corrector: Darcy face fluxes, the mass-conservative θ
///        update, moisture ghosts, and the final clamp.
///
/// Provenance: darcy_flux (subroutines.c:27-195), groundwater_flux
/// (groundwater.c:816-869), update_water_content (:903-930),
/// enforce_moisture_bc (:933-956), and the final θ check
/// (solve_groundwater:163-189).
///
/// The corrector re-evaluates its face conductivities from the *predicted*
/// head exactly as darcy_flux does — including the legacy quirks that the
/// x/y faces use the minus cell's saturated conductivity for both flanks
/// while z faces use the plus (lower) cell's, and that each flank's own vG
/// parameters enter its closure. Interface faces evaluate two-sided from
/// exchanged neighbor state (the legacy one-sided interface treatment was
/// rank-dependent; single-rank values are identical). The flux fields carry
/// the legacy face-area factor [m^3/s]; positive x/y flux points in the
/// -axis direction and positive z flux upward.

#include "core/Logger.hpp"
#include "gw/RichardsSolver.hpp"
#include "gw/VanGenuchten.hpp"

namespace frehg::gw {

namespace {

template <class V3>
KOKKOS_INLINE_FUNCTION VgSoil soilAt(const V3& vga, const V3& vgn, const V3& wcs, const V3& wcr,
                                     const V3& aev, int j, int i, int k) {
  return VgSoil{vga(j, i, k), vgn(j, i, k), wcs(j, i, k), wcr(j, i, k), aev(j, i, k)};
}

}  // namespace

void RichardsSolver::classifyCoupledTop(real_t dtg) {
  // Coupled top-exchange classification (amendment A13; SERGHEI
  // GwBC.h:277-288, the reference implementation's evolution of the legacy
  // exchange): a wet column whose saturated-Darcy demand the available
  // surface water can fund keeps the legacy Dirichlet-head top (mode 1); a
  // wet column it cannot fund becomes supply-limited (mode 2) — the
  // predictor and corrector both take exactly the remaining depth as a top
  // flux, so the predicted head stays physical and the consistency restore
  // cannot manufacture the shortfall (legacy applied the Dirichlet
  // unconditionally, which races the infiltration front far beyond the
  // supply; no legacy benchmark reaches this regime — b6 is always
  // capacity-limited). Dry columns are seepage faces (mode 0).
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();

  Field3<real_t> h = h_, ksz = ksz_;
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> mode = cplMode_;
  Field2<real_t> dept = cpl_.dept, avail = cpl_.avail;

  Kokkos::parallel_for(
      "gw_classify_coupled_top",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        int kTop = -1;
        for (int k = 0; k < nz; ++k) {
          if (gid(j, i, k) >= 0) {
            kTop = k;
            break;
          }
        }
        if (kTop < 0 || dept(j, i) <= 0.0) {
          mode(j, i) = 0;
          return;
        }
        // Saturated-conductivity demand estimate (SERGHEI GwBC.h:280;
        // identical to the legacy top-face Darcy form with kface = Ksz).
        const real_t ks = ksz(j, i, kTop);
        const real_t qEst =
            2.0 * ks * (h(j, i, kTop) - dept(j, i)) / dz3d(j, i, kTop) - ks;
        mode(j, i) = (qEst < 0.0 && -qEst * dtg > avail(j, i)) ? 2 : 1;
      });
}

void RichardsSolver::computeFluxes(real_t dtg) {
  const bool full3d = useFull3d_;
  const bool coupled = cpl_.active;
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t az = mesh_.areaZ();

  Field3<real_t> h = h_, wc = wc_, qx = qx_, qy = qy_, qzF = qzF_;
  Field3<real_t> ksx = ksx_, ksy = ksy_, ksz = ksz_;
  Field3<real_t> vga = vga_, vgn = vgn_, wcs = wcs_, wcr = wcr_, aev = aev_;
  Field3<real_t> rXp = rRhoXp_, rYp = rRhoYp_, rZp = rRhoZp_;
  Field3<real_t> vXp = rViscXp_, vYp = rViscYp_, vZp = rViscZp_;
  Field3<real_t> dz3d = mesh_.dz3d(), bot3d = mesh_.bot3d();
  Field3<real_t> ax = mesh_.areaX(), ay = mesh_.areaY();
  Field3<real_t> sinx = mesh_.sinX(), cosx = mesh_.cosX();
  Field3<real_t> siny = mesh_.sinY(), cosy = mesh_.cosY();
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> topCode = topCode_, botCode = botCode_;
  Field2<int> sideXm = sideCodeXm_, sideXp = sideCodeXp_;
  Field2<int> sideYm = sideCodeYm_, sideYp = sideCodeYp_;
  Field2<real_t> topValue = topValue_, botValue = botValue_;
  Field2<real_t> sideVXm = sideValueXm_, sideVXp = sideValueXp_;
  Field2<real_t> sideVYm = sideValueYm_, sideVYp = sideValueYp_;
  Field2<real_t> cplDept = coupled ? cpl_.dept : Field2<real_t>();
  Field2<real_t> cplAvail = coupled ? cpl_.avail : Field2<real_t>();
  Field2<int> cplMode = coupled ? cplMode_ : Field2<int>();

  // x faces (darcy_flux "x" branch, subroutines.c:32-46). Slot (j, i, k) is
  // the face between cells i and i+1; slot (j, 0, k) is the west
  // boundary/interface face.
  Kokkos::parallel_for(
      "gw_flux_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 0, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const int giMinus = i0 + i - 1;
        const bool westEdgeFace = (giMinus < 0);
        const bool eastEdgeFace = (giMinus == nxGlobal - 1);
        const real_t faceA = ax(j, i, k);
        real_t q = 0.0;
        if (westEdgeFace) {
          // Legacy back-x: half spacing, no terrain term, one-sided K from
          // the interior cell (subroutines.c:33-36 with defined ghosts).
          const real_t ks = ksx(j, 1, k);
          const real_t kface = conductivityFromHead(
              soilAt(vga, vgn, wcs, wcr, aev, j, 1, k), h(j, 1, k), ks);
          q = faceA * kface * vXp(j, 0, k) * cosx(j, 0, k) * (h(j, 1, k) - h(j, 0, k)) /
              (0.5 * dx);
        } else {
          // Normal face from the minus cell (subroutines.c:37-45): the
          // minus cell's Ksx enters both flanks; the plus flank collapses
          // one-sided at the east domain edge (half spacing there).
          const real_t ks = ksx(j, i, k);
          const real_t kc = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                                 h(j, i, k), ks);
          real_t kface;
          real_t delta = dx;
          int sign = (bot3d(j, i + 1, k) > bot3d(j, i, k)) ? 1 : -1;
          if (eastEdgeFace) {
            kface = kc;
            delta = 0.5 * dx;
          } else {
            const real_t kn = conductivityFromHead(
                soilAt(vga, vgn, wcs, wcr, aev, j, i + 1, k), h(j, i + 1, k), ks);
            kface = 0.5 * (kc + kn);
            const bool minusInDomain = (giMinus >= 0);
            const bool plusInDomain = (giMinus + 1 <= nxGlobal - 1);
            if (minusInDomain && plusInDomain &&
                (gid(j, i, k) < 0 || gid(j, i + 1, k) < 0)) {
              kface = 0.0;  // subroutines.c:99
            }
          }
          q = faceA * (kface * vXp(j, i, k) * cosx(j, i, k) * (h(j, i + 1, k) - h(j, i, k)) /
                           delta +
                       static_cast<real_t>(sign) * kface * sinx(j, i, k));
          // use_full3d = false: column mode also suppresses the corrector's
          // lateral flux at faces the predictor sealed. Legacy only zeroed
          // the matrix conductivities (no benchmark can observe the
          // difference: use_full3d = 0 appears only on b2's single column);
          // letting the corrector leak laterally would contradict the
          // documented 1D-column semantics (plan §2.1).
          if (!full3d && ((gid(j, i, k) >= 0 && wc(j, i, k) < wcs(j, i, k)) ||
                          (gid(j, i + 1, k) >= 0 && wc(j, i + 1, k) < wcs(j, i + 1, k)))) {
            q = 0.0;
          }
        }
        // Domain-edge overrides: prescribed flux (positive into the domain)
        // or a sealed edge. Prescribed-head edges keep the Darcy value.
        if (westEdgeFace) {
          const int code = sideXm(j, 1);
          if (code == static_cast<int>(GwBcCode::NoFlux)) {
            q = 0.0;
          } else if (code == static_cast<int>(GwBcCode::Flux)) {
            q = -sideVXm(j, 1) * faceA;
          }
        } else if (eastEdgeFace) {
          const int code = sideXp(j, i);
          if (code == static_cast<int>(GwBcCode::NoFlux)) {
            q = 0.0;
          } else if (code == static_cast<int>(GwBcCode::Flux)) {
            q = sideVXp(j, i) * faceA;
          }
        }
        qx(j, i, k) = q;
      });

  // y faces (darcy_flux "y" branch, subroutines.c:47-69).
  Kokkos::parallel_for(
      "gw_flux_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {0, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const int gjMinus = j0 + j - 1;
        const bool southEdgeFace = (gjMinus < 0);
        const bool northEdgeFace = (gjMinus == nyGlobal - 1);
        const real_t faceA = ay(j, i, k);
        real_t q = 0.0;
        if (southEdgeFace) {
          const real_t ks = ksy(1, i, k);
          const real_t kface = conductivityFromHead(
              soilAt(vga, vgn, wcs, wcr, aev, 1, i, k), h(1, i, k), ks);
          q = faceA * kface * vYp(0, i, k) * cosy(0, i, k) * (h(1, i, k) - h(0, i, k)) /
              (0.5 * dy);
        } else {
          const real_t ks = ksy(j, i, k);
          const real_t kc = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                                 h(j, i, k), ks);
          real_t kface;
          real_t delta = dy;
          int sign = (bot3d(j + 1, i, k) > bot3d(j, i, k)) ? 1 : -1;
          if (northEdgeFace) {
            kface = kc;
            delta = 0.5 * dy;
          } else {
            const real_t kn = conductivityFromHead(
                soilAt(vga, vgn, wcs, wcr, aev, j + 1, i, k), h(j + 1, i, k), ks);
            kface = 0.5 * (kc + kn);
            const bool minusInDomain = (gjMinus >= 0);
            const bool plusInDomain = (gjMinus + 1 <= nyGlobal - 1);
            if (minusInDomain && plusInDomain &&
                (gid(j, i, k) < 0 || gid(j + 1, i, k) < 0)) {
              kface = 0.0;
            }
          }
          q = faceA * (kface * vYp(j, i, k) * cosy(j, i, k) * (h(j + 1, i, k) - h(j, i, k)) /
                           delta +
                       static_cast<real_t>(sign) * kface * siny(j, i, k));
          if (!full3d && ((gid(j, i, k) >= 0 && wc(j, i, k) < wcs(j, i, k)) ||
                          (gid(j + 1, i, k) >= 0 && wc(j + 1, i, k) < wcs(j + 1, i, k)))) {
            q = 0.0;
          }
        }
        if (southEdgeFace) {
          const int code = sideYm(1, i);
          if (code == static_cast<int>(GwBcCode::NoFlux)) {
            q = 0.0;
          } else if (code == static_cast<int>(GwBcCode::Flux)) {
            q = -sideVYm(1, i) * faceA;
          }
        } else if (northEdgeFace) {
          const int code = sideYp(j, i);
          if (code == static_cast<int>(GwBcCode::NoFlux)) {
            q = 0.0;
          } else if (code == static_cast<int>(GwBcCode::Flux)) {
            q = sideVYp(j, i) * faceA;
          }
        }
        qy(j, i, k) = q;
      });

  // z faces (darcy_flux "z" branch, subroutines.c:70-105 and the boundary
  // overrides :160-191). Plane k is the face above cell k; plane nz the
  // bottom boundary face. Positive flux points upward.
  Kokkos::parallel_for(
      "gw_flux_z",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz + 1}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const bool belowActive = (k < nz) && (gid(j, i, k) >= 0);
        const bool aboveActive = (k >= 1) && (gid(j, i, k - 1) >= 0);
        const bool isTopFace = belowActive && !aboveActive;
        real_t q = 0.0;
        if (isTopFace) {
          // The face above a column's top cell (darcy_flux back-z,
          // subroutines.c:71-75, with the groundwater-only overrides
          // :182-185). Legacy computed the wrong face for masked columns
          // (ktop > 0) — fixed on port; ktop = 0 behavior is identical.
          if (coupled) {
            // Coupled exchange face (darcy_flux :87-105 with the coupled
            // overrides :160-180 and the amendment-A13 supply mode): the
            // ghost head is the surface depth; the face conductivity is the
            // mean of the two flanks' closures with the top cell's
            // saturated Ksz — replaced by the pure saturated value under a
            // wet surface (:102-103).
            const real_t ks = ksz(j, i, k);
            const VgSoil soil = soilAt(vga, vgn, wcs, wcr, aev, j, i, k);
            const real_t depthHere = cplDept(j, i);
            const int modeHere = cplMode(j, i);
            if (modeHere == 2) {
              // Supply-limited: the exchange takes the remaining surface
              // water (SERGHEI GwBC.h:408, q = -h_ghost / dt), bounded by
              // the top cell's pore room up to the legacy 0.9999 theta_s
              // mark (see the matching predictor branch); the rest of the
              // pond waits for the next substep, by which the filled top
              // cell reclassifies to mode 1.
              const real_t roomTop = Kokkos::fmax(
                  0.0, (0.9999 * wcs(j, i, k) - wc(j, i, k)) * dz3d(j, i, k));
              const real_t take = Kokkos::fmin(cplAvail(j, i), roomTop);
              q = -take * az / dtg;
              cplAvail(j, i) -= take;
            } else if (modeHere == 1) {
              const real_t kface = ks;  // saturated top, :102-103
              const real_t delta = 0.5 * dz3d(j, i, k);
              q = az * kface * vZp(j, i, k) *
                  ((h(j, i, k) - depthHere) / delta - rZp(j, i, k));
              // Infiltration cannot exceed the surface water available
              // (:164-169). The legacy bound carried a spurious porosity
              // factor (vseep = |q| dtg wcs) — dropped per the plan §5.7
              // unit resolution (amendment A9); the budget is the window's
              // remaining depth, debited across subcycles.
              if (q < 0.0) {
                const real_t need = -q * dtg / az;
                if (need > cplAvail(j, i)) {
                  q = -cplAvail(j, i) * az / dtg;
                  cplAvail(j, i) = 0.0;
                } else {
                  cplAvail(j, i) -= need;
                }
              }
            } else {
              // Dry surface: a seepage face — infiltration from dry land is
              // prohibited (:171-173); a configured top flux adds the
              // legacy qtop source with its moisture guards (:174-178).
              const real_t kGhost = conductivityFromHead(soil, depthHere, ks);
              const real_t kCell = conductivityFromHead(soil, h(j, i, k), ks);
              const real_t kface = 0.5 * (kGhost + kCell);
              const real_t delta = 0.5 * dz3d(j, i, k);
              q = az * kface * vZp(j, i, k) *
                  ((h(j, i, k) - depthHere) / delta - rZp(j, i, k));
              if (q <= 0.0) {
                q = 0.0;
              }
              if (topCode(j, i) == static_cast<int>(GwBcCode::Flux)) {
                const real_t v = topValue(j, i);
                if ((v > 0.0 && wc(j, i, k) > wcr(j, i, k)) || v < 0.0) {
                  q += v * az;
                }
              }
            }
          } else {
            const int code = topCode(j, i);
            if (code == static_cast<int>(GwBcCode::Flux)) {
              q = topValue(j, i) * az;
            } else if (code == static_cast<int>(GwBcCode::Head)) {
              const real_t kface = ksz(j, i, k);  // saturated top, :102-105
              const real_t delta = 0.5 * dz3d(j, i, k);
              q = az * kface * vZp(j, i, k) *
                  ((h(j, i, k) - topValue(j, i)) / delta - rZp(j, i, k));
            } else {
              q = 0.0;
            }
          }
        } else if (k == nz && aboveActive) {
          // Bottom boundary face (subroutines.c:76-86, overrides :187-191
          // with the bctype index fix noted in Predictor.cpp).
          const int code = botCode(j, i);
          const real_t ks = ksz(j, i, k - 1);
          const real_t kface = conductivityFromHead(
              soilAt(vga, vgn, wcs, wcr, aev, j, i, k - 1), h(j, i, k - 1), ks);
          if (code == static_cast<int>(GwBcCode::Flux)) {
            q = botValue(j, i) * az;
          } else if (code == static_cast<int>(GwBcCode::Gravity)) {
            q = -kface * vZp(j, i, k) * rZp(j, i, k) * az;
          } else if (code == static_cast<int>(GwBcCode::Head)) {
            const real_t delta = 0.5 * dz3d(j, i, k - 1);
            q = az * kface * vZp(j, i, k) *
                ((botValue(j, i) - h(j, i, k - 1)) / delta - rZp(j, i, k));
          } else {
            q = 0.0;
          }
        } else if (aboveActive && belowActive) {
          // Interior face: the lower cell's saturated conductivity enters
          // both flanks (subroutines.c:77, :90-98).
          const real_t ks = ksz(j, i, k);
          const real_t ka = conductivityFromHead(
              soilAt(vga, vgn, wcs, wcr, aev, j, i, k - 1), h(j, i, k - 1), ks);
          const real_t kb = conductivityFromHead(soilAt(vga, vgn, wcs, wcr, aev, j, i, k),
                                                 h(j, i, k), ks);
          const real_t kface = 0.5 * (ka + kb);
          const real_t delta = 0.5 * (dz3d(j, i, k - 1) + dz3d(j, i, k));
          q = az * kface * vZp(j, i, k) *
              ((h(j, i, k) - h(j, i, k - 1)) / delta - rZp(j, i, k));
        }
        qzF(j, i, k) = q;
      });
}

void RichardsSolver::applyCoupledTopBookkeeping(real_t dtg) {
  // Coupled top-face bookkeeping (legacy groundwater_flux,
  // groundwater.c:835-868), run right after the fluxes so the corrector's
  // divergence still sees the limited flux:
  //  - the saturation bounce-back (:844-851): downward flux into a top cell
  //    without pore room under a dry surface returns to the surface at once;
  //  - the seepage accumulation (:852-857): the top-face exchange joins the
  //    per-column accumulator the coupler applies to eta at the end of the
  //    surface step. Legacy accumulated the *rate* and multiplied by the
  //    surface dt on application, which over-draws when dtg < dt; the
  //    accumulator carries the exchanged volume per unit area instead
  //    (amendment A9) — identical in sync mode where dtg = dt.
  //  - the evaporation correction (:858-863): a positive configured top flux
  //    (evaporation) leaving a moist, surface-dry column does not pond.
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t az = mesh_.areaZ();
  const real_t minDepth = cpl_.minDepth;

  Field3<real_t> wc = wc_, wcs = wcs_, wcr = wcr_, qzF = qzF_;
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> topCode = topCode_;
  Field2<real_t> topValue = topValue_;
  Field2<real_t> eta = cpl_.eta, dept = cpl_.dept, avail = cpl_.avail;
  Field2<real_t> seepAccum = cpl_.seepAccum, gain = cpl_.gain;

  real_t exchanged = 0.0;
  real_t bounced = 0.0;
  real_t evaporated = 0.0;
  Kokkos::parallel_reduce(
      "gw_coupled_top",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumE, real_t& sumB, real_t& sumV) {
        int kTop = -1;
        for (int k = 0; k < nz; ++k) {
          if (gid(j, i, k) >= 0) {
            kTop = k;
            break;
          }
        }
        if (kTop < 0) {
          return;
        }
        const real_t q = qzF(j, i, kTop);
        sumE += q * dtg;
        if (q < 0.0 && dept(j, i) <= minDepth) {
          const real_t volume = az * dz3d(j, i, kTop);
          if ((wcs(j, i, kTop) - wc(j, i, kTop)) < -q * dtg / volume) {
            const real_t d = -q * dtg / az;
            eta(j, i) += d;
            dept(j, i) += d;
            avail(j, i) += d;
            gain(j, i) += d;
            sumB += d * az;
          }
        }
        seepAccum(j, i) += q * dtg / az;
        if (topCode(j, i) == static_cast<int>(GwBcCode::Flux) && topValue(j, i) > 0.0 &&
            wc(j, i, kTop) > wcr(j, i, kTop) && dept(j, i) <= minDepth) {
          seepAccum(j, i) -= topValue(j, i) * dtg;
          sumV += topValue(j, i) * dtg * az;
        }
      },
      exchanged, bounced, evaporated);
  audit_.cplExchanged = exchanged;
  audit_.cplBounce = bounced;
  audit_.cplEvap = evaporated;
}

void RichardsSolver::updateWaterContent(real_t dtg) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t az = mesh_.areaZ();
  const real_t ss = ss_;

  Field3<real_t> h = h_, hn = hn_, wc = wc_, wcn = wcn_;
  Field3<real_t> qx = qx_, qy = qy_, qzF = qzF_;
  Field3<real_t> wcs = wcs_;
  Field3<real_t> rRho = rRho_, rXp = rRhoXp_, rYp = rRhoYp_, rZp = rRhoZp_;
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<PetscInt> gid = grid_.gid3();

  // θ from the flux divergence with the Ss compressibility factor
  // (update_water_content:908-925). The reduction tracks the volume moved
  // into compressible storage for the mass audit.
  real_t ssVolume = 0.0;
  Kokkos::parallel_reduce(
      "gw_update_wc",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sum) {
        if (gid(j, i, k) < 0) {
          return;
        }
        const real_t volume = az * dz3d(j, i, k);
        const real_t coeff =
            volume * (rRho(j, i, k) +
                      rRho(j, i, k) * ss * (h(j, i, k) - hn(j, i, k)) / wcs(j, i, k));
        const real_t dqx = dtg * (qx(j, i, k) * rXp(j, i, k) -
                                  qx(j, i - 1, k) * rXp(j, i - 1, k));
        const real_t dqy = dtg * (qy(j, i, k) * rYp(j, i, k) -
                                  qy(j - 1, i, k) * rYp(j - 1, i, k));
        const real_t dqz = dtg * (qzF(j, i, k + 1) * rZp(j, i, k + 1) -
                                  qzF(j, i, k) * rZp(j, i, k));
        const real_t updated =
            (wcn(j, i, k) * rRho(j, i, k) * volume + dqx + dqy + dqz) / coeff;
        wc(j, i, k) = updated;
        sum += updated * volume * rRho(j, i, k) * ss * (h(j, i, k) - hn(j, i, k)) /
               wcs(j, i, k);
      },
      ssVolume);
  audit_.ssStorage = ssVolume;
}

void RichardsSolver::enforceMoistureBc() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const bool west = (grid_.rankWest() == MPI_PROC_NULL);
  const bool east = (grid_.rankEast() == MPI_PROC_NULL);
  const bool south = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool north = (grid_.rankNorth() == MPI_PROC_NULL);

  Field3<real_t> wc = wc_, wcs = wcs_;
  Field2<int> sideXm = sideCodeXm_, sideXp = sideCodeXp_;
  Field2<int> sideYm = sideCodeYm_, sideYp = sideCodeYp_;

  // Domain-edge moisture ghosts (enforce_moisture_bc:933-956): zero-gradient
  // copies; prescribed-head edges are saturated ghosts (:947-948).
  Kokkos::parallel_for(
      "gw_wc_ghost_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nyl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int k) {
        const auto isHead = [](int code) {
          return code == static_cast<int>(GwBcCode::Head) ||
                 code == static_cast<int>(GwBcCode::HeadHydrostatic);
        };
        if (west) {
          wc(j, 0, k) = isHead(sideXm(j, 1)) ? wcs(j, 1, k) : wc(j, 1, k);
        }
        if (east) {
          wc(j, nxl + 1, k) = isHead(sideXp(j, nxl)) ? wcs(j, nxl, k) : wc(j, nxl, k);
        }
      });
  Kokkos::parallel_for(
      "gw_wc_ghost_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nxl + 1, nz}),
      KOKKOS_LAMBDA(const int i, const int k) {
        const auto isHead = [](int code) {
          return code == static_cast<int>(GwBcCode::Head) ||
                 code == static_cast<int>(GwBcCode::HeadHydrostatic);
        };
        if (south) {
          wc(0, i, k) = isHead(sideYm(1, i)) ? wcs(1, i, k) : wc(1, i, k);
        }
        if (north) {
          wc(nyl + 1, i, k) = isHead(sideYp(nyl, i)) ? wcs(nyl, i, k) : wc(nyl, i, k);
        }
      });
}

void RichardsSolver::finalizeWaterContent() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t az = mesh_.areaZ();

  Field3<real_t> wc = wc_, wcs = wcs_, wcr = wcr_, vloss = vloss_, h = h_;
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<PetscInt> gid = grid_.gid3();

  // The final θ clamp with volume-loss accounting (solve_groundwater
  // :165-186); the reduction feeds the mass audit.
  real_t lost = 0.0;
  Kokkos::parallel_reduce(
      "gw_finalize_wc",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sum) {
        if (gid(j, i, k) < 0) {
          return;
        }
        const real_t volume = az * dz3d(j, i, k);
        if (wc(j, i, k) > wcs(j, i, k)) {
          const real_t excess = (wc(j, i, k) - wcs(j, i, k)) * volume;
          vloss(j, i, k) += excess;
          sum += excess;
          wc(j, i, k) = wcs(j, i, k);
        } else if (wc(j, i, k) < wcr(j, i, k)) {
          const real_t deficit = (wcr(j, i, k) - wc(j, i, k)) * volume;
          vloss(j, i, k) -= deficit;
          sum -= deficit;
          wc(j, i, k) = wcr(j, i, k);
        }
      },
      lost);
  audit_.vloss = lost;

#ifndef NDEBUG
  // Debug-build invariant (plan §10 P2 exit): θ within [θr, θs] and finite
  // head on every active cell at the end of the step.
  long violations = 0;
  Kokkos::parallel_reduce(
      "gw_debug_bounds",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, long& count) {
        if (gid(j, i, k) < 0) {
          return;
        }
        const bool ok = wc(j, i, k) >= wcr(j, i, k) && wc(j, i, k) <= wcs(j, i, k) &&
                        Kokkos::isfinite(h(j, i, k));
        if (!ok) {
          ++count;
        }
      },
      violations);
  if (violations > 0) {
    log::fatal(log::msg() << "groundwater state invariant violated on " << violations
                          << " cell(s): theta outside [theta_r, theta_s] or non-finite head");
  }
#endif
}

}  // namespace frehg::gw
