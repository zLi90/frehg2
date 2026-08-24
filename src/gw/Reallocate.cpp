/// \file Reallocate.cpp
/// \brief The PCA post-allocation step (legacy reallocate_water_content,
///        check_adj_sat, check_head_gradient, allocate_send;
///        groundwater.c:959-1424) — always on, per plan Appendix A
///        ("post-allocation always on"; the post_allocate knob is dropped).
///
/// The step restores θ/h consistency after the corrector and redistributes
/// over-saturation excess:
///  - **over-saturated cells** (θ ≥ θs) send their excess along the head
///    gradient (check_head_gradient's six-direction split, :1074-1150):
///    the vertical fractions walk the column up/down into available pore
///    room exactly as allocate_send does (:1153-1291), venting through a
///    prescribed-head top when the column fills; the cell clamps to θs;
///  - **unsaturated cells adjacent to a saturated cell** (or to a
///    prescribed-head boundary) take θ from the retention curve at the
///    predicted head, θ := θ(h) (:1001-1036);
///  - **isolated unsaturated cells** re-derive the head from the conserved
///    moisture, h := h(θ), unless effectively saturated (0.9999 θs,
///    :995-999).
///
/// Deviations from the legacy sweep, all recorded in the plan amendment log
/// and docs/theory/groundwater.md:
///  - cells are classified against the pre-reallocation state and the
///    vertical walks run as one deterministic sequential sweep per column
///    (the plan §5.2 design statement: "the Reallocate sweep is sequential
///    per column"). The legacy flat-order sweep interleaved sends with
///    classification, which made results depend on the rank decomposition.
///  - pore room is evaluated fresh (θs - θ) when the sweep runs; legacy
///    check_room measured it before the corrector (groundwater.c:139), and
///    the stale value could overfill receivers past θs — violating the θ
///    invariant the plan's P2 exit criterion asserts.
///  - the lateral (x/y) fractions of the split are not delivered to
///    neighbor columns; they are dropped and recorded in the audit. Legacy
///    delivered them order-dependently and dropped any remainder its
///    single-neighbor probe could not place (:1293-1400, retry disabled at
///    :983).
///  - check_head_gradient's up-face spacing added a cell index to a
///    thickness (:1110, the §2.1-flagged index bug); the intended
///    0.5 (dz(k) + dz(k-1)) is used.
/// Everything each cell publishes in the parallel phase touches only that
/// cell; the column sweep touches only its own column — race-free and
/// decomposition-invariant by construction.

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

void RichardsSolver::reallocateWaterContent() {
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
  const bool full3d = useFull3d_;
  const bool redistribute = surplusRedistribute_;
  const real_t dtg = dtgCurrent_;

  Field3<real_t> h = h_, wc = wc_;
  Field3<real_t> kx = kx_, ky = ky_, kzF = kzF_, qzF = qzF_;
  Field3<real_t> vga = vga_, vgn = vgn_, wcs = wcs_, wcr = wcr_, aev = aev_;
  Field3<real_t> dz3d = mesh_.dz3d();
  Field3<real_t> room = room_, sendUp = sendUp_, sendDown = sendDown_;
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> topCode = topCode_;
  Field2<real_t> topValue = topValue_;
  const bool coupled = cpl_.active;
  const real_t minDepth = cpl_.minDepth;
  Field2<real_t> cplEta = coupled ? cpl_.eta : Field2<real_t>();
  Field2<real_t> cplDept = coupled ? cpl_.dept : Field2<real_t>();
  Field2<real_t> cplAvail = coupled ? cpl_.avail : Field2<real_t>();
  Field2<real_t> cplGain = coupled ? cpl_.gain : Field2<real_t>();

  // Phase 1 (cell-parallel): classification, the θ/h consistency restore,
  // fresh pore room, and each over-saturated cell's gradient split. The
  // reductions track every θ change and, separately, the volume discarded
  // outright (lateral fractions; the adjacent-cell surplus in drop mode).
  real_t adjusted = 0.0;
  real_t discarded = 0.0;
  Kokkos::parallel_reduce(
      "gw_reallocate_classify",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sum, real_t& drop) {
        if (gid(j, i, k) < 0) {
          room(j, i, k) = 0.0;
          sendUp(j, i, k) = 0.0;
          sendDown(j, i, k) = 0.0;
          return;
        }
        const VgSoil soil = soilAt(vga, vgn, wcs, wcr, aev, j, i, k);
        const real_t volume = az * dz3d(j, i, k);
        const real_t before = wc(j, i, k);
        sendUp(j, i, k) = 0.0;
        sendDown(j, i, k) = 0.0;

        // The six-direction gradient split of check_head_gradient
        // (:1074-1150), shared by both send paths. Returns the up/down
        // fractions; lateral fractions are the remainder.
        struct SplitResult {
          real_t up;
          real_t down;
        };
        const auto gradientSplit = [&](int jc, int ic, int kc) -> SplitResult {
          real_t dh6[6];
          const auto active = [&](int jj, int ii, int kk) {
            return kk >= 0 && kk < nz && gid(jj, ii, kk) >= 0;
          };
          dh6[0] = (active(jc, ic + 1, kc) && kx(jc, ic, kc) > 0.0)
                       ? (h(jc, ic + 1, kc) - h(jc, ic, kc)) / dx
                       : 0.0;
          dh6[1] = (active(jc, ic - 1, kc) && kx(jc, ic - 1, kc) > 0.0)
                       ? (h(jc, ic, kc) - h(jc, ic - 1, kc)) / dx
                       : 0.0;
          dh6[2] = (active(jc + 1, ic, kc) && ky(jc, ic, kc) > 0.0)
                       ? (h(jc + 1, ic, kc) - h(jc, ic, kc)) / dy
                       : 0.0;
          dh6[3] = (active(jc - 1, ic, kc) && ky(jc - 1, ic, kc) > 0.0)
                       ? (h(jc, ic, kc) - h(jc - 1, ic, kc)) / dy
                       : 0.0;
          if (active(jc, ic, kc + 1) && kzF(jc, ic, kc + 1) > 0.0) {
            const real_t dzp = 0.5 * (dz3d(jc, ic, kc) + dz3d(jc, ic, kc + 1));
            dh6[4] = (h(jc, ic, kc + 1) - h(jc, ic, kc)) / dzp - 1.0;
          } else {
            dh6[4] = 0.0;
          }
          if (active(jc, ic, kc - 1) && kzF(jc, ic, kc) > 0.0) {
            const real_t dzm = 0.5 * (dz3d(jc, ic, kc) + dz3d(jc, ic, kc - 1));
            dh6[5] = (h(jc, ic, kc) - h(jc, ic, kc - 1)) / dzm - 1.0;
          } else {
            dh6[5] = 0.0;
          }
          if (!full3d) {
            dh6[0] = dh6[1] = dh6[2] = dh6[3] = 0.0;
          }
          real_t gradTot = 0.0;
          real_t gradx, grady, gradz;
          if (dh6[0] * dh6[1] >= 0.0) {
            gradx = 0.5 * (dh6[0] + dh6[1]);
            gradTot += Kokkos::fabs(gradx);
          } else {
            gradx = 0.0;
            gradTot += Kokkos::fabs(dh6[0]) + Kokkos::fabs(dh6[1]);
          }
          if (dh6[2] * dh6[3] >= 0.0) {
            grady = 0.5 * (dh6[2] + dh6[3]);
            gradTot += Kokkos::fabs(grady);
          } else {
            grady = 0.0;
            gradTot += Kokkos::fabs(dh6[2]) + Kokkos::fabs(dh6[3]);
          }
          if (dh6[4] * dh6[5] >= 0.0) {
            gradz = 0.5 * (dh6[4] + dh6[5]);
            gradTot += Kokkos::fabs(gradz);
          } else {
            gradz = 0.0;
            gradTot += Kokkos::fabs(dh6[4]) + Kokkos::fabs(dh6[5]);
          }
          SplitResult split{0.0, 0.0};
          if (gradTot > 0.0) {
            if (gradz > 0.0) {
              split.up = gradz / gradTot;
            } else if (gradz < 0.0) {
              split.down = -gradz / gradTot;
            } else {
              split.down = Kokkos::fabs(dh6[4]) / gradTot;
              split.up = Kokkos::fabs(dh6[5]) / gradTot;
            }
          }
          return split;
        };

        if (before >= soil.thetaS) {
          // Over-saturated: split the excess along the head gradient
          // (check_head_gradient, :1074-1150) and clamp to θs.
          const real_t excess = (before - soil.thetaS) * volume;
          const SplitResult split = gradientSplit(j, i, k);
          sendUp(j, i, k) = excess * split.up;
          sendDown(j, i, k) = excess * split.down;
          // Lateral fractions (and a zero-gradient split's whole excess)
          // are dropped and audited; see the file comment.
          sum -= excess - sendUp(j, i, k) - sendDown(j, i, k);
          drop += excess - sendUp(j, i, k) - sendDown(j, i, k);
          wc(j, i, k) = soil.thetaS;
          room(j, i, k) = 0.0;
          return;
        }

        room(j, i, k) = (soil.thetaS - before) * volume;

        // check_adj_sat (groundwater.c:1044-1071): saturated neighbor
        // through a conducting face, or a wet/prescribed-head top.
        const auto sat = [&](int jj, int ii, int kk) {
          return wc(jj, ii, kk) >= wcs(jj, ii, kk);
        };
        bool adjSat = false;
        if (k + 1 < nz && gid(j, i, k + 1) >= 0 && sat(j, i, k + 1)) {
          adjSat = true;
        } else if (kx(j, i, k) != 0.0 && sat(j, i + 1, k)) {
          adjSat = true;
        } else if (kx(j, i - 1, k) != 0.0 && sat(j, i - 1, k)) {
          adjSat = true;
        } else if (ky(j, i, k) != 0.0 && sat(j + 1, i, k)) {
          adjSat = true;
        } else if (ky(j - 1, i, k) != 0.0 && sat(j - 1, i, k)) {
          adjSat = true;
        } else {
          const bool isTop = (k == 0) || (gid(j, i, k - 1) < 0);
          if (isTop) {
            // A wet surface counts as saturated contact in coupled runs
            // (check_adj_sat:1061); groundwater-only runs test a
            // non-negative prescribed top head (:1062).
            if (coupled) {
              if (cplDept(j, i) > 0.0) {
                adjSat = true;
              }
            } else if (topCode(j, i) == static_cast<int>(GwBcCode::Head) &&
                       topValue(j, i) >= 0.0) {
              adjSat = true;
            }
          } else if (gid(j, i, k - 1) >= 0 && sat(j, i, k - 1)) {
            adjSat = true;
          }
        }

        if (!adjSat) {
          // Isolated unsaturated cell: restore h from the conserved θ
          // (:995-999).
          if (before < 0.9999 * soil.thetaS) {
            h(j, i, k) = headFromWaterContent(soil, before);
          }
        } else {
          // Saturated contact: θ follows the predicted head (:1034). When
          // the conserved θ exceeds θ(h), the surplus is sent down the head
          // gradient — the legacy alloc_type4 path computes exactly this dV
          // and gradient split (:1020-1025); its allocate_send call is
          // re-enabled here because discarding the surplus starves the b3
          // saturation bulb (amendment A7). The deficit direction (θ(h)
          // above the conserved θ) stays exactly as legacy left it: no
          // neighbor withdrawal (:1013-1014 disabled it for stability), the
          // created volume is audited.
          const real_t after = waterContentFromHead(soil, h(j, i, k));
          wc(j, i, k) = after;
          room(j, i, k) = (soil.thetaS - after) * volume;
          if (redistribute && after < before) {
            const real_t surplus = (before - after) * volume;
            const SplitResult split = gradientSplit(j, i, k);
            sendUp(j, i, k) = surplus * split.up;
            sendDown(j, i, k) = surplus * split.down;
            // Lateral fractions dropped and audited, as in the oversat
            // branch.
            sum -= surplus - sendUp(j, i, k) - sendDown(j, i, k);
            drop += surplus - sendUp(j, i, k) - sendDown(j, i, k);
          } else {
            sum += (after - before) * volume;
            if (after < before) {
              drop += (before - after) * volume;  // drop-mode surplus
            }
          }
        }
      },
      adjusted, discarded);

  // Phase 2 (column-parallel): the vertical send walks (allocate_send,
  // :1158-1291), k ascending exactly as the legacy flat sweep visited
  // cells. Deposits and vents are mirrored into the z-flux field the same
  // way legacy maintains qz, so the audit and adaptive controller see the
  // moved volume.
  real_t dropped = 0.0;
  real_t vented = 0.0;
  Kokkos::parallel_reduce(
      "gw_reallocate_send",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& lost, real_t& vent) {
        (void)i0;
        (void)j0;
        (void)nxGlobal;
        (void)nyGlobal;
        for (int k = 0; k < nz; ++k) {
          if (gid(j, i, k) < 0) {
            continue;
          }
          const bool topIsHead = (topCode(j, i) == static_cast<int>(GwBcCode::Head));
          // Send up (:1158-1249).
          real_t dV = sendUp(j, i, k);
          if (dV > 0.0) {
            int ll = k;
            while (ll > 0 && gid(j, i, ll - 1) >= 0 && dV > 0.0) {
              --ll;
              if (room(j, i, ll) > 0.0) {
                const real_t take = (room(j, i, ll) > dV) ? dV : room(j, i, ll);
                wc(j, i, ll) += take / (az * dz3d(j, i, ll));
                room(j, i, ll) -= take;
                qzF(j, i, ll + 1) += take / dtg;
                dV -= take;
              }
            }
            if (dV > 0.0) {
              if (coupled) {
                // Coupled columns vent onto the surface (allocate_send
                // :1185-1207): any amount when the surface is wet, only
                // above min_depth when dry — smaller dry-cell vents are
                // discarded outright (:1199-1206), audited here.
                const real_t d = dV / az;
                if (cplDept(j, i) > 0.0 || d > minDepth) {
                  cplEta(j, i) += d;
                  cplDept(j, i) += d;
                  cplAvail(j, i) += d;
                  cplGain(j, i) += d;
                  qzF(j, i, ll) += dV / dtg;
                  vent += dV;
                } else {
                  lost -= dV;
                }
                dV = 0.0;
              } else if (topIsHead) {
                // Release through the prescribed-head surface (:1210-1214).
                qzF(j, i, ll) += dV / dtg;
                dV = 0.0;
              } else {
                // No-flux/flux top: send the remainder back down (:1216-1240).
                while (ll + 1 < nz && dV > 0.0) {
                  ++ll;
                  if (room(j, i, ll) > 0.0) {
                    const real_t take = (room(j, i, ll) > dV) ? dV : room(j, i, ll);
                    wc(j, i, ll) += take / (az * dz3d(j, i, ll));
                    room(j, i, ll) -= take;
                    qzF(j, i, ll) -= take / dtg;
                    dV -= take;
                  }
                }
                lost -= dV;  // remainder dropped (legacy :1243-1248)
                dV = 0.0;
              }
            }
          }
          // Send down (:1252-1291).
          dV = sendDown(j, i, k);
          if (dV > 0.0) {
            int ll = k;
            while (ll + 1 < nz && dV > 0.0) {
              ++ll;
              if (room(j, i, ll) > 0.0) {
                const real_t take = (room(j, i, ll) > dV) ? dV : room(j, i, ll);
                wc(j, i, ll) += take / (az * dz3d(j, i, ll));
                room(j, i, ll) -= take;
                qzF(j, i, ll) -= take / dtg;
                dV -= take;
              }
            }
            lost -= dV;  // remainder at the bottom dropped (:1279-1291)
          }
        }
      },
      dropped, vented);
  audit_.reallocAdjust = adjusted + dropped;
  audit_.reallocDropped = discarded - dropped;  // dropped accumulates negative
  audit_.cplVent = vented;
}

}  // namespace frehg::gw
