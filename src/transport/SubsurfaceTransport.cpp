/// \file SubsurfaceTransport.cpp
/// \brief Subsurface scalar transport (scalar_groundwater, scalar.c:303-498,
///        with advective_flux :501-716 and dispersive_flux :719-847).
///
/// The step rebuilds the dispersion tensor from the window's last-substep
/// fluxes, advances the scalar mass s Vgn by the advective and dispersive
/// face fluxes over the surface step dt, removes the surface-exchange mass
/// (the sseepage the surface step computed), divides by the flux volume
/// Vgflux (volume_by_flux_subs, groundwater.c:1634-1644 — the pre-step pore
/// volume plus the last substep's flux divergence), and applies the local
/// min/max limiter, the configured bounds, and the ghost rules.
///
/// Legacy structure preserved exactly where defined; the deviations are all
/// of the determinism/decomposition class (docs/theory/transport.md):
///  - the superbee far-neighbor and cross-term guards are evaluated against
///    global domain edges with staged far values and corner-exchanged
///    scalars at rank interfaces (legacy degraded interface stencils to
///    upwind and read unexchanged corner ghosts — rank-count-dependent;
///    the serial b6 goldens never see either);
///  - the dead in-branch boundary doublings of dispersive_flux (:786, :831
///    — unreachable behind the ghost actv == 0 guards) are dropped; the
///    live jjp/jjm edge doublings (:348-350) are kept;
///  - the legacy top-face advective leak of uncoupled superbee runs with a
///    head-condition top reduces to the donor value (tvd of an equal pair),
///    which is what this port evaluates directly.

#include "transport/Limiters.hpp"
#include "transport/ScalarSolver.hpp"

#include "core/Logger.hpp"

namespace frehg::transport {

namespace {

/// Boundary codes as staged by the groundwater module's marker fields (the
/// ints decode the legacy bctype_GW values — gw::GwBcCode; transport reads
/// them through the wiring without including the module's headers, plan §4).
constexpr int kBcNoFlux = 0;
constexpr int kBcFlux = 2;

/// View bundle for the face-flux device helpers.
struct SubViews {
  Field3<real_t> s, dxx, dyy, dzz, dxy, dxz, dyz;
  Field3<real_t> ax, ay, cosx, cosy, dz3d;
  Field3<real_t> kx, ky, kzLower;
  Field3<PetscInt> gid;
  Field2<int> kTop;
  Field2<real_t> sSurf;
  real_t az = 0.0, dx = 0.0, dy = 0.0;
  int nz = 0;
  bool coupled = false;
};

/// Scalar "above" cell (j, i, k): the in-column neighbor, or the kM ghost
/// at a column top — the surface scalar in coupled runs, the top cell
/// itself otherwise (enforce_scalar_bc, scalar.c:897-909).
KOKKOS_INLINE_FUNCTION real_t scalarAbove(const SubViews& v, int j, int i, int k) {
  const int kTop = v.kTop(j, i);
  if (kTop >= 0 && k == kTop) {
    return v.coupled ? v.sSurf(j, i) : v.s(j, i, k);
  }
  return (k > 0) ? v.s(j, i, k - 1) : v.s(j, i, k);
}

/// Scalar "below" cell (j, i, k): zero-gradient past the box bottom.
KOKKOS_INLINE_FUNCTION real_t scalarBelow(const SubViews& v, int j, int i, int k) {
  return (k + 1 < v.nz) ? v.s(j, i, k + 1) : v.s(j, i, k);
}

/// Dispersive flux through the x+ face of cell (j, i, k) — the legacy
/// dispersive_flux(icell, "x") (scalar.c:727-763): the straight Dxx term
/// plus the Dxy/Dxz cross terms, every coefficient at the minus cell.
KOKKOS_INLINE_FUNCTION real_t dispFluxX(const SubViews& v, int j, int i, int k) {
  if (v.gid(j, i, k) < 0) {
    return 0.0;
  }
  real_t fx = 0.0;
  real_t fy = 0.0;
  real_t fz = 0.0;
  if (v.kx(j, i, k) > 0.0 && v.gid(j, i + 1, k) >= 0) {
    const real_t dist = v.dx / v.cosx(j, i, k);
    fx = v.dxx(j, i, k) * v.ax(j, i, k) * (v.s(j, i + 1, k) - v.s(j, i, k)) / dist;
  }
  if (v.ky(j, i, k) > 0.0 && v.gid(j + 1, i, k) >= 0) {
    const real_t syp = 0.25 * (v.s(j, i, k) + v.s(j, i + 1, k) + v.s(j + 1, i, k) +
                               v.s(j + 1, i + 1, k));
    const real_t sym = 0.25 * (v.s(j, i, k) + v.s(j, i + 1, k) + v.s(j - 1, i, k) +
                               v.s(j - 1, i + 1, k));
    const real_t dist = v.dy / v.cosy(j, i, k);
    fy = v.dxy(j, i, k) * v.ay(j, i, k) * (syp - sym) / dist;
  }
  if (v.kzLower(j, i, k) > 0.0 && k + 1 < v.nz && v.gid(j, i, k + 1) >= 0) {
    const real_t szp = 0.25 * (v.s(j, i, k) + v.s(j, i + 1, k) + v.s(j, i, k + 1) +
                               v.s(j, i + 1, k + 1));
    const real_t szm = 0.25 * (v.s(j, i, k) + v.s(j, i + 1, k) + scalarAbove(v, j, i, k) +
                               scalarAbove(v, j, i + 1, k));
    fz = v.dxz(j, i, k) * v.az * (szp - szm) / v.dz3d(j, i, k);
  }
  return fx + fy + fz;
}

/// Dispersive flux through the y+ face of cell (j, i, k) — the legacy
/// dispersive_flux(icell, "y") (scalar.c:764-801).
KOKKOS_INLINE_FUNCTION real_t dispFluxY(const SubViews& v, int j, int i, int k) {
  if (v.gid(j, i, k) < 0) {
    return 0.0;
  }
  real_t fx = 0.0;
  real_t fy = 0.0;
  real_t fz = 0.0;
  if (v.kx(j, i, k) > 0.0 && v.gid(j, i + 1, k) >= 0) {
    const real_t sxp = 0.25 * (v.s(j, i, k) + v.s(j, i + 1, k) + v.s(j + 1, i, k) +
                               v.s(j + 1, i + 1, k));
    const real_t sxm = 0.25 * (v.s(j, i, k) + v.s(j, i - 1, k) + v.s(j + 1, i, k) +
                               v.s(j + 1, i - 1, k));
    const real_t dist = v.dx / v.cosx(j, i, k);
    fx = v.dxy(j, i, k) * v.ay(j, i, k) * (sxp - sxm) / dist;
  }
  if (v.ky(j, i, k) > 0.0 && v.gid(j + 1, i, k) >= 0) {
    const real_t dist = v.dy / v.cosy(j, i, k);
    fy = v.dyy(j, i, k) * v.ay(j, i, k) * (v.s(j + 1, i, k) - v.s(j, i, k)) / dist;
  }
  if (v.kzLower(j, i, k) > 0.0 && k + 1 < v.nz && v.gid(j, i, k + 1) >= 0) {
    const real_t szp = 0.25 * (v.s(j, i, k) + v.s(j + 1, i, k) + v.s(j, i, k + 1) +
                               v.s(j + 1, i, k + 1));
    const real_t szm = 0.25 * (v.s(j, i, k) + v.s(j + 1, i, k) + scalarAbove(v, j, i, k) +
                               scalarAbove(v, j + 1, i, k));
    fz = v.dyz(j, i, k) * v.az * (szp - szm) / v.dz3d(j, i, k);
  }
  return fx + fy + fz;
}

/// Dispersive flux through the *lower* face of cell (j, i, k) — the legacy
/// dispersive_flux(icell, "z") for an active cell (scalar.c:802-833).
KOKKOS_INLINE_FUNCTION real_t dispFluxZLower(const SubViews& v, int j, int i, int k) {
  if (v.gid(j, i, k) < 0) {
    return 0.0;
  }
  real_t fx = 0.0;
  real_t fy = 0.0;
  real_t fz = 0.0;
  const bool belowActive = (k + 1 < v.nz) && (v.gid(j, i, k + 1) >= 0);
  if (v.kx(j, i, k) > 0.0 && v.gid(j, i + 1, k) >= 0 && belowActive) {
    const real_t sxp = 0.25 * (v.s(j, i, k) + v.s(j, i + 1, k) + v.s(j, i, k + 1) +
                               v.s(j, i + 1, k + 1));
    const real_t sxm = 0.25 * (v.s(j, i, k) + v.s(j, i - 1, k) + v.s(j, i, k + 1) +
                               v.s(j, i - 1, k + 1));
    const real_t dist = v.dx / v.cosx(j, i, k);
    fx = v.dxz(j, i, k) * v.ax(j, i, k) * (sxp - sxm) / dist;
  }
  if (v.ky(j, i, k) > 0.0 && v.gid(j + 1, i, k) >= 0 && belowActive) {
    const real_t syp = 0.25 * (v.s(j, i, k) + v.s(j, i, k + 1) + v.s(j + 1, i, k) +
                               v.s(j + 1, i, k + 1));
    const real_t sym = 0.25 * (v.s(j, i, k) + v.s(j, i, k + 1) + v.s(j - 1, i, k) +
                               v.s(j - 1, i, k + 1));
    const real_t dist = v.dy / v.cosy(j, i, k);
    fy = v.dyz(j, i, k) * v.ay(j, i, k) * (syp - sym) / dist;
  }
  if (v.kzLower(j, i, k) > 0.0 && belowActive) {
    fz = v.dzz(j, i, k) * v.az * (v.s(j, i, k + 1) - v.s(j, i, k)) / v.dz3d(j, i, k);
  }
  return fx + fy + fz;
}

}  // namespace

void ScalarSolver::stageSubsurfaceFarNeighbors() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  Field3<real_t> s = sSubs_;
  Field3<real_t> fxm = sFarXm3_, fxp = sFarXp3_, fym = sFarYm3_, fyp = sFarYp3_;
  Field3<real_t> kzLower = kzLower_, kzF = subs_.kzF;
  Kokkos::parallel_for(
      "transport_far_3d",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        fxm(j, i, k) = s(j, i - 1, k);
        fxp(j, i, k) = s(j, i + 1, k);
        fym(j, i, k) = s(j - 1, i, k);
        fyp(j, i, k) = s(j + 1, i, k);
        kzLower(j, i, k) = kzF(j, i, k + 1);
      });
  halo_.exchange({"s_far_xm_3d", "s_far_xp_3d", "s_far_ym_3d", "s_far_yp_3d", "s_kz_lower",
                  "s_gw_kx", "s_gw_ky"});
}

void ScalarSolver::stepSubsurface(real_t t, real_t dt, real_t dtgLast) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const bool superbee = superbee_;
  const bool coupled = surf_.active && cpl_.active;
  const real_t limHi = boundMax_;
  const real_t limLo = boundMin_;
  const bool hasMax = hasBoundMax_;

  stageSubsurfaceFarNeighbors();
  updateDispersionTensor();

  SubViews v;
  v.s = sSubs_;
  v.dxx = dxx_;
  v.dyy = dyy_;
  v.dzz = dzz_;
  v.dxy = dxy_;
  v.dxz = dxz_;
  v.dyz = dyz_;
  v.ax = subs_.ax;
  v.ay = subs_.ay;
  v.cosx = subs_.cosx;
  v.cosy = subs_.cosy;
  v.dz3d = subs_.dz3d;
  v.kx = subs_.kx;
  v.ky = subs_.ky;
  v.kzLower = kzLower_;
  v.gid = grid_.gid3();
  v.kTop = kTop_;
  v.sSurf = coupled ? sSurf_ : Field2<real_t>();
  v.az = subs_.az;
  v.dx = grid_.dx();
  v.dy = grid_.dy();
  v.nz = nz;
  v.coupled = coupled;

  Field3<real_t> s = sSubs_, sm = smSubs_, sMin = sMin3_, sMax = sMax3_;
  Field3<real_t> qx = subs_.qx, qy = subs_.qy, qzF = subs_.qzF;
  Field3<real_t> kx = subs_.kx, ky = subs_.ky, kzLower = kzLower_, kzFace = subs_.kzF;
  Field3<real_t> wcn = subs_.wcn, wcs = subs_.wcs;
  Field3<real_t> dz3d = subs_.dz3d, ax = subs_.ax, ay = subs_.ay;
  Field3<real_t> fxm = sFarXm3_, fxp = sFarXp3_, fym = sFarYm3_, fyp = sFarYp3_;
  Field3<PetscInt> gid = grid_.gid3();
  Field2<int> kTop = kTop_;
  Field2<int> sideYp = subs_.sideCodeYp;
  Field2<int> topCode = subs_.topCode;
  Field2<real_t> topValue = subs_.topValue;
  Field2<real_t> ssee = coupled ? sseepage_ : Field2<real_t>();
  Field2<real_t> dept = coupled ? surf_.dept : Field2<real_t>();
  Field2<real_t> sSurf = coupled ? sSurf_ : Field2<real_t>();
  const real_t az = subs_.az;
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t molecular = dispMol_;
  const real_t lon = dispLon_;

  // Increment pass: scalar mass, the advective and dispersive face fluxes,
  // the surface exchange, and the limiter's local extrema (scalar.c:318-435).
  real_t boundaryMass = 0.0;
  Kokkos::parallel_reduce(
      "transport_subsurface_increment",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sumBoundary) {
        if (gid(j, i, k) < 0) {
          sMin(j, i, k) = limHi;
          sMax(j, i, k) = limLo;
          return;
        }
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        const bool westCell = (gi == 0);
        const bool eastCell = (gi == nxGlobal - 1);
        const bool southCell = (gj == 0);
        const bool northCell = (gj == nyGlobal - 1);
        const int columnTop = kTop(j, i);
        const bool isTop = (k == columnTop);
        const real_t volume = az * dz3d(j, i, k);

        real_t mass = s(j, i, k) * wcn(j, i, k) * volume;

        // Far neighbors for the superbee stencils.
        const real_t sIp2 = (i + 2 <= nxl + 1) ? s(j, i + 2, k) : fxp(j, i + 1, k);
        const real_t sIm2 = (i - 2 >= 0) ? s(j, i - 2, k) : fxm(j, i - 1, k);
        const real_t sJp2 = (j + 2 <= nyl + 1) ? s(j + 2, i, k) : fyp(j + 1, i, k);
        const real_t sJm2 = (j - 2 >= 0) ? s(j - 2, i, k) : fym(j - 1, i, k);

        // x advection (advective_flux :505-562). Positive qx points toward
        // -x: the x+ face flux is an inflow from i+1 when positive.
        const real_t qxp = qx(j, i, k);
        const real_t qxm = qx(j, i - 1, k);
        real_t sip = 0.0;
        if (superbee) {
          if (qxp > 0.0) {
            sip = eastCell ? s(j, i + 1, k)
                           : tvdSuperbee(s(j, i, k), s(j, i + 1, k), sIp2,
                                         qxp / ax(j, i, k), dx, dt);
          } else if (qxp < 0.0) {
            sip = tvdSuperbee(s(j, i + 1, k), s(j, i, k), s(j, i - 1, k), qxp / ax(j, i, k),
                              dx, dt);
          }
        } else {
          if (qxp < 0.0) {
            sip = s(j, i, k);
          } else if (qxp > 0.0) {
            sip = s(j, i + 1, k);
          }
        }
        real_t sim = 0.0;
        if (superbee) {
          if (qxm > 0.0) {
            sim = tvdSuperbee(s(j, i - 1, k), s(j, i, k), s(j, i + 1, k),
                              qxm / ax(j, i - 1, k), dx, dt);
          } else if (qxm < 0.0) {
            sim = westCell ? s(j, i - 1, k)
                           : tvdSuperbee(s(j, i, k), s(j, i - 1, k), sIm2,
                                         qxm / ax(j, i - 1, k), dx, dt);
          }
        } else {
          if (qxm < 0.0) {
            sim = s(j, i - 1, k);
          } else if (qxm > 0.0) {
            sim = s(j, i, k);
          }
        }

        // y advection (:564-623).
        const real_t qyp = qy(j, i, k);
        const real_t qym = qy(j - 1, i, k);
        real_t sjp = 0.0;
        if (superbee) {
          if (qyp > 0.0) {
            sjp = northCell ? s(j + 1, i, k)
                            : tvdSuperbee(s(j, i, k), s(j + 1, i, k), sJp2,
                                          qyp / ay(j, i, k), dy, dt);
          } else if (qyp < 0.0) {
            sjp = tvdSuperbee(s(j + 1, i, k), s(j, i, k), s(j - 1, i, k), qyp / ay(j, i, k),
                              dy, dt);
          }
        } else {
          if (qyp < 0.0) {
            sjp = s(j, i, k);
          } else if (qyp > 0.0) {
            sjp = s(j + 1, i, k);
          }
        }
        real_t sjm = 0.0;
        if (superbee) {
          if (qym > 0.0) {
            sjm = tvdSuperbee(s(j - 1, i, k), s(j, i, k), s(j + 1, i, k),
                              qym / ay(j - 1, i, k), dy, dt);
          } else if (qym < 0.0) {
            sjm = southCell ? s(j - 1, i, k)
                            : tvdSuperbee(s(j, i, k), s(j - 1, i, k), sJm2,
                                          qym / ay(j - 1, i, k), dy, dt);
          }
        } else {
          if (qym < 0.0) {
            sjm = s(j - 1, i, k);
          } else if (qym > 0.0) {
            sjm = s(j, i, k);
          }
        }

        // z advection (:625-711). The face below carries qzF(k+1), the face
        // above qzF(k); positive flux points upward.
        const real_t qzp = qzF(j, i, k + 1);
        const real_t qzm = qzF(j, i, k);
        const real_t sBelow1 = scalarBelow(v, j, i, k);
        real_t skp = 0.0;
        if (superbee) {
          if (qzp > 0.0) {
            if (k == nz - 1) {
              skp = sBelow1;
            } else {
              const real_t sBelow2 = scalarBelow(v, j, i, k + 1);
              skp = tvdSuperbee(s(j, i, k), s(j, i, k + 1), sBelow2, qzp / az,
                                dz3d(j, i, k), dt);
            }
          } else if (qzp < 0.0) {
            skp = tvdSuperbee(sBelow1, s(j, i, k), scalarAbove(v, j, i, k), qzp / az,
                              dz3d(j, i, k), dt);
          }
        } else {
          if (qzp > 0.0) {
            skp = sBelow1;
          } else if (qzp < 0.0) {
            skp = s(j, i, k);
          }
        }
        real_t skm = 0.0;
        if (isTop) {
          // Top face: coupled exchange advects through the seepage term
          // instead (scalar.c:700-702); uncoupled runs keep the legacy
          // head-condition donor value and zero flux-condition faces
          // (:672-697 with the branch outcomes derived in the file comment).
          const bool topFlux = (topCode(j, i) == kBcFlux);
          if (coupled) {
            skm = 0.0;
          } else if (superbee) {
            skm = topFlux ? 0.0 : s(j, i, k);
          } else {
            skm = (qzm > 0.0 && !topFlux) ? s(j, i, k) : 0.0;
          }
        } else if (superbee) {
          if (qzm > 0.0) {
            skm = tvdSuperbee(s(j, i, k - 1), s(j, i, k), sBelow1, qzm / az, dz3d(j, i, k),
                              dt);
          } else if (qzm < 0.0) {
            skm = tvdSuperbee(s(j, i, k), s(j, i, k - 1), scalarAbove(v, j, i, k - 1),
                              qzm / az, dz3d(j, i, k), dt);
          }
        } else {
          if (qzm > 0.0) {
            skm = s(j, i, k);
          } else if (qzm < 0.0) {
            skm = s(j, i, k - 1);
          }
        }

        mass += dt * (qxp * sip - qxm * sim + qyp * sjp - qym * sjm + qzp * skp - qzm * skm);

        // Dispersive fluxes with the live edge doublings (scalar.c:343-365).
        const real_t jip = dispFluxX(v, j, i, k);
        const real_t jim = dispFluxX(v, j, i - 1, k);
        real_t jjp = dispFluxY(v, j, i, k);
        if (northCell) {
          jjp *= 2.0;
        }
        real_t jjm = dispFluxY(v, j - 1, i, k);
        if (southCell) {
          jjm *= 2.0;
        }
        const real_t jkp = dispFluxZLower(v, j, i, k);
        real_t jkm = 0.0;
        if (isTop) {
          // The surface-interface dispersive flux (dispersive_flux :834-845,
          // the ghost-cell else branch): live only under a wet coupled
          // surface. The ghost coefficient is the legacy tensor at the kM
          // ghost, whose only nonzero flux slot is the top exchange face —
          // molecular part plus the longitudinal term of |q_top|.
          if (coupled && dept(j, i) > 0.0) {
            const real_t dzzGhost =
                molecular * wcs(j, i, k) + lon * Kokkos::fabs(qzF(j, i, k));
            jkm = dzzGhost * az * (s(j, i, k) - sSurf(j, i)) / dz3d(j, i, k);
          }
        } else {
          jkm = dispFluxZLower(v, j, i, k - 1);
        }

        mass += dt * ((jip - jim) + (jjp - jjm) + (jkp - jkm));

        // Surface exchange (scalar.c:371-375).
        if (coupled && isTop) {
          mass -= dt * az * ssee(j, i);
        }

        sm(j, i, k) = mass;

        // Local extrema of the conductive stencil (scalar.c:378-434).
        real_t lo = limHi;
        real_t hi = limLo;
        const auto extend = [&](real_t value) {
          if (value > hi) {
            hi = value;
          }
          if (value < lo) {
            lo = value;
          }
        };
        if (kx(j, i, k) > 0.0 && gid(j, i + 1, k) >= 0) {
          extend(s(j, i + 1, k));
        }
        if (kx(j, i - 1, k) > 0.0 && gid(j, i - 1, k) >= 0) {
          extend(s(j, i - 1, k));
        }
        if (ky(j, i, k) > 0.0) {
          if (gid(j + 1, i, k) >= 0) {
            extend(s(j + 1, i, k));
          } else if (northCell && sideYp(j, i) != kBcNoFlux) {
            // The prescribed side ghost feeds the limiter (scalar.c:397-400)
            // — the sea-side Dirichlet salinity may raise the local bound.
            extend(s(j + 1, i, k));
          }
        }
        if (ky(j - 1, i, k) > 0.0 && gid(j - 1, i, k) >= 0) {
          extend(s(j - 1, i, k));
        }
        if (kzLower(j, i, k) > 0.0 && k + 1 < nz && gid(j, i, k + 1) >= 0) {
          extend(s(j, i, k + 1));
        }
        // Upper-face extension (scalar.c:412-434): the cell above, or the
        // surface scalar at a wet coupled top.
        if (kzFace(j, i, k) > 0.0) {
          if (isTop) {
            if (coupled && dept(j, i) > 0.0) {
              extend(sSurf(j, i));
            }
          } else if (gid(j, i, k - 1) >= 0) {
            extend(s(j, i, k - 1));
          }
        }
        sMin(j, i, k) = lo;
        sMax(j, i, k) = hi;

        // Boundary bookkeeping: scalar mass through global-edge, top, and
        // bottom faces (positive into the domain).
        real_t in = 0.0;
        if (westCell) {
          in -= dt * (qxm * sim + jim);
        }
        if (eastCell) {
          in += dt * (qxp * sip + jip);
        }
        if (southCell) {
          in -= dt * (qym * sjm + jjm);
        }
        if (northCell) {
          in += dt * (qyp * sjp + jjp);
        }
        if (k == nz - 1) {
          in += dt * (qzp * skp + jkp);
        }
        if (isTop && !coupled) {
          in -= dt * qzm * skm;
        }
        if (isTop && coupled) {
          // The one-sided interface dispersive gain (see TransportAudit).
          in -= dt * jkm;
        }
        sumBoundary += in;
      },
      boundaryMass);
  audit_.subsBoundary = boundaryMass;

  // Update pass (scalar.c:437-485): concentration from the flux volume
  // (volume_by_flux_subs, groundwater.c:1634-1644, evaluated over the last
  // substep), the qtop quirk, the limiter, and the bounds.
  Field3<real_t> wcNow = subs_.wc;
  real_t adjustMass = 0.0;
  real_t anchorMass = 0.0;
  long negatives = 0;
  Kokkos::parallel_reduce(
      "transport_subsurface_update",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k, real_t& sumAdjust,
                    real_t& sumAnchor, long& badCount) {
        if (gid(j, i, k) < 0) {
          s(j, i, k) = 0.0;
          return;
        }
        const real_t volume = az * dz3d(j, i, k);
        const real_t vgflux =
            wcn(j, i, k) * volume +
            dtgLast * (qx(j, i, k) - qx(j, i - 1, k) + qy(j, i, k) - qy(j - 1, i, k) +
                       qzF(j, i, k + 1) - qzF(j, i, k));
        real_t value = (vgflux > 0.0) ? sm(j, i, k) / vgflux : 0.0;
        const real_t raw = value;

        real_t hi = sMax(j, i, k);
        if (coupled && k == kTop(j, i) && topCode(j, i) == kBcFlux &&
            topValue(j, i) > 0.0 && dept(j, i) <= 0.0) {
          hi += 0.01;  // evaporative-concentration allowance (scalar.c:452-458)
        }
        if (value > hi && hi < limHi) {
          value = hi;
        } else if (value < sMin(j, i, k) && sMin(j, i, k) > limLo) {
          value = sMin(j, i, k);
        }
        if (hasMax && value > limHi) {
          value = limHi;
        } else if (value < limLo) {
          if (value < limLo - 0.01) {
            ++badCount;
          }
          value = limLo;
        }
        if (vgflux > 0.0) {
          sumAdjust += (value - raw) * vgflux;
          // Ledger re-anchor: the mass now lives on θ_new V, which differs
          // from the flux volume by the reallocation/clamp adjustments.
          sumAnchor += value * (wcNow(j, i, k) * volume - vgflux);
        } else {
          sumAdjust -= sm(j, i, k);
          sumAnchor += value * wcNow(j, i, k) * volume;
        }
        s(j, i, k) = value;
      },
      adjustMass, anchorMass, negatives);
  audit_.subsAdjust = adjustMass;
  audit_.subsAnchor = anchorMass;
  if (negatives > 0) {
    log::warn(log::msg() << "transport: " << negatives
                         << " subsurface cell(s) clipped from below the scalar bound");
  }

  enforceSubsurfaceBc(t);

  // Snapshot the top-cell Dzz for the next surface step's exchange term
  // (the tensor is rebuilt each step; see the file comment).
  Field3<real_t> dzz = dzz_;
  Field2<real_t> dzzTop = dzzTop_;
  Kokkos::parallel_for(
      "transport_dzz_top",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        const int columnTop = kTop(j, i);
        dzzTop(j, i) = (columnTop >= 0) ? dzz(j, i, columnTop) : 0.0;
      });

  halo_.exchangeWithCorners({"s_subs"});
}

void ScalarSolver::enforceSubsurfaceBc(real_t t) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const bool west = (grid_.rankWest() == MPI_PROC_NULL);
  const bool east = (grid_.rankEast() == MPI_PROC_NULL);
  const bool south = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool north = (grid_.rankNorth() == MPI_PROC_NULL);
  const bool coupled = surf_.active && cpl_.active;

  Field3<real_t> s = sSubs_;
  Field2<int> kTop = kTop_;
  Field2<real_t> sKp = sSurfKp_;

  // Zero-gradient domain-edge ghosts (enforce_scalar_bc, scalar.c:855-916).
  Kokkos::parallel_for(
      "transport_subs_ghost_x",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nyl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int k) {
        if (west) {
          s(j, 0, k) = s(j, 1, k);
        }
        if (east) {
          s(j, nxl + 1, k) = s(j, nxl, k);
        }
      });
  Kokkos::parallel_for(
      "transport_subs_ghost_y",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 0}, {nxl + 1, nz}),
      KOKKOS_LAMBDA(const int i, const int k) {
        if (south) {
          s(0, i, k) = s(1, i, k);
        }
        if (north) {
          s(nyl + 1, i, k) = s(nyl, i, k);
        }
      });

  // Prescribed side ghosts (the legacy s_yp/s_ym Dirichlet, :869-895,
  // generalized to all four sides): the configured value replaces the
  // zero-gradient copy along the member edge cells' ghost column.
  for (const ScalarBcList& list : sideGhost_) {
    const real_t value = list.bc->value(t);
    auto lj = list.j;
    auto li = list.i;
    auto lface = list.face;
    Kokkos::parallel_for(
        "transport_subs_side_dirichlet",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {0, 0}, {static_cast<int>(lj.extent(0)), nz}),
        KOKKOS_LAMBDA(const int m, const int k) {
          const int j = lj(static_cast<std::size_t>(m));
          const int i = li(static_cast<std::size_t>(m));
          switch (static_cast<BcFace>(lface(static_cast<std::size_t>(m)))) {
            case BcFace::XMinus:
              s(j, i - 1, k) = value;
              break;
            case BcFace::XPlus:
              s(j, i + 1, k) = value;
              break;
            case BcFace::YMinus:
              s(j - 1, i, k) = value;
              break;
            case BcFace::YPlus:
              s(j + 1, i, k) = value;
              break;
            default:
              break;
          }
        });
  }

  // Export the top-cell scalar for the next surface exchange (scalar.c:903).
  if (coupled) {
    Kokkos::parallel_for(
        "transport_subs_export_top",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i) {
          const int columnTop = kTop(j, i);
          sKp(j, i) = (columnTop >= 0) ? s(j, i, columnTop) : 0.0;
        });
  }
}

}  // namespace frehg::transport
