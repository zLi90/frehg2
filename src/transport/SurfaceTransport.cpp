/// \file SurfaceTransport.cpp
/// \brief Surface scalar transport (scalar_shallowwater, scalar.c:25-298).
///
/// The step advances the scalar mass s Vsn by the explicit advective
/// (current-step flow rates), diffusive, exchange, and inflow increments,
/// divides by the flux volume (volume_by_flux, shallowwater.c:996-1032 —
/// evaluated with the *previous* step's flow rates exactly as legacy calls
/// it before the velocity update), applies the local min/max limiter and
/// the configured bounds, dilutes by the rain/evaporation volume, and
/// enforces the ghost/tide/dry rules in the legacy order.
///
/// Deviations, all recorded in docs/theory/transport.md:
///  - the exchange-diffusion denominator is the top cell's thickness; the
///    legacy smap->dz it divided by is allocated and never written
///    (map.c:49) — an uninitialized read of the §2.1 hazard class;
///  - the flux-volume exchange term drops the legacy porosity factor
///    (qss dt wcs Asz, shallowwater.c:1008-1015) per the plan §5.7/A9
///    resolution — no porosity factor anywhere in the exchange;
///  - rain dilution uses the per-cell masked rain (the legacy global rate
///    ignored its own exclusion rows);
///  - the superbee far-neighbor guards are evaluated against the *global*
///    domain edges with staged far values at rank interfaces (legacy
///    guarded rank-local edges, degrading interface faces to upwind — a
///    rank-count-dependent stencil; the serial goldens never see it).

#include "transport/Limiters.hpp"
#include "transport/ScalarSolver.hpp"

namespace frehg::transport {

void ScalarSolver::stageSurfaceFarNeighbors() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  Field2<real_t> s = sSurf_;
  Field2<real_t> fxm = sFarXm2_, fxp = sFarXp2_, fym = sFarYm2_, fyp = sFarYp2_;
  Kokkos::parallel_for(
      "transport_far_2d",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        fxm(j, i) = s(j, i - 1);
        fxp(j, i) = s(j, i + 1);
        fym(j, i) = s(j - 1, i);
        fyp(j, i) = s(j + 1, i);
      });
  halo_.exchange({"s_far_xm_2d", "s_far_xp_2d", "s_far_ym_2d", "s_far_yp_2d"});
}

void ScalarSolver::stepSurface(real_t t, real_t dt, real_t rain, real_t evap) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const int nxGlobal = grid_.nx();
  const int nyGlobal = grid_.ny();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t cellArea = dx * dy;
  const bool superbee = superbee_;
  const bool coupled = subs_.active && cpl_.active;
  const real_t minDepth = surf_.minDepth;
  const real_t limHi = boundMax_;
  const real_t limLo = boundMin_;
  const bool hasMax = hasBoundMax_;
  const real_t difux = difuX_, difuy = difuY_;

  stageSurfaceFarNeighbors();

  Field2<real_t> s = sSurf_, sm = smSurf_, sMin = sMinS_, sMax = sMaxS_;
  Field2<real_t> vsn = vsn_, vflux = vflux_, fuOld = fuOld_, fvOld = fvOld_;
  Field2<real_t> sKp = sSurfKp_, ssee = sseepage_, dzzTop = dzzTop_;
  Field2<real_t> fxm = sFarXm2_, fxp = sFarXp2_, fym = sFarYm2_, fyp = sFarYp2_;
  Field2<real_t> dept = surf_.dept, etan = surf_.etan, bottom = surf_.bottom;
  Field2<real_t> uu = surf_.uu, vv = surf_.vv, Fu = surf_.fu, Fv = surf_.fv;
  Field2<real_t> Asx = surf_.asx, Asy = surf_.asy, rainMask = surf_.rainMask;
  Field2<real_t> qss = coupled ? cpl_.qss : Field2<real_t>();
  Field3<real_t> ksz = coupled ? subs_.ksz : Field3<real_t>();
  Field3<real_t> dz3d = coupled ? subs_.dz3d : Field3<real_t>();
  Field3<PetscInt> gid = grid_.gid3();

  // Increment pass (scalar.c:38-196): scalar mass, advection with the
  // current flow rates, diffusion, the seepage exchange, the limiter's
  // local extrema, and the legacy flux volume from the previous step's
  // flow rates.
  real_t exchangeMass = 0.0;
  real_t boundaryMass = 0.0;
  Kokkos::parallel_reduce(
      "transport_surface_increment",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumExchange, real_t& sumBoundary) {
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        const bool westCell = (gi == 0);
        const bool eastCell = (gi == nxGlobal - 1);
        const bool southCell = (gj == 0);
        const bool northCell = (gj == nyGlobal - 1);
        // Far neighbors (in-field where the halo covers them, staged copies
        // at the interface slots).
        const real_t sIp2 = (i + 2 <= nxl + 1) ? s(j, i + 2) : fxp(j, i + 1);
        const real_t sIm2 = (i - 2 >= 0) ? s(j, i - 2) : fxm(j, i - 1);
        const real_t sJp2 = (j + 2 <= nyl + 1) ? s(j + 2, i) : fyp(j + 1, i);
        const real_t sJm2 = (j - 2 >= 0) ? s(j - 2, i) : fym(j - 1, i);

        real_t mass = s(j, i) * vsn(j, i);

        // x advection (scalar.c:45-77).
        real_t sip;
        if (Fu(j, i) > 0.0) {
          sip = s(j, i);
          if (superbee) {
            sip = tvdSuperbee(s(j, i + 1), s(j, i), s(j, i - 1), uu(j, i), dx, dt);
          }
        } else {
          sip = s(j, i + 1);
          if (superbee && !eastCell) {
            sip = tvdSuperbee(s(j, i), s(j, i + 1), sIp2, uu(j, i), dx, dt);
          }
        }
        real_t sim;
        if (Fu(j, i - 1) > 0.0) {
          sim = s(j, i - 1);
          if (superbee && !westCell) {
            sim = tvdSuperbee(s(j, i), s(j, i - 1), sIm2, uu(j, i - 1), dx, dt);
          }
        } else {
          sim = s(j, i);
          if (superbee) {
            sim = tvdSuperbee(s(j, i - 1), s(j, i), s(j, i + 1), uu(j, i - 1), dx, dt);
          }
        }
        // y advection (scalar.c:78-110).
        real_t sjp;
        if (Fv(j, i) > 0.0) {
          sjp = s(j, i);
          if (superbee) {
            sjp = tvdSuperbee(s(j + 1, i), s(j, i), s(j - 1, i), vv(j, i), dy, dt);
          }
        } else {
          sjp = s(j + 1, i);
          if (superbee && !northCell) {
            sjp = tvdSuperbee(s(j, i), s(j + 1, i), sJp2, vv(j, i), dy, dt);
          }
        }
        real_t sjm;
        if (Fv(j - 1, i) > 0.0) {
          sjm = s(j - 1, i);
          if (superbee && !southCell) {
            sjm = tvdSuperbee(s(j, i), s(j - 1, i), sJm2, vv(j - 1, i), dy, dt);
          }
        } else {
          sjm = s(j, i);
          if (superbee) {
            sjm = tvdSuperbee(s(j - 1, i), s(j, i), s(j + 1, i), vv(j - 1, i), dy, dt);
          }
        }
        mass += dt * (-Fu(j, i) * sip + Fu(j, i - 1) * sim - Fv(j, i) * sjp +
                      Fv(j - 1, i) * sjm);

        // Advective mass through the global-edge faces (see TransportAudit):
        // the closed-edge fold keeps the water in, but the boundary-face
        // velocities are live and the scalar crosses them.
        if (westCell) {
          sumBoundary += dt * Fu(j, i - 1) * sim;
        }
        if (eastCell) {
          sumBoundary -= dt * Fu(j, i) * sip;
        }
        if (southCell) {
          sumBoundary += dt * Fv(j - 1, i) * sjm;
        }
        if (northCell) {
          sumBoundary -= dt * Fv(j, i) * sjp;
        }

        // Diffusion (scalar.c:114-119).
        mass += dt * ((difux * Asx(j, i) / dx) * (s(j, i + 1) - s(j, i)) -
                      (difux * Asx(j, i - 1) / dx) * (s(j, i) - s(j, i - 1)) +
                      (difuy * Asy(j, i) / dy) * (s(j + 1, i) - s(j, i)) -
                      (difuy * Asy(j - 1, i) / dy) * (s(j, i) - s(j - 1, i)));

        // Subsurface exchange (scalar.c:121-152): the seepage carries the
        // upwind concentration plus a two-point diffusive term against the
        // top-cell scalar; scalar does not leave a dry surface.
        real_t topKsz = 0.0;
        if (coupled) {
          int kTop = -1;
          for (int k = 0; k < nz; ++k) {
            if (gid(j, i, k) >= 0) {
              kTop = k;
              break;
            }
          }
          real_t rate = 0.0;
          if (kTop >= 0) {
            topKsz = ksz(j, i, kTop);
            const real_t diffusive = 2.0 * dzzTop(j, i) / dz3d(j, i, kTop);
            if (qss(j, i) > 0.0) {
              if (dept(j, i) > 0.0) {
                rate = qss(j, i) * sKp(j, i) + diffusive * (sKp(j, i) - s(j, i));
              }
            } else if (qss(j, i) < 0.0) {
              rate = qss(j, i) * s(j, i) + diffusive * (sKp(j, i) - s(j, i));
            }
          }
          ssee(j, i) = rate;
          mass += rate * cellArea * dt;
          sumExchange += rate * cellArea * dt;
        } else {
          ssee(j, i) = 0.0;
        }

        sm(j, i) = mass;

        // Local extrema of the wet stencil (scalar.c:154-179): neighbors
        // with a live face and a wet start-of-step stage, plus the top-cell
        // scalar where the exchange is conductive.
        real_t lo = s(j, i);
        real_t hi = s(j, i);
        const auto extend = [&](real_t v) {
          if (v > hi) {
            hi = v;
          }
          if (v < lo) {
            lo = v;
          }
        };
        if (Asx(j, i) > 0.0 && etan(j, i + 1) > bottom(j, i + 1)) {
          extend(s(j, i + 1));
        }
        if (Asx(j, i - 1) > 0.0 && etan(j, i - 1) > bottom(j, i - 1)) {
          extend(s(j, i - 1));
        }
        if (Asy(j, i) > 0.0 && etan(j + 1, i) > bottom(j + 1, i)) {
          extend(s(j + 1, i));
        }
        if (Asy(j - 1, i) > 0.0 && etan(j - 1, i) > bottom(j - 1, i)) {
          extend(s(j - 1, i));
        }
        if (coupled && topKsz > 0.0) {
          extend(sKp(j, i));
        }
        sMin(j, i) = lo;
        sMax(j, i) = hi;

        // Flux volume (volume_by_flux, shallowwater.c:996-1032) from the
        // previous step's flow rates; the coupled term is the applied
        // exchange volume without the legacy porosity factor (A9).
        real_t v = vsn(j, i) +
                   dt * (fuOld(j, i - 1) - fuOld(j, i) + fvOld(j - 1, i) - fvOld(j, i));
        if (coupled) {
          const real_t area = (dept(j, i) > 0.0) ? cellArea : 0.0;
          const real_t dv = qss(j, i) * dt * area;
          if (qss(j, i) < 0.0 && v > -dv) {
            v += dv;
          } else if (qss(j, i) > 0.0) {
            if (dept(j, i) > 0.0) {
              v += dv;
            }
          }
        }
        vflux(j, i) = v;
      },
      exchangeMass, boundaryMass);
  audit_.exchange = exchangeMass;
  audit_.surfBoundary = boundaryMass;

  // Inflow sources (scalar.c:180-195): mass injected at discharge-condition
  // cells with the paired concentration; the limiter is released there.
  real_t inflowMass = 0.0;
  for (const ScalarBcList& list : surfaceInflow_) {
    const real_t q = list.discharge->value(t);
    if (q <= 0.0) {
      continue;
    }
    const real_t sIn = list.bc->value(t);
    const real_t share = q * dt * sIn / static_cast<real_t>(list.dischargeCells);
    const real_t qShare = q * dt / static_cast<real_t>(list.dischargeCells);
    auto lj = list.j;
    auto li = list.i;
    real_t added = 0.0;
    Kokkos::parallel_reduce(
        "transport_surface_inflow", Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, lj.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m, real_t& sum) {
          const int j = lj(m);
          const int i = li(m);
          sm(j, i) += share;
          // The discharge volume itself joins the flux volume
          // (shallowwater.c:1020-1030).
          vflux(j, i) += qShare;
          sMax(j, i) = limHi;
          sMin(j, i) = limLo - 1.0;
          sum += share;
        },
        added);
    inflowMass += added;
  }
  audit_.surfSource = inflowMass;

  // Update pass (scalar.c:197-257): concentration from the flux volume,
  // the limiter with its sentinel guards, the configured bounds, and the
  // rain/evaporation dilution.
  const real_t rainRate = rain;
  const real_t evapRate = evap;
  real_t adjustMass = 0.0;
  real_t anchorMass = 0.0;
  Kokkos::parallel_reduce(
      "transport_surface_update",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumAdjust, real_t& sumAnchor) {
        const bool wet = (vflux(j, i) > 0.0 && dept(j, i) > 0.0);
        real_t value = wet ? sm(j, i) / vflux(j, i) : 0.0;
        const real_t raw = value;
        if (wet) {
          const bool moving = (uu(j, i) != 0.0 || uu(j, i - 1) != 0.0 || vv(j, i) != 0.0 ||
                               vv(j - 1, i) != 0.0);
          bool applyLimiter = moving;
          if (!moving && coupled) {
            // Legacy gates the still-water clamp on a conductive top
            // (param->Ksz > 0, scalar.c:217); the staged extrema already
            // encode the per-column conductivity, so reuse the wet gate.
            int kTop = -1;
            for (int k = 0; k < nz; ++k) {
              if (gid(j, i, k) >= 0) {
                kTop = k;
                break;
              }
            }
            applyLimiter = (kTop >= 0 && ksz(j, i, kTop) > 0.0);
          }
          if (applyLimiter) {
            if (value > sMax(j, i) && sMax(j, i) != limHi && sMax(j, i) != limLo) {
              value = sMax(j, i);
            } else if (value < sMin(j, i) && sMin(j, i) != limHi) {
              value = sMin(j, i);
            }
          }
          if (hasMax && value >= limHi) {
            value = limHi;
          } else if (value < limLo) {
            value = limLo;
          }
          sumAdjust += (value - raw) * vflux(j, i);
        } else {
          sumAdjust -= sm(j, i);
        }
        // Rain/evaporation dilution (scalar.c:243-254).
        real_t vre = 0.0;
        if (vflux(j, i) > 0.0) {
          const real_t area = (dept(j, i) > 0.0) ? cellArea : 0.0;
          const real_t cellRain = rainRate * rainMask(j, i);
          if (dept(j, i) > minDepth && cellRain > 0.0) {
            vre += cellRain * area * dt;
          } else {
            vre -= evapRate * area * dt;
          }
          value = value * vflux(j, i) / (vflux(j, i) + vre);
        }
        s(j, i) = value;
        // Ledger re-anchor: the mass now lives on the actual volume, not
        // the (lagged) flux volume the update divided by.
        if (wet) {
          sumAnchor += value * (dept(j, i) * cellArea - (vflux(j, i) + vre));
        }
      },
      adjustMass, anchorMass);
  audit_.surfAdjust = adjustMass;
  audit_.surfAnchor = anchorMass;

  // Ghost, tide, and dry rules in the legacy order (scalar.c:263-288):
  // zero-gradient edge ghosts, then the Dirichlet stage salinity, then the
  // dry zeroing.
  const bool west = (grid_.rankWest() == MPI_PROC_NULL);
  const bool east = (grid_.rankEast() == MPI_PROC_NULL);
  const bool south = (grid_.rankSouth() == MPI_PROC_NULL);
  const bool north = (grid_.rankNorth() == MPI_PROC_NULL);
  Kokkos::parallel_for(
      "transport_surface_ghost_x", Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nyl + 1),
      KOKKOS_LAMBDA(const int j) {
        if (west) {
          s(j, 0) = s(j, 1);
        }
        if (east) {
          s(j, nxl + 1) = s(j, nxl);
        }
      });
  Kokkos::parallel_for(
      "transport_surface_ghost_y", Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>(1, nxl + 1),
      KOKKOS_LAMBDA(const int i) {
        if (south) {
          s(0, i) = s(1, i);
        }
        if (north) {
          s(nyl + 1, i) = s(nyl, i);
        }
      });

  real_t resetMass = 0.0;
  for (const ScalarBcList& list : surfaceDirichlet_) {
    const real_t value = list.bc->value(t);
    auto lj = list.j;
    auto li = list.i;
    real_t delta = 0.0;
    Kokkos::parallel_reduce(
        "transport_surface_dirichlet",
        Kokkos::RangePolicy<ExecSpace>(std::size_t{0}, lj.extent(0)),
        KOKKOS_LAMBDA(const std::size_t m, real_t& sum) {
          const int j = lj(m);
          const int i = li(m);
          sum += (value - s(j, i)) * dept(j, i) * cellArea;
          s(j, i) = value;
        },
        delta);
    resetMass += delta;
  }
  audit_.surfSource += resetMass;

  Kokkos::parallel_for(
      "transport_surface_dry",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i) {
        if (dept(j, i) <= 0.0) {
          s(j, i) = 0.0;
        }
      });

  halo_.exchange({"s_surf"});
}

}  // namespace frehg::transport
