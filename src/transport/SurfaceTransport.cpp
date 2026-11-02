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
///
/// Opt-in correction (transport.surface_update: consistent): the flux
/// volume uses the current step's flow rates, matching the advective
/// increments, and the local min/max limiter applies in every wet cell.
/// The legacy pairing (lagged volume, moving-cell-only limiter) concentrates
/// thin films whose outflow the wetting limiter has just closed above every
/// input value and books the volume mismatch into surf_anchor; under
/// strongly unsteady wetting/drying that term can rival the scalar mass.

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
  const std::string p = spec_.prefix + "_";
  halo_.exchange({p + "far_xm_2d", p + "far_xp_2d", p + "far_ym_2d", p + "far_yp_2d"});
}

void ScalarSolver::stepSurface(real_t t, real_t dt, real_t rain, real_t evap) {
  if (conservativeSurface_) {
    stepSurfaceImpl<true>(t, dt, rain, evap);
  } else {
    stepSurfaceImpl<false>(t, dt, rain, evap);
  }
}

template <bool QVOL>
void ScalarSolver::stepSurfaceImpl(real_t t, real_t dt, real_t rain, real_t evap) {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const Grid::EdgeTopology topo = grid_.edgeTopology();
  const int nz = grid_.nz();
  const int i0 = grid_.i0();
  const int j0 = grid_.j0();
  const real_t dx = grid_.dx();
  const real_t dy = grid_.dy();
  const real_t cellArea = dx * dy;
  const bool superbee = superbee_;
  const bool coupled = subs_.active && cpl_.active;
  const bool consistent = consistentSurface_;
  // transport.subsurface_update: consistent — the exchange side of the
  // consistent subsurface update (see the exchange block below).
  const bool consistentExchange = coupled && consistentSubsurface_;
  const real_t minDepth = surf_.minDepth;
  const real_t limHi = boundMax_;
  const real_t limLo = boundMin_;
  const bool hasMax = hasBoundMax_;
  const bool hasMin = hasBoundMin_;
  const real_t difux = difuX_, difuy = difuY_;
  // v2 Q5 thermal switches (V2-A17): temperature skips the rain/evap
  // dilution (rain enters and evaporation leaves at the cell's own
  // temperature — the volume change carries no temperature change; the
  // rain-heat approximation lands in the anchor column), and may carry a
  // surface heat-exchange source, applied inside the update pass on the
  // final volume basis and audited into surf_atmos.
  const bool dilution = !spec_.isTemperature;
  const bool atmosEquilibrium =
      spec_.exchange == SurfaceExchangeConfig::Mode::Equilibrium;
  const bool atmosBulk = spec_.exchange == SurfaceExchangeConfig::Mode::Bulk;
  const real_t equilibriumT = spec_.equilibriumT;
  const real_t equilibriumK = spec_.equilibriumK;
  const real_t rhoCw = spec_.heatCapacityWater;
  atm::MetSample metNow;
  if (atmosBulk) {
    metNow = met_.sample(t);
  }
  const real_t metAirT = metNow.airTemperatureC;
  const real_t metPressure = metNow.pressureKpa;
  const real_t metWind = metNow.windSpeed;
  const real_t metQAir = metNow.airSpecificHumidity;
  const real_t metShortwave = metNow.shortwave;
  const real_t metLongwaveIn = metNow.longwaveIn;

  stageSurfaceFarNeighbors();

  Field2<real_t> s = sSurf_, sm = smSurf_, sMin = sMinS_, sMax = sMaxS_;
  Field2<real_t> vsn = vsn_, vflux = vflux_, fuOld = fuOld_, fvOld = fvOld_;
  Field2<real_t> sKp = sSurfKp_, ssee = sseepage_, dzzTop = dzzTop_;
  Field2<real_t> fxm = sFarXm2_, fxp = sFarXp2_, fym = sFarYm2_, fyp = sFarYp2_;
  Field2<real_t> dept = surf_.dept, etan = surf_.etan, bottom = surf_.bottom;
  Field2<real_t> uu = surf_.uu, vv = surf_.vv, Fu = surf_.fu, Fv = surf_.fv;
  Field2<real_t> Asx = surf_.asx, Asy = surf_.asy, rainMask = surf_.rainMask;
  Field2<real_t> qss = coupled ? cpl_.qss : Field2<real_t>();
  Field2<real_t> vent = consistentExchange ? cpl_.ventVolume : Field2<real_t>();
  // Stage-3 T10 (plan §12B C3): the coupler's exact applied-deposit book —
  // used for the flux volume under QVOL in place of the qss re-derivation.
  Field2<real_t> cplGain = (QVOL && coupled) ? cpl_.gain : Field2<real_t>();
  Field2<real_t> cplBounce = (QVOL && coupled) ? cpl_.bounce : Field2<real_t>();
  // Stage-3 T5: the salt store (W5.7), the ET per-cell dilution volumes
  // (D8), and the held-seepage salt split (W5.9).
  const bool store = store_;
  Field2<real_t> sStore = store ? sStore_ : Field2<real_t>();
  // Stage-3 T10 (W10.2, plan §12B C2): under surface_update: conservative
  // the advection, upwind signs, superbee weights and the flux volume all
  // come from the free-surface solve's own continuity face volumes — the
  // water that actually moved. legacy/consistent read Fu/uu through the
  // same locals with identical values (no arithmetic: bitwise-safe).
  Field2<real_t> qSolX = QVOL ? surf_.qxSolve : Field2<real_t>();
  Field2<real_t> qSolY = QVOL ? surf_.qySolve : Field2<real_t>();
  Field2<real_t> stageMask = QVOL ? surf_.isEtaBc : Field2<real_t>();
  const bool etSurf = surf_.etActive && !spec_.isTemperature;
  Field2<real_t> etRain = etSurf ? surf_.rainAdded : Field2<real_t>();
  Field2<real_t> etEvap = etSurf ? surf_.evapRemoved : Field2<real_t>();
  const bool held = coupled && consistentExchange && store && cpl_.heldSeepage;
  Field2<real_t> exfNew = held ? cpl_.exfNew : Field2<real_t>();
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
        // nvcc: an extended __host__ __device__ lambda may not FIRST capture
        // a variable inside an if-constexpr branch (CUDA C++ guide, "Extended
        // Lambda Restrictions"), so the QVOL branch's captures are named
        // here first. Discarded expressions: no code, no arithmetic.
        (void)etan;
        (void)bottom;
        (void)cplGain;
        (void)cplBounce;
        (void)stageMask;
        const int gi = i0 + i - 1;
        const int gj = j0 + j - 1;
        // Stage-3 T3: all false on a periodic axis (no edge: the superbee
        // far neighbors are real, no boundary mass is booked, W2.5's edge
        // diffusion has no face).
        const bool westCell = topo.west(gi);
        const bool eastCell = topo.east(gi);
        const bool southCell = topo.south(gj);
        const bool northCell = topo.north(gj);
        // Far neighbors (in-field where the halo covers them, staged copies
        // at the interface slots).
        const real_t sIp2 = (i + 2 <= nxl + 1) ? s(j, i + 2) : fxp(j, i + 1);
        const real_t sIm2 = (i - 2 >= 0) ? s(j, i - 2) : fxm(j, i - 1);
        const real_t sJp2 = (j + 2 <= nyl + 1) ? s(j + 2, i) : fyp(j + 1, i);
        const real_t sJm2 = (j - 2 >= 0) ? s(j - 2, i) : fym(j - 1, i);

        real_t mass = s(j, i) * vsn(j, i);

        const real_t flowE = QVOL ? qSolX(j, i) : Fu(j, i);
        const real_t flowW = QVOL ? qSolX(j, i - 1) : Fu(j, i - 1);
        const real_t flowN = QVOL ? qSolY(j, i) : Fv(j, i);
        const real_t flowS = QVOL ? qSolY(j - 1, i) : Fv(j - 1, i);
        const real_t velE =
            QVOL ? ((Asx(j, i) > 0.0) ? qSolX(j, i) / (dt * Asx(j, i)) : 0.0)
                 : uu(j, i);
        const real_t velW =
            QVOL ? ((Asx(j, i - 1) > 0.0) ? qSolX(j, i - 1) / (dt * Asx(j, i - 1)) : 0.0)
                 : uu(j, i - 1);
        const real_t velN =
            QVOL ? ((Asy(j, i) > 0.0) ? qSolY(j, i) / (dt * Asy(j, i)) : 0.0)
                 : vv(j, i);
        const real_t velS =
            QVOL ? ((Asy(j - 1, i) > 0.0) ? qSolY(j - 1, i) / (dt * Asy(j - 1, i)) : 0.0)
                 : vv(j - 1, i);

        // T10 (QVOL): infiltration concentration for a start-of-step-dry
        // cell (set in the pass-through block below; unused otherwise).
        real_t qvolInfilC = 0.0;
        bool qvolDryStart = false;
        (void)qvolInfilC;
        (void)qvolDryStart;

        // x advection (scalar.c:45-77).
        real_t sip;
        if (flowE > 0.0) {
          sip = s(j, i);
          if (superbee) {
            sip = tvdSuperbee(s(j, i + 1), s(j, i), s(j, i - 1), velE, dx, dt);
          }
        } else {
          sip = s(j, i + 1);
          if (superbee && !eastCell) {
            sip = tvdSuperbee(s(j, i), s(j, i + 1), sIp2, velE, dx, dt);
          }
        }
        real_t sim;
        if (flowW > 0.0) {
          sim = s(j, i - 1);
          if (superbee && !westCell) {
            sim = tvdSuperbee(s(j, i), s(j, i - 1), sIm2, velW, dx, dt);
          }
        } else {
          sim = s(j, i);
          if (superbee) {
            sim = tvdSuperbee(s(j, i - 1), s(j, i), s(j, i + 1), velW, dx, dt);
          }
        }
        // y advection (scalar.c:78-110).
        real_t sjp;
        if (flowN > 0.0) {
          sjp = s(j, i);
          if (superbee) {
            sjp = tvdSuperbee(s(j + 1, i), s(j, i), s(j - 1, i), velN, dy, dt);
          }
        } else {
          sjp = s(j + 1, i);
          if (superbee && !northCell) {
            sjp = tvdSuperbee(s(j, i), s(j + 1, i), sJp2, velN, dy, dt);
          }
        }
        real_t sjm;
        if (flowS > 0.0) {
          sjm = s(j - 1, i);
          if (superbee && !southCell) {
            sjm = tvdSuperbee(s(j, i), s(j - 1, i), sJm2, velS, dy, dt);
          }
        } else {
          sjm = s(j, i);
          if (superbee) {
            sjm = tvdSuperbee(s(j - 1, i), s(j, i), s(j + 1, i), velS, dy, dt);
          }
        }
        if constexpr (QVOL) {
          // C2 (pass-through front): the implicit solve can fill a cell
          // that started the step dry AND push water onward in the same
          // step; the outflow face's upwind value is then the stale dry
          // zero, so the water would move on without its salt (the cell
          // concentrates, downstream dilutes — the J1-breaking seed, found
          // on sx1 at the swash front). The water passing through a
          // start-of-step-dry cell IS this step's upstream water, so the
          // face carries the same-axis far-upwind value instead — one
          // level, both adjacent iterations compute the identical
          // substitute (the distance-2 reads use the superbee far staging,
          // with its edge guards).
          const bool dry0C = etan(j, i) - bottom(j, i) <= 0.0;
          if (flowE > 0.0) {
            if (dry0C) {
              sip = s(j, i - 1);
            }
          } else if (flowE < 0.0 && !eastCell &&
                     etan(j, i + 1) - bottom(j, i + 1) <= 0.0) {
            sip = sIp2;
          }
          if (flowW > 0.0) {
            if (!westCell && etan(j, i - 1) - bottom(j, i - 1) <= 0.0) {
              sim = sIm2;
            }
          } else if (flowW < 0.0 && dry0C) {
            sim = s(j, i + 1);
          }
          if (flowN > 0.0) {
            if (dry0C) {
              sjp = s(j - 1, i);
            }
          } else if (flowN < 0.0 && !northCell &&
                     etan(j + 1, i) - bottom(j + 1, i) <= 0.0) {
            sjp = sJp2;
          }
          if (flowS > 0.0) {
            if (!southCell && etan(j - 1, i) - bottom(j - 1, i) <= 0.0) {
              sjm = sJm2;
            }
          } else if (flowS < 0.0 && dry0C) {
            sjm = s(j + 1, i);
          }
          // C2 (bed face of a pass-through front): the same stale-zero
          // mechanism drains through the bed when the cell wets and
          // infiltrates in one step — the exchange term below then carries
          // the inflow mix instead of s(j, i).
          const real_t inVol = Kokkos::fmax(flowW, 0.0) + Kokkos::fmax(-flowE, 0.0) +
                               Kokkos::fmax(flowS, 0.0) + Kokkos::fmax(-flowN, 0.0);
          if (dry0C && inVol > 0.0) {
            qvolInfilC = (Kokkos::fmax(flowW, 0.0) * sim + Kokkos::fmax(-flowE, 0.0) * sip +
                          Kokkos::fmax(flowS, 0.0) * sjm + Kokkos::fmax(-flowN, 0.0) * sjp) /
                         inVol;
            qvolDryStart = true;
          }
          // C2: the salt moves with the water the solve moved.
          mass += -flowE * sip + flowW * sim - flowN * sjp + flowS * sjm;
        } else {
          mass += dt * (-Fu(j, i) * sip + Fu(j, i - 1) * sim - Fv(j, i) * sjp +
                        Fv(j - 1, i) * sjm);
        }

        // Advective mass through the global-edge faces (see TransportAudit):
        // the closed-edge fold keeps the water in, but the boundary-face
        // velocities are live and the scalar crosses them.
        if (westCell) {
          sumBoundary += QVOL ? flowW * sim : dt * Fu(j, i - 1) * sim;
        }
        if (eastCell) {
          sumBoundary -= QVOL ? flowE * sip : dt * Fu(j, i) * sip;
        }
        if (southCell) {
          sumBoundary += QVOL ? flowS * sjm : dt * Fv(j - 1, i) * sjm;
        }
        if (northCell) {
          sumBoundary -= QVOL ? flowN * sjp : dt * Fv(j, i) * sjp;
        }

        // Diffusion (scalar.c:114-119).
        mass += dt * ((difux * Asx(j, i) / dx) * (s(j, i + 1) - s(j, i)) -
                      (difux * Asx(j, i - 1) / dx) * (s(j, i) - s(j, i - 1)) +
                      (difuy * Asy(j, i) / dy) * (s(j + 1, i) - s(j, i)) -
                      (difuy * Asy(j - 1, i) / dy) * (s(j, i) - s(j - 1, i)));
        // Stage-3 T2 (W2.5, the 2026-10-26 Fix-B finding): diffusive mass
        // through the global-edge faces, booked with the advective edge
        // mass. The edge ghosts are zero-gradient copies taken before the
        // Dirichlet (tide) reset, so a reset stage cell diffuses against
        // its own pre-reset value; that exchange used to fall into the
        // surface ledger residual. Audit only: the update is unchanged.
        if (westCell) {
          sumBoundary -= dt * (difux * Asx(j, i - 1) / dx) * (s(j, i) - s(j, i - 1));
        }
        if (eastCell) {
          sumBoundary += dt * (difux * Asx(j, i) / dx) * (s(j, i + 1) - s(j, i));
        }
        if (southCell) {
          sumBoundary -= dt * (difuy * Asy(j - 1, i) / dy) * (s(j, i) - s(j - 1, i));
        }
        if (northCell) {
          sumBoundary += dt * (difuy * Asy(j, i) / dy) * (s(j + 1, i) - s(j, i));
        }

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
            if (held) {
              // W5.9 (D9): the exfiltrated water's salt moves AS THE WATER
              // LEAVES the subsurface — this step's fresh exfiltration
              // (exfNew) carries sKp whether the exchange applied or is
              // held; a held cell's share goes to the salt store (below,
              // via the dry branch: sm carries it and the cell is dry).
              // The APPLIED volume's held backlog carries no exchange salt
              // (its salt already sits in the store, which the wet branch
              // releases this very step).
              if (exfNew(j, i) > 0.0) {
                rate = exfNew(j, i) / dt * sKp(j, i);
              } else if (qss(j, i) < 0.0) {
                real_t cInf = s(j, i);
                if constexpr (QVOL) {
                  if (qvolDryStart) {
                    cInf = qvolInfilC;
                  }
                }
                rate = qss(j, i) * cInf;
              }
              if (dept(j, i) > 0.0) {
                rate += diffusive * (sKp(j, i) - s(j, i));
              }
              rate += vent(j, i) / (cellArea * dt) * sKp(j, i);
            } else if (consistentExchange) {
              // The interface diffusion acts on every wet step, not only on
              // steps whose seepage rate is nonzero (subcycled windows
              // without a substep have qss = 0), and the post-allocation's
              // vented water arrives at the top-cell concentration.
              if (qss(j, i) > 0.0 && dept(j, i) > 0.0) {
                rate = qss(j, i) * sKp(j, i);
              } else if (qss(j, i) < 0.0) {
                real_t cInf = s(j, i);
                if constexpr (QVOL) {
                  if (qvolDryStart) {
                    cInf = qvolInfilC;
                  }
                }
                rate = qss(j, i) * cInf;
              }
              if (dept(j, i) > 0.0) {
                rate += diffusive * (sKp(j, i) - s(j, i));
              }
              rate += vent(j, i) / (cellArea * dt) * sKp(j, i);
            } else if (qss(j, i) > 0.0) {
              if (dept(j, i) > 0.0) {
                rate = qss(j, i) * sKp(j, i) + diffusive * (sKp(j, i) - s(j, i));
              }
            } else if (qss(j, i) < 0.0) {
              rate = qss(j, i) * s(j, i) + diffusive * (sKp(j, i) - s(j, i));
            }
          }
          if constexpr (QVOL) {
            // T10 C3: the bounce-back deposits are in the flux volume (via
            // gain); their salt arrives with them at the top-cell
            // concentration — without this the bounced water dilutes the
            // film it lands on (the sub-35 seed measured on sx1's beach).
            rate += cplBounce(j, i) / dt * sKp(j, i);
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
        // T10: a face the solve moved water across admits its neighbour
        // into the window even when the post-step geometry closed it —
        // the arrived salt must not read as a minted extremum (J4).
        if ((Asx(j, i) > 0.0 && etan(j, i + 1) > bottom(j, i + 1)) ||
            (QVOL && flowE != 0.0)) {
          extend(s(j, i + 1));
        }
        if ((Asx(j, i - 1) > 0.0 && etan(j, i - 1) > bottom(j, i - 1)) ||
            (QVOL && flowW != 0.0)) {
          extend(s(j, i - 1));
        }
        if ((Asy(j, i) > 0.0 && etan(j + 1, i) > bottom(j + 1, i)) ||
            (QVOL && flowN != 0.0)) {
          extend(s(j + 1, i));
        }
        if ((Asy(j - 1, i) > 0.0 && etan(j - 1, i) > bottom(j - 1, i)) ||
            (QVOL && flowS != 0.0)) {
          extend(s(j - 1, i));
        }
        if (coupled && topKsz > 0.0) {
          extend(sKp(j, i));
        }
        sMin(j, i) = lo;
        sMax(j, i) = hi;

        // Flux volume (volume_by_flux, shallowwater.c:996-1032) from the
        // previous step's flow rates; the coupled term is the applied
        // exchange volume without the legacy porosity factor (A9). Under
        // transport.surface_update: consistent the volume uses the current
        // flow rates instead — the rates that moved the mass above — so the
        // division is the conservative flux-form update (a uniform scalar
        // stays uniform, and a film whose outflow the wetting limiter just
        // closed is not concentrated by last step's outflow).
        const real_t qxW = consistent ? Fu(j, i - 1) : fuOld(j, i - 1);
        const real_t qxE = consistent ? Fu(j, i) : fuOld(j, i);
        const real_t qyS = consistent ? Fv(j - 1, i) : fvOld(j - 1, i);
        const real_t qyN = consistent ? Fv(j, i) : fvOld(j, i);
        real_t v = vsn(j, i) + dt * (qxW - qxE + qyS - qyN);
        if constexpr (QVOL) {
          // C2: the post-solve, pre-manipulation water volume.
          v = vsn(j, i) + (flowW - flowE + flowS - flowN);
        }
        if constexpr (QVOL) {
          if (coupled) {
            // C3: the coupled deposits as they were actually applied to eta
            // this step (seepage post-clamp/post-hold, bounce-back, vents)
            // — the gain field keeps the manipulation's own books, so no
            // re-derivation from qss (whose clamp/hold approximation is
            // what the anchor used to absorb).
            v += cplGain(j, i) * cellArea;
          }
          if (stageMask(j, i) != 0.0) {
            // C3 (stage rows): a prescribed-stage cell's volume comes from
            // the identity row and the closure's reservoir routing, not
            // from its face volumes — its flux volume IS the actual water
            // it holds. The Dirichlet reset books the salt difference into
            // surf_source, so the anchor term vanishes here instead of
            // accumulating the reservoir's unbooked volume.
            v = dept(j, i) * cellArea;
          }
        } else if (coupled) {
          const real_t area = (dept(j, i) > 0.0) ? cellArea : 0.0;
          const real_t dv = qss(j, i) * dt * area;
          if (qss(j, i) < 0.0 && v > -dv) {
            v += dv;
          } else if (qss(j, i) > 0.0) {
            if (dept(j, i) > 0.0) {
              v += dv;
            }
          }
          if (consistentExchange) {
            v += vent(j, i);  // the vented water is already in the depth
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
  real_t atmosMass = 0.0;
  real_t storeInMass = 0.0;
  real_t storeOutMass = 0.0;
  // Stage-3 T9 (W9.3, plan §12A): the paired store's per-step surface
  // books. R4 (salt kept without its water) is salinity-only; the
  // temperature instance parks heat with its water and nothing else (R9).
  Field2<real_t> parkedV = store ? surf_.parkedVol : Field2<real_t>();
  Field2<real_t> releasedV = store ? surf_.releasedVol : Field2<real_t>();
  Field2<real_t> evapV = store ? surf_.evapStep : Field2<real_t>();
  const bool crust = store && !spec_.isTemperature;
  const bool heldKeep = crust && held;
  Field2<real_t> seepAcc = heldKeep ? cpl_.seepAccum : Field2<real_t>();
  Kokkos::parallel_reduce(
      "transport_surface_update",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
          {1, 1}, {nyl + 1, nxl + 1}),
      KOKKOS_LAMBDA(const int j, const int i, real_t& sumAdjust, real_t& sumAnchor,
                    real_t& sumAtmos, real_t& sumStoreIn, real_t& sumStoreOut) {
        // nvcc: an extended __host__ __device__ lambda may not FIRST capture
        // a variable inside an if-constexpr branch (CUDA C++ guide, "Extended
        // Lambda Restrictions"), so the QVOL branch's captures are named
        // here first. Discarded expressions: no code, no arithmetic.
        (void)stageMask;
        (void)s;
        const bool wet = (vflux(j, i) > 0.0 && dept(j, i) > 0.0);
        // T9 (W9.3, R8): under the paired store "dry" means no surface
        // water after the surface step (dept <= 0 — the W9.2 sweep made
        // that exact); a cell holding water but no flux volume (wetted from
        // dry by rain or a store release) keeps its mass in that water.
        const bool fromDry = store && !wet && dept(j, i) > 0.0;
        real_t value = wet ? sm(j, i) / vflux(j, i) : 0.0;
        const real_t raw = value;
        if (wet) {
          const bool moving = (uu(j, i) != 0.0 || uu(j, i - 1) != 0.0 || vv(j, i) != 0.0 ||
                               vv(j - 1, i) != 0.0);
          // Legacy clamps moving cells only; the consistent update clamps
          // every wet cell (a still film can still receive mass/volume
          // mismatches from the exchange or a closing face).
          bool applyLimiter = moving || consistent;
          if (!applyLimiter && coupled) {
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
          if constexpr (QVOL) {
            if (stageMask(j, i) != 0.0) {
              // C3 (stage rows): the stencil limiter would mint salt onto
              // the reservoir's fresh volume (raw < window on a rising
              // tide) and book it as adjust; the Dirichlet reset below is
              // the stage cell's one honest book (surf_source), so the
              // prescribed cell skips the limiter.
              applyLimiter = false;
            }
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
          } else if (hasMin && value < limLo) {
            value = limLo;
          }
          sumAdjust += (value - raw) * vflux(j, i);
        } else if (store) {
          if (!fromDry) {
            // T9 (W9.3, R3): the water parked this step (P1-P3 and the
            // sweep) takes its salt along at the cell's limited post-source
            // concentration, never above the limiter-window maximum (the
            // flux-volume basis when the cell had one, else the parked
            // volume itself). Salt is never created: at most the dry-branch
            // mass parks, so a parcel the transport books salt-free (an
            // isolated cell's solve-delivered water) parks fresh.
            const real_t dv = parkedV(j, i);
            real_t c = 0.0;
            if (vflux(j, i) > 0.0) {
              c = sm(j, i) / vflux(j, i);
            } else if (dv > 0.0) {
              c = sm(j, i) / dv;
            }
            if (c > sMax(j, i) && sMax(j, i) != limHi && sMax(j, i) != limLo) {
              c = sMax(j, i);
            }
            if (hasMax && c > limHi) {
              c = limHi;
            }
            real_t parkM = Kokkos::fmin(Kokkos::fmax(sm(j, i), 0.0),
                                        Kokkos::fmax(c, 0.0) * dv);
            if constexpr (QVOL) {
              // C4 (redesigned under the C-SX contingency; S3-A15): an
              // overdrawn cell — the solve drained it below its bed, the
              // below-bed lift restored dept (booked in `clamped`), and
              // the face volumes net negative here — parked the lifted
              // water. A donor-side outflow cap was measured first and
              // REJECTED: it starves the real receiving cells of their
              // water's salt (8.2 psu J1 breach on sx1). Instead the
              // lifted water parks PAIRED, at the cell's start-of-step
              // concentration bounded by the limiter window; the ledger
              // difference (the salt the faces over-exported plus this
              // pairing) stays on the negative-residual path into
              // surf_adjust — one honest book for the solve's minted
              // volume, and the store's water–salt pairing (J1 at the
              // release) survives.
              if (vflux(j, i) < 0.0 && dv > 0.0) {
                real_t cPair = s(j, i);
                if (cPair > sMax(j, i) && sMax(j, i) != limHi && sMax(j, i) != limLo) {
                  cPair = sMax(j, i);
                }
                if (hasMax && cPair > limHi) {
                  cPair = limHi;
                }
                parkM = Kokkos::fmax(cPair, 0.0) * dv;
              }
            }
            real_t keep = parkM;
            real_t rest = sm(j, i) - parkM;
            // R4: salt left without its water stays only as (a) the crust
            // of this step's evaporation — at most the salt the evaporated
            // water carried, c V_evap at the same limited concentration —
            // or (b) the W5.9 held-seepage salt, at most this step's held
            // exchange credit while the coupler holds its water. Any other
            // unpaired or negative remainder follows the discard path into
            // surf_adjust (an advective dry-front artifact never becomes a
            // store, S3-A13).
            if (rest > 0.0 && crust && evapV(j, i) > 0.0) {
              const real_t k = Kokkos::fmin(rest, Kokkos::fmax(c, 0.0) * evapV(j, i));
              keep += k;
              rest -= k;
            }
            if (rest > 0.0 && heldKeep && seepAcc(j, i) > 0.0) {
              const real_t k = Kokkos::fmin(rest, Kokkos::fmax(ssee(j, i) * cellArea * dt, 0.0));
              keep += k;
              rest -= k;
            }
            sumAdjust -= rest;
            sStore(j, i) += keep;
            sumStoreIn += keep;
          }
        } else {
          sumAdjust -= sm(j, i);
        }
        // Rain/evaporation dilution (scalar.c:243-254; salinity only —
        // temperature is invariant under volume exchange at the cell's own
        // temperature, V2-A17).
        real_t vre = 0.0;
        if (dilution && vflux(j, i) > 0.0) {
          const real_t area = (dept(j, i) > 0.0) ? cellArea : 0.0;
          if (etSurf) {
            // W5.7(c)/D8: the actual per-cell volumes — rain AND
            // evaporation together, masks already honoured by the fields
            // (legacy: either/or on a global threshold test).
            vre = (etRain(j, i) - etEvap(j, i)) * area;
          } else {
            const real_t cellRain = rainRate * rainMask(j, i);
            if (dept(j, i) > minDepth && cellRain > 0.0) {
              vre += cellRain * area * dt;
            } else {
              vre -= evapRate * area * dt;
            }
          }
          if (vflux(j, i) > 0.0) {
            value = value * vflux(j, i) / (vflux(j, i) + vre);
          }
        }
        // Surface heat exchange (v2 Q5, plan §4.1): an explicit source on
        // the final volume basis, after the limiter — an atmospheric flux
        // must be able to cross the local advective extrema (relaxation
        // toward a T_e below the wet-stencil minimum would otherwise
        // stall). Booked into surf_atmos; the anchor below then closes the
        // ledger on the final value.
        // T9 (R8): a wetted-from-dry cell's basis is the water it holds
        // apart from what the store returned this step (rain or nothing):
        // its mass sm lives there; the release below mixes in the store.
        real_t dVs = 0.0;
        if (store && dept(j, i) > 0.0) {
          dVs = Kokkos::fmax(releasedV(j, i) - parkedV(j, i), 0.0);
        }
        real_t vBasis = vflux(j, i) + vre;
        if (fromDry) {
          vBasis = dept(j, i) * cellArea - dVs;
          if (vBasis > 0.0) {
            value = sm(j, i) / vBasis;
          } else {
            vBasis = 0.0;
            sumAdjust -= sm(j, i);
          }
        }
        if ((atmosEquilibrium || atmosBulk) && dept(j, i) > 0.0 && vBasis > 0.0) {
          const real_t qnet =
              atmosEquilibrium
                  ? -equilibriumK * (value - equilibriumT)
                  : atm::netHeatFlux(value, metAirT, metPressure, metWind, metQAir,
                                     metShortwave, metLongwaveIn);
          const real_t dTv = dt * cellArea * qnet / rhoCw;  // [K m^3]
          value += dTv / vBasis;
          sumAtmos += dTv;
        }
        // T9 (W9.3, R5): the release is mixing. The cell holds water, so
        // the surface returned the whole water store to eta this step (net
        // dVs beyond this step's own parks); the salt store joins on the
        // combined basis, c = (M + M_s) / (V + V_s) — after the limiter,
        // bounds and dilution, so the released salt is never clipped, and
        // the anchor stays the basis-lag correction (never a salt source).
        // A crust (O4, V_s = 0) mixes the same way (dVs = 0).
        bool rewet = false;
        real_t basisOut = vBasis;
        if (store && dept(j, i) > 0.0 && (sStore(j, i) > 0.0 || dVs > 0.0) &&
            vBasis + dVs > 0.0) {
          rewet = true;
          basisOut = vBasis + dVs;
          value = (value * vBasis + sStore(j, i)) / basisOut;
          sumStoreOut += sStore(j, i);
          sStore(j, i) = 0.0;
        }
        s(j, i) = value;
        // Ledger re-anchor: the mass now lives on the actual volume, not
        // the (lagged) flux volume the update divided by.
        if (wet || rewet || fromDry) {
          sumAnchor += value * (dept(j, i) * cellArea - basisOut);
        }
      },
      adjustMass, anchorMass, atmosMass, storeInMass, storeOutMass);
  audit_.surfAdjust = adjustMass;
  audit_.surfAnchor = anchorMass;
  audit_.surfAtmos = atmosMass;
  audit_.storeIn = storeInMass;
  audit_.storeOut = storeOutMass;

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

  halo_.exchange({spec_.prefix + "_surf"});
}

}  // namespace frehg::transport
