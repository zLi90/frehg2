/// \file Dispersion.cpp
/// \brief The anisotropic dispersion tensor (dispersion_tensor,
///        scalar.c:958-1003).
///
/// D = θs D_molecular I + tensor(α_L, α_T, q) with q the *volumetric* face
/// fluxes [m^3/s] exactly as legacy evaluates it — the legacy tensor never
/// divides the face-area factor out of qx/qy/qz, so the dispersivities are
/// calibrated against face-scaled fluxes, and the b6 goldens embed that
/// scale (docs/theory/transport.md records the quirk). Each cell uses its
/// own plus-side face fluxes (qx/qy at the x+/y+ faces, qz at the lower
/// face). Off-diagonal entries are clamped at zero from below exactly as
/// legacy does (:996-1001). The symmetric pairs (Dyx = Dxy, Dzx = Dxz,
/// Dzy = Dyz — legacy stores both members and always assigns them the same
/// value) are stored once.
///
/// Interior cells are computed here and the six fields halo-exchanged: the
/// dispersive face fluxes evaluate their coefficient at each face's minus
/// cell, which is a neighbor-rank cell at interface faces (legacy read the
/// exchanged ghost D the same way).

#include "transport/ScalarSolver.hpp"

namespace frehg::transport {

void ScalarSolver::updateDispersionTensor() {
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  const int nz = grid_.nz();
  const real_t molecular = dispMol_;
  const real_t lon = dispLon_;
  const real_t lat = dispLat_;

  Field3<real_t> dxx = dxx_, dyy = dyy_, dzz = dzz_;
  Field3<real_t> dxy = dxy_, dxz = dxz_, dyz = dyz_;
  Field3<real_t> qx = subs_.qx, qy = subs_.qy, qzF = subs_.qzF;
  Field3<real_t> wcs = subs_.wcs;

  Kokkos::parallel_for(
      "transport_dispersion_tensor",
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
          {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
      KOKKOS_LAMBDA(const int j, const int i, const int k) {
        const real_t ux = qx(j, i, k);
        const real_t uy = qy(j, i, k);
        // Legacy qz[ii] is the flux at the face below the cell — plane k+1.
        const real_t uz = qzF(j, i, k + 1);
        real_t xx = molecular * wcs(j, i, k);
        real_t yy = xx;
        real_t zz = xx;
        real_t xy = 0.0;
        real_t xz = 0.0;
        real_t yz = 0.0;
        const real_t qAbs = Kokkos::sqrt(ux * ux + uy * uy + uz * uz);
        if (qAbs > 0.0) {
          xx += (lon * ux * ux + lat * uy * uy + lat * uz * uz) / qAbs;
          yy += (lon * uy * uy + lat * ux * ux + lat * uz * uz) / qAbs;
          zz += (lon * uz * uz + lat * uy * uy + lat * ux * ux) / qAbs;
          xy = (lon - lat) * ux * uy / qAbs;
          xz = (lon - lat) * ux * uz / qAbs;
          yz = (lon - lat) * uy * uz / qAbs;
        }
        if (xy < 0.0) {
          xy = 0.0;
        }
        if (xz < 0.0) {
          xz = 0.0;
        }
        if (yz < 0.0) {
          yz = 0.0;
        }
        dxx(j, i, k) = xx;
        dyy(j, i, k) = yy;
        dzz(j, i, k) = zz;
        dxy(j, i, k) = xy;
        dxz(j, i, k) = xz;
        dyz(j, i, k) = yz;
      });

  halo_.exchange({"s_dxx", "s_dyy", "s_dzz", "s_dxy", "s_dxz", "s_dyz"});
}

}  // namespace frehg::transport
