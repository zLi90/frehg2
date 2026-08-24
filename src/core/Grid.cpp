/// \file Grid.cpp
/// \brief Implementation of the structured grid and compressed global ids.

#include "core/Grid.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace frehg {

namespace {

/// MPI datatype matching PetscInt (32- or 64-bit builds).
MPI_Datatype petscIntMpiType() {
  return sizeof(PetscInt) == 8 ? MPI_LONG_LONG : MPI_INT;
}

}  // namespace

int Grid::blockSize(int N, int P, int p) { return N / P + (p < N % P ? 1 : 0); }

int Grid::blockStart(int N, int P, int p) { return p * (N / P) + std::min(p, N % P); }

Grid::Grid(MPI_Comm comm, const DomainConfig& domain)
    : comm_(comm),
      nx_(domain.nx),
      ny_(domain.ny),
      nz_(domain.nz),
      dx_(domain.dx),
      dy_(domain.dy) {
  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &size_);

  if (nx_ < 1 || ny_ < 1 || nz_ < 1) {
    log::fatal(log::msg() << "Grid: extents must be positive (nx=" << nx_ << ", ny=" << ny_
                          << ", nz=" << nz_ << ")");
  }

  // Vertical spacing: geometric stretch downward from the top layer.
  dz_.resize(static_cast<std::size_t>(nz_));
  zCenter_.resize(static_cast<std::size_t>(nz_));
  real_t depth = 0;
  for (int k = 0; k < nz_; ++k) {
    dz_[static_cast<std::size_t>(k)] = domain.dz * std::pow(domain.dzStretch, k);
    zCenter_[static_cast<std::size_t>(k)] = -(depth + 0.5 * dz_[static_cast<std::size_t>(k)]);
    depth += dz_[static_cast<std::size_t>(k)];
  }
  totalDepth_ = depth;

  chooseDecomposition(domain.decomposition);

  pi_ = rank_ % px_;
  pj_ = rank_ / px_;
  nxLocal_ = blockSize(nx_, px_, pi_);
  nyLocal_ = blockSize(ny_, py_, pj_);
  i0_ = blockStart(nx_, px_, pi_);
  j0_ = blockStart(ny_, py_, pj_);
  if (nxLocal_ < 1 || nyLocal_ < 1) {
    log::fatal(log::msg() << "Grid: decomposition " << px_ << " x " << py_
                          << " leaves rank " << rank_ << " with an empty block ("
                          << nxLocal_ << " x " << nyLocal_ << " cells); use fewer ranks");
  }

  rankWest_ = (pi_ > 0) ? rank_ - 1 : MPI_PROC_NULL;
  rankEast_ = (pi_ + 1 < px_) ? rank_ + 1 : MPI_PROC_NULL;
  rankSouth_ = (pj_ > 0) ? rank_ - px_ : MPI_PROC_NULL;
  rankNorth_ = (pj_ + 1 < py_) ? rank_ + px_ : MPI_PROC_NULL;

  // Default mask: everything active.
  HostField2<int> allActive("ktop_default", static_cast<std::size_t>(nyLocal_),
                            static_cast<std::size_t>(nxLocal_));
  Kokkos::deep_copy(allActive, 0);
  buildGlobalIds(allActive);
}

void Grid::chooseDecomposition(const DecompositionConfig& request) {
  int reqX = request.mpiNx;
  int reqY = request.mpiNy;
  if (reqX > 0 && reqY > 0) {
    if (reqX * reqY != size_) {
      log::fatal(log::msg() << "Grid: decomposition " << reqX << " x " << reqY << " needs "
                            << reqX * reqY << " ranks but the run has " << size_);
    }
    px_ = reqX;
    py_ = reqY;
    return;
  }
  if (reqX > 0) {
    if (size_ % reqX != 0) {
      log::fatal(log::msg() << "Grid: mpi_nx = " << reqX << " does not divide " << size_
                            << " ranks");
    }
    px_ = reqX;
    py_ = size_ / reqX;
    return;
  }
  if (reqY > 0) {
    if (size_ % reqY != 0) {
      log::fatal(log::msg() << "Grid: mpi_ny = " << reqY << " does not divide " << size_
                            << " ranks");
    }
    py_ = reqY;
    px_ = size_ / reqY;
    return;
  }
  std::array<int, 2> dims = {0, 0};
  MPI_Dims_create(size_, 2, dims.data());
  // Give the larger rank count to the larger extent for squarer blocks.
  const int large = std::max(dims[0], dims[1]);
  const int small = std::min(dims[0], dims[1]);
  if (nx_ >= ny_) {
    px_ = large;
    py_ = small;
  } else {
    px_ = small;
    py_ = large;
  }
  // MPI_Dims_create ignores the grid extents; when its layout would leave a
  // dimension over-decomposed (empty blocks — e.g. 2 x 2 ranks on a 1 x N
  // grid), fall back to the most balanced factorization that fits.
  if (px_ > nx_ || py_ > ny_) {
    int bestX = 0;
    long bestScore = -1;
    for (int candidate = 1; candidate <= size_; ++candidate) {
      if (size_ % candidate != 0 || candidate > nx_ || size_ / candidate > ny_) {
        continue;
      }
      // Prefer squarer blocks: maximize the smaller block extent.
      const long score = static_cast<long>(std::min(nx_ / candidate, ny_ / (size_ / candidate)));
      if (score > bestScore) {
        bestScore = score;
        bestX = candidate;
      }
    }
    if (bestX > 0) {
      px_ = bestX;
      py_ = size_ / bestX;
    }
    // No valid factorization: keep the invalid layout so the empty-block
    // check reports the failure with its actionable message.
  }
}

void Grid::buildGlobalIds(const HostField2<int>& ktop) {
  if (ktop.extent(0) != static_cast<std::size_t>(nyLocal_) ||
      ktop.extent(1) != static_cast<std::size_t>(nxLocal_)) {
    log::fatal(log::msg() << "Grid::buildGlobalIds: ktop is " << ktop.extent(0) << " x "
                          << ktop.extent(1) << ", expected " << nyLocal_ << " x " << nxLocal_);
  }

  ktopHost_ = HostField2<int>("ktop", static_cast<std::size_t>(nyLocal_),
                              static_cast<std::size_t>(nxLocal_));
  Kokkos::deep_copy(ktopHost_, ktop);

  long long count2 = 0;
  long long count3 = 0;
  for (int j = 0; j < nyLocal_; ++j) {
    for (int i = 0; i < nxLocal_; ++i) {
      const int kt = ktopHost_(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
      if (kt < 0 || kt > nz_) {
        log::fatal(log::msg() << "Grid::buildGlobalIds: ktop(" << j << "," << i << ") = " << kt
                              << " outside [0, nz=" << nz_ << "]");
      }
      if (kt < nz_) {
        count2 += 1;
        count3 += nz_ - kt;
      }
    }
  }

  long long scan2 = 0;
  long long scan3 = 0;
  MPI_Exscan(&count2, &scan2, 1, MPI_LONG_LONG, MPI_SUM, comm_);
  MPI_Exscan(&count3, &scan3, 1, MPI_LONG_LONG, MPI_SUM, comm_);
  if (rank_ == 0) {
    scan2 = 0;
    scan3 = 0;
  }
  long long total2 = 0;
  long long total3 = 0;
  MPI_Allreduce(&count2, &total2, 1, MPI_LONG_LONG, MPI_SUM, comm_);
  MPI_Allreduce(&count3, &total3, 1, MPI_LONG_LONG, MPI_SUM, comm_);

  constexpr long long kPetscIntMax = static_cast<long long>(std::numeric_limits<PetscInt>::max());
  if (total3 > kPetscIntMax || total2 > kPetscIntMax) {
    log::fatal(log::msg() << "Grid: " << total3
                          << " active cells overflow PetscInt; rebuild PETSc with 64-bit indices");
  }

  active2Local_ = static_cast<PetscInt>(count2);
  active3Local_ = static_cast<PetscInt>(count3);
  active2Global_ = static_cast<PetscInt>(total2);
  active3Global_ = static_cast<PetscInt>(total3);
  offset2_ = static_cast<PetscInt>(scan2);
  offset3_ = static_cast<PetscInt>(scan3);

  // Owned interior ids in (j, i, k) order, k innermost — matching both the
  // output flattening (j*NX + i)*NZ + k and PETSc's contiguous row ownership.
  const std::size_t nyH = static_cast<std::size_t>(nyLocal_) + 2;
  const std::size_t nxH = static_cast<std::size_t>(nxLocal_) + 2;
  gid2Host_ = HostField2<PetscInt>("gid2", nyH, nxH);
  gid3Host_ = HostField3<PetscInt>("gid3", nyH, nxH, static_cast<std::size_t>(nz_));
  Kokkos::deep_copy(gid2Host_, PetscInt{-1});
  Kokkos::deep_copy(gid3Host_, PetscInt{-1});

  PetscInt next2 = offset2_;
  PetscInt next3 = offset3_;
  for (int j = 0; j < nyLocal_; ++j) {
    for (int i = 0; i < nxLocal_; ++i) {
      const int kt = ktopHost_(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
      if (kt >= nz_) {
        continue;
      }
      const std::size_t jh = static_cast<std::size_t>(j) + 1;
      const std::size_t ih = static_cast<std::size_t>(i) + 1;
      gid2Host_(jh, ih) = next2++;
      for (int k = kt; k < nz_; ++k) {
        gid3Host_(jh, ih, static_cast<std::size_t>(k)) = next3++;
      }
    }
  }

  exchangeGidHalos();

  gid2_ = Field2<PetscInt>("gid2_device", nyH, nxH);
  gid3_ = Field3<PetscInt>("gid3_device", nyH, nxH, static_cast<std::size_t>(nz_));
  Kokkos::deep_copy(gid2_, gid2Host_);
  Kokkos::deep_copy(gid3_, gid3Host_);
}

void Grid::exchangeGidHalos() {
  const MPI_Datatype dtype = petscIntMpiType();
  const std::size_t nzS = static_cast<std::size_t>(nz_);
  const std::size_t nyl = static_cast<std::size_t>(nyLocal_);
  const std::size_t nxl = static_cast<std::size_t>(nxLocal_);

  // West/East halos: exchange full interior columns (2D + 3D coalesced).
  {
    const std::size_t planeCount = nyl * (1 + nzS);
    std::vector<PetscInt> sendW(planeCount), sendE(planeCount);
    std::vector<PetscInt> recvW(planeCount), recvE(planeCount);
    auto packColumn = [&](std::size_t ih, std::vector<PetscInt>& buf) {
      std::size_t n = 0;
      for (std::size_t j = 1; j <= nyl; ++j) {
        buf[n++] = gid2Host_(j, ih);
        for (std::size_t k = 0; k < nzS; ++k) {
          buf[n++] = gid3Host_(j, ih, k);
        }
      }
    };
    auto unpackColumn = [&](std::size_t ih, const std::vector<PetscInt>& buf) {
      std::size_t n = 0;
      for (std::size_t j = 1; j <= nyl; ++j) {
        gid2Host_(j, ih) = buf[n++];
        for (std::size_t k = 0; k < nzS; ++k) {
          gid3Host_(j, ih, k) = buf[n++];
        }
      }
    };
    packColumn(1, sendW);
    packColumn(nxl, sendE);
    const int cnt = static_cast<int>(planeCount);
    MPI_Sendrecv(sendW.data(), cnt, dtype, rankWest_, 100, recvE.data(), cnt, dtype, rankEast_,
                 100, comm_, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendE.data(), cnt, dtype, rankEast_, 101, recvW.data(), cnt, dtype, rankWest_,
                 101, comm_, MPI_STATUS_IGNORE);
    if (rankEast_ != MPI_PROC_NULL) {
      unpackColumn(nxl + 1, recvE);
    }
    if (rankWest_ != MPI_PROC_NULL) {
      unpackColumn(0, recvW);
    }
  }

  // South/North halos: exchange full interior rows.
  {
    const std::size_t planeCount = nxl * (1 + nzS);
    std::vector<PetscInt> sendS(planeCount), sendN(planeCount);
    std::vector<PetscInt> recvS(planeCount), recvN(planeCount);
    auto packRow = [&](std::size_t jh, std::vector<PetscInt>& buf) {
      std::size_t n = 0;
      for (std::size_t i = 1; i <= nxl; ++i) {
        buf[n++] = gid2Host_(jh, i);
        for (std::size_t k = 0; k < nzS; ++k) {
          buf[n++] = gid3Host_(jh, i, k);
        }
      }
    };
    auto unpackRow = [&](std::size_t jh, const std::vector<PetscInt>& buf) {
      std::size_t n = 0;
      for (std::size_t i = 1; i <= nxl; ++i) {
        gid2Host_(jh, i) = buf[n++];
        for (std::size_t k = 0; k < nzS; ++k) {
          gid3Host_(jh, i, k) = buf[n++];
        }
      }
    };
    packRow(1, sendS);
    packRow(nyl, sendN);
    const int cnt = static_cast<int>(planeCount);
    MPI_Sendrecv(sendS.data(), cnt, dtype, rankSouth_, 102, recvN.data(), cnt, dtype, rankNorth_,
                 102, comm_, MPI_STATUS_IGNORE);
    MPI_Sendrecv(sendN.data(), cnt, dtype, rankNorth_, 103, recvS.data(), cnt, dtype, rankSouth_,
                 103, comm_, MPI_STATUS_IGNORE);
    if (rankNorth_ != MPI_PROC_NULL) {
      unpackRow(nyl + 1, recvN);
    }
    if (rankSouth_ != MPI_PROC_NULL) {
      unpackRow(0, recvS);
    }
  }
}

}  // namespace frehg
