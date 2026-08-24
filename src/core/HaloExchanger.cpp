/// \file HaloExchanger.cpp
/// \brief Implementation of the coalescing halo exchange.

#include "core/HaloExchanger.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <numeric>

namespace frehg {

namespace {

/// Interior index of the plane sent in a direction, and the halo index of
/// the plane received from it, along the direction's axis.
struct PlaneIndices {
  std::size_t send;
  std::size_t recv;
};

PlaneIndices planeIndices(int direction, std::size_t nxl, std::size_t nyl) {
  switch (direction) {
    case 0:  return {1, 0};             // west: send i=1, recv into i=0
    case 1:  return {nxl, nxl + 1};     // east
    case 2:  return {1, 0};             // south: send j=1, recv into j=0
    default: return {nyl, nyl + 1};     // north
  }
}

}  // namespace

HaloExchanger::HaloExchanger(const Grid& grid, bool gpuAwareMpi)
    : grid_(grid), gpuAware_(gpuAwareMpi) {
  neighbor_ = {grid_.rankWest(), grid_.rankEast(), grid_.rankSouth(), grid_.rankNorth()};
}

std::size_t HaloExchanger::planeCount(const Entry& entry, int direction, bool wideRows) const {
  const bool alongX = (direction == kWest || direction == kEast);
  std::size_t rowCells = alongX ? static_cast<std::size_t>(grid_.nyLocal())
                                : static_cast<std::size_t>(grid_.nxLocal());
  if (wideRows && !alongX) {
    rowCells += 2;  // include the i-halo columns so corner values travel
  }
  return entry.is3d ? rowCells * static_cast<std::size_t>(grid_.nz()) : rowCells;
}

void HaloExchanger::add(const std::string& name, const Field2<real_t>& field) {
  for (const Entry& e : entries_) {
    if (e.name == name) {
      log::fatal(log::msg() << "HaloExchanger: field '" << name << "' registered twice");
    }
  }
  if (field.extent(0) != static_cast<std::size_t>(grid_.nyLocal()) + 2 ||
      field.extent(1) != static_cast<std::size_t>(grid_.nxLocal()) + 2) {
    log::fatal(log::msg() << "HaloExchanger: field '" << name << "' is " << field.extent(0)
                          << " x " << field.extent(1) << ", expected "
                          << grid_.nyLocal() + 2 << " x " << grid_.nxLocal() + 2);
  }
  Entry entry;
  entry.name = name;
  entry.is3d = false;
  entry.f2 = field;
  entries_.push_back(entry);
  ensureCapacity();
}

void HaloExchanger::add(const std::string& name, const Field3<real_t>& field) {
  for (const Entry& e : entries_) {
    if (e.name == name) {
      log::fatal(log::msg() << "HaloExchanger: field '" << name << "' registered twice");
    }
  }
  if (field.extent(0) != static_cast<std::size_t>(grid_.nyLocal()) + 2 ||
      field.extent(1) != static_cast<std::size_t>(grid_.nxLocal()) + 2 ||
      field.extent(2) != static_cast<std::size_t>(grid_.nz())) {
    log::fatal(log::msg() << "HaloExchanger: field '" << name << "' is " << field.extent(0)
                          << " x " << field.extent(1) << " x " << field.extent(2) << ", expected "
                          << grid_.nyLocal() + 2 << " x " << grid_.nxLocal() + 2 << " x "
                          << grid_.nz());
  }
  Entry entry;
  entry.name = name;
  entry.is3d = true;
  entry.f3 = field;
  entries_.push_back(entry);
  ensureCapacity();
}

void HaloExchanger::ensureCapacity() {
  for (int dir = 0; dir < 4; ++dir) {
    std::size_t needed = 0;
    for (const Entry& e : entries_) {
      needed += planeCount(e, dir, true);  // size for the widest (corner) mode
    }
    const std::size_t d = static_cast<std::size_t>(dir);
    if (sendBuf_[d].extent(0) < needed) {
      sendBuf_[d] = Buffer("halo_send", needed);
      recvBuf_[d] = Buffer("halo_recv", needed);
      sendHost_[d] = HostBuffer("halo_send_host", needed);
      recvHost_[d] = HostBuffer("halo_recv_host", needed);
    }
  }
}

void HaloExchanger::packEntry(const Entry& entry, int direction, std::size_t offset,
                              bool wideRows) {
  const std::size_t nxl = static_cast<std::size_t>(grid_.nxLocal());
  const std::size_t nyl = static_cast<std::size_t>(grid_.nyLocal());
  const std::size_t nzS = static_cast<std::size_t>(grid_.nz());
  const bool alongX = (direction == kWest || direction == kEast);
  const bool wide = wideRows && !alongX;
  const std::size_t rowCells = alongX ? nyl : (wide ? nxl + 2 : nxl);
  const std::size_t rowStart = wide ? 0 : 1;
  const std::size_t plane = planeIndices(direction, nxl, nyl).send;
  Buffer buf = sendBuf_[static_cast<std::size_t>(direction)];

  if (entry.is3d) {
    Field3<real_t> f = entry.f3;
    Kokkos::parallel_for(
        "halo_pack3", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0, 0}, {rowCells, nzS}),
        KOKKOS_LAMBDA(const std::size_t r, const std::size_t k) {
          const std::size_t jj = alongX ? r + rowStart : plane;
          const std::size_t ii = alongX ? plane : r + rowStart;
          buf(offset + r * nzS + k) = f(jj, ii, k);
        });
  } else {
    Field2<real_t> f = entry.f2;
    Kokkos::parallel_for(
        "halo_pack2", Kokkos::RangePolicy<ExecSpace>(0, rowCells),
        KOKKOS_LAMBDA(const std::size_t r) {
          const std::size_t jj = alongX ? r + rowStart : plane;
          const std::size_t ii = alongX ? plane : r + rowStart;
          buf(offset + r) = f(jj, ii);
        });
  }
}

void HaloExchanger::unpackEntry(const Entry& entry, int direction, std::size_t offset,
                                bool wideRows) {
  const std::size_t nxl = static_cast<std::size_t>(grid_.nxLocal());
  const std::size_t nyl = static_cast<std::size_t>(grid_.nyLocal());
  const std::size_t nzS = static_cast<std::size_t>(grid_.nz());
  const bool alongX = (direction == kWest || direction == kEast);
  const bool wide = wideRows && !alongX;
  const std::size_t rowCells = alongX ? nyl : (wide ? nxl + 2 : nxl);
  const std::size_t rowStart = wide ? 0 : 1;
  const std::size_t plane = planeIndices(direction, nxl, nyl).recv;
  Buffer buf = recvBuf_[static_cast<std::size_t>(direction)];

  if (entry.is3d) {
    Field3<real_t> f = entry.f3;
    Kokkos::parallel_for(
        "halo_unpack3", Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>>({0, 0}, {rowCells, nzS}),
        KOKKOS_LAMBDA(const std::size_t r, const std::size_t k) {
          const std::size_t jj = alongX ? r + rowStart : plane;
          const std::size_t ii = alongX ? plane : r + rowStart;
          f(jj, ii, k) = buf(offset + r * nzS + k);
        });
  } else {
    Field2<real_t> f = entry.f2;
    Kokkos::parallel_for(
        "halo_unpack2", Kokkos::RangePolicy<ExecSpace>(0, rowCells),
        KOKKOS_LAMBDA(const std::size_t r) {
          const std::size_t jj = alongX ? r + rowStart : plane;
          const std::size_t ii = alongX ? plane : r + rowStart;
          f(jj, ii) = buf(offset + r);
        });
  }
}

void HaloExchanger::exchangeAll() {
  std::vector<std::size_t> all(entries_.size());
  std::iota(all.begin(), all.end(), std::size_t{0});
  exchangeSelected(all, 0, 4, false);
}

std::vector<std::size_t> HaloExchanger::resolveNames(
    const std::vector<std::string>& names) const {
  std::vector<std::size_t> selected;
  selected.reserve(names.size());
  for (const std::string& name : names) {
    const auto it = std::find_if(entries_.begin(), entries_.end(),
                                 [&name](const Entry& e) { return e.name == name; });
    if (it == entries_.end()) {
      log::fatal(log::msg() << "HaloExchanger::exchange: field '" << name
                            << "' is not registered");
    }
    selected.push_back(static_cast<std::size_t>(it - entries_.begin()));
  }
  return selected;
}

void HaloExchanger::exchange(const std::vector<std::string>& names) {
  exchangeSelected(resolveNames(names), 0, 4, false);
}

void HaloExchanger::exchangeWithCorners(const std::vector<std::string>& names) {
  const std::vector<std::size_t> selected = resolveNames(names);
  // Phase 1 settles the i-halo columns; phase 2 sends full-width rows so the
  // corner cells travel inside the regular south/north messages. The
  // machinery is rank-generic: P1's uy/vx interpolation uses it for 2D
  // fields (amendment A1) and P4's dispersion cross terms for the 3D
  // subsurface scalar.
  exchangeSelected(selected, kWest, kEast + 1, false);
  exchangeSelected(selected, kSouth, kNorth + 1, true);
}

void HaloExchanger::exchangeSelected(const std::vector<std::size_t>& selected, int dirBegin,
                                     int dirEnd, bool wideRows) {
  std::array<std::size_t, 4> messageCount = {0, 0, 0, 0};
  for (int dir = dirBegin; dir < dirEnd; ++dir) {
    for (const std::size_t idx : selected) {
      messageCount[static_cast<std::size_t>(dir)] += planeCount(entries_[idx], dir, wideRows);
    }
  }

  // Post receives first (device buffers when GPU-aware, host mirrors else).
  std::array<MPI_Request, 8> requests;
  requests.fill(MPI_REQUEST_NULL);
  for (int dir = dirBegin; dir < dirEnd; ++dir) {
    const std::size_t d = static_cast<std::size_t>(dir);
    real_t* recvPtr = gpuAware_ ? recvBuf_[d].data() : recvHost_[d].data();
    MPI_Irecv(recvPtr, static_cast<int>(messageCount[d]), MPI_DOUBLE, neighbor_[d], 200 + dir,
              grid_.comm(), &requests[d]);
  }

  // Pack the selected fields for the active directions, then fence once.
  for (int dir = dirBegin; dir < dirEnd; ++dir) {
    std::size_t offset = 0;
    for (const std::size_t idx : selected) {
      packEntry(entries_[idx], dir, offset, wideRows);
      offset += planeCount(entries_[idx], dir, wideRows);
    }
  }
  Kokkos::fence();

  for (int dir = dirBegin; dir < dirEnd; ++dir) {
    const std::size_t d = static_cast<std::size_t>(dir);
    real_t* sendPtr = sendBuf_[d].data();
    if (!gpuAware_) {
      auto devSlice = Kokkos::subview(sendBuf_[d], std::make_pair(std::size_t{0}, messageCount[d]));
      auto hostSlice =
          Kokkos::subview(sendHost_[d], std::make_pair(std::size_t{0}, messageCount[d]));
      Kokkos::deep_copy(hostSlice, devSlice);
      sendPtr = sendHost_[d].data();
    }
    // Messages cross between opposite directions: what we send west arrives
    // at the neighbor as its east-received message, so tags pair (0<->1),
    // (2<->3).
    const int pairTag = 200 + (dir ^ 1);
    MPI_Isend(sendPtr, static_cast<int>(messageCount[d]), MPI_DOUBLE, neighbor_[d], pairTag,
              grid_.comm(), &requests[4 + d]);
  }

  MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);

  for (int dir = dirBegin; dir < dirEnd; ++dir) {
    const std::size_t d = static_cast<std::size_t>(dir);
    if (neighbor_[d] == MPI_PROC_NULL) {
      continue;  // domain edge: halos belong to the BC system
    }
    if (!gpuAware_) {
      auto devSlice = Kokkos::subview(recvBuf_[d], std::make_pair(std::size_t{0}, messageCount[d]));
      auto hostSlice =
          Kokkos::subview(recvHost_[d], std::make_pair(std::size_t{0}, messageCount[d]));
      Kokkos::deep_copy(devSlice, hostSlice);
    }
    std::size_t offset = 0;
    for (const std::size_t idx : selected) {
      unpackEntry(entries_[idx], dir, offset, wideRows);
      offset += planeCount(entries_[idx], dir, wideRows);
    }
  }
  Kokkos::fence();
}

}  // namespace frehg
