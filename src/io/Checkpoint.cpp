/// \file Checkpoint.cpp
/// \brief Implementation of checkpoint write/read.

#include "io/Checkpoint.hpp"

#include "core/Logger.hpp"

namespace frehg::io {

namespace {

/// Copy a field interior (halo layout) into a compact host block.
HostField2<real_t> interiorOf(const Grid& grid, const Field2<real_t>& field) {
  const std::size_t nyl = static_cast<std::size_t>(grid.nyLocal());
  const std::size_t nxl = static_cast<std::size_t>(grid.nxLocal());
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  HostField2<real_t> interior("ckpt_interior2", nyl, nxl);
  for (std::size_t j = 0; j < nyl; ++j) {
    for (std::size_t i = 0; i < nxl; ++i) {
      interior(j, i) = host(j + 1, i + 1);
    }
  }
  return interior;
}

/// \copydoc interiorOf
HostField3<real_t> interiorOf3(const Grid& grid, const Field3<real_t>& field) {
  const std::size_t nyl = static_cast<std::size_t>(grid.nyLocal());
  const std::size_t nxl = static_cast<std::size_t>(grid.nxLocal());
  const std::size_t nzg = static_cast<std::size_t>(grid.nz());
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);
  HostField3<real_t> interior("ckpt_interior3", nyl, nxl, nzg);
  for (std::size_t j = 0; j < nyl; ++j) {
    for (std::size_t i = 0; i < nxl; ++i) {
      for (std::size_t k = 0; k < nzg; ++k) {
        interior(j, i, k) = host(j + 1, i + 1, k);
      }
    }
  }
  return interior;
}

}  // namespace

void Checkpoint::write(real_t t, long step, const Scalars& scalars, const Fields2& fields2,
                       const Fields3& fields3, real_t labelTime) {
  const Grid& grid = out_.grid();
  const std::string groupPath =
      "/checkpoint/" + Hdf5Output::timeKey(labelTime >= 0.0 ? labelTime : t);
  out_.ensureGroup(groupPath);
  out_.writeDoubleAttribute(groupPath, "t", t);
  out_.writeLongAttribute(groupPath, "step", step);
  for (const auto& [name, value] : scalars) {
    out_.writeDoubleAttribute(groupPath, "scalar_" + name, value);
  }
  {
    std::string names;
    for (const auto& [name, value] : scalars) {
      names += names.empty() ? name : "," + name;
    }
    out_.writeStringAttribute(groupPath, "scalar_names", names);
  }

  for (const auto& [name, field] : fields2) {
    out_.writeDistributed2(groupPath + "/" + name, interiorOf(grid, field));
  }
  for (const auto& [name, field] : fields3) {
    out_.writeDistributed3(groupPath + "/" + name, interiorOf3(grid, field));
  }
  out_.flush();
}

bool Checkpoint::has(real_t t) const {
  return out_.exists("/checkpoint/" + Hdf5Output::timeKey(t));
}

Checkpoint::Header Checkpoint::read(real_t t, const Fields2& fields2, const Fields3& fields3) {
  const Grid& grid = out_.grid();
  const std::string groupPath = "/checkpoint/" + Hdf5Output::timeKey(t);
  if (!out_.exists(groupPath)) {
    log::fatal(log::msg() << "Checkpoint: no checkpoint at t = " << t << " (group '" << groupPath
                          << "') in the output file");
  }

  Header header;
  header.t = out_.readDoubleAttribute(groupPath, "t");
  header.step = out_.readLongAttribute(groupPath, "step");
  {
    const std::string names = out_.readStringAttribute(groupPath, "scalar_names");
    std::size_t begin = 0;
    while (begin < names.size()) {
      std::size_t end = names.find(',', begin);
      if (end == std::string::npos) {
        end = names.size();
      }
      const std::string name = names.substr(begin, end - begin);
      if (!name.empty()) {
        header.scalars[name] = out_.readDoubleAttribute(groupPath, "scalar_" + name);
      }
      begin = end + 1;
    }
  }

  const std::size_t nyl = static_cast<std::size_t>(grid.nyLocal());
  const std::size_t nxl = static_cast<std::size_t>(grid.nxLocal());
  const std::size_t nzg = static_cast<std::size_t>(grid.nz());

  for (const auto& [name, field] : fields2) {
    HostField2<real_t> interior("ckpt_read2", nyl, nxl);
    out_.readDistributed2(groupPath + "/" + name, interior);
    auto host = Kokkos::create_mirror_view(field);
    Kokkos::deep_copy(host, field);
    for (std::size_t j = 0; j < nyl; ++j) {
      for (std::size_t i = 0; i < nxl; ++i) {
        host(j + 1, i + 1) = interior(j, i);
      }
    }
    Kokkos::deep_copy(field, host);
  }
  for (const auto& [name, field] : fields3) {
    HostField3<real_t> interior("ckpt_read3", nyl, nxl, nzg);
    out_.readDistributed3(groupPath + "/" + name, interior);
    auto host = Kokkos::create_mirror_view(field);
    Kokkos::deep_copy(host, field);
    for (std::size_t j = 0; j < nyl; ++j) {
      for (std::size_t i = 0; i < nxl; ++i) {
        for (std::size_t k = 0; k < nzg; ++k) {
          host(j + 1, i + 1, k) = interior(j, i, k);
        }
      }
    }
    Kokkos::deep_copy(field, host);
  }
  Kokkos::fence();
  return header;
}

}  // namespace frehg::io
