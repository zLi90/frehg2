/// \file Monitor.cpp
/// \brief Implementation of point-probe monitors.

#include "io/Monitor.hpp"

#include "core/Logger.hpp"

#include <hdf5.h>

#include <sstream>

/// Every HDF5 return code is checked; failures are fatal (plan §11.1 rule 5).
#define FREHG_H5_CHECK(call)                                                          \
  do {                                                                                \
    if ((call) < 0) {                                                                 \
      ::frehg::log::fatal(::frehg::log::msg() << "HDF5 error from " << #call);        \
    }                                                                                 \
  } while (0)

namespace frehg::io {

namespace {

hid_t checkedId(hid_t id, const char* what) {
  if (id < 0) {
    log::fatal(log::msg() << "HDF5 error from " << what);
  }
  return id;
}

}  // namespace

Monitor::Monitor(Hdf5Output& out, const MonitorConfig& config) : out_(out), config_(config) {
  const Grid& grid = out_.grid();
  if (config_.i < 0 || config_.i >= grid.nx() || config_.j < 0 || config_.j >= grid.ny()) {
    log::fatal(log::msg() << "Monitor '" << config_.name << "': cell (" << config_.i << ", "
                          << config_.j << ") outside the " << grid.nx() << " x " << grid.ny()
                          << " domain");
  }
  if (config_.variables.empty()) {
    log::fatal(log::msg() << "Monitor '" << config_.name << "': no variables configured");
  }
  columns_ = 1 + config_.variables.size();
  owns_ = config_.i >= grid.i0() && config_.i < grid.i0() + grid.nxLocal() &&
          config_.j >= grid.j0() && config_.j < grid.j0() + grid.nyLocal();

  // Create the extendable table collectively.
  out_.ensureGroup("/monitor");
  const std::string path = "/monitor/" + config_.name;
  const hsize_t dims[2] = {0, static_cast<hsize_t>(columns_)};
  const hsize_t maxDims[2] = {H5S_UNLIMITED, static_cast<hsize_t>(columns_)};
  const hid_t space = checkedId(H5Screate_simple(2, dims, maxDims), "H5Screate_simple");
  const hid_t dcpl = checkedId(H5Pcreate(H5P_DATASET_CREATE), "H5Pcreate(dcpl)");
  const hsize_t chunk[2] = {256, static_cast<hsize_t>(columns_)};
  FREHG_H5_CHECK(H5Pset_chunk(dcpl, 2, chunk));
  const hid_t dset = checkedId(H5Dcreate2(out_.fileId(), path.c_str(), H5T_NATIVE_DOUBLE, space,
                                          H5P_DEFAULT, dcpl, H5P_DEFAULT),
                               "H5Dcreate2");
  FREHG_H5_CHECK(H5Pclose(dcpl));
  FREHG_H5_CHECK(H5Sclose(space));
  FREHG_H5_CHECK(H5Dclose(dset));

  std::ostringstream columns;
  columns << "time";
  for (const std::string& var : config_.variables) {
    columns << "," << var;
  }
  out_.writeStringAttribute(path, "columns", columns.str());
  out_.writeLongAttribute(path, "i", config_.i);
  out_.writeLongAttribute(path, "j", config_.j);
}

void Monitor::record(real_t t, const std::vector<real_t>& values) {
  if (!owns_) {
    return;
  }
  if (values.size() != config_.variables.size()) {
    log::fatal(log::msg() << "Monitor '" << config_.name << "': " << values.size()
                          << " values recorded for " << config_.variables.size()
                          << " configured variables");
  }
  buffer_.push_back(t);
  buffer_.insert(buffer_.end(), values.begin(), values.end());
}

void Monitor::flush() {
  const Grid& grid = out_.grid();

  // Agree on the number of fresh rows: the owner knows, everyone extends.
  long localRows = owns_ ? static_cast<long>(buffer_.size() / columns_) : 0;
  long fresh = 0;
  MPI_Allreduce(&localRows, &fresh, 1, MPI_LONG, MPI_MAX, grid.comm());
  if (fresh == 0) {
    return;
  }

  const std::string path = "/monitor/" + config_.name;
  const hid_t dset =
      checkedId(H5Dopen2(out_.fileId(), path.c_str(), H5P_DEFAULT), "H5Dopen2");
  const hsize_t newDims[2] = {static_cast<hsize_t>(rowsWritten_ + fresh),
                              static_cast<hsize_t>(columns_)};
  FREHG_H5_CHECK(H5Dset_extent(dset, newDims));

  // Row data: written by the owner alone with an independent transfer; the
  // other ranks select nothing.
  const hid_t fileSpace = checkedId(H5Dget_space(dset), "H5Dget_space");
  if (owns_) {
    const hsize_t start[2] = {static_cast<hsize_t>(rowsWritten_), 0};
    const hsize_t count[2] = {static_cast<hsize_t>(fresh), static_cast<hsize_t>(columns_)};
    FREHG_H5_CHECK(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, nullptr, count, nullptr));
    const hid_t memSpace = checkedId(H5Screate_simple(2, count, nullptr), "H5Screate_simple");
    FREHG_H5_CHECK(
        H5Dwrite(dset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, buffer_.data()));
    FREHG_H5_CHECK(H5Sclose(memSpace));
  }
  FREHG_H5_CHECK(H5Sclose(fileSpace));
  FREHG_H5_CHECK(H5Dclose(dset));

  buffer_.clear();
  rowsWritten_ += fresh;
}

}  // namespace frehg::io
