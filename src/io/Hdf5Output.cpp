/// \file Hdf5Output.cpp
/// \brief Implementation of the parallel HDF5 output file.

#include "io/Hdf5Output.hpp"

#include "core/Logger.hpp"

#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <limits>
#include <sstream>
#include <vector>

/// Every HDF5 return code is checked; failures are fatal (plan §11.1 rule 5).
#define FREHG_H5_CHECK(call)                                                          \
  do {                                                                                \
    if ((call) < 0) {                                                                 \
      ::frehg::log::fatal(::frehg::log::msg() << "HDF5 error from " << #call);        \
    }                                                                                 \
  } while (0)

namespace frehg::io {

namespace {

/// Check an hid_t-returning call and pass the id through.
hid_t checkedId(hid_t id, const char* what) {
  if (id < 0) {
    log::fatal(log::msg() << "HDF5 error from " << what);
  }
  return id;
}

/// Property list for collective transfers.
hid_t collectiveXfer() {
  const hid_t xfer = checkedId(H5Pcreate(H5P_DATASET_XFER), "H5Pcreate(xfer)");
  FREHG_H5_CHECK(H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE));
  return xfer;
}

std::string isoTimestamp() {
  const std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm tmValue = {};
  gmtime_r(&now, &tmValue);
  std::ostringstream out;
  out << (tmValue.tm_year + 1900) << "-";
  out.width(2);
  out.fill('0');
  out << (tmValue.tm_mon + 1) << "-";
  out.width(2);
  out << tmValue.tm_mday << "T";
  out.width(2);
  out << tmValue.tm_hour << ":";
  out.width(2);
  out << tmValue.tm_min << ":";
  out.width(2);
  out << tmValue.tm_sec << "Z";
  return out.str();
}

}  // namespace

std::string Hdf5Output::timeKey(real_t t) {
  const long long seconds = std::llround(t);
  if (std::abs(t - static_cast<real_t>(seconds)) > 1.0e-6) {
    log::fatal(log::msg() << "Hdf5Output: output time " << t
                          << " s does not land on a whole second; §7 keys times by integer "
                             "seconds and the configuration validator enforces integral "
                             "output intervals");
  }
  return std::to_string(seconds);
}

Hdf5Output::Hdf5Output(const Grid& grid, const std::string& filename,
                       const std::string& configText)
    : grid_(grid), filename_(filename) {
  const std::filesystem::path parent = std::filesystem::path(filename).parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  const hid_t fapl = checkedId(H5Pcreate(H5P_FILE_ACCESS), "H5Pcreate(fapl)");
  FREHG_H5_CHECK(H5Pset_fapl_mpio(fapl, grid_.comm(), MPI_INFO_NULL));
  file_ = checkedId(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl),
                    "H5Fcreate");
  FREHG_H5_CHECK(H5Pclose(fapl));
  writeRootAttributes(configText);
}

Hdf5Output::Hdf5Output(const Grid& grid, const std::string& filename)
    : grid_(grid), filename_(filename) {
  const hid_t fapl = checkedId(H5Pcreate(H5P_FILE_ACCESS), "H5Pcreate(fapl)");
  FREHG_H5_CHECK(H5Pset_fapl_mpio(fapl, grid_.comm(), MPI_INFO_NULL));
  file_ = checkedId(H5Fopen(filename.c_str(), H5F_ACC_RDWR, fapl), "H5Fopen");
  FREHG_H5_CHECK(H5Pclose(fapl));
}

Hdf5Output::~Hdf5Output() {
  if (file_ >= 0) {
    H5Fclose(file_);
  }
}

void Hdf5Output::writeRootAttributes(const std::string& configText) {
  ensureGroup("/frehg2");
  writeStringAttribute("/frehg2", "version", FREHG_VERSION);
  writeStringAttribute("/frehg2", "git_sha", FREHG_GIT_SHA);
  writeStringAttribute("/frehg2", "config", configText);
  writeStringAttribute("/frehg2", "created", isoTimestamp());
}

void Hdf5Output::ensureGroup(const std::string& path) {
  std::string sofar;
  std::istringstream parts(path);
  std::string part;
  while (std::getline(parts, part, '/')) {
    if (part.empty()) {
      continue;
    }
    sofar += "/" + part;
    const htri_t present = H5Lexists(file_, sofar.c_str(), H5P_DEFAULT);
    FREHG_H5_CHECK(present);
    if (present == 0) {
      const hid_t group = checkedId(
          H5Gcreate2(file_, sofar.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), "H5Gcreate2");
      FREHG_H5_CHECK(H5Gclose(group));
    }
  }
}

bool Hdf5Output::exists(const std::string& path) const {
  // Check every intermediate link so a missing parent reads as "absent"
  // rather than raising an HDF5 error.
  std::string sofar;
  std::istringstream parts(path);
  std::string part;
  while (std::getline(parts, part, '/')) {
    if (part.empty()) {
      continue;
    }
    sofar += "/" + part;
    const htri_t present = H5Lexists(file_, sofar.c_str(), H5P_DEFAULT);
    FREHG_H5_CHECK(present);
    if (present == 0) {
      return false;
    }
  }
  return true;
}

hid_t Hdf5Output::createDataset1D(const std::string& datasetPath, hsize_t globalSize) {
  const std::size_t slash = datasetPath.find_last_of('/');
  if (slash != std::string::npos && slash > 0) {
    ensureGroup(datasetPath.substr(0, slash));
  }
  const hid_t space = checkedId(H5Screate_simple(1, &globalSize, nullptr), "H5Screate_simple");
  const hid_t dset = checkedId(H5Dcreate2(file_, datasetPath.c_str(), H5T_NATIVE_DOUBLE, space,
                                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                               "H5Dcreate2");
  FREHG_H5_CHECK(H5Sclose(space));
  return dset;
}

void Hdf5Output::writeDistributed2(const std::string& datasetPath,
                                   const HostField2<real_t>& interior) {
  const hsize_t nxl = static_cast<hsize_t>(grid_.nxLocal());
  const hsize_t nyl = static_cast<hsize_t>(grid_.nyLocal());
  if (interior.extent(0) != nyl || interior.extent(1) != nxl) {
    log::fatal(log::msg() << "writeDistributed2('" << datasetPath << "'): block is "
                          << interior.extent(0) << " x " << interior.extent(1) << ", expected "
                          << nyl << " x " << nxl);
  }
  const hsize_t globalSize =
      static_cast<hsize_t>(grid_.nx()) * static_cast<hsize_t>(grid_.ny());
  const hid_t dset = createDataset1D(datasetPath, globalSize);

  // Regular 1D hyperslab: this rank's rows.
  const hsize_t start = static_cast<hsize_t>(grid_.j0()) * static_cast<hsize_t>(grid_.nx()) +
                        static_cast<hsize_t>(grid_.i0());
  const hsize_t stride = static_cast<hsize_t>(grid_.nx());
  const hsize_t count = nyl;
  const hsize_t block = nxl;
  const hid_t fileSpace = checkedId(H5Dget_space(dset), "H5Dget_space");
  FREHG_H5_CHECK(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &start, &stride, &count, &block));
  const hsize_t localSize = nyl * nxl;
  const hid_t memSpace = checkedId(H5Screate_simple(1, &localSize, nullptr), "H5Screate_simple");

  std::vector<double> buffer(static_cast<std::size_t>(localSize));
  for (hsize_t j = 0; j < nyl; ++j) {
    for (hsize_t i = 0; i < nxl; ++i) {
      buffer[static_cast<std::size_t>(j * nxl + i)] =
          interior(static_cast<std::size_t>(j), static_cast<std::size_t>(i));
    }
  }
  const hid_t xfer = collectiveXfer();
  FREHG_H5_CHECK(H5Dwrite(dset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, xfer, buffer.data()));
  FREHG_H5_CHECK(H5Pclose(xfer));
  FREHG_H5_CHECK(H5Sclose(memSpace));
  FREHG_H5_CHECK(H5Sclose(fileSpace));
  FREHG_H5_CHECK(H5Dclose(dset));
}

void Hdf5Output::writeDistributed3(const std::string& datasetPath,
                                   const HostField3<real_t>& interior) {
  const hsize_t nxl = static_cast<hsize_t>(grid_.nxLocal());
  const hsize_t nyl = static_cast<hsize_t>(grid_.nyLocal());
  const hsize_t nzg = static_cast<hsize_t>(grid_.nz());
  if (interior.extent(0) != nyl || interior.extent(1) != nxl || interior.extent(2) != nzg) {
    log::fatal(log::msg() << "writeDistributed3('" << datasetPath << "'): block is "
                          << interior.extent(0) << " x " << interior.extent(1) << " x "
                          << interior.extent(2) << ", expected " << nyl << " x " << nxl << " x "
                          << nzg);
  }
  const hsize_t globalSize =
      static_cast<hsize_t>(grid_.nx()) * static_cast<hsize_t>(grid_.ny()) * nzg;
  const hid_t dset = createDataset1D(datasetPath, globalSize);

  const hsize_t start = (static_cast<hsize_t>(grid_.j0()) * static_cast<hsize_t>(grid_.nx()) +
                         static_cast<hsize_t>(grid_.i0())) *
                        nzg;
  const hsize_t stride = static_cast<hsize_t>(grid_.nx()) * nzg;
  const hsize_t count = nyl;
  const hsize_t block = nxl * nzg;
  const hid_t fileSpace = checkedId(H5Dget_space(dset), "H5Dget_space");
  FREHG_H5_CHECK(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &start, &stride, &count, &block));
  const hsize_t localSize = nyl * nxl * nzg;
  const hid_t memSpace = checkedId(H5Screate_simple(1, &localSize, nullptr), "H5Screate_simple");

  std::vector<double> buffer(static_cast<std::size_t>(localSize));
  std::size_t n = 0;
  for (hsize_t j = 0; j < nyl; ++j) {
    for (hsize_t i = 0; i < nxl; ++i) {
      for (hsize_t k = 0; k < nzg; ++k) {
        buffer[n++] = interior(static_cast<std::size_t>(j), static_cast<std::size_t>(i),
                               static_cast<std::size_t>(k));
      }
    }
  }
  const hid_t xfer = collectiveXfer();
  FREHG_H5_CHECK(H5Dwrite(dset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, xfer, buffer.data()));
  FREHG_H5_CHECK(H5Pclose(xfer));
  FREHG_H5_CHECK(H5Sclose(memSpace));
  FREHG_H5_CHECK(H5Sclose(fileSpace));
  FREHG_H5_CHECK(H5Dclose(dset));
}

namespace {

/// Shared hyperslab reader used by both readDistributed variants.
void readBlock(hid_t file, const Grid& grid, const std::string& datasetPath, hsize_t perCell,
               std::vector<double>& buffer) {
  const hsize_t nxl = static_cast<hsize_t>(grid.nxLocal());
  const hsize_t nyl = static_cast<hsize_t>(grid.nyLocal());
  const hid_t dset = checkedId(H5Dopen2(file, datasetPath.c_str(), H5P_DEFAULT), "H5Dopen2");
  const hsize_t start = (static_cast<hsize_t>(grid.j0()) * static_cast<hsize_t>(grid.nx()) +
                         static_cast<hsize_t>(grid.i0())) *
                        perCell;
  const hsize_t stride = static_cast<hsize_t>(grid.nx()) * perCell;
  const hsize_t count = nyl;
  const hsize_t block = nxl * perCell;
  const hid_t fileSpace = checkedId(H5Dget_space(dset), "H5Dget_space");
  FREHG_H5_CHECK(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &start, &stride, &count, &block));
  const hsize_t localSize = nyl * nxl * perCell;
  const hid_t memSpace = checkedId(H5Screate_simple(1, &localSize, nullptr), "H5Screate_simple");
  buffer.resize(static_cast<std::size_t>(localSize));
  const hid_t xfer = collectiveXfer();
  FREHG_H5_CHECK(H5Dread(dset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, xfer, buffer.data()));
  FREHG_H5_CHECK(H5Pclose(xfer));
  FREHG_H5_CHECK(H5Sclose(memSpace));
  FREHG_H5_CHECK(H5Sclose(fileSpace));
  FREHG_H5_CHECK(H5Dclose(dset));
}

}  // namespace

void Hdf5Output::readDistributed2(const std::string& datasetPath,
                                  HostField2<real_t>& interior) const {
  const std::size_t nxl = static_cast<std::size_t>(grid_.nxLocal());
  const std::size_t nyl = static_cast<std::size_t>(grid_.nyLocal());
  if (interior.extent(0) != nyl || interior.extent(1) != nxl) {
    log::fatal(log::msg() << "readDistributed2('" << datasetPath << "'): block is "
                          << interior.extent(0) << " x " << interior.extent(1) << ", expected "
                          << nyl << " x " << nxl);
  }
  std::vector<double> buffer;
  readBlock(file_, grid_, datasetPath, 1, buffer);
  for (std::size_t j = 0; j < nyl; ++j) {
    for (std::size_t i = 0; i < nxl; ++i) {
      interior(j, i) = buffer[j * nxl + i];
    }
  }
}

void Hdf5Output::readDistributed3(const std::string& datasetPath,
                                  HostField3<real_t>& interior) const {
  const std::size_t nxl = static_cast<std::size_t>(grid_.nxLocal());
  const std::size_t nyl = static_cast<std::size_t>(grid_.nyLocal());
  const std::size_t nzg = static_cast<std::size_t>(grid_.nz());
  if (interior.extent(0) != nyl || interior.extent(1) != nxl || interior.extent(2) != nzg) {
    log::fatal(log::msg() << "readDistributed3('" << datasetPath << "'): block is "
                          << interior.extent(0) << " x " << interior.extent(1) << " x "
                          << interior.extent(2) << ", expected " << nyl << " x " << nxl << " x "
                          << nzg);
  }
  std::vector<double> buffer;
  readBlock(file_, grid_, datasetPath, static_cast<hsize_t>(nzg), buffer);
  std::size_t n = 0;
  for (std::size_t j = 0; j < nyl; ++j) {
    for (std::size_t i = 0; i < nxl; ++i) {
      for (std::size_t k = 0; k < nzg; ++k) {
        interior(j, i, k) = buffer[n++];
      }
    }
  }
}

void Hdf5Output::writeField2(const std::string& group, const std::string& var, real_t t,
                             const Field2<real_t>& field, const VarMeta& meta) {
  const std::size_t nxl = static_cast<std::size_t>(grid_.nxLocal());
  const std::size_t nyl = static_cast<std::size_t>(grid_.nyLocal());
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);

  HostField2<real_t> interior("field2_interior", nyl, nxl);
  const auto& ktop = grid_.ktopHost();
  for (std::size_t j = 0; j < nyl; ++j) {
    for (std::size_t i = 0; i < nxl; ++i) {
      const bool inactive = ktop(j, i) >= grid_.nz();
      interior(j, i) =
          inactive ? std::numeric_limits<real_t>::quiet_NaN() : host(j + 1, i + 1);
    }
  }
  const std::string datasetPath = "/" + group + "/" + var + "/" + timeKey(t);
  writeDistributed2(datasetPath, interior);
  writeStringAttribute(datasetPath, "units", meta.units);
  writeStringAttribute(datasetPath, "long_name", meta.longName);
  writeDoubleAttribute(datasetPath, "time", t);
}

void Hdf5Output::writeField3(const std::string& group, const std::string& var, real_t t,
                             const Field3<real_t>& field, const VarMeta& meta) {
  const std::size_t nxl = static_cast<std::size_t>(grid_.nxLocal());
  const std::size_t nyl = static_cast<std::size_t>(grid_.nyLocal());
  const std::size_t nzg = static_cast<std::size_t>(grid_.nz());
  auto host = Kokkos::create_mirror_view(field);
  Kokkos::deep_copy(host, field);

  HostField3<real_t> interior("field3_interior", nyl, nxl, nzg);
  const auto& ktop = grid_.ktopHost();
  for (std::size_t j = 0; j < nyl; ++j) {
    for (std::size_t i = 0; i < nxl; ++i) {
      const std::size_t kt = static_cast<std::size_t>(ktop(j, i));
      for (std::size_t k = 0; k < nzg; ++k) {
        interior(j, i, k) =
            (k < kt) ? std::numeric_limits<real_t>::quiet_NaN() : host(j + 1, i + 1, k);
      }
    }
  }
  const std::string datasetPath = "/" + group + "/" + var + "/" + timeKey(t);
  writeDistributed3(datasetPath, interior);
  writeStringAttribute(datasetPath, "units", meta.units);
  writeStringAttribute(datasetPath, "long_name", meta.longName);
  writeDoubleAttribute(datasetPath, "time", t);
}

void Hdf5Output::writeGridMeta(const HostField2<real_t>& bottomInterior) {
  // Uniform-layer zcell: every column carries the same level centers.
  HostField3<real_t> zcell("gridmeta_zcell", static_cast<std::size_t>(grid_.nyLocal()),
                           static_cast<std::size_t>(grid_.nxLocal()),
                           static_cast<std::size_t>(grid_.nz()));
  for (std::size_t j = 0; j < zcell.extent(0); ++j) {
    for (std::size_t i = 0; i < zcell.extent(1); ++i) {
      for (int k = 0; k < grid_.nz(); ++k) {
        zcell(j, i, static_cast<std::size_t>(k)) = grid_.zCenter(k);
      }
    }
  }
  writeGridMeta(bottomInterior, zcell);
}

void Hdf5Output::writeGridMeta(const HostField2<real_t>& bottomInterior,
                               const HostField3<real_t>& zcellInterior) {
  ensureGroup("/grid");

  // Coordinates are global 1D axes; write them from the collective
  // distributed path used everywhere else would shard them, so build the
  // full axes on every rank and write with all ranks selecting the whole
  // dataset (identical data, collective).
  {
    std::vector<double> xc(static_cast<std::size_t>(grid_.nx()));
    for (int i = 0; i < grid_.nx(); ++i) {
      xc[static_cast<std::size_t>(i)] = grid_.xCenter(i);
    }
    const hsize_t size = static_cast<hsize_t>(xc.size());
    const hid_t dset = createDataset1D("/grid/x_center", size);
    const hid_t xfer = collectiveXfer();
    FREHG_H5_CHECK(H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xfer, xc.data()));
    FREHG_H5_CHECK(H5Pclose(xfer));
    FREHG_H5_CHECK(H5Dclose(dset));
  }
  {
    std::vector<double> yc(static_cast<std::size_t>(grid_.ny()));
    for (int j = 0; j < grid_.ny(); ++j) {
      yc[static_cast<std::size_t>(j)] = grid_.yCenter(j);
    }
    const hsize_t size = static_cast<hsize_t>(yc.size());
    const hid_t dset = createDataset1D("/grid/y_center", size);
    const hid_t xfer = collectiveXfer();
    FREHG_H5_CHECK(H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xfer, yc.data()));
    FREHG_H5_CHECK(H5Pclose(xfer));
    FREHG_H5_CHECK(H5Dclose(dset));
  }
  {
    std::vector<double> zc(static_cast<std::size_t>(grid_.nz()));
    std::vector<double> dz(static_cast<std::size_t>(grid_.nz()));
    for (int k = 0; k < grid_.nz(); ++k) {
      zc[static_cast<std::size_t>(k)] = grid_.zCenter(k);
      dz[static_cast<std::size_t>(k)] = grid_.dzK(k);
    }
    const hsize_t size = static_cast<hsize_t>(zc.size());
    for (const auto& [name, data] :
         {std::pair<const char*, std::vector<double>*>{"/grid/z_center", &zc},
          std::pair<const char*, std::vector<double>*>{"/grid/dz", &dz}}) {
      const hid_t dset = createDataset1D(name, size);
      const hid_t xfer = collectiveXfer();
      FREHG_H5_CHECK(H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xfer, data->data()));
      FREHG_H5_CHECK(H5Pclose(xfer));
      FREHG_H5_CHECK(H5Dclose(dset));
    }
  }

  // Per-cell center elevations in the (j*NX + i)*NZ + k layout the b3
  // analysis script reads (plan §7 contract).
  writeDistributed3("/groundwater/zcell/0", zcellInterior);

  writeDistributed2("/grid/bottom", bottomInterior);

  HostField2<real_t> ktopReal("ktop_real", static_cast<std::size_t>(grid_.nyLocal()),
                              static_cast<std::size_t>(grid_.nxLocal()));
  const auto& ktop = grid_.ktopHost();
  for (std::size_t j = 0; j < ktopReal.extent(0); ++j) {
    for (std::size_t i = 0; i < ktopReal.extent(1); ++i) {
      ktopReal(j, i) = static_cast<real_t>(ktop(j, i));
    }
  }
  writeDistributed2("/grid/ktop", ktopReal);
}

std::vector<real_t> Hdf5Output::readAll1D(const std::string& path) const {
  const hid_t dset = checkedId(H5Dopen2(file_, path.c_str(), H5P_DEFAULT), "H5Dopen2");
  const hid_t space = checkedId(H5Dget_space(dset), "H5Dget_space");
  hsize_t size = 0;
  FREHG_H5_CHECK(H5Sget_simple_extent_dims(space, &size, nullptr));
  std::vector<real_t> data(static_cast<std::size_t>(size));
  const hid_t xfer = collectiveXfer();
  FREHG_H5_CHECK(H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xfer, data.data()));
  FREHG_H5_CHECK(H5Pclose(xfer));
  FREHG_H5_CHECK(H5Sclose(space));
  FREHG_H5_CHECK(H5Dclose(dset));
  return data;
}

void Hdf5Output::writeStringAttribute(const std::string& objectPath, const std::string& name,
                                      const std::string& value) {
  const hid_t obj = checkedId(H5Oopen(file_, objectPath.c_str(), H5P_DEFAULT), "H5Oopen");
  const hid_t type = checkedId(H5Tcopy(H5T_C_S1), "H5Tcopy");
  FREHG_H5_CHECK(H5Tset_size(type, value.empty() ? 1 : value.size()));
  FREHG_H5_CHECK(H5Tset_strpad(type, H5T_STR_NULLPAD));
  const hid_t space = checkedId(H5Screate(H5S_SCALAR), "H5Screate");
  if (H5Aexists(obj, name.c_str()) > 0) {
    FREHG_H5_CHECK(H5Adelete(obj, name.c_str()));
  }
  const hid_t attr = checkedId(
      H5Acreate2(obj, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT), "H5Acreate2");
  const std::string padded = value.empty() ? std::string(1, '\0') : value;
  FREHG_H5_CHECK(H5Awrite(attr, type, padded.data()));
  FREHG_H5_CHECK(H5Aclose(attr));
  FREHG_H5_CHECK(H5Sclose(space));
  FREHG_H5_CHECK(H5Tclose(type));
  FREHG_H5_CHECK(H5Oclose(obj));
}

void Hdf5Output::writeDoubleAttribute(const std::string& objectPath, const std::string& name,
                                      double value) {
  const hid_t obj = checkedId(H5Oopen(file_, objectPath.c_str(), H5P_DEFAULT), "H5Oopen");
  const hid_t space = checkedId(H5Screate(H5S_SCALAR), "H5Screate");
  if (H5Aexists(obj, name.c_str()) > 0) {
    FREHG_H5_CHECK(H5Adelete(obj, name.c_str()));
  }
  const hid_t attr = checkedId(
      H5Acreate2(obj, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT),
      "H5Acreate2");
  FREHG_H5_CHECK(H5Awrite(attr, H5T_NATIVE_DOUBLE, &value));
  FREHG_H5_CHECK(H5Aclose(attr));
  FREHG_H5_CHECK(H5Sclose(space));
  FREHG_H5_CHECK(H5Oclose(obj));
}

void Hdf5Output::writeLongAttribute(const std::string& objectPath, const std::string& name,
                                    long value) {
  const hid_t obj = checkedId(H5Oopen(file_, objectPath.c_str(), H5P_DEFAULT), "H5Oopen");
  const hid_t space = checkedId(H5Screate(H5S_SCALAR), "H5Screate");
  if (H5Aexists(obj, name.c_str()) > 0) {
    FREHG_H5_CHECK(H5Adelete(obj, name.c_str()));
  }
  const hid_t attr = checkedId(
      H5Acreate2(obj, name.c_str(), H5T_NATIVE_LONG, space, H5P_DEFAULT, H5P_DEFAULT),
      "H5Acreate2");
  FREHG_H5_CHECK(H5Awrite(attr, H5T_NATIVE_LONG, &value));
  FREHG_H5_CHECK(H5Aclose(attr));
  FREHG_H5_CHECK(H5Sclose(space));
  FREHG_H5_CHECK(H5Oclose(obj));
}

std::string Hdf5Output::readStringAttribute(const std::string& objectPath,
                                            const std::string& name) const {
  const hid_t obj = checkedId(H5Oopen(file_, objectPath.c_str(), H5P_DEFAULT), "H5Oopen");
  const hid_t attr = checkedId(H5Aopen(obj, name.c_str(), H5P_DEFAULT), "H5Aopen");
  const hid_t type = checkedId(H5Aget_type(attr), "H5Aget_type");
  const std::size_t size = H5Tget_size(type);
  std::vector<char> data(size);
  const hid_t memType = checkedId(H5Tcopy(H5T_C_S1), "H5Tcopy");
  FREHG_H5_CHECK(H5Tset_size(memType, size));
  FREHG_H5_CHECK(H5Tset_strpad(memType, H5T_STR_NULLPAD));
  FREHG_H5_CHECK(H5Aread(attr, memType, data.data()));
  FREHG_H5_CHECK(H5Tclose(memType));
  FREHG_H5_CHECK(H5Tclose(type));
  FREHG_H5_CHECK(H5Aclose(attr));
  FREHG_H5_CHECK(H5Oclose(obj));
  std::string value(data.begin(), data.end());
  while (!value.empty() && value.back() == '\0') {
    value.pop_back();
  }
  return value;
}

double Hdf5Output::readDoubleAttribute(const std::string& objectPath,
                                       const std::string& name) const {
  const hid_t obj = checkedId(H5Oopen(file_, objectPath.c_str(), H5P_DEFAULT), "H5Oopen");
  const hid_t attr = checkedId(H5Aopen(obj, name.c_str(), H5P_DEFAULT), "H5Aopen");
  double value = 0;
  FREHG_H5_CHECK(H5Aread(attr, H5T_NATIVE_DOUBLE, &value));
  FREHG_H5_CHECK(H5Aclose(attr));
  FREHG_H5_CHECK(H5Oclose(obj));
  return value;
}

long Hdf5Output::readLongAttribute(const std::string& objectPath, const std::string& name) const {
  const hid_t obj = checkedId(H5Oopen(file_, objectPath.c_str(), H5P_DEFAULT), "H5Oopen");
  const hid_t attr = checkedId(H5Aopen(obj, name.c_str(), H5P_DEFAULT), "H5Aopen");
  long value = 0;
  FREHG_H5_CHECK(H5Aread(attr, H5T_NATIVE_LONG, &value));
  FREHG_H5_CHECK(H5Aclose(attr));
  FREHG_H5_CHECK(H5Oclose(obj));
  return value;
}

void Hdf5Output::flush() { FREHG_H5_CHECK(H5Fflush(file_, H5F_SCOPE_GLOBAL)); }

}  // namespace frehg::io
