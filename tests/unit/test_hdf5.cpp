/// \file test_hdf5.cpp
/// \brief HDF5 output tests: 0-ULP write/read round-trip, §7 layout-index
///        contract, NaN masking of inactive cells, attributes, monitors,
///        and checkpoint round-trip (plan §8.1).

#include "core/Grid.hpp"
#include "core/Types.hpp"
#include "io/Checkpoint.hpp"
#include "io/Hdf5Output.hpp"
#include "io/Monitor.hpp"

#include <gtest/gtest.h>
#include <hdf5.h>
#include <mpi.h>

#include <cmath>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

namespace {

using frehg::real_t;

/// Bit-exact double comparison (0 ULP).
bool bitEqual(real_t a, real_t b) { return std::memcmp(&a, &b, sizeof(real_t)) == 0; }

frehg::DomainConfig testDomain() {
  frehg::DomainConfig dom;
  dom.nx = 4;
  dom.ny = 3;
  dom.nz = 2;
  dom.dx = 0.5;
  dom.dy = 0.25;
  dom.dz = 0.1;
  return dom;
}

/// Deterministic, ULP-hostile cell values (irrational-ish fractions).
real_t pattern2(int j, int i) { return (1.0 / 3.0) * (j + 1) + (1.0 / 7.0) * (i + 1); }
real_t pattern3(int j, int i, int k) { return pattern2(j, i) + (1.0 / 11.0) * (k + 1); }

class Hdf5Test : public ::testing::Test {
 protected:
  std::filesystem::path dir_;

  void SetUp() override {
    dir_ = std::filesystem::path(::testing::TempDir()) / "frehg_hdf5_test";
    std::filesystem::create_directories(dir_);
  }
  void TearDown() override { std::filesystem::remove_all(dir_); }

  std::string file(const std::string& name) { return (dir_ / name).string(); }
};

TEST_F(Hdf5Test, TimeKeyFormat) {
  EXPECT_EQ(frehg::io::Hdf5Output::timeKey(1800.0), "1800");
  EXPECT_EQ(frehg::io::Hdf5Output::timeKey(0.0), "0");
  EXPECT_EQ(frehg::io::Hdf5Output::timeKey(36000.0000001), "36000");
  EXPECT_THROW(frehg::io::Hdf5Output::timeKey(10.5), frehg::FatalError);
}

TEST_F(Hdf5Test, Field2RoundTripLayoutAndMask) {
  frehg::Grid grid(MPI_COMM_SELF, testDomain());
  frehg::HostField2<int> ktop("ktop", 3, 4);
  Kokkos::deep_copy(ktop, 0);
  ktop(1, 2) = 2;  // fully inactive column -> NaN in the output
  grid.buildGlobalIds(ktop);

  frehg::Field2<real_t> eta("eta", 5, 6);
  auto etaH = Kokkos::create_mirror_view(eta);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 4; ++i) {
      etaH(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1) = pattern2(j, i);
    }
  }
  Kokkos::deep_copy(eta, etaH);

  {
    frehg::io::Hdf5Output out(grid, file("field2.h5"), "config-text");
    out.writeField2("surface", "eta", 1800.0, eta, {"m", "free-surface elevation"});
    EXPECT_TRUE(out.exists("/surface/eta/1800"));
    EXPECT_FALSE(out.exists("/surface/eta/3600"));
    EXPECT_EQ(out.readStringAttribute("/surface/eta/1800", "units"), "m");
    EXPECT_EQ(out.readStringAttribute("/surface/eta/1800", "long_name"),
              "free-surface elevation");
    EXPECT_DOUBLE_EQ(out.readDoubleAttribute("/surface/eta/1800", "time"), 1800.0);
    EXPECT_EQ(out.readStringAttribute("/frehg2", "config"), "config-text");

    // Layout contract: flattened index j*NX + i (plan §7).
    const std::vector<real_t> flat = out.readAll1D("/surface/eta/1800");
    ASSERT_EQ(flat.size(), 12u);
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 4; ++i) {
        const real_t stored = flat[static_cast<std::size_t>(j * 4 + i)];
        if (j == 1 && i == 2) {
          EXPECT_TRUE(std::isnan(stored)) << "inactive column must be NaN";
        } else {
          EXPECT_TRUE(bitEqual(stored, pattern2(j, i)))
              << "0-ULP mismatch at (" << j << ", " << i << ")";
        }
      }
    }
  }
}

TEST_F(Hdf5Test, Field3RoundTripLayoutAndKtopMask) {
  frehg::Grid grid(MPI_COMM_SELF, testDomain());
  frehg::HostField2<int> ktop("ktop", 3, 4);
  Kokkos::deep_copy(ktop, 0);
  ktop(0, 1) = 1;  // top layer inactive at (j=0, i=1)
  grid.buildGlobalIds(ktop);

  frehg::Field3<real_t> head("head", 5, 6, 2);
  auto headH = Kokkos::create_mirror_view(head);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 4; ++i) {
      for (int k = 0; k < 2; ++k) {
        headH(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1,
              static_cast<std::size_t>(k)) = pattern3(j, i, k);
      }
    }
  }
  Kokkos::deep_copy(head, headH);

  frehg::io::Hdf5Output out(grid, file("field3.h5"), "config");
  out.writeField3("groundwater", "hydraulic_head", 600.0, head, {"m", "hydraulic head"});

  // Layout contract: flattened index (j*NX + i)*NZ + k (plan §7, pinned by
  // b3-kirkland/makeplot.py).
  const std::vector<real_t> flat = out.readAll1D("/groundwater/hydraulic_head/600");
  ASSERT_EQ(flat.size(), 24u);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 4; ++i) {
      for (int k = 0; k < 2; ++k) {
        const real_t stored = flat[static_cast<std::size_t>((j * 4 + i) * 2 + k)];
        if (j == 0 && i == 1 && k == 0) {
          EXPECT_TRUE(std::isnan(stored)) << "cell above ktop must be NaN";
        } else {
          EXPECT_TRUE(bitEqual(stored, pattern3(j, i, k)))
              << "0-ULP mismatch at (" << j << ", " << i << ", " << k << ")";
        }
      }
    }
  }
}

TEST_F(Hdf5Test, GridMetaDatasets) {
  frehg::Grid grid(MPI_COMM_SELF, testDomain());
  frehg::io::Hdf5Output out(grid, file("gridmeta.h5"), "config");

  frehg::HostField2<real_t> bottom("bottom", 3, 4);
  for (std::size_t j = 0; j < 3; ++j) {
    for (std::size_t i = 0; i < 4; ++i) {
      bottom(j, i) = -1.0 - static_cast<real_t>(j) * 0.1;
    }
  }
  out.writeGridMeta(bottom);

  const auto xc = out.readAll1D("/grid/x_center");
  ASSERT_EQ(xc.size(), 4u);
  EXPECT_DOUBLE_EQ(xc[0], 0.25);
  EXPECT_DOUBLE_EQ(xc[3], 1.75);
  const auto zc = out.readAll1D("/grid/z_center");
  ASSERT_EQ(zc.size(), 2u);
  EXPECT_DOUBLE_EQ(zc[0], -0.05);
  const auto dz = out.readAll1D("/grid/dz");
  EXPECT_DOUBLE_EQ(dz[1], 0.1);
  // b3 script compatibility dataset: per-cell centers in the (j*NX+i)*NZ+k
  // layout (the P2 correction to the plan §7 contract — b3's makeplot.py
  // reads NX*NY*NZ values). The uniform-layer overload tiles the z axis.
  const auto zcell = out.readAll1D("/groundwater/zcell/0");
  ASSERT_EQ(zcell.size(), 24u);
  EXPECT_DOUBLE_EQ(zcell[0], -0.05);
  EXPECT_DOUBLE_EQ(zcell[1], -0.15);
  EXPECT_DOUBLE_EQ(zcell[2], -0.05);
  const auto bottomBack = out.readAll1D("/grid/bottom");
  ASSERT_EQ(bottomBack.size(), 12u);
  EXPECT_DOUBLE_EQ(bottomBack[4], -1.1);
  const auto ktopBack = out.readAll1D("/grid/ktop");
  EXPECT_DOUBLE_EQ(ktopBack[0], 0.0);
}

TEST_F(Hdf5Test, MonitorAppendsRows) {
  frehg::Grid grid(MPI_COMM_SELF, testDomain());
  frehg::io::Hdf5Output out(grid, file("monitor.h5"), "config");

  frehg::MonitorConfig cfg;
  cfg.name = "outlet_q";
  cfg.i = 0;
  cfg.j = 1;
  cfg.variables = {"depth", "vv"};
  frehg::io::Monitor monitor(out, cfg);
  EXPECT_TRUE(monitor.ownsPoint());

  monitor.record(0.0, {0.1, -0.2});
  monitor.record(5.0, {0.3, -0.4});
  monitor.flush();
  monitor.record(10.0, {0.5, -0.6});
  monitor.flush();
  monitor.flush();  // nothing buffered: no-op
  EXPECT_EQ(monitor.rowsWritten(), 3);

  // Read the 2D table back directly through the HDF5 C API.
  const hid_t dset = H5Dopen2(out.fileId(), "/monitor/outlet_q", H5P_DEFAULT);
  ASSERT_GE(dset, 0);
  const hid_t space = H5Dget_space(dset);
  hsize_t dims[2] = {0, 0};
  H5Sget_simple_extent_dims(space, dims, nullptr);
  EXPECT_EQ(dims[0], 3u);
  EXPECT_EQ(dims[1], 3u);
  std::vector<double> table(9);
  ASSERT_GE(H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, table.data()), 0);
  H5Sclose(space);
  H5Dclose(dset);
  EXPECT_DOUBLE_EQ(table[0], 0.0);
  EXPECT_DOUBLE_EQ(table[1], 0.1);
  EXPECT_DOUBLE_EQ(table[2], -0.2);
  EXPECT_DOUBLE_EQ(table[6], 10.0);
  EXPECT_DOUBLE_EQ(table[8], -0.6);
  EXPECT_EQ(out.readStringAttribute("/monitor/outlet_q", "columns"), "time,depth,vv");

  // Misuse: wrong value count and out-of-domain point are fatal.
  EXPECT_THROW(monitor.record(1.0, {0.1}), frehg::FatalError);
  frehg::MonitorConfig bad = cfg;
  bad.name = "bad";
  bad.i = 100;
  EXPECT_THROW(frehg::io::Monitor(out, bad), frehg::FatalError);
}

TEST_F(Hdf5Test, CheckpointRoundTripZeroUlp) {
  frehg::Grid grid(MPI_COMM_SELF, testDomain());
  frehg::io::Hdf5Output out(grid, file("checkpoint.h5"), "config");
  frehg::io::Checkpoint checkpoint(out);

  frehg::Field2<real_t> eta("eta", 5, 6);
  frehg::Field3<real_t> theta("theta", 5, 6, 2);
  auto etaH = Kokkos::create_mirror_view(eta);
  auto thetaH = Kokkos::create_mirror_view(theta);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 4; ++i) {
      etaH(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1) =
          pattern2(j, i) * 1.0e-17 + 3.0e8;  // stress the mantissa
      for (int k = 0; k < 2; ++k) {
        thetaH(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1,
               static_cast<std::size_t>(k)) = pattern3(j, i, k);
      }
    }
  }
  Kokkos::deep_copy(eta, etaH);
  Kokkos::deep_copy(theta, thetaH);

  checkpoint.write(3600.0, 720, {{"dtg", 0.03125}, {"cumulative_rain", 1.0 / 3.0}},
                   {{"eta", eta}}, {{"theta", theta}});
  EXPECT_TRUE(checkpoint.has(3600.0));
  EXPECT_FALSE(checkpoint.has(7200.0));

  // Wipe and restore.
  Kokkos::deep_copy(eta, 0.0);
  Kokkos::deep_copy(theta, 0.0);
  const auto header = checkpoint.read(3600.0, {{"eta", eta}}, {{"theta", theta}});
  EXPECT_DOUBLE_EQ(header.t, 3600.0);
  EXPECT_EQ(header.step, 720);
  ASSERT_EQ(header.scalars.size(), 2u);
  EXPECT_TRUE(bitEqual(header.scalars.at("dtg"), 0.03125));
  EXPECT_TRUE(bitEqual(header.scalars.at("cumulative_rain"), 1.0 / 3.0));

  auto etaBack = Kokkos::create_mirror_view(eta);
  auto thetaBack = Kokkos::create_mirror_view(theta);
  Kokkos::deep_copy(etaBack, eta);
  Kokkos::deep_copy(thetaBack, theta);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 4; ++i) {
      EXPECT_TRUE(bitEqual(
          etaBack(static_cast<std::size_t>(j) + 1, static_cast<std::size_t>(i) + 1),
          pattern2(j, i) * 1.0e-17 + 3.0e8))
          << "checkpoint eta not bit-exact at (" << j << ", " << i << ")";
      for (int k = 0; k < 2; ++k) {
        EXPECT_TRUE(bitEqual(thetaBack(static_cast<std::size_t>(j) + 1,
                                       static_cast<std::size_t>(i) + 1,
                                       static_cast<std::size_t>(k)),
                             pattern3(j, i, k)));
      }
    }
  }

  EXPECT_THROW(checkpoint.read(7200.0, {}, {}), frehg::FatalError);
}

}  // namespace
