/// \file test_griddatareader.cpp
/// \brief Gridded-input reader tests: flat lists, headered rasters, nodata,
///        id maps, 3D fields, and malformed-input rejection.

#include "core/Types.hpp"
#include "io/GridDataReader.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

class GridDataReaderTest : public ::testing::Test {
 protected:
  std::filesystem::path dir_;

  void SetUp() override {
    dir_ = std::filesystem::path(::testing::TempDir()) / "frehg_griddata_test";
    std::filesystem::create_directories(dir_);
  }
  void TearDown() override { std::filesystem::remove_all(dir_); }

  std::string write(const std::string& name, const std::string& content) {
    const std::filesystem::path path = dir_ / name;
    std::ofstream out(path);
    out << content;
    return path.string();
  }
};

TEST_F(GridDataReaderTest, FlatList2D) {
  // One value per line (the legacy format), order j*nx + i.
  const std::string path = write("flat.dat", "1.5\n2.5\n3.5\n4.5\n5.5\n6.5\n");
  const auto values = frehg::io::readRaster2D(path, 3, 2);
  ASSERT_EQ(values.size(), 6u);
  EXPECT_DOUBLE_EQ(values[0 * 3 + 0], 1.5);
  EXPECT_DOUBLE_EQ(values[0 * 3 + 2], 3.5);
  EXPECT_DOUBLE_EQ(values[1 * 3 + 0], 4.5);
  EXPECT_DOUBLE_EQ(values[1 * 3 + 2], 6.5);
}

TEST_F(GridDataReaderTest, HeaderedRasterWithNodata) {
  const std::string path = write("raster.asc",
                                 "ncols 3\nnrows 2\nxllcorner 0.0\nyllcorner 0.0\n"
                                 "cellsize 1.0\nnodata_value -9999\n"
                                 "1 2 -9999\n4 5 6\n");
  const auto values = frehg::io::readRaster2D(path, 3, 2);
  ASSERT_EQ(values.size(), 6u);
  // First data row is j = 0 (file order, plan-fidelity convention).
  EXPECT_DOUBLE_EQ(values[0], 1.0);
  EXPECT_TRUE(std::isnan(values[2]));
  EXPECT_DOUBLE_EQ(values[5], 6.0);
}

TEST_F(GridDataReaderTest, HeaderExtentMismatchIsFatal) {
  const std::string path = write("wrong.asc", "ncols 3\nnrows 2\n1 2 3\n4 5 6\n");
  EXPECT_THROW(frehg::io::readRaster2D(path, 4, 2), frehg::FatalError);
}

TEST_F(GridDataReaderTest, ValueCountMismatchIsFatal) {
  const std::string path = write("short.dat", "1 2 3 4 5\n");
  EXPECT_THROW(frehg::io::readRaster2D(path, 3, 2), frehg::FatalError);
}

TEST_F(GridDataReaderTest, BadNumberIsFatal) {
  const std::string path = write("bad.dat", "1 2\nx 4\n5 6\n");
  EXPECT_THROW(frehg::io::readRaster2D(path, 2, 3), frehg::FatalError);
}

TEST_F(GridDataReaderTest, UnknownHeaderKeyIsFatal) {
  const std::string path = write("odd.asc", "n_soil 2\n1 2 3\n4 5 6\n");
  EXPECT_THROW(frehg::io::readRaster2D(path, 3, 2), frehg::FatalError);
}

TEST_F(GridDataReaderTest, MissingFileIsFatal) {
  EXPECT_THROW(frehg::io::readRaster2D((dir_ / "absent.dat").string(), 2, 2),
               frehg::FatalError);
}

TEST_F(GridDataReaderTest, IntRaster) {
  const std::string path = write("ids.dat", "0 1 1\n0 0 1\n");
  const auto ids = frehg::io::readIntRaster2D(path, 3, 2);
  ASSERT_EQ(ids.size(), 6u);
  EXPECT_EQ(ids[1], 1);
  EXPECT_EQ(ids[3], 0);

  const std::string frac = write("frac.dat", "0 1 1\n0 0.5 1\n");
  EXPECT_THROW(frehg::io::readIntRaster2D(frac, 3, 2), frehg::FatalError);
}

TEST_F(GridDataReaderTest, Field3DFlatOrder) {
  // 2 x 1 x 3: values (j*nx + i)*nz + k with comments interleaved.
  const std::string path = write("field3.dat", "# 3D field\n1 2 3\n4 5 6 # cell (0,1)\n");
  const auto values = frehg::io::readField3D(path, 2, 1, 3);
  ASSERT_EQ(values.size(), 6u);
  EXPECT_DOUBLE_EQ(values[(0 * 2 + 0) * 3 + 0], 1.0);
  EXPECT_DOUBLE_EQ(values[(0 * 2 + 0) * 3 + 2], 3.0);
  EXPECT_DOUBLE_EQ(values[(0 * 2 + 1) * 3 + 1], 5.0);
}

TEST_F(GridDataReaderTest, HeaderedFile3DIsFatal) {
  const std::string path = write("hdr3.dat", "ncols 2\nnrows 1\n1 2 3 4 5 6\n");
  EXPECT_THROW(frehg::io::readField3D(path, 2, 1, 3), frehg::FatalError);
}

}  // namespace
