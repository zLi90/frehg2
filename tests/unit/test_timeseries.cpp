/// \file test_timeseries.cpp
/// \brief TimeSeries tests: interpolation, clamping, monotonicity rejection
///        (plan §8.1).

#include "core/TimeSeries.hpp"
#include "core/Types.hpp"

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

namespace {

using frehg::TimeSeries;

TEST(TimeSeries, InteriorLinearInterpolation) {
  const TimeSeries ts({0.0, 10.0, 20.0}, {1.0, 3.0, -1.0});
  EXPECT_DOUBLE_EQ(ts.value(0.0), 1.0);
  EXPECT_DOUBLE_EQ(ts.value(5.0), 2.0);
  EXPECT_DOUBLE_EQ(ts.value(10.0), 3.0);
  EXPECT_DOUBLE_EQ(ts.value(15.0), 1.0);
  EXPECT_DOUBLE_EQ(ts.value(20.0), -1.0);
}

TEST(TimeSeries, EndClamping) {
  const TimeSeries ts({100.0, 200.0}, {5.0, 7.0});
  EXPECT_DOUBLE_EQ(ts.value(-1000.0), 5.0);
  EXPECT_DOUBLE_EQ(ts.value(99.0), 5.0);
  EXPECT_DOUBLE_EQ(ts.value(201.0), 7.0);
  EXPECT_DOUBLE_EQ(ts.value(1.0e9), 7.0);
  EXPECT_DOUBLE_EQ(ts.tMin(), 100.0);
  EXPECT_DOUBLE_EQ(ts.tMax(), 200.0);
}

TEST(TimeSeries, SinglePointIsConstant) {
  const TimeSeries ts({50.0}, {3.5});
  EXPECT_DOUBLE_EQ(ts.value(0.0), 3.5);
  EXPECT_DOUBLE_EQ(ts.value(50.0), 3.5);
  EXPECT_DOUBLE_EQ(ts.value(100.0), 3.5);
}

TEST(TimeSeries, NonMonotonicTimesRejected) {
  EXPECT_THROW(TimeSeries({0.0, 10.0, 10.0}, {1.0, 2.0, 3.0}), frehg::FatalError);
  EXPECT_THROW(TimeSeries({0.0, 10.0, 5.0}, {1.0, 2.0, 3.0}), frehg::FatalError);
}

TEST(TimeSeries, EmptyAndMismatchedRejected) {
  EXPECT_THROW(TimeSeries({}, {}), frehg::FatalError);
  EXPECT_THROW(TimeSeries({0.0, 1.0}, {1.0}), frehg::FatalError);
  EXPECT_THROW(TimeSeries().value(0.0), frehg::FatalError);
}

class TimeSeriesFileTest : public ::testing::Test {
 protected:
  std::string path_ = ::testing::TempDir() + "frehg_test_series.dat";
  void write(const std::string& content) {
    std::ofstream out(path_);
    out << content;
  }
  void TearDown() override { std::remove(path_.c_str()); }
};

TEST_F(TimeSeriesFileTest, ReadsTwoColumnFile) {
  write("# rain series\n0.0 1.0\n12000.0 1.0\n\n12001.0  0.0 # ramp off\n");
  const TimeSeries ts = TimeSeries::fromFile(path_);
  EXPECT_EQ(ts.size(), 3u);
  EXPECT_DOUBLE_EQ(ts.value(6000.0), 1.0);
  EXPECT_DOUBLE_EQ(ts.value(12000.5), 0.5);
  EXPECT_DOUBLE_EQ(ts.value(20000.0), 0.0);
}

TEST_F(TimeSeriesFileTest, MissingFileIsFatal) {
  EXPECT_THROW(TimeSeries::fromFile(path_ + ".does_not_exist"), frehg::FatalError);
}

TEST_F(TimeSeriesFileTest, OneColumnLineIsFatal) {
  write("0.0 1.0\n5.0\n");
  EXPECT_THROW(TimeSeries::fromFile(path_), frehg::FatalError);
}

TEST_F(TimeSeriesFileTest, ThreeColumnLineIsFatal) {
  write("0.0 1.0 2.0\n");
  EXPECT_THROW(TimeSeries::fromFile(path_), frehg::FatalError);
}

TEST_F(TimeSeriesFileTest, EmptyFileIsFatal) {
  write("# only comments\n\n");
  EXPECT_THROW(TimeSeries::fromFile(path_), frehg::FatalError);
}

TEST_F(TimeSeriesFileTest, NonMonotonicFileIsFatal) {
  write("0.0 1.0\n10.0 2.0\n10.0 3.0\n");
  EXPECT_THROW(TimeSeries::fromFile(path_), frehg::FatalError);
}

}  // namespace
