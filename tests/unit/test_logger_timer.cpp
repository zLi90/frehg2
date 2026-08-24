/// \file test_logger_timer.cpp
/// \brief Logger fatal-path and hierarchical timer tests.

#include "core/Logger.hpp"
#include "core/Timer.hpp"
#include "core/Types.hpp"

#include <gtest/gtest.h>
#include <mpi.h>

#include <chrono>
#include <thread>

namespace {

TEST(Logger, FatalThrowsFatalError) {
  EXPECT_THROW(frehg::log::fatal("intentional test failure"), frehg::FatalError);
}

TEST(Logger, MsgBuilderComposes) {
  const std::string text = frehg::log::msg() << "nx=" << 42 << ", dx=" << 0.5;
  EXPECT_EQ(text, "nx=42, dx=0.5");
}

TEST(Timer, AccumulatesAndCounts) {
  frehg::Timer::reset();
  for (int rep = 0; rep < 3; ++rep) {
    frehg::Timer::start("outer");
    frehg::Timer::start("inner");
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
    frehg::Timer::stop("inner");
    frehg::Timer::stop("outer");
  }
  EXPECT_EQ(frehg::Timer::count("outer"), 3);
  EXPECT_EQ(frehg::Timer::count("outer/inner"), 3);
  EXPECT_GT(frehg::Timer::elapsed("outer/inner"), 0.0);
  EXPECT_GE(frehg::Timer::elapsed("outer"), frehg::Timer::elapsed("outer/inner"));
  EXPECT_EQ(frehg::Timer::count("unknown"), 0);
  EXPECT_EQ(frehg::Timer::elapsed("unknown"), 0.0);
}

TEST(Timer, NestingBuildsPaths) {
  frehg::Timer::reset();
  {
    frehg::Timer::Scoped solve("solve");
    frehg::Timer::Scoped assembly("assembly");
  }
  EXPECT_EQ(frehg::Timer::count("solve/assembly"), 1);
  EXPECT_EQ(frehg::Timer::count("solve"), 1);
  EXPECT_EQ(frehg::Timer::count("assembly"), 0);
}

TEST(Timer, MismatchedStopIsFatal) {
  frehg::Timer::reset();
  frehg::Timer::start("a");
  EXPECT_THROW(frehg::Timer::stop("b"), frehg::FatalError);
  frehg::Timer::stop("a");
  EXPECT_THROW(frehg::Timer::stop("a"), frehg::FatalError);  // stack now empty
}

TEST(Timer, ReportMergesAcrossComm) {
  frehg::Timer::reset();
  frehg::Timer::start("reported");
  std::this_thread::sleep_for(std::chrono::milliseconds(1));
  frehg::Timer::stop("reported");
  const std::string report = frehg::Timer::report(MPI_COMM_SELF);
  EXPECT_NE(report.find("reported"), std::string::npos);
  EXPECT_NE(report.find("timer report"), std::string::npos);
}

}  // namespace
