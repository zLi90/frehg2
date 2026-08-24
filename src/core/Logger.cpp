/// \file Logger.cpp
/// \brief Implementation of the rank-aware logger.

#include "core/Logger.hpp"

#include "core/Types.hpp"

#include <iostream>

namespace frehg::log {

namespace {

struct LoggerState {
  int rank = 0;
  bool haveComm = false;
  Level minLevel = Level::Info;
};

LoggerState& state() {
  static LoggerState s;
  return s;
}

const char* levelTag(Level lvl) {
  switch (lvl) {
    case Level::Debug: return "DEBUG";
    case Level::Info:  return "INFO ";
    case Level::Warn:  return "WARN ";
    case Level::Error: return "ERROR";
  }
  return "?    ";
}

void emit(Level lvl, const std::string& message, bool rankZeroOnly) {
  LoggerState& s = state();
  if (lvl < s.minLevel) {
    return;
  }
  if (rankZeroOnly && s.rank != 0) {
    return;
  }
  std::ostream& out = (lvl >= Level::Warn) ? std::cerr : std::cout;
  out << "[frehg:" << s.rank << "][" << levelTag(lvl) << "] " << message << "\n";
  out.flush();
}

}  // namespace

void init(MPI_Comm comm, Level minLevel) {
  LoggerState& s = state();
  MPI_Comm_rank(comm, &s.rank);
  s.haveComm = true;
  s.minLevel = minLevel;
}

void setLevel(Level minLevel) { state().minLevel = minLevel; }

Level level() { return state().minLevel; }

void debug(const std::string& message) { emit(Level::Debug, message, false); }

void info(const std::string& message) { emit(Level::Info, message, true); }

void warn(const std::string& message) { emit(Level::Warn, message, false); }

void error(const std::string& message) { emit(Level::Error, message, false); }

void fatal(const std::string& message) {
  emit(Level::Error, "FATAL: " + message, false);
  throw FatalError(message);
}

}  // namespace frehg::log
