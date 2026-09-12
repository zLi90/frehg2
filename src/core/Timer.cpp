/// \file Timer.cpp
/// \brief Implementation of the hierarchical MPI-reduced timers.

#include "core/Timer.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <map>
#include <sstream>
#include <utility>
#include <vector>

namespace frehg {

namespace {

using Clock = std::chrono::steady_clock;

struct Section {
  double seconds = 0.0;
  long cycles = 0;
};

struct ActiveFrame {
  std::string name;   ///< bare section name as passed to start()
  std::string path;   ///< full nested path
  Clock::time_point began;
};

struct TimerState {
  std::map<std::string, Section> sections;  // keyed by full path
  std::vector<ActiveFrame> stack;
};

TimerState& state() {
  static TimerState s;
  return s;
}

}  // namespace

void Timer::start(const std::string& name) {
  TimerState& s = state();
  std::string path = s.stack.empty() ? name : s.stack.back().path + "/" + name;
  s.stack.push_back(ActiveFrame{name, std::move(path), Clock::now()});
}

void Timer::stop(const std::string& name) {
  TimerState& s = state();
  if (s.stack.empty()) {
    log::fatal(log::msg() << "Timer::stop(\"" << name << "\") with no active section");
  }
  const ActiveFrame& top = s.stack.back();
  if (top.name != name) {
    log::fatal(log::msg() << "Timer::stop(\"" << name << "\") does not match the active section \""
                          << top.name << "\" (sections must nest)");
  }
  const double dt = std::chrono::duration<double>(Clock::now() - top.began).count();
  Section& sec = s.sections[top.path];
  sec.seconds += dt;
  sec.cycles += 1;
  s.stack.pop_back();
}

double Timer::elapsed(const std::string& path) {
  const TimerState& s = state();
  auto it = s.sections.find(path);
  return it == s.sections.end() ? 0.0 : it->second.seconds;
}

long Timer::count(const std::string& path) {
  const TimerState& s = state();
  auto it = s.sections.find(path);
  return it == s.sections.end() ? 0L : it->second.cycles;
}

void Timer::reset() {
  state().sections.clear();
  state().stack.clear();
}

std::vector<TimerSection> Timer::merged(MPI_Comm comm) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // Serialize this rank's sections as "path\tseconds\tcycles\n" lines.
  std::ostringstream ser;
  ser << std::setprecision(17);
  for (const auto& [path, sec] : state().sections) {
    ser << path << "\t" << sec.seconds << "\t" << sec.cycles << "\n";
  }
  const std::string mine = ser.str();
  const int myLen = static_cast<int>(mine.size());

  std::vector<int> lengths(static_cast<std::size_t>(size), 0);
  MPI_Gather(&myLen, 1, MPI_INT, lengths.data(), 1, MPI_INT, 0, comm);

  std::vector<int> offsets(static_cast<std::size_t>(size), 0);
  int total = 0;
  if (rank == 0) {
    for (std::size_t r = 0; r < lengths.size(); ++r) {
      offsets[r] = total;
      total += lengths[r];
    }
  }
  std::vector<char> gathered(static_cast<std::size_t>(std::max(total, 1)));
  MPI_Gatherv(mine.data(), myLen, MPI_CHAR, gathered.data(), lengths.data(), offsets.data(),
              MPI_CHAR, 0, comm);

  if (rank != 0) {
    return {};
  }

  struct Merged {
    double minSec = 0.0;
    double maxSec = 0.0;
    double sumSec = 0.0;
    long cycles = 0;
    int nRanks = 0;
  };
  std::map<std::string, Merged> merged;
  std::istringstream all(std::string(gathered.data(), static_cast<std::size_t>(total)));
  std::string line;
  while (std::getline(all, line)) {
    const std::size_t tab1 = line.find('\t');
    const std::size_t tab2 = line.find('\t', tab1 + 1);
    if (tab1 == std::string::npos || tab2 == std::string::npos) {
      continue;
    }
    const std::string path = line.substr(0, tab1);
    const double sec = std::stod(line.substr(tab1 + 1, tab2 - tab1 - 1));
    const long cyc = std::stol(line.substr(tab2 + 1));
    Merged& m = merged[path];
    if (m.nRanks == 0) {
      m.minSec = sec;
      m.maxSec = sec;
    } else {
      m.minSec = std::min(m.minSec, sec);
      m.maxSec = std::max(m.maxSec, sec);
    }
    m.sumSec += sec;
    m.cycles = std::max(m.cycles, cyc);
    m.nRanks += 1;
  }

  std::vector<TimerSection> out;
  out.reserve(merged.size());
  for (const auto& [path, m] : merged) {
    TimerSection section;
    section.path = path;
    section.count = m.cycles;
    section.minSeconds = m.minSec;
    section.meanSeconds = m.sumSec / m.nRanks;
    section.maxSeconds = m.maxSec;
    out.push_back(std::move(section));
  }
  return out;
}

std::string Timer::report(MPI_Comm comm) {
  int size = 1;
  MPI_Comm_size(comm, &size);
  const std::vector<TimerSection> sections = merged(comm);
  if (sections.empty()) {
    return std::string();  // non-root ranks (and the no-sections edge case)
  }

  std::ostringstream out;
  out << "timer report (" << size << " rank" << (size > 1 ? "s" : "") << ")\n";
  out << "  " << std::left << std::setw(40) << "section" << std::right << std::setw(10) << "count"
      << std::setw(12) << "min[s]" << std::setw(12) << "mean[s]" << std::setw(12) << "max[s]"
      << "\n";
  out << std::fixed << std::setprecision(4);
  for (const TimerSection& s : sections) {
    out << "  " << std::left << std::setw(40) << s.path << std::right << std::setw(10) << s.count
        << std::setw(12) << s.minSeconds << std::setw(12) << s.meanSeconds << std::setw(12)
        << s.maxSeconds << "\n";
  }
  return out.str();
}

Timer::Scoped::Scoped(std::string name) : name_(std::move(name)) { Timer::start(name_); }

Timer::Scoped::~Scoped() { Timer::stop(name_); }

}  // namespace frehg
