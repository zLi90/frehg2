/// \file TimeSeries.cpp
/// \brief Implementation of the piecewise-linear time series.

#include "core/TimeSeries.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <utility>

namespace frehg {

TimeSeries::TimeSeries(std::vector<real_t> times, std::vector<real_t> values)
    : times_(std::move(times)), values_(std::move(values)) {
  if (times_.size() != values_.size()) {
    log::fatal(log::msg() << "TimeSeries: " << times_.size() << " times vs " << values_.size()
                          << " values");
  }
  if (times_.empty()) {
    log::fatal("TimeSeries: at least one (time, value) sample is required");
  }
  for (std::size_t n = 1; n < times_.size(); ++n) {
    if (!(times_[n] > times_[n - 1])) {
      log::fatal(log::msg() << "TimeSeries: times must be strictly increasing (t[" << n - 1
                            << "]=" << times_[n - 1] << ", t[" << n << "]=" << times_[n] << ")");
    }
  }
}

TimeSeries TimeSeries::fromFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    log::fatal(log::msg() << "TimeSeries: cannot open '" << path << "'");
  }
  std::vector<real_t> times;
  std::vector<real_t> values;
  std::string line;
  std::size_t lineNo = 0;
  while (std::getline(in, line)) {
    ++lineNo;
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) {
      line.erase(hash);
    }
    std::istringstream ls(line);
    real_t t = 0;
    real_t v = 0;
    if (!(ls >> t)) {
      continue;  // blank or comment-only line
    }
    if (!(ls >> v)) {
      log::fatal(log::msg() << "TimeSeries: '" << path << "' line " << lineNo
                            << ": expected 'time value', got '" << line << "'");
    }
    real_t extra = 0;
    if (ls >> extra) {
      log::fatal(log::msg() << "TimeSeries: '" << path << "' line " << lineNo
                            << ": more than two columns");
    }
    times.push_back(t);
    values.push_back(v);
  }
  if (times.empty()) {
    log::fatal(log::msg() << "TimeSeries: '" << path << "' contains no samples");
  }
  return TimeSeries(std::move(times), std::move(values));
}

real_t TimeSeries::value(real_t t) const {
  if (times_.empty()) {
    log::fatal("TimeSeries::value on an empty series");
  }
  if (t <= times_.front()) {
    return values_.front();
  }
  if (t >= times_.back()) {
    return values_.back();
  }
  const auto upper = std::upper_bound(times_.begin(), times_.end(), t);
  const std::size_t hi = static_cast<std::size_t>(upper - times_.begin());
  const std::size_t lo = hi - 1;
  const real_t w = (t - times_[lo]) / (times_[hi] - times_[lo]);
  return values_[lo] + w * (values_[hi] - values_[lo]);
}

real_t TimeSeries::tMin() const {
  if (times_.empty()) {
    log::fatal("TimeSeries::tMin on an empty series");
  }
  return times_.front();
}

real_t TimeSeries::tMax() const {
  if (times_.empty()) {
    log::fatal("TimeSeries::tMax on an empty series");
  }
  return times_.back();
}

}  // namespace frehg
