/// \file TimeSeries.hpp
/// \brief Piecewise-linear time series with end clamping (plan §5.6, §8.1).
///
/// Series drive boundary conditions and atmospheric forcing. Evaluation is
/// stateless (binary search), so restart determinism requires no cursor
/// state: re-evaluating at the restart time reproduces the interrupted run.

#ifndef FREHG_CORE_TIMESERIES_HPP
#define FREHG_CORE_TIMESERIES_HPP

#include "core/Types.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace frehg {

/// Piecewise-linear function of time built from (t, v) samples.
class TimeSeries {
 public:
  /// Empty series; value() on it is fatal. Exists so members can be
  /// default-constructed before assignment.
  TimeSeries() = default;

  /// Build from parallel sample arrays.
  /// Times must be strictly increasing and at least one sample is required;
  /// violations are fatal (plan §8.1 non-monotonic-time rejection).
  TimeSeries(std::vector<real_t> times, std::vector<real_t> values);

  /// Load a series from a text file: one "time value" pair per
  /// whitespace-separated line, '#' starts a comment. Units are the caller's
  /// contract (seconds + SI value in all Frehg2 configs).
  static TimeSeries fromFile(const std::string& path);

  /// Evaluate at time \p t: linear interpolation between samples, clamped to
  /// the first/last value outside the sampled range.
  real_t value(real_t t) const;

  /// \return number of samples.
  std::size_t size() const { return times_.size(); }

  /// \return true if the series has no samples.
  bool empty() const { return times_.empty(); }

  /// \return first sample time (fatal on empty series).
  real_t tMin() const;

  /// \return last sample time (fatal on empty series).
  real_t tMax() const;

 private:
  std::vector<real_t> times_;
  std::vector<real_t> values_;
};

}  // namespace frehg

#endif  // FREHG_CORE_TIMESERIES_HPP
