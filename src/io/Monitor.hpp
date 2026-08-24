/// \file Monitor.hpp
/// \brief Point-probe time series appended to /monitor/\<name\> (plan §7).
///
/// A monitor samples named variables at one global cell every time it is
/// recorded and appends rows (time, value...) to an extendable HDF5 table.
/// Rows are buffered on the rank that owns the cell and written on flush();
/// flush is collective (dataset extension is a metadata operation), while
/// the row data itself is written independently by the owner.

#ifndef FREHG_IO_MONITOR_HPP
#define FREHG_IO_MONITOR_HPP

#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "io/Hdf5Output.hpp"

#include <string>
#include <vector>

namespace frehg::io {

/// One point monitor.
class Monitor {
 public:
  /// Create the extendable dataset /monitor/\<name\> (collective).
  /// The monitored cell (config.i, config.j) uses global indices; the
  /// constructor determines the owning rank from the grid decomposition.
  Monitor(Hdf5Output& out, const MonitorConfig& config);

  /// \return true when this rank owns the monitored cell.
  bool ownsPoint() const { return owns_; }

  /// Buffer one row. \p values must match the configured variables in
  /// number and order. Only the owning rank stores the row; other ranks
  /// return immediately, so record() may be called unconditionally.
  void record(real_t t, const std::vector<real_t>& values);

  /// Append all buffered rows to the file. Collective on the grid
  /// communicator.
  void flush();

  /// \return the monitor configuration.
  const MonitorConfig& config() const { return config_; }

  /// \return rows written to the file so far (excludes buffered rows).
  long rowsWritten() const { return rowsWritten_; }

 private:
  Hdf5Output& out_;
  MonitorConfig config_;
  bool owns_ = false;
  std::size_t columns_ = 0;
  std::vector<real_t> buffer_;  ///< row-major (time, values...) rows
  long rowsWritten_ = 0;
};

}  // namespace frehg::io

#endif  // FREHG_IO_MONITOR_HPP
