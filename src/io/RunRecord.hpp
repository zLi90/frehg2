/// \file RunRecord.hpp
/// \brief Persistent run provenance and timing record (v2 plan §2A).
///
/// Every simulation writes `run-record.yaml` into the directory of the HDF5
/// output: provenance (version, hosts, ranks/threads, wall clock, input
/// hash), the fully resolved configuration, module/BC summaries, the merged
/// hierarchical timer tree, per-system solver telemetry, and the closing
/// mass-audit budgets. The record is rewritten at every output flush and
/// finalized at exit, so an aborted run still leaves a truthful partial
/// record (provenance.finished = false).

#ifndef FREHG_IO_RUNRECORD_HPP
#define FREHG_IO_RUNRECORD_HPP

#include "core/Config.hpp"

#include <mpi.h>

#include <chrono>
#include <string>

namespace YAML {
class Node;
}

namespace frehg {
class BoundarySet;
}

namespace frehg::io {

/// Builder and writer of the run record. Rank 0 writes; flush() is
/// collective (it merges the timer tree across ranks).
class RunRecord {
 public:
  /// Capture the static sections: provenance (start time, launch geometry,
  /// input path + SHA-256), the resolved configuration, modules, and the
  /// boundary-condition summary. \p px, \p py is the realized rank
  /// decomposition; \p boundaries supplies the rasterized global member
  /// counts.
  RunRecord(MPI_Comm comm, const FrehgConfig& config, const std::string& inputPath,
            const BoundarySet& boundaries, int px, int py);
  ~RunRecord();

  /// Replace the per-system solver telemetry section (rank-0 values;
  /// callers reduce times across ranks first — the driver's summary
  /// reduction). \p yamlText is a YAML map, one key per system.
  void setSolver(const std::string& yamlText);

  /// Replace the closure section (final mass-audit budgets; rank-0 values).
  void setClosure(const std::string& yamlText);

  /// Write the record. Collective on the construction communicator (merges
  /// the timer tree). \p finished marks a completed run; mid-run flushes
  /// pass false so a killed run's record stays truthful.
  void flush(bool finished);

 private:
  MPI_Comm comm_;
  int rank_ = 0;
  std::string path_;        ///< run-record.yaml beside the HDF5 output
  std::string staticText_;  ///< provenance-through-BC YAML (constant part)
  std::string solverText_;  ///< solver section body (YAML map or empty)
  std::string closureText_; ///< closure section body (YAML map or empty)
  std::chrono::steady_clock::time_point began_;
};

}  // namespace frehg::io

#endif  // FREHG_IO_RUNRECORD_HPP
