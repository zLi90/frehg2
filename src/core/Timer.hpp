/// \file Timer.hpp
/// \brief Hierarchical, MPI-reduced wall-clock timers (plan §3.3).
///
/// Sections are identified by name and nest automatically: starting a section
/// while another is active records it under the active section's path, e.g.
/// starting "assembly" inside "solve" accumulates into "solve/assembly". The
/// report merges timings across ranks (min/mean/max) so load imbalance is
/// visible per section.

#ifndef FREHG_CORE_TIMER_HPP
#define FREHG_CORE_TIMER_HPP

#include <mpi.h>

#include <string>
#include <vector>

namespace frehg {

/// One timer section merged across ranks (v2 plan §2A: the run record and
/// the text report share this so they cannot disagree).
struct TimerSection {
  std::string path;    ///< full nested path, e.g. "simulation/swe/velocity"
  long count = 0;      ///< completed cycles (max over ranks)
  double minSeconds = 0.0;
  double meanSeconds = 0.0;
  double maxSeconds = 0.0;
};

/// Static registry of named, nesting wall-clock timers.
class Timer {
 public:
  /// Start (or resume) the named section, nested under the active section.
  static void start(const std::string& name);

  /// Stop the named section. The name must match the innermost active
  /// section; a mismatch is a programming error and is fatal.
  static void stop(const std::string& name);

  /// \return accumulated wall seconds of a section by full path
  ///         (e.g. "solve/assembly"); 0 if the path was never timed.
  static double elapsed(const std::string& path);

  /// \return number of completed start/stop cycles for a path.
  static long count(const std::string& path);

  /// Merge all sections across \p comm: min/mean/max seconds per path.
  /// Collective on \p comm; the result is non-empty on rank 0 only.
  /// Sections still running are excluded.
  static std::vector<TimerSection> merged(MPI_Comm comm);

  /// Merge all sections across \p comm and return a formatted table
  /// (non-empty on rank 0 only). Collective on \p comm. Sections still
  /// running are excluded from the report.
  static std::string report(MPI_Comm comm);

  /// Clear all sections and the nesting stack (used by tests).
  static void reset();

  /// RAII helper timing the lifetime of a scope.
  class Scoped {
   public:
    /// Starts section \p name; the destructor stops it.
    explicit Scoped(std::string name);
    ~Scoped();
    Scoped(const Scoped&) = delete;
    Scoped& operator=(const Scoped&) = delete;

   private:
    std::string name_;
  };
};

}  // namespace frehg

#endif  // FREHG_CORE_TIMER_HPP
