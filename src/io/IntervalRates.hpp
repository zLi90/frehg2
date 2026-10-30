/// \file IntervalRates.hpp
/// \brief Interval-mean rate outputs (Stage-3 plan §8 W4.1).
///
/// A registry of named per-cell depth accumulators [m]. Contributors add
/// each step's per-cell deposit into the accumulator (addInto); on each
/// output write the registry writes acc / (t - t_prev) — an interval-MEAN
/// rate in m/s whose integral over the interval recovers the accumulated
/// volume exactly — then resets. This is the mechanism behind
/// `seepage_rate` (W4.2) and the T5-T7 `atmosphere` rates; snapshot
/// observables like `seepage` (qss) remain per-step samples and must not
/// be integrated (docs/user-guide/output.md). The accumulators are
/// prognostic across restarts: checkpoint fields `acc_<name>` plus the
/// scalar `acc_t_prev`. Only requested variables are allocated. The
/// Field3 (1/s) leg of the plan's registry arrives with its first
/// consumer, T6's `root_uptake` (W6.6) — nothing in T4/T5 reads it, and
/// authored-unexercised code is this project's documented failure class
/// (V2-A9..A12; amendment S3-A7).

#ifndef FREHG_IO_INTERVALRATES_HPP
#define FREHG_IO_INTERVALRATES_HPP

#include "core/Grid.hpp"
#include "core/Types.hpp"
#include "io/Checkpoint.hpp"
#include "io/Hdf5Output.hpp"

#include <string>
#include <vector>

namespace frehg::io {

/// The interval-rate registry (one per Simulation).
class IntervalRates {
 public:
  /// \param grid the decomposition (field extents, cell geometry).
  /// \param t0 the start of the first accumulation interval.
  IntervalRates(const Grid& grid, real_t t0) : grid_(grid), tPrev_(t0) {}

  /// Register a 2-D rate output: HDF5 dataset /\p group/\p name, written
  /// in m/s with \p meta. The accumulator starts zeroed.
  void add(const std::string& group, const std::string& name,
           const VarMeta& meta);

  /// Stage-3 T6 (the S3-A7 deferred Field3 leg, first consumer
  /// root_uptake): register a 3-D rate output. Contributors add the
  /// DIMENSIONLESS per-cell increment (sink·dtg, [(1/s)·s]); the write
  /// emits acc/(t − t_prev) [1/s].
  void add3(const std::string& group, const std::string& name,
            const VarMeta& meta);
  /// Add one step's per-cell increment into 3-D accumulator \p name.
  void addInto3(const std::string& name, const Field3<real_t>& increment);

  /// \return true when nothing is registered (all plumbing is skipped).
  bool empty() const { return entries_.empty() && entries3_.empty(); }

  /// \return true when \p name is a registered accumulator.
  bool has(const std::string& name) const {
    for (const Entry& entry : entries_) {
      if (entry.name == name) {
        return true;
      }
    }
    return false;
  }

  /// Add one step's per-cell deposit \p depth [m] into accumulator
  /// \p name (interior cells; both fields are (ny+2, nx+2)).
  void addInto(const std::string& name, const Field2<real_t>& depth);

  /// At an output write: write acc / (t - t_prev) for every entry under
  /// dataset key \p label, zero the accumulators, and start the next
  /// interval at \p t (the TRUE state time; an adaptive loop's label can
  /// trail it). At t == t_prev (the initial write) the frames are the
  /// zeros the accumulators hold.
  void writeAndReset(Hdf5Output& out, real_t label, real_t t);

  /// Checkpoint plumbing: the accumulators as fields `acc_<name>` and the
  /// interval start as the scalar `acc_t_prev` (appended by the caller).
  Checkpoint::Fields2 checkpointFields() const;
  /// The 3-D accumulators (T6), same naming convention.
  Checkpoint::Fields3 checkpointFields3() const;
  real_t tPrev() const { return tPrev_; }      ///< current interval start [s]
  void restore(real_t tPrev) { tPrev_ = tPrev; }  ///< restart: set the interval start

 private:
  struct Entry {
    std::string group;      ///< HDF5 group ("surface", "atmosphere", ...)
    std::string name;       ///< dataset / accumulator name
    VarMeta meta;           ///< units ("m s-1") and description
    Field2<real_t> acc;     ///< per-cell accumulated depth [m]
    Field2<real_t> rate;    ///< write buffer, acc / (t - t_prev) [m/s]
  };

  struct Entry3 {
    std::string group;      ///< HDF5 group ("groundwater")
    std::string name;       ///< dataset / accumulator name
    VarMeta meta;           ///< units ("s-1") and description
    Field3<real_t> acc;     ///< per-cell accumulated increment [-]
    Field3<real_t> rate;    ///< write buffer, acc / (t - t_prev) [1/s]
  };

  const Grid& grid_;
  real_t tPrev_;
  std::vector<Entry> entries_;
  std::vector<Entry3> entries3_;
};

}  // namespace frehg::io

#endif  // FREHG_IO_INTERVALRATES_HPP
