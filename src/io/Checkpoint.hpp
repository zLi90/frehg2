/// \file Checkpoint.hpp
/// \brief Checkpoint/restart state under /checkpoint/\<t\>/ (plan §7).
///
/// A checkpoint stores the full prognostic state exactly (no NaN masking,
/// values byte-identical on read-back) together with the simulation time,
/// step counter, and named restart scalars such as the adaptive subsurface
/// step dtg. Restart determinism (run - checkpoint - restart equals an
/// uninterrupted run to 1e-12) is gated in P3; this component provides the
/// exact-round-trip storage layer it builds on.

#ifndef FREHG_IO_CHECKPOINT_HPP
#define FREHG_IO_CHECKPOINT_HPP

#include "core/Grid.hpp"
#include "core/Types.hpp"
#include "io/Hdf5Output.hpp"

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace frehg::io {

/// Writer/reader for checkpoint groups inside an Hdf5Output file.
class Checkpoint {
 public:
  /// \param out the output file to write checkpoints into / read from.
  explicit Checkpoint(Hdf5Output& out) : out_(out) {}

  /// Named 2D prognostic fields (with halos; interiors are stored).
  using Fields2 = std::vector<std::pair<std::string, Field2<real_t>>>;
  /// Named 3D prognostic fields.
  using Fields3 = std::vector<std::pair<std::string, Field3<real_t>>>;
  /// Named restart scalars (e.g. {"dtg", ...}).
  using Scalars = std::map<std::string, real_t>;

  /// Collectively write /checkpoint/\<t\>/ with attrs t, step, and one
  /// attribute per scalar. \p labelTime, when non-negative, keys the group
  /// (it must land on a whole second per the §7 contract) while \p t stays
  /// the exact simulation time stored in the header — adaptive-dtg runs
  /// checkpoint at the first step past a whole-second boundary, so the key
  /// and the state time differ by less than one dtg.
  void write(real_t t, long step, const Scalars& scalars, const Fields2& fields2,
             const Fields3& fields3, real_t labelTime = -1.0);

  /// Header data returned by read().
  struct Header {
    real_t t = 0;     ///< checkpoint simulation time [s]
    long step = 0;    ///< checkpoint step counter
    Scalars scalars;  ///< named restart scalars (e.g. dtg)
  };

  /// Read /checkpoint/\<t\>/ back into the passed fields (interiors are
  /// overwritten; halos are the caller's to refresh via halo exchange and
  /// boundary conditions). Missing checkpoint groups or fields are fatal.
  Header read(real_t t, const Fields2& fields2, const Fields3& fields3);

  /// \return true when a checkpoint exists for time \p t.
  bool has(real_t t) const;

 private:
  Hdf5Output& out_;
};

}  // namespace frehg::io

#endif  // FREHG_IO_CHECKPOINT_HPP
