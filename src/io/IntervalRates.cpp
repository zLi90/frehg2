/// \file IntervalRates.cpp
/// \brief Interval-mean rate accumulation and write-out (plan §8 W4.1).

#include "io/IntervalRates.hpp"

#include "core/Logger.hpp"

#include <Kokkos_Core.hpp>

namespace frehg::io {

void IntervalRates::add(const std::string& group, const std::string& name,
                        const VarMeta& meta) {
  const std::size_t ny2 = static_cast<std::size_t>(grid_.nyLocal()) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(grid_.nxLocal()) + 2;
  Entry entry{group, name, meta, Field2<real_t>("rate_acc_" + name, ny2, nx2),
              Field2<real_t>("rate_out_" + name, ny2, nx2)};
  Kokkos::deep_copy(entry.acc, 0.0);
  Kokkos::deep_copy(entry.rate, 0.0);
  entries_.push_back(std::move(entry));
}

void IntervalRates::addInto(const std::string& name, const Field2<real_t>& depth) {
  for (Entry& entry : entries_) {
    if (entry.name != name) {
      continue;
    }
    const int nyl = grid_.nyLocal();
    const int nxl = grid_.nxLocal();
    Field2<real_t> acc = entry.acc;
    Field2<real_t> src = depth;
    Kokkos::parallel_for(
        "rate_accumulate",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i) { acc(j, i) += src(j, i); });
    return;
  }
  log::fatal(log::msg() << "interval rate '" << name << "' is not registered");
}

void IntervalRates::writeAndReset(Hdf5Output& out, real_t label, real_t t) {
  const real_t elapsed = t - tPrev_;
  const int nyl = grid_.nyLocal();
  const int nxl = grid_.nxLocal();
  for (Entry& entry : entries_) {
    Field2<real_t> acc = entry.acc;
    Field2<real_t> rate = entry.rate;
    // The initial write (and a zero-length interval after a forced write)
    // has nothing accumulated: the frame is zeros, not 0/0.
    const real_t inv = (elapsed > 0.0) ? 1.0 / elapsed : 0.0;
    Kokkos::parallel_for(
        "rate_write_reset",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, Kokkos::IndexType<int>>(
            {1, 1}, {nyl + 1, nxl + 1}),
        KOKKOS_LAMBDA(const int j, const int i) {
          rate(j, i) = acc(j, i) * inv;
          acc(j, i) = 0.0;
        });
    out.writeField2(entry.group, entry.name, label, entry.rate, entry.meta);
  }
  for (Entry3& entry : entries3_) {
    Field3<real_t> acc = entry.acc;
    Field3<real_t> rate = entry.rate;
    const real_t inv = (elapsed > 0.0) ? 1.0 / elapsed : 0.0;
    const int nz = grid_.nz();
    Kokkos::parallel_for(
        "rate3_write_reset",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
            {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
        KOKKOS_LAMBDA(const int j, const int i, const int k) {
          rate(j, i, k) = acc(j, i, k) * inv;
          acc(j, i, k) = 0.0;
        });
    out.writeField3(entry.group, entry.name, label, entry.rate, entry.meta);
  }
  tPrev_ = t;
}

void IntervalRates::add3(const std::string& group, const std::string& name,
                         const VarMeta& meta) {
  const std::size_t ny2 = static_cast<std::size_t>(grid_.nyLocal()) + 2;
  const std::size_t nx2 = static_cast<std::size_t>(grid_.nxLocal()) + 2;
  const std::size_t nz = static_cast<std::size_t>(grid_.nz());
  Entry3 entry{group, name, meta,
               Field3<real_t>("rate3_acc_" + name, ny2, nx2, nz),
               Field3<real_t>("rate3_out_" + name, ny2, nx2, nz)};
  Kokkos::deep_copy(entry.acc, 0.0);
  Kokkos::deep_copy(entry.rate, 0.0);
  entries3_.push_back(std::move(entry));
}

void IntervalRates::addInto3(const std::string& name, const Field3<real_t>& increment) {
  for (Entry3& entry : entries3_) {
    if (entry.name != name) {
      continue;
    }
    const int nyl = grid_.nyLocal();
    const int nxl = grid_.nxLocal();
    const int nz = grid_.nz();
    Field3<real_t> acc = entry.acc;
    Field3<real_t> src = increment;
    Kokkos::parallel_for(
        "rate3_accumulate",
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>, Kokkos::IndexType<int>>(
            {1, 1, 0}, {nyl + 1, nxl + 1, nz}),
        KOKKOS_LAMBDA(const int j, const int i, const int k) {
          acc(j, i, k) += src(j, i, k);
        });
    return;
  }
  log::fatal(log::msg() << "interval rate '" << name << "' is not registered (3-D)");
}

Checkpoint::Fields3 IntervalRates::checkpointFields3() const {
  Checkpoint::Fields3 fields;
  for (const Entry3& entry : entries3_) {
    fields.push_back({"acc_" + entry.name, entry.acc});
  }
  return fields;
}

Checkpoint::Fields2 IntervalRates::checkpointFields() const {
  Checkpoint::Fields2 fields;
  for (const Entry& entry : entries_) {
    fields.push_back({"acc_" + entry.name, entry.acc});
  }
  return fields;
}

}  // namespace frehg::io
