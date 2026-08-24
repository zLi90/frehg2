/// \file Polygon.cpp
/// \brief Ray-cast point-in-polygon with on-edge = inside.

#include "bc/Polygon.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace frehg {

Polygon::Polygon(std::vector<std::array<real_t, 2>> vertices) : vertices_(std::move(vertices)) {
  if (vertices_.size() < 3) {
    log::fatal(log::msg() << "Polygon: needs at least 3 vertices, got " << vertices_.size());
  }
  bbox_ = {vertices_[0][0], vertices_[0][1], vertices_[0][0], vertices_[0][1]};
  for (const auto& v : vertices_) {
    bbox_[0] = std::min(bbox_[0], v[0]);
    bbox_[1] = std::min(bbox_[1], v[1]);
    bbox_[2] = std::max(bbox_[2], v[0]);
    bbox_[3] = std::max(bbox_[3], v[1]);
  }
}

namespace {

/// True when point p lies on the closed segment [a, b], within a tolerance
/// scaled to the segment size.
bool onSegment(const std::array<real_t, 2>& a, const std::array<real_t, 2>& b, real_t x,
               real_t y) {
  const real_t cross = (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]);
  const real_t scale = std::max({std::abs(a[0]), std::abs(a[1]), std::abs(b[0]), std::abs(b[1]),
                                 std::abs(x), std::abs(y), real_t(1)});
  if (std::abs(cross) > 1.0e-12 * scale * scale) {
    return false;
  }
  const real_t dot = (x - a[0]) * (b[0] - a[0]) + (y - a[1]) * (b[1] - a[1]);
  const real_t len2 = (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]);
  return dot >= -1.0e-12 * scale * scale && dot <= len2 + 1.0e-12 * scale * scale;
}

}  // namespace

bool Polygon::contains(real_t x, real_t y) const {
  if (x < bbox_[0] - 1.0e-12 || x > bbox_[2] + 1.0e-12 || y < bbox_[1] - 1.0e-12 ||
      y > bbox_[3] + 1.0e-12) {
    return false;
  }

  // On-edge (including vertices) counts as inside.
  const std::size_t n = vertices_.size();
  for (std::size_t v = 0; v < n; ++v) {
    if (onSegment(vertices_[v], vertices_[(v + 1) % n], x, y)) {
      return true;
    }
  }

  // Standard crossing-number ray cast with the half-open rule, which makes
  // vertex crossings count exactly once.
  bool inside = false;
  for (std::size_t v = 0; v < n; ++v) {
    const auto& a = vertices_[v];
    const auto& b = vertices_[(v + 1) % n];
    const bool crossesY = (a[1] > y) != (b[1] > y);
    if (crossesY) {
      const real_t xCross = a[0] + (y - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
      if (x < xCross) {
        inside = !inside;
      }
    }
  }
  return inside;
}

}  // namespace frehg
