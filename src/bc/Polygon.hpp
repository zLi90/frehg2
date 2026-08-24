/// \file Polygon.hpp
/// \brief Point-in-polygon test used to rasterize boundary regions
///        (plan §5.6, adopted from SERGHEI's polygon BC system).

#ifndef FREHG_BC_POLYGON_HPP
#define FREHG_BC_POLYGON_HPP

#include "core/Types.hpp"

#include <array>
#include <vector>

namespace frehg {

/// A simple (non-self-intersecting) polygon in the horizontal plane.
class Polygon {
 public:
  /// \param vertices at least three (x, y) points; the closing edge from the
  ///        last vertex back to the first is implicit. Fewer than three
  ///        vertices is fatal.
  explicit Polygon(std::vector<std::array<real_t, 2>> vertices);

  /// Ray-cast containment test. Points exactly on an edge or vertex count as
  /// inside — the fixed convention of plan §5.6, unit-tested.
  bool contains(real_t x, real_t y) const;

  /// \return axis-aligned bounding box {xmin, ymin, xmax, ymax}.
  const std::array<real_t, 4>& boundingBox() const { return bbox_; }

  /// \return the vertex list.
  const std::vector<std::array<real_t, 2>>& vertices() const { return vertices_; }

 private:
  std::vector<std::array<real_t, 2>> vertices_;
  std::array<real_t, 4> bbox_ = {0, 0, 0, 0};
};

}  // namespace frehg

#endif  // FREHG_BC_POLYGON_HPP
