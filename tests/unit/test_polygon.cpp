/// \file test_polygon.cpp
/// \brief Point-in-polygon tests: convex, concave, on-edge convention, and
///        BoundarySet rasterization on a single rank (plan §8.1).

#include "bc/BoundarySet.hpp"
#include "bc/Polygon.hpp"
#include "core/Types.hpp"

#include <gtest/gtest.h>
#include <mpi.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace {

using frehg::Polygon;

TEST(Polygon, ConvexSquare) {
  const Polygon square({{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}});
  EXPECT_TRUE(square.contains(1.0, 1.0));
  EXPECT_TRUE(square.contains(0.1, 1.9));
  EXPECT_FALSE(square.contains(2.1, 1.0));
  EXPECT_FALSE(square.contains(-0.1, 1.0));
  EXPECT_FALSE(square.contains(1.0, -5.0));
}

TEST(Polygon, OnEdgeAndVertexAreInside) {
  const Polygon square({{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}});
  EXPECT_TRUE(square.contains(1.0, 0.0));  // mid bottom edge
  EXPECT_TRUE(square.contains(2.0, 1.0));  // mid right edge
  EXPECT_TRUE(square.contains(0.0, 0.0));  // vertex
  EXPECT_TRUE(square.contains(2.0, 2.0));  // vertex
}

TEST(Polygon, ConcaveLShape) {
  // L-shape: unit notch removed from the top-right of a 2 x 2 square.
  const Polygon lShape(
      {{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}, {0.0, 2.0}});
  EXPECT_TRUE(lShape.contains(0.5, 0.5));
  EXPECT_TRUE(lShape.contains(1.5, 0.5));
  EXPECT_TRUE(lShape.contains(0.5, 1.5));
  EXPECT_FALSE(lShape.contains(1.5, 1.5));  // inside the notch
  EXPECT_TRUE(lShape.contains(1.0, 1.5));   // on the notch edge: inside
  EXPECT_TRUE(lShape.contains(1.5, 1.0));   // on the notch edge: inside
}

TEST(Polygon, TriangleAndDegenerateRejection) {
  const Polygon triangle({{0.0, 0.0}, {4.0, 0.0}, {0.0, 3.0}});
  EXPECT_TRUE(triangle.contains(1.0, 1.0));
  EXPECT_FALSE(triangle.contains(3.0, 2.0));
  EXPECT_TRUE(triangle.contains(2.0, 1.5));  // on the hypotenuse
  EXPECT_THROW(Polygon({{0.0, 0.0}, {1.0, 1.0}}), frehg::FatalError);
}

TEST(Polygon, BoundingBox) {
  const Polygon triangle({{-1.0, 2.0}, {4.0, 0.0}, {0.0, 3.0}});
  const auto& bbox = triangle.boundingBox();
  EXPECT_DOUBLE_EQ(bbox[0], -1.0);
  EXPECT_DOUBLE_EQ(bbox[1], 0.0);
  EXPECT_DOUBLE_EQ(bbox[2], 4.0);
  EXPECT_DOUBLE_EQ(bbox[3], 3.0);
}

// ---------------------------------------------------------------------------
// BoundarySet rasterization (single rank; rank-spanning cases in tests/mpi)
// ---------------------------------------------------------------------------

frehg::DomainConfig b5LikeDomain() {
  frehg::DomainConfig dom;
  dom.nx = 101;
  dom.ny = 55;
  dom.nz = 2;
  dom.dx = 1.0;
  dom.dy = 1.0;
  dom.dz = 0.2;
  return dom;
}

frehg::BoundaryConditionConfig outletConfig() {
  frehg::BoundaryConditionConfig bc;
  bc.name = "outlet";
  // b5's outlet polygon from the plan §6 example.
  bc.polygon = {{-0.1, 50.1}, {0.9, 50.1}, {0.9, 55.1}, {-0.1, 55.1}};
  bc.target = frehg::BcTarget::Surface;
  bc.kind = frehg::BcKind::Eta;
  bc.value.form = frehg::BcValueConfig::Form::Constant;
  bc.value.constant = -0.5;
  return bc;
}

TEST(BoundarySet, RasterizesB5OutletFootprint) {
  const frehg::Grid grid(MPI_COMM_SELF, b5LikeDomain());
  const frehg::BoundarySet set(grid, {outletConfig()}, ".");
  const frehg::BoundaryCondition& outlet = set.byName("outlet");

  // Cell centers x = 0.5 (i = 0), y = 50.5 .. 54.5 (j = 50..54): 5 cells.
  EXPECT_EQ(outlet.globalCellCount(), 5);
  ASSERT_EQ(outlet.cells().size(), 5u);
  for (const frehg::BcCell& cell : outlet.cells()) {
    EXPECT_EQ(cell.i, 1);               // local interior index of global i = 0
    EXPECT_GE(cell.j, 51);              // local interior index of global j = 50
    EXPECT_LE(cell.j, 55);
    EXPECT_EQ(cell.face, frehg::BcFace::Top);
  }
  EXPECT_TRUE(outlet.active());
  EXPECT_NE(outlet.subComm(), MPI_COMM_NULL);
  EXPECT_DOUBLE_EQ(outlet.value(0.0), -0.5);
  EXPECT_EQ(outlet.kind(), frehg::BcKind::Eta);
  EXPECT_EQ(outlet.target(), frehg::BcTarget::Surface);
}

TEST(BoundarySet, SideTargetSelectsEdgeFaces) {
  frehg::DomainConfig dom = b5LikeDomain();
  dom.nx = 4;
  dom.ny = 3;
  const frehg::Grid grid(MPI_COMM_SELF, dom);

  frehg::BoundaryConditionConfig bc;
  bc.name = "west-side";
  // Covers the full west column of cell centers (x = 0.5) plus interior
  // cells at x = 1.5, which must NOT contribute (side BCs rasterize onto
  // domain-edge cells only).
  bc.polygon = {{-0.1, -0.1}, {1.9, -0.1}, {1.9, 3.1}, {-0.1, 3.1}};
  bc.target = frehg::BcTarget::GroundwaterSide;
  bc.kind = frehg::BcKind::Flux;
  bc.value.form = frehg::BcValueConfig::Form::Gravity;

  const frehg::BoundarySet set(grid, {bc}, ".");
  const frehg::BoundaryCondition& side = set.byName("west-side");

  // 3 west-edge cells, plus the south/north edge faces of the two polygon
  // columns: j = 0 and j = ny-1 rows at i = 0, 1.
  int westFaces = 0;
  int southFaces = 0;
  int northFaces = 0;
  for (const frehg::BcCell& cell : side.cells()) {
    if (cell.face == frehg::BcFace::XMinus) {
      ++westFaces;
      EXPECT_EQ(cell.i, 1);
    }
    if (cell.face == frehg::BcFace::YMinus) {
      ++southFaces;
    }
    if (cell.face == frehg::BcFace::YPlus) {
      ++northFaces;
    }
  }
  EXPECT_EQ(westFaces, 3);
  EXPECT_EQ(southFaces, 2);
  EXPECT_EQ(northFaces, 2);
  EXPECT_DOUBLE_EQ(side.value(123.0), 0.0);  // gravity form carries no value
}

/// A side condition whose polygon is a small box around one cell center
/// (the v2 Q7 §8.2 single-cell rasterization rows).
frehg::BoundaryConditionConfig sideCellConfig(
    const std::string& name, std::vector<std::array<frehg::real_t, 2>> polygon) {
  frehg::BoundaryConditionConfig bc;
  bc.name = name;
  bc.polygon = std::move(polygon);
  bc.target = frehg::BcTarget::GroundwaterSide;
  bc.kind = frehg::BcKind::Flux;
  bc.value.form = frehg::BcValueConfig::Form::Gravity;
  return bc;
}

int faceCount(const frehg::BoundaryCondition& bc, frehg::BcFace face) {
  int count = 0;
  for (const frehg::BcCell& cell : bc.cells()) {
    if (cell.face == face) {
      ++count;
    }
  }
  return count;
}

TEST(BoundarySet, SingleCellPolygonOnEachEdge) {
  frehg::DomainConfig dom = b5LikeDomain();
  dom.nx = 4;
  dom.ny = 3;
  const frehg::Grid grid(MPI_COMM_SELF, dom);

  // Four polygons, each a +-0.2 box around exactly one NON-corner edge
  // cell's center of the 4x3 domain (centers x = 0.5..3.5, y = 0.5..2.5):
  // each must select exactly that cell with exactly that edge's face.
  const frehg::BoundarySet set(
      grid,
      {sideCellConfig("west", {{0.3, 1.3}, {0.7, 1.3}, {0.7, 1.7}, {0.3, 1.7}}),    // (0.5, 1.5)
       sideCellConfig("east", {{3.3, 1.3}, {3.7, 1.3}, {3.7, 1.7}, {3.3, 1.7}}),    // (3.5, 1.5)
       sideCellConfig("south", {{2.3, 0.3}, {2.7, 0.3}, {2.7, 0.7}, {2.3, 0.7}}),   // (2.5, 0.5)
       sideCellConfig("north", {{1.3, 2.3}, {1.7, 2.3}, {1.7, 2.7}, {1.3, 2.7}})},  // (1.5, 2.5)
      ".");

  struct Expected {
    const char* name;
    frehg::BcFace face;
    int j, i;  // local interior indices (halo layout)
  };
  const Expected rows[] = {
      {"west", frehg::BcFace::XMinus, 2, 1},
      {"east", frehg::BcFace::XPlus, 2, 4},
      {"south", frehg::BcFace::YMinus, 1, 3},
      {"north", frehg::BcFace::YPlus, 3, 2},
  };
  for (const Expected& e : rows) {
    const frehg::BoundaryCondition& bc = set.byName(e.name);
    EXPECT_EQ(bc.globalCellCount(), 1) << e.name;
    ASSERT_EQ(bc.cells().size(), 1u) << e.name;
    EXPECT_EQ(bc.cells()[0].face, e.face) << e.name;
    EXPECT_EQ(bc.cells()[0].j, e.j) << e.name;
    EXPECT_EQ(bc.cells()[0].i, e.i) << e.name;
  }
}

TEST(BoundarySet, CornerCellPolygonSpillsBothFaces) {
  // A side polygon holding exactly one CORNER cell's center hands the
  // condition BOTH domain-edge faces of that cell, one member per face —
  // the documented corner spill b6's one-cell-wide tank relies on
  // (docs/theory/temperature.md, "Side-ghost limiter admission").
  frehg::DomainConfig dom = b5LikeDomain();
  dom.nx = 4;
  dom.ny = 3;
  const frehg::Grid grid(MPI_COMM_SELF, dom);

  const frehg::BoundarySet set(
      grid,
      {sideCellConfig("sw", {{0.3, 0.3}, {0.7, 0.3}, {0.7, 0.7}, {0.3, 0.7}}),   // (0.5, 0.5)
       sideCellConfig("se", {{3.3, 0.3}, {3.7, 0.3}, {3.7, 0.7}, {3.3, 0.7}}),   // (3.5, 0.5)
       sideCellConfig("nw", {{0.3, 2.3}, {0.7, 2.3}, {0.7, 2.7}, {0.3, 2.7}}),   // (0.5, 2.5)
       sideCellConfig("ne", {{3.3, 2.3}, {3.7, 2.3}, {3.7, 2.7}, {3.3, 2.7}})},  // (3.5, 2.5)
      ".");

  struct Expected {
    const char* name;
    frehg::BcFace faceA, faceB;
    int j, i;  // local interior indices of the corner cell
  };
  const Expected rows[] = {
      {"sw", frehg::BcFace::XMinus, frehg::BcFace::YMinus, 1, 1},
      {"se", frehg::BcFace::XPlus, frehg::BcFace::YMinus, 1, 4},
      {"nw", frehg::BcFace::XMinus, frehg::BcFace::YPlus, 3, 1},
      {"ne", frehg::BcFace::XPlus, frehg::BcFace::YPlus, 3, 4},
  };
  for (const Expected& e : rows) {
    const frehg::BoundaryCondition& bc = set.byName(e.name);
    // Exactly the two edge faces of the corner cell — nothing else.
    EXPECT_EQ(bc.globalCellCount(), 2) << e.name;
    ASSERT_EQ(bc.cells().size(), 2u) << e.name;
    EXPECT_EQ(faceCount(bc, e.faceA), 1) << e.name;
    EXPECT_EQ(faceCount(bc, e.faceB), 1) << e.name;
    for (const frehg::BcCell& cell : bc.cells()) {
      EXPECT_EQ(cell.j, e.j) << e.name;
      EXPECT_EQ(cell.i, e.i) << e.name;
    }
  }
}

TEST(BoundarySet, EmptyRegionIsFatal) {
  const frehg::Grid grid(MPI_COMM_SELF, b5LikeDomain());
  frehg::BoundaryConditionConfig bc = outletConfig();
  bc.polygon = {{500.0, 500.0}, {501.0, 500.0}, {501.0, 501.0}};  // off-domain
  EXPECT_THROW(frehg::BoundarySet(grid, {bc}, "."), frehg::FatalError);
}

TEST(BoundarySet, UnknownNameIsFatal) {
  const frehg::Grid grid(MPI_COMM_SELF, b5LikeDomain());
  const frehg::BoundarySet set(grid, {outletConfig()}, ".");
  EXPECT_THROW(set.byName("inlet"), frehg::FatalError);
}

}  // namespace
