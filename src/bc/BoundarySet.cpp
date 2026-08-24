/// \file BoundarySet.cpp
/// \brief Rasterization of polygon boundary regions and sub-communicator
///        creation.

#include "bc/BoundarySet.hpp"

#include "core/Logger.hpp"

#include <utility>

namespace frehg {

real_t BoundaryCondition::value(real_t t) const {
  switch (form_) {
    case BcValueConfig::Form::Constant:
      return constant_;
    case BcValueConfig::Form::Series:
      return series_.value(t);
    case BcValueConfig::Form::Gravity:
      return 0.0;
    case BcValueConfig::Form::Hydrostatic:
      return hydrostaticEta_;
  }
  return constant_;
}

namespace {

/// Append the domain-edge faces of an owned cell that lie inside \p poly.
void appendSideFaces(const Grid& grid, const Polygon& poly, int jLoc, int iLoc,
                     std::vector<BcCell>& cells) {
  const int iGlob = grid.i0() + iLoc - 1;
  const int jGlob = grid.j0() + jLoc - 1;
  const real_t xc = grid.xCenter(iGlob);
  const real_t yc = grid.yCenter(jGlob);
  if (!poly.contains(xc, yc)) {
    return;
  }
  if (iGlob == 0) {
    cells.push_back(BcCell{jLoc, iLoc, BcFace::XMinus});
  }
  if (iGlob == grid.nx() - 1) {
    cells.push_back(BcCell{jLoc, iLoc, BcFace::XPlus});
  }
  if (jGlob == 0) {
    cells.push_back(BcCell{jLoc, iLoc, BcFace::YMinus});
  }
  if (jGlob == grid.ny() - 1) {
    cells.push_back(BcCell{jLoc, iLoc, BcFace::YPlus});
  }
}

}  // namespace

BoundarySet::BoundarySet(const Grid& grid, const std::vector<BoundaryConditionConfig>& configs,
                         const std::string& configDir) {
  int worldRank = 0;
  MPI_Comm_rank(grid.comm(), &worldRank);

  for (const BoundaryConditionConfig& cfg : configs) {
    BoundaryCondition bc;
    bc.name_ = cfg.name;
    bc.target_ = cfg.target;
    bc.kind_ = cfg.kind;
    bc.form_ = cfg.value.form;
    bc.constant_ = cfg.value.constant;
    bc.hydrostaticEta_ = cfg.value.hydrostaticEta;
    if (cfg.value.form == BcValueConfig::Form::Series) {
      const std::string seriesPath = configDir.empty()
                                         ? cfg.value.seriesFile
                                         : configDir + "/" + cfg.value.seriesFile;
      bc.series_ = TimeSeries::fromFile(seriesPath);
    }

    const Polygon poly(cfg.polygon);
    const BcFace footprintFace = (cfg.target == BcTarget::GroundwaterBottom)
                                     ? BcFace::Bottom
                                     : BcFace::Top;
    // Side targets, surface velocity, and surface outflow conditions act on
    // domain-edge faces; everything else selects the (j, i) footprint.
    const bool edgeFaces =
        cfg.target == BcTarget::GroundwaterSide ||
        (cfg.target == BcTarget::Surface &&
         (cfg.kind == BcKind::Velocity || cfg.kind == BcKind::Outflow));
    for (int jLoc = 1; jLoc <= grid.nyLocal(); ++jLoc) {
      for (int iLoc = 1; iLoc <= grid.nxLocal(); ++iLoc) {
        if (edgeFaces) {
          appendSideFaces(grid, poly, jLoc, iLoc, bc.cells_);
          continue;
        }
        const int iGlob = grid.i0() + iLoc - 1;
        const int jGlob = grid.j0() + jLoc - 1;
        if (poly.contains(grid.xCenter(iGlob), grid.yCenter(jGlob))) {
          bc.cells_.push_back(BcCell{jLoc, iLoc, footprintFace});
        }
      }
    }

    // Sub-communicator over ranks owning members (plan §5.6).
    const int color = bc.cells_.empty() ? MPI_UNDEFINED : 1;
    MPI_Comm sub = MPI_COMM_NULL;
    MPI_Comm_split(grid.comm(), color, worldRank, &sub);
    bc.subComm_ = sub;

    const long localCells = static_cast<long>(bc.cells_.size());
    long globalCells = 0;
    MPI_Allreduce(&localCells, &globalCells, 1, MPI_LONG, MPI_SUM, grid.comm());
    bc.globalCells_ = globalCells;
    if (globalCells == 0) {
      log::fatal(log::msg() << "boundary condition '" << cfg.name
                            << "': its polygon selects no cell anywhere in the domain; "
                               "check the region coordinates against cell centers");
    }

    conditions_.push_back(std::move(bc));
  }
}

BoundarySet::~BoundarySet() {
  for (BoundaryCondition& bc : conditions_) {
    if (bc.subComm_ != MPI_COMM_NULL) {
      MPI_Comm_free(&bc.subComm_);
    }
  }
}

const BoundaryCondition& BoundarySet::byName(const std::string& name) const {
  for (const BoundaryCondition& bc : conditions_) {
    if (bc.name() == name) {
      return bc;
    }
  }
  log::fatal(log::msg() << "BoundarySet: no boundary condition named '" << name << "'");
}

}  // namespace frehg
