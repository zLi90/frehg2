/// \file BoundarySet.hpp
/// \brief Rasterized boundary conditions with per-BC sub-communicators
///        (plan §5.6, adopted from SERGHEI).
///
/// A BoundarySet turns the configuration's polygon regions into per-rank
/// member-cell lists at construction:
///  - target surface / groundwater_top / groundwater_bottom: every owned
///    cell whose center lies inside the polygon (the (j, i) footprint);
///  - target groundwater_side and surface velocity conditions: owned cells
///    on a domain edge whose center lies inside the polygon, one member per
///    edge face the cell touches.
///
/// Each condition gets an MPI sub-communicator spanning the ranks that own
/// at least one member cell, for cross-section reductions such as inflow
/// distribution and monitor fluxes.

#ifndef FREHG_BC_BOUNDARYSET_HPP
#define FREHG_BC_BOUNDARYSET_HPP

#include "bc/BoundaryKinds.hpp"
#include "bc/Polygon.hpp"
#include "core/Config.hpp"
#include "core/Grid.hpp"
#include "core/TimeSeries.hpp"

#include <mpi.h>

#include <string>
#include <vector>

namespace frehg {

/// One boundary condition after rasterization onto this rank.
class BoundaryCondition {
 public:
  /// \name Identity and classification
  ///@{
  /// The condition's unique configuration name.
  const std::string& name() const { return name_; }
  /// The sub-boundary this condition applies to.
  BcTarget target() const { return target_; }
  /// The prescribed quantity.
  BcKind kind() const { return kind_; }
  /// Value form (constant / series / gravity / hydrostatic).
  BcValueConfig::Form valueForm() const { return form_; }
  ///@}

  /// Evaluate the prescribed value at time \p t. For Gravity the value is 0
  /// (the flux is computed from local conductivity by the physics); for
  /// Hydrostatic it is the reference surface elevation.
  real_t value(real_t t) const;

  /// Member cells owned by this rank.
  const std::vector<BcCell>& cells() const { return cells_; }

  /// Sub-communicator over ranks with members (MPI_COMM_NULL elsewhere).
  MPI_Comm subComm() const { return subComm_; }

  /// \return true when this rank owns at least one member cell.
  bool active() const { return !cells_.empty(); }

  /// Total member cells across all ranks (counted at construction).
  long globalCellCount() const { return globalCells_; }

 private:
  friend class BoundarySet;
  std::string name_;
  BcTarget target_ = BcTarget::Surface;
  BcKind kind_ = BcKind::Eta;
  BcValueConfig::Form form_ = BcValueConfig::Form::Constant;
  real_t constant_ = 0.0;
  real_t hydrostaticEta_ = 0.0;
  TimeSeries series_;
  std::vector<BcCell> cells_;
  MPI_Comm subComm_ = MPI_COMM_NULL;
  long globalCells_ = 0;
};

/// All boundary conditions of a run, rasterized against a grid.
class BoundarySet {
 public:
  /// Rasterize the configured regions onto the grid. Collective on
  /// grid.comm().
  /// \param grid decomposition to rasterize against (kept by reference).
  /// \param configs the configuration's boundary_conditions list.
  /// \param configDir directory for resolving series files.
  /// A region that selects no cell on any rank is fatal (a mis-placed
  /// polygon would otherwise silently disable the condition).
  BoundarySet(const Grid& grid, const std::vector<BoundaryConditionConfig>& configs,
              const std::string& configDir);

  /// Frees the sub-communicators.
  ~BoundarySet();

  BoundarySet(const BoundarySet&) = delete;
  BoundarySet& operator=(const BoundarySet&) = delete;

  /// \return all conditions, in configuration order.
  const std::vector<BoundaryCondition>& all() const { return conditions_; }

  /// \return the condition with \p name (fatal if absent).
  const BoundaryCondition& byName(const std::string& name) const;

 private:
  std::vector<BoundaryCondition> conditions_;
};

}  // namespace frehg

#endif  // FREHG_BC_BOUNDARYSET_HPP
