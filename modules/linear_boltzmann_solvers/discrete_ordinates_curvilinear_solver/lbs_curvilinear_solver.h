// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

namespace opensn
{

/**
 * A neutral particle transport solver in point-symmetric and axial-symmetric curvilinear
 * coordinates.
 */
class DiscreteOrdinatesCurvilinearSolver : public DiscreteOrdinatesSolver
{
private:
  /// Coordinate system type.
  CoordinateSystemType coord_system_type_;
  /** Discretisation pointer to matrices of the secondary cell view  (matrices of the primary cell
   * view forwarded to the base class).
   */
  std::shared_ptr<opensn::SpatialDiscretization> discretization_secondary_;
  std::vector<UnitCellMatrices> secondary_unit_cell_matrices_;

public:
  explicit DiscreteOrdinatesCurvilinearSolver(const InputParameters& params);

  DiscreteOrdinatesCurvilinearSolver(const DiscreteOrdinatesCurvilinearSolver&) = delete;
  DiscreteOrdinatesCurvilinearSolver& operator=(const DiscreteOrdinatesCurvilinearSolver&) = delete;

protected:
  void PerformInputChecks() override;
  void InitializeSpatialDiscretization() override;
  void ComputeSecondaryUnitIntegrals();

private:
  std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset) override;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<DiscreteOrdinatesCurvilinearSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
