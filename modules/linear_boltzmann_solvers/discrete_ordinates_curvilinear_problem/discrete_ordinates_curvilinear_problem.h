// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

/**
 * A neutral particle transport solver in point-symmetric and axial-symmetric curvilinear
 * coordinates.
 */
class DiscreteOrdinatesCurvilinearProblem : public DiscreteOrdinatesProblem
{
private:
  /** Discretisation pointer to matrices of the secondary cell view  (matrices of the primary cell
   * view forwarded to the base class).
   */
  std::shared_ptr<opensn::SpatialDiscretization> discretization_secondary_;
  std::vector<UnitCellMatrices> secondary_unit_cell_matrices_;

public:
  explicit DiscreteOrdinatesCurvilinearProblem(const InputParameters& params);

  DiscreteOrdinatesCurvilinearProblem(const DiscreteOrdinatesCurvilinearProblem&) = delete;
  DiscreteOrdinatesCurvilinearProblem&
  operator=(const DiscreteOrdinatesCurvilinearProblem&) = delete;

protected:
  void PerformInputChecks() override;
  void InitializeSpatialDiscretization() override;
  void ComputeSecondaryUnitIntegrals();

private:
  std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset) override;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<DiscreteOrdinatesCurvilinearProblem> Create(const ParameterBlock& params);
};

} // namespace opensn
