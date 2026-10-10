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
 * coordinates (i.e., spherical and cylindrical coordinates, respectively).
 */
class DiscreteOrdinatesCurvilinearProblem : public DiscreteOrdinatesProblem
{
public:
  DiscreteOrdinatesCurvilinearProblem(const DiscreteOrdinatesCurvilinearProblem&) = delete;
  DiscreteOrdinatesCurvilinearProblem&
  operator=(const DiscreteOrdinatesCurvilinearProblem&) = delete;

protected:
  /// Curvilinear-geometry rules.
  void CheckConfigurationErrors(std::vector<std::string>& errors) const override;
  void InitializeSpatialDiscretization() override;
  void ComputeSecondaryUnitIntegrals();
  bool SupportsTimeDependentMode() const override { return false; }
  bool IsCurvilinear() const override { return true; }
  std::shared_ptr<SweepChunk> SetSweepChunk(LBSGroupset& groupset) override;

private:
  /// Cell matrices of the angular-derivative terms (one power of r less than the primary ones).
  std::vector<UnitCellMatrices> secondary_unit_cell_matrices_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<DiscreteOrdinatesCurvilinearProblem> Create(const ParameterBlock& params);

  const std::vector<UnitCellMatrices>& GetSecondaryUnitCellMatrices() const;

private:
  /// Factory-only constructor.
  explicit DiscreteOrdinatesCurvilinearProblem(const InputParameters& params);
};

} // namespace opensn
