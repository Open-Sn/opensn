// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/solvers/discrete_ordinates_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"

namespace opensn
{

/// Steady-state source solver that drives the across-groupset (AGS) solver.
class SteadyStateSourceSolver : public DiscreteOrdinatesSolver
{
public:
  explicit SteadyStateSourceSolver(const InputParameters& params);

  BalanceTable ComputeBalanceTable() const;

protected:
  void CheckRequirements(std::vector<std::string>& errors) const override;
  void InitializeSolver() override;
  void ExecuteSolver() override;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<SteadyStateSourceSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
