// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/compute/discrete_ordinates_compute.h"

namespace opensn
{

class DiscreteOrdinatesProblem;

class SteadyStateSourceSolver : public Solver
{
public:
  explicit SteadyStateSourceSolver(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

  BalanceTable ComputeBalanceTable() const;

protected:
  std::shared_ptr<DiscreteOrdinatesProblem> do_problem_;
  bool initialized_ = false;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<SteadyStateSourceSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
