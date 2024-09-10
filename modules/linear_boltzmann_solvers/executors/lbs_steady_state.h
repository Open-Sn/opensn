// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"

namespace opensn
{

class SteadyStateSolver : public opensn::Solver
{
protected:
  std::shared_ptr<LBSSolver> lbs_solver_;

public:
  explicit SteadyStateSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<SteadyStateSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
