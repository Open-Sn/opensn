// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn
{

class LBSSteadyStateSolver : public opensn::Solver
{
protected:
  std::shared_ptr<LBSProblem> lbs_problem_;

public:
  explicit LBSSteadyStateSolver(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<LBSSteadyStateSolver> Create(const ParameterBlock& params);

private:
  bool ReadRestartData();

  bool WriteRestartData();
};

} // namespace opensn
