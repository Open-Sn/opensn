// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"

namespace opensn
{

class LBSProblem;

class SteadyStateSourceSolver : public Solver
{
protected:
  std::shared_ptr<LBSProblem> lbs_problem_;

public:
  explicit SteadyStateSourceSolver(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<SteadyStateSourceSolver> Create(const ParameterBlock& params);

private:
  bool ReadRestartData();

  bool WriteRestartData();
};

} // namespace opensn
