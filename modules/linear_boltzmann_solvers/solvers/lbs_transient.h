// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"

namespace opensn
{
class TimeIntegration;

class TransientSolver : public opensn::Solver
{
protected:
  std::shared_ptr<LBSProblem> lbs_problem_;
  std::shared_ptr<TimeIntegration> time_integration_;

public:
  static InputParameters GetInputParameters();
  explicit TransientSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
  void Step() override;
  void Advance() override;
};

} // namespace opensn
