// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include <memory>

namespace opensn
{
class LBSProblem;

/**
 * A solver that drives a discrete ordinates problem through a time loop.
 */
class TimeDependentSourceSolver : public Solver
{
public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<TimeDependentSourceSolver> Create(const ParameterBlock& params);

  explicit TimeDependentSourceSolver(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

private:
  std::shared_ptr<LBSProblem> lbs_problem_;
  double dt_ = 1.0;
  double theta_ = 1.0;
  double stop_time_ = 1.0;
};

} // namespace opensn
