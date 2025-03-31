// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/physics/solver.h"

namespace opensn
{
class TimeIntegration;
class LBSProblem;

class TransientSolver : public Solver
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
