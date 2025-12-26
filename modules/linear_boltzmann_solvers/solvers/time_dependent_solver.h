// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include <cstddef>
#include <functional>
#include <memory>

namespace opensn
{
class LBSProblem;
class AGSLinearSolver;

/**
 * A solver that drives a discrete ordinates problem through a time loop.
 */
class TimeDependentSourceSolver : public Solver
{
public:
  explicit TimeDependentSourceSolver(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

  void Advance() override;

  void SetPreAdvanceCallback(std::function<void()> callback);

  void SetPreAdvanceCallback(std::nullptr_t);

  void SetPostAdvanceCallback(std::function<void()> callback);

  void SetPostAdvanceCallback(std::nullptr_t);

  void SetTimeStep(double dt);

  void SetTheta(double theta);

  double GetCurrentTime() const { return current_time_; }

  unsigned int GetStep() const { return step_; }

private:
  std::shared_ptr<LBSProblem> lbs_problem_;
  std::shared_ptr<AGSLinearSolver> ags_solver_;
  double stop_time_ = 1.0;
  double current_time_ = 0.0;
  unsigned int step_ = 0;
  bool verbose_ = true;
  std::function<void()> pre_advance_callback_;
  std::function<void()> post_advance_callback_;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<TimeDependentSourceSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
