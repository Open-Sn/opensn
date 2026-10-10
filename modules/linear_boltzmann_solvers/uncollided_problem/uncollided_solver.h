// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include <memory>
#include <string>

namespace opensn
{

class UncollidedProblem;

/// Solver that generates uncollided flux moments for first-collision transport.
class UncollidedSolver : public Solver
{
public:
  explicit UncollidedSolver(const InputParameters& params);

protected:
  /// Validates the problem's configuration; the solver has no further requirements.
  void ValidateState() const override;
  void InitializeSolver() override;
  void ExecuteSolver() override;

  std::shared_ptr<UncollidedProblem> problem_;
  std::string file_name_;
  unsigned int progress_interval_ = 5;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<UncollidedSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
