// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/solver.h"
#include <memory>
#include <string>

namespace opensn
{

class UncollidedProblem;

class UncollidedSolver : public Solver
{
public:
  explicit UncollidedSolver(const InputParameters& params);

  void Initialize() override;

  void Execute() override;

protected:
  std::shared_ptr<UncollidedProblem> problem_;
  std::string file_name_;
  unsigned int progress_interval_ = 5;
  bool initialized_ = false;

public:
  static InputParameters GetInputParameters();

  static std::shared_ptr<UncollidedSolver> Create(const ParameterBlock& params);
};

} // namespace opensn
