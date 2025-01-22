// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/post_processors/post_processor.h"
#include <memory>

namespace opensn
{
class Solver;

class SolverInfoPostProcessor : public PostProcessor
{
public:
  explicit SolverInfoPostProcessor(const InputParameters& params);

  void Execute(const Event& event_context) override;

private:
  std::shared_ptr<Solver> solver_;
  const ParameterBlock info_;

public:
  static InputParameters GetInputParameters();
  static std::shared_ptr<SolverInfoPostProcessor> Create(const ParameterBlock& params);
};

} // namespace opensn
