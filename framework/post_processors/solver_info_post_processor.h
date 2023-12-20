#pragma once

#include "framework/post_processors/post_processor.h"

namespace opensn
{
class Solver;

class SolverInfoPostProcessor : public PostProcessor
{
public:
  static InputParameters GetInputParameters();
  explicit SolverInfoPostProcessor(const InputParameters& params);

  void Execute(const Event& event_context) override;

private:
  const Solver& solver_;
  const ParameterBlock info_;
};

} // namespace opensn
