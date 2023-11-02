#ifndef CHITECH_SOLVERINFOPOSTPROCESSOR_H
#define CHITECH_SOLVERINFOPOSTPROCESSOR_H

#include "PostProcessor.h"

namespace chi_physics
{
class Solver;
}

namespace chi
{

class SolverInfoPostProcessor : public PostProcessor
{
public:
  static InputParameters GetInputParameters();
  explicit SolverInfoPostProcessor(const InputParameters& params);

  void Execute(const Event& event_context) override;

private:
  const chi_physics::Solver& solver_;
  const ParameterBlock info_;
};

}

#endif // CHITECH_SOLVERINFOPOSTPROCESSOR_H
