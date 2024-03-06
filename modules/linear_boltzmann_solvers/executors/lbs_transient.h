#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"

namespace opensn
{
class TimeIntegration;

namespace lbs
{

class TransientSolver : public opensn::Solver
{
protected:
  LBSSolver& lbs_solver_;
  std::shared_ptr<TimeIntegration> time_integration_;

public:
  static InputParameters GetInputParameters();
  explicit TransientSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
  void Step() override;
  void Advance() override;
};

} // namespace lbs
} // namespace opensn
