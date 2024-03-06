#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"

namespace opensn
{
namespace lbs
{

class SteadyStateSolver : public opensn::Solver
{
protected:
  LBSSolver& lbs_solver_;

public:
  static InputParameters GetInputParameters();

  explicit SteadyStateSolver(const InputParameters& params);

  void Initialize() override;
  void Execute() override;
};

} // namespace lbs
} // namespace opensn
