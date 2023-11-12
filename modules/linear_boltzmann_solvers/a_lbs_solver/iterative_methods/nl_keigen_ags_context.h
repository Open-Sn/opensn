#pragma once

#include "framework/math/nonlinear_solver/nonlinear_solver_context.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/snes_k_residual_func_context.h"

#include <vector>
#include <cstdint>

namespace lbs
{
class LBSSolver;
}

namespace lbs
{

struct NLKEigenAGSContext : public chi_math::NonLinearSolverContext
{
  LBSSolver& lbs_solver_;
  KResidualFunctionContext kresid_func_context_;

  std::vector<int> groupset_ids;

  explicit NLKEigenAGSContext(LBSSolver& lbs_solver)
    : lbs_solver_(lbs_solver), kresid_func_context_({lbs_solver.TextName(), 1.0})
  {
  }

  ~NLKEigenAGSContext() override = default;
};

} // namespace lbs
