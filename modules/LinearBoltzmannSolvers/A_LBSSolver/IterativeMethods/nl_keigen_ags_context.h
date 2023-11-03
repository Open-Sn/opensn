#pragma once

#include "math/NonLinearSolver/NonLinearSolverContext.h"

#include "A_LBSSolver/lbs_solver.h"
#include "A_LBSSolver/IterativeMethods/snes_k_residual_func_context.h"

#include <vector>
#include <cstdint>

namespace lbs
{
class LBSSolver;
}

namespace lbs
{

template <class VecType, class SolverType>
struct NLKEigenAGSContext : public chi_math::NonLinearSolverContext<VecType, SolverType>
{
  LBSSolver& lbs_solver_;
  KResidualFunctionContext kresid_func_context_;

  std::vector<int> groupset_ids;

  explicit NLKEigenAGSContext(LBSSolver& lbs_solver)
    : lbs_solver_(lbs_solver), kresid_func_context_({lbs_solver.TextName(), 1.0})
  {
  }

  virtual ~NLKEigenAGSContext() override = default;
};

} // namespace lbs


