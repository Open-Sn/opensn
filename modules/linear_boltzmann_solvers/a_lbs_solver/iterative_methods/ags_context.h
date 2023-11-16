#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"
#include "framework/math/linear_solver/linear_solver.h"
#include <vector>
#include <memory>
#include <petscksp.h>

namespace opensn
{
namespace lbs
{
class LBSSolver;

struct AGSContext : public LinearSolverContext
{
  typedef std::shared_ptr<LinearSolver> LinSolveBaseTypePtr;
  LBSSolver& lbs_solver_;
  std::vector<LinSolveBaseTypePtr> sub_solvers_list_;

  AGSContext(LBSSolver& lbs_solver, std::vector<LinSolveBaseTypePtr> sub_solvers_list)
    : lbs_solver_(lbs_solver), sub_solvers_list_(std::move(sub_solvers_list))
  {
  }

  std::pair<int64_t, int64_t> SystemSize();

  virtual void SetPreconditioner(KSP& solver);

  int MatrixAction(Mat& matrix, Vec& vector, Vec& action) override;
};

} // namespace lbs
} // namespace opensn
