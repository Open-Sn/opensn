// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/snes_k_residual_func_context.h"
#include "framework/math/nonlinear_solver/nonlinear_solver_context.h"
#include <vector>
#include <cstdint>

namespace opensn
{

class LBSProblem;

struct NLKEigenAGSContext : public NonLinearSolverContext
{
  explicit NLKEigenAGSContext(const std::shared_ptr<LBSProblem>& lbs_problem)
    : lbs_problem(lbs_problem), kresid_func_context({lbs_problem->GetName(), 1.0})
  {
  }

  ~NLKEigenAGSContext() override = default;

  std::shared_ptr<LBSProblem> lbs_problem;
  KResidualFunctionContext kresid_func_context;

  std::vector<int> groupset_ids;
};

} // namespace opensn
