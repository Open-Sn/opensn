// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/snes_k_residual_func_context.h"
#include "framework/math/nonlinear_solver/nonlinear_solver_context.h"
#include <vector>
#include <cstdint>

namespace opensn
{

class LBSSolver;

struct NLKEigenAGSContext : public NonLinearSolverContext
{
  std::shared_ptr<LBSSolver> lbs_solver;
  KResidualFunctionContext kresid_func_context;

  std::vector<int> groupset_ids;

  explicit NLKEigenAGSContext(std::shared_ptr<LBSSolver> lbs_solver)
    : lbs_solver(lbs_solver), kresid_func_context({lbs_solver->GetName(), 1.0})
  {
  }

  ~NLKEigenAGSContext() override = default;
};

} // namespace opensn
