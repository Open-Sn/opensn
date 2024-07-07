// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"
#include "framework/math/linear_solver/linear_solver.h"
#include <vector>
#include <memory>
#include <petscksp.h>

namespace opensn
{

class LBSSolver;

struct AGSContext : public LinearSolverContext
{
  LBSSolver& lbs_solver_;
  std::vector<std::shared_ptr<LinearSolver>> wgs_solvers_;

  AGSContext(LBSSolver& lbs_solver, std::vector<std::shared_ptr<LinearSolver>> wgs_solvers)
    : lbs_solver_(lbs_solver), wgs_solvers_(std::move(wgs_solvers))
  {
  }
};

} // namespace opensn
