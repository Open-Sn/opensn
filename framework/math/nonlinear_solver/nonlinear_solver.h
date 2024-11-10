// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/nonlinear_solver/nonlinear_solver_context.h"
#include <string>
#include <memory>

namespace opensn
{

/// Implementation of a general non-linear solver.
class NonLinearSolver
{
public:
  explicit NonLinearSolver(std::string solver_name,
                           std::shared_ptr<NonLinearSolverContext> context_ptr)
    : solver_name_(solver_name), context_ptr_(context_ptr)
  {
  }

  virtual ~NonLinearSolver() {}

  std::shared_ptr<NonLinearSolverContext> GetContext() { return context_ptr_; }

  virtual void Setup() {}

  virtual void Solve() = 0;

protected:
  const std::string solver_name_;
  std::shared_ptr<NonLinearSolverContext> context_ptr_;
};

} // namespace opensn
