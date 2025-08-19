// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver_context.h"
#include <memory>

namespace opensn
{

struct LinearSolverContext;

class LinearSolver
{
public:
  explicit LinearSolver(std::shared_ptr<LinearSolverContext> context_ptr)
    : context_ptr_(std::move(context_ptr))
  {
  }

  virtual ~LinearSolver() = default;

  /// Set up the linaer solver
  virtual void Setup() {}

  /// Solve the system
  virtual void Solve() = 0;

  std::shared_ptr<LinearSolverContext> GetContext() { return context_ptr_; }

protected:
  std::shared_ptr<LinearSolverContext> context_ptr_;
};

} // namespace opensn
