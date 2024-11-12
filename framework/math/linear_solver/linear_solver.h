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
  enum class IterativeMethod : int
  {
    NONE = 0,
    CLASSIC_RICHARDSON = 1, ///< Classic Richardson (source iteration)
    PETSC_RICHARDSON = 3,   ///< PETSc Richardson iteration
    PETSC_GMRES = 2,        ///< PETSc GMRES iterative algorithm
    PETSC_BICGSTAB = 4,     ///< PETSc BiCGStab iterative algorithm
  };

  LinearSolver(IterativeMethod iterative_method, std::shared_ptr<LinearSolverContext> context_ptr)
    : iterative_method_(iterative_method), context_ptr_(context_ptr)
  {
  }

  virtual ~LinearSolver() {}

  std::shared_ptr<LinearSolverContext> GetContext() { return context_ptr_; }

  /// Set up the linaer solver
  virtual void Setup() {}

  /// Solve the system
  virtual void Solve() = 0;

protected:
  const IterativeMethod iterative_method_;
  std::shared_ptr<LinearSolverContext> context_ptr_;

public:
  static std::string IterativeMethodName(IterativeMethod iterative_method);
};

} // namespace opensn
