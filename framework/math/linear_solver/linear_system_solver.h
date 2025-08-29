// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver.h"
#include "framework/math/linear_solver/linear_solver_context.h"

namespace opensn
{

class LinearSystemSolver : public LinearSolver
{
public:
  enum class IterativeMethod : int
  {
    NONE,
    CLASSIC_RICHARDSON,
    PETSC_RICHARDSON,
    PETSC_GMRES,
    PETSC_BICGSTAB
  };

  LinearSystemSolver(IterativeMethod method, const std::shared_ptr<LinearSystemContext>& ctx)
    : LinearSolver(ctx), method_(method)
  {
  }

  std::string GetIterativeMethodName()
  {
    switch (method_)
    {
      case IterativeMethod::NONE:
        return "NONE";
      case IterativeMethod::CLASSIC_RICHARDSON:
        return "CLASSIC_RICHARDSON";
      case IterativeMethod::PETSC_RICHARDSON:
        return "PETSC_RICHARDSON";
      case IterativeMethod::PETSC_GMRES:
        return "PETSC_GMRES";
      case IterativeMethod::PETSC_BICGSTAB:
        return "PETSC_BICGSTAB";
      default:
        return "UNKNOWN ITERATIVE METHOD";
    }
  }

protected:
  IterativeMethod method_;
};

} // namespace opensn
