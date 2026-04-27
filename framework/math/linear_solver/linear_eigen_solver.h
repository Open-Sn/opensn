// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/linear_solver.h"
#include "framework/math/linear_solver/linear_solver_context.h"

namespace opensn
{

class LinearEigenSolver : public LinearSolver
{
public:
  enum class IterativeMethod : int
  {
    KRYLOV_SCHUR,
    POWER,
    ARNOLDI,
    SUBSPACE
  };

  LinearEigenSolver(IterativeMethod method, std::shared_ptr<LinearEigenContext> ctx)
    : LinearSolver(ctx), method_(method)
  {
  }

  std::string GetIterativeMethodName()
  {
    switch (method_)
    {
      case IterativeMethod::KRYLOV_SCHUR:
        return "KRYLOV-SCHUR";
      case IterativeMethod::POWER:
        return "POWER ITERATION";
      case IterativeMethod::ARNOLDI:
        return "ARNOLDI";
      case IterativeMethod::SUBSPACE:
        return "SUBSPACE";
      default:
        return "UNKNOWN ITERATIVE METHOD";
    }
  }

protected:
  IterativeMethod method_;
};

} // namespace opensn
