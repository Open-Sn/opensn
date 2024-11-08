// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/linear_solver/linear_solver.h"

namespace opensn
{

std::string
LinearSolver::IterativeMethodName(LinearSolver::IterativeMethod iterative_method)
{
  switch (iterative_method)
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
      throw std::runtime_error("Unrecognized iterative method.");
  }
}

} // namespace opensn
