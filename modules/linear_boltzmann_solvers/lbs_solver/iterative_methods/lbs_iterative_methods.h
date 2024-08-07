// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>

namespace opensn
{

enum class IterativeMethod : int
{
  NONE = 0,
  CLASSIC_RICHARDSON = 1, ///< Classic Richardson (source iteration)
  KRYLOV_GMRES = 2,       ///< PETSc GMRES iterative algorithm
  KRYLOV_RICHARDSON = 3,  ///< PETSc Richardson iteration
  KRYLOV_BICGSTAB = 4,    ///< PETSc BiCGStab iterative algorithm
};

inline std::string
IterativeMethodPETScName(IterativeMethod it_method)
{
  switch (it_method)
  {
    case IterativeMethod::NONE:
      return "preonly";
    case IterativeMethod::KRYLOV_RICHARDSON:
      return "richardson";
    case IterativeMethod::KRYLOV_GMRES:
      return "gmres";
    case IterativeMethod::KRYLOV_BICGSTAB:
      return "bcgs";
    case IterativeMethod::CLASSIC_RICHARDSON:
      return "";
  }
  return "";
}

} // namespace opensn
