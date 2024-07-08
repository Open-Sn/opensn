// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <string>

namespace opensn
{

enum class IterativeMethod : int
{
  NONE = 0,
  CLASSIC_RICHARDSON = 1,        ///< Otherwise known as Source Iteration
  CLASSIC_RICHARDSON_CYCLES = 2, ///< Source Iteration with Cycles support
  GMRES = 3,                     ///< GMRES iterative algorithm
  GMRES_CYCLES = 4,              ///< GMRES with Cycles support
  KRYLOV_RICHARDSON = 5,         ///< Richardson iteration
  KRYLOV_RICHARDSON_CYCLES = 6,  ///< Richardson iteration with cycles support
  KRYLOV_GMRES = 7,              ///< GMRES iterative algorithm
  KRYLOV_GMRES_CYCLES = 8,       ///< GMRES with Cycles support
  KRYLOV_BICGSTAB = 9,           ///< BiCGStab iterative algorithm
  KRYLOV_BICGSTAB_CYCLES = 10,   ///< BiCGStab with Cycles support
};

inline std::string
IterativeMethodPETScName(IterativeMethod it_method)
{
  switch (it_method)
  {
    case IterativeMethod::NONE:
      return "preonly";
    case IterativeMethod::KRYLOV_RICHARDSON:
    case IterativeMethod::KRYLOV_RICHARDSON_CYCLES:
      return "richardson";
    case IterativeMethod::GMRES:
    case IterativeMethod::GMRES_CYCLES:
    case IterativeMethod::KRYLOV_GMRES:
    case IterativeMethod::KRYLOV_GMRES_CYCLES:
      return "gmres";
    case IterativeMethod::KRYLOV_BICGSTAB:
    case IterativeMethod::KRYLOV_BICGSTAB_CYCLES:
      return "bcgs";
    case IterativeMethod::CLASSIC_RICHARDSON:
    case IterativeMethod::CLASSIC_RICHARDSON_CYCLES:
      return "";
  }
  return "";
}

} // namespace opensn
