// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/physics/physics_namespace.h"

#include <sstream>

namespace opensn
{

std::string
GetPETScConvergedReasonstring(KSPConvergedReason reason)
{
  std::stringstream ostr;
  switch (reason)
  {
    case KSP_CONVERGED_RTOL_NORMAL:
      ostr << "KSP_CONVERGED_RTOL_NORMAL";
      break;
    case KSP_CONVERGED_ATOL_NORMAL:
      ostr << "KSP_CONVERGED_ATOL_NORMAL";
      break;
    case KSP_CONVERGED_RTOL:
      ostr << "KSP_CONVERGED_RTOL";
      break;
    case KSP_CONVERGED_ATOL:
      ostr << "KSP_CONVERGED_ATOL";
      break;
    case KSP_CONVERGED_ITS:
      ostr << "KSP_CONVERGED_ITS";
      break;
#if PETSC_VERSION_LT(3, 19, 0)
    case KSP_CONVERGED_CG_NEG_CURVE:
      ostr << "KSP_CONVERGED_CG_NEG_CURVE";
      break;
#else
    case KSP_CONVERGED_NEG_CURVE:
      ostr << "KSP_CONVERGED_NEG_CURVE";
      break;
#endif
#if PETSC_VERSION_LT(3, 19, 0)
    case KSP_CONVERGED_CG_CONSTRAINED:
      ostr << "KSP_CONVERGED_CG_CONSTRAINED";
      break;
#endif
    case KSP_CONVERGED_STEP_LENGTH:
      ostr << "KSP_CONVERGED_STEP_LENGTH";
      break;
    case KSP_CONVERGED_HAPPY_BREAKDOWN:
      ostr << "KSP_CONVERGED_HAPPY_BREAKDOWN";
      break;
      /* diverged */
    case KSP_DIVERGED_NULL:
      ostr << "KSP_DIVERGED_NULL";
      break;
    case KSP_DIVERGED_ITS:
      ostr << "KSP_DIVERGED_ITS";
      break;
    case KSP_DIVERGED_DTOL:
      ostr << "KSP_DIVERGED_DTOL";
      break;
    case KSP_DIVERGED_BREAKDOWN:
      ostr << "KSP_DIVERGED_BREAKDOWN";
      break;
    case KSP_DIVERGED_BREAKDOWN_BICG:
      ostr << "KSP_DIVERGED_BREAKDOWN_BICG";
      break;
    case KSP_DIVERGED_NONSYMMETRIC:
      ostr << "KSP_DIVERGED_NONSYMMETRIC";
      break;
    case KSP_DIVERGED_INDEFINITE_PC:
      ostr << "KSP_DIVERGED_INDEFINITE_PC";
      break;
    case KSP_DIVERGED_NANORINF:
      ostr << "KSP_DIVERGED_NANORINF";
      break;
    case KSP_DIVERGED_INDEFINITE_MAT:
      ostr << "KSP_DIVERGED_INDEFINITE_MAT";
      break;

    default:
      ostr << "Unknown convergence reason.";
  }

  return ostr.str();
}

} // namespace opensn
