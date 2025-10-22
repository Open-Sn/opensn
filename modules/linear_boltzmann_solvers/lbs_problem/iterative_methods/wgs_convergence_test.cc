// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_convergence_test.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <iomanip>

namespace opensn
{

/**Customized convergence test.*/
PetscErrorCode
GSConvergenceTest(
  KSP ksp, PetscInt n, PetscReal rnorm, KSPConvergedReason* convergedReason, void* /*unused*/)
{
  CALI_CXX_MARK_FUNCTION;

  // Get data context
  WGSContext* context = nullptr;
  KSPGetApplicationContext(ksp, static_cast<void*>(&context));

  // Set rhs norm
  double residual_scale = 1.0;
  switch (context->residual_scale_type)
  {
    case ResidualScaleType::NONE:
      residual_scale = 1.0;
      break;
    case ResidualScaleType::RHS_NORM:
      if (context->rhs_norm > 1.0e-25)
        residual_scale = 1.0 / context->rhs_norm;
      break;
    case ResidualScaleType::RHS_PRECONDITIONED_NORM:
      if (context->rhs_preconditioned_norm > 1.0e-25)
        residual_scale = 1.0 / context->rhs_preconditioned_norm;
      break;
    case ResidualScaleType::CUSTOM_SCALE:
      if (context->custom_residual_scale > 1.0e-25)
        residual_scale = 1.0 / context->custom_residual_scale;
      break;
  }

  // Compute test criterion
  double tol = 0.0;
  int64_t maxIts = 0;
  KSPGetTolerances(ksp, nullptr, &tol, nullptr, &maxIts);

  double scaled_residual = rnorm * residual_scale;

  // Print iteration information
  std::string offset;
  if (context->groupset.apply_wgdsa or context->groupset.apply_tgdsa)
    offset = std::string("    ");

  std::stringstream iter_info;
  iter_info << program_timer.GetTimeString() << " " << offset << "WGS groups ["
            << context->groupset.groups.front().id << "-" << context->groupset.groups.back().id
            << "]"
            << " Iteration " << std::setw(5) << n << " Residual " << std::setw(9)
            << scaled_residual;

  if (scaled_residual < tol)
  {
    *convergedReason = KSP_CONVERGED_RTOL;
    iter_info << " CONVERGED\n";
  }

  if (context->log_info)
    log.Log() << iter_info.str() << std::endl;

  return KSP_CONVERGED_ITERATING;
}

} // namespace opensn
