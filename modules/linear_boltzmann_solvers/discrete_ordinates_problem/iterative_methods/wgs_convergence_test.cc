// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_convergence_test.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

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
  PetscErrorCode ierr = KSPGetApplicationContext(ksp, static_cast<void*>(&context));
  if (ierr != PETSC_SUCCESS)
    return ierr;

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
  ierr = KSPGetTolerances(ksp, nullptr, &tol, nullptr, &maxIts);
  if (ierr != PETSC_SUCCESS)
    return ierr;

  double scaled_residual = rnorm * residual_scale;

  // Print iteration information
  std::stringstream iter_info;
  iter_info << program_timer.GetTimeString() << " WGS groups [" << context->groupset.first_group
            << "-" << context->groupset.last_group << "]"
            << " iteration = " << n;
  AppendNumericField(iter_info, "residual", scaled_residual, Scientific(6));

  if (scaled_residual < tol)
  {
    *convergedReason = KSP_CONVERGED_ATOL;
    iter_info << ", status = " << IterationStatusName(KSPReasonToIterationStatus(*convergedReason));
  }

  if (context->log_info)
    log.Log() << iter_info.str();

  return PETSC_SUCCESS;
}

} // namespace opensn
