// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/snes_k_monitor.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/iteration_logging.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/snes_k_residual_func_context.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <petscsnes.h>

namespace opensn
{

PetscErrorCode
KEigenSNESMonitor(SNES snes, PetscInt iter, PetscReal rnorm, void* ctx)
{
  auto& residual_context = *static_cast<KResidualFunctionContext*>(ctx);

  double k_eff = residual_context.k_eff;
  double reactivity = (k_eff - 1.0) / k_eff;

  KSP ksp = nullptr;
  OpenSnPETScCall(SNESGetKSP(snes, &ksp));
  KSPConvergedReason linear_reason = KSP_CONVERGED_ITERATING;
  OpenSnPETScCall(KSPGetConvergedReason(ksp, &linear_reason));
  const auto inner_status = KSPReasonToIterationStatus(linear_reason);

  std::stringstream out;
  out << "NLKE outer iteration = " << iter;
  AppendNumericField(out, "residual", rnorm, Scientific(6));
  AppendNumericField(out, "k_eff", k_eff, Fixed(7));
  AppendNumericField(out, "rho_pcm", reactivity * 1.0e5, Fixed(2));
  if (inner_status != IterationStatus::NONE)
  {
    out << FormatNestedStatusSummary("NLKE inners", inner_status);
    if (log.GetVerbosity() >= 2)
      out << ", detail = " << GetPETScConvergedReasonstring(linear_reason);
  }
  log.Log() << program_timer.GetTimeString() << " " << out.str();

  return 0;
}

PetscErrorCode
KEigenKSPMonitor(KSP /*ksp*/, PetscInt iter, PetscReal rnorm, void* ctx)
{
  auto& residual_context = *static_cast<KResidualFunctionContext*>(ctx);

  std::stringstream iter_info;
  iter_info << program_timer.GetTimeString() << " NLKE inner"
            << " iteration = " << iter;
  AppendNumericField(iter_info, "residual", rnorm, Scientific(6));
  log.Log() << iter_info.str();

  return 0;
}

} // namespace opensn
