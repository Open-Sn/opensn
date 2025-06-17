// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/solvers/iterative_methods/snes_k_residual_func_context.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <petscsnes.h>
#include <iomanip>

namespace opensn
{

PetscErrorCode
KEigenSNESMonitor(SNES, PetscInt iter, PetscReal rnorm, void* ctx)
{
  auto& residual_context = *(KResidualFunctionContext*)ctx;

  double k_eff = residual_context.k_eff;
  double reactivity = (k_eff - 1.0) / k_eff;

  std::stringstream iter_info;
  iter_info << program_timer.GetTimeString() << " " << residual_context.solver_name
            << "_NonLinearK_Outer"
            << " Iteration " << std::setw(5) << iter << " Residual " << std::setw(11) << rnorm
            << " k_eff " << std::fixed << std::setw(10) << std::setprecision(7) << k_eff
            << std::setprecision(2) << "  reactivity " << std::setw(10) << reactivity * 1e5;

  log.Log() << iter_info.str();

  return 0;
}

PetscErrorCode
KEigenKSPMonitor(KSP ksp, PetscInt iter, PetscReal rnorm, void* ctx)
{
  auto& residual_context = *(KResidualFunctionContext*)ctx;

  std::stringstream iter_info;
  iter_info << "      " << program_timer.GetTimeString() << " " << residual_context.solver_name
            << "_NonLinearK_Inner"
            << " Iteration " << std::setw(5) << iter << " Residual " << std::setw(11) << rnorm;

  log.Log() << iter_info.str();

  return 0;
}

} // namespace opensn
