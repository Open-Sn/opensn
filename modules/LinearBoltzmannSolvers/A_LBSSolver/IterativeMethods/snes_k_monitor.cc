#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/snes_k_residual_func_context.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/utils/chi_timer.h"

#include <petscsnes.h>
#include <iomanip>

namespace lbs
{

PetscErrorCode
KEigenSNESMonitor(SNES, PetscInt iter, PetscReal rnorm, void* ctx)
{
  auto& residual_context = *(KResidualFunctionContext*)ctx;

  double k_eff = residual_context.k_eff;
  double reactivity = (k_eff - 1.0) / k_eff;

  std::stringstream iter_info;
  iter_info << Chi::program_timer.GetTimeString() << " " << residual_context.solver_name
            << "_NonLinearK_Outer"
            << " Iteration " << std::setw(5) << iter << " Residual " << std::setw(11) << rnorm
            << " k_eff " << std::fixed << std::setw(10) << std::setprecision(7) << k_eff
            << std::setprecision(2) << "  reactivity " << std::setw(10) << reactivity * 1e5;

  Chi::log.Log() << iter_info.str();

  return 0;
}

PetscErrorCode
KEigenKSPMonitor(KSP ksp, PetscInt iter, PetscReal rnorm, void* ctx)
{
  auto& residual_context = *(KResidualFunctionContext*)ctx;

  std::stringstream iter_info;
  iter_info << "      " << Chi::program_timer.GetTimeString() << " " << residual_context.solver_name
            << "_NonLinearK_Inner"
            << " Iteration " << std::setw(5) << iter << " Residual " << std::setw(11) << rnorm;

  Chi::log.Log() << iter_info.str();

  return 0;
}

} // namespace lbs
