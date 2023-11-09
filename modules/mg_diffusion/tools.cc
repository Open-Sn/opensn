#include "modules/mg_diffusion/tools.h"
#include "framework/logging/log.h"
#include "modules/mg_diffusion/mg_diffusion_solver.h"

#include <iomanip>

#include "framework/runtime.h"

PetscErrorCode
mg_diffusion::MGKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void*)
{
  Vec Rhs;
  KSPGetRhs(ksp, &Rhs);
  double rhs_norm;
  VecNorm(Rhs, NORM_2, &rhs_norm);
  if (rhs_norm < 1.0e-12) rhs_norm = 1.0;

  // Get solver name
  const char* ksp_name;
  KSPGetOptionsPrefix(ksp, &ksp_name);

  // Default to this if ksp_name is NULL
  const char NONAME_SOLVER[] = "NoName-Solver\0";

  if (ksp_name == nullptr) ksp_name = NONAME_SOLVER;

  mg_diffusion::KSPAppContext* my_app_context;
  KSPGetApplicationContext(ksp, &my_app_context);

  // Print message
  if (my_app_context->verbose == PETSC_TRUE)
  {
    std::stringstream buff;
    buff << ksp_name << " iteration " << std::setw(4) << n << " - Residual " << std::scientific
         << std::setprecision(7) << rnorm / rhs_norm << std::endl;

    Chi::log.Log() << buff.str();
  }

  return 0;
}
