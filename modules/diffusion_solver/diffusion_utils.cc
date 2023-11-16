#include "diffusion.h"

#include <iomanip>

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/mpi/mpi.h"

namespace opensn
{
namespace diffusion
{

PetscErrorCode
KSPMonitorAChiTech(KSP ksp, PetscInt n, PetscReal rnorm, void* monitordestroy)
{

  Vec Rhs;
  KSPGetRhs(ksp, &Rhs);
  double rhs_norm;
  VecNorm(Rhs, NORM_2, &rhs_norm);
  if (rhs_norm < 1.0e-25) rhs_norm = 1.0;

  if (opensn::mpi.location_id == 0)
  {
    const auto ksp_name = "Diffusion";

    std::stringstream buff;
    buff << ksp_name << " iteration " << std::setw(4) << n << " - Residual " << std::scientific
         << std::setprecision(7) << rnorm / rhs_norm << std::endl;

    Chi::log.Log() << buff.str();
  }
  return 0;
}

PetscErrorCode
DiffusionConvergenceTestNPT(
  KSP ksp, PetscInt n, PetscReal rnorm, KSPConvergedReason* convergedReason, void*)
{
  // Compute rhs norm
  Vec Rhs;
  KSPGetRhs(ksp, &Rhs);
  double rhs_norm;
  VecNorm(Rhs, NORM_2, &rhs_norm);
  if (rhs_norm < 1.0e-25) rhs_norm = 1.0;

  // Compute test
  // criterion
  double tol;
  int64_t maxIts;
  KSPGetTolerances(ksp, NULL, &tol, NULL, &maxIts);

  double relative_residual = rnorm / rhs_norm;

  Chi::log.Log() << "Iteration " << n << " Residual " << rnorm / rhs_norm;

  if (relative_residual < tol) *convergedReason = KSP_CONVERGED_RTOL;

  return KSP_CONVERGED_ITERATING;
}

} // namespace diffusion
} // namespace opensn
