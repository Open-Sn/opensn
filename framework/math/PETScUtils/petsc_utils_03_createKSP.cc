#include "petsc_utils.h"

#include <iomanip>

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Creates a common Krylov-solver setup.
 *
This is a macro for:
\code
PETScSolverSetup setup;

KSPCreate(PETSC_COMM_WORLD,&setup.ksp);
KSPSetOperators(setup.ksp,ref_matrix,ref_matrix);
KSPSetType(setup.ksp,in_solver_type.c_str());

KSPSetOptionsPrefix(setup.ksp,in_solver_name.c_str());

KSPGetPC(setup.ksp,&setup.pc);
PCSetType(setup.pc,in_preconditioner_type.c_str());

KSPSetTolerances(setup.ksp,1.e-50,
                 in_relative_residual_tolerance,1.0e50,
                 in_maximum_iterations);
KSPSetInitialGuessNonzero(setup.ksp,PETSC_TRUE);

KSPMonitorSet(setup.ksp,&KSPMonitorRelativeToRHS,NULL,NULL);
KSPSetConvergenceTest(setup.ksp,&RelativeResidualConvergenceTest,NULL,NULL);

return setup;
\endcode*/
chi_math::PETScUtils::PETScSolverSetup
chi_math::PETScUtils::CreateCommonKrylovSolverSetup(
  Mat ref_matrix,
  const std::string& in_solver_name,
  const std::string& in_solver_type,
  const std::string& in_preconditioner_type,
  double in_relative_residual_tolerance,
  int64_t in_maximum_iterations)
{
  PETScSolverSetup setup;

  KSPCreate(PETSC_COMM_WORLD,&setup.ksp);
  KSPSetOperators(setup.ksp,ref_matrix,ref_matrix);
  KSPSetType(setup.ksp,in_solver_type.c_str());

  KSPSetOptionsPrefix(setup.ksp,in_solver_name.c_str());

  KSPGetPC(setup.ksp,&setup.pc);
  PCSetType(setup.pc,in_preconditioner_type.c_str());

  KSPSetTolerances(setup.ksp,1.e-50,
                   in_relative_residual_tolerance,1.0e50,
                   in_maximum_iterations);
  KSPSetInitialGuessNonzero(setup.ksp,PETSC_TRUE);

  KSPSetConvergenceTest(setup.ksp,&RelativeResidualConvergenceTest,
                        nullptr, nullptr);
  KSPSetFromOptions(setup.ksp);

  KSPMonitorSet(setup.ksp, &chi_math::PETScUtils::KSPMonitorRelativeToRHS,
                nullptr, nullptr);

  return setup;
}

//###################################################################
/**Relative Residual Convergence test. The test uses
 * the L2-norm of the residual (\f$ ||b-Ax||_2 \f$), divided by
 * by the L2-norm of the right hand side (\f$ ||b||_2 \f$) compared
 * to a tolerance, \f$ \epsilon \f$.

\f[
 \frac{||b-Ax||_2}{||b||_2} < \epsilon
\f]
 *
 * */
PetscErrorCode
chi_math::PETScUtils::RelativeResidualConvergenceTest(
  KSP ksp, PetscInt,
  PetscReal rnorm,
  KSPConvergedReason* convergedReason,
  void* )
{
  //============================== Compute rhs norm
  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-12)
    rhs_norm = 1.0;

  //============================== Compute test criterion
  double tol;
  int64_t maxIts;
  KSPGetTolerances(ksp, nullptr, &tol, nullptr, &maxIts);

  double relative_residual = rnorm/rhs_norm;

  if (relative_residual < tol)
    *convergedReason = KSP_CONVERGED_RTOL;

  return KSP_CONVERGED_ITERATING;
}

//###################################################################
/**General monitor that print the residual norm relative to the
 * right-hand side norm.*/
PetscErrorCode chi_math::PETScUtils::KSPMonitorRelativeToRHS(
  KSP ksp, PetscInt n, PetscReal rnorm, void*)
{
  Vec Rhs;
  KSPGetRhs(ksp,&Rhs);
  double rhs_norm;
  VecNorm(Rhs,NORM_2,&rhs_norm);
  if (rhs_norm < 1.0e-12)
    rhs_norm = 1.0;

  //Get solver name
  const char* ksp_name;
  KSPGetOptionsPrefix(ksp,&ksp_name);

  //Default to this if ksp_name is NULL
  const char NONAME_SOLVER[] = "NoName-Solver\0";

  if (ksp_name == nullptr)
    ksp_name = NONAME_SOLVER;

  //Print message
  std::stringstream buff;
  buff
    << ksp_name
    << " iteration "
    << std::setw(4) << n
    << " - Residual "
    << std::scientific << std::setprecision(7) << rnorm / rhs_norm
    << std::endl;

  Chi::log.Log() << buff.str();

  return 0;
}

//###################################################################
/**General monitor that print the residual norm relative to the
 * right-hand side norm.*/
PetscErrorCode chi_math::PETScUtils::KSPMonitorStraight(
  KSP ksp, PetscInt n, PetscReal rnorm, void*)
{
  //Get solver name
  const char* ksp_name;
  KSPGetOptionsPrefix(ksp,&ksp_name);

  //Default to this if ksp_name is NULL
  const char NONAME_SOLVER[] = "NoName-Solver\0";

  if (ksp_name == nullptr)
    ksp_name = NONAME_SOLVER;

  //Print message
  std::stringstream buff;
  buff
    << ksp_name
    << " iteration "
    << std::setw(4) << n
    << " - Residual "
    << std::scientific << std::setprecision(7) << rnorm
    << std::endl;

  Chi::log.Log() << buff.str();

  return 0;
}
