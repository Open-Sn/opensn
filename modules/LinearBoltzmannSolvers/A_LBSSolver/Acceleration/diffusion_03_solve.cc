#include "diffusion.h"

#include "math/PETScUtils/petsc_utils.h"
#include "math/SpatialDiscretization/SpatialDiscretization.h"

#include "physics/chi_physics_namespace.h"

#include "chi_runtime.h"
#include "chi_log.h"

// ###################################################################
/**Solves the system and stores the local solution in the vector provide.
 *
 * \param solution Vector in to which the solution will be parsed.
 * \param use_initial_guess bool [Default:False] Flag, when set, will
 *                 use the values of the output solution as initial guess.*/
void lbs::acceleration::DiffusionSolver::Solve(
  std::vector<double>& solution, bool use_initial_guess /*=false*/)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::Solve";
  Vec x;
  VecDuplicate(rhs_, &x);
  VecSet(x, 0.0);

  if (not use_initial_guess) KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);
  else
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);

  KSPSetTolerances(ksp_,
                   options.residual_tolerance,
                   options.residual_tolerance,
                   1.0e50,
                   options.max_iters);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(A_, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    using namespace chi_math::PETScUtils;
    KSPSetConvergenceTest(
      ksp_, &RelativeResidualConvergenceTest, nullptr, nullptr);

    KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr);

    double rhs_norm;
    VecNorm(rhs_, NORM_2, &rhs_norm);
    Chi::log.Log() << "RHS-norm " << rhs_norm;
  }

  if (use_initial_guess)
  {
    double* x_raw;
    VecGetArray(x, &x_raw);
    size_t k = 0;
    for (const auto& value : solution)
      x_raw[k++] = value;
    VecRestoreArray(x, &x_raw);
  }

  //============================================= Solve
  KSPSolve(ksp_, rhs_, x);

  //============================================= Print convergence info
  if (options.verbose)
  {
    double sol_norm;
    VecNorm(x, NORM_2, &sol_norm);
    Chi::log.Log() << "Solution-norm " << sol_norm;

    using namespace chi_physics;
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);

    Chi::log.Log() << "Convergence Reason: "
                   << GetPETScConvergedReasonstring(reason);
  }

  //============================================= Transfer petsc solution to
  //                                              vector
  if (requires_ghosts_)
  {
    chi_math::PETScUtils::CommunicateGhostEntries(x);
    sdm_.LocalizePETScVectorWithGhosts(x, solution, uk_man_);
  }
  else
    sdm_.LocalizePETScVector(x, solution, uk_man_);

  //============================================= Cleanup x
  VecDestroy(&x);
}

// ###################################################################
/**Solves the system and stores the local solution in the vector provide.
 *
 * \param petsc_solution Vector in to which the solution will be parsed.
 * \param use_initial_guess bool [Default:False] Flag, when set, will
 *                 use the values of the output solution as initial guess.*/
void lbs::acceleration::DiffusionSolver::Solve(
  Vec petsc_solution, bool use_initial_guess /*=false*/)
{
  const std::string fname = "lbs::acceleration::DiffusionMIPSolver::Solve";
  Vec x;
  VecDuplicate(rhs_, &x);
  VecSet(x, 0.0);

  if (not use_initial_guess) KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);
  else
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);

  KSPSetTolerances(ksp_,
                   options.residual_tolerance,
                   options.residual_tolerance,
                   1.0e50,
                   options.max_iters);

  if (options.perform_symmetry_check)
  {
    PetscBool symmetry = PETSC_FALSE;
    MatIsSymmetric(A_, 1.0e-6, &symmetry);
    if (symmetry == PETSC_FALSE)
      throw std::logic_error(fname + ":Symmetry check failed");
  }

  if (options.verbose)
  {
    using namespace chi_math::PETScUtils;
    KSPSetConvergenceTest(
      ksp_, &RelativeResidualConvergenceTest, nullptr, nullptr);

    KSPMonitorSet(ksp_, &KSPMonitorRelativeToRHS, nullptr, nullptr);

    double rhs_norm;
    VecNorm(rhs_, NORM_2, &rhs_norm);
    Chi::log.Log() << "RHS-norm " << rhs_norm;
  }

  if (use_initial_guess) { VecCopy(petsc_solution, x); }

  //============================================= Solve
  KSPSolve(ksp_, rhs_, x);

  //============================================= Print convergence info
  if (options.verbose)
  {
    double sol_norm;
    VecNorm(x, NORM_2, &sol_norm);
    Chi::log.Log() << "Solution-norm " << sol_norm;

    using namespace chi_physics;
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);

    Chi::log.Log() << "Convergence Reason: "
                   << GetPETScConvergedReasonstring(reason);
  }

  //============================================= Transfer petsc solution to
  //                                              vector
  VecCopy(x, petsc_solution);

  //============================================= Cleanup x
  VecDestroy(&x);
}