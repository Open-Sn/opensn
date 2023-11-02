#include "mg_diffusion_solver.h"
#include "chi_runtime.h"
#include "chi_log.h"

//========================================================== Solve 1g problem
void mg_diffusion::Solver::SolveOneGroupProblem(const unsigned int g,
                                                const int64_t verbose)
{
  if (verbose > 1) Chi::log.Log() << "Solving group: " << g;

  KSPSetOperators(petsc_solver_.ksp, A_[g], A_[g]);
  KSPSolve(petsc_solver_.ksp, b_, x_[g]);

  // this is required to compute the inscattering RHS correctly in parallel
  chi_math::PETScUtils::CommunicateGhostEntries(x_[g]);

  if (verbose > 1) Chi::log.Log() << "Done solving group " << g;
}

//  cout << "FLUX ###################################################### FLUX\n";
//  cout << "FLUX ###################################################### FLUX\n";
//  VecView(x[g], PETSC_VIEWER_STDERR_WORLD);

