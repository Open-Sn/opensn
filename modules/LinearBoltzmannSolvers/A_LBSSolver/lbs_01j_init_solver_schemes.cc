#include "lbs_solver.h"

#include "A_LBSSolver/IterativeMethods/ags_context.h"
#include "A_LBSSolver/IterativeMethods/ags_linear_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

void lbs::LBSSolver::InitializeSolverSchemes()
{
  Chi::log.Log() << "Initializing Solver schemes";

  InitializeWGSSolvers();

  /*This default behavior covers the situation when no Across-GroupSet (AGS)
   * solvers have been created for this solver.*/
  ags_solvers_.clear();
  //=========================================== Default AGS scheme
  if (options_.ags_scheme.empty())
  {
    auto ags_context = std::make_shared<AGSContext<Mat,Vec,KSP>>(
      *this, wgs_solvers_);

    auto ags_solver = std::make_shared<AGSLinearSolver<Mat,Vec,KSP>>(
      "richardson", ags_context,
      groupsets_.front().id_, groupsets_.back().id_);
    ags_solver->ToleranceOptions().maximum_iterations = 1;
    ags_solver->SetVerbosity(options_.verbose_ags_iterations);

    ags_solvers_.push_back(ags_solver);

    primary_ags_solver_ = ags_solvers_.front();
  }
}