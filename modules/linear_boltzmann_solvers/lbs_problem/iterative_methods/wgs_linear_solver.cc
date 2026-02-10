// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_convergence_test.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/sweep_wgs_context.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <petscksp.h>
#include "caliper/cali.h"
#include <memory>
#include <iomanip>

namespace opensn
{

WGSLinearSolver::WGSLinearSolver(const std::shared_ptr<WGSContext>& gs_context_ptr)
  : PETScLinearSolver(gs_context_ptr->groupset.iterative_method, gs_context_ptr)
{
  auto& groupset = gs_context_ptr->groupset;
  auto& solver_tol_options = this->GetToleranceOptions();
  solver_tol_options.residual_absolute = groupset.residual_tolerance;
  solver_tol_options.maximum_iterations = static_cast<PetscInt>(groupset.max_iterations);
  solver_tol_options.gmres_restart_interval = static_cast<PetscInt>(groupset.gmres_restart_intvl);
}

WGSLinearSolver::~WGSLinearSolver()
{
  MatDestroy(&A_);
}

void
WGSLinearSolver::PreSetupCallback()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->PreSetupCallback();
}

void
WGSLinearSolver::SetConvergenceTest()
{
  KSPSetConvergenceTest(ksp_, &GSConvergenceTest, nullptr, nullptr);
}

void
WGSLinearSolver::SetSystemSize()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  const auto sizes = gs_context_ptr->GetSystemSize();
  num_local_dofs_ = sizes.first;
  num_global_dofs_ = sizes.second;
}

void
WGSLinearSolver::SetSystem()
{
  if (IsSystemSet())
    return;

  x_ = CreateVector(num_local_dofs_, num_global_dofs_);

  VecSet(x_, 0.0);
  VecDuplicate(x_, &b_);

  // Create the matrix-shell
  MatCreateShell(opensn::mpi_comm,
                 num_local_dofs_,
                 num_local_dofs_,
                 num_global_dofs_,
                 num_global_dofs_,
                 &(*context_ptr_),
                 &A_);

  // Set the action-operator
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast,modernize-redundant-void-arg)
  MatShellSetOperation(A_, MATOP_MULT, (void (*)(void))LinearSolverMatrixAction);

  // Set solver operators
  KSPSetOperators(ksp_, A_, A_);
  KSPSetUp(ksp_);
}

void
WGSLinearSolver::SetPreconditioner()
{
  if (IsSystemSet())
    return;
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->SetPreconditioner(ksp_);
}

void
WGSLinearSolver::PostSetupCallback()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  gs_context_ptr->PostSetupCallback();
}

void
WGSLinearSolver::PreSolveCallback()
{
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;
  if (do_problem.GetOptions().verbose_inner_iterations)
  {
    log.Log() << "Solving groupset " << groupset.id << " with " << this->GetIterativeMethodName()
              << " (groups " << groupset.first_group << "-" << groupset.last_group << ", "
              << groupset.quadrature->abscissae.size() << " angles)\n";
  }
  gs_context_ptr->PreSolveCallback();
}

void
WGSLinearSolver::SetInitialGuess()
{
  // If the norm of the initial guess is large enough, the initial guess will be used, otherwise it
  // is assumed to be zero.

  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;

  LBSVecOps::SetGSPETScVecFromPrimarySTLvector(do_problem, groupset, x_, PhiSTLOption::PHI_OLD);

  double init_guess_norm = 0.0;
  VecNorm(x_, NORM_2, &init_guess_norm);

  if (init_guess_norm > 1.0e-10)
  {
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
    if (gs_context_ptr->log_info)
      log.Log() << "Using phi_old as initial guess.";
  }
}

void
WGSLinearSolver::SetRHS()
{
  CALI_CXX_MARK_SCOPE("WGSLinearSolver::SetRHS");

  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;

  if (gs_context_ptr->log_info)
    log.Log() << program_timer.GetTimeString() << " Computing b";

  // SetSource for RHS
  saved_q_moments_local_ = do_problem.GetQMomentsLocal();

  const bool single_richardson =
    groupset.iterative_method == LinearSystemSolver::IterativeMethod::PETSC_RICHARDSON and
    tolerance_options.maximum_iterations == 1;

  if (not single_richardson)
  {
    const auto scope = gs_context_ptr->rhs_src_scope | ZERO_INCOMING_DELAYED_PSI;
    gs_context_ptr->set_source_function(
      groupset, do_problem.GetQMomentsLocal(), do_problem.GetPhiOldLocal(), scope);

    // Enable RHS time (tau*psi^n)
    auto sweep_ctx = std::dynamic_pointer_cast<SweepWGSContext>(gs_context_ptr);
    if (sweep_ctx && sweep_ctx->sweep_chunk->IsTimeDependent())
      sweep_ctx->sweep_chunk->IncludeRHSTimeTerm(true);

    // Apply transport operator
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    // Assemble PETSc vector
    LBSVecOps::SetGSPETScVecFromPrimarySTLvector(do_problem, groupset, b_, PhiSTLOption::PHI_NEW);

    // Compute RHS norm
    VecNorm(b_, NORM_2, &gs_context_ptr->rhs_norm);

    // Compute precondition RHS norm
    PC pc = nullptr;
    KSPGetPC(ksp_, &pc);
    Vec temp_vec = nullptr;
    VecDuplicate(b_, &temp_vec);
    PCApply(pc, b_, temp_vec);
    VecNorm(temp_vec, NORM_2, &gs_context_ptr->rhs_preconditioned_norm);
    VecDestroy(&temp_vec);
  }
  // If we have a single richardson iteration then the user probably wants
  // only a single sweep. Therefore, we are going to combine the scattering
  // source (normally included in the lhs_src_scope) into the sweep for the
  // RHS, and just suppress the kspsolve part.
  else
  {
    const auto scope = gs_context_ptr->rhs_src_scope | gs_context_ptr->lhs_src_scope;
    gs_context_ptr->set_source_function(
      groupset, do_problem.GetQMomentsLocal(), do_problem.GetPhiOldLocal(), scope);

    // Apply transport operator
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    // Assemble PETSc vector
    LBSVecOps::SetGSPETScVecFromPrimarySTLvector(do_problem, groupset, x_, PhiSTLOption::PHI_NEW);

    // Compute RHS norm
    VecNorm(x_, NORM_2, &gs_context_ptr->rhs_norm);

    // Compute precondition RHS norm
    PC pc = nullptr;
    KSPGetPC(ksp_, &pc);
    Vec temp_vec = nullptr;
    VecDuplicate(x_, &temp_vec);
    PCApply(pc, x_, temp_vec);
    VecNorm(temp_vec, NORM_2, &gs_context_ptr->rhs_preconditioned_norm);
    VecDestroy(&temp_vec);

    SetKSPSolveSuppressionFlag(true);
  }
}

void
WGSLinearSolver::PostSolveCallback()
{
  // Get convergence reason
  if (not GetKSPSolveSuppressionFlag())
  {
    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    KSPGetConvergedReason(ksp_, &reason);
    PetscInt its = 0;
    KSPGetIterationNumber(ksp_, &its);
    if (reason < 0)
      log.Log0Warning() << "Krylov solver diverged. "
                        << "Reason: " << GetPETScConvergedReasonstring(reason);
    else if (reason == KSP_CONVERGED_RTOL)
    {
      auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
      if (gs_context_ptr && gs_context_ptr->log_info && its == 0)
        log.Log() << program_timer.GetTimeString() << " CONVERGED (relative tolerance)";
    }
    else if (reason == KSP_CONVERGED_ATOL)
    {
      auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
      if (gs_context_ptr && gs_context_ptr->log_info && its == 0)
        log.Log() << program_timer.GetTimeString() << " CONVERGED (absolute tolerance)";
    }
    else if (reason == KSP_DIVERGED_ITS)
      log.Log0Warning() << "Krylov solver reached iteration limit.";
  }

  // Copy x to local solution
  auto gs_context_ptr = std::dynamic_pointer_cast<WGSContext>(context_ptr_);
  auto& groupset = gs_context_ptr->groupset;
  auto& do_problem = gs_context_ptr->do_problem;
  LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(do_problem, groupset, x_, PhiSTLOption::PHI_NEW);
  LBSVecOps::SetPrimarySTLvectorFromGSPETScVec(do_problem, groupset, x_, PhiSTLOption::PHI_OLD);

  // Restore saved q_moms
  do_problem.GetQMomentsLocal() = saved_q_moments_local_;

  // Context specific callback
  gs_context_ptr->PostSolveCallback();
}

} // namespace opensn
