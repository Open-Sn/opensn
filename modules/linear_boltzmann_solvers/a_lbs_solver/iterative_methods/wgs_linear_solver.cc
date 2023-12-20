#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/wgs_linear_solver.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/iterative_methods/wgs_convergence_test.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/math/linear_solver/linear_matrix_action_Ax.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/utils/timer.h"

#include <petscksp.h>
#include <memory>
#include <iomanip>

#define sc_double static_cast<double>
#define sc_int64_t static_cast<int64_t>

#define GetGSContextPtr(x) std::dynamic_pointer_cast<WGSContext>(x)

namespace opensn
{
namespace lbs
{

WGSLinearSolver::WGSLinearSolver(WGSContextPtr gs_context_ptr)
  : LinearSolver(IterativeMethodPETScName(gs_context_ptr->groupset_.iterative_method_),
                 gs_context_ptr)
{
  auto& groupset = gs_context_ptr->groupset_;
  auto& solver_tol_options = this->ToleranceOptions();
  solver_tol_options.residual_absolute = groupset.residual_tolerance_;
  solver_tol_options.maximum_iterations = groupset.max_iterations_;
  solver_tol_options.gmres_restart_interval = groupset.gmres_restart_intvl_;
}

WGSLinearSolver::~WGSLinearSolver()
{
  MatDestroy(&A_);
}

void
WGSLinearSolver::PreSetupCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

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
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);
  const auto sizes = gs_context_ptr->SystemSize();

  num_local_dofs_ = sizes.first;
  num_global_dofs_ = sizes.second;
}

void
WGSLinearSolver::SetSystem()
{
  if (IsSystemSet()) return;

  x_ = CreateVector(sc_int64_t(num_local_dofs_), sc_int64_t(num_global_dofs_));

  VecSet(x_, 0.0);
  VecDuplicate(x_, &b_);

  // Create the matrix-shell
  MatCreateShell(PETSC_COMM_WORLD,
                 sc_int64_t(num_local_dofs_),
                 sc_int64_t(num_local_dofs_),
                 sc_int64_t(num_global_dofs_),
                 sc_int64_t(num_global_dofs_),
                 &(*context_ptr_),
                 &A_);

  // Set the action-operator
  MatShellSetOperation(A_, MATOP_MULT, (void (*)())LinearSolverMatrixAction);

  // Set solver operators
  KSPSetOperators(ksp_, A_, A_);
  KSPSetUp(ksp_);
}

void
WGSLinearSolver::SetPreconditioner()
{
  if (IsSystemSet()) return;
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->SetPreconditioner(ksp_);
}

void
WGSLinearSolver::PostSetupCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PostSetupCallback();
}

void
WGSLinearSolver::PreSolveCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PreSolveCallback();
}

void
WGSLinearSolver::SetInitialGuess()
{
  // If the norm of the initial guess is large enough, the initial guess will be used, otherwise it
  // is assumed to be zero.

  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  lbs_solver.SetGSPETScVecFromPrimarySTLvector(groupset, x_, PhiSTLOption::PHI_OLD);

  double init_guess_norm = 0.0;
  VecNorm(x_, NORM_2, &init_guess_norm);

  if (init_guess_norm > 1.0e-10)
  {
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
    if (gs_context_ptr->log_info_) Chi::log.Log() << "Using phi_old as initial guess.";
  }
}

void
WGSLinearSolver::SetRHS()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  if (gs_context_ptr->log_info_)
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " Computing b";

  // SetSource for RHS
  saved_q_moments_local_ = lbs_solver.QMomentsLocal();

  const bool single_richardson =
    iterative_method_ == "richardson" and tolerance_options_.maximum_iterations == 1;

  if (not single_richardson)
  {
    const auto scope = gs_context_ptr->rhs_src_scope_ | ZERO_INCOMING_DELAYED_PSI;
    gs_context_ptr->set_source_function_(
      groupset, lbs_solver.QMomentsLocal(), lbs_solver.PhiOldLocal(), scope);

    // Apply transport operator
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    // Assemble PETSc vector
    lbs_solver.SetGSPETScVecFromPrimarySTLvector(groupset, b_, PhiSTLOption::PHI_NEW);

    // Compute RHS norm
    VecNorm(b_, NORM_2, &context_ptr_->rhs_norm);

    // Compute precondition RHS norm
    PC pc;
    KSPGetPC(ksp_, &pc);
    Vec temp_vec;
    VecDuplicate(b_, &temp_vec);
    PCApply(pc, b_, temp_vec);
    VecNorm(temp_vec, NORM_2, &context_ptr_->rhs_preconditioned_norm);
    VecDestroy(&temp_vec);
  }
  // If we have a single richardson iteration then the user probably wants
  // only a single sweep. Therefore, we are going to combine the scattering
  // source (normally included in the lhs_src_scope) into the sweep for the
  // RHS, and just suppress the kspsolve part.
  else
  {
    const auto scope = gs_context_ptr->rhs_src_scope_ | gs_context_ptr->lhs_src_scope_;
    gs_context_ptr->set_source_function_(
      groupset, lbs_solver.QMomentsLocal(), lbs_solver.PhiOldLocal(), scope);

    // Apply transport operator
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    // Assemble PETSc vector
    lbs_solver.SetGSPETScVecFromPrimarySTLvector(groupset, x_, PhiSTLOption::PHI_NEW);

    // Compute RHS norm
    VecNorm(x_, NORM_2, &context_ptr_->rhs_norm);

    // Compute precondition RHS norm
    PC pc;
    KSPGetPC(ksp_, &pc);
    Vec temp_vec;
    VecDuplicate(x_, &temp_vec);
    PCApply(pc, x_, temp_vec);
    VecNorm(temp_vec, NORM_2, &context_ptr_->rhs_preconditioned_norm);
    VecDestroy(&temp_vec);

    SetKSPSolveSuppressionFlag(true);
  }
}

void
WGSLinearSolver::PostSolveCallback()
{
  // We simply restore the q_moments_local vector.

  // Get convergence reason
  if (not GetKSPSolveSuppressionFlag())
  {
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_, &reason);
    if (reason != KSP_CONVERGED_RTOL and reason != KSP_DIVERGED_ITS)
      Chi::log.Log0Warning() << "Krylov solver failed. "
                             << "Reason: " << GetPETScConvergedReasonstring(reason);
  }

  // Copy x to local solution
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(groupset, x_, PhiSTLOption::PHI_NEW);
  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(groupset, x_, PhiSTLOption::PHI_OLD);

  // Restore saved q_moms
  lbs_solver.QMomentsLocal() = saved_q_moments_local_;

  // Context specific callback
  gs_context_ptr->PostSolveCallback();
}

} // namespace lbs
} // namespace opensn
