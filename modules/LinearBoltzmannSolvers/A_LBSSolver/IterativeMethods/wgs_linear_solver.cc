#include "wgs_linear_solver.h"

#include "wgs_convergence_test.h"
#include "A_LBSSolver/lbs_solver.h"

#include "math/PETScUtils/petsc_utils.h"
#include "math/LinearSolver/linear_matrix_action_Ax.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

#include <petscksp.h>
#include <memory>
#include <iomanip>

#define sc_double static_cast<double>
#define sc_int64_t static_cast<int64_t>

#define GetGSContextPtr(x)                                                     \
  std::dynamic_pointer_cast<WGSContext<Mat, Vec, KSP>>(x)

namespace lbs
{

template <>
void WGSLinearSolver<Mat, Vec, KSP>::PreSetupCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PreSetupCallback();
}

template <>
void WGSLinearSolver<Mat, Vec, KSP>::SetConvergenceTest()
{
  KSPSetConvergenceTest(solver_, &GSConvergenceTest, nullptr, nullptr);
}

template <>
void WGSLinearSolver<Mat, Vec, KSP>::SetSystemSize()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);
  const auto sizes = gs_context_ptr->SystemSize();

  num_local_dofs_ = sizes.first;
  num_globl_dofs_ = sizes.second;
}

template <>
void WGSLinearSolver<Mat, Vec, KSP>::SetSystem()
{
  if (IsSystemSet()) return;

  x_ = chi_math::PETScUtils::CreateVector(sc_int64_t(num_local_dofs_),
                                          sc_int64_t(num_globl_dofs_));

  VecSet(x_, 0.0);
  VecDuplicate(x_, &b_);

  //============================================= Create the matrix-shell
  MatCreateShell(PETSC_COMM_WORLD,
                 sc_int64_t(num_local_dofs_),
                 sc_int64_t(num_local_dofs_),
                 sc_int64_t(num_globl_dofs_),
                 sc_int64_t(num_globl_dofs_),
                 &(*context_ptr_),
                 &A_);

  //============================================= Set the action-operator
  MatShellSetOperation(
    A_, MATOP_MULT, (void (*)())chi_math::LinearSolverMatrixAction<Mat, Vec>);

  //============================================= Set solver operators
  KSPSetOperators(solver_, A_, A_);
  KSPSetUp(solver_);
}

template <>
void WGSLinearSolver<Mat, Vec, KSP>::SetPreconditioner()
{
  if (IsSystemSet()) return;
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->SetPreconditioner(solver_);
}

template <>
void WGSLinearSolver<Mat, Vec, KSP>::PostSetupCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PostSetupCallback();
}

template <>
void WGSLinearSolver<Mat, Vec, KSP>::PreSolveCallback()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  gs_context_ptr->PreSolveCallback();
}

/**Sets the initial guess for a gs solver. If the initial guess's norm
 * is large enough the initial guess will be used, otherwise it is assumed
 * zero.*/
template <>
void WGSLinearSolver<Mat, Vec, KSP>::SetInitialGuess()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  lbs_solver.SetGSPETScVecFromPrimarySTLvector(
    groupset, x_, PhiSTLOption::PHI_OLD);

  double init_guess_norm = 0.0;
  VecNorm(x_, NORM_2, &init_guess_norm);

  if (init_guess_norm > 1.0e-10)
  {
    KSPSetInitialGuessNonzero(solver_, PETSC_TRUE);
    if (gs_context_ptr->log_info_)
      Chi::log.Log() << "Using phi_old as initial guess.";
  }
}

template <>
void WGSLinearSolver<Mat, Vec, KSP>::SetRHS()
{
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  if (gs_context_ptr->log_info_)
    Chi::log.Log() << Chi::program_timer.GetTimeString() << " Computing b";

  // SetSource for RHS
  saved_q_moments_local_ = lbs_solver.QMomentsLocal();

  const bool single_richardson = iterative_method_ == "richardson" and
                                 tolerance_options_.maximum_iterations == 1;

  if (not single_richardson)
  {
    const int scope =
      gs_context_ptr->rhs_src_scope_ | ZERO_INCOMING_DELAYED_PSI;
    gs_context_ptr->set_source_function_(
      groupset, lbs_solver.QMomentsLocal(), lbs_solver.PhiOldLocal(), scope);

    //=================================================== Apply transport
    //operator
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    //=================================================== Assemble PETSc vector
    lbs_solver.SetGSPETScVecFromPrimarySTLvector(
      groupset, b_, PhiSTLOption::PHI_NEW);

    //============================================= Compute RHS norm
    VecNorm(b_, NORM_2, &context_ptr_->rhs_norm);

    //============================================= Compute precondition RHS
    //norm
    PC pc;
    KSPGetPC(solver_, &pc);
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
    const int scope = gs_context_ptr->rhs_src_scope_ |
                      gs_context_ptr->lhs_src_scope_;
    gs_context_ptr->set_source_function_(
      groupset, lbs_solver.QMomentsLocal(), lbs_solver.PhiOldLocal(), scope);

    //=================================================== Apply transport
    //operator
    gs_context_ptr->ApplyInverseTransportOperator(scope);

    //=================================================== Assemble PETSc vector
    lbs_solver.SetGSPETScVecFromPrimarySTLvector(
      groupset, x_, PhiSTLOption::PHI_NEW);

    //============================================= Compute RHS norm
    VecNorm(x_, NORM_2, &context_ptr_->rhs_norm);

    //============================================= Compute precondition RHS
    //norm
    PC pc;
    KSPGetPC(solver_, &pc);
    Vec temp_vec;
    VecDuplicate(x_, &temp_vec);
    PCApply(pc, x_, temp_vec);
    VecNorm(temp_vec, NORM_2, &context_ptr_->rhs_preconditioned_norm);
    VecDestroy(&temp_vec);

    SetKSPSolveSuppressionFlag(true);
  }
}

/**For this callback we simply restore the q_moments_local vector.*/
template <>
void WGSLinearSolver<Mat, Vec, KSP>::PostSolveCallback()
{
  //============================================= Get convergence reason
  if (not GetKSPSolveSuppressionFlag())
  {
    KSPConvergedReason reason;
    KSPGetConvergedReason(solver_, &reason);
    if (reason != KSP_CONVERGED_RTOL and reason != KSP_DIVERGED_ITS)
      Chi::log.Log0Warning() << "Krylov solver failed. "
                             << "Reason: "
                             << chi_physics::GetPETScConvergedReasonstring(
                                  reason);
  }

  //============================================= Copy x to local solution
  auto gs_context_ptr = GetGSContextPtr(context_ptr_);

  auto& groupset = gs_context_ptr->groupset_;
  auto& lbs_solver = gs_context_ptr->lbs_solver_;

  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(
    groupset, x_, PhiSTLOption::PHI_NEW);
  lbs_solver.SetPrimarySTLvectorFromGSPETScVec(
    groupset, x_, PhiSTLOption::PHI_OLD);

  //============================================= Restore saved q_moms
  lbs_solver.QMomentsLocal() = saved_q_moments_local_;

  //============================================= Context specific callback
  gs_context_ptr->PostSolveCallback();
}

template <>
WGSLinearSolver<Mat, Vec, KSP>::~WGSLinearSolver()
{
  MatDestroy(&A_);
}
} // namespace lbs
