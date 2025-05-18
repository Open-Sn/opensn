// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/slepc_keigen_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/ags_linear_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_problem.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_compute.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_vecops.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/object_factory.h"
#include "framework/runtime.h"
#include <slepceps.h>

namespace opensn
{

namespace
{

// Shell matrix multiplication: y = A * x
PetscErrorCode
ShellMult(Mat M, Vec x, Vec y)
{
  void* raw_ctx = nullptr;
  MatShellGetContext(M, &raw_ctx);

  auto* ctx = static_cast<SLEPcLinearKEigenContext*>(raw_ctx);
  auto& lbs_problem = ctx->lbs_problem;
  auto& phi_old_local = lbs_problem->GetPhiOldLocal();
  auto& q_moments_local = lbs_problem->GetQMomentsLocal();
  auto& groupset_ids = ctx->groupset_ids;
  auto active_set_source_function = lbs_problem->GetActiveSetSourceFunction();
  auto ags_solver = lbs_problem->GetAGSSolver();

  LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(
    *lbs_problem, groupset_ids, x, PhiSTLOption::PHI_OLD);

  Set(q_moments_local, 0.0);
  for (auto& groupset : lbs_problem->GetGroupsets())
  {
    active_set_source_function(groupset,
                               q_moments_local,
                               phi_old_local,
                               APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  }

  ags_solver->Solve();

  LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector(
    *lbs_problem, groupset_ids, y, PhiSTLOption::PHI_NEW);

  return 0;
}

} // namespace

static PetscErrorCode
SLEPcLinearKEigenMonitor(EPS eps,
                         PetscInt its,
                         PetscInt nconv,
                         PetscScalar* eigr,
                         PetscScalar* eigi,
                         PetscReal* errest,
                         PetscInt nest,
                         void* ctx)
{
  std::stringstream iter_info;
  iter_info << program_timer.GetTimeString() << " SLEPc EPS Iteration " << std::setw(5) << its
            << " : ";
  if (nconv > 0)
    iter_info << " Residual " << std::setw(11) << errest[0];
  else
    iter_info << " No converged eigenpairs";
  log.Log() << iter_info.str();

  return 0;
}

void
SLEPcLinearKEigenSolver::PreSetupCallback()
{
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  ctx->groupset_ids.clear();
  for (auto& gs : ctx->lbs_problem->GetGroupsets())
    ctx->groupset_ids.push_back(gs.id);
}

void
SLEPcLinearKEigenSolver::SetMonitor()
{
}

void
SLEPcLinearKEigenSolver::SetSystemSize()
{
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  auto sizes = ctx->lbs_problem->GetNumPhiIterativeUnknowns();
  num_local_dofs_ = static_cast<int64_t>(sizes.first);
  num_global_dofs_ = static_cast<int64_t>(sizes.second);
}

void
SLEPcLinearKEigenSolver::SetSystem()
{
  // Create shell matrix A
  MatCreateShell(PETSC_COMM_WORLD,
                 (PetscInt)num_local_dofs_,
                 (PetscInt)num_local_dofs_,
                 (PetscInt)num_global_dofs_,
                 (PetscInt)num_global_dofs_,
                 context_ptr_.get(),
                 &A_);
  MatShellSetOperation(A_, MATOP_MULT, (void (*)(void))ShellMult);

  // Create x
  x_ = CreateVector(num_local_dofs_, num_global_dofs_);
}

void
SLEPcLinearKEigenSolver::SetInitialGuess()
{
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector(
    *ctx->lbs_problem, ctx->groupset_ids, x_, PhiSTLOption::PHI_OLD);
}

void
SLEPcLinearKEigenSolver::Solve()
{
  PreSolveCallback();

  EPSSetFromOptions(eps_);
  EPSSetOperators(eps_, A_, nullptr);
  EPSSetProblemType(eps_, EPS_NHEP);
  EPSSetWhichEigenpairs(eps_, EPS_LARGEST_MAGNITUDE);
  EPSSetTolerances(eps_, tolerance_options.residual_absolute, tolerance_options.maximum_iterations);
  EPSSetType(eps_, eps_type_.c_str());
  EPSSetInitialSpace(eps_, 1, &x_);
  EPSMonitorSet(eps_, SLEPcLinearKEigenMonitor, nullptr, nullptr);
  EPSSolve(eps_);

  // Check for convergence
  PetscInt nconv = 0;
  EPSGetConverged(eps_, &nconv);
  if (nconv == 0)
  {
    log.Log0Warning() << program_timer.GetTimeString() << " SLEPc EPS : No eigenpairs converged"
                      << std::endl;
    return;
  }

  PostSolveCallback();
}

void
SLEPcLinearKEigenSolver::PreSolveCallback()
{
  log.Log() << "Executing SLPEc Eigenvalue Problem Solver with iterative method "
            << GetIterativeMethodName() << std::endl;

  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  auto& lbs_problem = ctx->lbs_problem;

  for (auto& wgs_solver : lbs_problem->GetWGSSolvers())
  {
    auto context = wgs_solver->GetContext();
    auto wgs_context = std::dynamic_pointer_cast<WGSContext>(context);
    wgs_context->lhs_src_scope.Unset(APPLY_WGS_FISSION_SOURCES); // lhs_scope
    wgs_context->rhs_src_scope.Unset(APPLY_AGS_FISSION_SOURCES); // rhs_scope
  }
}

void
SLEPcLinearKEigenSolver::PostSolveCallback()
{
  auto ctx = std::dynamic_pointer_cast<SLEPcLinearKEigenContext>(context_ptr_);
  auto& lbs_problem = ctx->lbs_problem;
  auto& groupset_ids = ctx->groupset_ids;
  auto& phi_old_local = lbs_problem->GetPhiOldLocal();
  auto& phi_new_local = lbs_problem->GetPhiNewLocal();

  // Because the eigenvalue returned by EPS is internally scaled and does not correspond
  // to the k-eigenvalue of the transport operator, we extract the converged eigenvector
  // and compute k as the ratio of fission rates.

  // Get current eigenvector as phi_old and compute fission rate
  PetscScalar kr, ki;
  Vec x = CreateVector(num_local_dofs_, num_global_dofs_);
  EPSGetEigenpair(eps_, 0, &kr, &ki, x, NULL);
  LBSVecOps::SetPrimarySTLvectorFromMultiGSPETScVec(
    *lbs_problem, groupset_ids, x, PhiSTLOption::PHI_OLD);
  double F_prev = ComputeFissionProduction(*lbs_problem, phi_old_local);

  // Perform a transport solve, get the resulting eigenvector as phi_new, and compute fission rate
  Vec y = CreateVector(num_local_dofs_, num_global_dofs_);
  ShellMult(A_, x, y);
  LBSVecOps::SetMultiGSPETScVecFromPrimarySTLvector(
    *lbs_problem, ctx->groupset_ids, y, PhiSTLOption::PHI_NEW);
  double F = ComputeFissionProduction(*lbs_problem, phi_new_local);

  // Compute k-eigenvalue
  ctx->eigenvalue = F / F_prev;
  log.Log() << "Final k-eigenvalue: " << std::setprecision(7) << ctx->eigenvalue << "\n";

  VecDestroy(&x);
  VecDestroy(&y);
}

} // namespace opensn
