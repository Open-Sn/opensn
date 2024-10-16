// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/nonlinear_solver/nonlinear_solver.h"
#include "framework/logging/log.h"
#include "framework/logging/stringstream_color.h"

namespace opensn
{

NonLinearSolver::NonLinearSolver(std::shared_ptr<NonLinearSolverContext> context_ptr,
                                 const InputParameters& params)
  : solver_name_(params.GetParamValue<std::string>("name")),
    context_ptr_(context_ptr),
    options_(params)
{
}

NonLinearSolver::~NonLinearSolver()
{
  SNESDestroy(&nl_solver_);
  VecDestroy(&x_);
  VecDestroy(&r_);
  MatDestroy(&J_);
}

void
NonLinearSolver::ApplyToleranceOptions()
{
  SNESSetTolerances(nl_solver_,
                    options_.nl_abs_tol,
                    options_.nl_rel_tol,
                    options_.nl_sol_tol,
                    options_.nl_max_its,
                    options_.nl_max_r_evaluations);
  SNESSetMaxLinearSolveFailures(nl_solver_, options_.l_max_failed_iterations);
  KSP ksp;
  SNESGetKSP(nl_solver_, &ksp);
  KSPSetTolerances(
    ksp, options_.l_rel_tol, options_.l_abs_tol, options_.l_div_tol, options_.l_max_its);
  if (options_.l_method == "gmres")
  {
    KSPGMRESSetRestart(ksp, options_.l_gmres_restart_intvl);
    KSPGMRESSetBreakdownTolerance(ksp, options_.l_gmres_breakdown_tol);
  }
  KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
}

void
NonLinearSolver::PreSetupCallback()
{
}

void
NonLinearSolver::SetOptions()
{
}

void
NonLinearSolver::SetSolverContext()
{
  SNESSetApplicationContext(nl_solver_, &(*context_ptr_));
}

void
NonLinearSolver::SetConvergenceTest()
{
}

void
NonLinearSolver::SetMonitor()
{
}

void
NonLinearSolver::SetPreconditioner()
{
}

void
NonLinearSolver::PostSetupCallback()
{
}

void
NonLinearSolver::Setup()
{
  if (IsSystemSet())
    return;
  this->PreSetupCallback();

  SNESCreate(mpi_comm, &nl_solver_);

  SNESSetOptionsPrefix(nl_solver_, solver_name_.c_str());

  SNESSetType(nl_solver_, options_.petsc_snes_type.c_str());
  if (options_.nl_method == "LINEAR")
    SNESSetType(nl_solver_, SNESKSPONLY);
  SNESLineSearch linesearch;
  SNESGetLineSearch(nl_solver_, &linesearch);
  SNESLineSearchSetType(linesearch, SNESLINESEARCHBT);

  KSP ksp;
  SNESGetKSP(nl_solver_, &ksp);
  KSPSetType(ksp, options_.l_method.c_str());

  KSPSetOptionsPrefix(ksp, solver_name_.c_str());

  this->ApplyToleranceOptions();

  this->SetOptions();

  this->SetSolverContext();
  this->SetOptions();

  this->SetSolverContext();
  this->SetConvergenceTest();
  this->SetMonitor();

  this->SetSystemSize();
  this->SetSystem();

  this->SetFunction();
  this->SetJacobian();

  this->SetPreconditioner();

  this->PostSetupCallback();
  system_set_ = true;
}

void
NonLinearSolver::PreSolveCallback()
{
}

void
NonLinearSolver::PostSolveCallback()
{
}

void
NonLinearSolver::Solve()
{
  converged_ = false;
  converged_reason_string_ = "Reason not obtained";
  this->PreSolveCallback();
  this->SetInitialGuess();

  SNESSolve(nl_solver_, nullptr, x_);
  this->PostSolveCallback();

  SNESConvergedReason conv_reason;
  SNESGetConvergedReason(nl_solver_, &conv_reason);

  if (conv_reason > 0)
    converged_ = true;

  const char* strreason;
  SNESGetConvergedReasonString(nl_solver_, &strreason);

  converged_reason_string_ = std::string(strreason);
}

std::string
NonLinearSolver::GetConvergedReasonString() const
{
  std::stringstream outstr;
  if (converged_)
    outstr << StringStreamColor(FG_GREEN) << std::string(10, ' ') << "Converged "
           << converged_reason_string_ << StringStreamColor(RESET);
  else
    outstr << StringStreamColor(FG_RED) << std::string(10, ' ') << "Convergence failure "
           << converged_reason_string_ << StringStreamColor(RESET);

  return outstr.str();
}

} // namespace opensn
