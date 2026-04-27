// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/nonlinear_solver/petsc_nonlinear_solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/runtime.h"

namespace opensn
{

PETScNonLinearSolver::PETScNonLinearSolver(std::shared_ptr<NonLinearSolverContext> context_ptr,
                                           const InputParameters& params)
  : NonLinearSolver(params.GetParamValue<std::string>("name"), context_ptr),
    J_(nullptr),
    P_(nullptr),
    r_(nullptr),
    x_(nullptr),
    nl_solver_(nullptr),
    num_local_dofs_(0),
    num_global_dofs_(0),
    options_(params),
    system_set_(false),
    converged_(false)
{
}

PETScNonLinearSolver::~PETScNonLinearSolver()
{
  OpenSnPETScCall(SNESDestroy(&nl_solver_));
  OpenSnPETScCall(VecDestroy(&x_));
  OpenSnPETScCall(VecDestroy(&r_));
  OpenSnPETScCall(MatDestroy(&J_));
}

void
PETScNonLinearSolver::ApplyToleranceOptions()
{
  OpenSnPETScCall(SNESSetTolerances(nl_solver_,
                                    options_.nl_abs_tol,
                                    options_.nl_rel_tol,
                                    options_.nl_sol_tol,
                                    options_.nl_max_its,
                                    options_.nl_max_r_evaluations));
  OpenSnPETScCall(SNESSetMaxLinearSolveFailures(nl_solver_, options_.l_max_failed_iterations));
  KSP ksp = nullptr;
  OpenSnPETScCall(SNESGetKSP(nl_solver_, &ksp));
  OpenSnPETScCall(KSPSetTolerances(
    ksp, options_.l_rel_tol, options_.l_abs_tol, options_.l_div_tol, options_.l_max_its));
  if (options_.l_method == "gmres")
  {
    OpenSnPETScCall(KSPGMRESSetRestart(ksp, options_.l_gmres_restart_intvl));
    OpenSnPETScCall(KSPGMRESSetBreakdownTolerance(ksp, options_.l_gmres_breakdown_tol));
  }
  OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp, PETSC_FALSE));
}

void
PETScNonLinearSolver::PreSetupCallback()
{
}

void
PETScNonLinearSolver::SetOptions()
{
}

void
PETScNonLinearSolver::SetSolverContext()
{
  OpenSnPETScCall(SNESSetApplicationContext(nl_solver_, &(*context_ptr_)));
}

void
PETScNonLinearSolver::SetConvergenceTest()
{
}

void
PETScNonLinearSolver::SetMonitor()
{
}

void
PETScNonLinearSolver::SetPreconditioner()
{
}

void
PETScNonLinearSolver::PostSetupCallback()
{
}

void
PETScNonLinearSolver::Setup()
{
  if (IsSystemSet())
    return;
  this->PreSetupCallback();

  OpenSnPETScCall(SNESCreate(mpi_comm, &nl_solver_));

  OpenSnPETScCall(SNESSetOptionsPrefix(nl_solver_, solver_name_.c_str()));

  OpenSnPETScCall(SNESSetType(nl_solver_, options_.petsc_snes_type.c_str()));
  if (options_.nl_method == "LINEAR")
    OpenSnPETScCall(SNESSetType(nl_solver_, SNESKSPONLY));
  SNESLineSearch linesearch = nullptr;
  OpenSnPETScCall(SNESGetLineSearch(nl_solver_, &linesearch));
  OpenSnPETScCall(SNESLineSearchSetType(linesearch, SNESLINESEARCHBT));

  KSP ksp = nullptr;
  OpenSnPETScCall(SNESGetKSP(nl_solver_, &ksp));
  OpenSnPETScCall(KSPSetType(ksp, options_.l_method.c_str()));

  OpenSnPETScCall(KSPSetOptionsPrefix(ksp, solver_name_.c_str()));

  ApplyToleranceOptions();

  SetOptions();
  SetSolverContext();
  SetConvergenceTest();
  SetMonitor();

  SetSystemSize();
  SetSystem();

  SetFunction();
  SetJacobian();

  SetPreconditioner();

  PostSetupCallback();
  system_set_ = true;
}

void
PETScNonLinearSolver::PreSolveCallback()
{
}

void
PETScNonLinearSolver::PostSolveCallback()
{
}

void
PETScNonLinearSolver::Solve()
{
  converged_ = false;
  PreSolveCallback();
  SetInitialGuess();

  OpenSnPETScCall(SNESSolve(nl_solver_, nullptr, x_));
  PostSolveCallback();

  SNESConvergedReason conv_reason = SNES_CONVERGED_ITERATING;
  OpenSnPETScCall(SNESGetConvergedReason(nl_solver_, &conv_reason));

  if (SNESReasonToPETScSolverStatus(conv_reason) == PETScSolverStatus::CONVERGED)
    converged_ = true;
}

} // namespace opensn
