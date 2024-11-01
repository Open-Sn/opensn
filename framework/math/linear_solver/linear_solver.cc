// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/linear_solver/linear_solver.h"
#include "framework/runtime.h"

namespace opensn
{

LinearSolver::LinearSolver(const std::string& petsc_iterative_method,
                           std::shared_ptr<LinearSolverContext> context_ptr)
  : context_ptr_(context_ptr), petsc_iterative_method_(petsc_iterative_method)
{
}

LinearSolver::~LinearSolver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  KSPDestroy(&ksp_);
}

void
LinearSolver::ApplyToleranceOptions()
{
  KSPSetTolerances(ksp_,
                   tolerance_options.residual_relative,
                   tolerance_options.residual_absolute,
                   tolerance_options.residual_divergence,
                   tolerance_options.maximum_iterations);
}

void
LinearSolver::PreSetupCallback()
{
}

void
LinearSolver::SetOptions()
{
}

void
LinearSolver::SetSolverContext()
{
  KSPSetApplicationContext(ksp_, &(*context_ptr_));
}

void
LinearSolver::SetConvergenceTest()
{
  KSPSetConvergenceTest(ksp_, &KSPConvergedDefault, nullptr, nullptr);
}

void
LinearSolver::SetMonitor()
{
}

void
LinearSolver::SetPreconditioner()
{
}

void
LinearSolver::PostSetupCallback()
{
}

void
LinearSolver::Setup()
{
  if (IsSystemSet())
    return;
  PreSetupCallback();

  KSPCreate(opensn::mpi_comm, &ksp_);

  KSPSetType(ksp_, petsc_iterative_method_.c_str());

  ApplyToleranceOptions();

  if (petsc_iterative_method_ == "gmres")
  {
    KSPGMRESSetRestart(ksp_, tolerance_options.gmres_restart_interval);
    KSPGMRESSetBreakdownTolerance(ksp_, tolerance_options.gmres_breakdown_tolerance);
  }

  KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);

  SetOptions();

  SetSolverContext();
  SetConvergenceTest();
  SetMonitor();

  SetSystemSize();
  SetSystem();

  SetPreconditioner();

  PostSetupCallback();
  system_set_ = true;
}

void
LinearSolver::PreSolveCallback()
{
}

void
LinearSolver::PostSolveCallback()
{
}

void
LinearSolver::Solve()
{
  PreSolveCallback();
  SetInitialGuess();
  SetRHS();
  if (not suppress_kspsolve_)
    KSPSolve(ksp_, b_, x_);
  PostSolveCallback();
}

} // namespace opensn
