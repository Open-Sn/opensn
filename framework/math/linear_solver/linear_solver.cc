// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/linear_solver/linear_solver.h"
#include "framework/runtime.h"

namespace opensn
{

LinearSolver::LinearSolver(const std::string& iterative_method,
                           std::shared_ptr<LinearSolverContext> context_ptr)
  : solver_name_(iterative_method), iterative_method_(iterative_method), context_ptr_(context_ptr)
{
}

LinearSolver::LinearSolver(const std::string& solver_name,
                           const std::string& iterative_method,
                           std::shared_ptr<LinearSolverContext> context_ptr)
  : solver_name_(solver_name), iterative_method_(iterative_method), context_ptr_(context_ptr)
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

  // In OpenSn the PETSc version of Richardson iteration is referred to as krylov_richardson to
  // distinguish it from classic Richardson. At this point, we need to convert from
  // krylov_richardson to the correct PETSc algorithm name.
  if (iterative_method_ == "krylov_richardson")
    KSPSetType(ksp_, "richardson");
  else
    KSPSetType(ksp_, iterative_method_.c_str());

  ApplyToleranceOptions();

  if (iterative_method_ == "gmres")
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
