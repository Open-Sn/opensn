#include "framework/math/linear_solver/linear_solver.h"

namespace opensn
{

LinearSolver::LinearSolver(const std::string& iterative_method, LinSolveContextPtr context_ptr)
  : solver_name_(iterative_method), iterative_method_(iterative_method), context_ptr_(context_ptr)
{
}

LinearSolver::LinearSolver(const std::string& solver_name,
                           const std::string& iterative_method,
                           LinSolveContextPtr context_ptr)
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
                   tolerance_options_.residual_relative,
                   tolerance_options_.residual_absolute,
                   tolerance_options_.residual_divergence,
                   tolerance_options_.maximum_iterations);
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
  if (IsSystemSet()) return;
  PreSetupCallback();

  KSPCreate(PETSC_COMM_WORLD, &ksp_);
  KSPSetType(ksp_, iterative_method_.c_str());

  ApplyToleranceOptions();

  if (iterative_method_ == "gmres")
  {
    KSPGMRESSetRestart(ksp_, tolerance_options_.gmres_restart_interval);
    KSPGMRESSetBreakdownTolerance(ksp_, tolerance_options_.gmres_breakdown_tolerance);
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
  if (not suppress_kspsolve_) KSPSolve(ksp_, b_, x_);
  PostSolveCallback();
}

} // namespace opensn
