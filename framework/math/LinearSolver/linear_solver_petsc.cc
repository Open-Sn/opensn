#include "opensn/framework/math/LinearSolver/linear_solver.h"

#include <petscksp.h>

namespace chi_math
{

template <>
LinearSolver<Mat, Vec, KSP>::~LinearSolver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  KSPDestroy(&solver_);
}

template <>
void
LinearSolver<Mat, Vec, KSP>::ApplyToleranceOptions()
{
  KSPSetTolerances(solver_,
                   tolerance_options_.residual_relative,
                   tolerance_options_.residual_absolute,
                   tolerance_options_.residual_divergence,
                   tolerance_options_.maximum_iterations);
}

template <>
void
LinearSolver<Mat, Vec, KSP>::PreSetupCallback()
{
}

template <>
void
LinearSolver<Mat, Vec, KSP>::SetOptions()
{
}

template <>
void
LinearSolver<Mat, Vec, KSP>::SetSolverContext()
{
  KSPSetApplicationContext(solver_, &(*context_ptr_));
}

template <>
void
LinearSolver<Mat, Vec, KSP>::SetConvergenceTest()
{
  KSPSetConvergenceTest(solver_, &KSPConvergedDefault, nullptr, nullptr);
}

template <>
void
LinearSolver<Mat, Vec, KSP>::SetMonitor()
{
}

template <>
void
LinearSolver<Mat, Vec, KSP>::SetPreconditioner()
{
}

template <>
void
LinearSolver<Mat, Vec, KSP>::PostSetupCallback()
{
}

template <>
void
LinearSolver<Mat, Vec, KSP>::Setup()
{
  if (IsSystemSet()) return;
  this->PreSetupCallback();

  KSPCreate(PETSC_COMM_WORLD, &solver_);
  KSPSetType(solver_, iterative_method_.c_str());

  this->ApplyToleranceOptions();

  if (iterative_method_ == "gmres")
  {
    KSPGMRESSetRestart(solver_, tolerance_options_.gmres_restart_interval);
    KSPGMRESSetBreakdownTolerance(solver_, tolerance_options_.gmres_breakdown_tolerance);
  }

  KSPSetInitialGuessNonzero(solver_, PETSC_FALSE);

  this->SetOptions();

  this->SetSolverContext();
  this->SetConvergenceTest();
  this->SetMonitor();

  this->SetSystemSize();
  this->SetSystem();

  this->SetPreconditioner();

  this->PostSetupCallback();
  system_set_ = true;
}

template <>
void
LinearSolver<Mat, Vec, KSP>::PreSolveCallback()
{
}

template <>
void
LinearSolver<Mat, Vec, KSP>::PostSolveCallback()
{
}

template <>
void
LinearSolver<Mat, Vec, KSP>::Solve()
{
  this->PreSolveCallback();
  this->SetInitialGuess();
  this->SetRHS();

  if (not suppress_kspsolve_) KSPSolve(solver_, b_, x_);
  this->PostSolveCallback();
}

} // namespace chi_math
