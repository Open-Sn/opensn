// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/linear_solver/petsc_linear_solver.h"
#include "framework/runtime.h"

namespace opensn
{

PETScLinearSolver::PETScLinearSolver(IterativeMethod iterative_method,
                                     std::shared_ptr<LinearSolverContext> context_ptr)
  : LinearSolver(iterative_method, context_ptr),
    A_(nullptr),
    b_(nullptr),
    x_(nullptr),
    ksp_(nullptr),
    num_local_dofs_(0),
    num_global_dofs_(0),
    system_set_(false),
    suppress_kspsolve_(false)
{
}

PETScLinearSolver::~PETScLinearSolver()
{
  VecDestroy(&x_);
  VecDestroy(&b_);
  KSPDestroy(&ksp_);
}

void
PETScLinearSolver::ApplyToleranceOptions()
{
  KSPSetTolerances(ksp_,
                   tolerance_options.residual_relative,
                   tolerance_options.residual_absolute,
                   tolerance_options.residual_divergence,
                   tolerance_options.maximum_iterations);
}

void
PETScLinearSolver::PreSetupCallback()
{
}

void
PETScLinearSolver::SetOptions()
{
}

void
PETScLinearSolver::SetSolverContext()
{
  KSPSetApplicationContext(ksp_, &(*context_ptr_));
}

void
PETScLinearSolver::SetConvergenceTest()
{
  KSPSetConvergenceTest(ksp_, &KSPConvergedDefault, nullptr, nullptr);
}

void
PETScLinearSolver::SetMonitor()
{
}

void
PETScLinearSolver::SetPreconditioner()
{
}

void
PETScLinearSolver::PostSetupCallback()
{
}

void
PETScLinearSolver::Setup()
{
  if (IsSystemSet())
    return;

  PreSetupCallback();

  KSPCreate(opensn::mpi_comm, &ksp_);

  const auto petsc_iterative_method = PETScIterativeMethodName(iterative_method_);
  KSPSetType(ksp_, petsc_iterative_method.c_str());

  ApplyToleranceOptions();

  if (iterative_method_ == IterativeMethod::PETSC_GMRES)
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
PETScLinearSolver::PreSolveCallback()
{
}

void
PETScLinearSolver::PostSolveCallback()
{
}

void
PETScLinearSolver::Solve()
{
  PreSolveCallback();
  SetInitialGuess();
  SetRHS();
  if (not suppress_kspsolve_)
    KSPSolve(ksp_, b_, x_);
  PostSolveCallback();
}

std::string
PETScLinearSolver::PETScIterativeMethodName(opensn::LinearSolver::IterativeMethod iterative_method)
{
  switch (iterative_method)
  {
    case IterativeMethod::NONE:
      return "preonly";
    case IterativeMethod::PETSC_RICHARDSON:
      return "richardson";
    case IterativeMethod::PETSC_GMRES:
      return "gmres";
    case IterativeMethod::PETSC_BICGSTAB:
      return "bcgs";
    default:
      throw std::runtime_error("Cannot get a PETSc option name for a non-PETSc iterative method.");
  }
}

int
PETScLinearSolver::LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action)
{
  LinearSolverContext* context;
  MatShellGetContext(matrix, &context);

  context->MatrixAction(matrix, vector, action);

  return 0;
}

} // namespace opensn
