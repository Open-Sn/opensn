// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/linear_solver/petsc_linear_system_solver.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/runtime.h"

namespace opensn
{

PETScLinearSolver::PETScLinearSolver(IterativeMethod method,
                                     std::shared_ptr<LinearSystemContext> context_ptr)
  : LinearSystemSolver(method, context_ptr),
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
  OpenSnPETScCall(KSPDestroy(&ksp_));
  OpenSnPETScCall(MatDestroy(&A_));
  OpenSnPETScCall(VecDestroy(&x_));
  OpenSnPETScCall(VecDestroy(&b_));
}

void
PETScLinearSolver::ApplyToleranceOptions()
{
  OpenSnPETScCall(KSPSetTolerances(ksp_,
                                   tolerance_options.residual_relative,
                                   tolerance_options.residual_absolute,
                                   tolerance_options.residual_divergence,
                                   tolerance_options.maximum_iterations));
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
  OpenSnPETScCall(KSPSetApplicationContext(ksp_, &(*context_ptr_)));
}

void
PETScLinearSolver::SetConvergenceTest()
{
  OpenSnPETScCall(KSPSetConvergenceTest(ksp_, &KSPConvergedDefault, nullptr, nullptr));
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
  OpenSnPETScCall(KSPCreate(opensn::mpi_comm, &ksp_));
  const auto petsc_iterative_method = PETScIterativeMethodName();
  OpenSnPETScCall(KSPSetType(ksp_, petsc_iterative_method.c_str()));
  ApplyToleranceOptions();
  if (method_ == IterativeMethod::PETSC_GMRES)
  {
    OpenSnPETScCall(KSPGMRESSetRestart(ksp_, tolerance_options.gmres_restart_interval));
    OpenSnPETScCall(
      KSPGMRESSetBreakdownTolerance(ksp_, tolerance_options.gmres_breakdown_tolerance));
  }
  OpenSnPETScCall(KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE));
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
    OpenSnPETScCall(KSPSolve(ksp_, b_, x_));
  PostSolveCallback();
}

std::string
PETScLinearSolver::PETScIterativeMethodName()
{
  switch (method_)
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
      throw std::runtime_error("Unsupported PETSc iterative method.");
  }
}

PetscErrorCode
PETScLinearSolver::LinearSolverMatrixAction(Mat matrix, Vec vector, Vec action)
{
  LinearSystemContext* context = nullptr;
  const PetscErrorCode ierr = MatShellGetContext(matrix, static_cast<void*>(&context));
  if (ierr != PETSC_SUCCESS)
    return ierr;

  context->MatrixAction(matrix, vector, action);

  return PETSC_SUCCESS;
}

} // namespace opensn
