// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/linear_solver/slepc_linear_eigen_solver.h"
#include "framework/runtime.h"

namespace opensn
{

SLEPcLinearEigenSolver::SLEPcLinearEigenSolver(IterativeMethod method,
                                               std::shared_ptr<LinearEigenContext> context_ptr)
  : LinearEigenSolver(method, context_ptr),
    A_(nullptr),
    x_(nullptr),
    eps_(nullptr),
    eps_type_("krylovschur"),
    num_local_dofs_(0),
    num_global_dofs_(0),
    system_set_(false)
{
}

SLEPcLinearEigenSolver::~SLEPcLinearEigenSolver()
{
  VecDestroy(&x_);
  EPSDestroy(&eps_);
}

void
SLEPcLinearEigenSolver::PreSetupCallback()
{
}

void
SLEPcLinearEigenSolver::SetOptions()
{
}

void
SLEPcLinearEigenSolver::SetSolverContext()
{
}

void
SLEPcLinearEigenSolver::SetConvergenceTest()
{
}

void
SLEPcLinearEigenSolver::SetMonitor()
{
}

void
SLEPcLinearEigenSolver::SetPreconditioner()
{
}

void
SLEPcLinearEigenSolver::PostSetupCallback()
{
}

void
SLEPcLinearEigenSolver::Setup()
{
  if (IsSystemSet())
    return;

  PreSetupCallback();
  EPSCreate(opensn::mpi_comm, &eps_);
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
SLEPcLinearEigenSolver::PreSolveCallback()
{
}

void
SLEPcLinearEigenSolver::PostSolveCallback()
{
}

void
SLEPcLinearEigenSolver::Solve()
{
  PreSolveCallback();
  SetInitialGuess();
  PostSolveCallback();
}

} // namespace opensn
