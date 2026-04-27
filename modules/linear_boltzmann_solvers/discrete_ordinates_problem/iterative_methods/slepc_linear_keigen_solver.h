// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/math/linear_solver/slepc_linear_eigen_solver.h"
#include "framework/math/linear_solver/linear_solver_context.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include <memory>

namespace opensn
{

/// Context for the linear k-eigenvalue solve
struct SLEPcLinearKEigenContext : public LinearEigenContext
{
  std::shared_ptr<DiscreteOrdinatesProblem> do_problem;
  double last_eps_residual = 1.0;
  double last_operator_residual = 1.0;

  explicit SLEPcLinearKEigenContext(std::shared_ptr<DiscreteOrdinatesProblem> do_problem)
    : do_problem(std::move(do_problem))
  {
    eigenvalue = 1.0;
  }

  ~SLEPcLinearKEigenContext() override = default;
};

/// A SLEPc shell-based linear k-eigenvalue solver
class SLEPcLinearKEigenSolver : public SLEPcLinearEigenSolver
{
public:
  explicit SLEPcLinearKEigenSolver(std::shared_ptr<SLEPcLinearKEigenContext> context)
    : SLEPcLinearEigenSolver(IterativeMethod::KRYLOV_SCHUR, std::move(context))
  {
  }

  ~SLEPcLinearKEigenSolver() override = default;

  void Solve() override;

protected:
  void SetMonitor() override;
  void SetSystemSize() override;
  void SetSystem() override;
  void SetInitialGuess() override;
  void PostSolveCallback() override;
  void PreSolveCallback() override;
};

} // namespace opensn
