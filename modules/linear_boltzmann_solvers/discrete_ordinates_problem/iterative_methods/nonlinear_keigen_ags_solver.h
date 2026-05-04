// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/nonlinear_keigen_ags_context.h"
#include "framework/math/nonlinear_solver/petsc_nonlinear_solver.h"

namespace opensn
{

class NLKEigenvalueAGSSolver : public PETScNonLinearSolver
{
public:
  explicit NLKEigenvalueAGSSolver(std::shared_ptr<NLKEigenAGSContext> nlk_ags_context_ptr)
    : PETScNonLinearSolver(nlk_ags_context_ptr)
  {
  }

  ~NLKEigenvalueAGSSolver() override = default;

  void Solve() override;

protected:
  void PreSetupCallback() override;
  void SetMonitor() override;
  void SetSystemSize() override;
  void SetSystem() override;
  void SetFunction() override;
  void SetJacobian() override;
  void SetInitialGuess() override;
  void PostSolveCallback() override;
};

} // namespace opensn
