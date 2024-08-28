// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/nl_keigen_ags_context.h"
#include "framework/math/nonlinear_solver/nonlinear_solver.h"

namespace opensn
{

class NLKEigenvalueAGSSolver : public NonLinearSolver
{
public:
  explicit NLKEigenvalueAGSSolver(std::shared_ptr<NLKEigenAGSContext> nlk_ags_context_ptr)
    : NonLinearSolver(nlk_ags_context_ptr)
  {
  }

  ~NLKEigenvalueAGSSolver() override = default;

protected:
  void PreSetupCallback() override;
  void SetMonitor() override;

  void SetSystemSize() override;
  void SetSystem() override;
  void SetFunction() override;
  void SetJacobian() override;

protected:
  void SetInitialGuess() override;
  void PostSolveCallback() override;
};

} // namespace opensn
