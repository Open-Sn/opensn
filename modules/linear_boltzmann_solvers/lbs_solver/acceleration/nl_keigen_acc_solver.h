// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/acceleration/nl_keigen_acc_context.h"
#include "framework/math/nonlinear_solver/nonlinear_solver.h"
#include <petscsnes.h>

namespace opensn
{

class NLKEigenDiffSolver : public NonLinearSolver
{
public:
  explicit NLKEigenDiffSolver(std::shared_ptr<NLKEigenDiffContext> nlk_diff_context_ptr)
    : NonLinearSolver(nlk_diff_context_ptr)
  {
  }

  ~NLKEigenDiffSolver() override = default;

protected:
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