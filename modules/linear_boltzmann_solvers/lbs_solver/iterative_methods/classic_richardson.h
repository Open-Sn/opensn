// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_solver/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "framework/math/linear_solver/linear_solver.h"
#include <memory>
#include <vector>
#include <functional>

namespace opensn
{

/**
 * Linear Solver specialization for Within GroupSet (WGS) solves with classic
 * Richardson.
 */
class ClassicRichardson : public LinearSolver
{
public:
  /**
   * Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.
   */
  explicit ClassicRichardson(const std::shared_ptr<WGSContext>& gs_context_ptr);

  ~ClassicRichardson() override;

  void Setup() override {}

  void Solve() override;

protected:
  void PreSetupCallback() override {}

  void SetConvergenceTest() override {}

  void SetSystemSize() override {}

  void SetSystem() override {}

  void SetPreconditioner() override {}

  void PostSetupCallback() override {}

  void PreSolveCallback() override {}

  void SetRHS() override {}

  void SetInitialGuess() override {}

  void PostSolveCallback() override {}

private:
  std::vector<double> saved_q_moments_local_;
  std::vector<double> psi_new_, psi_old_;
};

} // namespace opensn
