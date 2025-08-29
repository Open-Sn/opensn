// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/iterative_methods/wgs_context.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/linear_solver/linear_system_solver.h"
#include <memory>
#include <vector>
#include <functional>

namespace opensn
{

/**
 * Linear Solver specialization for Within GroupSet (WGS) solves with classic
 * Richardson.
 */
class ClassicRichardson : public LinearSystemSolver
{
public:
  /**
   * Constructor.
   * \param gs_context_ptr Context Pointer to abstract context.
   */
  explicit ClassicRichardson(const std::shared_ptr<WGSContext>& gs_context_ptr, bool verbose);

  ~ClassicRichardson() override = default;

  void Solve() override;

private:
  bool verbose_;
  std::vector<double> saved_q_moments_local_;
  std::vector<double> psi_new_, psi_old_;
};

} // namespace opensn
