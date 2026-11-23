// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/iterative_methods/wgs_context.h"
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
   * \param verbose If true, print per-iteration convergence information.
   */
  explicit ClassicRichardson(const std::shared_ptr<WGSContext>& gs_context_ptr, bool verbose);

  ~ClassicRichardson() override = default;

  void Solve() override;

  /// Frees the vectors that Solve() allocates.
  void ReleaseWorkspace() override;

private:
  bool verbose_;
  std::vector<double> psi_new_, psi_old_;
  /// The groupset's block of phi_new - phi_old, which drives the DSA correction.
  std::vector<double> phi_change_;
  /// The groupset's block of the scalar flux before acceleration, then the DSA correction.
  std::vector<double> phi_unaccelerated_;

  void SyncLaggedStateToLatestIterate(WGSContext& gs_context);
};

} // namespace opensn
