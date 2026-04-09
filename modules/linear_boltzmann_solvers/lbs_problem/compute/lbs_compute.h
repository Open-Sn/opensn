// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <optional>
#include <vector>

namespace opensn
{

class LBSProblem;

struct BalanceTable
{
  double absorption_rate = 0.0;
  double production_rate = 0.0;
  double inflow_rate = 0.0;
  double outflow_rate = 0.0;
  double balance = 0.0;
  std::optional<double> initial_inventory;
  std::optional<double> final_inventory;
  std::optional<double> predicted_inventory_change;
  std::optional<double> actual_inventory_change;
  std::optional<double> inventory_residual;
};

/**
 * Compute the total fission production in the problem.
 */
double ComputeFissionProduction(const LBSProblem& lbs_problem, const std::vector<double>& phi);

/**
 * Computes the total fission rate in the problem.
 */
double ComputeFissionRate(const LBSProblem& lbs_problem, const std::vector<double>& phi);

/// Compute the steady state delayed neutron precursor concentrations.
void ComputePrecursors(LBSProblem& lbs_problem);

} // namespace opensn
