// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <map>
#include <optional>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

struct BalanceTable
{
  double absorption_rate = 0.0;
  double production_rate = 0.0;
  double inflow_rate = 0.0;
  double outflow_rate = 0.0;
  double balance = 0.0;
  std::optional<double> csda_charge_deposition_rate;
  std::optional<double> csda_particle_deposition_rate;
  /// Signed particle residual divided by production plus boundary inflow.
  /// Zero gain gives zero for a zero residual, otherwise signed infinity.
  std::optional<double> csda_particle_balance;
  /// Magnitude of csda_particle_balance.
  std::optional<double> csda_particle_relative_balance;
  /// Midpoint-weighted collision loss plus CSDA loss, including terminal cutoff.
  /// This is distinct from the imported CEPXS deposition response field integral.
  std::optional<double> csda_energy_deposition_rate;
  std::optional<double> csda_energy_collision_loss_rate;
  std::optional<double> csda_energy_continuous_loss_rate;
  std::optional<double> csda_energy_production_rate;
  std::optional<double> csda_energy_inflow_rate;
  std::optional<double> csda_energy_outflow_rate;
  /// Signed energy residual divided by energy production plus boundary inflow.
  /// Zero gain gives zero for a zero residual, otherwise signed infinity.
  std::optional<double> csda_energy_balance;
  /// Magnitude of csda_energy_balance.
  std::optional<double> csda_energy_relative_balance;
  std::optional<double> initial_inventory;
  std::optional<double> final_inventory;
  std::optional<double> predicted_inventory_change;
  std::optional<double> actual_inventory_change;
  std::optional<double> inventory_residual;
};

BalanceTable ComputeBalanceTable(DiscreteOrdinatesProblem& do_problem, double scaling_factor = 1.0);

/// Compute balance
void ComputeBalance(DiscreteOrdinatesProblem& do_problem, double scaling_factor = 1.0);

/**
 * Computes the angular flux based leakage from boundary surfaces.
 * \param do_problem The discrete ordinates problem supplying the current angular flux state.
 * \param groupset_id The groupset for which to compute the leakage.
 * \param boundary_id The boundary id for which to perform the integration.
 *
 * \return The leakage as a value.
 */
std::vector<double> ComputeLeakage(DiscreteOrdinatesProblem& do_problem,
                                   unsigned int groupset_id,
                                   uint64_t boundary_id);

/**
 * Computes the group-wise angular flux-based leakage from the specified boundaries.
 *
 * \param do_problem The discrete ordinates problem supplying the current angular flux state.
 * \param boundary_ids The boundary ids to compute leakages on.
 * \return A map of boundary ids to group-wise leakages.
 */
std::map<uint64_t, std::vector<double>> ComputeLeakage(DiscreteOrdinatesProblem& do_problem,
                                                       const std::vector<uint64_t>& boundary_ids);

} // namespace opensn
