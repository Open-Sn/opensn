// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <map>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

/// Compute balance
void ComputeBalance(DiscreteOrdinatesProblem& do_problem, double scaling_factor = 1.0);

/**
 * Computes the angular flux based leakage from boundary surfaces.
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
 * \param boundary_ids The boundary ids to compute leakages on.
 * \return A map of boundary ids to group-wise leakages.
 */
std::map<uint64_t, std::vector<double>> ComputeLeakage(DiscreteOrdinatesProblem& do_problem,
                                                       const std::vector<uint64_t>& boundary_ids);

} // namespace opensn
