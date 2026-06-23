// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

namespace opensn
{

class LBSProblem;
class MultiGroupXS;

/// Returns true when delayed neutron production contributes to transport sources.
bool UseDelayedNeutronProduction(const LBSProblem& lbs_problem);

/// Computes the delayed fission production from a source group to a destination group.
double ComputeDelayedFissionProduction(const MultiGroupXS& xs,
                                       unsigned int to_group,
                                       unsigned int from_group);

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
