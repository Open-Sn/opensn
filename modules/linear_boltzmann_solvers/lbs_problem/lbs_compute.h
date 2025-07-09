// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

namespace opensn
{

class LBSProblem;

/**
 * Compute the total fission production in the problem.
 */
double ComputeFissionProduction(LBSProblem& lbs_problem, const std::vector<double>& phi);

/// Compute the steady state delayed neutron precursor concentrations.
void ComputePrecursors(LBSProblem& lbs_problem);

} // namespace opensn
