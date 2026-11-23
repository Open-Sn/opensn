// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/lbs_problem/groupset_blocked_vector.h"
#include <span>
#include <vector>

namespace opensn
{

class LBSProblem;
class MultiGroupXS;

/**
 * Returns true when delayed neutron production contributes to fission sources: the cross
 * sections have delayed-neutron data and options.use_precursors is set. Every solve (steady
 * state, k-eigenvalue and transient) uses the same model.
 */
bool UseDelayedNeutronProduction(const LBSProblem& lbs_problem);

/// Computes the delayed fission production from a source group to a destination group.
double ComputeDelayedFissionProduction(const MultiGroupXS& xs,
                                       unsigned int to_group,
                                       unsigned int from_group);

/**
 * Compute the total fission production in the problem from the flux moments `phi`, given by
 * groupset block. Collective.
 */
double ComputeFissionProduction(const LBSProblem& lbs_problem,
                                const GroupsetBlockedSpan<const double>& phi);

/// ComputeFissionProduction of a whole flux-moment vector `phi` in the problem's layout.
double ComputeFissionProduction(const LBSProblem& lbs_problem, std::span<const double> phi);

/**
 * Computes the total fission rate in the problem from the flux moments `phi`, given by groupset
 * block. Collective.
 */
double ComputeFissionRate(const LBSProblem& lbs_problem,
                          const GroupsetBlockedSpan<const double>& phi);

/// ComputeFissionRate of a whole flux-moment vector `phi` in the problem's layout.
double ComputeFissionRate(const LBSProblem& lbs_problem, std::span<const double> phi);

/// Compute the steady state delayed neutron precursor concentrations (zero in adjoint mode).
void ComputePrecursors(LBSProblem& lbs_problem);

} // namespace opensn
