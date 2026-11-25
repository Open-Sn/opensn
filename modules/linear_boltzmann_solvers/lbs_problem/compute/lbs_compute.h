// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

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
 * Compute the total fission production in the problem.
 */
double ComputeFissionProduction(const LBSProblem& lbs_problem, const std::vector<double>& phi);

/**
 * Fission production of each fissionable local cell divided by the global fission production, in
 * local cell order. The sources of all ranks sum to one. Collective. Throws if the global
 * production is not positive and finite.
 */
std::vector<double> ComputeNormalizedCellFissionSources(const LBSProblem& lbs_problem,
                                                        const std::vector<double>& phi);

/**
 * Root-mean-square over all fissionable cells of the relative change from previous to current,
 * two results of ComputeNormalizedCellFissionSources for the same problem. A cell whose previous
 * source is not positive adds nothing to the sum but counts as a cell, as in OpenMOC. Returns zero
 * when no rank has a fissionable cell. Collective.
 */
double ComputeFissionSourceRMSChange(const std::vector<double>& previous,
                                     const std::vector<double>& current);

/**
 * Computes the total fission rate in the problem.
 */
double ComputeFissionRate(const LBSProblem& lbs_problem, const std::vector<double>& phi);

/// Compute the steady state delayed neutron precursor concentrations (zero in adjoint mode).
void ComputePrecursors(LBSProblem& lbs_problem);

} // namespace opensn
