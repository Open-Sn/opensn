// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <utility>
#include <vector>

namespace opensn
{
/// Global relative L2 balance residual, scaled before squaring to avoid underflow/overflow.
/// An absent source or nonfinite balance is not converged and returns infinity.
double CMFDRelativeBalanceResidual(const std::vector<std::pair<double, double>>& residual_rhs);

class CMFDCoarseMesh;
class DiscreteOrdinatesProblem;

/**
 * Restrict the scalar flux moments from the transport grid to coarse-cell volume averages.
 */
std::vector<double> CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                                           unsigned int first_group,
                                           unsigned int num_groups,
                                           unsigned int group_aggregation_size,
                                           const CMFDCoarseMesh& coarse_mesh,
                                           const std::vector<double>& phi);
std::vector<double> CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                                           unsigned int first_group,
                                           unsigned int num_groups,
                                           const CMFDCoarseMesh& coarse_mesh,
                                           const std::vector<double>& phi);

/**
 * Scale fine scalar-flux moments by coarse-cell, coarse-group ratios.
 */
void CMFDProlongateScalarFluxRatio(const DiscreteOrdinatesProblem& do_problem,
                                   unsigned int first_group,
                                   unsigned int num_groups,
                                   unsigned int group_aggregation_size,
                                   const CMFDCoarseMesh& coarse_mesh,
                                   const std::vector<double>& coarse_phi_old,
                                   const std::vector<double>& coarse_phi_new,
                                   std::vector<double>& phi);

} // namespace opensn
