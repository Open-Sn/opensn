// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <vector>

namespace opensn
{
class CMFDCoarseMesh;
class DiscreteOrdinatesProblem;
class LBSGroupset;

/**
 * Restrict the scalar flux moments from the transport grid to coarse-cell volume averages.
 */
std::vector<double> CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                                           const LBSGroupset& groupset,
                                           const CMFDCoarseMesh& coarse_mesh,
                                           const std::vector<double>& phi);
std::vector<double> CMFDRestrictScalarFlux(const DiscreteOrdinatesProblem& do_problem,
                                           unsigned int first_group,
                                           unsigned int num_groups,
                                           const CMFDCoarseMesh& coarse_mesh,
                                           const std::vector<double>& phi);

/**
 * Add a coarse-cell scalar-flux correction back to the transport scalar flux moments.
 */
void CMFDProlongateScalarFluxCorrection(const DiscreteOrdinatesProblem& do_problem,
                                        const LBSGroupset& groupset,
                                        const CMFDCoarseMesh& coarse_mesh,
                                        const std::vector<double>& coarse_delta_phi,
                                        std::vector<double>& phi);
void CMFDProlongateScalarFluxCorrection(const DiscreteOrdinatesProblem& do_problem,
                                        unsigned int first_group,
                                        unsigned int num_groups,
                                        const CMFDCoarseMesh& coarse_mesh,
                                        const std::vector<double>& coarse_delta_phi,
                                        std::vector<double>& phi);

} // namespace opensn
