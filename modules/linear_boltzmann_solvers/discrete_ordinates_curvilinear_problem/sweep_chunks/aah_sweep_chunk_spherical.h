// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include <vector>

namespace opensn
{

class LBSGroupset;
class DiscreteOrdinatesProblem;

/**
 * AAH sweep chunk for 1D spherical geometry.
 *
 * Solves, for each direction mu_n of a GLProductQuadrature1DSpherical, the PWLD discretization
 * (radial weight r^2) of
 *   mu d psi/dr + (1/r) [alpha_{n+1/2}/(w_n tau_n) + 2 mu_n] (psi_n - psi_{n-1/2})
 *     + sigma_t psi_n = q_n,
 * where psi_{n-1/2} is the angular-edge flux from the previous direction in order of increasing
 * mu and psi_{n+1/2} = (psi_n - (1 - tau_n) psi_{n-1/2}) / tau_n (weighted diamond difference in
 * angle). The recursion starts from the zero-weight direction mu = -1, which has no angular term.
 * The 1/r term uses mass matrices with weight r. Directions must reach each cell in order of
 * increasing mu: the angle aggregation groups them into an inward and an outward set and the inward
 * set precedes the outward set.
 *
 * A boundary face at r = 0 (the center of a solid sphere) has zero area and receives no incident
 * flux; any other boundary face, including the inner surface of a hollow sphere, uses the
 * boundary condition.
 */
class AAHSweepChunkSpherical : public SweepChunk
{
public:
  AAHSweepChunkSpherical(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  void Sweep(AngleSet& angle_set) override;

private:
  /// Mass matrices with radial weight r for the 1/r angular-redistribution term.
  const std::vector<UnitCellMatrices>& secondary_unit_cell_matrices_;
  /// Unknown manager for the angular-edge flux (one vector of groupset groups per node).
  UnknownManager unknown_manager_;
  /// Angular-edge flux psi_{n+1/2} of the most recently swept direction at each node.
  std::vector<double> psi_edge_;
};

} // namespace opensn
