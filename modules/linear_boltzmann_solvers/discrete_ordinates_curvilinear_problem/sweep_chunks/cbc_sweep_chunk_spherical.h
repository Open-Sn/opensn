// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_kernels.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include <vector>

namespace opensn
{

class LBSGroupset;
class DiscreteOrdinatesProblem;
class CurvilinearProductQuadrature;

/**
 * Host CBC sweep chunk for 1D spherical geometry.
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
 * set precedes the outward set. SINGLE aggregation chains individual directions in the same order.
 * The CBC scheduler invokes this chunk for one ready cell at a time. Angular-edge storage is local
 * to each owned cell and groupset; only spatial upwind fluxes cross rank boundaries. Outer-surface
 * reflection is delayed by the boundary machinery so that the angular sequence remains acyclic.
 * The chunk is used serially by the host scheduler and lives no longer than its problem/groupset.
 *
 * A boundary face at r = 0 (the center of a solid sphere) has zero area and receives no incident
 * flux; any other boundary face, including the inner surface of a hollow sphere, uses the
 * boundary condition.
 */
class CBCSweepChunkSpherical : public SweepChunk
{
public:
  CBCSweepChunkSpherical(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  /// Select the CBC flux storage and communicator before sweeping cells of an angle set.
  void SetAngleSet(AngleSet& angle_set) override;

  /// Sweep the cell selected by SetCell through the angle set in increasing mu.
  void Sweep(AngleSet& angle_set) override;

private:
  void Sweep_Generic(AngleSet& angle_set);

  /// Two-node radial cells use the same blocked SIMD solve as the RZ chunks.
  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);

  /// Validated spherical quadrature, owned by the groupset.
  const CurvilinearProductQuadrature& curvilinear_quadrature_;
  /// Number of groups solved in one block.
  const unsigned int group_block_size_;
  /// Mass matrices with radial weight r for the 1/r angular-redistribution term.
  const std::vector<UnitCellMatrices>& secondary_unit_cell_matrices_;
  /// Unknown manager for the angular-edge flux (one vector of groupset groups per node).
  UnknownManager unknown_manager_;
  /// Angular-edge flux psi_{n+1/2} of the most recently swept direction at each node.
  std::vector<double> psi_edge_;
  /// CBC flux storage and communicator for the active angle set (owned by that angle set).
  CBC_FLUDS* fluds_ = nullptr;
  CBC_AsynchronousCommunicator* async_comm_ = nullptr;
  /// Staged nonlocal outgoing face fluxes, reused between cell tasks.
  CBCSweepWorkspace workspace_;
  /// Dense cell solve buffers, reused between cell tasks.
  DenseMatrix<double> streaming_matrix_;
  DenseMatrix<double> system_matrix_;
  std::vector<Vector<double>> rhs_;
  std::vector<double> source_;
};

} // namespace opensn
