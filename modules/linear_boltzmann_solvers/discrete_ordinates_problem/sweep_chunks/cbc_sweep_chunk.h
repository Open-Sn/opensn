// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk_shared.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{

class CellMapping;
class DiscreteOrdinatesProblem;

/**
 * Implements the core sweep operation for a single cell within the
 * cell-by-cell (CBC) sweep algorithm.
 *
 * This class is responsible for performing the discrete ordinates transport
 * calculation on a given cell for all angles and groups managed by its
 * current AngleSet
 * It interacts with a CBC_FLUDS object to obtain upwind angular flux data
 * (from local neighbors, MPI remote buffers, or boundaries) and to store
 * outgoing angular flux data (to local neighbors or MPI send buffers)
 */
class CBCSweepChunk : public SweepChunk
{
public:
  CBCSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset);

  /// Set the current AngleSet
  void SetAngleSet(AngleSet& angle_set) override;

  /// Set the current cell to be swept
  void SetCell(Cell const* cell_ptr, AngleSet& angle_set) override;

  /**
   * Performs the discrete ordinates sweep calculation for the currently
   * set cell, for all angles and groups within the provided AngleSet.
   *
   * It:
   * - Assembles the local transport equation system for each angle and group
   * - Retrieves upwind angular fluxes from local neighbors, remote locations
   *   (via MPI data managed by CBC_FLUDS), or boundaries
   * - Solves the local system for the outgoing angular fluxes at the cell nodes
   * - Updates the global scalar flux moments
   * - If save_angular_flux is true, stores the computed angular fluxes into
   *   the global angular flux vector
   * - Propagates outgoing angular fluxes to local downwind neighbors or stages
   *   them for MPI transmission to remote downwind neighbors
   */
  void Sweep(AngleSet& angle_set) override;

protected:
  DiscreteOrdinatesProblem& problem_;
  CBCSweepChunkContext ctx_;
  unsigned int group_block_size_ = 0;

private:
  using SweepFunc = void (CBCSweepChunk::*)(AngleSet&);
  SweepFunc sweep_impl_ = nullptr;

  void Sweep_Generic(AngleSet& angle_set);
  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);
};

} // namespace opensn
