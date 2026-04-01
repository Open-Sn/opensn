// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/avx_sweep_chunk_utils.h"

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

private:
  using SweepFunc = void (CBCSweepChunk::*)(AngleSet&);
  SweepFunc sweep_impl_ = nullptr;

  void Sweep_Generic(AngleSet& angle_set);
  template <unsigned int NumNodes>
  void Sweep_FixedN(AngleSet& angle_set);

  unsigned int group_block_size_;

  CBC_FLUDS* fluds_;
  size_t gs_size_;
  unsigned int gs_gi_;
  size_t num_angles_in_as_;
  unsigned int group_stride_; // Stride for consecutive angles
  size_t group_angle_stride_; // Stride for consecutive spatial DOFs
  bool surface_source_active_;

  const Cell* cell_;
  std::uint32_t cell_local_id_;
  const CellMapping* cell_mapping_;
  CellLBSView* cell_transport_view_;
  size_t cell_num_faces_;
  size_t cell_num_nodes_;

  const DenseMatrix<Vector3>* G_;
  const DenseMatrix<double>* M_;
  const std::vector<DenseMatrix<double>>* M_surf_;
  const std::vector<Vector<double>>* IntS_shapeI_;
};

} // namespace opensn
