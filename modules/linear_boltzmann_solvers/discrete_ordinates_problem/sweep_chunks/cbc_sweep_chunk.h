// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"

namespace opensn
{
class CellMapping;

/**
 * @class CBCSweepChunk
 * @brief Implements the core sweep operation for a single cell within the
 *        cell-by-cell (CBC) sweep algorithm.
 *
 * This class is responsible for performing the discrete ordinates transport
 * calculation on a given cell for all angles and groups managed by its
 * current `AngleSet`. It interacts with a `CBC_FLUDS` object to obtain
 * upwind angular flux data (from local neighbors, MPI remote buffers, or boundaries)
 * and to store outgoing angular flux data (to local neighbors or MPI send buffers).
 */
class CBCSweepChunk : public SweepChunk
{
public:
  /**
   * @brief Constructs a CBCSweepChunk object.
   * @param destination_phi Reference to the global vector for storing scalar flux moments.
   * @param destination_psi Reference to the global vector for storing angular fluxes
   *                        (written to if `save_angular_flux` is enabled).
   *                        This refers to `LBSSolver::psi_new_local_[groupset_id]`.
   * @param grid Shared pointer to the mesh continuum.
   * @param discretization Reference to the spatial discretization manager.
   * @param unit_cell_matrices Vector of precomputed unit cell matrices.
   * @param cell_transport_views Vector of LBS cell views for transport properties.
   * @param densities Vector of cell material densities.
   * @param source_moments Vector of source moments.
   * @param groupset Reference to the LBS groupset being solved.
   * @param xs Map of material IDs to multi-group cross sections.
   * @param num_moments Number of flux moments.
   * @param max_num_cell_dofs Maximum number of degrees of freedom on any cell
   *                          (used for sizing local solver matrices).
   */
  CBCSweepChunk(std::vector<double>& destination_phi,
                std::vector<double>& destination_psi,
                const std::shared_ptr<MeshContinuum> grid,
                const SpatialDiscretization& discretization,
                const std::vector<UnitCellMatrices>& unit_cell_matrices,
                std::vector<CellLBSView>& cell_transport_views,
                const std::vector<double>& densities,
                const std::vector<double>& source_moments,
                const LBSGroupset& groupset,
                const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
                int num_moments,
                int max_num_cell_dofs);

  /**
   * @brief Sets the current AngleSet for this sweep chunk.
   *
   * This method initializes internal references (like `fluds_`) and calculates
   * necessary strides based on the properties of the given `AngleSet` and its
   * associated `LBSGroupset`. This includes strides for both the compact
   * `local_psi_data_` in `CBC_FLUDS` and for remote communication buffers.
   *
   * @param angle_set Reference to the current `AngleSet` being processed.
   */
  void SetAngleSet(AngleSet& angle_set) override;

  /**
   * @brief Sets the current cell to be processed by the sweep chunk.
   *
   * This method caches pointers and properties related to the `cell_ptr`,
   * such as its local ID, cell mapping, transport view, and precomputed
   * unit cell matrices, to optimize access during the `Sweep` operation.
   *
   * @param cell_ptr Pointer to the constant `Cell` object to be processed.
   * @param angle_set Reference to the current `AngleSet`
   */
  void SetCell(Cell const* cell_ptr, AngleSet& angle_set) override;

  /**
   * @brief Performs the discrete ordinates sweep calculation for the currently
   *        set cell, for all angles and groups within the provided `AngleSet`.
   *
   * This is the core computational method.
   * It:
   * - Assembles the local transport equation system for each angle and group.
   * - Retrieves upwind angular fluxes from local neighbors (via `CBC_FLUDS::local_psi_data_`),
   *   remote locations (via MPI data managed by `CBC_FLUDS`), or boundaries.
   * - Solves the local system for the outgoing angular fluxes at the cell nodes.
   * - Updates the global scalar flux moments (`destination_phi_`).
   * - If `save_angular_flux_` is true, stores the computed angular fluxes into
   *   the global angular flux vector (`destination_psi_`).
   * - Propagates outgoing angular fluxes to local downwind neighbors (by writing
   *   to `CBC_FLUDS::local_psi_data_`) or stages them for MPI transmission
   *   to remote downwind neighbors.
   *
   * @param angle_set Reference to the current `AngleSet` being swept.
   */
  void Sweep(AngleSet& angle_set) override;

private:
  CBC_FLUDS* fluds_; ///< Pointer to the CBC_FLUDS for the current AngleSet.
  size_t gs_size_;   ///< Number of energy groups in the parent LBSGroupset.
  int gs_gi_;        ///< Global starting group index of the parent LBSGroupset.

  /// --- Strides for accessing `local_psi_data_` in `CBC_FLUDS` ---

  /// Number of angular directions managed by the current AngleSet.
  size_t num_angles_in_set_local_;

  /// Stride to jump between data for different angles (local to AngleSet)
  /// for the same spatial DOF within the compact `local_psi_data_`.
  /// Equal to `gs_size_` (number of groups).
  size_t local_compact_angle_stride_;

  /// Stride to jump between data for different spatial DOFs (nodes) within
  /// the compact `local_psi_data_`.
  /// Equal to `num_angles_in_set_local_ * gs_size_`.
  size_t local_compact_node_stride_;

  /// --- Strides for MPI communication buffers (e.g., `psi_dnwnd_data` for outgoing) ---

  /// These buffers are structured to contain data for all angles in THIS AngleSet,
  /// for all nodes on a particular face.
  /// Layout: face spatial DOF major -> angle in set major -> group major.

  /// Number of angular directions in the current AngleSet (same as local, for clarity).
  size_t num_angles_in_set_remote_;

  /// Stride between angles within a single face-node's data block in an MPI buffer.
  /// Equal to `gs_size_`.
  size_t remote_angle_stride_;

  /// Stride between different face-nodes' data blocks in an MPI buffer.
  /// Equal to `num_angles_in_set_remote_ * gs_size_`.
  size_t remote_node_stride_;

  bool surface_source_active_;

  const Cell* cell_;
  uint64_t cell_local_id_;
  const CellMapping* cell_mapping_;
  CellLBSView* cell_transport_view_;
  size_t cell_num_faces_;
  size_t cell_num_nodes_;

  DenseMatrix<Vector3> G_;
  DenseMatrix<double> M_;
  std::vector<DenseMatrix<double>> M_surf_;
  std::vector<Vector<double>> IntS_shapeI_;
};

} // namespace opensn
