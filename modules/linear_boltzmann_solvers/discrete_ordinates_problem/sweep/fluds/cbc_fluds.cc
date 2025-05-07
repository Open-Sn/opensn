// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/math/unknown_manager/unknown_manager.h"

namespace opensn
{

/**
 * @brief Constructs a CBC_FLUDS object.
 *
 * The `local_psi_data_` vector is sized based on the number of angles specific
 * to this `AngleSet` instance, the number of groups in the `LBSGroupset`, and
 * the number of local spatial degrees of freedom.
 * Verification logging is performed to compare the optimized size against the
 * previous, larger sizing approach.
 */
CBC_FLUDS::CBC_FLUDS(size_t num_groups_in_angle_set,
                     size_t num_angles_in_angle_set,
                     const CBC_FLUDSCommonData& common_data,
                     const UnknownManager& lbs_groupset_psi_uk_man,
                     const SpatialDiscretization& sdm)
  : FLUDS(num_groups_in_angle_set, num_angles_in_angle_set, common_data.GetSPDS()),
    common_data_(common_data),
    psi_uk_man_(lbs_groupset_psi_uk_man),
    sdm_(sdm)
{
  /// --- Calculate required size for the optimized `local_psi_data_` ---

  /// 1. Get number of purely spatial DOFs on this MPI rank.
  ///    This is achieved by using a temporary UnknownManager representing a single
  ///    scalar unknown per spatial node.
  const UnknownManager temp_unitary_uk_man({std::make_pair(UnknownType::SCALAR, 0)},
                                           UnknownStorageType::NODAL);
  const size_t num_local_spatial_dofs = sdm_.GetNumLocalDOFs(temp_unitary_uk_man);

  /// 2. Calculate the required size for the `local_psi_data_` vector.
  ///    Layout: spatial DOF major -> angle in set major -> group major
  ///    @note `this->num_angles_` = Number of angles specific to this AngleSet
  ///    @note `this->num_groups_` = Number of groups specific to this AngleSet
  size_t local_psi_data_size = num_local_spatial_dofs * this->num_angles_ * this->num_groups_;

  /// 3. Allocate the local_psi_data_ vector.
  local_psi_data_.assign(local_psi_data_size, 0.0);

  /// --- Verification logging ---
  /// Calculate what the old size would have been for comparison,
  /// when `local_psi_data_` was sized to account for the number of angles
  /// in the groupset's qaudrature

  /// Number of angles in the parent LBSGroupset's full quadrature
  const size_t num_angles_in_gs_quadrature = psi_uk_man_.GetNumberOfUnknowns();

  /// Calculate estimated old size in number of doubles
  size_t size_before_doubles_calc =
    num_local_spatial_dofs * num_angles_in_gs_quadrature * this->num_groups_;

  /// For direct comparison, get size using the LBSGroupset's UnknownManager
  size_t size_before_doubles_direct = sdm.GetNumLocalDOFs(lbs_groupset_psi_uk_man);

  if (local_psi_data_.size() != local_psi_data_size) /// Sanity check allocation
  {
    log.Log0Warning() << "CBC_FLUDS Warning: Allocated local_psi_data_ size ("
                      << local_psi_data_.size() << ") does not match calculated size ("
                      << local_psi_data_size << ").";
  }

  const double mb_divisor = 1024.0 * 1024.0; /// For conversion to MBs
  auto size_before_mb = static_cast<double>(size_before_doubles_calc * sizeof(double)) / mb_divisor;
  auto size_before_direct_mb =
    static_cast<double>(size_before_doubles_direct * sizeof(double)) / mb_divisor;

  auto size_after_mb = static_cast<double>(local_psi_data_size * sizeof(double)) / mb_divisor;

  log.Log() << "CBC_FLUDS Size Comparison for AngleSet (Angles in Set: " << this->num_angles_
            << ", Angles in Groupset Quadrature: " << num_angles_in_gs_quadrature
            << ", Groups: " << this->num_groups_ << "):";
  log.Log() << "  Original estimated size: " << size_before_doubles_calc << " doubles ("
            << std::fixed << std::setprecision(3) << size_before_mb << " MB)";
  log.Log() << "  Original direct size (for verification): " << size_before_doubles_direct
            << " doubles (" << std::fixed << std::setprecision(3) << size_before_direct_mb
            << " MB)";
  log.Log() << "  Optimized size (current):  " << local_psi_data_size << " doubles (" << std::fixed
            << std::setprecision(3) << size_after_mb << " MB)";

  /// Avoid division by zero if original size was effectively 0
  if (size_before_mb > 1e-6)
  {
    double reduction = (size_before_mb - size_after_mb) / size_before_mb * 100.0;
    log.Log() << "  Memory Reduction for local_psi_data_: " << std::fixed << std::setprecision(1)
              << reduction << "%";
  }
}

/**
 * @brief Gets the common data shared among FLUDS instances.
 * @return Constant reference to the `CBC_FLUDSCommonData` object.
 */
const FLUDSCommonData&
CBC_FLUDS::GetCommonData() const
{
  return common_data_;
}

/**
 * @brief Returns a base pointer to an upwind neighbor cell's data block within
 *        the compact `local_psi_data_`.
 *
 * The method calculates the starting memory address for the `face_neighbor` cell's
 * angular flux data. This data block contains fluxes for all nodes of `face_neighbor`,
 * for all angles managed by this `AngleSet`, and for all groups.
 *
 * @param face_neighbor The upwind cell from which data is required.
 * @return `const double*` pointing to the beginning of `face_neighbor`'s data
 *         in `local_psi_data_`. The caller must add relative offsets for specific
 *         nodes and angles.
 * @throws std::runtime_error if the calculated offset is out of bounds.
 */
const double*
CBC_FLUDS::GetLocalUpwindPsi(const Cell& face_neighbor) const
{
  /// Stride to jump from one spatial DOF's full block of (angle-group) data to the next
  /// spatial DOF's block within the compact `local_psi_data_`.
  const size_t node_stride_compact = this->num_angles_ * this->num_groups_;

  /// Use a temporary unitary UnknownManager to get the "flat" spatial DOF index
  /// for the first node of the face_neighbor cell.
  const UnknownManager temp_unitary_uk_man({std::make_pair(UnknownType::SCALAR, 0)},
                                           UnknownStorageType::NODAL);

  /// Get the unique index (0 to num_local_spatial_dofs-1) for the first node of the neighbor cell.
  const int64_t node0_spatial_map =
    sdm_.MapDOFLocal(face_neighbor, 0, temp_unitary_uk_man, 0, 0); // Use local temp manager

  /// Offset to the start of the neighbor cell's data block in the compact local_psi_data_.
  const int64_t offset = node0_spatial_map * node_stride_compact;

  if (offset < 0 || static_cast<size_t>(offset) >= local_psi_data_.size())
  {
    std::ostringstream err_stream;
    err_stream << "CBC_FLUDS::GetLocalUpwindPsi: Offset out of bounds. "
               << "NeighborCell global_id = " << face_neighbor.global_id
               << ", Calculated Offset = " << offset
               << ", CompactVectorSize = " << local_psi_data_.size()
               << ", node0_spatial_map = " << node0_spatial_map
               << ", node_stride_compact = " << node_stride_compact;
    throw std::runtime_error(err_stream.str());
  }

  /// Returns pointer to start of neighbor cell's data block.
  return &local_psi_data_[offset];
}

/**
 * @brief Returns a base pointer to the current cell's data block within
 *        the compact `local_psi_data_` for writing.
 *
 * Similar to `GetLocalUpwindPsi`, but for the cell currently being processed. This
 * pointer allows `CbcSweepChunk` to write the newly computed outgoing angular fluxes.
 *
 * @param cell The current cell being swept.
 * @return `double*` pointing to the beginning of `cell`'s data in `local_psi_data_`.
 *         The caller must add relative offsets for specific nodes and angles.
 * @throws std::runtime_error if the calculated offset is out of bounds.
 */
double*
CBC_FLUDS::GetLocalDownwindPsi(const Cell& cell)
{
  const size_t node_stride_compact = this->num_angles_ * this->num_groups_;

  const UnknownManager temp_unitary_uk_man({std::make_pair(UnknownType::SCALAR, 0)},
                                           UnknownStorageType::NODAL);

  const int64_t node0_spatial_map = sdm_.MapDOFLocal(cell, 0, temp_unitary_uk_man, 0, 0);

  const int64_t offset = node0_spatial_map * node_stride_compact;

  if (offset < 0 || static_cast<size_t>(offset) >= local_psi_data_.size())
  {
    std::ostringstream err_stream;
    err_stream << "CBC_FLUDS::GetLocalDownwindPsi: Offset out of bounds. "
               << "CurrentCell global_id = " << cell.global_id << ", Calculated Offset =" << offset
               << ", CompactVectorSize = " << local_psi_data_.size()
               << ", node0_spatial_map = " << node0_spatial_map
               << ", node_stride_compact = " << node_stride_compact;
    throw std::runtime_error(err_stream.str());
  }
  return &local_psi_data_[offset];
}

/**
 * @brief Retrieves a data packet received via MPI for a specific face of a local cell.
 *
 * This packet contains angular fluxes from a remote upwind cell, covering all angles
 * managed by this `AngleSet` and all nodes on the specified face.
 *
 * @param cell_global_id Global ID of the local cell.
 * @param face_id Local index of the face on the cell.
 * @return Constant reference to the `std::vector<double>` data packet.
 */
const std::vector<double>&
CBC_FLUDS::GetNonLocalUpwindData(uint64_t cell_global_id, unsigned int face_id) const
{
  return deplocs_outgoing_messages_.at({cell_global_id, face_id});
}

/**
 * @brief Extracts a pointer to angular flux data for a specific angle and face node
 *        from a raw MPI data packet.
 *
 * The input `psi_data` is assumed to be a packet containing data for all angles
 * in this `AngleSet`, for all nodes on a face, received via MPI.
 * The layout within `psi_data` is face spatial DOF major -> angle in set major -> group major.
 *
 * @param psi_data The MPI data packet (vector of doubles).
 * @param face_node_mapped The 0-indexed node on the face.
 * @param angle_set_index The local 0-indexed angular direction within this `AngleSet`.
 * @return `const double*` pointing to the start of group data for the specified
 *         `face_node_mapped` and `angle_set_index` within the `psi_data` packet.
 * @throws std::runtime_error if calculated offsets are out of bounds or `psi_data` is empty.
 */
const double*
CBC_FLUDS::GetNonLocalUpwindPsi(const std::vector<double>& psi_data,
                                unsigned int face_node_mapped,
                                unsigned int angle_set_index)
{
  /// Stride to jump from one face node's data block to the next within `psi_data`.
  /// Each face node block contains data for all angles in this AngleSet and all groups.
  const size_t num_psi_per_face_node_for_set = this->num_angles_ * this->num_groups_;

  /// Stride to jump from one angle's data block to the next (for the same face node) within
  /// `psi_data`. Each angle block contains data for all groups.
  const size_t num_groups_stride = this->num_groups_;

  /// Calculate the 1D offset into psi_data.
  const size_t dof_map =
    face_node_mapped *
      num_psi_per_face_node_for_set +    /// Offset to the start of data for this face_node_mapped
    angle_set_index * num_groups_stride; /// Further offset to the start of data for this specific
                                         /// angle_set_index (for group 0)

  if (dof_map + num_groups_stride > psi_data.size() && num_groups_stride > 0)
  {
    std::ostringstream err_stream;
    err_stream << "CBC_FLUDS::GetNonLocalUpwindPsi: Offset out of bounds. "
               << "Calculated_dof_map=" << dof_map << ", num_groups_stride=" << num_groups_stride
               << ", psi_data.size()=" << psi_data.size()
               << ", face_node_mapped=" << face_node_mapped
               << ", angle_set_index=" << angle_set_index
               << ", this->num_angles_=" << this->num_angles_
               << ", this->num_groups_=" << this->num_groups_;
    throw std::runtime_error(err_stream.str());
  }
  if (psi_data.empty() && (face_node_mapped > 0 || angle_set_index > 0))
  {
    throw std::runtime_error(
      "CBC_FLUDS::GetNonLocalUpwindPsi: Accessing non-empty psi_data with non-zero indices.");
  }

  return &psi_data[dof_map];
}

} // namespace opensn
