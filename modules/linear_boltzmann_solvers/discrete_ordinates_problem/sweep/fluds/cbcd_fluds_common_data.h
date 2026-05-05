// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <cstdint>
#include <span>
#include <vector>

namespace opensn
{

class SpatialDiscretization;

/**
 * Shared CBCD FLUDS metadata.
 *
 * Builds and owns the flattened indexing tables used by every CBCD FLUDS instance
 * associated with one SPDS. The tables translate cell-face-node accesses into
 * compact local, boundary, and non-local storage indices on both the host and device.
 */
class CBCD_FLUDSCommonData : public FLUDSCommonData
{
public:
  /**
   * Construct the shared CBCD FLUDS metadata for one SPDS.
   *
   * \param spds Sweep plane data structure providing the CBC cell and face ordering.
   * \param grid_nodal_mappings Per-cell face-node mappings from the spatial discretization.
   * \param sdm Spatial discretization used to enumerate face nodes.
   */
  CBCD_FLUDSCommonData(const SPDS& spds,
                       const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
                       const SpatialDiscretization& sdm);

  ~CBCD_FLUDSCommonData() override;

  /// Get number of incoming boundary face nodes.
  std::size_t GetNumIncomingBoundaryNodes() const { return num_incoming_boundary_nodes_; }

  /// Get number of outgoing boundary face nodes.
  std::size_t GetNumOutgoingBoundaryNodes() const { return num_outgoing_boundary_nodes_; }

  /// Get number of incoming non-local face nodes.
  std::size_t GetNumIncomingNonlocalNodes() const { return num_incoming_nonlocal_nodes_; }

  /// Get number of outgoing non-local face nodes.
  std::size_t GetNumOutgoingNonlocalNodes() const { return num_outgoing_nonlocal_nodes_; }

  /// Get number of incoming non-local faces.
  std::size_t GetNumIncomingNonlocalFaces() const { return num_incoming_nonlocal_faces_; }

  /// Get number of outgoing non-local faces.
  std::size_t GetNumOutgoingNonlocalFaces() const { return num_outgoing_nonlocal_faces_; }

  /// Return grouped incoming-boundary faces.
  const std::vector<IncomingBoundaryFacePlan>& GetIncomingBoundaryFaces() const
  {
    return incoming_boundary_face_plans_;
  }

  /// Return the number of grouped incoming non-local faces from one source locality slot.
  std::size_t GetNumIncomingFacesFromSource(const std::size_t source_slot) const
  {
    return source_to_incoming_face_offsets_[source_slot + 1] -
           source_to_incoming_face_offsets_[source_slot];
  }

  /// Return outgoing-boundary nodes for one cell.
  std::span<const BoundaryNodeInfo> GetOutgoingBoundaryNodes(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_outgoing_boundary_node_offsets_[cell_local_id];
    const auto end = cell_to_outgoing_boundary_node_offsets_[cell_local_id + 1];
    return {outgoing_boundary_nodes_.data() + begin, end - begin};
  }

  /// Return grouped outgoing non-local faces for one cell.
  std::span<const GroupedOutgoingNonlocalFace>
  GetOutgoingNonlocalFaces(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_outgoing_nonlocal_face_offsets_[cell_local_id];
    const auto end = cell_to_outgoing_nonlocal_face_offsets_[cell_local_id + 1];
    return {outgoing_nonlocal_faces_.data() + begin, end - begin};
  }

  /// Return grouped incoming non-local faces for one cell.
  std::span<const GroupedIncomingNonlocalFace>
  GetIncomingNonlocalFaces(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_incoming_nonlocal_face_offsets_[cell_local_id];
    const auto end = cell_to_incoming_nonlocal_face_offsets_[cell_local_id + 1];
    return {incoming_nonlocal_faces_.data() + begin, end - begin};
  }

  /// Return the number of local cells represented in the grouped-face tables.
  std::size_t GetNumLocalCells() const
  {
    return cell_to_outgoing_nonlocal_face_offsets_.size() - 1;
  }

  /// Return the ordered outgoing-locality table used to build communicator queue indices.
  const std::vector<int>& GetOutgoingLocalities() const { return outgoing_localities_; }

  /// Return the ordered incoming source-locality table.
  const std::vector<int>& GetIncomingSourcePartitions() const
  {
    return incoming_source_partitions_;
  }

  /// Resolve one grouped incoming non-local face by source-slot-local face index.
  const GroupedIncomingNonlocalFace& GetIncomingNonlocalFace(std::uint32_t source_slot,
                                                             std::uint32_t source_face_index) const;

  /// Return the outgoing-node-copy descriptors for one grouped outgoing face.
  std::span<const OutgoingNodeCopy>
  GetOutgoingNodeCopies(const GroupedOutgoingNonlocalFace& face) const
  {
    return {outgoing_nonlocal_face_node_copies_.data() + face.node_copy_offset,
            face.num_node_copies};
  }

  /// Get pointer to cell-face-node map on device.
  const std::uint64_t* GetDeviceIndex() const { return device_cell_face_node_map_; }

private:
  /// Number of incoming boundary face nodes.
  size_t num_incoming_boundary_nodes_;
  /// Number of outgoing boundary face nodes.
  size_t num_outgoing_boundary_nodes_;
  /// Number of incoming non-local faces.
  size_t num_incoming_nonlocal_faces_;
  /// Number of incoming non-local face nodes.
  size_t num_incoming_nonlocal_nodes_;
  /// Number of outgoing non-local faces.
  size_t num_outgoing_nonlocal_faces_;
  /// Number of outgoing non-local face nodes.
  size_t num_outgoing_nonlocal_nodes_;
  /// Device pointer to cell-face-node map for angular flux buffer access.
  std::uint64_t* device_cell_face_node_map_;
  /// Flat grouped incoming-boundary face copy plans.
  std::vector<IncomingBoundaryFacePlan> incoming_boundary_face_plans_;
  /// Cell-to-outgoing-boundary-node offset table.
  std::vector<std::uint32_t> cell_to_outgoing_boundary_node_offsets_;
  /// Flat outgoing-boundary node list.
  std::vector<BoundaryNodeInfo> outgoing_boundary_nodes_;
  /// Cell-to-incoming-face offset table.
  std::vector<std::uint32_t> cell_to_incoming_nonlocal_face_offsets_;
  /// Cell-to-outgoing-face offset table.
  std::vector<std::uint32_t> cell_to_outgoing_nonlocal_face_offsets_;
  /// Flat grouped incoming nonlocal faces.
  std::vector<GroupedIncomingNonlocalFace> incoming_nonlocal_faces_;
  /// Flat grouped outgoing nonlocal faces.
  std::vector<GroupedOutgoingNonlocalFace> outgoing_nonlocal_faces_;
  /// Flat outgoing-node-copy metadata referenced by grouped outgoing faces.
  std::vector<OutgoingNodeCopy> outgoing_nonlocal_face_node_copies_;
  /// Ordered table of distinct outgoing localities.
  std::vector<int> outgoing_localities_;
  /// Ordered table of incoming source localities.
  std::vector<int> incoming_source_partitions_;
  /// Source-major incoming grouped-face spans.
  std::vector<std::uint32_t> source_to_incoming_face_offsets_;
  /// Source-major ordered incoming grouped-face indices.
  std::vector<std::uint32_t> incoming_face_indices_by_source_;

  /**
   * Build and upload the flattened cell-face-node index map.
   *
   * \param sdm Spatial discretization used to enumerate face nodes.
   */
  void CopyFlattenedNodeIndexToDevice(const SpatialDiscretization& sdm);
  /// Deallocate device memory for cell-face-node map.
  void DeallocateDeviceMemory();
};

} // namespace opensn
