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

/// Shared CBCD FLUDS metadata.
class CBCD_FLUDSCommonData : public FLUDSCommonData
{
public:
  /**
   * Construct immutable indexing and communication metadata shared by an SPDS's angle sets.
   *
   * \param spds Cell-by-cell sweep ordering.
   * \param grid_nodal_mappings Face-node mappings for every local cell.
   * \param sdm Spatial discretization defining face-node counts.
   */
  CBCD_FLUDSCommonData(const SPDS& spds,
                       const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
                       const SpatialDiscretization& sdm);

  /// Release the device cell-face-node index.
  ~CBCD_FLUDSCommonData() override;

  /// Return the number of incoming boundary nodes.
  std::size_t GetNumIncomingBoundaryNodes() const { return num_incoming_boundary_nodes_; }

  /// Return the number of outgoing boundary nodes.
  std::size_t GetNumOutgoingBoundaryNodes() const { return num_outgoing_boundary_nodes_; }

  /// Return the number of incoming nonlocal nodes.
  std::size_t GetNumIncomingNonlocalNodes() const { return num_incoming_nonlocal_nodes_; }

  /// Return the number of outgoing nonlocal nodes.
  std::size_t GetNumOutgoingNonlocalNodes() const { return num_outgoing_nonlocal_nodes_; }

  /// Return the number of incoming nonlocal faces.
  std::size_t GetNumIncomingNonlocalFaces() const { return incoming_nonlocal_faces_.size(); }

  /// Return the number of outgoing nonlocal faces.
  std::size_t GetNumOutgoingNonlocalFaces() const { return outgoing_nonlocal_faces_.size(); }

  /// Return grouped incoming-boundary copy plans.
  const std::vector<IncomingBoundaryFacePlan>& GetIncomingBoundaryFaces() const
  {
    return incoming_boundary_face_plans_;
  }

  /// Return outgoing boundary nodes for a local cell.
  std::span<const OutgoingBoundaryNode> GetOutgoingBoundaryNodes(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_outgoing_boundary_node_offsets_[cell_local_id];
    const auto end = cell_to_outgoing_boundary_node_offsets_[cell_local_id + 1];
    return {outgoing_boundary_nodes_.data() + begin, end - begin};
  }

  /// Return outgoing nonlocal faces for a local cell.
  std::span<const OutgoingNonlocalFace> GetOutgoingNonlocalFaces(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_outgoing_nonlocal_face_offsets_[cell_local_id];
    const auto end = cell_to_outgoing_nonlocal_face_offsets_[cell_local_id + 1];
    return {outgoing_nonlocal_faces_.data() + begin, end - begin};
  }

  /// Return incoming nonlocal faces for a local cell.
  std::span<const IncomingNonlocalFace> GetIncomingNonlocalFaces(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_incoming_nonlocal_face_offsets_[cell_local_id];
    const auto end = cell_to_incoming_nonlocal_face_offsets_[cell_local_id + 1];
    return {incoming_nonlocal_faces_.data() + begin, end - begin};
  }

  /// Return the number of local cells represented by the metadata.
  std::size_t GetNumLocalCells() const
  {
    return cell_to_outgoing_nonlocal_face_offsets_.size() - 1;
  }

  /// Return destination MPI ranks indexed by `OutgoingNonlocalFace::destination_index`.
  const std::vector<int>& GetDestinationRanks() const { return destination_ranks_; }

  /// Return source MPI partitions indexed by `IncomingNonlocalFace::source_partition_index`.
  const std::vector<int>& GetIncomingSourcePartitions() const
  {
    return incoming_source_partitions_;
  }

  /// Return an incoming face by source-partition index and source-local face index.
  const IncomingNonlocalFace& GetIncomingNonlocalFace(std::uint32_t source_partition_index,
                                                      std::uint32_t incoming_face_index) const;

  /// Return source-to-destination node mappings for an outgoing face.
  std::span<const OutgoingFaceNodeCopy>
  GetOutgoingFaceNodeCopies(const OutgoingNonlocalFace& face) const
  {
    return {outgoing_nonlocal_face_node_copies_.data() + face.node_copy_begin,
            face.num_node_copies};
  }

  /// Return the device cell-face-node index map.
  const std::uint64_t* GetDeviceIndex() const { return device_cell_face_node_map_; }

private:
  /// Boundary-node counts used to size mapped psi storage.
  std::size_t num_incoming_boundary_nodes_ = 0;
  std::size_t num_outgoing_boundary_nodes_ = 0;
  /// Nonlocal-node counts used to size mapped psi storage.
  std::size_t num_incoming_nonlocal_nodes_ = 0;
  std::size_t num_outgoing_nonlocal_nodes_ = 0;
  /// Device copy of the flattened cell-face-node index map.
  std::uint64_t* device_cell_face_node_map_ = nullptr;
  /// Incoming boundary copies grouped into contiguous face-node ranges.
  std::vector<IncomingBoundaryFacePlan> incoming_boundary_face_plans_;
  /// CSR offsets for outgoing boundary nodes indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_outgoing_boundary_node_offsets_;
  std::vector<OutgoingBoundaryNode> outgoing_boundary_nodes_;
  /// CSR offsets for incoming nonlocal faces indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_incoming_nonlocal_face_offsets_;
  /// CSR offsets for outgoing nonlocal faces indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_outgoing_nonlocal_face_offsets_;
  std::vector<IncomingNonlocalFace> incoming_nonlocal_faces_;
  std::vector<OutgoingNonlocalFace> outgoing_nonlocal_faces_;
  /// Flattened node permutations for outgoing nonlocal faces.
  std::vector<OutgoingFaceNodeCopy> outgoing_nonlocal_face_node_copies_;
  /// Peer arrays referenced by compact indices in face metadata.
  std::vector<int> destination_ranks_;
  std::vector<int> incoming_source_partitions_;
  /// CSR lookup from source-partition index and face index to incoming face metadata.
  std::vector<std::uint32_t> source_to_incoming_face_offsets_;
  std::vector<std::uint32_t> incoming_face_indices_by_source_;

  /// Build host metadata and copy the flattened node index to device storage.
  void BuildMetadataAndCopyNodeIndex(const SpatialDiscretization& sdm);
  /// Release the device node-index allocation.
  void DeallocateDeviceNodeIndex();
};

} // namespace opensn
