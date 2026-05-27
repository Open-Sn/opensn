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

/// Device CBC FLUDS common data.
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

  ~CBCD_FLUDSCommonData() override;

  CBCD_FLUDSCommonData(const CBCD_FLUDSCommonData&) = delete;
  CBCD_FLUDSCommonData& operator=(const CBCD_FLUDSCommonData&) = delete;

  /// Return the number of incoming boundary face nodes.
  std::size_t GetNumIncomingBoundaryNodes() const { return num_incoming_boundary_nodes_; }

  /// Return the number of outgoing boundary face nodes.
  std::size_t GetNumOutgoingBoundaryNodes() const { return num_outgoing_boundary_nodes_; }

  /// Return the number of incoming nonlocal face nodes.
  std::size_t GetNumIncomingNonlocalNodes() const { return num_incoming_nonlocal_nodes_; }

  /// Return the number of outgoing nonlocal face nodes.
  std::size_t GetNumOutgoingNonlocalNodes() const { return num_outgoing_nonlocal_nodes_; }

  /// Return the number of incoming nonlocal faces.
  std::size_t GetNumIncomingNonlocalFaces() const { return incoming_nonlocal_faces_.size(); }

  /// Return the exact incoming face count for one source-partition index.
  std::size_t GetNumIncomingNonlocalFaces(std::size_t source_partition_index) const
  {
    return source_to_incoming_face_offsets_[source_partition_index + 1] -
           source_to_incoming_face_offsets_[source_partition_index];
  }

  /// Return the number of outgoing nonlocal faces.
  std::size_t GetNumOutgoingNonlocalFaces() const { return outgoing_nonlocal_faces_.size(); }

  /// Return incoming boundary faces in local-cell and face order.
  const std::vector<CBCD_IncomingBoundaryFace>& GetIncomingBoundaryFaces() const
  {
    return incoming_boundary_faces_;
  }

  /// Return outgoing boundary faces for a local cell.
  std::span<const CBCD_OutgoingBoundaryFace>
  GetOutgoingBoundaryFaces(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_outgoing_boundary_face_offsets_[cell_local_id];
    const auto end = cell_to_outgoing_boundary_face_offsets_[cell_local_id + 1];
    return {outgoing_boundary_faces_.data() + begin, end - begin};
  }

  /// Return outgoing nonlocal faces for a local cell.
  std::span<const CBCD_OutgoingNonlocalFace>
  GetOutgoingNonlocalFaces(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_outgoing_nonlocal_face_offsets_[cell_local_id];
    const auto end = cell_to_outgoing_nonlocal_face_offsets_[cell_local_id + 1];
    return {outgoing_nonlocal_faces_.data() + begin, end - begin};
  }

  /// Return incoming nonlocal faces for a local cell.
  std::span<const CBCD_IncomingNonlocalFace>
  GetIncomingNonlocalFaces(std::uint64_t cell_local_id) const
  {
    const auto begin = cell_to_incoming_nonlocal_face_offsets_[cell_local_id];
    const auto end = cell_to_incoming_nonlocal_face_offsets_[cell_local_id + 1];
    return {incoming_nonlocal_faces_.data() + begin, end - begin};
  }

  /// Return the number of local cells.
  std::size_t GetNumLocalCells() const
  {
    return cell_to_outgoing_nonlocal_face_offsets_.size() - 1;
  }

  /// Return destination MPI ranks indexed by `CBCD_OutgoingNonlocalFace::destination_index`.
  const std::vector<int>& GetDestinationRanks() const { return destination_ranks_; }

  /// Return source MPI partitions indexed by `CBCD_IncomingNonlocalFace::source_partition_index`.
  const std::vector<int>& GetIncomingSourcePartitions() const
  {
    return incoming_source_partitions_;
  }

  /**
   * Return an incoming nonlocal face by source partition and position.
   * Faces from one source partition are ordered by downwind cell global ID and face index. This is
   * also the order in which the source numbers them.
   *
   * \param source_partition_index Index into the source-partition array.
   * \param incoming_face_index Index of the face among the faces from that source partition.
   * \return Metadata of the face.
   */
  const CBCD_IncomingNonlocalFace& GetIncomingNonlocalFace(std::uint32_t source_partition_index,
                                                           std::uint32_t incoming_face_index) const;

  /// Return the device cell-face-node index map.
  const std::uint64_t* GetDeviceIndex() const { return device_node_indexes_; }

private:
  /// Number of incoming boundary face nodes.
  std::size_t num_incoming_boundary_nodes_ = 0;
  /// Number of outgoing boundary face nodes.
  std::size_t num_outgoing_boundary_nodes_ = 0;
  /// Number of incoming nonlocal face nodes.
  std::size_t num_incoming_nonlocal_nodes_ = 0;
  /// Number of outgoing nonlocal face nodes.
  std::size_t num_outgoing_nonlocal_nodes_ = 0;
  /// Device cell-face-node index map.
  std::uint64_t* device_node_indexes_ = nullptr;
  /// Incoming boundary faces in local-cell and face order.
  std::vector<CBCD_IncomingBoundaryFace> incoming_boundary_faces_;
  /// CSR offsets of the outgoing boundary faces indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_outgoing_boundary_face_offsets_;
  /// Outgoing boundary faces in local-cell and face order.
  std::vector<CBCD_OutgoingBoundaryFace> outgoing_boundary_faces_;
  /// CSR offsets of the incoming nonlocal faces indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_incoming_nonlocal_face_offsets_;
  /// CSR offsets of the outgoing nonlocal faces indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_outgoing_nonlocal_face_offsets_;
  /// Incoming nonlocal faces in local-cell and face order.
  std::vector<CBCD_IncomingNonlocalFace> incoming_nonlocal_faces_;
  /// Outgoing nonlocal faces in local-cell and face order.
  std::vector<CBCD_OutgoingNonlocalFace> outgoing_nonlocal_faces_;
  /// Destination partitions indexed by `CBCD_OutgoingNonlocalFace::destination_index`.
  std::vector<int> destination_ranks_;
  /// Source partitions indexed by `CBCD_IncomingNonlocalFace::source_partition_index`.
  std::vector<int> incoming_source_partitions_;
  /// CSR offsets into incoming_face_indices_by_source_ indexed by source-partition index.
  std::vector<std::uint32_t> source_to_incoming_face_offsets_;
  /// Indices into incoming_nonlocal_faces_, which is grouped by source partition in source order.
  std::vector<std::uint32_t> incoming_face_indices_by_source_;

  /// Build host metadata and copy the flattened node index to device storage.
  void CopyFlattenNodeIndexToDevice(const SpatialDiscretization& sdm);
  /// Release the device node-index allocation.
  void DeallocateDeviceMemory();
};

} // namespace opensn
