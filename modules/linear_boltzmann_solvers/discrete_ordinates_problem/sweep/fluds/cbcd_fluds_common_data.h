// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include "caribou/main.hpp"
#include <cstdint>
#include <span>
#include <vector>

namespace crb = caribou;

namespace opensn
{

class SpatialDiscretization;

/// Immutable indexing and communication metadata for device CBC FLUDS.
/// The device map stores two uint64 entries per local cell (first face-node entry and count),
/// then CBCD_NodeIndex entries in cell, face and node order.
/// The SPDS and nodal mappings must outlive this object and its FLUDS.
class CBCD_FLUDSCommonData : public FLUDSCommonData
{
public:
  /// Build the metadata and upload the index map during collective runtime setup.
  /// The SPDS's angle sets share this data.
  CBCD_FLUDSCommonData(const SPDS& spds,
                       const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
                       const SpatialDiscretization& sdm);

  /// Return the number of incoming boundary face nodes.
  std::size_t GetNumIncomingBoundaryNodes() const { return num_incoming_boundary_nodes_; }

  /// Return the number of outgoing boundary face nodes.
  std::size_t GetNumOutgoingBoundaryNodes() const { return num_outgoing_boundary_nodes_; }

  /// Return the number of incoming nonlocal face nodes.
  std::size_t GetNumIncomingNonlocalNodes() const { return num_incoming_nonlocal_nodes_; }

  /// Return the number of incoming nonlocal faces.
  std::size_t GetNumIncomingNonlocalFaces() const { return incoming_nonlocal_faces_.size(); }

  /// Return the number of outgoing nonlocal faces.
  std::size_t GetNumOutgoingNonlocalFaces() const { return outgoing_nonlocal_faces_.size(); }

  /// Return the incoming boundary faces in local-cell and face order.
  const std::vector<CBCD_BoundaryFace>& GetIncomingBoundaryFaces() const
  {
    return incoming_boundary_faces_;
  }

  /// Return the outgoing boundary faces in local-cell and face order.
  const std::vector<CBCD_BoundaryFace>& GetOutgoingBoundaryFaces() const
  {
    return outgoing_boundary_faces_;
  }

  /// Return the outgoing nonlocal faces in local-cell and face order.
  const std::vector<CBCD_OutgoingNonlocalFace>& GetOutgoingNonlocalFaces() const
  {
    return outgoing_nonlocal_faces_;
  }

  /// Borrow a host span of a source partition's incoming faces until destruction.
  /// Both peers order these faces by downwind cell global ID and face index.
  std::span<const CBCD_IncomingNonlocalFace>
  GetIncomingNonlocalFacesFromSource(std::size_t source_partition_index) const
  {
    const auto begin = source_to_incoming_face_offsets_[source_partition_index];
    const auto end = source_to_incoming_face_offsets_[source_partition_index + 1];
    return {incoming_nonlocal_faces_.data() + begin, end - begin};
  }

  /// Return the first outgoing face index for a local cell in cell and face order.
  /// Pass the local cell count for the final offset.
  std::uint32_t GetOutgoingNonlocalFaceBegin(std::uint64_t cell_local_id) const
  {
    return cell_to_outgoing_nonlocal_face_offsets_[cell_local_id];
  }

  /// Return the first incoming face index for a local cell in cell and face order.
  /// Pass the local cell count for the final offset.
  std::uint32_t GetIncomingNonlocalFaceBegin(std::uint64_t cell_local_id) const
  {
    return cell_to_incoming_nonlocal_face_offsets_[cell_local_id];
  }

  std::size_t GetNumLocalCells() const
  {
    return cell_to_outgoing_nonlocal_face_offsets_.size() - 1;
  }

  /// Return destination partitions indexed by `CBCD_OutgoingNonlocalFace::destination_index`.
  const std::vector<int>& GetDestinationPartitions() const { return destination_partitions_; }

  const std::vector<int>& GetIncomingSourcePartitions() const
  {
    return incoming_source_partitions_;
  }

  /// Return the device cell-face-node index map.
  const std::uint64_t* GetDeviceIndex() const { return device_node_indexes_.get(); }

private:
  std::size_t num_incoming_boundary_nodes_ = 0;
  std::size_t num_outgoing_boundary_nodes_ = 0;
  std::size_t num_incoming_nonlocal_nodes_ = 0;
  crb::DeviceMemory<std::uint64_t> device_node_indexes_;
  std::vector<CBCD_BoundaryFace> incoming_boundary_faces_;
  std::vector<CBCD_BoundaryFace> outgoing_boundary_faces_;
  /// CSR offsets of the incoming nonlocal faces indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_incoming_nonlocal_face_offsets_;
  /// CSR offsets of the outgoing nonlocal faces indexed by local cell ID.
  std::vector<std::uint32_t> cell_to_outgoing_nonlocal_face_offsets_;
  /// Incoming nonlocal faces grouped by source partition, each group in source order.
  std::vector<CBCD_IncomingNonlocalFace> incoming_nonlocal_faces_;
  /// CSR offsets of incoming_nonlocal_faces_ indexed by source-partition index.
  std::vector<std::uint32_t> source_to_incoming_face_offsets_;
  std::vector<CBCD_OutgoingNonlocalFace> outgoing_nonlocal_faces_;
  std::vector<int> destination_partitions_;
  std::vector<int> incoming_source_partitions_;
};

} // namespace opensn
