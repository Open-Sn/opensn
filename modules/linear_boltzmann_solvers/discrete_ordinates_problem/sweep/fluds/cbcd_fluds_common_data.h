// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <cstdint>
#include <map>

namespace opensn
{

class SpatialDiscretization;

/// Device CBC FLUDS common data.
class CBCD_FLUDSCommonData : public FLUDSCommonData
{
public:
  CBCD_FLUDSCommonData(const SPDS& spds,
                       const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
                       const SpatialDiscretization& sdm);

  ~CBCD_FLUDSCommonData() override;

  /// Return the number of incoming boundary face nodes.
  std::size_t GetNumIncomingBoundaryNodes() const { return num_incoming_boundary_nodes_; }

  /// Return the number of outgoing boundary face nodes.
  std::size_t GetNumOutgoingBoundaryNodes() const { return num_outgoing_boundary_nodes_; }

  /// Return the number of incoming nonlocal face nodes.
  std::size_t GetNumIncomingNonlocalNodes() const { return num_incoming_nonlocal_nodes_; }

  /// Return the number of outgoing nonlocal face nodes.
  std::size_t GetNumOutgoingNonlocalNodes() const { return num_outgoing_nonlocal_nodes_; }

  /// Return the number of incoming nonlocal faces.
  std::size_t GetNumIncomingNonlocalFaces() const { return num_incoming_nonlocal_faces_; }

  /// Return the number of outgoing nonlocal faces.
  std::size_t GetNumOutgoingNonlocalFaces() const { return num_outgoing_nonlocal_faces_; }

  /// Return the incoming boundary node map.
  const std::vector<BoundaryNodeInfo>& GetIncomingBoundaryNodeMap() const
  {
    return incoming_boundary_node_map_;
  }

  /// Return the outgoing boundary node map.
  const std::map<std::uint64_t, std::vector<BoundaryNodeInfo>>& GetOutgoingBoundaryNodeMap() const
  {
    return cell_to_outgoing_boundary_nodes_;
  }

  /// Return the incoming nonlocal node map.
  const std::map<std::uint64_t, std::vector<NonlocalNodeInfo>>& GetIncomingNonlocalNodeMap() const
  {
    return cell_to_incoming_nonlocal_nodes_;
  }

  /// Return the outgoing nonlocal node map.
  const std::map<std::uint64_t, std::vector<NonlocalNodeInfo>>& GetOutgoingNonlocalNodeMap() const
  {
    return cell_to_outgoing_nonlocal_nodes_;
  }

  /// Return the device cell-face-node index map.
  const std::uint64_t* GetDeviceIndex() const { return device_cell_face_node_map_; }

private:
  std::size_t num_incoming_boundary_nodes_;
  std::size_t num_outgoing_boundary_nodes_;
  std::size_t num_incoming_nonlocal_faces_;
  std::size_t num_incoming_nonlocal_nodes_;
  std::size_t num_outgoing_nonlocal_faces_;
  std::size_t num_outgoing_nonlocal_nodes_;
  /// Device pointer to cell-face-node map for angular flux buffer access.
  std::uint64_t* device_cell_face_node_map_;
  /// Map from incoming face boundary node to indexing metadata.
  std::vector<BoundaryNodeInfo> incoming_boundary_node_map_;
  /// Map from cell to outgoing boundary nodes.
  std::map<std::uint64_t, std::vector<BoundaryNodeInfo>> cell_to_outgoing_boundary_nodes_;
  /// Map from cell to incoming nonlocal nodes.
  std::map<std::uint64_t, std::vector<NonlocalNodeInfo>> cell_to_incoming_nonlocal_nodes_;
  /// Map from cell to outgoing nonlocal nodes.
  std::map<std::uint64_t, std::vector<NonlocalNodeInfo>> cell_to_outgoing_nonlocal_nodes_;

  /// Compute device and host angular-flux indexing maps.
  void CopyFlattenedNodeIndexToDevice(const SpatialDiscretization& sdm);
  /// Deallocate device memory for cell-face-node map.
  void DeallocateDeviceMemory();
};

} // namespace opensn
