// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_structs.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds_common_data.h"
#include <cstdint>
#include <map>

namespace opensn
{

class SpatialDiscretization;

/// Common data for CBCD_FLUDS
class CBCD_FLUDSCommonData : public CBC_FLUDSCommonData
{
public:
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

  /// Get incoming boundary node map.
  const std::vector<BoundaryNodeInfo>& GetIncomingBoundaryNodeMap() const
  {
    return incoming_boundary_node_map_;
  }

  /// Get outgoing boundary node map.
  const std::map<std::uint64_t, std::vector<BoundaryNodeInfo>>& GetOutgoingBoundaryNodeMap() const
  {
    return cell_to_outgoing_boundary_nodes_;
  }

  /// Get incoming nonlocal node map.
  const std::map<std::uint64_t, std::vector<NonlocalNodeInfo>>& GetIncomingNonlocalNodeMap() const
  {
    return cell_to_incoming_nonlocal_nodes_;
  }

  /// Get outgoing nonlocal node map.
  const std::map<std::uint64_t, std::vector<NonlocalNodeInfo>>& GetOutgoingNonlocalNodeMap() const
  {
    return cell_to_outgoing_nonlocal_nodes_;
  }

  /// Get pointer to cell-face-node map on device.
  const std::uint64_t* GetDeviceIndex() const { return device_cell_face_node_map_; }

private:
  /// Number of incoming boundary face nodes.
  size_t num_incoming_boundary_nodes_;
  /// Number of outgoing boundary face nodes.
  size_t num_outgoing_boundary_nodes_;
  /// Number of incoming non-local face nodes.
  size_t num_incoming_nonlocal_nodes_;
  /// Number of outgoing non-local face nodes.
  size_t num_outgoing_nonlocal_nodes_;
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

  /**
   * Compute cell-face-node map for device angular flux buffer access, and
   * create auxiliary indexing maps for boundary and non-local nodes for host access.
   */
  void CopyFlattenedNodeIndexToDevice(const SpatialDiscretization& sdm);
  /// Deallocate device memory for cell-face-node map.
  void DeallocateDeviceMemory();
};

} // namespace opensn