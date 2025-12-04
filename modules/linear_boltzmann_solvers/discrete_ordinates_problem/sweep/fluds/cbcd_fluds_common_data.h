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

/**
 * \brief Get pointer to index data of a cell and the number of indexes for that cell.
 * \param device_cell_face_node_map Pointer to the flattened node index structure on device.
 * \param cell_local_idx Cell local index.
 * \return Pointer to the indexes of the cell and the number of indexes (total number of cell face
 * nodes) for the cell.
 */
constexpr std::pair<const std::uint64_t*, std::uint64_t>
GetCBCDCellFaceDataIndex(const std::uint64_t* device_cell_face_node_map,
                         const std::uint32_t& cell_local_idx)
{
  const std::uint64_t* cell_face_data =
    device_cell_face_node_map + static_cast<std::uint64_t>(2 * cell_local_idx);
  return {device_cell_face_node_map + cell_face_data[0], cell_face_data[1]};
}

/**
 * \brief Host-analogous version of GetCBCDCellFaceDataIndex for host copy of the map
 * \param host_cell_face_node_map Reference to the host copy of the flattened node index structure.
 * \param cell_local_idx Cell local index.
 * \return Pointer to the indexes of the cell and the number of indexes (total number of cell face
 * nodes) for the cell.
 */
constexpr std::pair<const std::uint64_t*, std::uint64_t>
GetCBCDCellFaceDataIndexHost(const std::vector<uint64_t>& host_cell_face_node_map,
                             const std::uint32_t& cell_local_idx)
{
  const std::uint64_t* cell_face_data =
    host_cell_face_node_map.data() + static_cast<std::uint64_t>(2 * cell_local_idx);
  return {host_cell_face_node_map.data() + cell_face_data[0], cell_face_data[1]};
}

/**
 * \brief Common data for CBCD_FLUDS (used for encoding/decoding face node indices on host/device)
 */
class CBCD_FLUDSCommonData : public FLUDSCommonData
{
public:
  /**
   * \brief Constructor from SPDS and grid nodal mappings.
   * \details Initializes index trackers for all face nodes of a grid.
   *	        Once created, the indices are utilized by the CBCD_FLUDS to access the correct
   *          location in the flux storage arrays.
   */
  CBCD_FLUDSCommonData(const SPDS& spds,
                       const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
                       const SpatialDiscretization& sdm);

  ~CBCD_FLUDSCommonData() override;

  /// Get constant reference to the spatial discretization.
  SpatialDiscretization const& GetSDM() const { return sdm_; }

  /// Get constant reference to the face node tracker map.
  const std::map<CBCD_FaceNode, CBCD_NodeIndex>& GetNodeTracker() const { return node_tracker_; }

  /// Compute the cell-face-node map for host/device angular flux buffer access.
  void ComputeCellFaceNodeMap();

  /// Get number of incoming boundary face nodes.
  std::size_t GetNumIncomingBoundaryNodes() const { return num_incoming_boundary_nodes_; }

  /// Get number of outgoing boundary face nodes.
  std::size_t GetNumOutgoingBoundaryNodes() const { return num_outgoing_boundary_nodes_; }

  /// Get number of incoming non-local face nodes.
  std::size_t GetNumIncomingNonlocalNodes() const { return num_incoming_nonlocal_nodes_; }

  /// Get number of outgoing non-local face nodes.
  std::size_t GetNumOutgoingNonlocalNodes() const { return num_outgoing_nonlocal_nodes_; }

  /// Get pointer to cell-face-node map on device.
  const std::uint64_t* GetDeviceCellFaceNodeMap() const { return device_cell_face_node_map_; }

  /// Host-analogous version of GetDeviceCellFaceNodeMap for host copy of the map
  const std::vector<uint64_t>& GetHostCellFaceNodeMap() const { return host_cell_face_node_map_; }

protected:
  /// Constant reference to the spatial discretization.
  const SpatialDiscretization& sdm_;

  /// Map from face node to associated index in the corresponding angular flux buffer.
  std::map<CBCD_FaceNode, CBCD_NodeIndex> node_tracker_;

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

  /// Host copy of cell-face-node map for angular flux buffer access.
  std::vector<std::uint64_t> host_cell_face_node_map_;

  /// Constructs flattened node index structure and copy to device.
  void CopyCellFaceNodeMapToDevice();

  /// Deallocate memory for flattened node indices on device.
  void DeallocateDeviceCellFaceNodeMap();
};

} // namespace opensn