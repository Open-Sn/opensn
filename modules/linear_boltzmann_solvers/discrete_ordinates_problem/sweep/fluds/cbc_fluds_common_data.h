// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <cstdint>
#include <cstddef>
#include <unordered_map>

namespace opensn
{

/**
 * Shared CBC FLUDS metadata.
 *
 * Owns the flat face-level lookup tables used to index compact local-face slots and
 * incoming/outgoing non-local face storage for one CBC sweep plane.
 */
class CBC_FLUDSCommonData : public FLUDSCommonData
{
public:
  /// Incoming-face key: `(cell_global_id, face_id)`.
  using CellFaceKey = std::pair<std::uint64_t, unsigned int>;

  /// Hash for CellFaceKey.
  struct CellFaceKeyHash
  {
    size_t operator()(const CellFaceKey& key) const noexcept
    {
      size_t h = std::hash<std::uint64_t>{}(key.first);
      h ^= std::hash<unsigned int>{}(key.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  /// Metadata for one incoming non-local face.
  struct IncomingNonlocalFaceInfo
  {
    /// Local ID of the cell owning this face.
    std::uint32_t cell_local_id = 0;
    /// Offset into the incoming non-local psi buffer for this face's node data.
    std::uint32_t face_node_offset = 0;
    /// Number of face nodes.
    std::uint32_t num_face_nodes = 0;
  };

  /// Metadata for one outgoing non-local face.
  struct OutgoingNonlocalFaceInfo
  {
    /// Destination MPI rank locality index.
    int locality = 0;
    /// Global ID of the destination cell.
    std::uint64_t cell_global_id = 0;
    /// Face index on the destination cell.
    unsigned int associated_face = 0;
    /// Number of face nodes.
    std::uint32_t num_face_nodes = 0;
  };

  /**
   * Construct common data from the SPDS and grid nodal mappings.
   *
   * \param spds Sweep-plane data structure providing face orientations.
   * \param grid_nodal_mappings Per-cell-face nodal mapping data.
   */
  CBC_FLUDSCommonData(const SPDS& spds,
                      const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

  /// Return the number of incoming non-local faces.
  size_t GetNumIncomingNonlocalFaces() const { return num_incoming_nonlocal_faces_; }

  /// Return the number of incoming non-local face nodes.
  size_t GetNumIncomingNonlocalFaceNodes() const { return num_incoming_nonlocal_face_nodes_; }

  /// Return the number of outgoing non-local faces.
  size_t GetNumOutgoingNonlocalFaces() const { return num_outgoing_nonlocal_faces_; }

  /// Return the number of local directed faces.
  size_t GetNumLocalFaces() const { return num_local_faces_; }

  /// Return the maximum local-face node count.
  size_t GetMaxLocalFaceNodeCount() const { return max_local_face_node_count_; }

  /// Return the number of reusable local-face slots.
  size_t GetNumLocalFaceSlots() const { return num_local_face_slots_; }

  /// Get number of outgoing non-local faces for dependent locality `deplocI`.
  size_t GetDeplocIFaceCount(std::size_t deplocI) const noexcept
  {
    return outgoing_nonlocal_face_counts_[deplocI];
  }

  /// Get number of outgoing non-local face nodes for dependent locality `deplocI`.
  size_t GetDeplocIFaceNodeCount(std::size_t deplocI) const noexcept
  {
    return outgoing_nonlocal_face_node_counts_[deplocI];
  }

  /// Look up incoming nonlocal face info by cell local ID and face index.
  const IncomingNonlocalFaceInfo& GetIncomingNonlocalFaceInfo(std::uint32_t cell_local_id,
                                                              unsigned int face_id) const noexcept;

  /// Look up incoming nonlocal face info by flat storage index.
  const IncomingNonlocalFaceInfo&
  GetIncomingNonlocalFaceInfoByStorageIndex(std::size_t storage_index) const noexcept;

  /// Resolve a (cell_global_id, face_id) pair to a flat storage index.
  std::size_t GetIncomingNonlocalFaceStorageIndexByKey(std::uint64_t cell_global_id,
                                                       unsigned int face_id) const noexcept;

  /// Total number of cell-face entries in the flat face table.
  std::size_t GetNumCellFaces() const noexcept { return cell_face_offsets_.back(); }

  /// Look up outgoing nonlocal face info by cell local ID and face index.
  const OutgoingNonlocalFaceInfo& GetOutgoingNonlocalFaceInfo(std::uint32_t cell_local_id,
                                                              unsigned int face_id) const noexcept;

  /// Look up the static local-face slot id by cell local ID and face index.
  std::uint32_t GetLocalFaceSlotID(std::uint32_t cell_local_id, unsigned int face_id) const noexcept
  {
    return local_face_slot_ids_[cell_face_offsets_[cell_local_id] + face_id];
  }

  /// Flat face-table offset for a given cell.
  size_t GetCellFaceOffset(std::uint32_t cell_local_id) const noexcept
  {
    return cell_face_offsets_[cell_local_id];
  }

  /// Flat face-table index for a given cell face.
  size_t GetFaceStorageIndex(std::uint32_t cell_local_id, unsigned int face_id) const noexcept
  {
    return cell_face_offsets_[cell_local_id] + face_id;
  }

  /// Return the flat cell-face offsets table.
  const std::vector<size_t>& GetCellFaceOffsets() const noexcept { return cell_face_offsets_; }

private:
  /// Number of incoming non-local faces.
  size_t num_incoming_nonlocal_faces_;
  /// Number of incoming non-local face nodes.
  size_t num_incoming_nonlocal_face_nodes_;
  /// Number of outgoing non-local faces.
  size_t num_outgoing_nonlocal_faces_;
  /// Number of local directed faces.
  size_t num_local_faces_;
  /// Maximum number of nodes on any local directed face.
  size_t max_local_face_node_count_;
  /// Number of reusable local-face storage slots.
  size_t num_local_face_slots_;
  /// Prefix-sum offsets into the flat face tables, indexed by cell local ID.
  std::vector<size_t> cell_face_offsets_;
  /// Flat local-face slot IDs, indexed by face storage index.
  std::vector<std::uint32_t> local_face_slot_ids_;
  /// Flat incoming non-local face metadata, indexed by face storage index.
  std::vector<IncomingNonlocalFaceInfo> incoming_nonlocal_face_info_;
  /// Flat outgoing non-local face metadata, indexed by face storage index.
  std::vector<OutgoingNonlocalFaceInfo> outgoing_nonlocal_face_info_;
  /// Per-dependent locality outgoing face counts.
  std::vector<size_t> outgoing_nonlocal_face_counts_;
  /// Per-dependent locality outgoing face node counts.
  std::vector<size_t> outgoing_nonlocal_face_node_counts_;
  /// Map from (cell_global_id, face_id) to flat storage index for incoming non-local faces.
  std::unordered_map<CellFaceKey, std::size_t, CellFaceKeyHash> incoming_nonlocal_face_info_by_key_;
};

} // namespace opensn
