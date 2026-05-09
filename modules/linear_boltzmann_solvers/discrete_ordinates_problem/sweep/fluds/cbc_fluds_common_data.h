// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace opensn
{

/// Host CBC FLUDS common data shared by angle-set subsets.
class CBC_FLUDSCommonData : public FLUDSCommonData
{
public:
  /// Slot metadata for a lagged local incoming face.
  struct DelayedLocalFaceInfo
  {
    /// Offset into the delayed-local face-node storage.
    std::size_t slot_address = 0;
    /// Number of nodes on the delayed local face.
    std::size_t num_face_nodes = 0;
  };

  /// Slot metadata for a lagged nonlocal incoming face.
  struct DelayedNonlocalFaceInfo
  {
    /// Delayed upstream-location index in the SPDS dependency list.
    std::size_t prelocI = 0;
    /// Offset into the delayed nonlocal face-node storage for `prelocI`.
    std::size_t slot_address = 0;
    /// Number of nodes on the delayed nonlocal face.
    std::size_t num_face_nodes = 0;
  };

  CBC_FLUDSCommonData(const SPDS& spds,
                      const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

  /// Return the number of non-delayed incoming nonlocal faces.
  std::size_t NumIncomingFaces() const { return num_incoming_faces_; }

  /// Return the number of delayed incoming nonlocal faces.
  std::size_t NumDelayedNonlocalFaces() const { return delayed_nonlocal_face_info_by_slot_.size(); }

  /// Return the total number of delayed local face nodes.
  std::size_t NumDelayedLocalFaceNodes() const { return num_delayed_local_face_nodes_; }

  /// Return the delayed nonlocal face-node count for a delayed upstream location.
  std::size_t DelayedPrelocIFaceNodeCount(std::size_t prelocI) const;

  /// Return the face-node count for a delayed nonlocal face slot.
  std::size_t DelayedNonlocalFaceNodeCount(std::size_t delayed_face_slot) const;

  /// Return whether a local face reads lagged local incoming face psi.
  bool IsDelayedLocalIncomingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return whether a local face writes lagged local outgoing face psi.
  bool IsDelayedLocalOutgoingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return whether a local face reads lagged nonlocal incoming face psi.
  bool IsDelayedNonlocalIncomingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return whether a local face writes lagged nonlocal outgoing face psi.
  bool IsDelayedNonlocalOutgoingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return delayed local face metadata for a local cell face.
  const DelayedLocalFaceInfo& DelayedLocalFace(std::uint32_t cell_local_id,
                                               unsigned int face_id) const;

  /// Return delayed nonlocal face metadata for a local cell face.
  const DelayedNonlocalFaceInfo& DelayedNonlocalFaceByLocalFace(std::uint32_t cell_local_id,
                                                                unsigned int face_id) const;

  /// Return delayed nonlocal face metadata for a delayed face slot.
  const DelayedNonlocalFaceInfo& DelayedNonlocalFaceBySlot(std::size_t delayed_face_slot) const;

  /**
   * Return the incoming nonlocal face slot for a local face.
   *
   * \param cell_local_id Local cell ID.
   * \param face_id Local face ID.
   * \return Incoming face slot or `INVALID_FACE_SLOT`.
   */
  std::size_t IncomingFaceSlot(std::uint32_t cell_local_id, unsigned int face_id) const;

  /**
   * Return the downstream incoming face slot for an outgoing local face.
   *
   * \param cell_local_id Local cell ID.
   * \param face_id Local face ID.
   * \return Downstream incoming face slot or `INVALID_FACE_SLOT`.
   */
  std::size_t OutgoingFaceSlot(std::uint32_t cell_local_id, unsigned int face_id) const;

  /**
   * Return the SPDS-successor peer index for an outgoing local face.
   *
   * \param cell_local_id Local cell ID.
   * \param face_id Local face ID.
   * \return Successor peer index or `INVALID_PEER_INDEX`.
   */
  std::size_t OutgoingPeerIndex(std::uint32_t cell_local_id, unsigned int face_id) const;

  /**
   * Return the destination location for an outgoing local face.
   *
   * \param cell_local_id Local cell ID.
   * \param face_id Local face ID.
   * \return Destination MPI location or -1.
   */
  int OutgoingFaceLocation(std::uint32_t cell_local_id, unsigned int face_id) const;

  /**
   * Return the local cell associated with an incoming face slot.
   *
   * \param incoming_face_slot Incoming face slot.
   * \return Local cell ID associated with the incoming face slot.
   */
  std::uint32_t IncomingFaceCell(std::size_t incoming_face_slot) const;

  /// Marker for cell faces without a nonlocal slot.
  static constexpr std::size_t INVALID_FACE_SLOT = std::numeric_limits<std::size_t>::max();

  /// Marker for cell faces without an outgoing nonlocal peer.
  static constexpr std::size_t INVALID_PEER_INDEX = std::numeric_limits<std::size_t>::max();

private:
  /// Number of non-delayed incoming nonlocal faces.
  std::size_t num_incoming_faces_;
  /// Total number of delayed local face nodes.
  std::size_t num_delayed_local_face_nodes_;
  /// Prefix offsets into local-face-indexed slot arrays.
  std::vector<std::size_t> face_offsets_;
  /// Local-face-indexed incoming slots.
  std::vector<std::size_t> incoming_face_slots_;
  /// Slot-indexed local cells associated with received face psi.
  std::vector<std::uint32_t> incoming_face_cells_;
  /// Local-face-indexed downstream incoming slots.
  std::vector<std::size_t> outgoing_face_slots_;
  /// Local-face-indexed SPDS-successor peer indices.
  std::vector<std::size_t> outgoing_peer_indices_;
  /// Local-face-indexed outgoing destination MPI locations.
  std::vector<int> outgoing_face_locations_;
  /// Local-face-indexed delayed local face metadata.
  std::vector<DelayedLocalFaceInfo> delayed_local_face_info_by_face_;
  /// Slot-indexed delayed nonlocal face metadata.
  std::vector<DelayedNonlocalFaceInfo> delayed_nonlocal_face_info_by_slot_;
  /// Local-face-indexed delayed nonlocal face metadata.
  std::vector<DelayedNonlocalFaceInfo> delayed_nonlocal_face_info_by_face_;
  /// Delayed nonlocal face-node counts by delayed upstream-location index.
  std::vector<std::size_t> delayed_prelocI_face_node_counts_;
  /// Local-face-indexed flags for delayed local incoming faces.
  std::vector<unsigned char> delayed_local_incoming_faces_;
  /// Local-face-indexed flags for delayed local outgoing faces.
  std::vector<unsigned char> delayed_local_outgoing_faces_;
  /// Local-face-indexed flags for delayed nonlocal incoming faces.
  std::vector<unsigned char> delayed_nonlocal_incoming_faces_;
  /// Local-face-indexed flags for delayed nonlocal outgoing faces.
  std::vector<unsigned char> delayed_nonlocal_outgoing_faces_;
};

} // namespace opensn
