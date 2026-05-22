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

/// Host CBC FLUDS common data.
class CBC_FLUDSCommonData : public FLUDSCommonData
{
public:
  CBC_FLUDSCommonData(const SPDS& spds,
                      const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

  std::size_t NumIncomingFaces() const { return num_incoming_faces_; }

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
  std::size_t num_incoming_faces_;
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
};

} // namespace opensn
