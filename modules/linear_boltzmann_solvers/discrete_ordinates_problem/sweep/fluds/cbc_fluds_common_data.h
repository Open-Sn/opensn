// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_common_data.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

namespace opensn
{

/**
 * Shared face-slot metadata for host CBC angular-flux storage.
 * The SPDS and mesh must outlive this object. Finish beta setup before constructing
 * FLUDS or querying nonlocal slots.
 */
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

  /// Construct alpha phase only (local work, no MPI). Call FinalizeBeta() afterwards.
  static std::unique_ptr<CBC_FLUDSCommonData>
  MakeAlpha(const SPDS& spds, const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

  /**
   * Exchange nonlocal face slots and finalize the metadata.
   *
   * Collective on mpi_comm. All ranks call in identical SPDS order after local
   * construction. A second call is rejected.
   *
   * \throw std::logic_error If already finalized, an exchange count is too large, or received
   * face-slot records are malformed or inconsistent with the sweep graph.
   */
  void FinalizeBeta();

  /**
   * Initialize beta elements from exchanged receiving-face slots, without MPI calls.
   *
   * Group records by downstream location using four uint64_t fields per face (receiving
   * cell global ID, face index, receiving slot and a 0 or 1 delayed flag). Include each
   * outgoing nonlocal face exactly once, with the SPDS delayed classification. A failed
   * validation leaves the object unfinalized.
   *
   * \param downstream_slot_records Complete receiving-face records from every downstream location.
   * \throw std::logic_error If already finalized, or a record count, endpoint, slot, flag or
   * delayed classification is invalid, missing or repeated.
   */
  void
  InitializeBetaElements(const std::map<int, std::vector<std::uint64_t>>& downstream_slot_records);

  /**
   * Require finalized face metadata before allocating angular values.
   *
   * \throw std::logic_error If beta setup has not completed.
   */
  void CheckFinalized() const
  {
    if (not finalized_)
      throw std::logic_error("CBC FLUDS common data must be finalized before use");
  }

  /// Return the number of non-delayed incoming nonlocal faces.
  std::size_t NumIncomingFaces() const { return num_incoming_faces_; }

  /// Return the number of delayed incoming nonlocal face slots.
  std::size_t GetNumDelayedIncomingFaces() const
  {
    return delayed_nonlocal_face_info_by_slot_.size();
  }

  /// Return the total number of delayed local face nodes.
  std::size_t GetNumDelayedLocalFaceNodes() const { return num_delayed_local_face_nodes_; }

  /// Return the delayed nonlocal face-node count for a delayed upstream location.
  std::size_t GetDelayedPrelocIFaceNodeCount(std::size_t prelocI) const;

  /// Return whether a local face reads lagged local incoming face psi.
  bool IsDelayedLocalIncomingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return whether a local face writes lagged local outgoing face psi.
  bool IsDelayedLocalOutgoingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return whether a local face reads lagged nonlocal incoming face psi.
  bool IsDelayedNonlocalIncomingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return whether a local face writes lagged nonlocal outgoing face psi.
  bool IsDelayedNonlocalOutgoingFace(std::uint32_t cell_local_id, unsigned int face_id) const;

  /// Return delayed local face metadata for a local cell face.
  const DelayedLocalFaceInfo& GetDelayedLocalFace(std::uint32_t cell_local_id,
                                                  unsigned int face_id) const;

  /// Return delayed nonlocal face metadata for a local cell face.
  const DelayedNonlocalFaceInfo& GetDelayedNonlocalFaceByLocalFace(std::uint32_t cell_local_id,
                                                                   unsigned int face_id) const;

  /// Return delayed nonlocal face metadata for a delayed face slot.
  const DelayedNonlocalFaceInfo& GetDelayedNonlocalFaceBySlot(std::size_t delayed_face_slot) const;

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
   * Return the ordinary or delayed SPDS-successor index for an outgoing nonlocal face.
   * Select the successor list using IsDelayedNonlocalOutgoingFace.
   *
   * \param cell_local_id Local cell ID.
   * \param face_id Local face ID.
   * \return Index in the selected successor list, or `INVALID_PEER_INDEX`.
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
  enum DelayedFaceFlag : std::uint8_t
  {
    LOCAL_INCOMING = 1U << 0U,
    LOCAL_OUTGOING = 1U << 1U,
    NONLOCAL_INCOMING = 1U << 2U,
    NONLOCAL_OUTGOING = 1U << 3U
  };

  CBC_FLUDSCommonData(const SPDS& spds,
                      const std::vector<CellFaceNodalMapping>& grid_nodal_mappings);

  /// Number of ordinary incoming nonlocal faces.
  std::size_t num_incoming_faces_;
  /// Number of outgoing nonlocal faces, including delayed faces.
  std::size_t num_outgoing_faces_;
  /// Whether the nonlocal face-slot exchange has completed.
  bool finalized_ = false;
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
  /// Local-face-indexed ordinary or delayed SPDS-successor peer indices.
  std::vector<std::size_t> outgoing_peer_indices_;
  /// Incoming face-slot records grouped by upstream location.
  std::map<int, std::vector<std::uint64_t>> incoming_slot_records_by_upstream_location_;
  /// Outgoing face counts used to size the setup exchange by downstream location.
  std::map<int, std::size_t> outgoing_face_counts_by_location_;
  /// Local-face-indexed delayed local face metadata.
  std::vector<DelayedLocalFaceInfo> delayed_local_face_info_by_face_;
  /// Slot-indexed delayed nonlocal face metadata.
  std::vector<DelayedNonlocalFaceInfo> delayed_nonlocal_face_info_by_slot_;
  /// Local-face-indexed delayed nonlocal face metadata.
  std::vector<DelayedNonlocalFaceInfo> delayed_nonlocal_face_info_by_face_;
  /// Delayed nonlocal face-node counts by delayed upstream-location index.
  std::vector<std::size_t> delayed_prelocI_face_node_counts_;
  /// Bit-packed delayed-face classifications indexed by local face.
  std::vector<std::uint8_t> delayed_face_flags_;
};

} // namespace opensn
