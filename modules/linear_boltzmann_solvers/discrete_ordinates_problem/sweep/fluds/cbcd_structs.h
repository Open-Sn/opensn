// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_structs.h"
#include <array>
#include <cstddef>
#include <functional>

namespace opensn
{

class SweepBoundary;

/**
 * Node index specific to CBCD FLUDS.
 * Does not support delayed nodes. Reclaims the delayed bit for indices.
 * - Bit 63: Incoming/outgoing bit.
 * - Bit 62: Boundary bit.
 * - Bit 61: Local bit.
 * - Bits 0-60: Index bits (capacity ~2.3e18).
 */
class CBCD_NodeIndex : public NodeIndex
{
public:
  /// Default constructor.
  constexpr CBCD_NodeIndex() = default;

  /// Direct assign core value.
  constexpr CBCD_NodeIndex(const std::uint64_t& value) : NodeIndex(value) {}

  /**
   * Construct a non-boundary node index.
   * \param index Index into the corresponding bank. Cannot exceed 2^61 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   * \param is_local Flag indicating if the index is in a local bank.
   */
  CBCD_NodeIndex(std::uint64_t index, bool is_outgoing, bool is_local)
  {
    SetInOut(is_outgoing);
    SetLocal(is_local);
    SetBoundary(false);
    SetIndex(index);
  }

  /**
   * Construct a boundary node index.
   * \param index Index into the corresponding bank. Cannot exceed 2^61 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   */
  CBCD_NodeIndex(std::uint64_t index, bool is_outgoing)
  {
    SetInOut(is_outgoing);
    SetLocal(true);
    SetBoundary(true);
    SetIndex(index);
  }

  /// Check if the current index corresponds to a local bank.
  constexpr bool IsLocal() const noexcept { return (value_ & local_bit_mask) != 0; }

  /// Get the index into the bank.
  constexpr std::uint64_t GetIndex() const noexcept { return value_ & index_bit_mask; }

private:
  /// \name Local bit
  /// \{
  /// Third bit mask (``001`` followed by 61 zeros) - Bit 61.
  static constexpr std::uint64_t local_bit_mask = std::uint64_t{1} << (64 - 3);
  /// Encode the value as local.
  constexpr void SetLocal(bool is_local) noexcept
  {
    if (is_local)
      value_ |= local_bit_mask;
    else
      value_ &= ~local_bit_mask;
  }
  /// \}

  /// \name Index bits
  /// \{
  /// Index bit mask (``1`` at the last 61 bits).
  static constexpr std::uint64_t index_bit_mask = (std::uint64_t{1} << (64 - 3)) - 1;
  /// Encode the index.
  constexpr void SetIndex(std::uint64_t index) noexcept
  {
    value_ &= ~index_bit_mask;
    value_ |= (index & index_bit_mask);
  }
  /// \}
};

/**
 * Set of device pointers to local, boundary, and non-local buffers for CBCD FLUDS.
 */
struct CBCD_FLUDSPointerSet : public FLUDSPointerSet
{
  /// Pointer to incoming boundary angular fluxes.
  double* __restrict__ incoming_boundary_psi = nullptr;
  /// Pointer to outgoing boundary angular fluxes.
  double* __restrict__ outgoing_boundary_psi = nullptr;

  /// Get pointer to the incoming angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetIncomingFluxPointer(const CBCD_NodeIndex& node_index,
                                           const unsigned int angle_group_idx) const noexcept
  {
    // Undefined case (corresponds to a parallel face)
    if (node_index.IsUndefined())
      return nullptr;

    // Outgoing case : nullptr
    if (node_index.IsOutgoing())
      return nullptr;

    // Incoming boundary case
    if (node_index.IsBoundary())
    {
      return incoming_boundary_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
    // Incoming local case
    if (node_index.IsLocal())
    {
      return local_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
    // Incoming non-local case
    else
    {
      return nonlocal_incoming_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
  }

  /// Get pointer to the outgoing angular flux (if the face is not outgoing, a nullptr is returned).
  constexpr double* GetOutgoingFluxPointer(const CBCD_NodeIndex& node_index,
                                           const unsigned int angle_group_idx) const noexcept
  {
    // Undefined case (corresponds to a parallel face)
    if (node_index.IsUndefined())
      return nullptr;

    // Incoming case : nullptr
    if (!node_index.IsOutgoing())
      return nullptr;

    // Outgoing boundary case
    if (node_index.IsBoundary())
    {
      return outgoing_boundary_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
    // Outgoing local case
    if (node_index.IsLocal())
    {
      return local_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
    // Outgoing non-local case
    else
    {
      return nonlocal_outgoing_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
  }
};

/// Host metadata for one outgoing boundary node.
struct OutgoingBoundaryNode
{
  /// Boundary identifier.
  std::uint64_t boundary_id = 0;
  /// Local cell identifier.
  std::uint32_t cell_local_id = 0;
  /// Cell-local face identifier.
  unsigned int face_id = 0;
  /// Node offset in outgoing boundary storage.
  std::uint32_t storage_offset = 0;
  /// Face-local node identifier.
  std::uint16_t face_node_index = 0;
};

/// Grouped incoming-boundary face copy plan.
struct IncomingBoundaryFacePlan
{
  /// Boundary identifier.
  std::uint64_t boundary_id = 0;
  /// Local cell identifier.
  std::uint32_t cell_local_id = 0;
  /// Cell-local face identifier.
  unsigned int face_id = 0;
  /// First face-local node in the contiguous group.
  std::uint16_t first_face_node_index = 0;
  /// Node offset in incoming boundary storage.
  std::uint32_t storage_offset = 0;
  /// Number of face nodes in the group.
  std::uint16_t num_face_nodes = 0;
};

/// Receive metadata for one incoming nonlocal face.
struct IncomingNonlocalFace
{
  /// Local downwind cell identifier.
  std::uint32_t cell_local_id = 0;
  /// Node offset in incoming nonlocal storage.
  std::uint32_t storage_offset = 0;
  /// Index into the source-partition array.
  std::uint32_t source_partition_index = 0;
  /// Number of nodes on the face.
  std::uint16_t num_face_nodes = 0;
};

/// Source-to-destination node mapping for an outgoing nonlocal face.
struct OutgoingFaceNodeCopy
{
  /// Node offset in outgoing nonlocal storage.
  std::uint32_t storage_offset = 0;
  /// Node identifier on the destination face.
  std::uint16_t destination_face_node_index = 0;
};

/// Send metadata for one outgoing nonlocal face.
struct OutgoingNonlocalFace
{
  /// Index into the destination-rank array.
  std::uint32_t destination_index = 0;
  /// Face index assigned by the destination for this source rank.
  std::uint32_t destination_face_index = 0;
  /// First entry in the outgoing node-copy array.
  std::uint32_t node_copy_begin = 0;
  /// Number of nodes in the serialized destination face.
  std::uint16_t num_face_nodes = 0;
  /// Number of source-to-destination node copies.
  std::uint16_t num_node_copies = 0;
};

/// Reflecting-boundary face copy plan.
struct ReflectingBoundaryFacePlan
{
  /// Reflecting boundary receiving the outgoing psi.
  SweepBoundary* boundary = nullptr;
  /// Local cell identifier.
  std::uint32_t cell_local_id = 0;
  /// Cell-local face identifier.
  unsigned int face_id = 0;
  /// First face-local node in the contiguous group.
  std::uint16_t first_face_node_index = 0;
  /// Value offset in outgoing boundary storage.
  std::size_t source_offset = 0;
  /// Number of face nodes in the group.
  std::uint16_t num_face_nodes = 0;
};

/// Contiguous psi copy for one outgoing face node.
struct OutgoingPsiCopy
{
  /// Value offset in outgoing nonlocal storage.
  std::size_t source_offset = 0;
  /// Value offset in the serialized destination face.
  std::size_t destination_offset = 0;
};

} // namespace opensn
