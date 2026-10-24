// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_structs.h"
#include <cstddef>
#include <stdexcept>

namespace opensn
{

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
    if (index >= (std::uint64_t{1} << 61) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^61.");
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
    if (index >= (std::uint64_t{1} << 61) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^61.");
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

/// Address of one incoming boundary node in the device boundary bank.
struct CBCD_BoundaryInput
{
  /// Offset of the node's values in the groupset's device boundary bank.
  std::uint64_t offset = 0;
  /// True for a reflecting boundary.
  bool reflecting = false;
  /// True when the boundary values depend on the angle.
  bool angle_dependent = false;
};

/**
 * Set of device pointers to local, boundary, and non-local buffers for CBCD FLUDS.
 */
struct CBCD_FLUDSPointerSet : public FLUDSPointerSet
{
  /// Offsets of the incoming nonlocal nodes' values in the received device packets.
  const std::size_t* __restrict__ incoming_packet_offsets = nullptr;
  /// Device boundary bank addresses indexed by incoming boundary node.
  const CBCD_BoundaryInput* __restrict__ incoming_boundary_map = nullptr;
  /// Device boundary bank offsets indexed by outgoing boundary node.
  const std::uint64_t* __restrict__ outgoing_boundary_map = nullptr;
  /// Reflecting flags indexed by outgoing boundary node (nullptr when no boundary reflects).
  const std::uint8_t* __restrict__ outgoing_boundary_reflecting = nullptr;

  /// Check whether an outgoing face node's angular flux has a downwind consumer.
  constexpr bool ShouldStoreOutgoingFlux(const CBCD_NodeIndex& node_index) const noexcept
  {
    if (node_index.IsUndefined() or not node_index.IsOutgoing())
      return false;
    return not node_index.IsBoundary() or
           (outgoing_boundary_reflecting != nullptr and
            outgoing_boundary_reflecting[node_index.GetIndex()] != 0);
  }

  /// Get pointer to the incoming angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetIncomingFluxPointer(const CBCD_NodeIndex& node_index,
                                           const unsigned int angle_group_idx,
                                           const unsigned int group_idx,
                                           double* boundary,
                                           const bool surface_source_active) const noexcept
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
      const auto& input = incoming_boundary_map[node_index.GetIndex()];
      if (input.reflecting)
        return boundary + input.offset + angle_group_idx;
      if (not surface_source_active)
        return boundary + group_idx;
      return boundary + input.offset + (input.angle_dependent ? angle_group_idx : group_idx);
    }
    // Incoming local case
    if (node_index.IsLocal())
    {
      return local_psi + node_index.GetIndex() * stride_size + angle_group_idx;
    }
    // Incoming non-local case
    else
    {
      const auto node = node_index.GetIndex();
      const auto offset = incoming_packet_offsets[node];
      return nonlocal_incoming_psi + offset + angle_group_idx;
    }
  }

  /// Get pointer to the outgoing angular flux (if the face is not outgoing, a nullptr is returned).
  constexpr double* GetOutgoingFluxPointer(const CBCD_NodeIndex& node_index,
                                           const unsigned int angle_group_idx,
                                           double* boundary) const noexcept
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
      return boundary + outgoing_boundary_map[node_index.GetIndex()] + angle_group_idx;
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

/// Host metadata for one outgoing boundary face.
struct CBCD_OutgoingBoundaryFace
{
  /// Boundary ID.
  std::uint64_t boundary_id = 0;
  /// Face index within the cell.
  unsigned int face_id = 0;
  /// Index of the face's first node among the outgoing boundary nodes.
  std::uint32_t storage_offset = 0;
  /// Number of face nodes.
  std::uint16_t num_face_nodes = 0;
};

/// Host metadata for one incoming boundary face.
struct CBCD_IncomingBoundaryFace
{
  /// Boundary ID.
  std::uint64_t boundary_id = 0;
  /// Local ID of the cell.
  std::uint32_t cell_local_id = 0;
  /// Face index within the cell.
  unsigned int face_id = 0;
  /// Index of the face's first node among the incoming boundary nodes.
  std::uint32_t storage_offset = 0;
  /// Number of face nodes.
  std::uint16_t num_face_nodes = 0;
};

/// Host metadata for one incoming nonlocal face.
struct CBCD_IncomingNonlocalFace
{
  /// Local ID of the downwind cell.
  std::uint32_t cell_local_id = 0;
  /// Index of the face's first node among the incoming nonlocal nodes.
  std::uint32_t storage_offset = 0;
  /// Index into the source-partition array.
  std::uint32_t source_partition_index = 0;
  /// Number of face nodes.
  std::uint16_t num_face_nodes = 0;
};

/// Host metadata for one outgoing nonlocal face.
struct CBCD_OutgoingNonlocalFace
{
  /// Index into the destination-rank array.
  std::uint32_t destination_index = 0;
  /// Index of the face among the faces sent to the destination.
  std::uint32_t destination_face_index = 0;
  /// Index of the face's first node among the outgoing nonlocal nodes, in destination-node order.
  std::uint32_t storage_offset = 0;
  /// Number of face nodes.
  std::uint16_t num_face_nodes = 0;
};

} // namespace opensn
