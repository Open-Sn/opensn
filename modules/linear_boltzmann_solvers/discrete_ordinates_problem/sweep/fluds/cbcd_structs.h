// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/fluds_structs.h"
#include <stdexcept>
#include <limits>

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
   * \param index Index into the corresponding bank. Must be less than 2^61 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   * \param is_local Flag indicating if the index is in a local bank.
   */
  CBCD_NodeIndex(std::uint64_t index, bool is_outgoing, bool is_local)
  {
    if (index >= (std::uint64_t{1} << 61) - 1)
      throw std::runtime_error("CBCD node index must be less than 2^61 - 1.");
    SetInOut(is_outgoing);
    SetLocal(is_local);
    SetBoundary(false);
    SetIndex(index);
  }

  /**
   * Construct a boundary node index.
   * \param index Index into the corresponding bank. Must be less than 2^61 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   */
  CBCD_NodeIndex(std::uint64_t index, bool is_outgoing)
  {
    if (index >= (std::uint64_t{1} << 61) - 1)
      throw std::runtime_error("CBCD node index must be less than 2^61 - 1.");
    SetInOut(is_outgoing);
    SetLocal(true);
    SetBoundary(true);
    SetIndex(index);
  }

  /// Check if the current index corresponds to a local bank.
  constexpr bool IsLocal() const noexcept { return (value_ & local_bit_mask) != 0; }

  /// Get the index into the bank.
  constexpr std::uint64_t GetIndex() const noexcept { return value_ & index_bit_mask; }

  /**
   * \name Nonlocal nodes
   * The node index packs the cell's nonlocal face index (bits 16–60) and the
   * face-node position (bits 0–15). Outgoing nodes use destination order.
   * \{
   */
  static constexpr unsigned int face_node_bits = std::numeric_limits<std::uint16_t>::digits;
  static_assert(face_node_bits == 16 and std::numeric_limits<std::uint64_t>::digits == 64,
                "Device CBC node indices require uint64 storage and uint16 face nodes.");

  /// Encode an incoming or outgoing node. The face-node position must fit in face_node_bits
  /// and the combined index in 61 bits.
  static CBCD_NodeIndex Nonlocal(std::uint64_t cell_face, std::uint64_t face_node, bool is_outgoing)
  {
    if (face_node >= (std::uint64_t{1} << face_node_bits))
      throw std::runtime_error("CBCD: a nonlocal face has more than 65536 nodes.");
    return {(cell_face << face_node_bits) | face_node, is_outgoing, false};
  }

  constexpr std::uint64_t GetCellFace() const noexcept { return GetIndex() >> face_node_bits; }

  constexpr std::uint64_t GetFaceNode() const noexcept
  {
    return GetIndex() & ((std::uint64_t{1} << face_node_bits) - 1);
  }
  /** \} */

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

/// Offset and boundary flags of one incoming node in the device boundary bank.
struct CBCD_BoundaryInput
{
  std::uint64_t offset = 0;
  bool reflecting = false;
  bool angle_dependent = false;
};

/**
 * Borrowed device views of CBC angular fluxes and per-cell face offsets.
 * Offsets count doubles. Group varies fastest within the angle/group stride.
 * Retain the FLUDS and banks through kernel completion. The kernel supplies the
 * per-cell offset arrays from its launch metadata.
 */
struct CBCD_FLUDSPointerSet : public FLUDSPointerSet
{
  /// Offsets into the incoming bank, indexed by the cell's incoming nonlocal face.
  const std::uint64_t* __restrict__ incoming_payload_offsets = nullptr;
  /// Offsets into the outgoing bank, indexed by the cell's outgoing nonlocal face.
  const std::uint64_t* __restrict__ outgoing_payload_offsets = nullptr;
  /// Device boundary bank offsets and flags indexed by incoming boundary node.
  const CBCD_BoundaryInput* __restrict__ incoming_boundary_map = nullptr;
  /// Device boundary bank offsets indexed by outgoing boundary node.
  const std::uint64_t* __restrict__ outgoing_boundary_map = nullptr;
  /// Reflecting flags indexed by outgoing boundary node (nullptr when no boundary reflects).
  const std::uint8_t* __restrict__ outgoing_boundary_reflecting = nullptr;

  /// Return true for outgoing local, nonlocal and reflecting boundary nodes.
  constexpr bool ShouldStoreOutgoingFlux(const CBCD_NodeIndex& node_index) const noexcept
  {
    if (node_index.IsUndefined() or not node_index.IsOutgoing())
      return false;
    return not node_index.IsBoundary() or
           (outgoing_boundary_reflecting != nullptr and
            outgoing_boundary_reflecting[node_index.GetIndex()] != 0);
  }

  /// Get a borrowed upwind flux pointer, or nullptr if no source is read.
  /// angle_group_idx uses the angle-set stride and group_idx is groupset-relative.
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
      const auto offset =
        incoming_payload_offsets[node_index.GetCellFace()] + node_index.GetFaceNode() * stride_size;
      return nonlocal_incoming_psi + offset + angle_group_idx;
    }
  }

  /// Get a borrowed downwind flux pointer if ShouldStoreOutgoingFlux is true for the node.
  /// angle_group_idx uses the angle-set stride.
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
      const auto offset =
        outgoing_payload_offsets[node_index.GetCellFace()] + node_index.GetFaceNode() * stride_size;
      return nonlocal_outgoing_psi + offset + angle_group_idx;
    }
  }
};

/// Host metadata for one incoming or outgoing boundary face.
struct CBCD_BoundaryFace
{
  std::uint64_t boundary_id = 0;
  std::uint32_t cell_local_id = 0;
  unsigned int face_id = 0;
  /// Index of the face's first node among the incoming or outgoing boundary nodes.
  std::uint32_t storage_offset = 0;
  std::uint16_t num_face_nodes = 0;
};

/// Host metadata for one incoming nonlocal face.
struct CBCD_IncomingNonlocalFace
{
  /// Index of the face among the incoming nonlocal faces in local-cell and face order.
  std::uint32_t face_index = 0;
  std::uint32_t cell_local_id = 0;
  std::uint16_t num_face_nodes = 0;
};

/// Host metadata for one outgoing nonlocal face.
struct CBCD_OutgoingNonlocalFace
{
  /// Index into the destination-partition array.
  std::uint32_t destination_index = 0;
  /// Index of the face among the faces sent to the destination.
  std::uint32_t destination_face_index = 0;
  std::uint16_t num_face_nodes = 0;
};

} // namespace opensn
