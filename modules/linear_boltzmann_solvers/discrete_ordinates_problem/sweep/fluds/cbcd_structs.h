// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>

namespace opensn
{

/**
 * \brief Represents a cell face node.
 * \details This class represents a cell face node. It stores the cell local index, the face index
 * within the cell, and the index of the node relative to the face as a compact 64-bit integer. This
 * class abstracts the node for usage with std::map and std::set.
 */
class CBCD_FaceNode
{
public:
  /// Default constructor.
  constexpr CBCD_FaceNode() = default;

  /// Member constructor.
  constexpr CBCD_FaceNode(std::uint32_t cell_idx,
                          std::uint16_t face_idx,
                          std::uint16_t face_node_idx)
    : value_(0)
  {
    value_ |= static_cast<std::uint64_t>(cell_idx) << 32;
    value_ |= static_cast<std::uint64_t>(face_idx) << 16;
    value_ |= static_cast<std::uint64_t>(face_node_idx);
  }

  /// Comparison operator for ordering.
  constexpr bool operator<(const CBCD_FaceNode& other) const { return value_ < other.value_; }

  /// Equality operator.
  constexpr bool operator==(const CBCD_FaceNode& other) const { return value_ == other.value_; }

  /// Get cell local index.
  constexpr std::uint32_t GetCellIndex() const { return static_cast<std::uint32_t>(value_ >> 32); }

  /// Get face index.
  constexpr std::uint16_t GetFaceIndex() const
  {
    return static_cast<std::uint16_t>((value_ >> 16) & 0xFFFFU);
  }

  /// Get face node index.
  constexpr std::uint16_t GetFaceNodeIndex() const
  {
    return static_cast<std::uint16_t>(value_ & 0xFFFFU);
  }

  /// Check if face node is initialized.
  constexpr bool IsInitialized() const
  {
    return value_ != std::numeric_limits<std::uint64_t>::max();
  }

private:
  /// Core value.
  std::uint64_t value_ = std::numeric_limits<std::uint64_t>::max();
};

/**
 * \brief Class that represents a 64-bit integer encoding the index and the host/device information
 *        of a face node.
 * \note The index follows the following mutually exclusive rules:
 * 	     - Boundary status is checked BEFORE local status. Although boundary nodes are treated as
 *         local nodes, their locations are set to an undefined value.
 * 	     - A node cannot be both a boundary node and a non-local node.
 */
class CBCD_NodeIndex
{
public:
  /// Default constructor
  CBCD_NodeIndex() = default;

  /// Direct assign core value
  constexpr CBCD_NodeIndex(const std::uint64_t& value) : value_(value) {}

  /**
   * Construct a non-boundary node index.
   * \param index Index into the corresponding bank. Cannot exceed 2^60 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   * \param is_local Flag indicating if the index is in a local bank.
   */
  CBCD_NodeIndex(std::uint64_t index, bool is_outgoing, bool is_local)
  {
    if (index >= (std::uint64_t(1) << 60) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^60.");
    SetInOut(is_outgoing);
    SetLocal(is_local);
    SetBoundary(false);
    SetIndex(index);
  }

  /**
   * Construct a boundary node index.
   * \param index Index into the corresponding bank. Cannot exceed 2^40 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   */
  CBCD_NodeIndex(std::uint64_t index, bool is_outgoing)
  {
    if (index >= (std::uint64_t(1) << 60) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^60.");
    SetInOut(is_outgoing);
    SetLocal(true);
    SetBoundary(true);
    SetIndex(index);
  }

  /// Check if the current node's index is undefined.
  constexpr bool IsUndefined() const noexcept
  {
    return value_ == std::numeric_limits<std::uint64_t>::max();
  }

  /// Check if the node corresponds to an outgoing face.
  constexpr bool IsOutgoing() const noexcept { return (value_ & inout_bit_mask) != 0; }

  /// Check if the node corresponds to a boundary face.
  constexpr bool IsBoundary() const noexcept { return (value_ & boundary_bit_mask) != 0; }

  /// Check if the node index corresponds to a local face.
  constexpr bool IsLocal() const noexcept { return (value_ & local_bit_mask) != 0; }

  /// Get the node index into the appropriate host/device buffer.
  constexpr std::uint64_t GetIndex() const noexcept { return value_ & index_bit_mask; }

  /// Get the node's core value (can be used with the node tracker map to retrieve a given node's
  /// index into the appropriate host/device buffer).
  constexpr std::uint64_t GetCoreValue() const noexcept { return value_; }

private:
  /// Incoming/outgoing bit
  /// First bit mask (``1`` followed by 63 zeros).
  static constexpr std::uint64_t inout_bit_mask = std::uint64_t(1) << (64 - 1);

  /// Encode the node as incoming or outgoing.
  constexpr void SetInOut(bool is_outgoing) noexcept
  {
    if (is_outgoing)
      value_ |= inout_bit_mask;
    else
      value_ &= ~inout_bit_mask;
  }

  /// Boundary bit
  /// Second bit mask (``01`` followed by 62 zeros).
  static constexpr std::uint64_t boundary_bit_mask = std::uint64_t(1) << (64 - 2);

  /// Encode the node as either a boundary node or non-boundary node.
  constexpr void SetBoundary(bool is_boundary) noexcept
  {
    if (is_boundary)
      value_ |= boundary_bit_mask;
    else
      value_ &= ~boundary_bit_mask;
  }

  /// Local bit
  /// Third bit mask (``001`` followed by 61 zeros).
  static constexpr std::uint64_t local_bit_mask = std::uint64_t(1) << (64 - 3);

  /// Encode the node as either a local or non-local node.
  constexpr void SetLocal(bool is_local) noexcept
  {
    if (is_local)
      value_ |= local_bit_mask;
    else
      value_ &= ~local_bit_mask;
  }

  // Index bits
  // Index bit mask (``1`` at the last 61 bits)
  static constexpr std::uint64_t index_bit_mask = (std::uint64_t(1) << (64 - 3)) - 1;

  /// Encode the node index.
  constexpr void SetIndex(const std::uint64_t& index) noexcept
  {
    value_ &= ~index_bit_mask;
    value_ |= (index & index_bit_mask);
  }

  /**
   * \brief Core value.
   * \details The value's bits contain the following information:
   *	- 1 bit for incoming/outgoing node
   *	- 1 bit for boundary or non-boundary node
   *	- 1 bit for local or non-local node
   *	- 61 bits for index into the appropriate host/device buffer
   */
  std::uint64_t value_ = std::numeric_limits<std::uint64_t>::max();
};

/**
 * \brief Set of device pointers to local, boundary, and non-local buffers for CBCD FLUDS.
 */
struct CBCD_FLUDSPointerSet
{
  /// Pointer to device local angular flux buffer.
  double* local_psi = nullptr;

  /// Pointer to device incoming boundary angular flux buffer.
  double* incoming_boundary_psi = nullptr;

  /// Pointer to device outgoing boundary angular flux buffer.
  double* outgoing_boundary_psi = nullptr;

  /// Pointer to device incoming non-local angular flux buffer.
  double* incoming_nonlocal_psi = nullptr;

  /// Pointer to device outgoing non-local angular flux buffer.
  double* outgoing_nonlocal_psi = nullptr;

  // Stride size (number of groups times number of angles).
  std::uint64_t stride_size = 0;

  /// Default constructor
  CBCD_FLUDSPointerSet() = default;

  /// Copy constructor
  CBCD_FLUDSPointerSet(const CBCD_FLUDSPointerSet&) = default;

  /// Move constructor
  CBCD_FLUDSPointerSet(CBCD_FLUDSPointerSet&&) = default;

  /// Copy and move assignment operator
  CBCD_FLUDSPointerSet& operator=(const CBCD_FLUDSPointerSet&) = default;

  /// Move assignment operator
  CBCD_FLUDSPointerSet& operator=(CBCD_FLUDSPointerSet&&) = default;

  /// Get pointer to the incoming angular flux (if the face is not incoming, a nullptr is returned).
  constexpr double* GetIncomingFluxPointer(const CBCD_NodeIndex& node_index) const noexcept
  {
    // Undefined case (corresponds to a parallel face)
    if (node_index.IsUndefined())
    {
      return nullptr;
    }
    // Outgoing case : nullptr
    if (node_index.IsOutgoing())
    {
      return nullptr;
    }
    // Incoming boundary case
    if (node_index.IsBoundary())
    {
      return incoming_boundary_psi + node_index.GetIndex() * stride_size;
    }
    // Incoming local case
    if (node_index.IsLocal())
    {
      return local_psi + node_index.GetIndex() * stride_size;
    }
    // Incoming non-local case
    else
    {
      return incoming_nonlocal_psi + node_index.GetIndex() * stride_size;
    }
  }

  /// Get pointer to the outgoing angular flux (if the face is not outgoing, a nullptr is returned).
  constexpr double* GetOutgoingFluxPointer(const CBCD_NodeIndex& node_index) const noexcept
  {
    // Undefined case (corresponds to a parallel face)
    if (node_index.IsUndefined())
    {
      return nullptr;
    }
    // Incoming case : nullptr
    if (!node_index.IsOutgoing())
    {
      return nullptr;
    }
    // Outgoing boundary case
    if (node_index.IsBoundary())
    {
      return outgoing_boundary_psi + node_index.GetIndex() * stride_size;
    }
    // Outgoing local case
    if (node_index.IsLocal())
    {
      return local_psi + node_index.GetIndex() * stride_size;
    }
    // Outgoing non-local case
    else
    {
      return outgoing_nonlocal_psi + node_index.GetIndex() * stride_size;
    }
  }
};

} // namespace opensn