// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <ostream>
#include <stdexcept>

namespace opensn
{

/**
 * \brief Node on a face, which is also an edge in the sweep graph.
 * \details This class represents a node on a face. It stores the cell local index, the face index
 * within the cell and the index of the node relative to the face as compact 64-bit integer. This
 * class abstracts the node for usage with std::map and std::set.
 */
class DeviceFaceNode
{
public:
  /// Default constructor
  constexpr DeviceFaceNode() = default;

  /// Member constructor
  constexpr DeviceFaceNode(std::uint32_t cell_idx,
                           std::uint16_t face_idx,
                           std::uint16_t face_node_idx)
    : value_(0)
  {
    // Encode the indices into the 64-bit data member
    value_ |= static_cast<std::uint64_t>(cell_idx) << 32;
    value_ |= static_cast<std::uint64_t>(face_idx) << 16;
    value_ |= static_cast<std::uint64_t>(face_node_idx);
  }

  /// Comparison operator for ordering.
  constexpr bool operator<(const DeviceFaceNode& other) const { return value_ < other.value_; }

  /// Equality operator.
  constexpr bool operator==(const DeviceFaceNode& other) const { return value_ == other.value_; }

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
 * \brief 64-bit integer encoding the index and the bank information of a face node.
 * \note The index must follow strict mutual-exclusion rules:
 *   - Always check boundary status BEFORE local status. Boundary instances are assigned
 *     as local, but their locations are actually set to an undefined value.
 *   - Boundary and delayed status are mutually exclusive. Both cannot be true at the same time.
 *   - Boundary and non-local status are mutually exclusive. Both cannot be true at the same time.
 *   - Non-local outgoing and delayed status are mutually exclusive. The sweep scheduler does not
 *     distinguish between delayed and non-delayed non-local outgoing faces. For non-local outgoing
 *     faces, delayed status is always set to false.
 */
class DeviceNodeIndex
{
public:
  /// Default constructor.
  constexpr DeviceNodeIndex() = default;

  /// Direct assign core value.
  constexpr DeviceNodeIndex(const std::uint64_t& value) : value_(value) {}

  /**
   * Contruct a non-boundary node index.
   * \param index Index into the corresponding bank. Cannot exceed 2^60 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   * \param is_local Flag indicating if the index is in a local bank.
   * \param is_delayed Flag indicating if the index is in a delayed bank (always false for CBC).
   * \param loc_index Index into delayed/non-delayed location dependencies/successors vector for
   * non-local nodes. Cannot exceed 2^21 - 1.
   */
  DeviceNodeIndex(std::uint64_t index, bool is_outgoing, bool is_local, bool is_delayed)
  {
    if (index >= (std::uint64_t(1) << 60) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^60.");
    SetInOut(is_outgoing);
    SetLocal(is_local);
    SetDelayed(is_delayed);
    SetBoundary(false);
    SetIndex(index);
    if (is_delayed && is_outgoing && !is_local)
      throw std::runtime_error(
        "Non-local outgoing nodes cannot be in delayed banks for AAHD FLUDS.");
  }

  /**
   * Contruct a boundary node index.
   * \param index Index into the corresponding bank. Cannot exceed 2^40 - 1.
   * \param is_outgoing Flag indicating if the node corresponds to an outgoing face.
   */
  DeviceNodeIndex(std::uint64_t index, bool is_outgoing)
  {
    if (index >= (std::uint64_t(1) << 60) - 1)
      throw std::runtime_error("Cannot hold an index greater than 2^60.");
    SetInOut(is_outgoing);
    SetDelayed(false);
    SetLocal(true);
    SetBoundary(true);
    SetIndex(index);
  }

  /// Check if the current index is undefined.
  constexpr bool IsUndefined() const noexcept
  {
    return value_ == std::numeric_limits<std::uint64_t>::max();
  }

  /// Check if the current index is incoming or outgoing.
  constexpr bool IsOutgoing() const noexcept { return (value_ & inout_bit_mask) != 0; }

  /// Check if the current index corresponds to a delayed bank.
  constexpr bool IsDelayed() const noexcept { return (value_ & delayed_bit_mask) != 0; }

  /// Check if the current index corresponds to boundary.
  constexpr bool IsBoundary() const noexcept { return (value_ & boundary_bit_mask) != 0; }

  /// Check if the current index corresponds to a local bank.
  constexpr bool IsLocal() const noexcept { return (value_ & local_bit_mask) != 0; }

  /// Get the index into the bank.
  constexpr std::uint64_t GetIndex() const noexcept { return value_ & index_bit_mask; }

  /// Get the core value.
  constexpr std::uint64_t GetCoreValue() const noexcept { return value_; }

private:
  /// \name Incoming/outgoing bit
  /// \{
  /// First bit mask (``1`` followed by 63 zeros).
  static constexpr std::uint64_t inout_bit_mask = std::uint64_t(1) << (64 - 1);
  /// Encode the value as incoming or outgoing.
  constexpr void SetInOut(bool is_outgoing) noexcept
  {
    if (is_outgoing)
      value_ |= inout_bit_mask;
    else
      value_ &= ~inout_bit_mask;
  }
  /// \}

  /// \name Boundary bit
  /// \{
  /// Second bit mask (``01`` followed by 62 zeros).
  static constexpr std::uint64_t boundary_bit_mask = std::uint64_t(1) << (64 - 2);
  /// Encode the value as boundary.
  constexpr void SetBoundary(bool is_boundary) noexcept
  {
    if (is_boundary)
      value_ |= boundary_bit_mask;
    else
      value_ &= ~boundary_bit_mask;
  }
  /// \}

  /// \name Delayed bit
  /// \{
  /// Third bit mask (``001`` followed by 61 zeros).
  static constexpr std::uint64_t delayed_bit_mask = std::uint64_t(1) << (64 - 3);
  /// Encode the value as delayed.
  constexpr void SetDelayed(bool is_delayed) noexcept
  {
    if (is_delayed)
      value_ |= delayed_bit_mask;
    else
      value_ &= ~delayed_bit_mask;
  }
  /// \}

  /// \name Local bit
  /// \{
  /// Third bit mask (``0001`` followed by 60 zeros).
  static constexpr std::uint64_t local_bit_mask = std::uint64_t(1) << (64 - 4);
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
  /// Index bit mask (``1`` at the last 60 bits).
  static constexpr std::uint64_t index_bit_mask = (std::uint64_t(1) << (64 - 4)) - 1;
  /// Encode the index.
  constexpr void SetIndex(std::uint64_t index) noexcept
  {
    value_ &= ~index_bit_mask;
    value_ |= (index & index_bit_mask);
  }
  /// \}

  /**
   * \brief Core value.
   * \details 1 bit for incoming/outgoing; 1 bit for delayed or non delayed; 1 bit for boundary or
   * not; 21 bit (~2.1M) for MPI rank (all ones if local); 40 bits (~4T) for index.
   */
  std::uint64_t value_ = std::numeric_limits<std::uint64_t>::max();
};

} // namespace opensn

namespace std
{

/// Print face node.
ostream& operator<<(ostream& out, const opensn::DeviceFaceNode& n);

/// Print node index.
ostream& operator<<(ostream& out, const opensn::DeviceNodeIndex& e);

} // namespace std