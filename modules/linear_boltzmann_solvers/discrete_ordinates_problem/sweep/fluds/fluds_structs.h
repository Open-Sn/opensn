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
 * Base class for 64-bit integer encoding the index and bank information of a face node.
 * This base class manages the common bits:
 * - Bit 63: Incoming/Outgoing
 * - Bit 62: Boundary
 *
 * Derived classes (AAHD_NodeIndex, CBCD_NodeIndex) manage the remaining bits (Local, Delayed,
 * Index) as their layouts differ.
 */
class NodeIndex
{
public:
  /// Default constructor.
  constexpr NodeIndex() = default;

  /// Direct assign core value.
  constexpr NodeIndex(const std::uint64_t& value) : value_(value) {}

  /// Check if the current index is undefined.
  constexpr bool IsUndefined() const noexcept
  {
    return value_ == std::numeric_limits<std::uint64_t>::max();
  }

  /// Check if the current index is incoming or outgoing.
  constexpr bool IsOutgoing() const noexcept { return (value_ & inout_bit_mask) != 0; }

  /// Check if the current index corresponds to boundary.
  constexpr bool IsBoundary() const noexcept { return (value_ & boundary_bit_mask) != 0; }

  /// Get the core value.
  constexpr std::uint64_t GetCoreValue() const noexcept { return value_; }

protected:
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

  /// Core value.
  std::uint64_t value_ = std::numeric_limits<std::uint64_t>::max();
};

/**
 * Base pointer set for FLUDs on the device.
 * Holds the common pointers and stride information shared between
 * the device AAH and CBC FLUDS implementations.
 */
struct FLUDSPointerSet
{
  /// Pointer to local angular fluxes.
  double* __restrict__ local_psi = nullptr;
  /// Pointer to non-local incoming non-local angular fluxes.
  double* __restrict__ nonlocal_incoming_psi = nullptr;
  /// Pointer to non-local outgoing non-local angular fluxes.
  double* __restrict__ nonlocal_outgoing_psi = nullptr;
  /// Stride size (number of groups times number of angles).
  std::uint64_t stride_size = 0;
};

} // namespace opensn