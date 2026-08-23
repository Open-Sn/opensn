// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace opensn
{
class Cell;

namespace detail
{

/// Checked addition in the existing 32-bit CBC task/face index domain. Setup only.
inline std::uint32_t
CheckedCBCTaskCount(std::uint32_t count, std::size_t increment)
{
  if (increment > std::numeric_limits<std::uint32_t>::max() - count)
    throw std::length_error("CBC host task count exceeds its 32-bit storage range.");
  return count + static_cast<std::uint32_t>(increment);
}

/**
 * Host-only immutable CSR task graph, indexed by mesh-local cell ID.
 *
 * Each outgoing non-delayed local FACE contributes one successor, in cell-face order;
 * parallel faces are not merged. Offsets have C+1 entries, including {0} on an empty
 * rank. Cell pointers borrow the SPDS-owned mesh. All arrays are replaced on a task
 * rebuild; angle sets must reset before borrowing them again. No sweep mutates them.
 */
struct CBCHostTaskGraph
{
  std::vector<const Cell*> cells;
  std::vector<std::uint32_t> successor_offsets;
  std::vector<std::uint32_t> successors;
};

} // namespace detail
} // namespace opensn
