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

/**
 * Add a setup count within the 32-bit CBC task and face index range.
 *
 * \param count Current count.
 * \param increment Number of entries to add.
 * \return Sum representable as a 32-bit count.
 * \throw std::length_error If the sum exceeds the 32-bit range.
 */
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
 * Each outgoing non-delayed local face contributes one successor, in cell-face order;
 * parallel faces are not merged. Offsets have C+1 entries, including {0} on an empty
 * rank. Cell pointers borrow the SPDS-owned mesh. All arrays are replaced on a task
 * rebuild; angle sets must reset before borrowing them again. No sweep mutates them.
 */
struct CBCHostTaskGraph
{
  /// Borrowed cell pointers indexed by local cell ID.
  std::vector<const Cell*> cells;
  /// Prefix offsets into successors, including a final end offset.
  std::vector<std::uint32_t> successor_offsets;
  /// Local successor cell IDs, with one entry per ordinary outgoing local face.
  std::vector<std::uint32_t> successors;
};

} // namespace detail
} // namespace opensn
