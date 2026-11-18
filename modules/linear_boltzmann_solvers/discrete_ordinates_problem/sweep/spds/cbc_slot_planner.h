// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace opensn::detail
{

/// Local-face slot assignment size and fallback diagnostics.
struct LocalFaceSlotPlanResult
{
  /// Number of slots in the returned plan.
  std::size_t num_slots = 0;
  /// Whether conservative 32-bit network bounds require returning one slot per face.
  bool used_identity_fallback = false;
  /// Whether a mixed-face plan uses strict precedence to bound its bank size.
  bool used_strict_fallback = false;
};

/**
 * Assign storage slots to ordinary local angular-flux faces in a host CBC sweep.
 *
 * Each task is a cell on this rank, indexed by topological rank after removing delayed
 * dependencies. Each CSR successor entry represents a directed face from its row's producer
 * to the named consumer. Parallel faces have separate entries. Delayed faces, nonlocal
 * faces and reflecting boundaries use separate storage. Curvilinear RZ and spherical
 * sweeps use AAH.
 *
 * A slot can pass from face A to face B when a nonempty dependency path leads from A's
 * consumer to B's producer. The scheduler finishes each task before releasing its
 * successors. The plan uses the minimum number of slots under this order
 * unless the network exceeds the 32-bit index bounds. The assignment is deterministic for
 * fixed input and Boost version. Workspace is local to the call and uses O(T + F) space
 * for T tasks and F faces, including the flow network and decomposition, excluding input.
 *
 * \param successor_rank_offsets Nondecreasing CSR offsets, starting at zero and ending
 *        at the face count. Empty offsets are allowed for an empty successor array.
 * \param successor_ranks Consumer indices. Each must be in range and greater than its producer.
 * \param face_slot_ids Output slot IDs in face order. May alias either CSR input.
 * \return Slot count and whether the index bounds required one slot per face.
 * \throw std::invalid_argument If the CSR offsets are inconsistent.
 * \throw std::length_error If the face count exceeds the 32-bit range.
 * \throw std::logic_error If a consumer does not follow its producer or decomposition fails.
 */
LocalFaceSlotPlanResult
BuildMinimumLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                              const std::vector<std::uint32_t>& successor_ranks,
                              std::vector<std::uint32_t>& face_slot_ids);

/**
 * Assign ordinary local face slots with read-before-write reuse.
 *
 * This overload also permits reuse when A's consumer is B's producer. Directions
 * and groups occupy distinct entries.
 *
 * \note For a given cell, direction and group, CBC kernels must read all incoming
 * local values before the first outgoing local write. An interleaved kernel must
 * stage its incoming values before writing or use strict-precedence slots.
 *
 * A slot reserves the largest face-node count assigned to it. For mixed face sizes the
 * planner compares three decompositions of the minimum-slot flow by total slot nodes.
 * If all exceed the strict plan's bank size, it returns that plan and sets
 * used_strict_fallback. The minimum slot count applies under the returned order unless
 * the index bounds require one slot per face. The bank is bounded by the strict plan's
 * size for the same input and Boost version. The size comparison is a heuristic.
 *
 * The CSR, output aliasing and workspace contracts match the strict overload.
 * For nonempty input within the index bounds, the planner uses one flow solve for
 * uniform face sizes and two sequential solves sharing one workspace for mixed sizes.
 *
 * \param successor_rank_offsets CSR offsets indexed by producer task.
 * \param successor_ranks Consumer indices, one per directed face.
 * \param face_slot_ids Output slot IDs in face order.
 * \param face_nodes Positive node counts, one per directed face.
 * \return Slot count and diagnostics for the strict and index-bound fallbacks.
 * \throw std::invalid_argument If CSR offsets or face-node counts are inconsistent or a count is
 * zero.
 * \throw std::length_error If the face count exceeds the 32-bit range.
 * \throw std::logic_error If a consumer does not follow its producer or decomposition fails.
 */
LocalFaceSlotPlanResult
BuildMinimumLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                              const std::vector<std::uint32_t>& successor_ranks,
                              std::vector<std::uint32_t>& face_slot_ids,
                              std::span<const std::uint16_t> face_nodes);

} // namespace opensn::detail
