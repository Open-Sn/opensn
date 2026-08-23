// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace opensn::detail
{

/// Result of constructing an exact local-face slot plan.
struct LocalFaceSlotPlanResult
{
  /// Number of slots in the returned plan.
  std::size_t num_slots = 0;
  /// Whether an internal consistency check required the identity plan.
  bool used_identity_fallback = false;
};

/**
 * Compute an exact minimum safe local-face slot assignment.
 *
 * For directed local faces A and B, define A < B when the consumer task of A strictly
 * precedes the producer task of B in the reduced task DAG. Faces in a chain of this
 * strict partial order have nonoverlapping lifetimes and may share one slot. By
 * Dilworth's theorem, the minimum number of chains in a partition equals the poset
 * width, which is the maximum number of pairwise incomparable faces.
 *
 * The planner obtains this minimum without materializing the face-poset transitive
 * closure. Its sparse flow network contains split task vertices. Unit-capacity source
 * and sink arcs identify predecessor and successor faces, whereas internal task-DAG arcs
 * have capacity at least the maximum possible flow. Because every capacity is integral,
 * the maximum flow is integral and is equivalent to a maximum matching in the implicit
 * face-reuse bipartite graph. Thus, if F is the number of directed local faces and M is
 * the maximum-flow value, the minimum slot count is F - M.
 *
 * An internal flow-extraction or chain-consistency failure returns the conservative identity
 * assignment, in which every face has its own slot.
 *
 * \param successor_rank_offsets CSR offsets for the reduced task DAG.
 * \param successor_ranks Task-DAG successors in topological-rank coordinates.
 * \param face_producer_ranks Producer task rank for each directed local face.
 * \param face_consumer_ranks Consumer task rank for each directed local face.
 * \param producer_face_offsets CSR offsets for faces grouped by producer task.
 * \param face_slot_ids Output slot assignment indexed by directed local-face ID.
 * \return The number of slots and whether the identity fallback was used.
 */
LocalFaceSlotPlanResult
BuildMinimumLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                              const std::vector<std::uint32_t>& successor_ranks,
                              const std::vector<std::uint32_t>& face_producer_ranks,
                              const std::vector<std::uint32_t>& face_consumer_ranks,
                              const std::vector<std::uint32_t>& producer_face_offsets,
                              std::vector<std::uint32_t>& face_slot_ids);

} // namespace opensn::detail
