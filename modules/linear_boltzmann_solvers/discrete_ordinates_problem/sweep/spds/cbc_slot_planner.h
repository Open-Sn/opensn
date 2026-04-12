// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace opensn::detail
{

/// Result of an exact local-face slot-planning solve.
struct SlotSolveResult
{
  /// Exact number of reusable slots required by the computed chain cover.
  std::size_t slot_count = 0;
  /// Flag indicating that the post-solve verifier rejected the computed assignment.
  bool verifier_rejected = false;
};

/**
 * Compute the exact minimum safe local-face slot assignment.
 *
 * Let `F` denote the local directed faces, and define `u < v` when the consumer cell of
 * face `u` reaches the producer cell of face `v` in the local CBC task DAG. This is the
 * safe reuse relation: if `u < v`, then every admissible CBC or CBCD sweep consumes the
 * angular flux stored for `u` before `v` may overwrite the same slot.
 *
 * Computing the minimum number of reusable cell-face slots is equivalent to
 * the minimum chain-cover problem for the induced face poset.
 * A chain is one statically reusable slot. The minimum number of slots equals the poset
 * width (i.e. the maximum cardinality of any antichain of pairwise incomparable
 * faces). By Dilworth's theorem, this is the minimum chain-cover cardinality.
 *
 * The implementation uses the standard bipartite split-graph reduction. The reuse relation
 * defines the bipartite edges, Hopcroft-Karp computes a maximum cardinality matching, and
 * the matching induces a minimum chain cover of size `|F| - |M|`. Koenig's theorem provides
 * the matching-cover duality for the bipartite graph. Consequently, the returned slot count
 * is exact.
 *
 * Algorithm flow:
 * 1. Build the reflexive transitive closure of the local CBC task DAG in topological-rank
 * space.
 * 2. Group local directed faces by consumer-cell rank and cache the reachable producer-cell
 * ranks that define the reuse graph rows.
 * 3. Run Hopcroft-Karp on the implicit bipartite reuse graph:
 *    a. greedy seeding,
 *    b. BFS layer construction,
 *    c. iterative DFS augmentation.
 * 4. Extract one slot chain per unmatched right-side face.
 * 5. Verify the extracted assignment and report whether the caller should fall back to the
 * identity assignment.
 *
 * After chain extraction, the assignment is verified by checking each consecutive reuse
 * handoff in face enumeration order. If the verifier rejects the result, the caller may
 * conservatively fall back to the identity assignment.
 *
 * \param successor_rank_offsets Offsets into the flat successor-rank adjacency list of the
 * local CBC task DAG.
 * \param successor_ranks Flat successor-rank adjacency list of the local CBC task DAG.
 * \param face_producer_ranks Producer-cell topological rank for each local directed face.
 * \param face_consumer_ranks Consumer-cell topological rank for each local directed face.
 * \param producer_cell_face_offsets Offsets grouping local faces by producer-cell topological
 * rank.
 * \param face_slot_ids Output slot assignment keyed by local face rank.
 * \return Exact slot count and verifier status for the computed assignment.
 */
SlotSolveResult
ComputeLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                         const std::vector<std::uint32_t>& successor_ranks,
                         const std::vector<std::uint32_t>& face_producer_ranks,
                         const std::vector<std::uint32_t>& face_consumer_ranks,
                         const std::vector<std::uint32_t>& producer_cell_face_offsets,
                         std::vector<std::uint32_t>& face_slot_ids);

} // namespace opensn::detail
