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
  /// Whether conservative 32-bit network bounds require returning one slot per face.
  bool used_identity_fallback = false;
};

/**
 * Compute a safe local-face slot assignment with the minimum number of slots.
 *
 * The tasks are the local cells, numbered by topological rank in the local task DAG after
 * removing delayed edges. Every arc p -> c of this DAG is one directed local face, produced by
 * task p and consumed by task c; parallel arcs are distinct faces. Face k is entry k of the DAG's
 * CSR arrays: its producer is the row of k and its consumer is `successor_ranks[k]`.
 *
 * For faces A and B, define A < B when the consumer of A strictly precedes the producer of B in
 * the DAG. Every schedule that respects the DAG then finishes reading A before B is written, so
 * the faces of a chain of this strict partial order may share one slot. By Dilworth's theorem,
 * the minimum number of chains that cover the faces equals the width of the order, F - M, where F
 * is the number of faces and M is a maximum matching in the bipartite graph with an edge (A, B)
 * for every A < B (Fulkerson's reduction).
 *
 * The planner computes M as a maximum flow without building the transitive closure, by
 * boost::push_relabel_max_flow (highest-label push-relabel with gap and global-relabel
 * heuristics) on a Boost.Graph view of the network. Every task t is split into entry(t) -> exit(t),
 * and every DAG arc u -> v becomes exit(u) -> entry(v), all with capacity F. Face A adds a unit arc
 * from the source to exit(consumer of A), and face B a unit arc from entry(producer of B) to the
 * sink. An augmenting path enters at an exit node and leaves at an entry node, so it uses at least
 * one DAG arc: the split encodes strict precedence. Internal arcs may carry several paths, so the
 * integral maximum flow equals M. Equivalent source/sink arcs are aggregated by consumer/producer
 * task with capacity equal to the number of corresponding faces; duplicate DAG arcs collapse to
 * one capacity-F arc. Summing parallel flows preserves feasibility, and an integral aggregated
 * flow can be expanded by assigning distinct faces within each endpoint group. Thus these
 * reductions preserve M. The chains follow from a decomposition of the flow into paths, assigning
 * endpoint faces in their original order.
 *
 * The objective is the number of slots, not the number of nodes: a slot holds the largest of its
 * faces, so for faces of unequal sizes the slot bank is not the smallest possible. Among the
 * plans with the minimum number of slots, the plan is a deterministic function of the input for
 * a given Boost version. Another Boost version or maximum-flow algorithm may return a different
 * minimum plan and change the slot-bank bytes.
 *
 * References: A. V. Goldberg and R. E. Tarjan, "A new approach to the maximum-flow problem",
 * J. ACM 35(4), 921-940, 1988; B. V. Cherkassky and A. V. Goldberg, "On implementing the
 * push-relabel method for the maximum flow problem", Algorithmica 19(4), 390-410, 1997.
 * Boost implements the highest-label method with gap and global relabeling described there.
 *
 * \param successor_rank_offsets CSR offsets of the task successors, indexed by topological rank
 *        (one more entry than tasks).
 * \param successor_ranks Topological ranks of the task successors; one entry per face.
 * \param face_slot_ids Output slot ID of each face, indexed as `successor_ranks`.
 * \return The number of slots and whether one slot per face was returned, because the flow network
 *         exceeds the conservative 32-bit bounds for the unaggregated network.
 * \throw std::length_error If there are more faces than 32-bit face IDs.
 * \throw std::logic_error If a successor rank is not larger than its task's rank, or the flow
 *         decomposition is inconsistent (an internal error).
 */
LocalFaceSlotPlanResult
BuildMinimumLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                              const std::vector<std::uint32_t>& successor_ranks,
                              std::vector<std::uint32_t>& face_slot_ids);

} // namespace opensn::detail
