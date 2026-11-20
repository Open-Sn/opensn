// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_slot_planner.h"
#include <boost/concept/assert.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/property_map/property_map.hpp>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <span>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace opensn::detail
{
namespace
{

constexpr std::uint32_t INVALID_INDEX = std::numeric_limits<std::uint32_t>::max();

void
ValidateSuccessorOffsets(const std::vector<std::uint32_t>& offsets, std::size_t num_successors)
{
  if (offsets.empty())
  {
    if (num_successors != 0)
      throw std::invalid_argument("CBC slot planner requires offsets for nonempty successors.");
    return;
  }
  if (offsets.front() != 0 or offsets.back() != num_successors or
      not std::is_sorted(offsets.begin(), offsets.end()))
    throw std::invalid_argument(
      "CBC slot planner requires nondecreasing offsets from zero to the successor count.");
}

struct FlowArc
{
  std::uint32_t to = 0;
  std::uint32_t reverse = 0;
  std::uint32_t residual = 0;
};
static_assert(sizeof(FlowArc) == 3 * sizeof(std::uint32_t));

// Boost.Graph view borrowing the CSR offsets and residual arcs.
struct SplitNetwork
{
  using vertex_descriptor = std::uint32_t;
  using edge_descriptor = std::uint32_t;
  using directed_category = boost::directed_tag;
  using edge_parallel_category = boost::allow_parallel_edge_tag;
  struct Traversal : boost::incidence_graph_tag, boost::vertex_list_graph_tag
  {
  };
  using traversal_category = Traversal;
  using out_edge_iterator = boost::counting_iterator<std::uint32_t>;
  using vertex_iterator = boost::counting_iterator<std::uint32_t>;
  // Boost uses these types for n, m and its relabel budget 6*n + m, not just indices.
  using degree_size_type = std::size_t;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  static vertex_descriptor null_vertex() { return INVALID_INDEX; }

  const std::uint32_t* offsets; // num_nodes + 1 entries
  FlowArc* arcs;                // offsets[num_nodes] entries
  std::uint32_t num_nodes;
};

using Edge = SplitNetwork::edge_descriptor;
using OutEdgeIterator = boost::graph_traits<SplitNetwork>::out_edge_iterator;
using VertexIterator = boost::graph_traits<SplitNetwork>::vertex_iterator;

// IncidenceGraph and VertexListGraph (found by argument-dependent lookup).
std::pair<OutEdgeIterator, OutEdgeIterator>
out_edges(std::uint32_t u, const SplitNetwork& g)
{
  return {OutEdgeIterator(g.offsets[u]), OutEdgeIterator(g.offsets[u + 1])};
}
std::size_t
out_degree(std::uint32_t u, const SplitNetwork& g)
{
  return g.offsets[u + 1] - g.offsets[u];
}
std::uint32_t
target(Edge e, const SplitNetwork& g)
{
  return g.arcs[e].to;
}
std::uint32_t
source(Edge e, const SplitNetwork& g)
{
  return g.arcs[g.arcs[e].reverse].to;
}
std::pair<VertexIterator, VertexIterator>
vertices(const SplitNetwork& g)
{
  return {VertexIterator(0), VertexIterator(g.num_nodes)};
}
std::size_t
num_vertices(const SplitNetwork& g)
{
  return g.num_nodes;
}

std::size_t
num_edges(const SplitNetwork& g)
{
  return g.offsets[g.num_nodes];
}

// Capacities depend only on immutable endpoints and input counts, including during
// Boost initialization, flow conversion and assertion checks. Residuals are never read here.
struct CapacityMap
{
  using key_type = Edge;
  using value_type = std::uint32_t;
  using reference = value_type;
  using category = boost::readable_property_map_tag;
  const FlowArc* arcs;
  const std::uint32_t* producer_offsets;
  const std::uint32_t* incoming;
  std::uint32_t source;
  std::uint32_t sink;
  std::uint32_t num_faces;
  value_type operator[](Edge e) const
  {
    const auto from = arcs[arcs[e].reverse].to;
    const auto to = arcs[e].to;
    if (from == source)
      return incoming[to / 2];
    if (to == sink)
      return producer_offsets[from / 2 + 1] - producer_offsets[from / 2];
    return from < to and to < source ? num_faces : 0;
  }
};
CapacityMap::reference
get(const CapacityMap& m, Edge e)
{
  return m[e];
}

struct ResidualCapacityMap
{
  using key_type = Edge;
  using value_type = std::uint32_t;
  using reference = std::uint32_t&;
  using category = boost::lvalue_property_map_tag;
  FlowArc* arcs;
  reference operator[](Edge e) const { return arcs[e].residual; }
};
ResidualCapacityMap::reference
get(const ResidualCapacityMap& m, Edge e)
{
  return m[e];
}
void
put(const ResidualCapacityMap& m, Edge e, std::uint32_t value)
{
  m[e] = value;
}

struct ReverseEdgeMap
{
  using key_type = Edge;
  using value_type = Edge;
  using reference = const Edge&;
  using category = boost::lvalue_property_map_tag;
  const FlowArc* arcs;
  reference operator[](Edge e) const { return arcs[e].reverse; }
};
ReverseEdgeMap::reference
get(const ReverseEdgeMap& m, Edge e)
{
  return m[e];
}

using VertexIndexMap = boost::typed_identity_property_map<std::uint32_t>;

BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<SplitNetwork>));
BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<SplitNetwork>));
BOOST_CONCEPT_ASSERT((boost::ReadablePropertyMapConcept<CapacityMap, Edge>));
BOOST_CONCEPT_ASSERT((boost::Mutable_LvaluePropertyMapConcept<ResidualCapacityMap, Edge>));
BOOST_CONCEPT_ASSERT((boost::LvaluePropertyMapConcept<ReverseEdgeMap, Edge>));
BOOST_CONCEPT_ASSERT((boost::ReadablePropertyMapConcept<VertexIndexMap, std::uint32_t>));

LocalFaceSlotPlanResult
UseIdentityAssignment(std::size_t num_faces, std::vector<std::uint32_t>& face_slot_ids)
{
  face_slot_ids.resize(num_faces);
  std::iota(face_slot_ids.begin(), face_slot_ids.end(), std::uint32_t{0});
  return {num_faces, true};
}

bool
RequiresIdentityAssignment(std::size_t num_tasks, std::size_t num_faces)
{
  if (num_faces > INVALID_INDEX)
    throw std::length_error("CBC local face IDs exceed their 32-bit range.");

  // Before aggregation the network has 2T + 2 vertices and at most 2(3F + T) residual arcs.
  // Bound the largest vertex index by INT_MAX and the CSR arc count by UINT32_MAX.
  return num_tasks > (static_cast<std::size_t>(std::numeric_limits<int>::max()) - 1) / 2 or
         2 * (std::uint64_t{3} * num_faces + num_tasks) > INVALID_INDEX;
}

std::uint64_t
BankNodes(std::span<const std::uint32_t> slots,
          std::size_t count,
          std::span<const std::uint16_t> face_nodes)
{
  std::vector<std::uint16_t> maxima(count, 0);
  for (std::size_t face = 0; face < slots.size(); ++face)
    maxima[slots[face]] = std::max(maxima[slots[face]], face_nodes[face]);
  return std::accumulate(maxima.begin(), maxima.end(), std::uint64_t{0});
}

struct SlotPlannerWorkspace
{
  std::vector<std::uint32_t> offsets;
  std::vector<FlowArc> arcs;
  // Next edge indices for network assembly and flow decomposition.
  std::vector<std::uint32_t> next_edge_indices;
  // Build: incoming face counts, distinct-successor CSR and producer task indices.
  // Solve: marks maps consumers to source arcs. Extract: incoming holds next producer-face indices,
  // and unique_ranks holds consumer IDs, preserving input/output aliases before output mutation.
  std::vector<std::uint32_t> incoming, unique_offsets, unique_ranks, marks;
};

class LocalFaceSlotPlanner
{
public:
  LocalFaceSlotPlanner(const std::vector<std::uint32_t>& successor_rank_offsets,
                       const std::vector<std::uint32_t>& successor_ranks,
                       SlotPlannerWorkspace& workspace,
                       bool in_place = false,
                       std::span<const std::uint16_t> face_nodes = {})
    : successor_rank_offsets_(successor_rank_offsets),
      successor_ranks_(successor_ranks),
      workspace_(workspace),
      in_place_(in_place),
      face_nodes_(face_nodes),
      num_tasks_(static_cast<std::uint32_t>(successor_rank_offsets.size() - 1)),
      num_faces_(static_cast<std::uint32_t>(successor_ranks.size())),
      source_(2 * num_tasks_),
      sink_(source_ + 1),
      num_nodes_(static_cast<std::size_t>(source_) + 2)
  {
    BuildNetwork();
  }

  // Extraction consumes the flow, so solve only once per network.
  // clang-format off
  // R. P. Dilworth, "A decomposition theorem for partially ordered sets", Ann. Math., 1950.
  // D. R. Fulkerson, "Note on Dilworth's decomposition theorem for partially ordered sets", Proc. Amer. Math. Soc., 1956.
  // A. V. Goldberg and R. E. Tarjan, "A new approach to the maximum-flow problem", J. ACM, 1988.
  // B. V. Cherkassky and A. V. Goldberg, "On implementing the push-relabel method for the maximum flow problem", Algorithmica, 1997.
  // clang-format on
  LocalFaceSlotPlanResult Solve(std::vector<std::uint32_t>& face_slot_ids)
  {
    SplitNetwork network{
      workspace_.offsets.data(), workspace_.arcs.data(), static_cast<std::uint32_t>(num_nodes_)};
    const auto matching_size =
      boost::push_relabel_max_flow(network,
                                   source_,
                                   sink_,
                                   CapacityMap{workspace_.arcs.data(),
                                               successor_rank_offsets_.data(),
                                               workspace_.incoming.data(),
                                               source_,
                                               sink_,
                                               num_faces_},
                                   ResidualCapacityMap{workspace_.arcs.data()},
                                   ReverseEdgeMap{workspace_.arcs.data()},
                                   VertexIndexMap{});
    return ExtractPlan(matching_size, face_slot_ids);
  }

private:
  static constexpr std::uint32_t TaskEntry(const std::uint32_t task) noexcept { return 2 * task; }

  static constexpr std::uint32_t TaskExit(const std::uint32_t task) noexcept
  {
    return 2 * task + 1;
  }

  std::uint32_t SourceTarget(std::uint32_t task) const noexcept
  {
    return in_place_ ? TaskEntry(task) : TaskExit(task);
  }

  void AddEdge(const std::uint32_t from, const std::uint32_t to, const std::uint32_t capacity)
  {
    const auto forward = workspace_.next_edge_indices[from]++;
    const auto reverse = workspace_.next_edge_indices[to]++;
    workspace_.arcs[forward] = {to, reverse, capacity};
    workspace_.arcs[reverse] = {from, forward, 0};
  }

  // Split each task into entry and exit. A distinct dependency links the upstream exit
  // to the downstream entry. Internal arcs have capacity F, the number of faces.
  // Source arcs carry consumer face counts and sink arcs carry producer face counts.
  // A source arc to a consumer's entry permits same-cell reuse. Entering at its exit
  // forces the path to traverse at least one dependency before reaching a producer's entry.
  // Each unit of flow links two faces, giving F minus the flow value as the chain count.
  void BuildNetwork()
  {
    auto& workspace = workspace_;
    workspace.incoming.assign(num_tasks_, 0);
    workspace.marks.assign(num_tasks_, INVALID_INDEX);
    workspace.unique_offsets.assign(num_tasks_ + 1, 0);
    workspace.unique_ranks.clear();
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
    {
      for (auto face = successor_rank_offsets_[task]; face < successor_rank_offsets_[task + 1];
           ++face)
      {
        const auto consumer = successor_ranks_[face];
        if (consumer <= task)
          throw std::logic_error("CBC local-face slot planner: task " + std::to_string(task) +
                                 " has non-topological successor " + std::to_string(consumer));
        ++workspace.incoming[consumer];
        if (workspace.marks[consumer] != task)
        {
          workspace.marks[consumer] = task;
          workspace.unique_ranks.push_back(consumer);
        }
      }
      workspace.unique_offsets[task + 1] =
        static_cast<std::uint32_t>(workspace.unique_ranks.size());
    }
    workspace.offsets.assign(num_nodes_ + 1, 0);
    auto Count = [&](auto from, auto to)
    {
      ++workspace.offsets[from + 1];
      ++workspace.offsets[to + 1];
    };
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
    {
      if (workspace.incoming[task])
        Count(source_, SourceTarget(task));
      if (successor_rank_offsets_[task + 1] != successor_rank_offsets_[task])
        Count(TaskEntry(task), sink_);
      Count(TaskEntry(task), TaskExit(task));
      for (auto edge = workspace.unique_offsets[task]; edge < workspace.unique_offsets[task + 1];
           ++edge)
        Count(TaskExit(task), TaskEntry(workspace.unique_ranks[edge]));
    }
    std::partial_sum(workspace.offsets.begin(), workspace.offsets.end(), workspace.offsets.begin());
    workspace.arcs.resize(workspace.offsets.back());
    workspace.next_edge_indices.assign(workspace.offsets.begin(), workspace.offsets.end() - 1);
    workspace.marks.assign(num_tasks_, INVALID_INDEX);
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
      if (workspace.incoming[task])
      {
        workspace.marks[task] = workspace.next_edge_indices[source_];
        AddEdge(source_, SourceTarget(task), workspace.incoming[task]);
      }
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
      if (auto num_produced = successor_rank_offsets_[task + 1] - successor_rank_offsets_[task])
        AddEdge(TaskEntry(task), sink_, num_produced);
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
    {
      AddEdge(TaskEntry(task), TaskExit(task), num_faces_);
      for (auto edge = workspace.unique_offsets[task]; edge < workspace.unique_offsets[task + 1];
           ++edge)
        AddEdge(TaskExit(task), TaskEntry(workspace.unique_ranks[edge]), num_faces_);
    }
  }

  std::uint32_t FlowOnArc(const std::uint32_t arc_index) const noexcept
  {
    return workspace_.arcs[workspace_.arcs[arc_index].reverse].residual;
  }

  void ConsumeFlowUnit(const std::uint32_t arc_index) noexcept
  {
    // The solver has finished. Decrease only the flow read by extraction.
    --workspace_.arcs[workspace_.arcs[arc_index].reverse].residual;
  }

  enum class DecompositionOrder
  {
    ORIGINAL,
    LARGER_FACE_TIE,
    SMALLER_FACE_TIE
  };

  LocalFaceSlotPlanResult ExtractPlan(std::size_t matching_size,
                                      std::vector<std::uint32_t>& face_slot_ids)
  {
    const bool mixed =
      not face_nodes_.empty() and
      std::adjacent_find(face_nodes_.begin(), face_nodes_.end(), std::not_equal_to<>{}) !=
        face_nodes_.end();
    if (not mixed)
    {
      ExtractChains(face_slot_ids, DecompositionOrder::ORIGINAL);
      return ExtractSlotAssignment(matching_size, face_slot_ids);
    }

    std::vector<std::uint32_t> residuals(workspace_.arcs.size());
    for (std::size_t edge = 0; edge < residuals.size(); ++edge)
      residuals[edge] = workspace_.arcs[edge].residual;
    std::vector<std::uint32_t> best;
    ExtractChains(best, DecompositionOrder::ORIGINAL);
    const auto result = ExtractSlotAssignment(matching_size, best);
    auto best_nodes = BankNodes(best, result.num_slots, face_nodes_);
    for (const auto order :
         {DecompositionOrder::LARGER_FACE_TIE, DecompositionOrder::SMALLER_FACE_TIE})
    {
      for (std::size_t edge = 0; edge < residuals.size(); ++edge)
        workspace_.arcs[edge].residual = residuals[edge];
      std::vector<std::uint32_t> candidate;
      ExtractChains(candidate, order);
      ExtractSlotAssignment(matching_size, candidate);
      const auto nodes = BankNodes(candidate, result.num_slots, face_nodes_);
      if (nodes < best_nodes)
      {
        best_nodes = nodes;
        best = std::move(candidate);
      }
    }
    // Output may alias either input. All candidates finish reading them before this assignment.
    face_slot_ids = std::move(best);
    return result;
  }

  // Decompose flow into paths along forward arcs. Flow conservation and task order
  // guarantee that each path reaches the sink. Exhausted arcs can be skipped permanently.
  void ExtractChains(std::vector<std::uint32_t>& previous_face, DecompositionOrder order)
  {
    // Preserve both CSR inputs before writing output, which may alias either input.
    // The distinct-successor list and incoming counts have finished their build-phase use.
    workspace_.unique_ranks.assign(successor_ranks_.begin(), successor_ranks_.end());
    workspace_.incoming.assign(successor_rank_offsets_.begin(), successor_rank_offsets_.end() - 1);
    previous_face.assign(num_faces_, INVALID_INDEX);
    workspace_.next_edge_indices.assign(workspace_.offsets.begin(),
                                        workspace_.offsets.begin() + source_);
    std::vector<std::uint32_t> pool;
    std::vector<std::uint16_t> chain_max;
    if (order != DecompositionOrder::ORIGINAL)
    {
      pool.resize(num_faces_);
      std::iota(pool.begin(), pool.end(), std::uint32_t{0});
      chain_max.assign(face_nodes_.begin(), face_nodes_.end());
    }

    for (std::uint32_t u_face = 0; u_face < num_faces_; ++u_face)
    {
      // Faces with the same consumer use successive units of its aggregated source arc.
      const auto source_arc = workspace_.marks[workspace_.unique_ranks[u_face]];
      if (FlowOnArc(source_arc) == 0)
        continue;

      ConsumeFlowUnit(source_arc);
      auto node = workspace_.arcs[source_arc].to;
      while (true)
      {
        auto& next_edge = workspace_.next_edge_indices[node];
        const auto arc_end = workspace_.offsets[node + 1];
        while (next_edge < arc_end and
               (workspace_.arcs[next_edge].to <= node or workspace_.arcs[next_edge].to == source_ or
                FlowOnArc(next_edge) == 0))
          ++next_edge;
        if (next_edge == arc_end)
          throw std::logic_error("CBC local-face slot planner: flow decomposition found no arc.");

        const auto arc_index = next_edge;
        ConsumeFlowUnit(arc_index);
        const auto to = workspace_.arcs[arc_index].to;
        if (to == sink_)
        {
          const auto begin = workspace_.incoming[node / 2]++;
          auto v_face = begin;
          if (order != DecompositionOrder::ORIGINAL)
          {
            auto chosen = begin;
            const auto cap = chain_max[u_face];
            const auto end = successor_rank_offsets_[node / 2 + 1];
            for (auto position = begin + 1; position < end; ++position)
            {
              const auto weight = face_nodes_[pool[position]];
              const auto chosen_weight = face_nodes_[pool[chosen]];
              // Merging a singleton face into this chain saves min(chain maximum, face nodes).
              const auto saving = std::min(cap, weight);
              const auto chosen_saving = std::min(cap, chosen_weight);
              const bool better_tie = order == DecompositionOrder::LARGER_FACE_TIE
                                        ? weight > chosen_weight
                                        : weight < chosen_weight;
              if (saving > chosen_saving or
                  (saving == chosen_saving and
                   (better_tie or (weight == chosen_weight and pool[position] < pool[chosen]))))
                chosen = position;
            }
            std::swap(pool[begin], pool[chosen]);
            v_face = pool[begin];
            chain_max[v_face] = std::max(chain_max[u_face], face_nodes_[v_face]);
          }
          previous_face[v_face] = u_face;
          break;
        }
        node = to;
      }
    }
  }

  LocalFaceSlotPlanResult ExtractSlotAssignment(const std::size_t matching_size,
                                                std::vector<std::uint32_t>& face_slot_ids) const
  {
    // A predecessor's producer precedes this face's producer, so its slot ID is available.
    std::uint32_t next_slot = 0;
    for (auto& slot : face_slot_ids)
      slot = slot == INVALID_INDEX ? next_slot++ : face_slot_ids[slot];

    const auto slot_count = static_cast<std::size_t>(num_faces_) - matching_size;
    if (next_slot != slot_count)
      throw std::logic_error("CBC local-face slot planner: the chains do not match the flow.");
    return {slot_count, false};
  }

  const std::vector<std::uint32_t>& successor_rank_offsets_;
  const std::vector<std::uint32_t>& successor_ranks_;
  SlotPlannerWorkspace& workspace_;
  bool in_place_;
  std::span<const std::uint16_t> face_nodes_;
  std::uint32_t num_tasks_ = 0;
  std::uint32_t num_faces_ = 0;
  std::uint32_t source_ = 0;
  std::uint32_t sink_ = 0;
  std::size_t num_nodes_ = 0;
};

} // namespace

LocalFaceSlotPlanResult
BuildMinimumLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                              const std::vector<std::uint32_t>& successor_ranks,
                              std::vector<std::uint32_t>& face_slot_ids)
{
  ValidateSuccessorOffsets(successor_rank_offsets, successor_ranks.size());
  if (successor_ranks.empty())
  {
    face_slot_ids.clear();
    return {};
  }

  const auto num_tasks = successor_rank_offsets.size() - 1;
  const auto num_faces = successor_ranks.size();
  if (RequiresIdentityAssignment(num_tasks, num_faces))
    return UseIdentityAssignment(num_faces, face_slot_ids);

  SlotPlannerWorkspace workspace;
  LocalFaceSlotPlanner planner(successor_rank_offsets, successor_ranks, workspace);
  return planner.Solve(face_slot_ids);
}

LocalFaceSlotPlanResult
BuildMinimumLocalFaceSlotPlan(const std::vector<std::uint32_t>& successor_rank_offsets,
                              const std::vector<std::uint32_t>& successor_ranks,
                              std::vector<std::uint32_t>& face_slot_ids,
                              std::span<const std::uint16_t> face_nodes)
{
  ValidateSuccessorOffsets(successor_rank_offsets, successor_ranks.size());
  if (face_nodes.size() != successor_ranks.size() or
      std::ranges::any_of(face_nodes, [](auto count) { return count == 0; }))
    throw std::invalid_argument("CBC slot planner requires one positive node count per face.");
  if (successor_ranks.empty())
  {
    face_slot_ids.clear();
    return {};
  }
  const auto num_tasks = successor_rank_offsets.size() - 1;
  const auto num_faces = successor_ranks.size();
  if (RequiresIdentityAssignment(num_tasks, num_faces))
    return UseIdentityAssignment(num_faces, face_slot_ids);

  const bool mixed =
    std::adjacent_find(face_nodes.begin(), face_nodes.end(), std::not_equal_to<>{}) !=
    face_nodes.end();
  std::vector<std::uint32_t> strict;
  LocalFaceSlotPlanResult strict_result;
  SlotPlannerWorkspace workspace;
  if (mixed)
  {
    LocalFaceSlotPlanner strict_planner(successor_rank_offsets, successor_ranks, workspace);
    strict_result = strict_planner.Solve(strict);
  }

  std::vector<std::uint32_t> candidate;
  LocalFaceSlotPlanResult result;
  {
    LocalFaceSlotPlanner planner(
      successor_rank_offsets, successor_ranks, workspace, true, face_nodes);
    result = planner.Solve(candidate);
  }
  if (mixed and BankNodes(strict, strict_result.num_slots, face_nodes) <
                  BankNodes(candidate, result.num_slots, face_nodes))
  {
    candidate = std::move(strict);
    result = strict_result;
    result.used_strict_fallback = true;
  }
  // Both input arrays may alias the output; publish only after every input read is complete.
  face_slot_ids = std::move(candidate);
  return result;
}

} // namespace opensn::detail
