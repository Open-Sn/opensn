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
#include <limits>
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

struct FlowArc
{
  std::uint32_t to = 0;
  std::uint32_t reverse = 0;
  std::uint32_t capacity = 0;
  std::uint32_t residual = 0;
};

/**
 * Non-owning Boost.Graph view of the split-task network in CSR form. Vertices and edges are
 * 32-bit indices. The view and maps remain valid while their arrays are neither resized nor
 * freed. Internal linkage lets the compiler inline the view operations even in shared builds.
 */
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

/// Stored initial capacities remain valid through Boost's flow/optimality assertion checks.
struct CapacityMap
{
  using key_type = Edge;
  using value_type = std::uint32_t;
  using reference = const std::uint32_t&;
  using category = boost::lvalue_property_map_tag;
  const FlowArc* arcs;
  reference operator[](Edge e) const { return arcs[e].capacity; }
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

// Boost's CSR vertex_index map is typed_identity_property_map, and push-relabel reads it only
// through get(index, v). Follow that convention despite the documentation's Lvalue wording.
using VertexIndexMap = boost::typed_identity_property_map<std::uint32_t>;

BOOST_CONCEPT_ASSERT((boost::IncidenceGraphConcept<SplitNetwork>));
BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<SplitNetwork>));
BOOST_CONCEPT_ASSERT((boost::LvaluePropertyMapConcept<CapacityMap, Edge>));
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

struct SlotPlannerWorkspace
{
  std::vector<std::uint32_t> offsets;
  std::vector<FlowArc> arcs;
  // Network write offsets, then decomposition cursors. Boost owns its per-call workspace.
  std::vector<std::uint32_t> cursors;
  // Build: incoming multiplicities, distinct-successor CSR and deduplication stamps.
  // Solve: marks maps consumers to source arcs. Extract: incoming holds producer-face cursors,
  // and unique_ranks holds consumer IDs, preserving input/output aliases before output mutation.
  std::vector<std::uint32_t> incoming, unique_offsets, unique_ranks, marks;
};

class LocalFaceSlotPlanner
{
public:
  LocalFaceSlotPlanner(const std::vector<std::uint32_t>& successor_rank_offsets,
                       const std::vector<std::uint32_t>& successor_ranks,
                       SlotPlannerWorkspace& workspace)
    : successor_rank_offsets_(successor_rank_offsets),
      successor_ranks_(successor_ranks),
      workspace_(workspace),
      num_tasks_(static_cast<std::uint32_t>(successor_rank_offsets.size() - 1)),
      num_faces_(static_cast<std::uint32_t>(successor_ranks.size())),
      source_(2 * num_tasks_),
      sink_(source_ + 1),
      num_nodes_(static_cast<std::size_t>(source_) + 2)
  {
    BuildNetwork();
  }

  /// Solve once per constructed planner: extraction consumes the network's flow. Every public
  /// call constructs a new planner and rebuilds all arcs in the reusable workspace.
  LocalFaceSlotPlanResult Solve(std::vector<std::uint32_t>& face_slot_ids)
  {
    SplitNetwork network{
      workspace_.offsets.data(), workspace_.arcs.data(), static_cast<std::uint32_t>(num_nodes_)};
    const auto matching_size =
      boost::push_relabel_max_flow(network,
                                   source_,
                                   sink_,
                                   CapacityMap{workspace_.arcs.data()},
                                   ResidualCapacityMap{workspace_.arcs.data()},
                                   ReverseEdgeMap{workspace_.arcs.data()},
                                   VertexIndexMap{});
    ExtractChains(face_slot_ids);
    return ExtractSlotAssignment(matching_size, face_slot_ids);
  }

private:
  static constexpr std::uint32_t TaskEntry(const std::uint32_t task) noexcept { return 2 * task; }

  static constexpr std::uint32_t TaskExit(const std::uint32_t task) noexcept
  {
    return 2 * task + 1;
  }

  void AddEdge(const std::uint32_t from, const std::uint32_t to, const std::uint32_t capacity)
  {
    const auto forward = workspace_.cursors[from]++;
    const auto reverse = workspace_.cursors[to]++;
    workspace_.arcs[forward] = {to, reverse, capacity, capacity};
    workspace_.arcs[reverse] = {from, forward, 0, 0};
  }

  // Aggregate equivalent face endpoints and duplicate DAG arcs. A source arc carries the
  // number of faces consumed by its task, a sink arc the number produced by its task.
  // Internal arcs retain capacity F, so only strict task reachability limits a matching.
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
        Count(source_, TaskExit(task));
      if (successor_rank_offsets_[task + 1] != successor_rank_offsets_[task])
        Count(TaskEntry(task), sink_);
      Count(TaskEntry(task), TaskExit(task));
      for (auto edge = workspace.unique_offsets[task]; edge < workspace.unique_offsets[task + 1];
           ++edge)
        Count(TaskExit(task), TaskEntry(workspace.unique_ranks[edge]));
    }
    std::partial_sum(workspace.offsets.begin(), workspace.offsets.end(), workspace.offsets.begin());
    workspace.arcs.resize(workspace.offsets.back());
    workspace.cursors.assign(workspace.offsets.begin(), workspace.offsets.end() - 1);
    workspace.marks.assign(num_tasks_, INVALID_INDEX);
    for (std::uint32_t task = 0; task < num_tasks_; ++task)
      if (workspace.incoming[task])
      {
        workspace.marks[task] = workspace.cursors[source_];
        AddEdge(source_, TaskExit(task), workspace.incoming[task]);
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
    const auto& arc = workspace_.arcs[arc_index];
    return arc.capacity - arc.residual;
  }

  void ConsumeFlowUnit(const std::uint32_t arc_index) noexcept
  {
    // Extraction consumes only forward flow; reverse residuals are no longer used.
    ++workspace_.arcs[arc_index].residual;
  }

  /**
   * Boost converts its preflow to a flow before returning. Positive-capacity arcs form a DAG
   * (source -> split tasks in topological order -> sink), so that flow has no cycles.
   * Decompose it into source-sink paths, one unit at a time, and link the face of each
   * path's source arc to the face of its sink arc. From a node that one unit has entered, flow
   * conservation leaves a forward arc with positive flow, and forward arcs form a DAG, so every
   * walk ends at the sink. Each node keeps one cursor over its arcs, which only moves past arcs
   * whose flow is 0; flows only decrease here, so they stay 0. Store the predecessor face in
   * the output vector; ExtractSlotAssignment replaces these links with slot IDs.
   */
  void ExtractChains(std::vector<std::uint32_t>& previous_face)
  {
    // Preserve both CSR inputs before writing output, which may alias either input.
    // The distinct-successor list and incoming counts have finished their build-phase use.
    workspace_.unique_ranks.assign(successor_ranks_.begin(), successor_ranks_.end());
    workspace_.incoming.assign(successor_rank_offsets_.begin(), successor_rank_offsets_.end() - 1);
    previous_face.assign(num_faces_, INVALID_INDEX);
    workspace_.cursors.assign(workspace_.offsets.begin(), workspace_.offsets.begin() + source_);

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
        auto& cursor = workspace_.cursors[node];
        const auto arc_end = workspace_.offsets[node + 1];
        while (cursor < arc_end and
               (workspace_.arcs[cursor].capacity == 0 or FlowOnArc(cursor) == 0))
          ++cursor;
        if (cursor == arc_end)
          throw std::logic_error("CBC local-face slot planner: flow decomposition found no arc.");

        const auto arc_index = cursor;
        ConsumeFlowUnit(arc_index);
        const auto to = workspace_.arcs[arc_index].to;
        if (to == sink_)
        {
          // Allocate the next unused face of this producer, in original face order.
          const auto v_face = workspace_.incoming[node / 2]++;
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
    // Every predecessor has an earlier producer task and therefore a smaller face index.
    // Its link has already been replaced by its slot ID when this face is processed. Roots
    // still receive slots in face order, exactly as in an explicit traversal of each chain.
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
  if (successor_ranks.empty())
  {
    face_slot_ids.clear();
    return {};
  }

  const auto num_tasks = successor_rank_offsets.size() - 1;
  const auto num_faces = successor_ranks.size();
  if (num_faces > INVALID_INDEX)
    throw std::length_error("CBC local face IDs exceed their 32-bit range.");

  // The unmerged network has n = 2 T + 2 nodes and 2 m = 2 (3 F + T) arcs; merging only
  // reduces the arc count, so retain these conservative fallback bounds. They keep
  // descriptors/CSR counts in uint32_t and source excess/residual pairs <= F. Boost's n, m,
  // 6*n+m budget and vector-size distance labels fit size_t on LP64; 32-bit graph count types
  // would overflow that budget. Boost uses long for relabel work and cumulative statistics:
  // the latter have no bound implied by these input limits and remain an upstream scale limit.
  // Narrower size_t/long platforms need a separate audit of Boost's budgets and counters.
  if (num_tasks > (static_cast<std::size_t>(std::numeric_limits<int>::max()) - 1) / 2 or
      2 * (std::uint64_t{3} * num_faces + num_tasks) > INVALID_INDEX)
  {
    return UseIdentityAssignment(num_faces, face_slot_ids);
  }

  static thread_local SlotPlannerWorkspace workspace;
  LocalFaceSlotPlanner planner(successor_rank_offsets, successor_ranks, workspace);
  return planner.Solve(face_slot_ids);
}

} // namespace opensn::detail
