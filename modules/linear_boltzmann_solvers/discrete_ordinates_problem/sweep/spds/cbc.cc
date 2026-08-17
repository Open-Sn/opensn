// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <unordered_map>

namespace opensn
{

namespace
{

// Bounds exact SCC optimization to 2^16 states (roughly 0.6 MiB of DP storage).
constexpr std::size_t MAX_EXACT_FAS_STATES = std::size_t{1} << 16;

std::uint64_t
PackDirectedEdge(const std::uint32_t upstream_id, const std::uint32_t downstream_id) noexcept
{
  return (static_cast<std::uint64_t>(upstream_id) << 32) |
         static_cast<std::uint64_t>(downstream_id);
}

std::vector<std::pair<Vertex, Vertex>>
FindBackwardEdges(const Graph& graph,
                  const std::vector<Vertex>& component,
                  const std::vector<Vertex>& ordering)
{
  using boost::make_iterator_range;

  constexpr auto invalid_position = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> position(boost::num_vertices(graph), invalid_position);
  for (std::size_t i = 0; i < ordering.size(); ++i)
    position[ordering[i]] = i;

  std::vector<std::pair<Vertex, Vertex>> feedback_arc_set;
  for (const auto upstream : component)
    for (const auto edge : make_iterator_range(boost::out_edges(upstream, graph)))
    {
      const auto downstream = boost::target(edge, graph);
      if (position[downstream] != invalid_position and position[downstream] <= position[upstream])
        feedback_arc_set.emplace_back(upstream, downstream);
    }
  return feedback_arc_set;
}

} // namespace

CBC_SPDS::CBC_SPDS(int id,
                   const Vector3& omega,
                   const std::shared_ptr<MeshContinuum>& grid,
                   const SPDSFaceNeighborInfoVec& face_neighbor_info,
                   bool allow_cycles)
  : SPDS(omega, grid), id_(id), allow_cycles_(allow_cycles)
{

  const auto num_loc_cells = grid->local_cells.size();

  // Populate cell relationships
  std::vector<std::vector<std::pair<std::uint32_t, double>>> cell_successors(num_loc_cells);
  std::vector<int> location_successors;
  std::vector<int> location_dependencies;

  PopulateCellRelationships(
    omega, face_neighbor_info, location_dependencies, location_successors, cell_successors);

  location_successors_ = std::move(location_successors);
  location_dependencies_ = std::move(location_dependencies);

  // Build local cell graph
  Graph local_DG(num_loc_cells);

  /*
   * One lagged local edge stores one face's angular flux. Weighting by face-node count makes the
   * FAS objective proportional to the lagged storage and copy work for this SPDS.
   */
  for (const auto& cell : grid->local_cells)
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if (cell_face_orientations_[cell.local_id][f] == FaceOrientation::OUTGOING and
          face.has_neighbor and face.IsNeighborLocal(grid.get()))
        boost::add_edge(cell.local_id,
                        face.GetNeighborLocalID(grid.get()),
                        static_cast<double>(face.vertex_ids.size()),
                        local_DG);
    }

  if (allow_cycles) // NOLINT
  {
    auto edges_to_remove = RemoveCyclicDependencies(local_DG);
    for (const auto& edge_to_remove : edges_to_remove)
    {
      const auto upwind = static_cast<std::uint32_t>(edge_to_remove.first);
      const auto downwind = static_cast<std::uint32_t>(edge_to_remove.second);
      local_sweep_fas_.emplace_back(upwind, downwind);
      delayed_local_dependency_set_.insert(PackDirectedEdge(upwind, downwind));
    }
  }

  // Generate topological sorting
  spls_.clear();
  boost::topological_sort(local_DG, std::back_inserter(spls_)); // NOLINT
  std::reverse(spls_.begin(), spls_.end());

  delayed_location_dependency_flags_.assign(static_cast<std::size_t>(opensn::mpi_comm.size()), 0);
  BuildTaskList();
}

std::vector<std::pair<Vertex, Vertex>>
CBC_SPDS::FindApproxMinimumFAS(const Graph& graph, const std::vector<Vertex>& component)
{
  using boost::make_iterator_range;

  const auto num_vertices = boost::num_vertices(graph);
  std::vector<bool> active(num_vertices, false);
  for (const auto vertex : component)
    active[vertex] = true;

  const auto edge_weights = boost::get(boost::edge_weight, graph);
  const auto has_active_successor = [&](const Vertex vertex)
  {
    return std::ranges::any_of(make_iterator_range(boost::out_edges(vertex, graph)),
                               [&](const auto edge) { return active[boost::target(edge, graph)]; });
  };
  const auto has_active_predecessor = [&](const Vertex vertex)
  {
    return std::ranges::any_of(make_iterator_range(boost::in_edges(vertex, graph)),
                               [&](const auto edge) { return active[boost::source(edge, graph)]; });
  };
  const auto weighted_degree_difference = [&](const Vertex vertex)
  {
    double difference = 0.0;
    for (const auto edge : make_iterator_range(boost::out_edges(vertex, graph)))
      if (active[boost::target(edge, graph)])
        difference += edge_weights[edge];
    for (const auto edge : make_iterator_range(boost::in_edges(vertex, graph)))
      if (active[boost::source(edge, graph)])
        difference -= edge_weights[edge];
    return difference;
  };

  std::vector<Vertex> prefix;
  std::vector<Vertex> suffix;
  prefix.reserve(component.size());
  suffix.reserve(component.size());
  auto remaining = component.size();

  while (remaining != 0)
  {
    bool removed_vertex = true;
    while (removed_vertex)
    {
      removed_vertex = false;
      for (const auto vertex : component)
        if (active[vertex] and not has_active_successor(vertex))
        {
          active[vertex] = false;
          suffix.push_back(vertex);
          --remaining;
          removed_vertex = true;
        }
    }

    removed_vertex = true;
    while (removed_vertex)
    {
      removed_vertex = false;
      for (const auto vertex : component)
        if (active[vertex] and not has_active_predecessor(vertex))
        {
          active[vertex] = false;
          prefix.push_back(vertex);
          --remaining;
          removed_vertex = true;
        }
    }

    double maximum_difference = -std::numeric_limits<double>::infinity();
    auto selected = boost::graph_traits<Graph>::null_vertex();
    for (const auto vertex : component)
      if (active[vertex])
      {
        const double difference = weighted_degree_difference(vertex);
        if (difference > maximum_difference)
        {
          maximum_difference = difference;
          selected = vertex;
        }
      }

    if (selected != boost::graph_traits<Graph>::null_vertex())
    {
      active[selected] = false;
      prefix.push_back(selected);
      --remaining;
    }
  }

  std::vector<Vertex> ordering;
  ordering.reserve(component.size());
  ordering.insert(ordering.end(), prefix.begin(), prefix.end());
  ordering.insert(ordering.end(), suffix.rbegin(), suffix.rend());

  return FindBackwardEdges(graph, component, ordering);
}

std::vector<std::pair<Vertex, Vertex>>
CBC_SPDS::FindExactMinimumFAS(const Graph& graph, const std::vector<Vertex>& component)
{
  using boost::make_iterator_range;

  const auto num_component_vertices = component.size();
  const auto num_states = std::size_t{1} << num_component_vertices;
  std::vector<std::vector<double>> edge_weights(num_component_vertices,
                                                std::vector<double>(num_component_vertices, 0.0));
  std::vector<std::size_t> component_index(boost::num_vertices(graph), num_component_vertices);
  for (std::size_t i = 0; i < num_component_vertices; ++i)
    component_index[component[i]] = i;

  const auto weight_map = boost::get(boost::edge_weight, graph);
  for (std::size_t i = 0; i < num_component_vertices; ++i)
    for (const auto edge : make_iterator_range(boost::out_edges(component[i], graph)))
    {
      const auto j = component_index[boost::target(edge, graph)];
      if (j < num_component_vertices and i != j)
        edge_weights[i][j] += weight_map[edge];
    }

  std::vector<double> minimum_cost(num_states, std::numeric_limits<double>::infinity());
  std::vector<std::uint8_t> last_vertex(num_states, 0);
  minimum_cost[0] = 0.0;

  // If v is last in S, precisely the edges v -> u with u in S\{v} are backward.
  for (std::size_t state = 1; state < num_states; ++state)
    for (std::size_t v = 0; v < num_component_vertices; ++v)
    {
      const auto vertex_bit = std::size_t{1} << v;
      if ((state & vertex_bit) == 0)
        continue;

      const auto prefix = state ^ vertex_bit;
      double cost = minimum_cost[prefix];
      for (std::size_t u = 0; u < num_component_vertices; ++u)
        if (prefix & (std::size_t{1} << u))
          cost += edge_weights[v][u];

      if (cost < minimum_cost[state] or
          (cost == minimum_cost[state] and component[v] > component[last_vertex[state]]))
      {
        minimum_cost[state] = cost;
        last_vertex[state] = static_cast<std::uint8_t>(v);
      }
    }

  std::vector<Vertex> ordering(num_component_vertices);
  auto state = num_states - 1;
  for (std::size_t position = num_component_vertices; position-- > 0;)
  {
    const auto v = static_cast<std::size_t>(last_vertex[state]);
    ordering[position] = component[v];
    state ^= std::size_t{1} << v;
  }

  return FindBackwardEdges(graph, component, ordering);
}

std::vector<std::pair<Vertex, Vertex>>
CBC_SPDS::RemoveCyclicDependencies(Graph& graph)
{
  const auto num_vertices = boost::num_vertices(graph);
  std::vector<int> component_ids(num_vertices, -1);
  const int num_components =
    boost::strong_components(graph,
                             boost::make_iterator_property_map(
                               component_ids.begin(), boost::get(boost::vertex_index, graph)));

  std::vector<std::vector<Vertex>> components(num_components);
  for (Vertex vertex = 0; vertex < num_vertices; ++vertex)
    components[component_ids[vertex]].push_back(vertex);

  std::vector<std::pair<Vertex, Vertex>> feedback_arc_set;
  for (const auto& component : components)
  {
    if (component.size() == 1)
    {
      const auto vertex = component.front();
      if (boost::edge(vertex, vertex, graph).second)
        feedback_arc_set.emplace_back(vertex, vertex);
      continue;
    }

    const bool exact_state_space_fits =
      component.size() < std::numeric_limits<std::size_t>::digits and
      (std::size_t{1} << component.size()) <= MAX_EXACT_FAS_STATES;
    auto component_fas = exact_state_space_fits ? FindExactMinimumFAS(graph, component)
                                                : FindApproxMinimumFAS(graph, component);
    feedback_arc_set.insert(feedback_arc_set.end(), component_fas.begin(), component_fas.end());
  }

  std::ranges::sort(feedback_arc_set);
  const auto unique_end = std::ranges::unique(feedback_arc_set).begin();
  feedback_arc_set.erase(unique_end, feedback_arc_set.end());
  for (const auto& [upstream, downstream] : feedback_arc_set)
    boost::remove_edge(upstream, downstream, graph);

  return feedback_arc_set;
}

void
CBC_SPDS::BuildTaskList()
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::BuildTaskList");

  const auto& grid = *grid_;
  constexpr auto INCOMING = FaceOrientation::INCOMING;
  constexpr auto OUTGOING = FaceOrientation::OUTGOING;

  task_list_.assign(grid.local_cells.size(), Task{});
  initial_task_dependencies_.assign(task_list_.size(), 0);
  initial_ready_tasks_.clear();
  for (const auto& cell : grid.local_cells)
  {
    const auto num_faces = cell.faces.size();
    unsigned int num_dependencies = 0;
    std::vector<std::uint32_t> successors;
    successors.reserve(num_faces);

    for (std::size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const auto& orientation = cell_face_orientations_[cell.local_id][f];

      if (orientation == INCOMING)
      {
        if (face.has_neighbor)
        {
          if (face.IsNeighborLocal(&grid))
          {
            const auto upwind_local_id = grid.cells[face.neighbor_id].local_id;
            if (IsDelayedLocalDependency(upwind_local_id, cell.local_id))
              continue;
          }
          else
          {
            const auto neighbor_location = face.GetNeighborPartitionID(&grid);
            if (neighbor_location >= 0 and
                delayed_location_dependency_flags_[static_cast<std::size_t>(neighbor_location)])
              continue;
          }

          ++num_dependencies;
        }
      }
      else if (orientation == OUTGOING)
      {
        if (face.has_neighbor and face.IsNeighborLocal(&grid))
        {
          const auto successor_local_id = grid.cells[face.neighbor_id].local_id;
          if (IsDelayedLocalDependency(cell.local_id, successor_local_id))
            continue;

          successors.push_back(successor_local_id);
        }
      }
    }

    task_list_[cell.local_id] = {
      num_dependencies, std::move(successors), cell.local_id, &cell, false};
    initial_task_dependencies_[cell.local_id] = num_dependencies;
  }

  for (std::uint32_t task = 0; task < initial_task_dependencies_.size(); ++task)
    if (initial_task_dependencies_[task] == 0)
      initial_ready_tasks_.push_back(task);
}

const std::vector<Task>&
CBC_SPDS::GetTaskList() const
{
  return task_list_;
}

bool
CBC_SPDS::IsDelayedLocalDependency(const std::uint32_t upwind_local_id,
                                   const std::uint32_t downwind_local_id) const noexcept
{
  return delayed_local_dependency_set_.contains(
    PackDirectedEdge(upwind_local_id, downwind_local_id));
}

void
CBC_SPDS::SetGlobalDependencies(std::vector<std::vector<int>> global_dependencies)
{
  global_dependencies_ = std::move(global_dependencies);
}

void
CBC_SPDS::BuildGlobalSweepFAS()
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::BuildGlobalSweepFAS");

  const int comm_size = opensn::mpi_comm.size();
  Graph global_tdg(comm_size);

  for (int loc = 0; loc < comm_size; ++loc)
  {
    for (const auto dep : global_dependencies_[loc])
    {
      double weight = 1.0;
      if (not global_edge_weights_.empty())
      {
        const auto edge_weight = global_edge_weights_.find(
          PackDirectedEdge(static_cast<std::uint32_t>(dep), static_cast<std::uint32_t>(loc)));
        if (edge_weight != global_edge_weights_.end() and edge_weight->second > 0.0)
          weight = edge_weight->second;
      }
      boost::add_edge(dep, loc, weight, global_tdg);
    }
  }

  global_sweep_fas_.clear();
  if (allow_cycles_)
  {
    const auto edges_to_remove = RemoveCyclicDependencies(global_tdg);
    for (const auto& [upwind, downwind] : edges_to_remove)
    {
      global_sweep_fas_.push_back(static_cast<int>(upwind));
      global_sweep_fas_.push_back(static_cast<int>(downwind));
    }
  }
}

void
CBC_SPDS::ApplyGlobalSweepFAS()
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::ApplyGlobalSweepFAS");

  delayed_location_dependencies_.clear();
  delayed_location_successors_.clear();
  delayed_location_dependency_flags_.assign(static_cast<std::size_t>(opensn::mpi_comm.size()), 0);

  std::vector<std::pair<int, int>> edges_to_remove(global_sweep_fas_.size() / 2);
  int edge_i = 0;
  for (auto& edge : edges_to_remove)
  {
    edge.first = global_sweep_fas_[edge_i++];
    edge.second = global_sweep_fas_[edge_i++];
  }

  for (const auto& [pred_loc, succ_loc] : edges_to_remove)
  {
    if (succ_loc == opensn::mpi_comm.rank())
    {
      const auto it =
        std::find(location_dependencies_.begin(), location_dependencies_.end(), pred_loc);
      if (it != location_dependencies_.end())
        location_dependencies_.erase(it);
      delayed_location_dependencies_.push_back(pred_loc);
      delayed_location_dependency_flags_[static_cast<std::size_t>(pred_loc)] = 1;
    }

    if (pred_loc == opensn::mpi_comm.rank())
      delayed_location_successors_.push_back(succ_loc);
  }

  BuildTaskList();
}

std::vector<CBC_SPDS::LocationEdgeWeight>
CBC_SPDS::ComputeLocalLocationEdgeWeights() const
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::ComputeLocalLocationEdgeWeights");

  constexpr double tolerance = FACE_ORIENTATION_TOLERANCE;
  std::unordered_map<int, double> downstream_weights;

  for (const auto& cell : grid_->local_cells)
  {
    const auto& face_orientations = cell_face_orientations_[cell.local_id];
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if (face.has_neighbor and not face.IsNeighborLocal(grid_.get()) and
          face_orientations[f] == FaceOrientation::OUTGOING)
      {
        const double mu = omega_.Dot(face.normal);
        if (mu > tolerance)
        {
          const auto& adj_cell = grid_->cells[face.neighbor_id];
          downstream_weights[adj_cell.partition_id] += static_cast<double>(face.vertex_ids.size());
        }
      }
    }
  }

  std::vector<LocationEdgeWeight> edge_weights;
  edge_weights.reserve(downstream_weights.size());
  for (const auto& [downstream_location, weight] : downstream_weights)
    edge_weights.push_back({opensn::mpi_comm.rank(), downstream_location, weight});

  return edge_weights;
}

void
CBC_SPDS::SetGlobalEdgeWeights(std::span<const LocationEdgeWeight> edge_weights)
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::SetGlobalEdgeWeights");

  global_edge_weights_.clear();
  global_edge_weights_.reserve(edge_weights.size());
  for (const auto& edge_weight : edge_weights)
  {
    if (edge_weight.weight <= 0.0)
      continue;

    global_edge_weights_[PackDirectedEdge(
      static_cast<std::uint32_t>(edge_weight.upstream_location),
      static_cast<std::uint32_t>(edge_weight.downstream_location))] += edge_weight.weight;
  }
}

} // namespace opensn
