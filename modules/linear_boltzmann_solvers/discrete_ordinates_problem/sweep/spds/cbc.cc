// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_graph_validation.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "framework/graphs/feedback_arc_set.h"
#include "caliper/cali.h"
#include <boost/graph/topological_sort.hpp>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <set>
#include <unordered_map>
#include <utility>

namespace opensn
{

namespace
{

std::uint64_t
PackDirectedEdge(const std::uint32_t upstream_id, const std::uint32_t downstream_id) noexcept
{
  return (static_cast<std::uint64_t>(upstream_id) << 32) |
         static_cast<std::uint64_t>(downstream_id);
}

} // namespace

CBC_SPDS::CBC_SPDS(int id,
                   const Vector3& omega,
                   const std::shared_ptr<Mesh>& grid,
                   const SPDSFaceNeighborInfoVec& face_neighbor_info,
                   bool allow_cycles)
  : SPDS(omega, grid), id_(id), allow_cycles_(allow_cycles)
{

  const auto num_loc_cells = grid->GetLocalCellCount();

  // Populate cell relationships
  std::vector<std::vector<std::pair<std::uint32_t, double>>> cell_successors(num_loc_cells);
  std::vector<int> location_successors;
  std::vector<int> location_dependencies;

  PopulateCellRelationships(
    omega, face_neighbor_info, location_dependencies, location_successors, cell_successors);

  location_successors_ = std::move(location_successors);
  location_dependencies_ = std::move(location_dependencies);

  // Build local cell graph
  Graph local_cell_graph(num_loc_cells);

  // With cycles allowed, each outgoing local face is weighted by its node count, so the
  // weighted FAS favors smaller lagged storage and copy work (exact only for small SCCs).
  if (allow_cycles) // NOLINT
  {
    std::vector<LocalFaceEdge> face_edges;
    for (const auto& cell_ptr : grid->GetLocalCells())
    {
      const auto& cell = *cell_ptr;
      for (std::size_t f = 0; f < cell.faces.size(); ++f)
      {
        const auto& face = cell.faces[f];
        if (cell_face_orientations_[cell.local_id][f] == FaceOrientation::OUTGOING and
            face.has_neighbor and face.IsNeighborLocal(grid.get()))
          face_edges.push_back({static_cast<std::uint32_t>(cell.local_id),
                                static_cast<std::uint32_t>(face.GetNeighborLocalID(grid.get())),
                                static_cast<double>(face.vertex_ids.size())});
      }
    }

    local_sweep_fas_ = RemoveLocalCellCycles(local_cell_graph, face_edges);
    for (const auto& [upwind, downwind] : local_sweep_fas_)
      delayed_local_dependency_set_.insert(PackDirectedEdge(upwind, downwind));
  }
  else
    for (std::size_t cell_id = 0; cell_id < cell_successors.size(); ++cell_id)
      for (const auto& [successor_id, weight] : cell_successors[cell_id])
        boost::add_edge(cell_id, successor_id, weight, local_cell_graph);

  // Generate topological sorting
  spls_.clear();
  try
  {
    // An explicit color map avoids analyzer false positives in Boost's shared_array map.
    std::vector<boost::default_color_type> colors(num_loc_cells, boost::white_color);
    boost::topological_sort(
      local_cell_graph,
      std::back_inserter(spls_),
      boost::color_map(boost::make_iterator_property_map(
        colors.begin(), boost::get(boost::vertex_index, local_cell_graph)))); // NOLINT
  }
  catch (const boost::not_a_dag&)
  {
    throw std::logic_error("CBC_SPDS: Cyclic dependencies found in the local cell graph.\n"
                           "Cycles need to be allowed by the calling application.");
  }
  std::reverse(spls_.begin(), spls_.end());

  if (allow_cycles_)
    delayed_location_dependency_flags_.assign(static_cast<std::size_t>(opensn::mpi_comm.size()), 0);
  BuildTaskList();
}

std::vector<std::pair<std::uint32_t, std::uint32_t>>
CBC_SPDS::RemoveLocalCellCycles(Graph& graph, const std::vector<LocalFaceEdge>& face_edges)
{
  // The greedy FAS of a component with more than MAX_EXACT_FAS_SCC_SIZE vertices lists a pair once
  // per parallel edge, so parallel face edges are merged first. The summed weight is the one both
  // FAS methods use for parallel edges, so the removed pairs are the same.
  std::unordered_map<std::uint64_t, std::size_t> pair_index;
  std::vector<LocalFaceEdge> cell_edges;
  cell_edges.reserve(face_edges.size());
  for (const auto& edge : face_edges)
  {
    assert(edge.upwind != edge.downwind);
    const auto [it, inserted] =
      pair_index.try_emplace(PackDirectedEdge(edge.upwind, edge.downwind), cell_edges.size());
    if (inserted)
      cell_edges.push_back(edge);
    else
      cell_edges[it->second].weight += edge.weight;
  }
  for (const auto& edge : cell_edges)
    boost::add_edge(edge.upwind, edge.downwind, edge.weight, graph);

  const auto edges_to_remove = RemoveCyclicDependencies(graph);
  std::vector<std::pair<std::uint32_t, std::uint32_t>> removed_pairs;
  removed_pairs.reserve(edges_to_remove.size());
  for (const auto& [upwind, downwind] : edges_to_remove)
    removed_pairs.emplace_back(static_cast<std::uint32_t>(upwind),
                               static_cast<std::uint32_t>(downwind));
#ifndef NDEBUG
  const std::set<std::pair<std::uint32_t, std::uint32_t>> unique_pairs(removed_pairs.begin(),
                                                                       removed_pairs.end());
  assert(unique_pairs.size() == removed_pairs.size());
#endif
  return removed_pairs;
}

void
CBC_SPDS::BuildTaskList()
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::BuildTaskList");

  const auto& grid = *grid_;
  constexpr auto INCOMING = FaceOrientation::INCOMING;
  constexpr auto OUTGOING = FaceOrientation::OUTGOING;

  task_list_.assign(grid.GetLocalCellCount(), Task{});
  initial_task_dependency_counts_.assign(task_list_.size(), 0);
  initial_ready_tasks_.clear();
  for (const auto& cell_ptr : grid.GetLocalCells())
  {
    const auto& cell = *cell_ptr;
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
            const auto upwind_local_id = grid.GetGlobalCell(face.neighbor_id).local_id;
            if (allow_cycles_ and IsDelayedLocalDependency(upwind_local_id, cell.local_id))
              continue;
          }
          else
          {
            const auto neighbor_location = face.GetNeighborPartitionID(&grid);
            if (allow_cycles_ and neighbor_location >= 0 and
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
          const auto successor_local_id = grid.GetGlobalCell(face.neighbor_id).local_id;
          if (allow_cycles_ and IsDelayedLocalDependency(cell.local_id, successor_local_id))
            continue;

          successors.push_back(successor_local_id);
        }
      }
    }

    task_list_[cell.local_id] = {
      num_dependencies, std::move(successors), cell.local_id, &cell, false};
    initial_task_dependency_counts_[cell.local_id] = num_dependencies;
  }

  for (std::uint32_t task = 0; task < initial_task_dependency_counts_.size(); ++task)
    if (initial_task_dependency_counts_[task] == 0)
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

std::string
CBC_SPDS::CheckGlobalCellCycles() const
{
  if (allow_cycles_)
    throw std::logic_error("CBC global no-lag validation requires allow_cycles=False");
  const auto components =
    detail::GetCBCCyclicRankComponents(location_dependencies_, location_successors_);
  if (std::none_of(components.begin(), components.end(), [](int c) { return c >= 0; }))
    return {};
  // No cell-level data or traversal is needed outside the nontrivial quotient SCCs. Empty ranks
  // still participate in elimination, the verdict and the witness broadcasts.
  if (components[mpi_comm.rank()] < 0)
    return detail::CheckCBCGlobalCellGraph({}, {}, {}, {}, components);
  std::vector<std::uint64_t> ids(task_list_.size());
  std::vector<std::vector<std::uint32_t>> successors(task_list_.size());
  std::vector<detail::CBCBoundaryEdge> incoming, outgoing;
  for (const auto& task : task_list_)
  {
    const auto& cell = *task.cell_ptr;
    ids[cell.local_id] = cell.global_id;
    successors[cell.local_id] = task.successors;
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if (not face.has_neighbor or face.IsNeighborLocal(grid_.get()))
        continue;
      detail::CBCBoundaryEdge edge{
        cell.local_id, face.GetNeighborPartitionID(grid_.get()), face.neighbor_id};
      if (cell_face_orientations_[cell.local_id][f] == FaceOrientation::INCOMING)
        incoming.push_back(edge);
      else if (cell_face_orientations_[cell.local_id][f] == FaceOrientation::OUTGOING)
        outgoing.push_back(edge);
    }
  }
  return detail::CheckCBCGlobalCellGraph(ids, successors, incoming, outgoing, components);
}

void
CBC_SPDS::SetGlobalSweepFAS(std::vector<int> edges)
{
  if (edges.size() % 2 != 0)
    throw std::logic_error("CBC FAS requires pairs of rank endpoints");
  std::set<std::pair<int, int>> pairs;
  for (std::size_t i = 0; i < edges.size(); i += 2)
    if (edges[i] < 0 or edges[i + 1] < 0 or edges[i] >= mpi_comm.size() or
        edges[i + 1] >= mpi_comm.size() or edges[i] == edges[i + 1] or
        not pairs.emplace(edges[i], edges[i + 1]).second)
      throw std::logic_error("CBC FAS has invalid or duplicate rank endpoints");
  global_sweep_fas_ = std::move(edges);
}

void
CBC_SPDS::SetGlobalDependencies(std::vector<std::vector<int>> global_dependencies)
{
  if (global_dependencies.size() != static_cast<std::size_t>(mpi_comm.size()))
    throw std::logic_error("CBC global dependencies require one row per MPI rank");
  for (std::size_t rank = 0; rank < global_dependencies.size(); ++rank)
  {
    std::set<int> unique;
    for (const int dependency : global_dependencies[rank])
      if (dependency < 0 or dependency >= mpi_comm.size() or std::cmp_equal(dependency, rank) or
          not unique.insert(dependency).second)
        throw std::logic_error("CBC global dependencies have an invalid or duplicate rank");
  }
  global_dependencies_ = std::move(global_dependencies);
}

void
CBC_SPDS::BuildGlobalSweepFAS()
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::BuildGlobalSweepFAS");

  // One edge per (dependency, location) pair. Each rank's dependencies are distinct other ranks
  // (PopulateCellRelationships sorts and deduplicates them), so the graph has no parallel edges
  // and no self loops.
  const int comm_size = opensn::mpi_comm.size();
  if (global_dependencies_.size() != static_cast<std::size_t>(comm_size))
    throw std::logic_error("CBC global dependencies have not been initialized");
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

  // Validate against the original incoming set, including cuts from an earlier application.
  // Build a temporary set first so a rejected replacement leaves task state untouched.
  auto dependencies = location_dependencies_;
  dependencies.insert(dependencies.end(),
                      delayed_location_dependencies_.begin(),
                      delayed_location_dependencies_.end());
  std::ranges::sort(dependencies);
  dependencies.erase(std::unique(dependencies.begin(), dependencies.end()), dependencies.end());
  for (std::size_t i = 0; i < global_sweep_fas_.size(); i += 2)
  {
    const int pred = global_sweep_fas_[i], succ = global_sweep_fas_[i + 1];
    if (not global_dependencies_.empty())
    {
      const auto& row = global_dependencies_.at(succ);
      if (std::find(row.begin(), row.end(), pred) == row.end())
        throw std::logic_error("CBC FAS removes a rank edge absent from the global graph");
    }
    else if ((succ == mpi_comm.rank() and
              std::find(dependencies.begin(), dependencies.end(), pred) == dependencies.end()) or
             (pred == mpi_comm.rank() and
              std::find(location_successors_.begin(), location_successors_.end(), succ) ==
                location_successors_.end()))
      throw std::logic_error("CBC FAS removes a rank edge absent from this rank's graph");
  }
  location_dependencies_ = std::move(dependencies);
  delayed_location_dependencies_.clear();
  delayed_location_successors_.clear();
  delayed_location_dependency_flags_.assign(static_cast<std::size_t>(opensn::mpi_comm.size()), 0);

  for (std::size_t edge_i = 0; edge_i < global_sweep_fas_.size(); edge_i += 2)
  {
    const int pred_loc = global_sweep_fas_[edge_i];
    const int succ_loc = global_sweep_fas_[edge_i + 1];
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

  for (const auto& cell_ptr : grid_->GetLocalCells())
  {
    const auto& cell = *cell_ptr;
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
          const auto& adj_cell = grid_->GetGlobalCell(face.neighbor_id);
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

  for (const auto& edge : edge_weights)
    if (edge.upstream_location < 0 or edge.downstream_location < 0 or
        edge.upstream_location >= mpi_comm.size() or edge.downstream_location >= mpi_comm.size() or
        edge.upstream_location == edge.downstream_location or not std::isfinite(edge.weight))
      throw std::logic_error("CBC rank graph weight has an invalid endpoint/value");
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
