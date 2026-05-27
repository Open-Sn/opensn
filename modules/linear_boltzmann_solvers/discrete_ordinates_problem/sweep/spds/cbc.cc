// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/topological_sort.hpp>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>

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
                   const std::shared_ptr<MeshContinuum>& grid,
                   bool allow_cycles)
  : SPDS(omega, grid), id_(id), allow_cycles_(allow_cycles)
{

  const auto num_loc_cells = grid->local_cells.size();

  // Populate cell relationships
  std::vector<std::set<std::pair<std::uint32_t, double>>> cell_successors(num_loc_cells);
  std::set<int> location_successors;
  std::set<int> location_dependencies;

  PopulateCellRelationships(omega, location_dependencies, location_successors, cell_successors);

  location_successors_.reserve(location_successors.size());
  location_dependencies_.reserve(location_dependencies.size());

  for (const auto location : location_successors)
    location_successors_.push_back(location);

  for (const auto location : location_dependencies)
    location_dependencies_.push_back(location);

  // Build local cell graph
  Graph local_DG(num_loc_cells);

  // Create graph edges
  for (std::size_t c = 0; c < num_loc_cells; ++c) // NOLINT
    for (const auto& successor : cell_successors[c])
      boost::add_edge(c, successor.first, successor.second, local_DG);

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
  if (spls_.empty())
  {
    throw std::logic_error("CBC_SPDS: Cyclic dependencies found in the local cell graph.\n"
                           "Cycles need to be allowed by the calling application.");
  }

  global_dependencies_.resize(opensn::mpi_comm.size());
  CommunicateLocationDependencies(location_dependencies_, global_dependencies_);

  BuildTaskList();
}

void
CBC_SPDS::BuildTaskList()
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::BuildTaskList");

  const auto& grid = *grid_;
  constexpr auto INCOMING = FaceOrientation::INCOMING;
  constexpr auto OUTGOING = FaceOrientation::OUTGOING;

  task_list_.assign(grid.local_cells.size(), Task{});
  for (const auto& cell : grid.local_cells)
  {
    if (cell.local_id >= task_list_.size())
      throw std::logic_error("CBC_SPDS: local cell ID is outside the task-list bounds.");

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
          else if (std::find(delayed_location_dependencies_.begin(),
                             delayed_location_dependencies_.end(),
                             face.GetNeighborPartitionID(&grid)) !=
                   delayed_location_dependencies_.end())
            continue;

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
  }
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

  const int comm_size = opensn::mpi_comm.size();
  if (comm_size <= 0)
    return;

  Graph global_tdg(comm_size);
  for (int loc = 0; loc < comm_size; ++loc)
    for (const auto dep : global_dependencies_[loc])
      boost::add_edge(dep, loc, 1.0, global_tdg);

  std::vector<std::pair<int, int>> edges_to_remove(global_sweep_fas_.size() / 2);
  int edge_i = 0;
  for (auto& edge : edges_to_remove)
  {
    edge.first = global_sweep_fas_[edge_i++];
    edge.second = global_sweep_fas_[edge_i++];
  }

  for (const auto& [pred_loc, succ_loc] : edges_to_remove)
  {
    boost::remove_edge(pred_loc, succ_loc, global_tdg);

    if (succ_loc == opensn::mpi_comm.rank())
    {
      const auto it =
        std::find(location_dependencies_.begin(), location_dependencies_.end(), pred_loc);
      if (it != location_dependencies_.end())
        location_dependencies_.erase(it);
      delayed_location_dependencies_.push_back(pred_loc);
    }

    if (pred_loc == opensn::mpi_comm.rank())
      delayed_location_successors_.push_back(succ_loc);
  }

  std::vector<int> global_linear_sweep_order;
  boost::topological_sort(global_tdg, std::back_inserter(global_linear_sweep_order)); // NOLINT
  std::reverse(global_linear_sweep_order.begin(), global_linear_sweep_order.end());
  if (global_linear_sweep_order.empty())
  {
    throw std::logic_error("CBC_SPDS: Cyclic dependencies found in the global sweep graph.\n"
                           "Cycles need to be allowed by the calling application.");
  }

  BuildTaskList();
}

std::vector<CBC_SPDS::LocationEdgeWeight>
CBC_SPDS::ComputeLocalLocationEdgeWeights() const
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::ComputeLocalLocationEdgeWeights");

  constexpr double tolerance = 1.0e-16;
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
          downstream_weights[adj_cell.partition_id] += mu * mu * face.area;
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
