// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <boost/graph/topological_sort.hpp>
#include <algorithm>
#include <map>

namespace opensn
{

AAH_SPDS::AAH_SPDS(int id,
                   const Vector3& omega,
                   const std::shared_ptr<MeshContinuum> grid,
                   const SPDSFaceNeighborInfoVec& face_neighbor_info,
                   bool allow_cycles,
                   bool use_gpus)
  : SPDS(omega, grid), id_(id), allow_cycles_(allow_cycles)
{

  // Populate cell relationships
  size_t num_loc_cells = grid->local_cells.size();
  std::vector<std::vector<std::pair<std::uint32_t, double>>> cell_successors(num_loc_cells);
  std::vector<int> location_successors;
  std::vector<int> location_dependencies;

  PopulateCellRelationships(
    omega, face_neighbor_info, location_dependencies, location_successors, cell_successors);

  location_successors_ = std::move(location_successors);
  location_dependencies_ = std::move(location_dependencies);

  // Create local cell graph
  Graph local_cell_graph(num_loc_cells);

  for (size_t c = 0; c < num_loc_cells; ++c) // NOLINT
    for (const auto& successor : cell_successors[c])
      boost::add_edge(c, successor.first, successor.second, local_cell_graph);

  // Remove cycles
  if (allow_cycles) // NOLINT
  {
    auto edges_to_remove = RemoveCyclicDependencies(local_cell_graph);
    for (auto& edge_to_remove : edges_to_remove)
      local_sweep_fas_.emplace_back(edge_to_remove.first, edge_to_remove.second);
  }

  // Generate topological ordering
  spls_.clear();
  boost::topological_sort(local_cell_graph, std::back_inserter(spls_)); // NOLINT
  std::reverse(spls_.begin(), spls_.end());
  if (spls_.empty())
  {
    throw std::logic_error("AAH_SPDS: Cyclic dependencies found in the local cell graph.\n"
                           "Cycles need to be allowed by the calling application.");
  }

  // Generate levelized spls
  int max_level = 0;
  std::vector<int> levels(num_vertices(local_cell_graph), 0);
  for (auto& v : spls_)
  {
    for (auto ei = out_edges(v, local_cell_graph); ei.first != ei.second; ++ei.first)
    {
      auto successor = target(*ei.first, local_cell_graph);
      levels[successor] = std::max(levels[successor], levels[v] + 1);
      max_level = std::max(max_level, levels[successor]);
    }
  }
  levelized_spls_.resize(max_level + 1);
  for (auto v = 0; v < num_vertices(local_cell_graph); ++v)
    levelized_spls_[levels[v]].push_back(v);

  // Regenerate spls to match levelized spls
  spls_.clear();
  for (auto& level : levelized_spls_)
    for (auto& cell : level)
      spls_.push_back(cell);

  // Copy levelized spls data to GPU
  if (use_gpus)
    CopySPLSDataOnDevice();
}

AAH_SPDS::GlobalSweepMetadata
AAH_SPDS::BuildGlobalSweepMetadata()
{
  if (global_sweep_metadata_built_)
    throw std::logic_error("AAH_SPDS: Global sweep metadata has already been built.");
  if (not global_edges_initialized_)
    throw std::logic_error("AAH_SPDS: Global sweep edges are not initialized.");

  const int comm_size = opensn::mpi_comm.size();
  Graph global_tdg(comm_size);
  for (const auto& edge : global_edges_)
  {
    const auto dep = edge.dependency;
    const auto loc = edge.location;
    if (dep < 0 or dep >= comm_size or loc < 0 or loc >= comm_size)
      throw std::logic_error("AAH_SPDS: Invalid edge in the global sweep graph.");
    boost::add_edge(dep, loc, edge.weight, global_tdg);
  }
  std::vector<GlobalSweepEdge>().swap(global_edges_);

  std::vector<std::pair<Vertex, Vertex>> edges_to_remove;
  if (allow_cycles_)
    edges_to_remove = RemoveCyclicDependencies(global_tdg);

  std::vector<int> global_linear_sweep_order;
  boost::topological_sort(global_tdg, std::back_inserter(global_linear_sweep_order)); // NOLINT
  std::reverse(global_linear_sweep_order.begin(), global_linear_sweep_order.end());
  if (global_linear_sweep_order.size() != static_cast<std::size_t>(comm_size))
    throw std::logic_error("AAH_SPDS: Cyclic dependencies found in the global sweep graph.\n"
                           "Cycles need to be allowed by the calling application.");

  int max_level = 0;
  std::vector<int> levels(comm_size, 0);
  for (const int loc : global_linear_sweep_order)
  {
    for (auto [edge, edge_end] = boost::in_edges(loc, global_tdg); edge != edge_end; ++edge)
      levels[loc] = std::max(levels[loc], levels[boost::source(*edge, global_tdg)] + 1);
    max_level = std::max(max_level, levels[loc]);
  }

  GlobalSweepMetadata metadata;
  metadata.location_depths.resize(comm_size);
  metadata.delayed_dependencies.resize(comm_size);
  metadata.delayed_successors.resize(comm_size);
  for (int loc = 0; loc < comm_size; ++loc)
    metadata.location_depths[loc] = max_level - levels[loc] + 1;
  for (const auto& [dep_vertex, loc_vertex] : edges_to_remove)
  {
    if (std::cmp_greater_equal(dep_vertex, comm_size) or
        std::cmp_greater_equal(loc_vertex, comm_size))
      throw std::logic_error("AAH_SPDS: Invalid feedback edge in the global sweep graph.");
    const auto dep = static_cast<int>(dep_vertex);
    const auto loc = static_cast<int>(loc_vertex);
    metadata.delayed_dependencies[loc].push_back(dep);
    metadata.delayed_successors[dep].push_back(loc);
  }

  global_sweep_metadata_built_ = true;
  return metadata;
}

void
AAH_SPDS::SetGlobalSweepMetadata(int location_depth,
                                 std::vector<int> delayed_dependencies,
                                 std::vector<int> delayed_successors)
{
  if (global_sweep_metadata_set_)
    throw std::logic_error("AAH_SPDS: Global sweep metadata has already been set.");
  if (location_depth < 1)
    throw std::logic_error("AAH_SPDS: Invalid location depth.");

  for (const int dep : delayed_dependencies)
  {
    const auto it = std::find(location_dependencies_.begin(), location_dependencies_.end(), dep);
    if (it == location_dependencies_.end())
      throw std::logic_error("AAH_SPDS: Feedback edge is not a local dependency.");
    location_dependencies_.erase(it);
  }
  for (const int successor : delayed_successors)
    if (std::find(location_successors_.begin(), location_successors_.end(), successor) ==
        location_successors_.end())
      throw std::logic_error("AAH_SPDS: Feedback edge is not a local successor.");

  location_depth_ = location_depth;
  delayed_location_dependencies_ = std::move(delayed_dependencies);
  delayed_location_successors_ = std::move(delayed_successors);
  global_sweep_metadata_set_ = true;
}

std::vector<std::pair<int, double>>
AAH_SPDS::ComputeLocalLocationEdgeWeights() const
{
  std::map<int, double> row;
  constexpr double tolerance = FACE_ORIENTATION_TOLERANCE;

  for (const auto& cell : grid_->local_cells)
  {
    const auto& face_orientations = cell_face_orientations_[cell.local_id];
    std::size_t f = 0;
    for (const auto& face : cell.faces)
    {
      if (face.has_neighbor and not face.IsNeighborLocal(grid_.get()) and
          face_orientations[f] == FaceOrientation::OUTGOING)
      {
        const double mu = omega_.Dot(face.normal);
        if (mu > tolerance)
        {
          const auto& adj_cell = grid_->cells[face.neighbor_id];
          row[adj_cell.partition_id] += mu * mu * face.area;
        }
      }
      ++f;
    }
  }

  return {row.begin(), row.end()};
}

#ifndef __OPENSN_WITH_GPU__
void
AAH_SPDS::CopySPLSDataOnDevice()
{
}

void
AAH_SPDS::FreeDeviceData()
{
}
#endif

AAH_SPDS::~AAH_SPDS()
{
  FreeDeviceData();
}

} // namespace opensn
