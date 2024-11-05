// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/aah.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/topological_sort.hpp>
#include <algorithm>

namespace opensn
{

AAH_SPDS::AAH_SPDS(int id, const Vector3& omega, const MeshContinuum& grid, bool allow_cycles)
  : SPDS(omega, grid), id_(id), allow_cycles_(allow_cycles)
{
  CALI_CXX_MARK_SCOPE("AAH_SPDS::AAH_SPDS");

  // Populate Cell Relationships
  size_t num_loc_cells = grid.local_cells.size();
  std::vector<std::set<std::pair<int, double>>> cell_successors(num_loc_cells);
  std::set<int> location_successors;
  std::set<int> location_dependencies;

  PopulateCellRelationships(omega, location_dependencies, location_successors, cell_successors);

  location_successors_.reserve(location_successors.size());
  for (auto v : location_successors)
    location_successors_.push_back(v);

  location_dependencies_.reserve(location_dependencies.size());
  for (auto v : location_dependencies)
    location_dependencies_.push_back(v);

  // Create local cell graph
  Graph local_cell_graph(num_loc_cells);

  for (auto c = 0; c < num_loc_cells; ++c)
    for (auto& successor : cell_successors[c])
      boost::add_edge(c, successor.first, successor.second, local_cell_graph);

  // Remove cycles
  if (allow_cycles)
  {
    auto edges_to_remove = RemoveCyclicDependencies(local_cell_graph);
    for (auto& edge_to_remove : edges_to_remove)
      local_sweep_fas_.emplace_back(edge_to_remove.first, edge_to_remove.second);
  }

  // Generate topological ordering
  spls_.item_id.clear();
  boost::topological_sort(local_cell_graph, std::back_inserter(spls_.item_id));
  std::reverse(spls_.item_id.begin(), spls_.item_id.end());
  if (spls_.item_id.empty())
  {
    throw std::logic_error("AAH_SPDS: Cyclic dependencies found in the local cell graph.\n"
                           "Cycles need to be allowed by the calling application.");
  }

  // Generate location-to-location dependencies
  global_dependencies_.resize(opensn::mpi_comm.size());
  CommunicateLocationDependencies(location_dependencies_, global_dependencies_);
}

void
AAH_SPDS::BuildGlobalSweepFAS()
{
  assert(not global_dependencies_.empty());

  CALI_CXX_MARK_SCOPE("AAH_SPDS::BuildGlobalSweepFAS");

  // Create global sweep graph
  Graph global_tdg(opensn::mpi_comm.size());

  for (int loc = 0; loc < opensn::mpi_comm.size(); ++loc)
    for (int dep : global_dependencies_[loc])
      boost::add_edge(dep, loc, 1.0, global_tdg);

  // Remove cycles and generate the feedback arc set (FAS). The FAS is the list of edges that must
  // be removed from the graph to make it acyclic.
  if (allow_cycles_)
  {
    auto edges_to_remove = RemoveCyclicDependencies(global_tdg);
    for (const auto& [e0, e1] : edges_to_remove)
    {
      global_sweep_fas_.emplace_back(e0);
      global_sweep_fas_.emplace_back(e1);
    }
  }
}

void
AAH_SPDS::BuildGlobalSweepTDG()
{
  CALI_CXX_MARK_SCOPE("AAH_SPDS::BuildGlobalSweepTDG");

  // Create graph
  Graph global_tdg(opensn::mpi_comm.size());

  for (int loc = 0; loc < opensn::mpi_comm.size(); ++loc)
    for (int dep : global_dependencies_[loc])
      boost::add_edge(dep, loc, 1.0, global_tdg);

  // De-serialize edges
  std::vector<std::pair<int, int>> edges_to_remove;
  edges_to_remove.resize(global_sweep_fas_.size() / 2, std::pair<int, int>(0, 0));
  int i = 0;
  for (auto& edge : edges_to_remove)
  {
    edge.first = global_sweep_fas_[i++];
    edge.second = global_sweep_fas_[i++];
  }

  // Remove edges
  for (auto& edge_to_remove : edges_to_remove)
  {
    int rlocI = edge_to_remove.first;
    int locI = edge_to_remove.second;

    boost::remove_edge(rlocI, locI, global_tdg);

    if (locI == opensn::mpi_comm.rank())
    {
      auto dependent_location =
        std::find(location_dependencies_.begin(), location_dependencies_.end(), rlocI);
      location_dependencies_.erase(dependent_location);
      delayed_location_dependencies_.push_back(rlocI);
    }

    if (rlocI == opensn::mpi_comm.rank())
      delayed_location_successors_.push_back(locI);
  }

  // Generate topological ordering
  std::vector<size_t> global_linear_sweep_order;
  boost::topological_sort(global_tdg, std::back_inserter(global_linear_sweep_order));
  std::reverse(global_linear_sweep_order.begin(), global_linear_sweep_order.end());
  if (global_linear_sweep_order.empty())
  {
    throw std::logic_error("AAH_SPDS: Cyclic dependencies found in the global sweep graph.\n"
                           "Cycles need to be allowed by the calling application.");
  }

  // Rank to global_tdg id mapping
  std::vector<int> global_order_mapping(opensn::mpi_comm.size(), -1);
  for (int k = 0; k < opensn::mpi_comm.size(); ++k)
  {
    int loc = global_linear_sweep_order[k];
    global_order_mapping[loc] = k;
  }

  // Determine sweep order ranks
  int abs_max_rank = 0;
  std::vector<int> global_sweep_order_rank(opensn::mpi_comm.size(), -1);
  for (int k = 0; k < opensn::mpi_comm.size(); ++k)
  {
    int loc = global_linear_sweep_order[k];
    if (global_dependencies_[loc].empty())
      global_sweep_order_rank[k] = 0;
    else
    {
      int max_rank = -1;
      for (auto dep_loc : global_dependencies_[loc])
      {
        if (dep_loc < 0)
          continue;

        int dep_mapped_index = global_order_mapping[dep_loc];
        if (global_sweep_order_rank[dep_mapped_index] > max_rank)
          max_rank = global_sweep_order_rank[dep_mapped_index];
      }
      global_sweep_order_rank[k] = max_rank + 1;
      if ((max_rank + 1) > abs_max_rank)
        abs_max_rank = max_rank + 1;
    }
  }

  // Generate TDG
  for (int rank = 0; rank <= abs_max_rank; ++rank)
  {
    STDG stdg;
    for (int k = 0; k < opensn::mpi_comm.size(); ++k)
      if (global_sweep_order_rank[k] == rank)
        stdg.item_id.push_back(global_linear_sweep_order[k]);
    global_sweep_planes_.push_back(stdg);
  }
}

} // namespace opensn
