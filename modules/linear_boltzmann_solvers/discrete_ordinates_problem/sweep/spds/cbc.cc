// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/topological_sort.hpp>

namespace opensn
{

CBC_SPDS::CBC_SPDS(const Vector3& omega,
                   const std::shared_ptr<MeshContinuum>& grid,
                   bool allow_cycles)
  : SPDS(omega, grid)
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::CBC_SPDS");

  size_t num_loc_cells = grid->local_cells.size();

  // Populate Cell Relationships
  std::vector<std::set<std::pair<int, double>>> cell_successors(num_loc_cells);
  std::set<int> location_successors;
  std::set<int> location_dependencies;

  PopulateCellRelationships(omega, location_dependencies, location_successors, cell_successors);

  location_successors_.reserve(location_successors.size());
  location_dependencies_.reserve(location_dependencies.size());

  for (auto v : location_successors)
    location_successors_.push_back(v);

  for (auto v : location_dependencies)
    location_dependencies_.push_back(v);

  // Build local cell graph
  Graph local_DG(num_loc_cells);

  // Create graph edges
  for (size_t c = 0; c < num_loc_cells; ++c) // NOLINT
    for (const auto& successor : cell_successors[c])
      boost::add_edge(c, successor.first, successor.second, local_DG);

  if (allow_cycles) // NOLINT
  {
    auto edges_to_remove = RemoveCyclicDependencies(local_DG);
    for (auto& edge_to_remove : edges_to_remove)
      local_sweep_fas_.emplace_back(edge_to_remove.first, edge_to_remove.second);
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

  // Create task list
  std::vector<std::vector<int>> global_dependencies;
  global_dependencies.resize(opensn::mpi_comm.size());
  CommunicateLocationDependencies(location_dependencies_, global_dependencies);

  constexpr auto INCOMING = FaceOrientation::INCOMING;
  constexpr auto OUTGOING = FaceOrientation::OUTGOING;

  // For each local cell create a task
  for (const auto& cell : grid_->local_cells)
  {
    const size_t num_faces = cell.faces.size();
    unsigned int num_dependencies = 0;
    unsigned int num_satisfied_downwind_deps = 0;
    unsigned int num_remote_predecessors = 0;
    unsigned int num_remote_successors = 0;
    std::vector<uint64_t> local_predecessors;
    std::vector<uint64_t> local_successors;

    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const auto& cell_face_orientation = cell_face_orientations_[cell.local_id][f];

      if (cell_face_orientation == INCOMING)
      {
        if (face.has_neighbor)
        {
          ++num_dependencies;
          if (grid->IsCellLocal(face.neighbor_id))
            local_predecessors.push_back(grid->cells[face.neighbor_id].local_id);
          else
            ++num_remote_predecessors;
        }
      }
      else if (cell_face_orientation == OUTGOING)
      {
        if (face.has_neighbor)
        {
          if (grid->IsCellLocal(face.neighbor_id))
            local_successors.push_back(grid->cells[face.neighbor_id].local_id);
          else
            ++num_remote_successors;
        }
      }
    }

    task_list_.push_back({num_dependencies,
                          num_satisfied_downwind_deps,
                          num_remote_predecessors,
                          num_remote_successors,
                          local_predecessors,
                          local_successors,
                          cell.local_id,
                          &cell,
                          false});
  }

  min_num_pool_allocator_slots_ = std::min(SimulateLocalSweep(), num_loc_cells);
}

const std::vector<Task>&
CBC_SPDS::GetTaskList() const
{
  return task_list_;
}

size_t
CBC_SPDS::SimulateLocalSweep() const
{
  std::vector<Task> simulated_task_list = task_list_;

  size_t min_num_slots = 0;
  size_t current_num_slots = 0;

  size_t num_permanent_slots = 0;

  // For each task that has remote dependencies, set aside a permanent slot
  // in the pool
  for (const auto& task : simulated_task_list)
    if ((task.num_remote_predecessors > 0) or (task.num_remote_successors > 0))
      ++num_permanent_slots;

  // Assume all remote dependencies are satisfied at the start
  for (auto& task : simulated_task_list)
    if ((task.num_remote_predecessors > 0) and
        (task.num_dependencies >= task.num_remote_predecessors))
      task.num_dependencies -= task.num_remote_predecessors;

  bool a_task_executed = true;
  while (a_task_executed)
  {
    a_task_executed = false;

    for (auto& task : simulated_task_list)
    {
      if (task.num_dependencies == 0 and (not task.completed))
      {
        // Allocate a slot for this task if it has no remote dependencies
        if ((task.num_remote_predecessors == 0) and (task.num_remote_successors == 0))
          ++current_num_slots;

        min_num_slots = std::max(min_num_slots, current_num_slots);
        a_task_executed = true;

        for (const auto& local_task_num : task.local_successors)
          --simulated_task_list[local_task_num].num_dependencies;

        task.completed = true;

        // Update predecessor task downwind dependency consumption counts
        for (const auto& local_task_num : task.local_predecessors)
        {
          auto& predecessor_task = simulated_task_list[local_task_num];
          ++predecessor_task.num_satisfied_downwind_deps;

          // Deallocate a slot if the predecessor task has satisfied all its
          // downwind dependencies and has no remote dependencies
          if ((predecessor_task.num_satisfied_downwind_deps >=
               predecessor_task.local_successors.size()) and
              (predecessor_task.num_remote_predecessors == 0) and
              (predecessor_task.num_remote_successors == 0))
            --current_num_slots;
        }

        // Deallocate a slot for this task if it has no local or remote dependencies
        if ((task.local_successors.empty()) and (task.num_remote_predecessors == 0) and
            (task.num_remote_successors == 0))
          --current_num_slots;
      }
    }
  }

  return (min_num_slots + num_permanent_slots);
}

} // namespace opensn
