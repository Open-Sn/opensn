// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/runtime.h"
#include <boost/graph/topological_sort.hpp>
#include <algorithm>

namespace opensn
{

CBC_SPDS::CBC_SPDS(const Vector3& omega,
                   const std::shared_ptr<Mesh>& grid,
                   const SPDSFaceNeighborInfoVec& face_neighbor_info,
                   bool allow_cycles)
  : SPDS(omega, grid)
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
  Graph local_DG(num_loc_cells);

  // Create graph edges
  for (std::size_t c = 0; c < num_loc_cells; ++c) // NOLINT
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
  constexpr auto INCOMING = FaceOrientation::INCOMING;
  constexpr auto OUTGOING = FaceOrientation::OUTGOING;

  // For each local cell create a task
  task_list_.reserve(grid_->GetLocalCellCount());
  for (std::uint32_t cell_local_id = 0; cell_local_id < grid_->GetLocalCellCount(); ++cell_local_id)
  {
    const auto& cell = grid_->GetLocalCell(cell_local_id);
    const auto cell_faces = grid_->GetCellFaces(cell_local_id);
    const auto num_faces = grid_->GetCellFaceCount(cell_local_id);
    unsigned int num_dependencies = 0;
    std::vector<std::uint32_t> successors;
    successors.reserve(num_faces);

    for (std::size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell_faces[f];
      const auto& orientation = cell_face_orientations_[cell_local_id][f];

      if (orientation == INCOMING)
      {
        if (face.has_neighbor)
          ++num_dependencies;
      }
      else if (orientation == OUTGOING)
      {
        if (face.has_neighbor and grid->IsCellLocal(face.neighbor_id))
          successors.push_back(grid->MapCellGlobalID2LocalID(face.neighbor_id));
      }
    }

    task_list_.push_back({num_dependencies, successors, cell_local_id, false});
  }
}

const std::vector<Task>&
CBC_SPDS::GetTaskList() const
{
  return task_list_;
}

} // namespace opensn
