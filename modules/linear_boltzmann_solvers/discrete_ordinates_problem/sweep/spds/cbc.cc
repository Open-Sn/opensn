// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc_slot_planner.h"
#include "framework/logging/log.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <boost/graph/topological_sort.hpp>
#include <numeric>
#include <stdexcept>

namespace opensn
{

void
CBC_SPDS::BuildTaskGraph()
{
  constexpr auto INCOMING = FaceOrientation::INCOMING;
  constexpr auto OUTGOING = FaceOrientation::OUTGOING;

  const auto num_loc_cells = grid_->local_cells.size();
  task_list_.assign(num_loc_cells, Task{});
  task_successor_rank_offsets_.assign(num_loc_cells + 1, 0);
  task_successor_ranks_.clear();
  task_successor_ranks_.reserve(num_loc_cells * 4);

  for (std::size_t rank = 0; rank < topo_order_.size(); ++rank)
  {
    const auto& cell = grid_->local_cells[topo_order_[rank]];
    unsigned int num_dependencies = 0;
    std::vector<std::uint32_t> successors;

    successors.reserve(cell.faces.size());
    task_successor_rank_offsets_[rank] = static_cast<std::uint32_t>(task_successor_ranks_.size());
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      const auto& orientation = cell_face_orientations_[cell.local_id][f];

      if (orientation == INCOMING and face.has_neighbor)
        ++num_dependencies;
      else if ((orientation == OUTGOING) and (face.has_neighbor) and
               (face.IsNeighborLocal(grid_.get())))
      {
        const auto successor_local_id = grid_->cells[face.neighbor_id].local_id;
        successors.push_back(successor_local_id);
        task_successor_ranks_.push_back(topo_rank_by_cell_local_id_[successor_local_id]);
      }
    }

    task_list_[cell.local_id] =
      Task{num_dependencies, std::move(successors), cell.local_id, &cell, false};
  }
  task_successor_rank_offsets_.back() = static_cast<std::uint32_t>(task_successor_ranks_.size());
}

void
CBC_SPDS::BuildLocalFaceTaskGraph()
{
  // Each outgoing local face becomes one directed-face task.
  // The task is keyed by the producer-cell topological rank, the consumer-cell topological
  // rank, and the face-node count needed later when CBC/CBCD size the compact slot bank.
  const auto num_loc_cells = grid_->local_cells.size();
  cell_face_offsets_.assign(num_loc_cells + 1, 0);
  size_t total_num_faces = 0;
  for (const auto& cell : grid_->local_cells)
  {
    cell_face_offsets_[cell.local_id] = static_cast<std::uint32_t>(total_num_faces);
    total_num_faces += cell.faces.size();
  }
  cell_face_offsets_.back() = static_cast<std::uint32_t>(total_num_faces);
  outgoing_local_face_task_ids_.assign(total_num_faces, INVALID_LOCAL_FACE_TASK_ID);
  incoming_local_face_task_ids_.assign(total_num_faces, INVALID_LOCAL_FACE_TASK_ID);

  producer_cell_face_offsets_.assign(num_loc_cells + 1, 0);
  local_face_producer_ranks_.clear();
  local_face_consumer_ranks_.clear();
  local_face_node_counts_.clear();
  max_local_face_node_count_ = 0;

  for (std::size_t producer_rank = 0; producer_rank < topo_order_.size(); ++producer_rank)
  {
    producer_cell_face_offsets_[producer_rank] =
      static_cast<std::uint32_t>(local_face_producer_ranks_.size());

    const auto producer_cell_local_id = topo_order_[producer_rank];
    const auto& cell = grid_->local_cells[producer_cell_local_id];
    const auto& face_orientations = cell_face_orientations_[producer_cell_local_id];

    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      const auto& orientation = face_orientations[f];
      if ((orientation != FaceOrientation::OUTGOING) or (not face.IsNeighborLocal(grid_.get())))
        continue;

      const auto consumer_cell_local_id = face.GetNeighborLocalID(grid_.get());
      const auto consumer_face_id =
        static_cast<std::uint16_t>(face.GetNeighborAdjacentFaceIndex(grid_.get()));
      const auto num_face_nodes = static_cast<std::uint32_t>(face.vertex_ids.size());
      max_local_face_node_count_ =
        std::max(max_local_face_node_count_, static_cast<std::size_t>(num_face_nodes));

      const auto face_task_id = static_cast<std::uint32_t>(local_face_producer_ranks_.size());
      local_face_producer_ranks_.push_back(static_cast<std::uint32_t>(producer_rank));
      local_face_consumer_ranks_.push_back(topo_rank_by_cell_local_id_[consumer_cell_local_id]);
      local_face_node_counts_.push_back(static_cast<std::uint16_t>(num_face_nodes));
      outgoing_local_face_task_ids_[cell_face_offsets_[producer_cell_local_id] + f] = face_task_id;
      incoming_local_face_task_ids_[cell_face_offsets_[consumer_cell_local_id] + consumer_face_id] =
        face_task_id;
    }
  }

  producer_cell_face_offsets_.back() =
    static_cast<std::uint32_t>(local_face_producer_ranks_.size());
  local_face_slot_ids_.resize(local_face_producer_ranks_.size());
  std::iota(local_face_slot_ids_.begin(), local_face_slot_ids_.end(), std::uint32_t{0});
}

void
CBC_SPDS::UpdateLocalFaceSlotLayout()
{
  // The slot planner only decides which faces may share one slot. The physical storage bank is
  // then sized slot-by-slot by taking the maximum face-node extent over each slot chain.
  local_face_slot_node_counts_.assign(max_num_local_psi_slots_, std::uint16_t{0});
  local_face_slot_node_offsets_.assign(max_num_local_psi_slots_ + 1, std::uint32_t{0});
  total_local_face_slot_nodes_ = 0;

  bool is_identity_layout = max_num_local_psi_slots_ == local_face_slot_ids_.size();
  for (std::size_t face_task_id = 0;
       is_identity_layout and face_task_id < local_face_slot_ids_.size();
       ++face_task_id)
    is_identity_layout = local_face_slot_ids_[face_task_id] == face_task_id;

  if (is_identity_layout)
  {
    for (std::size_t slot_id = 0; slot_id < local_face_node_counts_.size(); ++slot_id)
    {
      local_face_slot_node_counts_[slot_id] = local_face_node_counts_[slot_id];
      local_face_slot_node_offsets_[slot_id] =
        static_cast<std::uint32_t>(total_local_face_slot_nodes_);
      total_local_face_slot_nodes_ += local_face_node_counts_[slot_id];
    }
    local_face_slot_node_offsets_.back() = static_cast<std::uint32_t>(total_local_face_slot_nodes_);
    return;
  }

  for (std::size_t face_task_id = 0; face_task_id < local_face_slot_ids_.size(); ++face_task_id)
  {
    const auto slot_id = local_face_slot_ids_[face_task_id];
    assert(slot_id < local_face_slot_node_counts_.size());
    local_face_slot_node_counts_[slot_id] =
      std::max(local_face_slot_node_counts_[slot_id], local_face_node_counts_[face_task_id]);
  }

  for (std::size_t slot_id = 0; slot_id < local_face_slot_node_counts_.size(); ++slot_id)
  {
    local_face_slot_node_offsets_[slot_id] =
      static_cast<std::uint32_t>(total_local_face_slot_nodes_);
    total_local_face_slot_nodes_ += local_face_slot_node_counts_[slot_id];
  }
  local_face_slot_node_offsets_.back() = static_cast<std::uint32_t>(total_local_face_slot_nodes_);
}

CBC_SPDS::CBC_SPDS(const Vector3& omega,
                   const std::shared_ptr<MeshContinuum>& grid,
                   bool allow_cycles)
  : SPDS(omega, grid)
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::CBC_SPDS");

  size_t num_loc_cells = grid->local_cells.size();

  std::vector<std::set<std::pair<std::uint32_t, double>>> cell_successors(num_loc_cells);
  std::set<int> location_successors;
  std::set<int> location_dependencies;

  PopulateCellRelationships(omega, location_dependencies, location_successors, cell_successors);

  location_successors_.reserve(location_successors.size());
  location_dependencies_.reserve(location_dependencies.size());

  for (auto v : location_successors)
    location_successors_.push_back(v);

  for (auto v : location_dependencies)
    location_dependencies_.push_back(v);

  Graph local_DG(num_loc_cells);

  for (size_t c = 0; c < num_loc_cells; ++c) // NOLINT
    for (const auto& successor : cell_successors[c])
      boost::add_edge(c, successor.first, successor.second, local_DG);

  if (allow_cycles) // NOLINT
  {
    auto edges_to_remove = RemoveCyclicDependencies(local_DG);
    for (const auto& [u, v] : edges_to_remove)
      local_sweep_fas_.emplace_back(u, v);
  }

  spls_.clear();
  boost::topological_sort(local_DG, std::back_inserter(spls_)); // NOLINT
  std::reverse(spls_.begin(), spls_.end());
  if (spls_.empty())
  {
    throw std::logic_error("CBC_SPDS: Cyclic dependencies found in the local cell graph.\n"
                           "Cycles need to be allowed by the calling application.");
  }

  topo_order_.assign(spls_.begin(), spls_.end());
  topo_rank_by_cell_local_id_.assign(num_loc_cells, 0);
  for (std::size_t rank = 0; rank < topo_order_.size(); ++rank)
    topo_rank_by_cell_local_id_[topo_order_[rank]] = static_cast<std::uint32_t>(rank);

  std::vector<std::vector<int>> global_dependencies(opensn::mpi_comm.size());
  CommunicateLocationDependencies(location_dependencies_, global_dependencies);
  BuildTaskGraph();
  BuildLocalFaceTaskGraph();

  max_num_local_psi_slots_ = local_face_producer_ranks_.size();
  UpdateLocalFaceSlotLayout();
}

const std::vector<Task>&
CBC_SPDS::GetTaskList() const noexcept
{
  return task_list_;
}

void
CBC_SPDS::ComputeMaxNumLocalPsiSlots()
{
  CALI_CXX_MARK_SCOPE("CBC_SPDS::ComputeMaxNumLocalPsiSlots");

  if (task_list_.empty())
  {
    max_num_local_psi_slots_ = 0;
    local_face_slot_ids_.clear();
    UpdateLocalFaceSlotLayout();
    return;
  }

  if (local_face_producer_ranks_.empty())
  {
    max_num_local_psi_slots_ = 0;
    local_face_slot_ids_.clear();
    UpdateLocalFaceSlotLayout();
    return;
  }

  // Solve the exact minimum chain cover of the local-face reuse poset, then turn that chain
  // decomposition into a static slot assignment and compact slot-bank layout.
  const auto result = detail::ComputeLocalFaceSlotPlan(task_successor_rank_offsets_,
                                                       task_successor_ranks_,
                                                       local_face_producer_ranks_,
                                                       local_face_consumer_ranks_,
                                                       producer_cell_face_offsets_,
                                                       local_face_slot_ids_);
  max_num_local_psi_slots_ = result.slot_count;
  UpdateLocalFaceSlotLayout();
  if (result.verifier_rejected)
    opensn::log.LogAllWarning()
      << "CBC_SPDS::ComputeMaxNumLocalPsiSlots: local cell-face slot assignment verifier rejected "
      << " the computed slot count; falling back to the identity assignment "
      << " (one slot per local directed face, no reuse).";
}

std::uint32_t
CBC_SPDS::GetOutgoingLocalFaceTaskID(const std::uint32_t cell_local_id,
                                     const unsigned int face_id) const noexcept
{
  return outgoing_local_face_task_ids_[cell_face_offsets_[cell_local_id] + face_id];
}

std::uint32_t
CBC_SPDS::GetIncomingLocalFaceTaskID(const std::uint32_t cell_local_id,
                                     const unsigned int face_id) const noexcept
{
  return incoming_local_face_task_ids_[cell_face_offsets_[cell_local_id] + face_id];
}

} // namespace opensn
