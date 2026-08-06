// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <algorithm>
#include <queue>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace opensn
{

AAHD_FLUDSCommonData::AAHD_FLUDSCommonData(
  const SPDS& spds,
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SpatialDiscretization& sdm)
  : FLUDSCommonData(spds, grid_nodal_mappings)
{
  ComputeNodeIndexForNonDelayedLocalFaces(sdm);
  ComputeNodeIndexForDelayedLocalFaces(sdm);
  ComputeNodeIndexForNonLocalFaces(sdm);
  ComputeNodeIndexForParallelFaces(sdm);
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForNonDelayedLocalFaces(const SpatialDiscretization& sdm)
{
  // get the topological levels
  const std::vector<std::vector<std::uint32_t>>& topological_levels =
    spds_.GetLevelizedLocalSubgrid();
  if (topological_levels.empty())
    return;
  // initialize the set of FAS edges
  std::set<std::pair<std::uint32_t, std::uint32_t>> fas_edges(spds_.GetLocalSweepFAS().begin(),
                                                              spds_.GetLocalSweepFAS().end());
  // for each level
  const MeshContinuum& grid = *(spds_.GetGrid());
  std::vector<AAHD_DirectedEdgeNode> local_node_stack;
  std::queue<std::uint64_t> idle_slots;

  // Store the stack slots that provide incoming data to each downstream cell. This lets us find a
  // cell's incoming slots directly instead of searching the entire stack. A slot remains valid
  // until its downstream cell is processed, because that is the only time the slot is released.
  std::unordered_map<std::uint32_t, std::vector<std::uint64_t>> downwind_index;
  downwind_index.reserve(grid.GetLocalCellCount());

  // Place a directed edge node at the stack_index and register it in downwind_index.
  auto RegisterNode = [&](std::uint64_t stack_index, const AAHD_DirectedEdgeNode& new_node)
  {
    downwind_index[new_node.downwind_node.GetCellIndex()].push_back(stack_index);
    node_tracker_.emplace(new_node.upwind_node,
                          AAHD_NodeIndex(stack_index, true, false, true, false));
    node_tracker_.emplace(new_node.downwind_node,
                          AAHD_NodeIndex(stack_index, false, false, true, false));
  };

  for (const std::vector<std::uint32_t>& level : topological_levels)
  {
    std::queue<std::uint64_t> level_idle_slots;
    // for each cell in the level
    for (const auto& cell_local_idx : level)
    {
      std::queue<std::uint64_t> cell_idle_slots;
      // Collect idle slots for this cell
      auto it = downwind_index.find(cell_local_idx);
      if (it != downwind_index.end())
      {
        for (std::uint64_t pos : it->second)
        {
          cell_idle_slots.push(pos);
          local_node_stack[pos] = AAHD_DirectedEdgeNode();
        }
        it->second.clear();
      }
      // Build a list of outgoing nodes for the current cell
      const Cell& cell = grid.GetLocalCell(cell_local_idx);
      std::vector<AAHD_DirectedEdgeNode> new_outgoing_nodes;
      for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
      {
        const CellFace& face = cell.faces[f];
        const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell_local_idx][f];
        const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell_local_idx][f];
        // Skip if the face is not outgoing or its neighbor is not local
        if (orientation != FaceOrientation::OUTGOING || !(face.IsNeighborLocal(&grid)))
          continue;
        // Check if the face is not in the FAS
        std::uint32_t neighbor_local_idx = face.GetNeighborLocalID(&grid);
        if (fas_edges.contains(
              std::pair<std::uint32_t, std::uint32_t>(cell_local_idx, neighbor_local_idx)))
          continue;
        std::uint32_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          FaceNode upwind(cell_local_idx, f, fnode);
          FaceNode downwind(neighbor_local_idx,
                            face_nodal_mapping.associated_face_,
                            face_nodal_mapping.face_node_mapping_.at(fnode));
          AAHD_DirectedEdgeNode new_node{upwind, downwind};
          new_outgoing_nodes.push_back(new_node);
        }
      }
      // Insert new outgoing nodes into the stack
      for (const AAHD_DirectedEdgeNode& new_node : new_outgoing_nodes)
      {
        std::uint64_t stack_index = std::numeric_limits<std::uint64_t>::max();
        if (!cell_idle_slots.empty())
        {
          stack_index = cell_idle_slots.front();
          cell_idle_slots.pop();
          local_node_stack[stack_index] = new_node;
        }
        else if (!idle_slots.empty())
        {
          stack_index = idle_slots.front();
          idle_slots.pop();
          local_node_stack[stack_index] = new_node;
        }
        else
        {
          stack_index = local_node_stack.size();
          local_node_stack.push_back(new_node);
        }
        RegisterNode(stack_index, new_node);
      }
      while (!cell_idle_slots.empty())
      {
        level_idle_slots.push(cell_idle_slots.front());
        cell_idle_slots.pop();
      }
    }
    while (!level_idle_slots.empty())
    {
      idle_slots.push(level_idle_slots.front());
      level_idle_slots.pop();
    }
  }
  local_node_stack_size_ = local_node_stack.size();
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForDelayedLocalFaces(const SpatialDiscretization& sdm)
{
  // get reference to the grid
  const MeshContinuum& grid = *(spds_.GetGrid());
  // record the FAS edges
  std::uint64_t fas_node_index = 0;
  const std::vector<std::pair<std::uint32_t, std::uint32_t>>& fas_edges = spds_.GetLocalSweepFAS();
  for (const auto& edge : fas_edges)
  {
    const Cell& upwind_cell = grid.GetLocalCell(edge.first);
    for (std::uint32_t f = 0; f < upwind_cell.faces.size(); ++f)
    {
      const CellFace& face = upwind_cell.faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[edge.first][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[edge.first][f];
      // check if the face is outgoing and its neighbor is the cell in the FAS edge
      if (orientation == FaceOrientation::OUTGOING && face.IsNeighborLocal(&grid))
      {
        std::uint32_t neighbor_local_idx = face.GetNeighborLocalID(&grid);
        if (std::cmp_equal(neighbor_local_idx, edge.second))
        {
          // record the address of all the nodes
          std::uint32_t num_face_nodes = sdm.GetCellMapping(upwind_cell).GetNumFaceNodes(f);
          for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
          {
            FaceNode upwind(edge.first, f, fnode);
            node_tracker_.emplace(upwind, AAHD_NodeIndex(fas_node_index, true, false, true, true));
            FaceNode downwind(edge.second,
                              face_nodal_mapping.associated_face_,
                              face_nodal_mapping.face_node_mapping_.at(fnode));
            node_tracker_.emplace(downwind,
                                  AAHD_NodeIndex(fas_node_index, false, false, true, true));
            fas_node_index++;
          }
        }
      }
    }
  }
  delayed_local_node_stack_size_ = fas_node_index;
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForNonLocalFaces(const SpatialDiscretization& sdm)
{
  const MeshContinuum& grid = *(spds_.GetGrid());
  // Sort each bank before assigning node indices so that index ordering is deterministic. Face
  // nodes are unique, so the banks do not require duplicate removal.
  std::vector<std::vector<AAHD_NonLocalFaceNode>> incoming_bank, delayed_incoming_bank,
    outgoing_bank;
  incoming_bank.resize(spds_.GetLocationDependencies().size());
  delayed_incoming_bank.resize(spds_.GetDelayedLocationDependencies().size());
  outgoing_bank.resize(spds_.GetLocationSuccessors().size());
  std::uint64_t incremental_boundary_index = 0;
  // Loop over cells and isolate the non-local neighbor or boundary faces
  for (const auto& cell : grid.GetLocalCells())
  {
    for (std::uint32_t f = 0; f < cell->faces.size(); ++f)
    {
      const CellFace& face = cell->faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell->local_id][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell->local_id][f];
      std::uint32_t num_face_nodes = sdm.GetCellMapping(*cell).GetNumFaceNodes(f);
      // Skip for local face and boundary face
      if (face.IsNeighborLocal(&grid) or not face.has_neighbor)
        continue;
      // Store non-local face nodes
      std::vector<AAHD_NonLocalFaceNode> nl_en_vec(num_face_nodes);
      for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
      {
        nl_en_vec[fnode] = AAHD_NonLocalFaceNode(cell->global_id,
                                                 cell->local_id,
                                                 f,
                                                 fnode,
                                                 face.neighbor_id,
                                                 face_nodal_mapping.associated_face_,
                                                 face_nodal_mapping.face_node_mapping_.at(fnode));
      }
      int loc = face.GetNeighborPartitionID(&grid);
      if (orientation == FaceOrientation::INCOMING)
      {
        int preloc = spds_.MapLocJToPrelocI(loc);
        // Non-delayed incoming
        if (preloc >= 0)
        {
          incoming_bank[preloc].insert(
            incoming_bank[preloc].end(), nl_en_vec.begin(), nl_en_vec.end());
        }
        // Delayed incoming
        else
        {
          preloc = -preloc - 1;
          delayed_incoming_bank[preloc].insert(
            delayed_incoming_bank[preloc].end(), nl_en_vec.begin(), nl_en_vec.end());
        }
      }
      // Append to outgoing bank
      else if (orientation == FaceOrientation::OUTGOING)
      {
        int deploc = spds_.MapLocJToDeplocI(loc);
        outgoing_bank[deploc].insert(
          outgoing_bank[deploc].end(), nl_en_vec.begin(), nl_en_vec.end());
      }
    }
  }
  // Sort each bank
  for (auto& bank : incoming_bank)
    std::sort(bank.begin(), bank.end());
  for (auto& bank : delayed_incoming_bank)
    std::sort(bank.begin(), bank.end());
  for (auto& bank : outgoing_bank)
    std::sort(bank.begin(), bank.end());
  // For each bank, loop in lexicographic order and retrieve the index
  std::uint64_t node_index = 0;
  for (std::uint32_t loc_idx = 0; loc_idx < incoming_bank.size(); ++loc_idx)
  {
    for (const AAHD_NonLocalFaceNode& nl_en : incoming_bank[loc_idx])
    {
      node_tracker_.emplace(nl_en.node, AAHD_NodeIndex(node_index, false, false, false, false));
      node_index++;
    }
  }
  node_index = 0;
  for (std::uint32_t loc_idx = 0; loc_idx < delayed_incoming_bank.size(); ++loc_idx)
  {
    for (const AAHD_NonLocalFaceNode& nl_en : delayed_incoming_bank[loc_idx])
    {
      node_tracker_.emplace(nl_en.node, AAHD_NodeIndex(node_index, false, false, false, true));
      node_index++;
    }
  }
  node_index = 0;
  for (std::uint32_t loc_idx = 0; loc_idx < outgoing_bank.size(); ++loc_idx)
  {
    for (const AAHD_NonLocalFaceNode& nl_en : outgoing_bank[loc_idx])
    {
      node_tracker_.emplace(nl_en.node, AAHD_NodeIndex(node_index, true, false, false, false));
      node_index++;
    }
  }
  // Store size
  std::size_t cumulative_offset = 0;
  nonlocal_incoming_node_offsets_ = {0};
  for (const auto& bank : incoming_bank)
  {
    nonlocal_incoming_node_sizes_.push_back(bank.size());
    cumulative_offset += bank.size();
    nonlocal_incoming_node_offsets_.push_back(cumulative_offset);
  }
  cumulative_offset = 0;
  nonlocal_delayed_incoming_node_offsets_ = {0};
  for (const auto& bank : delayed_incoming_bank)
  {
    nonlocal_delayed_incoming_node_sizes_.push_back(bank.size());
    cumulative_offset += bank.size();
    nonlocal_delayed_incoming_node_offsets_.push_back(cumulative_offset);
  }
  cumulative_offset = 0;
  nonlocal_outgoing_node_offsets_ = {0};
  for (const auto& bank : outgoing_bank)
  {
    nonlocal_outgoing_node_sizes_.push_back(bank.size());
    cumulative_offset += bank.size();
    nonlocal_outgoing_node_offsets_.push_back(cumulative_offset);
  }
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForParallelFaces(const SpatialDiscretization& sdm)
{
  // get reference to the mesh
  const MeshContinuum& grid = *(spds_.GetGrid());
  // loop for each cell and detect parallel faces
  for (const auto& cell : grid.GetLocalCells())
  {
    for (std::uint32_t f = 0; f < cell->faces.size(); ++f)
    {
      // get face data
      const CellFace& face = cell->faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell->local_id][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell->local_id][f];
      std::uint32_t num_face_nodes = sdm.GetCellMapping(*cell).GetNumFaceNodes(f);
      // skip for non-parallel face
      if (orientation != FaceOrientation::PARALLEL)
        continue;
      // construct index for parallel faces
      for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
      {
        node_tracker_.emplace(FaceNode(cell->local_id, f, fnode), AAHD_NodeIndex());
      }
    }
  }
}

#ifndef __OPENSN_WITH_GPU__
void
AAHD_FLUDSCommonData::CopyFlattenNodeIndexToDevice(const SpatialDiscretization& sdm)
{
}

void
AAHD_FLUDSCommonData::DeallocateDeviceMemory()
{
}
#endif

AAHD_FLUDSCommonData::~AAHD_FLUDSCommonData()
{
  DeallocateDeviceMemory();
}

} // namespace opensn
