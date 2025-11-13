// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include <algorithm>
#include <queue>
#include <set>
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
  ComputeNodeIndexForBoundaryAndNonLocalFaces(sdm);
  ComputeNodeIndexForParallelFaces(sdm);
  CopyFlattenNodeIndexToDevice(sdm);
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
  for (const std::vector<std::uint32_t>& level : topological_levels)
  {
    std::queue<std::uint64_t> level_idle_slots;
    // for each cell in the level
    for (const auto& cell_local_idx : level)
    {
      std::queue<std::uint64_t> cell_idle_slots;
      // push incoming face nodes in the cell's queue of idle slots
      for (std::uint32_t i_node = 0; i_node < local_node_stack.size(); ++i_node)
      {
        AAHD_DirectedEdgeNode& node = local_node_stack[i_node];
        if (node.IsInitialized() && node.downwind_node.GetCellIndex() == cell_local_idx)
        {
          cell_idle_slots.push(i_node);
          node = AAHD_DirectedEdgeNode();
        }
      }
      // build a list of outgoing nodes for the current cell
      const Cell& cell = grid.local_cells[cell_local_idx];
      std::vector<AAHD_DirectedEdgeNode> new_outgoing_nodes;
      for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
      {
        // get data of the face
        const CellFace& face = cell.faces[f];
        const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell_local_idx][f];
        const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell_local_idx][f];
        // skip if the face is not outgoing face or the neighbor cell is not a local one
        if (orientation != FaceOrientation::OUTGOING || !(face.IsNeighborLocal(&grid)))
          continue;
        // check if the face is not in the FAS
        std::uint32_t neighbor_local_idx = face.GetNeighborLocalID(&grid);
        if (fas_edges.contains(
              std::pair<std::uint32_t, std::uint32_t>(cell_local_idx, neighbor_local_idx)))
          continue;
        // for each node on the face, initialize face node
        std::uint32_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          AAHD_FaceNode upwind(cell_local_idx, f, fnode);
          AAHD_FaceNode downwind(neighbor_local_idx,
                                 face_nodal_mapping.associated_face_,
                                 face_nodal_mapping.face_node_mapping_.at(fnode));
          AAHD_DirectedEdgeNode new_node{upwind, downwind};
          new_outgoing_nodes.push_back(new_node);
        }
      }
      // insert new outgoing nodes into the stack
      for (const AAHD_DirectedEdgeNode& new_node : new_outgoing_nodes)
      {
        std::uint64_t stack_index = std::numeric_limits<std::uint64_t>::max();
        // check for cell idle slots first
        if (!cell_idle_slots.empty())
        {
          stack_index = cell_idle_slots.front();
          cell_idle_slots.pop();
          local_node_stack[stack_index] = new_node;
        }
        // then check for idle slots
        else if (!idle_slots.empty())
        {
          stack_index = idle_slots.front();
          idle_slots.pop();
          local_node_stack[stack_index] = new_node;
        }
        // otherwise, expand the stack
        else
        {
          stack_index = local_node_stack.size();
          local_node_stack.push_back(new_node);
        }
        // record the stack index
        node_tracker_.emplace(new_node.upwind_node, AAHD_NodeIndex(stack_index, true, true, false));
        node_tracker_.emplace(new_node.downwind_node,
                              AAHD_NodeIndex(stack_index, false, true, false));
        ++stack_index;
      }
      // merge cell idle slots to level idle slots
      while (!cell_idle_slots.empty())
      {
        level_idle_slots.push(cell_idle_slots.front());
        cell_idle_slots.pop();
      }
    }
    // merge level idle slots to idle slots
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
    const Cell& upwind_cell = grid.local_cells[edge.first];
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
            AAHD_FaceNode upwind(edge.first, f, fnode);
            node_tracker_.emplace(upwind, AAHD_NodeIndex(fas_node_index, true, true, true));
            AAHD_FaceNode downwind(edge.second,
                                   face_nodal_mapping.associated_face_,
                                   face_nodal_mapping.face_node_mapping_.at(fnode));
            node_tracker_.emplace(downwind, AAHD_NodeIndex(fas_node_index, false, true, true));
            fas_node_index++;
          }
        }
      }
    }
  }
  delayed_local_node_stack_size_ = fas_node_index;
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForBoundaryAndNonLocalFaces(const SpatialDiscretization& sdm)
{
  // get reference to the mesh
  const MeshContinuum& grid = *(spds_.GetGrid());
  // initialize index for boundary nodes and banks for tracking non-local nodes for each location
  std::vector<std::set<AAHD_NonLocalFaceNode>> incoming_bank, delayed_incoming_bank, outgoing_bank;
  incoming_bank.resize(spds_.GetLocationDependencies().size());
  delayed_incoming_bank.resize(spds_.GetDelayedLocationDependencies().size());
  outgoing_bank.resize(spds_.GetLocationSuccessors().size());
  std::uint64_t incremental_boundary_index = 0;
  // loop for each cell and isolate the non-local neighbor or boundary faces
  for (const Cell& cell : grid.local_cells)
  {
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
    {
      // get face data
      const CellFace& face = cell.faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell.local_id][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell.local_id][f];
      std::uint32_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
      // skip for local face
      if (face.IsNeighborLocal(&grid))
        continue;
      // construct index for boundary faces
      if (not face.has_neighbor && orientation != FaceOrientation::PARALLEL)
      {
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          node_tracker_.emplace(
            AAHD_FaceNode(cell.local_id, f, fnode),
            AAHD_NodeIndex(incremental_boundary_index, orientation == FaceOrientation::OUTGOING));
          ++incremental_boundary_index;
        }
      }
      // store non-local face nodes to the banks
      else if (face.has_neighbor)
      {
        // initialize non-local face nodes
        std::vector<AAHD_NonLocalFaceNode> nl_en_vec(num_face_nodes);
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          nl_en_vec[fnode] = AAHD_NonLocalFaceNode(cell.global_id,
                                                   cell.local_id,
                                                   f,
                                                   fnode,
                                                   face.neighbor_id,
                                                   face_nodal_mapping.associated_face_,
                                                   face_nodal_mapping.face_node_mapping_.at(fnode));
        }
        // get neighbor partition ID
        int loc = face.GetNeighborPartitionID(&grid);
        // append to incoming bank
        if (orientation == FaceOrientation::INCOMING)
        {
          int preloc = spds_.MapLocJToPrelocI(loc);
          // non-delayed incoming
          if (preloc >= 0)
          {
            incoming_bank[preloc].insert(nl_en_vec.begin(), nl_en_vec.end());
          }
          // delayed incoming
          else
          {
            preloc = -preloc - 1;
            delayed_incoming_bank[preloc].insert(nl_en_vec.begin(), nl_en_vec.end());
          }
        }
        // append to outgoing bank
        else if (orientation == FaceOrientation::OUTGOING)
        {
          int deploc = spds_.MapLocJToDeplocI(loc);
          outgoing_bank[deploc].insert(nl_en_vec.begin(), nl_en_vec.end());
        }
      }
    }
  }
  // for each bank, loop in lexicographic order and retrieve the index
  std::uint64_t node_index = 0;
  for (std::uint32_t loc_idx = 0; loc_idx < incoming_bank.size(); ++loc_idx)
  {
    for (const AAHD_NonLocalFaceNode& nl_en : incoming_bank[loc_idx])
    {
      node_tracker_.emplace(nl_en.node, AAHD_NodeIndex(node_index, false, false, false));
      node_index++;
    }
  }
  node_index = 0;
  for (std::uint32_t loc_idx = 0; loc_idx < delayed_incoming_bank.size(); ++loc_idx)
  {
    for (const AAHD_NonLocalFaceNode& nl_en : delayed_incoming_bank[loc_idx])
    {
      node_tracker_.emplace(nl_en.node, AAHD_NodeIndex(node_index, false, false, true));
      node_index++;
    }
  }
  node_index = 0;
  for (std::uint32_t loc_idx = 0; loc_idx < outgoing_bank.size(); ++loc_idx)
  {
    for (const AAHD_NonLocalFaceNode& nl_en : outgoing_bank[loc_idx])
    {
      node_tracker_.emplace(nl_en.node, AAHD_NodeIndex(node_index, true, false, false));
      node_index++;
    }
  }
  // store size
  boundary_node_size_ = incremental_boundary_index;
  std::size_t cumulative_offset = 0;
  nonlocal_incoming_node_offsets_ = {0};
  for (const std::set<AAHD_NonLocalFaceNode>& nl_en_set : incoming_bank)
  {
    nonlocal_incoming_node_sizes_.push_back(nl_en_set.size());
    cumulative_offset += nl_en_set.size();
    nonlocal_incoming_node_offsets_.push_back(cumulative_offset);
  }
  cumulative_offset = 0;
  nonlocal_delayed_incoming_node_offsets_ = {0};
  for (const std::set<AAHD_NonLocalFaceNode>& nl_en_set : delayed_incoming_bank)
  {
    nonlocal_delayed_incoming_node_sizes_.push_back(nl_en_set.size());
    cumulative_offset += nl_en_set.size();
    nonlocal_delayed_incoming_node_offsets_.push_back(cumulative_offset);
  }
  cumulative_offset = 0;
  nonlocal_outgoing_node_offsets_ = {0};
  for (const std::set<AAHD_NonLocalFaceNode>& nl_en_set : outgoing_bank)
  {
    nonlocal_outgoing_node_sizes_.push_back(nl_en_set.size());
    cumulative_offset += nl_en_set.size();
    nonlocal_outgoing_node_offsets_.push_back(cumulative_offset);
  }
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForParallelFaces(const SpatialDiscretization& sdm)
{
  // get reference to the mesh
  const MeshContinuum& grid = *(spds_.GetGrid());
  // loop for each cell and detect parallel faces
  for (const Cell& cell : grid.local_cells)
  {
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
    {
      // get face data
      const CellFace& face = cell.faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell.local_id][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell.local_id][f];
      std::uint32_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
      // skip for non-parallel face
      if (orientation != FaceOrientation::PARALLEL)
        continue;
      // construct index for parallel faces
      for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
      {
        node_tracker_.emplace(AAHD_FaceNode(cell.local_id, f, fnode), AAHD_NodeIndex());
      }
    }
  }
}

#ifndef __OPENSN_USE_CUDA__
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
