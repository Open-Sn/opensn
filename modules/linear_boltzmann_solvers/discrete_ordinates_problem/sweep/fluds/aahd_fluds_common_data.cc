// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
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
  const std::vector<std::vector<int>>& topological_levels = spds_.GetLevelizedLocalSubgrid();
  if (topological_levels.empty())
    return;
  std::int32_t max_level = topological_levels.size() - 1u;
  std::set<std::pair<int, int>> fas_edges(spds_.GetLocalSweepFAS().begin(),
                                          spds_.GetLocalSweepFAS().end());
  const MeshContinuum& grid = *(spds_.GetGrid());
  // initialize the local face node stack for non-delay local cells
  std::vector<AAHD_DirectedEdgeNode> non_delay_local_node_stack(1, AAHD_DirectedEdgeNode());
  for (std::int32_t level = max_level; level != 0; --level)
  {
    // sort cell local index in the level
    std::set<std::uint32_t> sorted_level(topological_levels[level].begin(),
                                         topological_levels[level].end());
    // compute new face nodes that going to the cells in the current level to the face node stack
    std::vector<AAHD_DirectedEdgeNode> new_non_delay_local_nodes;
    for (const std::uint32_t& cell_local_idx : sorted_level)
    {
      // check for each cell face
      const Cell& cell = grid.local_cells[cell_local_idx];
      for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
      {
        // get data of the face
        const CellFace& face = cell.faces[f];
        const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell_local_idx][f];
        const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell_local_idx][f];
        // skip if the face is not incoming face or the neighbor cell is not a local one
        if (orientation != FaceOrientation::INCOMING || !(face.IsNeighborLocal(&grid)))
          continue;
        // check if the face is not in the FAS
        std::uint32_t neighbor_local_idx = face.GetNeighborLocalID(&grid);
        if (fas_edges.contains(std::pair<int, int>(neighbor_local_idx, cell_local_idx)))
          continue;
        // for each node on the face, initialize face node
        std::uint32_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          AAHD_FaceNode downwind(cell_local_idx, f, fnode);
          AAHD_FaceNode upwind(neighbor_local_idx,
                               face_nodal_mapping.associated_face_,
                               face_nodal_mapping.face_node_mapping_.at(fnode));
          new_non_delay_local_nodes.push_back({upwind, downwind});
        }
      }
    }
    // append new nodes into the face node stack and record its index in the stack
    for (std::uint64_t stack_index = 0;
         const AAHD_DirectedEdgeNode& new_node : new_non_delay_local_nodes)
    {
      // skip if the element is already initialized
      while (stack_index != non_delay_local_node_stack.size() &&
             non_delay_local_node_stack[stack_index].IsInitialized())
      {
        ++stack_index;
      }
      // overwrite available slot, otherwise, push back
      if (stack_index != non_delay_local_node_stack.size())
      {
        non_delay_local_node_stack[stack_index] = new_node;
      }
      else
      {
        non_delay_local_node_stack.push_back(new_node);
      }
      // record the stack index
      node_tracker_.emplace(
        std::make_pair(new_node.upwind_node, AAHD_NodeIndex::CreateIndex(stack_index, true)));
      node_tracker_.emplace(
        std::make_pair(new_node.downwind_node, AAHD_NodeIndex::CreateIndex(stack_index, false)));
      ++stack_index;
    }
    // reset nodes from the node stack if the upwind cell is in the current level
    for (AAHD_DirectedEdgeNode& den : non_delay_local_node_stack)
    {
      if (den.IsComingFrom(sorted_level))
        den = AAHD_DirectedEdgeNode();
    }
  }
  local_node_stack_size_ = non_delay_local_node_stack.size();
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForDelayedLocalFaces(const SpatialDiscretization& sdm)
{
  // get reference to the grid
  const MeshContinuum& grid = *(spds_.GetGrid());
  // record the FAS edges
  std::uint64_t fas_node_index = 0;
  const std::vector<std::pair<int, int>>& fas_edges = spds_.GetLocalSweepFAS();
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
        if (neighbor_local_idx == edge.second)
        {
          // record the address of all the nodes
          std::uint32_t num_face_nodes = sdm.GetCellMapping(upwind_cell).GetNumFaceNodes(f);
          for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
          {
            AAHD_FaceNode upwind(edge.first, f, fnode);
            node_tracker_.emplace(
              std::make_pair(upwind, AAHD_NodeIndex::CreateIndex(fas_node_index, true, true)));
            AAHD_FaceNode downwind(edge.second,
                                   face_nodal_mapping.associated_face_,
                                   face_nodal_mapping.face_node_mapping_.at(fnode));
            node_tracker_.emplace(downwind,
                                  AAHD_NodeIndex::CreateIndex(fas_node_index, false, true));
            fas_node_index++;
          }
        }
      }
    }
  }
  delayed_local_node_stack_size_ = fas_node_index;
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForBoundaryAndNonLocalFaces(
  const SpatialDiscretization& sdm)
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
          node_tracker_.emplace(std::make_pair(
            AAHD_FaceNode(cell.local_id, f, fnode),
            AAHD_NodeIndex::CreateBoundaryIndex(incremental_boundary_index,
                                                orientation == FaceOrientation::OUTGOING)));
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
  for (std::uint32_t loc_idx = 0; loc_idx < incoming_bank.size(); ++loc_idx)
  {
    std::uint64_t node_index = 0;
    for (const AAHD_NonLocalFaceNode& nl_en : incoming_bank[loc_idx])
    {
      node_tracker_.emplace(
        std::make_pair(nl_en.node, AAHD_NodeIndex::CreateIndex(node_index, false, false, loc_idx)));
      node_index++;
    }
  }
  for (std::uint32_t loc_idx = 0; loc_idx < delayed_incoming_bank.size(); ++loc_idx)
  {
    std::uint64_t node_index = 0;
    for (const AAHD_NonLocalFaceNode& nl_en : delayed_incoming_bank[loc_idx])
    {
      node_tracker_.emplace(
        std::make_pair(nl_en.node, AAHD_NodeIndex::CreateIndex(node_index, false, true, loc_idx)));
      node_index++;
    }
  }
  for (std::uint32_t loc_idx = 0; loc_idx < outgoing_bank.size(); ++loc_idx)
  {
    std::uint64_t node_index = 0;
    for (const AAHD_NonLocalFaceNode& nl_en : outgoing_bank[loc_idx])
    {
      node_tracker_.emplace(
        std::make_pair(nl_en.node, AAHD_NodeIndex::CreateIndex(node_index, true, false, loc_idx)));
      node_index++;
    }
  }
  // store size
  boundary_node_size_ = incremental_boundary_index;
  for (const std::set<AAHD_NonLocalFaceNode>& nl_en_set : incoming_bank)
  {
    nonlocal_incoming_node_sizes_.push_back(nl_en_set.size());
  }
  for (const std::set<AAHD_NonLocalFaceNode>& nl_en_set : delayed_incoming_bank)
  {
    nonlocal_delayed_incoming_node_sizes_.push_back(nl_en_set.size());
  }
  for (const std::set<AAHD_NonLocalFaceNode>& nl_en_set : outgoing_bank)
  {
    nonlocal_outgoing_node_sizes_.push_back(nl_en_set.size());
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
        node_tracker_.emplace(
          std::make_pair(AAHD_FaceNode(cell.local_id, f, fnode), AAHD_NodeIndex()));
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
