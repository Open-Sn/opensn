// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caribou/main.hpp"
#include <cinttypes>

namespace crb = caribou;

namespace opensn
{

void
CBCD_FLUDSCommonData::CopyFlattenedNodeIndexToDevice(const SpatialDiscretization& sdm)
{
  const MeshContinuum& grid = *(spds_.GetGrid());
  const size_t num_local_cells = grid.local_cells.size();
  std::uint64_t total_face_nodes = 0;
  for (const auto& cell : grid.local_cells)
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
      total_face_nodes += sdm.GetCellMapping(cell).GetNumFaceNodes(f);
  std::vector<size_t> cell_spatial_dof_offsets(num_local_cells);
  size_t current_dof_offset = 0;
  for (const auto& cell : grid.local_cells)
  {
    cell_spatial_dof_offsets[cell.local_id] = current_dof_offset;
    current_dof_offset += sdm.GetCellMapping(cell).GetNumNodes();
  }
  const size_t offsets_size = 2 * num_local_cells;
  const size_t total_size = offsets_size + total_face_nodes;
  std::vector<std::uint64_t> local_map(total_size);
  std::uint64_t* cell_offsets_ptr = local_map.data();
  std::uint64_t* indices_ptr = local_map.data() + offsets_size;
  std::uint64_t current_index_offset = offsets_size;
  std::uint64_t local_indices_filled = 0;
  // Iterate over cells to fill the map and populate metadata structures
  for (const auto& cell : grid.local_cells)
  {
    cell_offsets_ptr[2 * cell.local_id] = current_index_offset;
    std::uint64_t num_cell_nodes = 0;
    for (size_t f = 0; f < cell.faces.size(); ++f)
    {
      const CellFace& face = cell.faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell.local_id][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell.local_id][f];
      const size_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
      const bool is_outgoing_face = (orientation == FaceOrientation::OUTGOING);
      const bool is_incoming_face = (orientation == FaceOrientation::INCOMING);
      const bool is_local_face = face.IsNeighborLocal(&grid);
      const bool is_boundary_face = not face.has_neighbor;
      for (size_t fn = 0; fn < num_face_nodes; ++fn)
      {
        CBCD_NodeIndex node_index;

        if (is_incoming_face)
        {
          if (is_local_face)
          {
            std::uint32_t nbr_local_idx = face.GetNeighborLocalID(&grid);
            std::uint32_t adj_cell_node = face_nodal_mapping.cell_node_mapping_[fn];
            const std::uint64_t index = cell_spatial_dof_offsets[nbr_local_idx] + adj_cell_node;
            node_index = CBCD_NodeIndex(index, is_outgoing_face, is_local_face);
          }
          else if (not is_boundary_face)
          {
            node_index =
              CBCD_NodeIndex(num_incoming_nonlocal_nodes_, is_outgoing_face, is_local_face);
            cell_to_incoming_nonlocal_nodes_[cell.local_id].emplace_back(
              NonlocalNodeInfo{cell.local_id,
                               cell.global_id,
                               static_cast<unsigned int>(f),
                               fn,
                               face_nodal_mapping.face_node_mapping_[fn],
                               static_cast<std::uint64_t>(num_incoming_nonlocal_nodes_)});
            ++num_incoming_nonlocal_nodes_;
          }
          else
          {
            node_index = CBCD_NodeIndex(num_incoming_boundary_nodes_, is_outgoing_face);
            incoming_boundary_node_map_.emplace_back(
              BoundaryNodeInfo{cell.local_id,
                               static_cast<unsigned int>(f),
                               fn,
                               static_cast<std::uint64_t>(num_incoming_boundary_nodes_),
                               face.neighbor_id});
            ++num_incoming_boundary_nodes_;
          }
        }
        else if (is_outgoing_face)
        {
          if (is_local_face)
          {
            const int cell_node = sdm.GetCellMapping(cell).MapFaceNode(f, fn);
            const std::uint64_t index = cell_spatial_dof_offsets[cell.local_id] + cell_node;
            node_index = CBCD_NodeIndex(index, is_outgoing_face, is_local_face);
          }
          else if (not is_boundary_face)
          {
            node_index =
              CBCD_NodeIndex(num_outgoing_nonlocal_nodes_, is_outgoing_face, is_local_face);
            cell_to_outgoing_nonlocal_nodes_[cell.local_id].emplace_back(
              NonlocalNodeInfo{cell.local_id,
                               cell.global_id,
                               static_cast<unsigned int>(f),
                               fn,
                               face_nodal_mapping.face_node_mapping_[fn],
                               static_cast<std::uint64_t>(num_outgoing_nonlocal_nodes_)});
            ++num_outgoing_nonlocal_nodes_;
          }
          else
          {
            node_index = CBCD_NodeIndex(num_outgoing_boundary_nodes_, is_outgoing_face);
            cell_to_outgoing_boundary_nodes_[cell.local_id].emplace_back(
              BoundaryNodeInfo{cell.local_id,
                               static_cast<unsigned int>(f),
                               fn,
                               static_cast<std::uint64_t>(num_outgoing_boundary_nodes_),
                               face.neighbor_id});
            ++num_outgoing_boundary_nodes_;
          }
        }
        else
        {
          node_index = CBCD_NodeIndex();
        }
        indices_ptr[local_indices_filled++] = node_index.GetCoreValue();
      }
      num_cell_nodes += num_face_nodes;
    }
    cell_offsets_ptr[2 * cell.local_id + 1] = num_cell_nodes;
    current_index_offset += num_cell_nodes;
  }
  if (local_map.empty())
    return;
  crb::HostVector<std::uint64_t> host_mem(local_map.begin(), local_map.end());
  crb::DeviceMemory<std::uint64_t> device_mem(local_map.size());
  crb::copy(device_mem, host_mem, host_mem.size());
  device_cell_face_node_map_ = device_mem.release();
}

void
CBCD_FLUDSCommonData::DeallocateDeviceMemory()
{
  if (device_cell_face_node_map_ != nullptr)
  {
    crb::DeviceMemory<std::uint64_t> device_cell_face_node_map(device_cell_face_node_map_);
    device_cell_face_node_map.release();
    device_cell_face_node_map_ = nullptr;
  }
}
} // namespace opensn