// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"

namespace opensn
{

CBCD_FLUDSCommonData::CBCD_FLUDSCommonData(
  const SPDS& spds,
  const std::vector<CellFaceNodalMapping>& grid_nodal_mappings,
  const SpatialDiscretization& sdm)
  : FLUDSCommonData(spds, grid_nodal_mappings),
    sdm_(sdm),
    num_incoming_boundary_nodes_(0),
    num_outgoing_boundary_nodes_(0),
    num_incoming_nonlocal_nodes_(0),
    num_outgoing_nonlocal_nodes_(0),
    device_cell_face_node_map_(nullptr),
    host_cell_face_node_map_()
{
  ComputeCellFaceNodeMap();
  CopyCellFaceNodeMapToDevice();
}

CBCD_FLUDSCommonData::~CBCD_FLUDSCommonData()
{
  DeallocateDeviceCellFaceNodeMap();
}

void
CBCD_FLUDSCommonData::ComputeCellFaceNodeMap()
{
  const MeshContinuum& grid = *(spds_.GetGrid());

  // Pre-compute spatial DOF offset for each cell
  std::vector<size_t> cell_spatial_dof_offsets(grid.local_cells.size());
  size_t current_offset = 0;
  for (const auto& cell : grid.local_cells)
  {
    cell_spatial_dof_offsets[cell.local_id] = current_offset;
    current_offset += sdm_.GetCellMapping(cell).GetNumNodes();
  }

  for (const auto& cell : grid.local_cells)
  {
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
    {
      const CellFace& face = cell.faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell.local_id][f];
      const FaceNodalMapping& face_nodal_mapping = grid_nodal_mappings_[cell.local_id][f];
      std::uint32_t num_face_nodes = sdm_.GetCellMapping(cell).GetNumFaceNodes(f);

      const bool is_outgoing_face = (orientation == FaceOrientation::OUTGOING);
      const bool is_incoming_face = (orientation == FaceOrientation::INCOMING);
      const bool is_local_face = face.IsNeighborLocal(&grid);
      const bool is_boundary_face = not face.has_neighbor;

      for (std::uint32_t fn = 0; fn < num_face_nodes; ++fn)
      {
        CBCD_FaceNode face_node(cell.local_id, f, fn);

        if (is_incoming_face)
        {
          if (is_local_face)
          {
            std::uint32_t nbr_local_idx = face.GetNeighborLocalID(&grid);
            std::uint32_t adj_cell_node = face_nodal_mapping.cell_node_mapping_[fn];
            const std::uint64_t index = cell_spatial_dof_offsets[nbr_local_idx] + adj_cell_node;
            node_tracker_.emplace(face_node,
                                  CBCD_NodeIndex(index, is_outgoing_face, is_local_face, false));
          }
          else if (not is_boundary_face)
          {
            node_tracker_.emplace(
              face_node,
              CBCD_NodeIndex(num_incoming_nonlocal_nodes_, is_outgoing_face, is_local_face, false));
            ++num_incoming_nonlocal_nodes_;
          }
          else
          {
            node_tracker_.emplace(face_node,
                                  CBCD_NodeIndex(num_incoming_boundary_nodes_, is_outgoing_face));
            ++num_incoming_boundary_nodes_;
          }
        }
        else if (is_outgoing_face)
        {
          if (is_local_face)
          {
            const int cell_node = sdm_.GetCellMapping(cell).MapFaceNode(f, fn);
            const std::uint64_t index = cell_spatial_dof_offsets[cell.local_id] + cell_node;
            node_tracker_.emplace(face_node,
                                  CBCD_NodeIndex(index, is_outgoing_face, is_local_face, false));
          }
          else if (not is_boundary_face)
          {
            node_tracker_.emplace(
              face_node,
              CBCD_NodeIndex(num_outgoing_nonlocal_nodes_, is_outgoing_face, is_local_face, false));
            ++num_outgoing_nonlocal_nodes_;
          }
          else
          {
            node_tracker_.emplace(face_node,
                                  CBCD_NodeIndex(num_outgoing_boundary_nodes_, is_outgoing_face));
            ++num_outgoing_boundary_nodes_;
          }
        }
        else
        {
          node_tracker_.emplace(face_node, CBCD_NodeIndex());
        }
      }
    }
  }
}

#ifndef __OPENSN_USE_CUDA__
void
CBCD_FLUDSCommonData::CopyCellFaceNodeMapToDevice()
{
}

void
CBCD_FLUDSCommonData::DeallocateDeviceCellFaceNodeMap()
{
}
#endif

} // namespace opensn