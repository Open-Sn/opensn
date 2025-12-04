// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caribou/caribou.h"
#include <cinttypes>

namespace crb = caribou;

namespace opensn
{

void
CBCD_FLUDSCommonData::CopyCellFaceNodeMapToDevice()
{
  crb::HostVector<std::uint64_t> cell_offsets;
  crb::HostVector<std::uint64_t> cell_face_node_index_core_values;

  const MeshContinuum& grid = *(spds_.GetGrid());

  std::uint64_t cell_offset = 2 * grid.local_cells.size();

  for (const auto& cell : grid.local_cells)
  {
    cell_offsets.push_back(cell_offset);

    std::uint64_t num_nodes = 0;
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
    {
      std::uint32_t num_face_nodes = sdm_.GetCellMapping(cell).GetNumFaceNodes(f);
      for (std::uint32_t fn = 0; fn < num_face_nodes; ++fn)
      {
        CBCD_FaceNode face_node(cell.local_id, f, fn);
        CBCD_NodeIndex index = node_tracker_.at(face_node);
        cell_face_node_index_core_values.push_back(index.GetCoreValue());
      }
      num_nodes += num_face_nodes;
      cell_offset += num_face_nodes;
    }

    cell_offsets.push_back(num_nodes);
  }

  crb::DeviceMemory<std::uint64_t> device_cell_face_node_map(
    cell_offsets.size() + cell_face_node_index_core_values.size());
  crb::copy(device_cell_face_node_map, cell_offsets, cell_offsets.size());
  crb::copy(device_cell_face_node_map,
            cell_face_node_index_core_values,
            cell_face_node_index_core_values.size(),
            0,
            cell_offsets.size());
  device_cell_face_node_map_ = device_cell_face_node_map.release();

  // Make a copy on host for performing asynchronous transfers of boundary and non-local angular
  // flux data
  host_cell_face_node_map_.resize(cell_offsets.size() + cell_face_node_index_core_values.size());
  std::copy(cell_offsets.begin(), cell_offsets.end(), host_cell_face_node_map_.begin());
  std::copy(cell_face_node_index_core_values.begin(),
            cell_face_node_index_core_values.end(),
            host_cell_face_node_map_.begin() + cell_offsets.size());
}

void
CBCD_FLUDSCommonData::DeallocateDeviceCellFaceNodeMap()
{
  if (device_cell_face_node_map_ != nullptr)
  {
    crb::DeviceMemory<std::uint64_t> device_cell_face_node_map(device_cell_face_node_map_);
    device_cell_face_node_map.release();
    device_cell_face_node_map_ = nullptr;
  }
}
} // namespace opensn