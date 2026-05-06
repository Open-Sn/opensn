// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds_common_data.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caribou/main.hpp"
#include <cinttypes>

namespace crb = caribou;

namespace opensn
{

void
AAHD_FLUDSCommonData::UpdateBoundaryAndSyncWithDevice(
  const SpatialDiscretization& sdm,
  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries)
{
  ComputeNodeIndexForBoundaryFaces(sdm, boundaries);
  CopyFlattenNodeIndexToDevice(sdm);
  for (auto angleset : associated_anglesets_)
  {
    angleset->CopyBoundaryOffsetToDevice();
  }
}

void
AAHD_FLUDSCommonData::ComputeNodeIndexForBoundaryFaces(
  const SpatialDiscretization& sdm,
  const std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries)
{
  std::uint64_t incremental_boundary_index = 0;
  const MeshContinuum& grid = *(spds_.GetGrid());
  for (const Cell& cell : grid.local_cells)
  {
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
    {
      const CellFace& face = cell.faces[f];
      const FaceOrientation& orientation = spds_.GetCellFaceOrientations()[cell.local_id][f];
      std::uint32_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);

      if (face.has_neighbor)
        continue;

      auto& boundary = *(boundaries.at(face.neighbor_id));
      if (orientation == FaceOrientation::INCOMING)
      {
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          FaceNode fn(cell.local_id, f, fnode);
          for (auto angleset : associated_anglesets_)
          {
            std::uint64_t offset = boundary.GetOffsetToAngleset(fn, *angleset, false);
            angleset->AddBoundaryOffset(offset);
          }
          node_tracker_.emplace(fn,
                                AAHD_NodeIndex(incremental_boundary_index,
                                               false,
                                               true,
                                               boundary.IsReflecting(),
                                               boundary.IsAngleDependent()));
          ++incremental_boundary_index;
        }
      }
      else if (orientation == FaceOrientation::OUTGOING && boundary.IsReflecting())
      {
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          FaceNode fn(cell.local_id, f, fnode);
          for (auto angleset : associated_anglesets_)
          {
            std::uint64_t offset = boundary.GetOffsetToAngleset(fn, *angleset, true);
            angleset->AddBoundaryOffset(offset);
          }
          node_tracker_.emplace(fn,
                                AAHD_NodeIndex(incremental_boundary_index, true, true, true, true));
          ++incremental_boundary_index;
        }
      }
      else
      {
        for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
        {
          FaceNode fn(cell.local_id, f, fnode);
          node_tracker_.emplace(FaceNode(cell.local_id, f, fnode), AAHD_NodeIndex());
        }
      }
    }
  }
}

void
AAHD_FLUDSCommonData::CopyFlattenNodeIndexToDevice(const SpatialDiscretization& sdm)
{
  // initialize flatten index and cell offset
  crb::HostVector<std::uint64_t> cell_offset;
  crb::HostVector<std::uint64_t> data;
  // loop for each cell
  const MeshContinuum& grid = *(spds_.GetGrid());
  std::uint64_t offset = 2 * grid.local_cells.size();
  for (const Cell& cell : grid.local_cells)
  {
    // record the offset
    cell_offset.push_back(offset);
    // record each face node to the data
    std::uint64_t num_nodes = 0;
    for (std::uint32_t f = 0; f < cell.faces.size(); ++f)
    {
      std::uint32_t num_face_nodes = sdm.GetCellMapping(cell).GetNumFaceNodes(f);
      for (std::uint32_t fnode = 0; fnode < num_face_nodes; ++fnode)
      {
        FaceNode node(cell.local_id, f, fnode);
        AAHD_NodeIndex index = node_tracker_.at(node);
        data.push_back(index.GetCoreValue());
        ++num_nodes;
        ++offset;
      }
    }
    cell_offset.push_back(num_nodes);
  }
  // copy data to GPU
  crb::DeviceMemory<std::uint64_t> device_mem_ptr(cell_offset.size() + data.size());
  crb::copy(device_mem_ptr, cell_offset, cell_offset.size());
  crb::copy(device_mem_ptr, data, data.size(), 0, cell_offset.size());
  device_node_indexes_ = device_mem_ptr.release();
}

void
AAHD_FLUDSCommonData::DeallocateDeviceMemory()
{
  if (device_node_indexes_)
  {
    crb::DeviceMemory<std::uint64_t> device_mem_ptr(device_node_indexes_);
    device_mem_ptr.release();
    device_node_indexes_ = nullptr;
  }
}

} // namespace opensn
