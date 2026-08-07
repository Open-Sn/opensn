// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_base_mapping.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

PieceWiseLinearBaseMapping::PieceWiseLinearBaseMapping(
  const std::shared_ptr<Mesh> grid,
  std::uint32_t cell_local_id,
  size_t num_nodes,
  std::vector<std::vector<int>> face_node_mappings)
  : CellMapping(grid,
                cell_local_id,
                num_nodes,
                GetVertexLocations(grid, cell_local_id),
                std::move(face_node_mappings))
{
}

std::vector<std::vector<int>>
PieceWiseLinearBaseMapping::MakeFaceNodeMapping(const std::shared_ptr<Mesh>& grid,
                                                std::uint32_t cell_local_id)
{
  const auto& cell = grid->GetLocalCell(cell_local_id);
  auto cell_vertex_ids = grid->GetCellConnectivity(cell_local_id);
  const size_t num_faces = cell.faces.size();
  std::vector<std::vector<int>> mappings;
  mappings.reserve(num_faces);
  for (std::uint32_t face_idx = 0; face_idx < cell.faces.size(); ++face_idx)
  {
    auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, face_idx);
    std::vector<int> face_dof_mapping;
    face_dof_mapping.reserve(face_vertex_ids.size());
    for (uint64_t fvid : face_vertex_ids)
    {
      int mapping = -1;
      for (size_t ci = 0; ci < cell_vertex_ids.size(); ++ci)
      {
        if (fvid == cell_vertex_ids[ci])
        {
          mapping = static_cast<int>(ci);
          break;
        }
      } // for cell i
      if (mapping < 0)
        throw std::runtime_error("Unknown face mapping encountered.");

      face_dof_mapping.push_back(mapping);
    } // for face i

    mappings.push_back(face_dof_mapping);
  }
  return mappings;
}

std::vector<Vector3>
PieceWiseLinearBaseMapping::GetVertexLocations(const std::shared_ptr<Mesh>& grid,
                                               std::uint32_t cell_local_id)
{
  auto cell_vertex_ids = grid->GetCellConnectivity(cell_local_id);
  std::vector<Vector3> verts;
  verts.reserve(cell_vertex_ids.size());

  for (const auto vid : cell_vertex_ids)
    verts.push_back(grid->GlobalVertex(vid));

  return verts;
}

std::size_t
PieceWiseLinearBaseMapping::GetNumberOfNodes(const std::shared_ptr<Mesh>& grid,
                                             std::uint32_t cell_local_id)
{
  auto cell_vertex_ids = grid->GetCellConnectivity(cell_local_id);
  return cell_vertex_ids.size();
}

} // namespace opensn
