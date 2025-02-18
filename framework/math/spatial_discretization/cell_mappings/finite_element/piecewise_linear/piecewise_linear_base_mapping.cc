// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/cell_mappings/finite_element/piecewise_linear/piecewise_linear_base_mapping.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

PieceWiseLinearBaseMapping::PieceWiseLinearBaseMapping(
  const std::shared_ptr<MeshContinuum> grid,
  const Cell& cell,
  size_t num_nodes,
  std::vector<std::vector<int>> face_node_mappings)
  : CellMapping(grid,
                cell,
                num_nodes,
                GetVertexLocations(grid, cell),
                std::move(face_node_mappings),
                &CellMapping::ComputeCellVolumeAndAreas)
{
}

std::vector<std::vector<int>>
PieceWiseLinearBaseMapping::MakeFaceNodeMapping(const Cell& cell)
{
  const size_t num_faces = cell.faces.size();
  std::vector<std::vector<int>> mappings;
  mappings.reserve(num_faces);
  for (auto& face : cell.faces)
  {
    std::vector<int> face_dof_mapping;
    face_dof_mapping.reserve(face.vertex_ids.size());
    for (uint64_t fvid : face.vertex_ids)
    {
      int mapping = -1;
      for (size_t ci = 0; ci < cell.vertex_ids.size(); ++ci)
      {
        if (fvid == cell.vertex_ids[ci])
        {
          mapping = static_cast<int>(ci);
          break;
        }
      } // for cell i
      if (mapping < 0)
      {
        log.LogAllError() << "Unknown face mapping encountered. "
                             "pwl_polyhedron.h";
        Exit(EXIT_FAILURE);
      }
      face_dof_mapping.push_back(mapping);
    } // for face i

    mappings.push_back(face_dof_mapping);
  }
  return mappings;
}

std::vector<Vector3>
PieceWiseLinearBaseMapping::GetVertexLocations(const std::shared_ptr<MeshContinuum> grid,
                                               const Cell& cell)
{
  std::vector<Vector3> verts;
  verts.reserve(cell.vertex_ids.size());

  for (const auto vid : cell.vertex_ids)
    verts.push_back(grid->vertices[vid]);

  return verts;
}

} // namespace opensn
