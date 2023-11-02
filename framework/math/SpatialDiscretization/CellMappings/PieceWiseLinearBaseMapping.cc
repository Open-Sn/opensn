#include "math/SpatialDiscretization/CellMappings/PieceWiseLinearBaseMapping.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace chi_math::cell_mapping
{

PieceWiseLinearBaseMapping::PieceWiseLinearBaseMapping(
  const chi_mesh::MeshContinuum& grid,
  const chi_mesh::Cell& cell,
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

/** This section just determines a mapping of face dofs
to cell dofs. This is pretty simple since we can
just loop over each face dof then subsequently
loop over cell dofs, if the face dof node index equals
the cell dof node index then the mapping is assigned.

This mapping is not used by any of the methods in
    this class but is used by methods requiring the
      surface integrals of the shape functions.*/
std::vector<std::vector<int>>
PieceWiseLinearBaseMapping::MakeFaceNodeMapping(const chi_mesh::Cell& cell)
{
  const size_t num_faces = cell.faces_.size();
  std::vector<std::vector<int>> mappings;
  mappings.reserve(num_faces);
  for (auto& face : cell.faces_)
  {
    std::vector<int> face_dof_mapping;
    face_dof_mapping.reserve(face.vertex_ids_.size());
    for (uint64_t fvid : face.vertex_ids_)
    {
      int mapping = -1;
      for (size_t ci = 0; ci < cell.vertex_ids_.size(); ci++)
      {
        if (fvid == cell.vertex_ids_[ci])
        {
          mapping = static_cast<int>(ci);
          break;
        }
      } // for cell i
      if (mapping < 0)
      {
        Chi::log.LogAllError() << "Unknown face mapping encountered. "
                                  "pwl_polyhedron.h";
        Chi::Exit(EXIT_FAILURE);
      }
      face_dof_mapping.push_back(mapping);
    } // for face i

    mappings.push_back(face_dof_mapping);
  }
  return mappings;
}

std::vector<chi_mesh::Vector3> PieceWiseLinearBaseMapping::GetVertexLocations(
  const chi_mesh::MeshContinuum& grid, const chi_mesh::Cell& cell)
{
  std::vector<chi_mesh::Vector3> verts;
  verts.reserve(cell.vertex_ids_.size());

  for (const auto vid : cell.vertex_ids_)
    verts.push_back(grid.vertices[vid]);

  return verts;
}

} // namespace chi_math::cell_mapping
