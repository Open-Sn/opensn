#include "mesh/VolumeMesher/chi_volumemesher.h"
#include "mesh/UnpartitionedMesh/chi_unpartitioned_mesh.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Creates 2D polygon cells for each face of an unpartitioned mesh.*/
void
chi_mesh::VolumeMesher::CreatePolygonCells(const chi_mesh::UnpartitionedMesh& umesh,
                                           chi_mesh::MeshContinuumPtr& grid)
{
  //=================================== Copy nodes
  {
    uint64_t id = 0;
    for (const auto& vertex : umesh.GetVertices())
      grid->vertices.Insert(id++, vertex);
  }

  size_t num_cells = 0;
  for (auto& raw_cell : umesh.GetRawCells())
  {
    // Check valid template cell
    if (raw_cell->type != chi_mesh::CellType::POLYGON)
    {
      Chi::log.LogAllError() << "chi_mesh::VolumeMesher::CreatePolygonCells "
                                "called with a cell not being of primary type"
                                " chi_mesh::CellType::POLYGON.";
      Chi::Exit(EXIT_FAILURE);
    }

    //====================================== Make cell
    auto cell = std::make_unique<chi_mesh::Cell>(CellType::POLYGON, raw_cell->sub_type);

    cell->global_id_ = num_cells;
    cell->local_id_ = num_cells;
    cell->partition_id_ = Chi::mpi.location_id;

    cell->centroid_ = raw_cell->centroid;
    cell->material_id_ = raw_cell->material_id;
    cell->vertex_ids_ = raw_cell->vertex_ids;

    // Copy faces + compute face centroid and normal
    const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);
    for (auto& raw_face : raw_cell->faces)
    {
      chi_mesh::CellFace new_face;

      new_face.vertex_ids_ = raw_face.vertex_ids;

      const auto& v0 = grid->vertices[new_face.vertex_ids_[0]];
      const auto& v1 = grid->vertices[new_face.vertex_ids_[1]];
      new_face.centroid_ = v0 * 0.5 + v1 * 0.5;

      chi_mesh::Vector3 va = v1 - v0;
      chi_mesh::Vector3 vn = va.Cross(khat);
      vn = vn / vn.Norm();
      new_face.normal_ = vn;

      new_face.has_neighbor_ = raw_face.has_neighbor;
      new_face.neighbor_id_ = raw_face.neighbor;

      cell->faces_.push_back(new_face);
    }

    //====================================== Push to grid
    grid->cells.push_back(std::move(cell));
    ++num_cells;
  } // for raw_cell
}