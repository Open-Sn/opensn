#include "volmesher_extruder.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_mpi.h"

//###################################################################
/**Extrude template cells into polygons.*/
void
chi_mesh::VolumeMesherExtruder::ExtrudeCells(chi_mesh::MeshContinuum& template_grid,
                                             chi_mesh::MeshContinuum& grid)
{
  const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);
  //================================================== Start extrusion
  size_t num_global_cells = 0;
  for (size_t iz = 0; iz < (vertex_layers_.size() - 1); iz++)
  {
    for (const auto& template_cell : template_grid.local_cells)
    {
      //========================================= Check template cell type
      if (template_cell.Type() != chi_mesh::CellType::POLYGON)
        throw std::logic_error("Extruder::ExtrudeCells: "
                               "Template cell error. Not of base type POLYGON");

      //========================================= Check cell not inverted
      {
        const auto& v0 = template_cell.centroid_;
        const auto& v1 = template_grid.vertices[template_cell.vertex_ids_[0]];
        const auto& v2 = template_grid.vertices[template_cell.vertex_ids_[1]];

        auto v01 = v1 - v0;
        auto v02 = v2 - v0;

        if (v01.Cross(v02).Dot(khat) < 0.0)
          throw std::logic_error("Extruder attempting to extrude a template"
                                 " cell with a normal pointing downward. This"
                                 " causes erratic behavior and needs to be"
                                 " corrected.");
      }

      auto projected_centroid = ProjectCentroidToLevel(template_cell.centroid_, iz);
      int pid = GetCellKBAPartitionIDFromCentroid(projected_centroid);

      bool has_local_scope = HasLocalScope(template_cell, template_grid, iz);

      if (has_local_scope)
      {
        auto cell = MakeExtrudedCell(
          template_cell, grid, iz, num_global_cells, pid, template_grid.local_cells.size());

        cell->material_id_ = template_cell.material_id_;

        grid.cells.push_back(std::move(cell));
      }
      ++num_global_cells;

    } // for template cell

  } // for iz
}
