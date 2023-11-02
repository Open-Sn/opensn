#include "chi_grid_vtk_utils.h"

#include "chi_meshcontinuum.h"

#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedIntArray.h>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include "chi_runtime.h"
#include "chi_log.h"

// ###################################################################
/**Uploads vertices and cells to an unstructured grid. This routine
 * also uploads cell material ids (sub-domain ids) and partition ids.*/
vtkNew<vtkUnstructuredGrid>
chi_mesh::PrepareVtkUnstructuredGrid(const chi_mesh::MeshContinuum& grid,
                                     bool discontinuous /*=true*/)
{
  //============================================= Instantiate VTK items
  vtkNew<vtkUnstructuredGrid> ugrid;
  vtkNew<vtkPoints> points;
  vtkNew<vtkIntArray> material_array;
  vtkNew<vtkUnsignedIntArray> partition_id_array;

  points->SetDataType(VTK_DOUBLE);

  //============================================= Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");

  std::vector<uint64_t> vertex_map;
  if (not discontinuous)
  {
    vertex_map.assign(grid.GetGlobalVertexCount(), 0);
    const size_t num_verts = grid.GetGlobalVertexCount();
    for (size_t v = 0; v < num_verts; ++v)
      vertex_map[v] = v;
  }

  // ############################################# Populate cell information
  int64_t node_count = 0;
  for (const auto& cell : grid.local_cells)
  {
    if (discontinuous)
      chi_mesh::UploadCellGeometryDiscontinuous(grid, cell, node_count, points, ugrid);
    else
    {
      for (uint64_t vid : cell.vertex_ids_)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x, vertex.y, vertex.z);
        vertex_map[vid] = node_count;
        ++node_count;
      }
      chi_mesh::UploadCellGeometryContinuous(cell, vertex_map, ugrid);
    }

    material_array->InsertNextValue(cell.material_id_);
    partition_id_array->InsertNextValue(cell.partition_id_);
  } // for local cells
  ugrid->SetPoints(points);

  ugrid->GetCellData()->AddArray(material_array);
  ugrid->GetCellData()->AddArray(partition_id_array);

  return ugrid;
}

// ###################################################################
/**Writes an unstructured grid to files (.pvtu and .vtu).*/
void
chi_mesh::WritePVTUFiles(vtkNew<vtkUnstructuredGrid>& ugrid, const std::string& file_base_name)
{
  //============================================= Construct file name
  std::string base_filename = std::string(file_base_name);
  std::string location_filename =
    base_filename + std::string("_") + std::to_string(Chi::mpi.location_id) + std::string(".vtu");

  //============================================= Write master file
  if (Chi::mpi.location_id == 0)
  {
    std::string pvtu_file_name = base_filename + std::string(".pvtu");

    auto pgrid_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

    pgrid_writer->EncodeAppendedDataOff();
    pgrid_writer->SetFileName(pvtu_file_name.c_str());
    pgrid_writer->SetNumberOfPieces(Chi::mpi.process_count);
    pgrid_writer->SetStartPiece(Chi::mpi.location_id);
    pgrid_writer->SetEndPiece(Chi::mpi.process_count - 1);
    pgrid_writer->SetInputData(ugrid);

    pgrid_writer->Write();
  }
  Chi::mpi.Barrier();

  //============================================= Serial output each piece
  auto grid_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();
}