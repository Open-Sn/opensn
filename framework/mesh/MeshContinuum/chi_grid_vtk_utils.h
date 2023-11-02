#ifndef CHITECH_CHI_GRID_VTK_UTILS_H
#define CHITECH_CHI_GRID_VTK_UTILS_H

#include <cstdint>
#include <string>
#include <vector>
#include <map>

template<class T>
class vtkNew;
class vtkPoints;
class vtkUnstructuredGrid;
template<class T>
class vtkSmartPointer;

namespace chi_mesh
{
  class MeshContinuum;
  class Cell;
  class CellFace;

  //00
  void UploadCellGeometryDiscontinuous(const chi_mesh::MeshContinuum& grid,
                                       const chi_mesh::Cell &cell,
                                       int64_t& node_counter,
                                       vtkNew<vtkPoints>& points,
                                       vtkNew<vtkUnstructuredGrid>& ugrid);
  void UploadCellGeometryContinuous(const chi_mesh::Cell &cell,
                                    const std::vector<uint64_t>& vertex_map,
                                    vtkNew<vtkUnstructuredGrid>& ugrid);
  void UploadFaceGeometry(const chi_mesh::CellFace& cell_face,
                          const std::vector<uint64_t>& vertex_map,
                          vtkNew<vtkUnstructuredGrid> &ugrid);

  //01 Utils for Reading
  typedef vtkSmartPointer<vtkUnstructuredGrid> vtkUGridPtr;
  typedef std::pair<vtkUGridPtr, std::string> vtkUGridPtrAndName;

  int FindHighestDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks);

  vtkUGridPtr
  ConsolidateGridBlocks(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                        const std::string& block_id_array_name = "BlockID");

  std::vector<vtkUGridPtrAndName>
    GetBlocksOfDesiredDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                                int desired_dimension);

  std::vector<uint64_t>
  BuildBlockCellExtents(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                        int desired_dimension);

  void SetBlockIDArrays(std::vector<vtkUGridPtrAndName>& ugrid_blocks);

  std::vector<int>
  BuildCellMaterialIDsFromField(vtkUGridPtr &ugrid,
                                const std::string& field_name,
                                const std::string& file_name);

  //04 Writing VTK files
  vtkNew<vtkUnstructuredGrid>
    PrepareVtkUnstructuredGrid(const chi_mesh::MeshContinuum& grid,
                              bool discontinuous = true);

  void WritePVTUFiles(vtkNew<vtkUnstructuredGrid> &ugrid,
                      const std::string &file_base_name);

}//namespace chi_mesh

#endif //CHITECH_CHI_GRID_VTK_UTILS_H
