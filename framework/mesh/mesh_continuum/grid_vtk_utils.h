// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <memory>

template <class T>
class vtkNew;
class vtkPoints;
class vtkUnstructuredGrid;
template <class T>
class vtkSmartPointer;

namespace opensn
{
class MeshContinuum;
class Cell;
class CellFace;

/**
 * Uploads vertices and cells to an unstructured grid.
 */
void UploadCellGeometryDiscontinuous(const std::shared_ptr<MeshContinuum> grid,
                                     const Cell& cell,
                                     int64_t& node_counter,
                                     vtkNew<vtkPoints>& points,
                                     vtkNew<vtkUnstructuredGrid>& ugrid);

/**
 * Uploads vertices and cells to an unstructured grid.
 */
void UploadCellGeometryContinuous(const Cell& cell,
                                  const std::vector<uint64_t>& vertex_map,
                                  vtkNew<vtkUnstructuredGrid>& ugrid);
/**
 * Uploads vertices and cells to an unstructured grid.
 */
void UploadFaceGeometry(const CellFace& cell_face,
                        const std::vector<uint64_t>& vertex_map,
                        vtkNew<vtkUnstructuredGrid>& ugrid);

using vtkUGridPtr = vtkSmartPointer<vtkUnstructuredGrid>;
using vtkUGridPtrAndName = std::pair<vtkUGridPtr, std::string>;

/**
 * Finds the highest dimension across all the grid blocks. This is useful when a vtk-read mesh
 * contains multiple blocks. Some of which are boundary faces.
 */
int FindHighestDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks);

/**
 * Consolidates all blocks containing cells with the desired dimension.
 * Thereafter it removes duplicate vertices.
 */
vtkUGridPtr ConsolidateGridBlocks(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                                  const std::string& block_id_array_name = "BlockID");

/**
 * Provides a map of the different grids that have the requested dimension.
 */
std::vector<vtkUGridPtrAndName>
GetBlocksOfDesiredDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks, int desired_dimension);

/**
 * Given several unstructured grid blocks, each denoting a material id, this function sets material
 * ids accordingly.
 */
std::vector<uint64_t> BuildBlockCellExtents(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                                            int desired_dimension);

/**
 * Given several unstructured grid blocks, each denoting a material id, this function creates a VTK
 * cell-data array called "BlockID" that holds this information.
 */
void SetBlockIDArrays(std::vector<vtkUGridPtrAndName>& ugrid_blocks);

/**
 * Retrieves material-ids from a field.
 */
std::vector<int> BuildCellMaterialIDsFromField(vtkUGridPtr& ugrid,
                                               const std::string& field_name,
                                               const std::string& file_name);

/**
 * Uploads vertices and cells to an unstructured grid. This routine also uploads cell material ids
 * (sub-domain ids) and partition ids.
 */
vtkNew<vtkUnstructuredGrid> PrepareVtkUnstructuredGrid(const std::shared_ptr<MeshContinuum> grid,
                                                       bool discontinuous = true);

/**
 * Writes an unstructured grid to files (.pvtu and .vtu).
 */
void WritePVTUFiles(vtkNew<vtkUnstructuredGrid>& ugrid, const std::string& file_base_name);

} // namespace opensn
