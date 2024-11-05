// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/grid_vtk_utils.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedIntArray.h>
#include <vtkAppendFilter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

namespace opensn
{

void
UploadCellGeometryDiscontinuous(const MeshContinuum& grid,
                                const Cell& cell,
                                int64_t& node_counter,
                                vtkNew<vtkPoints>& points,
                                vtkNew<vtkUnstructuredGrid>& ugrid)
{
  size_t num_verts = cell.vertex_ids.size();

  std::vector<vtkIdType> cell_vids(num_verts);
  for (size_t v = 0; v < num_verts; ++v)
  {
    uint64_t vgi = cell.vertex_ids[v];
    std::vector<double> d_node(3);
    d_node[0] = grid.vertices[vgi].x;
    d_node[1] = grid.vertices[vgi].y;
    d_node[2] = grid.vertices[vgi].z;

    points->InsertPoint(node_counter, d_node.data());
    cell_vids[v] = node_counter++;
  }

  if (cell.Type() == CellType::SLAB)
  {
    ugrid->InsertNextCell(VTK_LINE, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == CellType::POLYGON)
  {
    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYGON:
        vtk_subtype = VTK_POLYGON;
        break;
      case CellType::QUADRILATERAL:
        vtk_subtype = VTK_QUAD;
        break;
      case CellType::TRIANGLE:
        vtk_subtype = VTK_TRIANGLE;
        break;
      default:
        vtk_subtype = VTK_POLYGON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == CellType::POLYHEDRON)
  {
    // Build polyhedron faces
    std::vector<vtkIdType> faces_vids;

    size_t num_faces = cell.faces.size();
    for (auto& face : cell.faces)
    {
      size_t num_fverts = face.vertex_ids.size();
      std::vector<vtkIdType> face_info(num_fverts);
      for (size_t fv = 0; fv < num_fverts; ++fv)
      {
        size_t v = 0;
        for (size_t cv = 0; cv < num_verts; ++cv)
          if (cell.vertex_ids[cv] == face.vertex_ids[fv])
          {
            v = cv;
            break;
          }

        face_info[fv] = cell_vids[v];
      }

      faces_vids.push_back(static_cast<vtkIdType>(num_fverts));
      for (auto vid : face_info)
        faces_vids.push_back(vid);
    } // for f

    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYHEDRON:
        vtk_subtype = VTK_POLYHEDRON;
        break;
      case CellType::PYRAMID:
        vtk_subtype = VTK_PYRAMID;
        break;
      case CellType::WEDGE:
        vtk_subtype = VTK_WEDGE;
        break;
      case CellType::HEXAHEDRON:
        vtk_subtype = VTK_HEXAHEDRON;
        break;
      case CellType::TETRAHEDRON:
        vtk_subtype = VTK_TETRA;
        break;
      default:
        vtk_subtype = VTK_POLYHEDRON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype,
                          static_cast<vtkIdType>(num_verts),
                          cell_vids.data(),
                          static_cast<vtkIdType>(num_faces),
                          faces_vids.data());
  } // polyhedron
}

void
UploadCellGeometryContinuous(const Cell& cell,
                             const std::vector<uint64_t>& vertex_map,
                             vtkNew<vtkUnstructuredGrid>& ugrid)
{
  size_t num_verts = cell.vertex_ids.size();

  std::vector<vtkIdType> cell_vids(num_verts);
  for (size_t v = 0; v < num_verts; ++v)
    cell_vids[v] = static_cast<vtkIdType>(vertex_map[cell.vertex_ids[v]]);

  if (cell.Type() == CellType::SLAB)
  {
    ugrid->InsertNextCell(VTK_LINE, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == CellType::POLYGON)
  {
    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYGON:
        vtk_subtype = VTK_POLYGON;
        break;
      case CellType::QUADRILATERAL:
        vtk_subtype = VTK_QUAD;
        break;
      case CellType::TRIANGLE:
        vtk_subtype = VTK_TRIANGLE;
        break;
      default:
        vtk_subtype = VTK_POLYGON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (cell.Type() == CellType::POLYHEDRON)
  {
    // Build polyhedron faces
    std::vector<vtkIdType> faces_vids;

    int vtk_subtype;
    switch (cell.SubType())
    {
      case CellType::POLYHEDRON:
        vtk_subtype = VTK_POLYHEDRON;
        break;
      case CellType::PYRAMID:
        vtk_subtype = VTK_PYRAMID;
        break;
      case CellType::WEDGE:
        vtk_subtype = VTK_WEDGE;
        break;
      case CellType::HEXAHEDRON:
        vtk_subtype = VTK_HEXAHEDRON;
        break;
      case CellType::TETRAHEDRON:
        vtk_subtype = VTK_TETRA;
        break;
      default:
        vtk_subtype = VTK_POLYHEDRON;
        break;
    }

    switch (cell.SubType())
    {
      case CellType::POLYHEDRON:
      {
        size_t num_faces = cell.faces.size();
        for (auto& face : cell.faces)
        {
          size_t num_fverts = face.vertex_ids.size();
          std::vector<vtkIdType> face_info(num_fverts);
          for (size_t fv = 0; fv < num_fverts; ++fv)
          {
            size_t v = 0;
            for (size_t cv = 0; cv < num_verts; ++cv)
              if (cell.vertex_ids[cv] == face.vertex_ids[fv])
              {
                v = cv;
                break;
              }

            face_info[fv] = cell_vids[v];
          }

          faces_vids.push_back(static_cast<vtkIdType>(num_fverts));
          for (auto vid : face_info)
            faces_vids.push_back(vid);
        } // for f

        ugrid->InsertNextCell(vtk_subtype,
                              static_cast<vtkIdType>(num_verts),
                              cell_vids.data(),
                              static_cast<vtkIdType>(num_faces),
                              faces_vids.data());
        break;
      }
      default:
        ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
    }
  } // polyhedron
}

void
UploadFaceGeometry(const CellFace& cell_face,
                   const std::vector<uint64_t>& vertex_map,
                   vtkNew<vtkUnstructuredGrid>& ugrid)
{
  const size_t num_verts = cell_face.vertex_ids.size();

  std::vector<vtkIdType> cell_vids;
  for (uint64_t vid : cell_face.vertex_ids)
    cell_vids.push_back(static_cast<vtkIdType>(vertex_map[vid]));

  if (num_verts == 1)
  {
    ugrid->InsertNextCell(VTK_VERTEX, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (num_verts == 2)
  {
    ugrid->InsertNextCell(VTK_LINE, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
  if (num_verts >= 3)
  {
    int vtk_subtype;
    switch (num_verts)
    {
      case 3:
        vtk_subtype = VTK_TRIANGLE;
        break;
      case 4:
        vtk_subtype = VTK_QUAD;
        break;
      default:
        vtk_subtype = VTK_POLYGON;
        break;
    }

    ugrid->InsertNextCell(vtk_subtype, static_cast<vtkIdType>(num_verts), cell_vids.data());
  }
}

int
FindHighestDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks)
{
  int max_dim = 0;
  for (auto& ugrid : ugrid_blocks)
  {
    const int64_t num_cells = ugrid.first->GetNumberOfCells();
    for (int64_t c = 0; c < num_cells; ++c)
    {
      auto cell = ugrid.first->GetCell(c);
      OpenSnLogicalErrorIf(not cell, "Failed to obtain VTK-cell pointer");
      max_dim = std::max(max_dim, cell->GetCellDimension());
    }
  } // for ugrid in block

  return max_dim;
}

vtkUGridPtr
ConsolidateGridBlocks(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                      const std::string& block_id_array_name)
{
  const std::string fname = "ConsolidateGridBlocks";

  // Determine if all blocks have global-ids
  bool has_global_ids = true;
  for (auto& ugrid_name : ugrid_blocks)
  {
    auto& ugrid = ugrid_name.first;
    const bool has_cell_gids = ugrid->GetCellData()->GetGlobalIds();
    const bool has_pnts_gids = ugrid->GetPointData()->GetGlobalIds();
    const bool has_block_ids = ugrid->GetCellData()->GetArray(block_id_array_name.c_str());

    if ((not has_cell_gids) or (not has_pnts_gids))
      has_global_ids = false;

    if (not has_block_ids)
      throw std::logic_error(fname + ": Grid block " + ugrid_name.second + " does not have \"" +
                             block_id_array_name + "\" array.");
  } // for grid_name pairs

  if (has_global_ids)
    log.Log() << fname << ": blocks have global-id arrays";

  // Consolidate the blocks
  auto append = vtkSmartPointer<vtkAppendFilter>::New();
  for (auto& ugrid : ugrid_blocks)
    append->AddInputData(ugrid.first);

  append->MergePointsOn();
  append->Update();

  auto consolidated_ugrid =
    vtkSmartPointer<vtkUnstructuredGrid>(vtkUnstructuredGrid::SafeDownCast(append->GetOutput()));

  log.Log0Verbose1() << "Consolidated grid num cells and points: "
                     << consolidated_ugrid->GetNumberOfCells() << " "
                     << consolidated_ugrid->GetNumberOfPoints();

  if (has_global_ids)
  {
    const vtkIdType num_points = consolidated_ugrid->GetNumberOfPoints();
    vtkIdType min_id = num_points;
    vtkIdType max_id = 0;
    for (vtkIdType p = 0; p < num_points; ++p)
    {
      auto point_gids =
        vtkIdTypeArray::SafeDownCast(consolidated_ugrid->GetPointData()->GetGlobalIds());
      auto point_gid = point_gids->GetValue(p);

      min_id = std::min(min_id, point_gid);
      max_id = std::max(max_id, point_gid);
    }

    log.Log() << "Minimum and Maximum node-ids " << min_id << " " << max_id;
  }

  std::map<std::string, size_t> cell_type_count_map;
  const int64_t num_cells = consolidated_ugrid->GetNumberOfCells();
  for (int64_t c = 0; c < num_cells; ++c)
  {
    auto cell = consolidated_ugrid->GetCell(c);
    OpenSnLogicalErrorIf(not cell, "Failed to obtain VTK-cell pointer");

    auto cell_type_name = cell->GetClassName();
    cell_type_count_map[cell_type_name] += 1;
  }

  if (log.Verbosity() >= 1)
  {
    std::stringstream outstr;
    /**Lambda to right pad an entry.*/
    auto RightPad = [](std::string& entry, size_t width)
    {
      const size_t pad_size = width - entry.size();
      entry.append(std::string(pad_size, ' '));
    };

    outstr << "Block statistictics:\n";
    for (const auto& [type_name, count] : cell_type_count_map)
    {
      auto temp_name = type_name;
      RightPad(temp_name, 20);
      outstr << "  " << temp_name << " " << count << "\n";
    }
    log.Log0Verbose1() << outstr.str();
  }

  return consolidated_ugrid;
}

std::vector<vtkUGridPtrAndName>
BlocksOfDesiredDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks, int desired_dimension)
{
  std::vector<vtkUGridPtrAndName> desired_blocks;
  for (auto& ugrid : ugrid_blocks)
  {
    if (ugrid.first->GetNumberOfCells() == 0)
      continue;

    std::vector<vtkUGridPtrAndName> single_grid = {ugrid};
    int block_dimension = FindHighestDimension(single_grid);

    if (block_dimension == desired_dimension)
      desired_blocks.push_back(ugrid);
  }

  return desired_blocks;
}

std::vector<uint64_t>
BuildBlockCellExtents(std::vector<vtkUGridPtrAndName>& ugrid_blocks, const int desired_dimension)
{
  std::vector<uint64_t> block_mat_ids;
  size_t total_cells = 0;

  for (auto& ugrid : ugrid_blocks)
  {
    uint64_t num_cells = ugrid.first->GetNumberOfCells();

    if (num_cells == 0)
      continue;

    if (ugrid.first->GetCell(0)->GetCellDimension() == desired_dimension)
    {
      total_cells += num_cells;
      block_mat_ids.push_back(total_cells);
    }
  }
  return block_mat_ids;
}

void
SetBlockIDArrays(std::vector<vtkUGridPtrAndName>& ugrid_blocks)
{
  int block_id = 0;
  for (auto& ugrid : ugrid_blocks)
  {
    const vtkIdType num_cells = ugrid.first->GetNumberOfCells();

    if (num_cells == 0)
      continue;

    vtkNew<vtkIntArray> block_id_list;
    block_id_list->SetName("BlockID");

    for (vtkIdType c = 0; c < num_cells; ++c)
      block_id_list->InsertNextValue(block_id);

    auto arr = ugrid.first->GetCellData()->GetArray("BlockID");
    if (not arr)
      ugrid.first->GetCellData()->RemoveArray("BlockID");

    ugrid.first->GetCellData()->AddArray(block_id_list);
    ++block_id;
  }
}

std::vector<int>
BuildCellMaterialIDsFromField(vtkUGridPtr& ugrid,
                              const std::string& field_name,
                              const std::string& file_name)
{
  const size_t total_cell_count = ugrid->GetNumberOfCells();
  std::vector<int> material_ids(total_cell_count, -1);

  // Determine if reading cell identifiers
  vtkDataArray* cell_id_array_ptr;
  if (field_name.empty())
  {
    log.Log0Warning() << "A user-supplied field name from which to recover material "
                         "identifiers "
                      << "has not been found. Material-ids will be left unassigned.";
    goto end_error_checks;
  }
  else
  {
    auto cell_data = ugrid->GetCellData();
    const auto vtk_abstract_array_ptr = cell_data->GetAbstractArray(field_name.c_str());

    if (not vtk_abstract_array_ptr)
    {
      log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
                        << "does not contain a vtkCellData field of name : \"" << field_name
                        << "\". Material-ids will be left unassigned.";
      goto end_error_checks;
    }

    cell_id_array_ptr = vtkArrayDownCast<vtkDataArray>(vtk_abstract_array_ptr);
    if (not cell_id_array_ptr)
    {
      log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
                        << "with vtkCellData field of name : \"" << field_name << "\" "
                        << "cannot be downcast to vtkDataArray. Material-ids will be left "
                           "unassigned.";
      goto end_error_checks;
    }

    const auto cell_id_n_tup = cell_id_array_ptr->GetNumberOfTuples();
    if (cell_id_n_tup != total_cell_count)
    {
      log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
                        << "with vtkCellData field of name : \"" << field_name
                        << "\" has n. tuples : " << cell_id_n_tup
                        << ", but differs from the value expected : " << total_cell_count
                        << ". Material-ids will be left unassigned.";
      goto end_error_checks;
    }

    const auto cell_id_n_val = cell_id_array_ptr->GetNumberOfValues();
    if (cell_id_n_val != total_cell_count)
    {
      log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
                        << "with vtkCellData field of name : \"" << field_name
                        << "\" has n. values : " << cell_id_n_val
                        << ", but differs from the value expected : " << total_cell_count
                        << ". Material-ids will be left unassigned.";
      goto end_error_checks;
    }
  }

  //  apply cell identifier
  for (size_t c = 0; c < total_cell_count; ++c)
  {
    std::vector<double> cell_id_vec(1);
    cell_id_array_ptr->GetTuple(static_cast<vtkIdType>(c), cell_id_vec.data());
    const auto mat_id = (int)cell_id_vec.front();

    material_ids[c] = mat_id;
  }

end_error_checks:
  return material_ids;
}

vtkNew<vtkUnstructuredGrid>
PrepareVtkUnstructuredGrid(const MeshContinuum& grid, bool discontinuous)
{
  // Instantiate VTK items
  vtkNew<vtkUnstructuredGrid> ugrid;
  vtkNew<vtkPoints> points;
  vtkNew<vtkIntArray> material_array;
  vtkNew<vtkUnsignedIntArray> partition_id_array;

  points->SetDataType(VTK_DOUBLE);

  // Set names
  material_array->SetName("Material");
  partition_id_array->SetName("Partition");

  std::vector<uint64_t> vertex_map;
  if (not discontinuous)
  {
    vertex_map.assign(grid.GlobalVertexCount(), 0);
    const size_t num_verts = grid.GlobalVertexCount();
    for (size_t v = 0; v < num_verts; ++v)
      vertex_map[v] = v;
  }

  // Populate cell information
  int64_t node_count = 0;
  for (const auto& cell : grid.local_cells)
  {
    if (discontinuous)
      UploadCellGeometryDiscontinuous(grid, cell, node_count, points, ugrid);
    else
    {
      for (uint64_t vid : cell.vertex_ids)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x, vertex.y, vertex.z);
        vertex_map[vid] = node_count;
        ++node_count;
      }
      UploadCellGeometryContinuous(cell, vertex_map, ugrid);
    }

    material_array->InsertNextValue(cell.material_id);
    partition_id_array->InsertNextValue(cell.partition_id);
  } // for local cells
  ugrid->SetPoints(points);

  ugrid->GetCellData()->AddArray(material_array);
  ugrid->GetCellData()->AddArray(partition_id_array);

  return ugrid;
}

void
WritePVTUFiles(vtkNew<vtkUnstructuredGrid>& ugrid, const std::string& file_base_name)
{
  // Construct file name
  std::string base_filename = std::string(file_base_name);
  std::string location_filename = base_filename + std::string("_") +
                                  std::to_string(opensn::mpi_comm.rank()) + std::string(".vtu");

  // Write master file
  if (opensn::mpi_comm.rank() == 0)
  {
    std::string pvtu_file_name = base_filename + std::string(".pvtu");

    auto pgrid_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

    pgrid_writer->EncodeAppendedDataOff();
    pgrid_writer->SetFileName(pvtu_file_name.c_str());
    pgrid_writer->SetNumberOfPieces(opensn::mpi_comm.size());
    pgrid_writer->SetStartPiece(opensn::mpi_comm.rank());
    pgrid_writer->SetEndPiece(opensn::mpi_comm.size() - 1);
    pgrid_writer->SetInputData(ugrid);

    pgrid_writer->Write();
  }
  opensn::mpi_comm.barrier();

  // Serial output each piece
  auto grid_writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();
}

} // namespace opensn
