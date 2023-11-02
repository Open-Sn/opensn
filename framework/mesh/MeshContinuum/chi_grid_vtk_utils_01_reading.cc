#include "chi_grid_vtk_utils.h"

#include <vtkUnstructuredGrid.h>
#include <vtkAppendFilter.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

#include "chi_runtime.h"
#include "chi_log.h"

// ###################################################################
/**Finds the highest dimension across all the grid blocks. This
 * is useful when a vtk-read mesh contains multiple blocks. Some of which
 * are boundary faces.*/
int
chi_mesh::FindHighestDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks)
{
  int max_dim = 0;
  for (auto& ugrid : ugrid_blocks)
  {
    const int64_t num_cells = ugrid.first->GetNumberOfCells();
    for (int64_t c = 0; c < num_cells; ++c)
    {
      auto cell = ugrid.first->GetCell(c);
      ChiLogicalErrorIf(not cell, "Failed to obtain VTK-cell pointer");
      max_dim = std::max(max_dim, cell->GetCellDimension());
    }
  } // for ugrid in block

  return max_dim;
}

// ###################################################################
/**Consolidates all blocks containing cells with the desired dimension.
 * Thereafter it removes duplicate vertices.*/
chi_mesh::vtkUGridPtr
chi_mesh::ConsolidateGridBlocks(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                                const std::string& block_id_array_name /*="BlockID"*/)
{
  const std::string fname = "chi_mesh::ConsolidateGridBlocks";

  //======================================== Determine if all blocks have
  //                                         global-ids
  bool has_global_ids = true;
  for (auto& ugrid_name : ugrid_blocks)
  {
    auto& ugrid = ugrid_name.first;
    const bool has_cell_gids = ugrid->GetCellData()->GetGlobalIds();
    const bool has_pnts_gids = ugrid->GetPointData()->GetGlobalIds();
    const bool has_block_ids = ugrid->GetCellData()->GetArray(block_id_array_name.c_str());

    if ((not has_cell_gids) or (not has_pnts_gids)) has_global_ids = false;

    if (not has_block_ids)
      throw std::logic_error(fname + ": Grid block " + ugrid_name.second + " does not have \"" +
                             block_id_array_name + "\" array.");
  } // for grid_name pairs

  if (has_global_ids) Chi::log.Log() << fname << ": blocks have global-id arrays";

  //======================================== Consolidate the blocks
  auto append = vtkSmartPointer<vtkAppendFilter>::New();
  for (auto& ugrid : ugrid_blocks)
    append->AddInputData(ugrid.first);

  append->MergePointsOn();
  append->Update();

  auto consolidated_ugrid =
    vtkSmartPointer<vtkUnstructuredGrid>(vtkUnstructuredGrid::SafeDownCast(append->GetOutput()));

  Chi::log.Log0Verbose1() << "Consolidated grid num cells and points: "
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

    Chi::log.Log() << "Minimum and Maximum node-ids " << min_id << " " << max_id;
  }

  std::map<std::string, size_t> cell_type_count_map;
  const int64_t num_cells = consolidated_ugrid->GetNumberOfCells();
  for (int64_t c = 0; c < num_cells; ++c)
  {
    auto cell = consolidated_ugrid->GetCell(c);
    ChiLogicalErrorIf(not cell, "Failed to obtain VTK-cell pointer");

    auto cell_type_name = cell->GetClassName();
    cell_type_count_map[cell_type_name] += 1;
  }

  if (Chi::log.GetVerbosity() >= 1)
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
    Chi::log.Log0Verbose1() << outstr.str();
  }

  return consolidated_ugrid;
}

// ###################################################################
/**Provides a map of the different grids that have the
 * requested dimension.*/
std::vector<chi_mesh::vtkUGridPtrAndName>
chi_mesh::GetBlocksOfDesiredDimension(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                                      int desired_dimension)
{
  std::vector<chi_mesh::vtkUGridPtrAndName> desired_blocks;
  for (auto& ugrid : ugrid_blocks)
  {
    if (ugrid.first->GetNumberOfCells() == 0) continue;

    std::vector<vtkUGridPtrAndName> single_grid = {ugrid};
    int block_dimension = chi_mesh::FindHighestDimension(single_grid);

    if (block_dimension == desired_dimension) desired_blocks.push_back(ugrid);
  }

  return desired_blocks;
}

// ###################################################################
/**Given several unstructured grid blocks, each denoting a material id,
 * this function sets material ids accordingly.*/
std::vector<uint64_t>
chi_mesh::BuildBlockCellExtents(std::vector<vtkUGridPtrAndName>& ugrid_blocks,
                                const int desired_dimension)
{
  std::vector<uint64_t> block_mat_ids;
  size_t total_cells = 0;

  for (auto& ugrid : ugrid_blocks)
  {
    uint64_t num_cells = ugrid.first->GetNumberOfCells();

    if (num_cells == 0) continue;

    if (ugrid.first->GetCell(0)->GetCellDimension() == desired_dimension)
    {
      total_cells += num_cells;
      block_mat_ids.push_back(total_cells);
    }
  }
  return block_mat_ids;
}

// ###################################################################
/**Given several unstructured grid blocks, each denoting a material id,
 * this function creates a VTK cell-data array called "BlockID" that holds
 * this information.*/
void
chi_mesh::SetBlockIDArrays(std::vector<vtkUGridPtrAndName>& ugrid_blocks)
{
  int block_id = 0;
  for (auto& ugrid : ugrid_blocks)
  {
    const vtkIdType num_cells = ugrid.first->GetNumberOfCells();

    if (num_cells == 0) continue;

    vtkNew<vtkIntArray> block_id_list;
    block_id_list->SetName("BlockID");

    for (vtkIdType c = 0; c < num_cells; ++c)
      block_id_list->InsertNextValue(block_id);

    auto arr = ugrid.first->GetCellData()->GetArray("BlockID");
    if (not arr) ugrid.first->GetCellData()->RemoveArray("BlockID");

    ugrid.first->GetCellData()->AddArray(block_id_list);
    ++block_id;
  }
}

// ###################################################################
/**Retrieves material-ids from a field.*/
std::vector<int>
chi_mesh::BuildCellMaterialIDsFromField(vtkUGridPtr& ugrid,
                                        const std::string& field_name,
                                        const std::string& file_name)
{
  const size_t total_cell_count = ugrid->GetNumberOfCells();
  std::vector<int> material_ids(total_cell_count, -1);

  //======================================== Determine if reading
  //                                         cell identifiers
  vtkDataArray* cell_id_array_ptr;
  if (field_name.empty())
  {
    Chi::log.Log0Warning() << "A user-supplied field name from which to recover material "
                              "identifiers "
                           << "has not been found. Material-ids will be left unassigned.";
    goto end_error_checks;
  }
  else
  {
    auto cell_data = ugrid->GetCellData();
    const auto vtk_abstract_array_ptr = cell_data->GetAbstractArray(field_name.c_str());

    if (!vtk_abstract_array_ptr)
    {
      Chi::log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
                             << "does not contain a vtkCellData field of name : \"" << field_name
                             << "\". Material-ids will be left unassigned.";
      goto end_error_checks;
    }

    cell_id_array_ptr = vtkArrayDownCast<vtkDataArray>(vtk_abstract_array_ptr);
    if (!cell_id_array_ptr)
    {
      Chi::log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
                             << "with vtkCellData field of name : \"" << field_name << "\" "
                             << "cannot be downcast to vtkDataArray. Material-ids will be left "
                                "unassigned.";
      goto end_error_checks;
    }

    const auto cell_id_n_tup = cell_id_array_ptr->GetNumberOfTuples();
    if (cell_id_n_tup != total_cell_count)
    {
      Chi::log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
                             << "with vtkCellData field of name : \"" << field_name
                             << "\" has n. tuples : " << cell_id_n_tup
                             << ", but differs from the value expected : " << total_cell_count
                             << ". Material-ids will be left unassigned.";
      goto end_error_checks;
    }

    const auto cell_id_n_val = cell_id_array_ptr->GetNumberOfValues();
    if (cell_id_n_val != total_cell_count)
    {
      Chi::log.Log0Warning() << "The VTU file : \"" << file_name << "\" "
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