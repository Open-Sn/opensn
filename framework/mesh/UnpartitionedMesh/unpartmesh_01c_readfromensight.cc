#include "chi_unpartitioned_mesh.h"

#include "mesh/MeshContinuum/chi_grid_vtk_utils.h"

#include "utils/chi_utils.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkEnSightGoldBinaryReader.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkInformation.h>

#include "chi_runtime.h"
#include "chi_log.h"

#define ErrorReadingFile(fname) \
std::runtime_error("Failed to open file: " + options.file_name + \
" in call to " + #fname + ".")

//###################################################################
/**Reads an Ensight-Gold unstructured mesh.*/
void chi_mesh::UnpartitionedMesh::
  ReadFromEnsightGold(const chi_mesh::UnpartitionedMesh::Options &options)
{
  Chi::log.Log() << "Reading Ensight-Gold file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromEnsightGold);
  file.close();

  //======================================== Read the file
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkEnSightGoldBinaryReader>::New();
  reader->SetCaseFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");
  reader->UpdateInformation();
  reader->Update();

  //======================================== Get all the grid blocks
  auto multiblock = reader->GetOutput();

  std::vector<vtkUGridPtrAndName> grid_blocks;
  auto iter_a = multiblock->NewIterator();
  iter_a->GoToFirstItem();
  while (not iter_a->IsDoneWithTraversal())
  {
    auto block_a = iter_a->GetCurrentDataObject();

    const std::string block_name = chi::StringTrim(
      iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()));

    if (block_a->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    {
      grid_blocks.emplace_back(
        vtkUnstructuredGrid::SafeDownCast(block_a),
        chi::StringTrim(
          iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME())));

      Chi::log.Log()
        << "Reading block " << block_name
        << " Number of cells: " << grid_blocks.back().first->GetNumberOfCells()
        << " Number of points: " << grid_blocks.back().first->GetNumberOfPoints();
    }

    iter_a->GoToNextItem();
  }

  //======================================== Get the main + bndry blocks
  const int max_dimension = chi_mesh::FindHighestDimension(grid_blocks);
  Chi::log.Log0Verbose1() << "Maximum dimension : " << max_dimension << "\n";
  std::vector<vtkUGridPtrAndName> domain_grid_blocks =
    chi_mesh::GetBlocksOfDesiredDimension(grid_blocks, max_dimension);
  std::vector<vtkUGridPtrAndName> bndry_grid_blocks =
    chi_mesh::GetBlocksOfDesiredDimension(grid_blocks, max_dimension-1);



  //======================================== Process blocks
  chi_mesh::SetBlockIDArrays(domain_grid_blocks);
  auto ugrid = chi_mesh::ConsolidateGridBlocks(domain_grid_blocks);

  //======================================== Copy Data
  // Material-IDs will get set form block-id arrays
  CopyUGridCellsAndPoints(*ugrid, options.scale, max_dimension);

  //======================================== Always do this
  chi_mesh::MeshAttributes dimension = NONE;
  switch (max_dimension)
  {
    case 1: dimension = DIMENSION_1; break;
    case 2: dimension = DIMENSION_2; break;
    case 3: dimension = DIMENSION_3; break;
    default: break;
  }

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  //======================================== Set boundary ids
  SetBoundaryIDsFromBlocks(bndry_grid_blocks);

  Chi::log.Log() << "Done reading Ensight-Gold file: "
                 << options.file_name << ".";
}