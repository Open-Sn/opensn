#include "chi_unpartitioned_mesh.h"

#include "mesh/MeshContinuum/chi_grid_vtk_utils.h"

#include "utils/chi_utils.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkExodusIIReader.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkInformation.h>

#include "chi_runtime.h"
#include "chi_log.h"

#define ErrorReadingFile(fname)                                                                    \
  std::runtime_error("Failed to open file: " + options.file_name + " in call to " + #fname + ".")

//###################################################################
/**Reads an Exodus unstructured mesh.*/
void
chi_mesh::UnpartitionedMesh::ReadFromExodus(const chi_mesh::UnpartitionedMesh::Options& options)
{
  Chi::log.Log() << "Reading Exodus file: " << options.file_name << ".";

  //======================================== Attempt to open file
  std::ifstream file;
  file.open(options.file_name);
  if (!file.is_open()) throw ErrorReadingFile(ReadFromExodus);
  file.close();

  //======================================== Read the file
  mesh_options_ = options;
  auto reader = vtkSmartPointer<vtkExodusIIReader>::New();
  reader->SetFileName(options.file_name.c_str());

  if (not reader->CanReadFile(options.file_name.c_str()))
    throw std::logic_error("Unable to read file-type with this routine");

  reader->UpdateInformation();
  // Exodus ships boundary-ids via SideSets and NodeSets. This allows
  // it to be read from the file. Here we have to enable the reader
  // to process this because it does not do it by default.
  reader->SetAllArrayStatus(reader->NODE_SET, 1);
  reader->SetAllArrayStatus(reader->NODE_SET_CONN, 1);
  reader->SetAllArrayStatus(reader->SIDE_SET, 1);
  reader->SetAllArrayStatus(reader->SIDE_SET_CONN, 1);

  // The exodusII file format ships blocks of elements
  // together with their points/vertices as self-contained (localized)
  // unstructured meshes. To relate these localized blocks to the original
  // mesh, where are the blocks formed a whole, we need to know the mapping
  // from block-local ids to the original ids. This information can
  // be derived from the GlobalNodeID arrays loaded onto point-data and
  // cell-data. Again, this information is not read by default so we have to
  // turn this on.
  reader->SetGenerateGlobalNodeIdArray(true);
  reader->SetGenerateGlobalElementIdArray(true);
  reader->Update();

  //======================================== Get all the grid blocks
  // This part was quite difficult. I eventually found how to do
  // this from this post:
  // https://public.kitware.com/pipermail/vtkusers/2010-December/064921.html
  // It indicated that all exodus formats are read with a single
  // top level vtkMultiBlockDataSet (Level 1). Each of the blocks in
  // level 2 is also of type vtkMultiBlockDataSet. The level 2 block that has
  // the actual elements is also split into blocks but these, level 3,
  // blocks each contain a structure castable to vtkUnstructuredGrid.
  auto multiblock = reader->GetOutput();

  std::vector<vtkUGridPtrAndName> grid_blocks;
  auto iter_a = multiblock->NewIterator();
  iter_a->GoToFirstItem();
  while (not iter_a->IsDoneWithTraversal())
  {
    auto block_a = iter_a->GetCurrentDataObject();

    const std::string block_name =
      chi::StringTrim(iter_a->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME()));

    if (block_a->GetDataObjectType() == VTK_UNSTRUCTURED_GRID)
    {
      grid_blocks.emplace_back(vtkUnstructuredGrid::SafeDownCast(block_a), block_name);

      Chi::log.Log() << "Reading block " << block_name
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
    chi_mesh::GetBlocksOfDesiredDimension(grid_blocks, max_dimension - 1);

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
    case 1:
      dimension = DIMENSION_1;
      break;
    case 2:
      dimension = DIMENSION_2;
      break;
    case 3:
      dimension = DIMENSION_3;
      break;
    default:
      break;
  }

  attributes_ = dimension | UNSTRUCTURED;

  ComputeCentroidsAndCheckQuality();
  BuildMeshConnectivity();

  //======================================== Set boundary ids
  SetBoundaryIDsFromBlocks(bndry_grid_blocks);

  Chi::log.Log() << "Done reading Exodus file: " << options.file_name << ".";
}