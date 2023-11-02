#include "chi_meshcontinuum.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/VolumeMesher/chi_volumemesher.h"
#include "mesh/MeshContinuum/chi_grid_vtk_utils.h"

#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkInformation.h>
#include <vtkModelMetadata.h>

#include <vtkExodusIIWriter.h>

#include <vtkPointData.h>
#include <vtkCellData.h>

#include "chi_runtime.h"
#include "chi_log.h"

#define scvtkid static_cast<vtkIdType>

//###################################################################
/**Exports just the portion of the mesh to ExodusII format.*/
void chi_mesh::MeshContinuum::
  ExportCellsToExodus(const std::string& file_base_name,
                      bool suppress_node_sets/*= false*/,
                      bool suppress_side_sets/*= false*/) const
{
  const std::string fname = "chi_mesh::MeshContinuum::ExportCellsToExodus";
  Chi::log.Log() << "Exporting mesh to Exodus file with base " << file_base_name;

  if (Chi::mpi.process_count != 1)
    throw std::logic_error(fname + ": Currently this routine is only allowed "
                                   "in serial.");

  const auto& grid = *this;

  //============================================= Check block consistency
  std::map<int, chi_mesh::CellType> block_id_map;
  for (const auto& cell : local_cells)
  {
    const int mat_id = cell.material_id_;
    if (block_id_map.count(mat_id) == 0)
      block_id_map[mat_id] = cell.SubType();
    else
    {
      if (cell.SubType() != block_id_map.at(mat_id))
        throw std::logic_error(fname + ": Material id " +
                               std::to_string(mat_id) + " appearing for more "
                               "than one cell type.");
    }
  }

  //============================================= Create unstructured meshes
  //                                              for each material-type pair
  vtkNew<vtkMultiBlockDataSet> grid_blocks;
  int max_dimension = 0;
  {
    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkPoints> points;

    points->SetDataType(VTK_DOUBLE);

    vtkNew<vtkIdTypeArray> global_node_id_list;
    global_node_id_list->SetName("GlobalNodeId");

    vtkNew<vtkIdTypeArray> global_elem_id_list;
    global_elem_id_list->SetName("GlobalElementId");

    vtkNew<vtkIntArray> block_id_list;
    block_id_list->SetName("BlockID");

    //============================ Load vertices
    std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
    const size_t num_verts = grid.GetGlobalVertexCount();
    for (size_t v=0; v<num_verts; ++v)
    {
      vertex_map[v] = v;
      const auto& vertex = grid.vertices[v];
      points->InsertNextPoint(vertex.x,vertex.y,vertex.z);

      //Exodus node- and cell indices are 1-based
      //therefore we add a 1 here.
      global_node_id_list->InsertNextValue(scvtkid(v+1));
    }

    //============================ Load cells
    for (const auto& cell : grid.local_cells)
    {
      if (cell.SubType() == CellType::POLYGON or
          cell.SubType() == CellType::POLYHEDRON)
        throw std::logic_error(fname + ": Cell-subtype \"" +
          chi_mesh::CellTypeName(cell.SubType()) + "\" encountered that is not"
          "supported by exodus.");
      chi_mesh::UploadCellGeometryContinuous(cell, vertex_map, ugrid);
      block_id_list->InsertNextValue(cell.material_id_);
      max_dimension =
        std::max(max_dimension,MeshContinuum::GetCellDimension(cell));

      //Exodus node- and cell indices are 1-based
      //therefore we add a 1 here.
      global_elem_id_list->InsertNextValue(scvtkid(cell.global_id_+1));
    }//for local cells

    //============================ Set arrays
    ugrid->SetPoints(points);
    ugrid->GetPointData()->AddArray(global_node_id_list);
    ugrid->GetCellData()->AddArray(global_elem_id_list);
    ugrid->GetCellData()->AddArray(block_id_list);

    ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");
    ugrid->GetCellData()->SetActiveGlobalIds("GlobalElementId");

    //============================ Set block
    grid_blocks->SetBlock(0, ugrid);

    Chi::log.Log()
      << "Writing grid block "
      << " Number of cells: " << ugrid->GetNumberOfCells()
      << " Number of points: " << ugrid->GetNumberOfPoints();
  }//end of grid_blocks assignment

  //============================================= Separate faces by boundary
  //                                              id
  struct FaceInfo
  {
    chi_mesh::CellFace const* face_ptr;
    uint64_t source_cell_id;
    int source_face_id;
  };
  typedef std::vector<FaceInfo> ListOfFaces;
  std::map<uint64_t, ListOfFaces> boundary_id_faces_map;
  for (const auto& cell : grid.local_cells)
  {
    //Here we build a face mapping because ChiTech's face orientation
    //for prisms (wedges) and hexahedrons differ from that of VTK.
    //ChiTech's orientation for prisms and hexes actually matches that
    //of Exodus but VTK assumes the incoming mesh to be conformant to
    //VTK and therefore, internally performs a mapping. Fortunately,
    //the only relevant cell-types, for which a special mapping is
    //required, are the prisms and hexes.
    const size_t num_faces = cell.faces_.size();
    std::vector<int> face_mapping(num_faces, 0);
    if (cell.SubType() == CellType::WEDGE)
      face_mapping = { 2, 3, 4, 0, 1 };
    else if (cell.SubType() == CellType::HEXAHEDRON)
      face_mapping = { 2, 1, 3, 0, 4, 5 };
    else
    {
      for (size_t f=0; f<cell.faces_.size(); ++f)
        face_mapping[f] = static_cast<int>(f);
    }

    //Here we store face information as a triplet, i.e., a face pointer,
    //the id of the cell owning it, and the local face index (relative to
    //the cell) of the face.
    int f=0;
    for (const auto& face : cell.faces_)
    {
      if (not face.has_neighbor_)
        boundary_id_faces_map[face.neighbor_id_].push_back(
          {&face,cell.global_id_,face_mapping[f]});
      ++f;
    }
  }

  //============================================= Make NodeSets and/or SideSets
  vtkNew<vtkMultiBlockDataSet> nodesets_blocks;
  vtkNew<vtkMultiBlockDataSet> sidesets_blocks;
  for (const auto& [bndry_id, face_list] : boundary_id_faces_map)
  {
    const std::string block_name = grid.GetBoundaryIDMap().at(bndry_id);
    Chi::log.Log0Verbose1() << "bid: " + std::to_string(bndry_id) +
                               " name=\"" + block_name + "\"";

    //====================================== NodeSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      points->SetDataType(VTK_DOUBLE);

      vtkNew<vtkIdTypeArray> node_global_ids;
      node_global_ids->SetName("GlobalNodeId");

      //========================== Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids_)
          vid_set.insert(vid);

      //========================== Build vertex map
      std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      //========================== Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x,vertex.y,vertex.z);

        //Exodus node- and cell indices are 1-based
        //therefore we add a 1 here.
        node_global_ids->InsertNextValue(scvtkid(vid+1));
      }

      //========================== Load cells
      for (uint64_t vid : vid_set)
      {
        std::vector<vtkIdType> cell_vids = {scvtkid(vertex_map[vid])};
        ugrid->InsertNextCell(VTK_VERTEX,
                              scvtkid(1),
                              cell_vids.data());
      }

      ugrid->SetPoints(points);
      ugrid->GetPointData()->AddArray(node_global_ids);

      ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");

      nodesets_blocks->SetBlock(bndry_id, ugrid);
      nodesets_blocks->GetMetaData(bndry_id)->
        Set(vtkCompositeDataSet::NAME(), block_name);

      Chi::log.Log()
      << "Writing nodeset block " << block_name
      << " Number of cells: " << ugrid->GetNumberOfCells()
      << " Number of points: " << ugrid->GetNumberOfPoints();
    }

    //====================================== SideSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      points->SetDataType(VTK_DOUBLE);

      vtkNew<vtkIdTypeArray> src_cell_global_ids;
      vtkNew<vtkIntArray>    src_cell_face_id;
      src_cell_global_ids->SetName("SourceElementId");
      src_cell_face_id->SetName("SourceElementSide");

      //========================== Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids_)
          vid_set.insert(vid);

      //========================== Build vertex map
      std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      //========================== Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x,vertex.y,vertex.z);
      }

      //========================== Load faces
      for (const auto& face_info : face_list)
      {
        UploadFaceGeometry(*face_info.face_ptr, vertex_map, ugrid);
        src_cell_global_ids->InsertNextValue(scvtkid(face_info.source_cell_id));
        src_cell_face_id->InsertNextValue(face_info.source_face_id);
      }

      ugrid->SetPoints(points);
      ugrid->GetCellData()->AddArray(src_cell_global_ids);
      ugrid->GetCellData()->AddArray(src_cell_face_id);

      sidesets_blocks->SetBlock(bndry_id, ugrid);
      sidesets_blocks->GetMetaData(bndry_id)->
        Set(vtkCompositeDataSet::NAME(), block_name);

      Chi::log.Log()
        << "Writing sideset block " << block_name
        << " Number of cells: " << ugrid->GetNumberOfCells()
        << " Number of points: " << ugrid->GetNumberOfPoints();
    }//End of side-set
  }

  //============================================= Write the file
  unsigned int next_block = 0;
  vtkNew<vtkMultiBlockDataSet> main_block;
  main_block->SetBlock(next_block++, grid_blocks);
  if (not suppress_node_sets)
  {
    Chi::log.Log0Verbose1() << "Exporting nodeset";
    main_block->SetBlock(next_block, nodesets_blocks);
    main_block->GetMetaData(next_block++)->Set(vtkCompositeDataSet::NAME(), "Node Sets");
  }
  if (not suppress_side_sets)
  {
    Chi::log.Log0Verbose1() << "Exporting sideset";
    main_block->SetBlock(next_block, sidesets_blocks);
    main_block->GetMetaData(next_block++)->Set(vtkCompositeDataSet::NAME(), "Side Sets");
  }

  vtkNew<vtkExodusIIWriter> writer;
  writer->SetBlockIdArrayName("BlockID");

  writer->SetFileName((file_base_name + ".e").c_str());
  writer->SetStoreDoubles(1);

  writer->SetInputData(main_block);

  writer->WriteOutGlobalNodeIdArrayOff();
  writer->WriteOutGlobalElementIdArrayOff();
  writer->WriteOutBlockIdArrayOff();

// The code below requires a VTK patch
//
//  {
//    auto em_in = vtkModelMetadata::New();
//
//    char **dimNames = new char *[3];
//    dimNames[0] = new char[]{"X"};
//    dimNames[1] = new char[]{"Y"};
//    dimNames[2] = new char[]{"Z"};
//
//    max_dimension = std::min(max_dimension, 3);
//    chi::log.Log() << "Max dimension set to " << max_dimension;
//
//    em_in->SetCoordinateNames(max_dimension, dimNames);
//
//    writer->SetModelMetadata(em_in);
//  }

  writer->Write();

  auto em = writer->GetModelMetadata();

  Chi::log.Log() << "Num Blocks   :  " << em->GetNumberOfBlocks();
  Chi::log.Log() << "Num Node Sets:  " << em->GetNumberOfNodeSets();
  Chi::log.Log() << "Num Side Sets:  " << em->GetNumberOfSideSets();
  Chi::log.Log() << "Dimension    :  " << em->GetDimension();

  //writer->PrintSelf(std::cout, vtkIndent());

  Chi::log.Log() << "Done exporting mesh to exodus.";
  Chi::mpi.Barrier();
}