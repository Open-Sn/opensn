#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi_comm_set.h"
#include "framework/mpi/mpi.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "framework/mesh/mesh_continuum/grid_vtk_utils.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/mesh/mesh_continuum/grid_face_histogram.h"
#include "framework/data_types/ndarray.h"
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkModelMetadata.h>
#include <vtkExodusIIWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <algorithm>

#define scvtkid static_cast<vtkIdType>

namespace opensn
{

std::shared_ptr<ChiMPICommunicatorSet>
MeshContinuum::MakeMPILocalCommunicatorSet() const
{
  // Build the communicator
  Chi::log.Log0Verbose1() << "Building communicator.";
  std::set<int> local_graph_edges;

  // Loop over local cells
  // Populate local_graph_edges
  local_graph_edges.insert(Chi::mpi.location_id); // add current location
  for (auto& cell : local_cells)
  {
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_)
        if (not face.IsNeighborLocal(*this))
          local_graph_edges.insert(face.GetNeighborPartitionID(*this));
    } // for f
  }   // for local cells

  // Convert set to vector
  // This is just done for convenience because MPI
  // needs a contiguous array
  std::vector<int> local_connections(local_graph_edges.begin(), local_graph_edges.end());

  // Broadcast local connection size
  Chi::log.Log0Verbose1() << "Communicating local connections.";

  std::vector<std::vector<int>> global_graph(Chi::mpi.process_count, std::vector<int>());
  for (int locI = 0; locI < Chi::mpi.process_count; locI++)
  {
    int locI_num_connections = static_cast<int>(local_connections.size());

    // If chi::mpi.location_id == locI then this call will
    // act like a send instead of receive. Otherwise
    // It receives the count.
    MPI_Bcast(&locI_num_connections, 1, MPI_INT, locI, Chi::mpi.comm);

    if (Chi::mpi.location_id != locI) { global_graph[locI].resize(locI_num_connections, -1); }
    else
    {
      std::copy(
        local_connections.begin(), local_connections.end(), std::back_inserter(global_graph[locI]));
    }
  }

  // Broadcast local connections
  for (int locI = 0; locI < Chi::mpi.process_count; locI++)
  {
    // If chi::mpi.location_id == locI then this call will
    // act like a send instead of receive. Otherwise
    // It receives the count.
    MPI_Bcast(global_graph[locI].data(),
              static_cast<int>(global_graph[locI].size()),
              MPI_INT,
              locI,
              Chi::mpi.comm);
  }

  Chi::log.Log0Verbose1() << "Done communicating local connections.";

  // Build groups
  MPI_Group world_group;
  MPI_Comm_group(Chi::mpi.comm, &world_group);

  std::vector<MPI_Group> location_groups;
  location_groups.resize(Chi::mpi.process_count, MPI_Group());

  for (int locI = 0; locI < Chi::mpi.process_count; locI++)
  {
    MPI_Group_incl(world_group,
                   static_cast<int>(global_graph[locI].size()),
                   global_graph[locI].data(),
                   &location_groups[locI]);
  }

  // Build communicators
  std::vector<MPI_Comm> communicators;
  Chi::log.Log0Verbose1() << "Building communicators.";
  communicators.resize(Chi::mpi.process_count, MPI_Comm());

  for (int locI = 0; locI < Chi::mpi.process_count; locI++)
  {
    int err = MPI_Comm_create_group(Chi::mpi.comm, location_groups[locI], 0, &communicators[locI]);

    if (err != MPI_SUCCESS) { Chi::log.Log0Verbose1() << "Communicator creation failed."; }
  }

  Chi::log.Log0Verbose1() << "Done building communicators.";

  return std::make_shared<ChiMPICommunicatorSet>(communicators, location_groups, world_group);
}

void
MeshContinuum::ExportCellsToExodus(const std::string& file_base_name,
                                   bool suppress_node_sets,
                                   bool suppress_side_sets) const
{
  const std::string fname = "MeshContinuum::ExportCellsToExodus";
  Chi::log.Log() << "Exporting mesh to Exodus file with base " << file_base_name;

  if (Chi::mpi.process_count != 1)
    throw std::logic_error(fname + ": Currently this routine is only allowed "
                                   "in serial.");

  const auto& grid = *this;

  // Check block consistency
  std::map<int, CellType> block_id_map;
  for (const auto& cell : local_cells)
  {
    const int mat_id = cell.material_id_;
    if (block_id_map.count(mat_id) == 0) block_id_map[mat_id] = cell.SubType();
    else
    {
      if (cell.SubType() != block_id_map.at(mat_id))
        throw std::logic_error(fname + ": Material id " + std::to_string(mat_id) +
                               " appearing for more "
                               "than one cell type.");
    }
  }

  // Create unstructured meshes for each material-type pair
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

    // Load vertices
    std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
    const size_t num_verts = grid.GetGlobalVertexCount();
    for (size_t v = 0; v < num_verts; ++v)
    {
      vertex_map[v] = v;
      const auto& vertex = grid.vertices[v];
      points->InsertNextPoint(vertex.x, vertex.y, vertex.z);

      // Exodus node- and cell indices are 1-based
      // therefore we add a 1 here.
      global_node_id_list->InsertNextValue(scvtkid(v + 1));
    }

    // Load cells
    for (const auto& cell : grid.local_cells)
    {
      if (cell.SubType() == CellType::POLYGON or cell.SubType() == CellType::POLYHEDRON)
        throw std::logic_error(fname + ": Cell-subtype \"" + CellTypeName(cell.SubType()) +
                               "\" encountered that is not"
                               "supported by exodus.");
      UploadCellGeometryContinuous(cell, vertex_map, ugrid);
      block_id_list->InsertNextValue(cell.material_id_);
      max_dimension = std::max(max_dimension, MeshContinuum::GetCellDimension(cell));

      // Exodus node- and cell indices are 1-based
      // therefore we add a 1 here.
      global_elem_id_list->InsertNextValue(scvtkid(cell.global_id_ + 1));
    } // for local cells

    // Set arrays
    ugrid->SetPoints(points);
    ugrid->GetPointData()->AddArray(global_node_id_list);
    ugrid->GetCellData()->AddArray(global_elem_id_list);
    ugrid->GetCellData()->AddArray(block_id_list);

    ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");
    ugrid->GetCellData()->SetActiveGlobalIds("GlobalElementId");

    // Set block
    grid_blocks->SetBlock(0, ugrid);

    Chi::log.Log() << "Writing grid block "
                   << " Number of cells: " << ugrid->GetNumberOfCells()
                   << " Number of points: " << ugrid->GetNumberOfPoints();
  } // end of grid_blocks assignment

  // Separate faces by boundary id
  struct FaceInfo
  {
    CellFace const* face_ptr;
    uint64_t source_cell_id;
    int source_face_id;
  };
  typedef std::vector<FaceInfo> ListOfFaces;
  std::map<uint64_t, ListOfFaces> boundary_id_faces_map;
  for (const auto& cell : grid.local_cells)
  {
    // Here we build a face mapping because ChiTech's face orientation
    // for prisms (wedges) and hexahedrons differ from that of VTK.
    // ChiTech's orientation for prisms and hexes actually matches that
    // of Exodus but VTK assumes the incoming mesh to be conformant to
    // VTK and therefore, internally performs a mapping. Fortunately,
    // the only relevant cell-types, for which a special mapping is
    // required, are the prisms and hexes.
    const size_t num_faces = cell.faces_.size();
    std::vector<int> face_mapping(num_faces, 0);
    if (cell.SubType() == CellType::WEDGE) face_mapping = {2, 3, 4, 0, 1};
    else if (cell.SubType() == CellType::HEXAHEDRON)
      face_mapping = {2, 1, 3, 0, 4, 5};
    else
    {
      for (size_t f = 0; f < cell.faces_.size(); ++f)
        face_mapping[f] = static_cast<int>(f);
    }

    // Here we store face information as a triplet, i.e., a face pointer,
    // the id of the cell owning it, and the local face index (relative to
    // the cell) of the face.
    int f = 0;
    for (const auto& face : cell.faces_)
    {
      if (not face.has_neighbor_)
        boundary_id_faces_map[face.neighbor_id_].push_back(
          {&face, cell.global_id_, face_mapping[f]});
      ++f;
    }
  }

  // Make NodeSets and/or SideSets
  vtkNew<vtkMultiBlockDataSet> nodesets_blocks;
  vtkNew<vtkMultiBlockDataSet> sidesets_blocks;
  for (const auto& [bndry_id, face_list] : boundary_id_faces_map)
  {
    const std::string block_name = grid.GetBoundaryIDMap().at(bndry_id);
    Chi::log.Log0Verbose1() << "bid: " + std::to_string(bndry_id) + " name=\"" + block_name + "\"";

    // NodeSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      points->SetDataType(VTK_DOUBLE);

      vtkNew<vtkIdTypeArray> node_global_ids;
      node_global_ids->SetName("GlobalNodeId");

      // Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids_)
          vid_set.insert(vid);

      // Build vertex map
      std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      // Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x, vertex.y, vertex.z);

        // Exodus node- and cell indices are 1-based
        // therefore we add a 1 here.
        node_global_ids->InsertNextValue(scvtkid(vid + 1));
      }

      // Load cells
      for (uint64_t vid : vid_set)
      {
        std::vector<vtkIdType> cell_vids = {scvtkid(vertex_map[vid])};
        ugrid->InsertNextCell(VTK_VERTEX, scvtkid(1), cell_vids.data());
      }

      ugrid->SetPoints(points);
      ugrid->GetPointData()->AddArray(node_global_ids);

      ugrid->GetPointData()->SetActiveGlobalIds("GlobalNodeId");

      nodesets_blocks->SetBlock(bndry_id, ugrid);
      nodesets_blocks->GetMetaData(bndry_id)->Set(vtkCompositeDataSet::NAME(), block_name);

      Chi::log.Log() << "Writing nodeset block " << block_name
                     << " Number of cells: " << ugrid->GetNumberOfCells()
                     << " Number of points: " << ugrid->GetNumberOfPoints();
    }

    // SideSet
    {
      vtkNew<vtkUnstructuredGrid> ugrid;
      vtkNew<vtkPoints> points;

      points->SetDataType(VTK_DOUBLE);

      vtkNew<vtkIdTypeArray> src_cell_global_ids;
      vtkNew<vtkIntArray> src_cell_face_id;
      src_cell_global_ids->SetName("SourceElementId");
      src_cell_face_id->SetName("SourceElementSide");

      // Build vertex set
      std::set<uint64_t> vid_set;
      for (const auto& face_info : face_list)
        for (uint64_t vid : face_info.face_ptr->vertex_ids_)
          vid_set.insert(vid);

      // Build vertex map
      std::vector<uint64_t> vertex_map(grid.GetGlobalVertexCount(), 0);
      {
        uint64_t mapped_id = 0;
        for (uint64_t vid : vid_set)
          vertex_map[vid] = mapped_id++;
      }

      // Load vertices
      for (uint64_t vid : vid_set)
      {
        const auto& vertex = grid.vertices[vid];
        points->InsertNextPoint(vertex.x, vertex.y, vertex.z);
      }

      // Load faces
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
      sidesets_blocks->GetMetaData(bndry_id)->Set(vtkCompositeDataSet::NAME(), block_name);

      Chi::log.Log() << "Writing sideset block " << block_name
                     << " Number of cells: " << ugrid->GetNumberOfCells()
                     << " Number of points: " << ugrid->GetNumberOfPoints();
    } // End of side-set
  }

  // Write the file
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

  // writer->PrintSelf(std::cout, vtkIndent());

  Chi::log.Log() << "Done exporting mesh to exodus.";
  Chi::mpi.Barrier();
}

void
MeshContinuum::ExportCellsToObj(const char* fileName, bool per_material, int options) const
{
  if (!per_material)
  {
    FILE* of = fopen(fileName, "w");

    if (of == nullptr)
    {
      Chi::log.LogAllError() << "Could not open file: " << std::string(fileName);
      Chi::Exit(EXIT_FAILURE);
    }

    // Develop list of faces and nodes
    std::set<int> nodes_set;
    std::vector<CellFace> faces_to_export;
    for (auto& cell : local_cells)
    {
      if (cell.Type() == CellType::POLYHEDRON)
      {
        for (auto& face : cell.faces_)
        {
          if (not face.has_neighbor_)
          {
            faces_to_export.push_back(face);

            for (int vid : face.vertex_ids_)
              nodes_set.insert(vid);
          } // if boundary
        }   // for face
      }     // if polyhedron
    }       // for local cell

    // Write header
    fprintf(of, "# Exported mesh file from Extrusion script\n");
    std::string str_file_name(fileName);
    std::string file_base_name = str_file_name.substr(0, str_file_name.find('.'));
    fprintf(of, "o %s\n", file_base_name.c_str());

    // Develop node mapping and write them
    std::vector<int> node_mapping(GetGlobalVertexCount(), -1);

    int node_counter = 0;
    for (auto node : nodes_set)
    {
      node_counter++;
      int node_g_index = node;
      node_mapping[node_g_index] = node_counter;

      Vertex cur_v = vertices[node_g_index];

      fprintf(of, "v %9.6f %9.6f %9.6f\n", cur_v.x, cur_v.y, cur_v.z);
    }

    // Write face normals
    for (const auto& face : faces_to_export)
    {
      fprintf(of, "vn %.4f %.4f %.4f\n", face.normal_.x, face.normal_.y, face.normal_.z);
    }

    // Write faces
    int normal_counter = 0;
    for (const auto& face : faces_to_export)
    {
      normal_counter++;
      fprintf(of, "f");

      for (auto v_g_index : face.vertex_ids_)
        fprintf(of, " %d//%d", node_mapping[v_g_index], normal_counter);

      fprintf(of, "\n");
    }

    fclose(of);

    Chi::log.Log() << "Exported Volume mesh to " << str_file_name;
  } // Whole mesh
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PER MATERIAL
  else
  {
    // Get base name
    std::string str_file_name(fileName);
    std::string file_base_name = str_file_name.substr(0, str_file_name.find('.'));

    if (Chi::material_stack.empty())
    {
      Chi::log.Log0Warning() << "ExportCellsToObj: No mesh will be exported because there "
                             << "are no physics materials present";
    }

    for (int mat = 0; mat < Chi::material_stack.size(); mat++)
    {
      std::string mat_base_name = file_base_name + std::string("_m") + std::to_string(mat);
      std::string mat_file_name = mat_base_name + std::string(".obj");
      FILE* of = fopen(mat_file_name.c_str(), "w");

      if (of == nullptr)
      {
        Chi::log.LogAllError() << "Could not open file: " << mat_file_name;
        Chi::Exit(EXIT_FAILURE);
      }

      // Develop list of faces and nodes
      std::set<int> nodes_set;
      std::vector<CellFace> faces_to_export;
      for (const auto& cell : local_cells)
      {
        if (cell.Type() == CellType::POLYHEDRON)
        {
          if (cell.material_id_ != mat) continue;

          for (const auto& face : cell.faces_)
          {
            int adjcell_glob_index = face.neighbor_id_;

            if (adjcell_glob_index < 0)
            {
              faces_to_export.push_back(face);

              for (auto vid : face.vertex_ids_)
                nodes_set.insert(vid);
            } // if boundary
            else
            {
              auto& adj_cell = cells[adjcell_glob_index];

              if (adj_cell.material_id_ != mat)
              {
                faces_to_export.push_back(face);

                for (auto vid : face.vertex_ids_)
                  nodes_set.insert(vid);
              } // if material missmatch
            }   // if neigbor cell
          }     // for face
        }       // if polyhedron
      }         // for local cell

      // Write header
      fprintf(of, "# Exported mesh file from Extrusion script\n");
      fprintf(of, "o %s\n", mat_base_name.c_str());

      // Develop node mapping and write
      // them
      std::vector<int> node_mapping(GetGlobalVertexCount(), -1);

      int node_counter = 0;
      for (auto node : nodes_set)
      {
        node_counter++;
        int node_g_index = node;
        node_mapping[node_g_index] = node_counter;

        Vertex cur_v = vertices[node_g_index];

        fprintf(of, "v %9.6f %9.6f %9.6f\n", cur_v.x, cur_v.y, cur_v.z);
      }

      // Write face normals
      for (const auto& face : faces_to_export)
      {
        fprintf(of, "vn %.4f %.4f %.4f\n", face.normal_.x, face.normal_.y, face.normal_.z);
      }

      // Write faces
      int normal_counter = 0;
      for (const auto& face : faces_to_export)
      {
        normal_counter++;
        fprintf(of, "f");

        for (auto v_g_index : face.vertex_ids_)
          fprintf(of, " %d//%d", node_mapping[v_g_index], normal_counter);

        fprintf(of, "\n");
      }

      fclose(of);

      Chi::log.Log() << "Exported Material Volume mesh to " << mat_file_name;
    } // for mat
  }   // if per material
}

void
MeshContinuum::ExportCellsToVTK(const std::string& file_base_name) const
{
  Chi::log.Log() << "Exporting mesh to VTK files with base " << file_base_name;

  const auto& grid = *this;

  auto ugrid = PrepareVtkUnstructuredGrid(grid, false);

  WritePVTUFiles(ugrid, file_base_name);

  Chi::log.Log() << "Done exporting mesh to VTK.";
}

std::vector<uint64_t>
MeshContinuum::GetDomainUniqueBoundaryIDs() const
{
  Chi::mpi.Barrier();
  Chi::log.Log() << "Identifying unique boundary-ids.";

  // Develop local bndry-id set
  std::set<uint64_t> local_bndry_ids_set;
  for (auto& cell : local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_) local_bndry_ids_set.insert(face.neighbor_id_);

  // Vectorify it and get local count
  std::vector<uint64_t> local_bndry_ids(local_bndry_ids_set.begin(), local_bndry_ids_set.end());
  int local_num_bndry_ids = (int)local_bndry_ids.size();

  // Everyone now tells everyone
  //                                       how many bndry-ids they have
  std::vector<int> locI_bndry_count(Chi::mpi.process_count, 0);

  MPI_Allgather(
    &local_num_bndry_ids, 1, MPI_INT, locI_bndry_count.data(), 1, MPI_INT, Chi::mpi.comm);

  // Build a displacement list, in prep for gathering all bndry-ids
  std::vector<int> locI_bndry_ids_displs(Chi::mpi.process_count, 0);
  size_t total_num_global_bndry_ids = locI_bndry_count[0];
  for (int locI = 1; locI < Chi::mpi.process_count; ++locI)
  {
    locI_bndry_ids_displs[locI] = locI_bndry_ids_displs[locI - 1] + locI_bndry_count[locI - 1];
    total_num_global_bndry_ids += locI_bndry_count[locI];
  }

  // Everyone now sends everyone their boundary-ids
  std::vector<uint64_t> globl_bndry_ids(total_num_global_bndry_ids);

  MPI_Allgatherv(local_bndry_ids.data(),
                 local_num_bndry_ids,
                 MPI_UNSIGNED_LONG_LONG,
                 globl_bndry_ids.data(),
                 locI_bndry_count.data(),
                 locI_bndry_ids_displs.data(),
                 MPI_UNSIGNED_LONG_LONG,
                 Chi::mpi.comm);

  std::set<uint64_t> globl_bndry_ids_set(globl_bndry_ids.begin(), globl_bndry_ids.end());

  std::vector<uint64_t> unique_bdnry_ids(globl_bndry_ids_set.begin(), globl_bndry_ids_set.end());
  return unique_bdnry_ids;
}

std::shared_ptr<GridFaceHistogram>
MeshContinuum::MakeGridFaceHistogram(double master_tolerance, double slave_tolerance) const
{
  std::vector<std::pair<size_t, size_t>> face_categories_list;
  // Fill histogram
  std::vector<size_t> face_size_histogram;
  for (const auto& cell : local_cells)
    for (const auto& face : cell.faces_)
      face_size_histogram.push_back(face.vertex_ids_.size());

  std::stable_sort(face_size_histogram.begin(), face_size_histogram.end());

  // Determine total face dofs
  size_t total_face_dofs_count = 0;
  for (auto face_size : face_size_histogram)
    total_face_dofs_count += face_size;

  // Compute average and ratio
  size_t smallest_face = face_size_histogram.front();
  size_t largest_face = face_size_histogram.back();
  size_t total_num_faces = face_size_histogram.size();
  double average_dofs_per_face = (double)total_face_dofs_count / (double)total_num_faces;

  std::stringstream outstr;
  outstr << "\nSmallest face = " << smallest_face;
  outstr << "\nLargest face = " << largest_face;
  outstr << "\nTotal face dofs = " << total_face_dofs_count;
  outstr << "\nTotal faces = " << face_size_histogram.size();
  outstr << "\nAverage dofs/face = " << average_dofs_per_face;
  outstr << "\nMax to avg ratio = " << (double)largest_face / average_dofs_per_face;
  Chi::log.LogAllVerbose2() << outstr.str();

  // Determine number of bins
  size_t last_bin_num_faces = total_num_faces;
  if (((double)largest_face / average_dofs_per_face) > master_tolerance)
  {
    Chi::log.LogAllVerbose2() << "The ratio of max face dofs to average face dofs "
                              << "is larger than " << master_tolerance
                              << ", therefore a binned histogram "
                              << "will be constructed.";

    // Build categories
    size_t running_total_face_dofs = 0;
    size_t running_face_count = 0;
    size_t running_face_size = face_size_histogram[0];

    double running_average = (double)face_size_histogram[0];

    for (size_t f = 0; f < total_num_faces; ++f)
    {
      if (((double)face_size_histogram[f] / running_average) > slave_tolerance)
      {
        face_categories_list.emplace_back(running_face_size, running_face_count);
        running_total_face_dofs = 0;
        running_face_count = 0;
      }

      running_face_size = face_size_histogram[f];
      running_total_face_dofs += face_size_histogram[f];
      running_face_count++;
      running_average = (double)running_total_face_dofs / double(running_face_count);
      last_bin_num_faces = running_face_count;
    }
  }
  face_categories_list.emplace_back(largest_face, last_bin_num_faces);

  // Verbose print bins
  outstr.str(std::string());
  outstr << "A total of " << face_categories_list.size() << " bins were created:\n";

  size_t bin_counter = -1;
  for (auto bins : face_categories_list)
  {
    outstr << "Bin " << ++bin_counter << ": " << bins.second << " faces with max face dofs "
           << bins.first << "\n";
  }

  Chi::log.LogAllVerbose2() << outstr.str();

  return std::make_shared<GridFaceHistogram>(face_categories_list);
}

bool
MeshContinuum::IsCellLocal(uint64_t cell_global_index) const
{
  auto native_index = global_cell_id_to_local_id_map_.find(cell_global_index);

  if (native_index != global_cell_id_to_local_id_map_.end()) return true;

  return false;
}

int
MeshContinuum::GetCellDimension(const Cell& cell)
{
  switch (cell.Type())
  {
    case CellType::POINT:
    case CellType::GHOST:
      return 0;
    case CellType::SLAB:
      return 1;
    case CellType::TRIANGLE:
    case CellType::QUADRILATERAL:
    case CellType::POLYGON:
      return 2;
    case CellType::TETRAHEDRON:
    case CellType::HEXAHEDRON:
    case CellType::WEDGE:
    case CellType::PYRAMID:
    case CellType::POLYHEDRON:
      return 3;
    default:
      throw std::logic_error("MeshContinuum::GetCellDimension: "
                             "Dimension mapping unavailable for cell type.");
  }
  return false;
}

void
MeshContinuum::FindAssociatedVertices(const CellFace& cur_face,
                                      std::vector<short>& dof_mapping) const
{
  const int associated_face = cur_face.GetNeighborAssociatedFace(*this);
  // Check face validity
  ChiLogicalErrorIf(not cur_face.has_neighbor_,
                    "Invalid cell index encountered in call to "
                    "MeshContinuum::FindAssociatedVertices. Index "
                    "points to a boundary");

  auto& adj_cell = cells[cur_face.neighbor_id_];

  dof_mapping.reserve(cur_face.vertex_ids_.size());

  const auto& adj_face = adj_cell.faces_[associated_face];

  for (auto cfvid : cur_face.vertex_ids_)
  {
    bool found = false;
    short afv = 0;
    for (auto afvid : adj_face.vertex_ids_)
    {
      if (cfvid == afvid)
      {
        dof_mapping.push_back((short)afv);
        found = true;
        break;
      }
      afv++;
    }

    if (!found)
    {
      Chi::log.LogAllError() << "Face DOF mapping failed in call to "
                             << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
                                "node."
                             << cur_face.neighbor_id_ << " " << cur_face.centroid_.PrintS();
      Chi::Exit(EXIT_FAILURE);
    }
  }
}

void
MeshContinuum::FindAssociatedCellVertices(const CellFace& cur_face,
                                          std::vector<short>& dof_mapping) const
{
  // Check face validity
  ChiLogicalErrorIf(not cur_face.has_neighbor_,
                    "Invalid cell index encountered in call to "
                    "MeshContinuum::FindAssociatedVertices. Index "
                    "points to a boundary");

  auto& adj_cell = cells[cur_face.neighbor_id_];

  dof_mapping.reserve(cur_face.vertex_ids_.size());

  for (auto cfvid : cur_face.vertex_ids_)
  {
    bool found = false;
    short acv = 0;
    for (auto acvid : adj_cell.vertex_ids_)
    {
      if (cfvid == acvid)
      {
        dof_mapping.push_back(acv);
        found = true;
        break;
      }
      ++acv;
    }

    if (!found)
    {
      Chi::log.LogAllError() << "Face DOF mapping failed in call to "
                             << "MeshContinuum::FindAssociatedVertices. Could not find a matching"
                                "node."
                             << cur_face.neighbor_id_ << " " << cur_face.centroid_.PrintS();
      Chi::Exit(EXIT_FAILURE);
    }
  }
}

size_t
MeshContinuum::MapCellFace(const Cell& cur_cell, const Cell& adj_cell, unsigned int f)
{
  const auto& ccface = cur_cell.faces_[f]; // current cell face
  std::set<uint64_t> ccface_vids;
  for (auto vid : ccface.vertex_ids_)
    ccface_vids.insert(vid);

  size_t fmap;
  bool map_found = false;
  for (size_t af = 0; af < adj_cell.faces_.size(); af++)
  {
    const auto& acface = adj_cell.faces_[af]; // adjacent cell face

    std::set<uint64_t> acface_vids;
    for (auto vid : acface.vertex_ids_)
      acface_vids.insert(vid);

    if (acface_vids == ccface_vids)
    {
      fmap = af;
      map_found = true;
      break;
    }
  } // for adj faces

  if (not map_found) throw std::logic_error("MeshContinuum::MapCellFace: Mapping failure.");

  return fmap;
}

size_t
MeshContinuum::MapCellGlobalID2LocalID(uint64_t global_id) const
{
  return global_cell_id_to_local_id_map_.at(global_id);
}

Vector3
MeshContinuum::ComputeCentroidFromListOfNodes(const std::vector<uint64_t>& list) const
{
  if (list.empty())
  {
    Chi::log.LogAllError() << "ComputeCentroidFromListOfNodes, empty list";
    Chi::Exit(EXIT_FAILURE);
  }
  Vector3 centroid;
  for (auto node_id : list)
    centroid = centroid + vertices[node_id];

  return centroid / double(list.size());
}

size_t
MeshContinuum::CountCellsInLogicalVolume(const LogicalVolume& log_vol) const
{
  size_t local_count = 0;
  for (const auto& cell : local_cells)
    if (log_vol.Inside(cell.centroid_)) ++local_count;

  size_t global_count = 0;

  MPI_Allreduce(&local_count, &global_count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chi::mpi.comm);

  return global_count;
}

bool
MeshContinuum::CheckPointInsideCell(const Cell& cell, const Vector3& point) const
{
  const auto& grid_ref = *this;
  typedef Vector3 Vec3;
  auto InsideTet = [](const Vec3& point, const Vec3& v0, const Vec3& v1, const Vec3& v2)
  {
    const auto& v01 = v1 - v0;
    const auto& v02 = v2 - v0;

    const auto n = v01.Cross(v02).Normalized();
    const auto c = (v0 + v1 + v2) / 3.0;

    const auto pc = point - c;

    if (pc.Dot(n) > 0.0) return true;
    else
      return false;
  };

  bool inside = true;
  if (cell.Type() == CellType::SLAB)
  {
    const auto& v0 = grid_ref.vertices[cell.vertex_ids_[0]];
    const auto& v1 = grid_ref.vertices[cell.vertex_ids_[1]];

    const auto v01 = v1 - v0;
    const auto v0p = point - v0;

    const double v0p_dot_v01 = v0p.Dot(v01);

    if (not(v0p_dot_v01 >= 0 and v0p_dot_v01 < v01.Norm())) inside = false;
  } // slab

  else if (cell.Type() == CellType::POLYGON)
  {
    for (const auto& face : cell.faces_)
    {
      const auto& vcp = point - face.centroid_;

      if (vcp.Dot(face.normal_) > 0)
      {
        inside = false;
        break;
      }
    } // for face
  }   // polygon

  else if (cell.Type() == CellType::POLYHEDRON)
  {
    inside = false;
    // form tetra hedrons
    const auto& vcc = cell.centroid_;
    for (const auto& face : cell.faces_)
    {
      const auto& vfc = face.centroid_;

      const size_t num_sides = face.vertex_ids_.size();
      for (size_t s = 0; s < num_sides; ++s)
      {
        const size_t sp1 = (s < (num_sides - 1)) ? s + 1 : 0;
        const auto& v0 = grid_ref.vertices[face.vertex_ids_[s]];
        const auto& v1 = vfc;
        const auto& v2 = grid_ref.vertices[face.vertex_ids_[sp1]];
        const auto& v3 = vcc;

        typedef std::tuple<Vec3, Vec3, Vec3> TetFace;

        std::vector<TetFace> tet_faces;
        tet_faces.emplace_back(v0, v1, v2);
        tet_faces.emplace_back(v0, v2, v3);
        tet_faces.emplace_back(v1, v3, v2);
        tet_faces.emplace_back(v0, v3, v1);

        bool inside_tet = true;
        for (const auto& tet_face : tet_faces)
        {
          if (not InsideTet(
                point, std::get<0>(tet_face), std::get<1>(tet_face), std::get<2>(tet_face)))
          {
            inside_tet = false;
            break;
          }
        } // for triangular tet_face
        if (inside_tet)
        {
          inside = true;
          break;
        }
      } // for side
      if (inside) break;
    } // for face
  }   // polyhedron
  else
    throw std::logic_error("MeshContinuum::CheckPointInsideCell: "
                           "Unsupported cell-type encountered.");

  return inside;
}

std::array<size_t, 3>
MeshContinuum::GetIJKInfo() const
{
  const std::string fname = "GetIJKInfo";
  if (not(this->Attributes() & MeshAttributes::ORTHOGONAL))
    throw std::logic_error(fname + " can only be run on orthogonal meshes.");

  return {ortho_attributes.Nx, ortho_attributes.Ny, ortho_attributes.Nz};
}

NDArray<uint64_t>
MeshContinuum::MakeIJKToGlobalIDMapping() const
{
  const std::string fname = "MakeIJKToGlobalIDMapping";
  if (not(this->Attributes() & MeshAttributes::ORTHOGONAL))
    throw std::logic_error(fname + " can only be run on orthogonal meshes.");

  const auto ijk_info = this->GetIJKInfo();
  const auto Nx = static_cast<int64_t>(ijk_info[0]);
  const auto Ny = static_cast<int64_t>(ijk_info[1]);
  const auto Nz = static_cast<int64_t>(ijk_info[2]);

  NDArray<uint64_t> m_ijk_to_i({Nx, Ny, Nz});
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k)
        m_ijk_to_i(i, j, k) = static_cast<uint64_t>(m_ijk_to_i.MapNDtoLin(i, j, k));

  return m_ijk_to_i;
}

std::vector<Vector3>
MeshContinuum::MakeCellOrthoSizes() const
{
  std::vector<Vector3> cell_ortho_sizes(local_cells.size());
  for (const auto& cell : local_cells)
  {
    Vector3 vmin = vertices[cell.vertex_ids_.front()];
    Vector3 vmax = vmin;

    for (const auto vid : cell.vertex_ids_)
    {
      const auto& vertex = vertices[vid];
      vmin.x = std::min(vertex.x, vmin.x);
      vmin.y = std::min(vertex.y, vmin.y);
      vmin.z = std::min(vertex.z, vmin.z);

      vmax.x = std::max(vertex.x, vmax.x);
      vmax.y = std::max(vertex.y, vmax.y);
      vmax.z = std::max(vertex.z, vmax.z);
    }

    cell_ortho_sizes[cell.local_id_] = vmax - vmin;
  } // for cell

  return cell_ortho_sizes;
}

uint64_t
MeshContinuum::MakeBoundaryID(const std::string& boundary_name) const
{
  if (boundary_id_map_.empty()) return 0;

  for (const auto& [id, name] : boundary_id_map_)
    if (boundary_name == name) return id;

  uint64_t max_id = 0;
  for (const auto& [id, name] : boundary_id_map_)
    max_id = std::max(id, max_id);

  return max_id + 1;
}

std::pair<Vector3, Vector3>
MeshContinuum::GetLocalBoundingBox() const
{
  Vector3 xyz_min;
  Vector3 xyz_max;

  auto Vec3Min = [](const Vector3& xyz_A, const Vector3& xyz_B)
  {
    return Vector3(
      std::min(xyz_A.x, xyz_B.x), std::min(xyz_A.y, xyz_B.y), std::min(xyz_A.z, xyz_B.z));
  };
  auto Vec3Max = [](const Vector3& xyz_A, const Vector3& xyz_B)
  {
    return Vector3(
      std::max(xyz_A.x, xyz_B.x), std::max(xyz_A.y, xyz_B.y), std::max(xyz_A.z, xyz_B.z));
  };

  bool initialized = false;
  for (const auto& cell : local_cells)
  {
    for (const uint64_t vid : cell.vertex_ids_)
    {
      const auto& vertex = vertices[vid];
      if (not initialized)
      {
        xyz_min = vertex;
        xyz_max = vertex;
        initialized = true;
      }
      xyz_min = Vec3Min(xyz_min, vertex);
      xyz_max = Vec3Max(xyz_max, vertex);
    }
  }
  return {xyz_min, xyz_max};
}

size_t
MeshContinuum::GetGlobalNumberOfCells() const
{
  size_t num_local_cells = local_cells_.size();
  size_t num_globl_cells = 0;

  MPI_Allreduce(
    &num_local_cells, &num_globl_cells, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, Chi::mpi.comm);

  return num_globl_cells;
}

} // namespace opensn
