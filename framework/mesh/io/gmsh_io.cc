// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <fstream>

namespace opensn
{

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromGmsh(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::FromFile";

  // Opening the file
  std::ifstream file;
  file.open(options.file_name);
  if (not file.is_open())
  {
    std::ostringstream oss;
    oss << fname + ": Failed to open file: " << options.file_name;
    throw std::runtime_error(oss.str());
  }

  log.Log() << "Making Unpartitioned mesh from msh format file " << options.file_name;
  opensn::mpi_comm.barrier();

  std::shared_ptr<UnpartitionedMesh> mesh = std::make_shared<UnpartitionedMesh>();

  // Declarations
  std::string file_line;
  std::istringstream iss;
  const std::string node_section_name = "$Nodes";
  const std::string elements_section_name = "$Elements";
  const std::string format_section_name = "$MeshFormat";

  // Check the format of this input
  // Role file forward until "$MeshFormat" line is encountered.
  while (std::getline(file, file_line))
    if (format_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  double format;
  if (not(iss >> format))
    throw std::logic_error(fname + ": Failed to read the file format.");
  else if (format != 2.2)
    throw std::logic_error(fname + ": Currently, only msh format 2.2 is supported.");

  // Find section with node information and then read the nodes
  file.seekg(0);
  while (std::getline(file, file_line))
    if (node_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_nodes;
  if (not(iss >> num_nodes))
    throw std::logic_error(fname + ": Failed while trying to read "
                                   "the number of nodes.");

  auto& vertices = mesh->Vertices();
  vertices.clear();
  vertices.resize(num_nodes);

  for (int n = 0; n < num_nodes; n++)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    int vert_index;
    if (not(iss >> vert_index))
      throw std::logic_error(fname + ": Failed to read vertex index.");

    if (not(iss >> vertices[vert_index - 1].x >> vertices[vert_index - 1].y >>
            vertices[vert_index - 1].z))
      throw std::logic_error(fname + ": Failed while reading the vertex "
                                     "coordinates.");
  }

  // Define utility lambdas
  /**Lambda for reading nodes.*/
  auto ReadNodes = [&iss, &fname](int num_nodes)
  {
    std::vector<int> raw_nodes(num_nodes, 0);
    for (int i = 0; i < num_nodes; ++i)
      if (not(iss >> raw_nodes[i]))
        throw std::logic_error(fname + ": Failed when reading element "
                                       "node index.");

    std::vector<uint64_t> nodes(num_nodes, 0);
    for (int i = 0; i < num_nodes; ++i)
      if ((raw_nodes[i] - 1) >= 0)
        nodes[i] = raw_nodes[i] - 1;
    return nodes;
  };

  /**Lamda for checking if an element is 1D.*/
  auto IsElementType1D = [](int element_type)
  {
    if (element_type == 1)
      return true;

    return false;
  };

  /**Lamda for checking if an element is 2D.*/
  auto IsElementType2D = [](int element_type)
  {
    if (element_type == 2 or element_type == 3)
      return true;

    return false;
  };

  /**Lamda for checking if an element is 3D.*/
  auto IsElementType3D = [](int element_type)
  {
    if (element_type >= 4 and element_type <= 7)
      return true;

    return false;
  };

  /**Lambda for checking supported elements.*/
  auto IsElementSupported = [](int element_type)
  {
    if (element_type >= 1 and element_type <= 7)
      return true;

    return false;
  };

  /**Lambda giving the cell subtype, given the MSH cell type.*/
  auto CellTypeFromMSHTypeID = [](int element_type)
  {
    CellType cell_type = CellType::GHOST;

    if (element_type == 1)
      cell_type = CellType::SLAB;
    else if (element_type == 2)
      cell_type = CellType::TRIANGLE;
    else if (element_type == 3)
      cell_type = CellType::QUADRILATERAL;
    else if (element_type == 4)
      cell_type = CellType::TETRAHEDRON;
    else if (element_type == 5)
      cell_type = CellType::HEXAHEDRON;
    else if (element_type == 6 or // Prism
             element_type == 7)   // Pyramid
      cell_type = CellType::POLYHEDRON;

    return cell_type;
  };

  // Determine mesh type 2D/3D
  // Only 2D and 3D meshes are supported. If the mesh
  // is 1D then no elements will be read but the state
  // would still be safe.
  // This section will run through all the elements
  // looking for a 3D element. It will not process
  // any elements.
  bool mesh_is_2D_assumption = true;
  file.seekg(0);
  while (std::getline(file, file_line))
    if (elements_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_elems;
  if (not(iss >> num_elems))
    throw std::logic_error(fname + ": Failed to read number of elements.");

  for (int n = 0; n < num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if (not(iss >> element_index >> elem_type >> num_tags))
      throw std::logic_error(fname + ": Failed while reading element index, "
                                     "element type, and number of tags.");

    if (not(iss >> physical_reg))
      throw std::logic_error(fname + ": Failed while reading physical region.");

    for (int i = 1; i < num_tags; i++)
      if (not(iss >> tag))
        throw std::logic_error(fname + ": Failed when reading tags.");

    if (IsElementType3D(elem_type))
    {
      mesh_is_2D_assumption = false;
      log.Log() << "Mesh identified as 3D.";
      break; // have the answer now leave loop
    }

    if (elem_type == 15) // skip point type element
      continue;

    if (not IsElementSupported(elem_type))
      throw std::logic_error(fname + ": Unsupported element encountered.");
  } // for n

  // Return to the element listing section
  // Now we will actually read the elements.
  file.seekg(0);
  while (std::getline(file, file_line))
    if (elements_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  if (not(iss >> num_elems))
    throw std::logic_error(fname + ": Failed to read number of elements.");

  for (int n = 0; n < num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if (not(iss >> element_index >> elem_type >> num_tags))
      throw std::logic_error(fname + ": Failed while reading element index, "
                                     "element type, and number of tags.");

    if (not(iss >> physical_reg))
      throw std::logic_error(fname + ": Failed while reading physical region.");

    for (int i = 1; i < num_tags; i++)
      if (not(iss >> tag))
        throw std::logic_error(fname + ": Failed when reading tags.");

    if (elem_type == 15) // skip point type elements
      continue;

    if (not IsElementSupported(elem_type))
      throw std::logic_error(fname + ": Unsupported element encountered.");

    log.Log0Verbose2() << "Reading element: " << file_line << " type: " << elem_type;

    int num_cell_nodes;
    if (elem_type == 1)
      num_cell_nodes = 2;
    else if (elem_type == 2) // 3-node triangle
      num_cell_nodes = 3;
    else if (elem_type == 3 or elem_type == 4) // 4-node quadrangle or tet
      num_cell_nodes = 4;
    else if (elem_type == 5) // 8-node hexahedron
      num_cell_nodes = 8;
    else
      continue;

    auto& raw_boundary_cells = mesh->RawBoundaryCells();
    auto& raw_cells = mesh->RawCells();
    // Make the cell on either the volume or the boundary
    UnpartitionedMesh::LightWeightCell* raw_cell = nullptr;
    if (mesh_is_2D_assumption)
    {
      if (IsElementType1D(elem_type))
      {
        raw_cell = new UnpartitionedMesh::LightWeightCell(CellType::SLAB, CellType::SLAB);
        raw_boundary_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType2D(elem_type))
      {
        raw_cell = new UnpartitionedMesh::LightWeightCell(CellType::POLYGON,
                                                          CellTypeFromMSHTypeID(elem_type));
        raw_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_cells.";
      }
    }
    else
    {
      if (IsElementType2D(elem_type))
      {
        raw_cell = new UnpartitionedMesh::LightWeightCell(CellType::POLYGON,
                                                          CellTypeFromMSHTypeID(elem_type));
        raw_boundary_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType3D(elem_type))
      {
        raw_cell = new UnpartitionedMesh::LightWeightCell(CellType::POLYHEDRON,
                                                          CellTypeFromMSHTypeID(elem_type));
        raw_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_cells.";
      }
    }

    if (raw_cell == nullptr)
      continue;

    auto& cell = *raw_cell;
    cell.material_id = physical_reg;
    cell.vertex_ids = ReadNodes(num_cell_nodes);

    // Populate faces
    if (elem_type == 1) // 2-node edge
    {
      UnpartitionedMesh::LightWeightFace face0;
      UnpartitionedMesh::LightWeightFace face1;

      face0.vertex_ids = {cell.vertex_ids.at(0)};
      face1.vertex_ids = {cell.vertex_ids.at(1)};

      cell.faces.push_back(face0);
      cell.faces.push_back(face1);
    }
    else if (elem_type == 2 or elem_type == 3) // 3-node triangle or 4-node quadrangle
    {
      size_t num_verts = cell.vertex_ids.size();
      for (size_t e = 0; e < num_verts; e++)
      {
        size_t ep1 = (e < (num_verts - 1)) ? e + 1 : 0;
        UnpartitionedMesh::LightWeightFace face;

        face.vertex_ids = {cell.vertex_ids[e], cell.vertex_ids[ep1]};

        cell.faces.push_back(std::move(face));
      }                      // for e
    }                        // if 2D elements
    else if (elem_type == 4) // 4-node tetrahedron
    {
      auto& v = cell.vertex_ids;
      std::vector<UnpartitionedMesh::LightWeightFace> lw_faces(4);
      lw_faces[0].vertex_ids = {v[0], v[2], v[1]}; // base-face
      lw_faces[1].vertex_ids = {v[0], v[3], v[2]};
      lw_faces[2].vertex_ids = {v[3], v[1], v[2]};
      lw_faces[3].vertex_ids = {v[3], v[0], v[1]};

      for (auto& lw_face : lw_faces)
        cell.faces.push_back(lw_face);
    }
    else if (elem_type == 5) // 8-node hexahedron
    {
      auto& v = cell.vertex_ids;
      std::vector<UnpartitionedMesh::LightWeightFace> lw_faces(6);
      lw_faces[0].vertex_ids = {v[5], v[1], v[2], v[6]}; // East face
      lw_faces[1].vertex_ids = {v[0], v[4], v[7], v[3]}; // West face
      lw_faces[2].vertex_ids = {v[0], v[3], v[2], v[1]}; // North face
      lw_faces[3].vertex_ids = {v[4], v[5], v[6], v[7]}; // South face
      lw_faces[4].vertex_ids = {v[2], v[3], v[7], v[6]}; // Top face
      lw_faces[5].vertex_ids = {v[0], v[1], v[5], v[4]}; // Bottom face

      for (auto& lw_face : lw_faces)
        cell.faces.push_back(lw_face);
    }
    else
      throw std::runtime_error(fname + ": Unsupported cell type");

  } // for elements

  file.close();

  // Remap material-ids
  std::set<int> material_ids_set_as_read;
  std::map<int, int> material_mapping;

  for (auto& cell : mesh->RawCells())
    material_ids_set_as_read.insert(cell->material_id);

  std::set<int> boundary_ids_set_as_read;
  std::map<int, int> boundary_mapping;

  for (auto& cell : mesh->RawBoundaryCells())
    boundary_ids_set_as_read.insert(cell->material_id);

  {
    int m = 0;
    for (const auto& mat_id : material_ids_set_as_read)
      material_mapping.insert(std::make_pair(mat_id, m++));

    int b = 0;
    for (const auto& bndry_id : boundary_ids_set_as_read)
      boundary_mapping.insert(std::make_pair(bndry_id, b++));
  }

  for (auto& cell : mesh->RawCells())
    cell->material_id = material_mapping[cell->material_id];

  for (auto& cell : mesh->RawBoundaryCells())
    cell->material_id = boundary_mapping[cell->material_id];

  // Always do this
  unsigned int dimension = 2;
  if (not mesh_is_2D_assumption)
    dimension = 3;

  mesh->SetDimension(dimension);
  mesh->SetType(UNSTRUCTURED);

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  log.Log() << "Done processing " << options.file_name << ".\n"
            << "Number of nodes read: " << mesh->Vertices().size() << "\n"
            << "Number of cells read: " << mesh->RawCells().size();

  return mesh;
}

} // namespace opensn
