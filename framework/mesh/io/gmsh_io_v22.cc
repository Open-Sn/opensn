// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <fstream>

namespace opensn
{

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromGmshV22(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::FromFile";
  const std::string node_section_name = "$Nodes";
  const std::string elements_section_name = "$Elements";
  const std::string format_section_name = "$MeshFormat";

  // Opening file
  std::ifstream file(options.file_name);
  if (not file.is_open())
    throw std::runtime_error(fname + ": Failed to open file " + options.file_name);

  // Check file format version
  std::string file_line;
  std::istringstream iss;
  while (std::getline(file, file_line))
    if (format_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  double gmsh_version;
  if (!(iss >> gmsh_version) || gmsh_version != 2.2)
    throw std::logic_error(fname + ": Only Gmsh version 4.1 format is supported.");

  auto mesh = std::make_shared<UnpartitionedMesh>();
  log.Log() << "Making unpartitioned mesh from Gmsh file " << options.file_name << " (format v2.2)";

  // Read node data
  file.seekg(0);
  while (std::getline(file, file_line))
    if (node_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_nodes;
  if (not(iss >> num_nodes))
    throw std::logic_error(fname + ": Failed to read the number of nodes.");

  auto& vertices = mesh->GetVertices();
  vertices.clear();
  vertices.resize(num_nodes);

  for (int n = 0; n < num_nodes; n++)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    int vert_index;
    if (not(iss >> vert_index))
      throw std::logic_error(fname + ": Failed to read vertex index.");

    double x, y, z;
    if (not(iss >> x >> y >> z))
      throw std::logic_error(fname + ": Failed while reading vertex coordinates.");
    vertices[vert_index - 1] = {x, y, z};
  }

  auto IsElementType1D = [](int type) { return type == 1; };
  auto IsElementType2D = [](int type) { return type == 2 || type == 3; };
  auto IsElementType3D = [](int type) { return type >= 4 && type <= 7; };
  auto IsElementSupported = [](int type) { return type >= 1 && type <= 7; };
  auto CellTypeFromMSHTypeID = [](int type)
  {
    switch (type)
    {
      case 1:
        return CellType::SLAB;
      case 2:
        return CellType::TRIANGLE;
      case 3:
        return CellType::QUADRILATERAL;
      case 4:
        return CellType::TETRAHEDRON;
      case 5:
        return CellType::HEXAHEDRON;
      case 6:
      case 7:
        return CellType::POLYHEDRON;
      default:
        return CellType::GHOST;
    }
  };
  auto ReadNodes = [&iss, &fname](int num_nodes)
  {
    std::vector<int> raw_nodes(num_nodes, 0);
    for (int i = 0; i < num_nodes; ++i)
      if (not(iss >> raw_nodes[i]))
        throw std::logic_error(fname + ": Failed reading element node index.");

    std::vector<uint64_t> nodes(num_nodes, 0);
    for (int i = 0; i < num_nodes; ++i)
      if ((raw_nodes[i] - 1) >= 0)
        nodes[i] = raw_nodes[i] - 1;
    return nodes;
  };

  // Determine dimension of mesh. Only 2D and 3D meshes are supported. If the mesh is 1D, no
  // elements will be read.
  bool mesh_is_2D = true;
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
    int element_type, num_tags, physical_reg, tag, element_index;

    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    if (not(iss >> element_index >> element_type >> num_tags))
      throw std::logic_error(fname + ": Failed reading element index, type, and number of tags.");

    if (not(iss >> physical_reg))
      throw std::logic_error(fname + ": Failed reading physical region.");

    for (int i = 1; i < num_tags; i++)
      if (not(iss >> tag))
        throw std::logic_error(fname + ": Failed reading element tags.");

    // Skip point type element
    if (element_type == 15)
      continue;

    if (IsElementType3D(element_type))
    {
      mesh_is_2D = false;
      log.Log() << "Mesh identified as 3D.";
      break;
    }

    if (not IsElementSupported(element_type))
      throw std::logic_error(fname + ": Found unsupported element type.");
  }

  // Read element data
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
    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    int element_index, element_type, num_tags;
    if (not(iss >> element_index >> element_type >> num_tags))
      throw std::logic_error(fname + ": Failed reading element index, type, and number of tags.");

    int physical_region;
    if (not(iss >> physical_region))
      throw std::logic_error(fname + ": Failed reading physical region.");

    int tag;
    for (int i = 1; i < num_tags; i++)
      if (not(iss >> tag))
        throw std::logic_error(fname + ": Failed reading tags.");

    // Skip point type elements
    if (element_type == 15)
      continue;

    int num_cell_nodes = 0;
    if (element_type == 1) // 2-node edge
      num_cell_nodes = 2;
    else if (element_type == 2) // 3-node triangle
      num_cell_nodes = 3;
    else if (element_type == 3 or element_type == 4) // 4-node quadrangle or tet
      num_cell_nodes = 4;
    else if (element_type == 5) // 8-node hexahedron
      num_cell_nodes = 8;

    auto& raw_boundary_cells = mesh->GetRawBoundaryCells();
    auto& raw_cells = mesh->GetRawCells();

    // Make the cell on either the volume or the boundary
    std::shared_ptr<UnpartitionedMesh::LightWeightCell> raw_cell;
    if (mesh_is_2D)
    {
      if (IsElementType1D(element_type))
      {
        raw_cell =
          std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::SLAB, CellType::SLAB);
        raw_boundary_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType2D(element_type))
      {
        raw_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(
          CellType::POLYGON, CellTypeFromMSHTypeID(element_type));
        raw_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_cells.";
      }
    }
    else
    {
      if (IsElementType2D(element_type))
      {
        raw_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(
          CellType::POLYGON, CellTypeFromMSHTypeID(element_type));
        raw_boundary_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType3D(element_type))
      {
        raw_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(
          CellType::POLYHEDRON, CellTypeFromMSHTypeID(element_type));
        raw_cells.push_back(raw_cell);
        log.Log0Verbose2() << "Added to raw_cells.";
      }
    }

    if (raw_cell == nullptr)
      continue;

    auto& cell = *raw_cell;
    cell.block_id = physical_region;
    cell.vertex_ids = ReadNodes(num_cell_nodes);

    // Populate faces
    if (element_type == 1) // 2-node edge
    {
      UnpartitionedMesh::LightWeightFace face0;
      UnpartitionedMesh::LightWeightFace face1;

      face0.vertex_ids = {cell.vertex_ids.at(0)};
      face1.vertex_ids = {cell.vertex_ids.at(1)};

      cell.faces.push_back(face0);
      cell.faces.push_back(face1);
    }
    else if (element_type == 2 or element_type == 3) // 3-node triangle or 4-node quadrangle
    {
      size_t num_verts = cell.vertex_ids.size();
      for (size_t e = 0; e < num_verts; e++)
      {
        size_t ep1 = (e < (num_verts - 1)) ? e + 1 : 0;
        UnpartitionedMesh::LightWeightFace face;

        face.vertex_ids = {cell.vertex_ids[e], cell.vertex_ids[ep1]};

        cell.faces.push_back(std::move(face));
      }
    }
    else if (element_type == 4) // 4-node tetrahedron
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
    else if (element_type == 5) // 8-node hexahedron
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

  unsigned int dimension = (mesh_is_2D) ? 2 : 3;
  mesh->SetDimension(dimension);
  mesh->SetType(UNSTRUCTURED);
  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  log.Log() << "Done processing " << options.file_name << ".\n"
            << "Number of nodes read: " << mesh->GetVertices().size() << "\n"
            << "Number of cells read: " << mesh->GetRawCells().size();

  return mesh;
}

} // namespace opensn
