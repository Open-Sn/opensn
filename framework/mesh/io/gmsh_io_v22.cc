// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <cstdint>
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
  const std::string physical_names_section_name = "$PhysicalNames";

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
  double gmsh_version = 0.0;
  if (!(iss >> gmsh_version) or gmsh_version != 2.2)
    throw std::logic_error(fname + ": Only Gmsh version 4.1 format is supported.");

  auto mesh = std::make_shared<UnpartitionedMesh>();
  log.Log() << "Making unpartitioned mesh from Gmsh file " << options.file_name << " (format v2.2)";

  // Read physical section name
  file.seekg(0);
  bool has_physical_names_section = false;
  while (std::getline(file, file_line))
    if (physical_names_section_name == file_line)
    {
      has_physical_names_section = true;
      break;
    }

  // map from physical name ID to [dim, physical name]
  std::map<int, std::tuple<int, std::string>> physical_names;
  if (has_physical_names_section)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);
    int num_physical_names = 0;
    if (not(iss >> num_physical_names))
      throw std::logic_error(fname + ": Failed to read the number of nodes.");

    for (int n = 0; n < num_physical_names; n++)
    {
      std::getline(file, file_line);
      iss = std::istringstream(file_line);

      int dim = 0;
      if (not(iss >> dim))
        throw std::logic_error(fname + ": Failed to read physical name dimension.");

      int phys_name_id = 0;
      if (not(iss >> phys_name_id))
        throw std::logic_error(fname + ": Failed to read physical name ID.");

      std::string name;
      if (not(iss >> name))
        throw std::logic_error(fname + ": Failed to read physical name.");
      name = name.substr(1, name.length() - 2);

      physical_names[phys_name_id] = {dim, name};
    }
  }

  // Read node data
  file.clear();
  file.seekg(0);
  while (std::getline(file, file_line))
    if (node_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  int num_nodes = 0;
  if (not(iss >> num_nodes))
    throw std::logic_error(fname + ": Failed to read the number of nodes.");

  auto& vertices = mesh->GetVertices();
  vertices.clear();
  vertices.resize(num_nodes);

  for (int n = 0; n < num_nodes; n++)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    int vert_index = 0;
    if (not(iss >> vert_index))
      throw std::logic_error(fname + ": Failed to read vertex index.");

    double x = 0.0, y = 0.0, z = 0.0;
    if (not(iss >> x >> y >> z))
      throw std::logic_error(fname + ": Failed while reading vertex coordinates.");
    vertices[vert_index - 1] = {x, y, z};
  }

  auto IsElementType1D = [](int type) { return type == 1; };
  auto IsElementType2D = [](int type) { return type == 2 or type == 3; };
  auto IsElementType3D = [](int type) { return type >= 4 and type <= 7; };
  auto IsElementSupported = [](int type) { return type >= 1 and type <= 7; };
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
  int num_elems = 0;
  if (not(iss >> num_elems))
    throw std::logic_error(fname + ": Failed to read number of elements.");

  for (int n = 0; n < num_elems; n++)
  {
    int element_type = 0, num_tags = 0, physical_reg = 0, tag = 0, element_index = 0;

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

  std::vector<Cell> raw_cells;
  raw_cells.reserve(num_elems);
  std::vector<std::vector<std::uint64_t>> cell_connect;
  cell_connect.reserve(num_elems);
  std::vector<Cell> raw_boundary_cells;
  std::vector<std::vector<std::uint64_t>> bnd_cell_connect;
  std::vector<std::vector<std::vector<std::uint64_t>>> cell_face_connect;
  cell_face_connect.reserve(num_elems);

  for (int n = 0; n < num_elems; n++)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);

    int element_index = 0, element_type = 0, num_tags = 0;
    if (not(iss >> element_index >> element_type >> num_tags))
      throw std::logic_error(fname + ": Failed reading element index, type, and number of tags.");

    unsigned int physical_region = 0;
    if (not(iss >> physical_region))
      throw std::logic_error(fname + ": Failed reading physical region.");

    int tag = 0;
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

    // Make the cell on either the volume or the boundary
    std::optional<Cell> raw_cell;
    bool is_boundary_cell = false;
    if (mesh_is_2D)
    {
      if (IsElementType1D(element_type))
      {
        raw_cell = Cell(CellType::SLAB, CellType::SLAB);
        is_boundary_cell = true;
        log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType2D(element_type))
      {
        raw_cell = Cell(CellType::POLYGON, CellTypeFromMSHTypeID(element_type));
        log.Log0Verbose2() << "Added to raw_cells.";
      }
    }
    else
    {
      if (IsElementType2D(element_type))
      {
        raw_cell = Cell(CellType::POLYGON, CellTypeFromMSHTypeID(element_type));
        is_boundary_cell = true;
        log.Log0Verbose2() << "Added to raw_boundary_cells.";
      }
      else if (IsElementType3D(element_type))
      {
        raw_cell = Cell(CellType::POLYHEDRON, CellTypeFromMSHTypeID(element_type));
        log.Log0Verbose2() << "Added to raw_cells.";
      }
    }

    if (not raw_cell.has_value())
      continue;

    auto& cell = raw_cell.value();
    cell.block_id = physical_region;
    auto cell_vertex_ids = ReadNodes(num_cell_nodes);

    std::vector<std::vector<std::uint64_t>> cell_face_vertex_ids;

    // Populate faces
    if (element_type == 1) // 2-node edge
    {
      CellFace face0;
      CellFace face1;

      std::vector<std::uint64_t> f0_vids = {cell_vertex_ids.at(0)};
      std::vector<std::uint64_t> f1_vids = {cell_vertex_ids.at(1)};

      cell.faces.push_back(face0);
      cell.faces.push_back(face1);
      cell_face_vertex_ids.push_back(std::move(f0_vids));
      cell_face_vertex_ids.push_back(std::move(f1_vids));
    }
    else if (element_type == 2 or element_type == 3) // 3-node triangle or 4-node quadrangle
    {
      size_t num_verts = cell_vertex_ids.size();
      for (size_t e = 0; e < num_verts; e++)
      {
        size_t ep1 = (e < (num_verts - 1)) ? e + 1 : 0;
        CellFace face;

        std::vector<std::uint64_t> f_vids = {cell_vertex_ids[e], cell_vertex_ids[ep1]};

        cell.faces.push_back(std::move(face));
        cell_face_vertex_ids.push_back(std::move(f_vids));
      }
    }
    else if (element_type == 4) // 4-node tetrahedron
    {
      const auto& v = cell_vertex_ids;
      std::vector<CellFace> lw_faces(4);
      cell_face_vertex_ids = {
        {v[0], v[2], v[1]},
        {v[0], v[3], v[2]},
        {v[3], v[1], v[2]},
        {v[3], v[0], v[1]}
      };

      for (auto& lw_face : lw_faces)
        cell.faces.push_back(lw_face);
    }
    else if (element_type == 5) // 8-node hexahedron
    {
      const auto& v = cell_vertex_ids;
      std::vector<CellFace> lw_faces(6);
      cell_face_vertex_ids = {
        {v[5], v[1], v[2], v[6]},
        {v[0], v[4], v[7], v[3]},
        {v[0], v[3], v[2], v[1]},
        {v[4], v[5], v[6], v[7]},
        {v[2], v[3], v[7], v[6]},
        {v[0], v[1], v[5], v[4]}
      };

      for (auto& lw_face : lw_faces)
        cell.faces.push_back(lw_face);
    }
    else
      throw std::runtime_error(fname + ": Unsupported cell type");

    if (is_boundary_cell)
    {
      raw_boundary_cells.emplace_back(cell);
      bnd_cell_connect.emplace_back(cell_vertex_ids);
    }
    else
    {
      raw_cells.emplace_back(cell);
      cell_connect.emplace_back(cell_vertex_ids);
      cell_face_connect.emplace_back(std::move(cell_face_vertex_ids));
    }
  } // for elements

  file.close();

  // create boundary names and IDs
  unsigned int dimension = mesh_is_2D ? 2 : 3;
  for (auto& [id, e] : physical_names)
  {
    if (std::get<0>(e) == dimension - 1)
    {
      auto name = std::get<1>(e);
      mesh->AddBoundary(id, name);
    }
  }

  mesh->SetDimension(dimension);
  mesh->SetType(UNSTRUCTURED);
  mesh->SetCells(std::move(raw_cells), cell_connect);
  mesh->SetCellFaces(cell_face_connect);
  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // remap boundary cells onto cell faces
  std::map<std::set<uint64_t>, unsigned int> bnd_cell_to_bnd_id_map;
  for (std::size_t i = 0; i < raw_boundary_cells.size(); ++i)
  {
    const auto& bnd_cell = raw_boundary_cells[i];
    const auto& bnd_cell_vertex_ids = bnd_cell_connect[i];
    std::set<uint64_t> key;
    for (const auto& vid : bnd_cell_vertex_ids)
      key.insert(vid);
    bnd_cell_to_bnd_id_map[key] = bnd_cell.block_id;
  }
  size_t cell_idx = 0;
  for (auto& cell : mesh->GetCells())
  {
    for (size_t f = 0; f < cell.faces.size(); ++f)
    {
      auto& face = cell.faces[f];
      if (not face.has_neighbor)
      {
        const auto& face_vids = cell_face_connect[cell_idx][f];
        std::set<uint64_t> key(face_vids.begin(), face_vids.end());

        auto it = bnd_cell_to_bnd_id_map.find(key);
        if (it != bnd_cell_to_bnd_id_map.end())
          face.neighbor_id = it->second;
      }
    }
    ++cell_idx;
  }

  log.Log() << "Done processing " << options.file_name << ".\n"
            << "Number of nodes read: " << mesh->GetVertices().size() << "\n"
            << "Number of cells read: " << mesh->GetCells().size();

  return mesh;
}

} // namespace opensn
