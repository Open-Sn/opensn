// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include <fstream>
#include <cstdint>
#include <limits>
#include <cctype>
#include <array>
#include <cstring>
#include <utility>

namespace opensn
{
namespace
{
template <typename T>
T
ReadBinary(std::ifstream& file)
{
  std::array<char, sizeof(T)> buffer{};
  file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  T value{};
  std::memcpy(&value, buffer.data(), buffer.size());
  return value;
}

uint32_t
ReadUInt32(std::ifstream& file)
{
  return ReadBinary<uint32_t>(file);
}

uint64_t
ReadSize(std::ifstream& file)
{
  return ReadBinary<uint64_t>(file);
}

double
ReadReal(std::ifstream& file, int data_size)
{
  if (data_size == 4)
  {
    const auto value = ReadBinary<float>(file);
    return static_cast<double>(value);
  }
  if (data_size == 8)
  {
    return ReadBinary<double>(file);
  }
  throw std::logic_error("Unsupported Gmsh data size.");
}

void
ExpectEndSection(std::ifstream& file, const std::string& end_name)
{
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::string line;
  std::getline(file, line);
  if (line != end_name)
    throw std::logic_error("Missing expected section terminator " + end_name);
}

std::string
ReadBinaryString(std::ifstream& file, int data_size)
{
  const auto length = static_cast<size_t>(ReadUInt32(file));
  std::string value(length, '\0');
  if (length > 0)
    file.read(value.data(), static_cast<std::streamsize>(length));
  return value;
}

bool
NextIsAsciiNumeric(std::ifstream& file)
{
  const int c = file.peek();
  if (c == '-' || c == '+')
    return true;
  return std::isdigit(c);
}

} // namespace

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromGmshV41(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::FromGmshV41";
  const std::string format_section_name = "$MeshFormat";

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
  int file_type = 0;
  int data_size = 0;
  if (!(iss >> gmsh_version >> file_type >> data_size) || gmsh_version != 4.1)
    throw std::logic_error(fname + ": Only Gmsh version 4.1 format is supported.");

  file.close();

  if (file_type == 1)
    return FromGmshV41Binary(options, data_size);
  if (file_type == 0)
    return FromGmshV41ASCII(options);

  throw std::logic_error(fname + ": Unsupported Gmsh file type.");
}

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromGmshV41ASCII(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::FromGmshV41ASCII";
  const std::string node_section_name = "$Nodes";
  const std::string element_section_name = "$Elements";
  const std::string entities_section_name = "$Entities";
  const std::string partitioned_entities_section_name = "$PartitionedEntities";
  const std::string physical_names_section_name = "$PhysicalNames";

  // Open file
  std::ifstream file(options.file_name);
  if (not file.is_open())
    throw std::runtime_error(fname + ": Failed to open file " + options.file_name);

  std::string file_line;
  std::istringstream iss;

  // Check for unsupported $PartitionedEntities section
  while (std::getline(file, file_line))
  {
    if (partitioned_entities_section_name == file_line)
    {
      log.Log0Warning() << "Found unsuported $PartitionedEntities section in " + options.file_name;
      break;
    }
  }
  file.clear();
  file.seekg(0);

  auto mesh = std::make_shared<UnpartitionedMesh>();
  log.Log() << "Making unpartitioned mesh from Gmsh file " << options.file_name << " (format v4.1)";

  // Read physical section name
  file.seekg(0);
  bool has_physical_names_section = false;
  while (std::getline(file, file_line))
    if (physical_names_section_name == file_line)
    {
      has_physical_names_section = true;
      break;
    }

  // Map from physical name ID to [dim, physical name]
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

  // Read $Entities section
  file.clear();
  file.seekg(0);
  while (std::getline(file, file_line))
    if (entities_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  size_t num_points = 0, num_curves = 0, num_surfaces = 0, num_volumes = 0;
  if (!(iss >> num_points >> num_curves >> num_surfaces >> num_volumes))
    throw std::logic_error(fname + ": Failed to read number of entitities.");

  // Skip point entities
  for (size_t i = 0; i < num_points; ++i)
    std::getline(file, file_line);

  // Read curve, surface, and volume entities
  auto read_entity = [&](std::map<int, int>& entity_map, size_t num_entities)
  {
    for (size_t i = 0; i < num_entities; ++i)
    {
      int entity_tag = 0, physical_tag = 0;
      size_t num_physical_tags = 0;
      double minx = 0.0, miny = 0.0, minz = 0.0, maxx = 0.0, maxy = 0.0, maxz = 0.0;
      std::getline(file, file_line);
      iss = std::istringstream(file_line);
      iss >> entity_tag >> minx >> miny >> minz >> maxx >> maxy >> maxz >> num_physical_tags >>
        physical_tag;
      entity_map[entity_tag] = physical_tag;
    }
  };
  std::map<int, int> curve_entities, surface_entities, volume_entities;
  read_entity(curve_entities, num_curves);
  read_entity(surface_entities, num_surfaces);
  read_entity(volume_entities, num_volumes);

  // Read node data
  file.seekg(0);
  while (std::getline(file, file_line))
    if (node_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  size_t num_entity_blocks = 0, num_nodes = 0, min_node_tag = 0, max_node_tag = 0;
  if (not(iss >> num_entity_blocks >> num_nodes >> min_node_tag >> max_node_tag))
    throw std::logic_error(fname + ": Failed to read the number of node entity blocks.");

  auto& vertices = mesh->GetVertices();
  vertices.clear();
  vertices.resize(num_nodes);

  for (size_t n = 0; n < num_entity_blocks; ++n)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);
    int entity_dim = 0, entity_tag = 0, parametric = 0;
    size_t num_nodes_in_block = 0;
    if (not(iss >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block))
    {
      throw std::logic_error(fname + ": Failed to read number of nodes in block " +
                             std::to_string(n) + ".");
    }

    std::vector<size_t> node_tags(num_nodes_in_block);
    for (size_t i = 0; i < num_nodes_in_block; ++i)
    {
      std::getline(file, file_line);
      iss = std::istringstream(file_line);
      iss >> node_tags[i];
    }

    for (size_t i = 0; i < num_nodes_in_block; ++i)
    {
      std::getline(file, file_line);
      iss = std::istringstream(file_line);
      double x = 0.0, y = 0.0, z = 0.0;
      iss >> x >> y >> z;
      vertices[node_tags[i] - 1] = {x, y, z};
    }
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

  // Determine dimension of mesh. Only 2D and 3D meshes are supported. If the mesh is 1D, no
  // elements will be read.
  bool mesh_is_2D = true;
  file.seekg(0);
  while (std::getline(file, file_line))
    if (element_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  size_t num_elements = 0, min_element_tag = 0, max_element_tag = 0;
  if (not(iss >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag))
    throw std::logic_error(fname + ": Failed to read number of element entity blocks.");

  for (size_t n = 0; n < num_entity_blocks; ++n)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);
    int entity_dim = 0, entity_tag = 0, element_type = 0;
    size_t num_elems_in_block = 0;
    if (not(iss >> entity_dim >> entity_tag >> element_type >> num_elems_in_block))
    {
      throw std::logic_error(fname + ": Failed to read number of elements in block " +
                             std::to_string(n));
    }

    // Skip point type element
    if (element_type == 15)
    {
      for (size_t i = 0; i < num_elems_in_block; ++i)
        std::getline(file, file_line);
      continue;
    }

    if (IsElementType3D(element_type))
    {
      mesh_is_2D = false;
      break;
    }

    if (not IsElementSupported(element_type))
      throw std::logic_error(fname + ": Found unsupported element type.");

    for (size_t i = 0; i < num_elems_in_block; ++i)
      std::getline(file, file_line);
  }

  // Read element data
  file.seekg(0);
  while (std::getline(file, file_line))
    if (element_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  if (not(iss >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag))
    throw std::logic_error(fname + ": Failed to read number of element entity blocks.");

  std::vector<std::tuple<size_t, int, unsigned int, std::vector<size_t>>> element_data;
  for (size_t n = 0; n < num_entity_blocks; ++n)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);
    int entity_dim = 0, entity_tag = 0, element_type = 0;
    size_t num_elems_in_block = 0;
    if (not(iss >> entity_dim >> entity_tag >> element_type >> num_elems_in_block))
    {
      throw std::logic_error(fname + ": Failed to read number of elements in block " +
                             std::to_string(n) + ".");
    }

    // Skip point type elements
    if (element_type == 15)
    {
      for (size_t i = 0; i < num_elems_in_block; ++i)
        std::getline(file, file_line);
      continue;
    }

    int physical_reg = -1;
    if (entity_dim == 1)
      physical_reg = curve_entities[entity_tag];
    else if (entity_dim == 2)
      physical_reg = surface_entities[entity_tag];
    else if (entity_dim == 3)
      physical_reg = volume_entities[entity_tag];
    if (physical_reg == -1)
    {
      throw std::logic_error(fname + ": Failed to map physical region for block " +
                             std::to_string(n) + ".");
    }

    int num_cell_nodes = 0;
    if (element_type == 1) // 2-node edge
      num_cell_nodes = 2;
    else if (element_type == 2) // 3-node triangle
      num_cell_nodes = 3;
    else if (element_type == 3 or element_type == 4) // 4-node quadrangle or tet
      num_cell_nodes = 4;
    else if (element_type == 5) // 8-node hexahedron
      num_cell_nodes = 8;

    for (size_t i = 0; i < num_elems_in_block; ++i)
    {
      std::getline(file, file_line);
      iss = std::istringstream(file_line);
      size_t element_tag = 0;
      iss >> element_tag;

      std::vector<size_t> node_tags(num_cell_nodes);
      for (int j = 0; j < num_cell_nodes; ++j)
        iss >> node_tags[j];

      element_data.emplace_back(element_tag, element_type, physical_reg, node_tags);
    }
  }

  for (const auto& [element_tag, element_type, physical_reg, node_tags] : element_data)
  {
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
    cell.block_id = physical_reg;
    std::vector<uint64_t> nodes(node_tags.size());
    for (size_t i = 0; i < node_tags.size(); ++i)
      nodes[i] = node_tags[i] - 1;
    cell.vertex_ids = nodes;

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
      lw_faces[0].vertex_ids = {v[0], v[2], v[1]}; // Base face
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
      throw std::runtime_error(fname + ": Unsupported cell type.");

  } // for elements

  file.close();

  // create boundary names and IDs
  unsigned int dimension = (mesh_is_2D) ? 2 : 3;
  mesh->SetDimension(dimension);
  for (auto& [id, e] : physical_names)
  {
    if (std::get<0>(e) == dimension - 1)
    {
      auto name = std::get<1>(e);
      mesh->AddBoundary(id, name);
    }
  }

  mesh->SetType(UNSTRUCTURED);
  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // remap boundary cells onto cell faces
  std::map<std::set<uint64_t>, unsigned int> bnd_cell_to_bnd_id_map;
  for (auto& bnd_cell : mesh->GetRawBoundaryCells())
  {
    std::set<uint64_t> key;
    for (auto& vid : bnd_cell->vertex_ids)
      key.insert(vid);
    bnd_cell_to_bnd_id_map[key] = bnd_cell->block_id;
  }
  auto& raw_cells = mesh->GetRawCells();
  for (auto& cell_ptr : raw_cells)
    for (auto& face : cell_ptr->faces)
      if (not face.has_neighbor)
      {
        std::set<uint64_t> key;
        for (auto& vid : face.vertex_ids)
          key.insert(vid);

        auto it = bnd_cell_to_bnd_id_map.find(key);
        if (it != bnd_cell_to_bnd_id_map.end())
          face.neighbor = it->second;
      }

  log.Log() << "Done processing " << options.file_name << ".\n"
            << "Number of nodes read: " << mesh->GetVertices().size() << "\n"
            << "Number of cells read: " << mesh->GetRawCells().size();

  return mesh;
}

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromGmshV41Binary(const UnpartitionedMesh::Options& options, int data_size)
{
  const std::string fname = "MeshIO::FromGmshV41Binary";
  const std::string node_section_name = "$Nodes";
  const std::string element_section_name = "$Elements";
  const std::string entities_section_name = "$Entities";
  const std::string physical_names_section_name = "$PhysicalNames";

  std::ifstream file(options.file_name, std::ios::binary);
  if (not file.is_open())
    throw std::runtime_error(fname + ": Failed to open file " + options.file_name);
  std::string line;
  std::istringstream iss;

  auto mesh = std::make_shared<UnpartitionedMesh>();
  log.Log() << "Making unpartitioned mesh from Gmsh file " << options.file_name
            << " (format v4.1 binary)";

  // Map from physical name ID to [dim, physical name]
  std::map<int, std::tuple<int, std::string>> physical_names;
  std::map<int, int> curve_entities, surface_entities, volume_entities;
  bool have_entities = false;
  bool mesh_is_2D = true;

  // Scan sections
  file.clear();
  file.seekg(0);
  while (std::getline(file, line))
  {
    if (physical_names_section_name == line)
    {
      if (NextIsAsciiNumeric(file))
      {
        std::getline(file, line);
        iss = std::istringstream(line);
        int num_physical_names = 0;
        if (not(iss >> num_physical_names))
          throw std::logic_error(fname + ": Failed to read number of physical names.");
        for (int n = 0; n < num_physical_names; n++)
        {
          std::getline(file, line);
          iss = std::istringstream(line);
          int dim = 0;
          int phys_name_id = 0;
          std::string name;
          if (not(iss >> dim >> phys_name_id >> name))
            throw std::logic_error(fname + ": Failed to read physical name.");
          name = name.substr(1, name.length() - 2);
          physical_names[phys_name_id] = {dim, name};
        }
      }
      else
      {
        const auto num_physical_names_bin = static_cast<size_t>(ReadUInt32(file));
        for (size_t n = 0; n < num_physical_names_bin; ++n)
        {
          const int dim = static_cast<int>(ReadUInt32(file));
          const int phys_name_id = static_cast<int>(ReadUInt32(file));
          std::string name = ReadBinaryString(file, data_size);
          physical_names[phys_name_id] = {dim, name};
        }
        ExpectEndSection(file, "$EndPhysicalNames");
      }
    }
    else if (entities_section_name == line)
    {
      have_entities = true;
      size_t num_points = 0, num_curves = 0, num_surfaces = 0, num_volumes = 0;
      if (NextIsAsciiNumeric(file))
      {
        std::getline(file, line);
        iss = std::istringstream(line);
        if (not(iss >> num_points >> num_curves >> num_surfaces >> num_volumes))
          throw std::logic_error(fname + ": Failed to read number of entities.");
        for (size_t i = 0; i < num_points; ++i)
          std::getline(file, line);

        auto read_entity = [&](std::map<int, int>& entity_map, size_t num_entities)
        {
          for (size_t i = 0; i < num_entities; ++i)
          {
            int entity_tag = 0, physical_tag = 0;
            size_t num_physical_tags = 0;
            double minx = 0.0, miny = 0.0, minz = 0.0, maxx = 0.0, maxy = 0.0, maxz = 0.0;
            std::getline(file, line);
            iss = std::istringstream(line);
            iss >> entity_tag >> minx >> miny >> minz >> maxx >> maxy >> maxz >>
              num_physical_tags >> physical_tag;
            entity_map[entity_tag] = physical_tag;
          }
        };

        read_entity(curve_entities, num_curves);
        read_entity(surface_entities, num_surfaces);
        read_entity(volume_entities, num_volumes);
      }
      else
      {
        num_points = ReadSize(file);
        num_curves = ReadSize(file);
        num_surfaces = ReadSize(file);
        num_volumes = ReadSize(file);
        const size_t max_reasonable = 10'000'000;
        if (num_points > max_reasonable || num_curves > max_reasonable ||
            num_surfaces > max_reasonable || num_volumes > max_reasonable)
          throw std::logic_error(fname + ": $Entities binary counts are invalid.");

        auto read_entity_bin = [&](std::map<int, int>& entity_map, size_t num_entities)
        {
          for (size_t i = 0; i < num_entities; ++i)
          {
            const int entity_tag = static_cast<int>(ReadUInt32(file));
            ReadReal(file, data_size);
            ReadReal(file, data_size);
            ReadReal(file, data_size);
            ReadReal(file, data_size);
            ReadReal(file, data_size);
            ReadReal(file, data_size);
            const size_t num_physical_tags = ReadSize(file);
            int physical_tag = 0;
            if (num_physical_tags > 0)
              physical_tag = static_cast<int>(ReadUInt32(file));
            for (size_t t = 1; t < num_physical_tags; ++t)
              ReadUInt32(file);
            const size_t num_bounding_tags = ReadSize(file);
            for (size_t t = 0; t < num_bounding_tags; ++t)
              ReadUInt32(file);
            entity_map[entity_tag] = physical_tag;
          }
        };

        for (size_t i = 0; i < num_points; ++i)
        {
          (void)ReadUInt32(file);
          (void)ReadReal(file, data_size);
          (void)ReadReal(file, data_size);
          (void)ReadReal(file, data_size);
          const size_t num_physical_tags = ReadSize(file);
          const size_t max_tags = 1'000'000;
          if (num_physical_tags > max_tags)
            throw std::logic_error(fname + ": $Entities point has too many physical tags.");
          if (num_physical_tags > 0)
            ReadUInt32(file);
          for (size_t t = 1; t < num_physical_tags; ++t)
            ReadUInt32(file);
        }
        read_entity_bin(curve_entities, num_curves);
        read_entity_bin(surface_entities, num_surfaces);
        read_entity_bin(volume_entities, num_volumes);
        ExpectEndSection(file, "$EndEntities");
      }
    }
    else if (node_section_name == line)
    {
      // Read node data (binary)
      const uint64_t num_entity_blocks = ReadSize(file);
      const uint64_t num_nodes = ReadSize(file);
      ReadSize(file); // min_node_tag
      ReadSize(file); // max_node_tag

      auto& vertices = mesh->GetVertices();
      vertices.clear();
      vertices.resize(num_nodes);
      for (uint64_t n = 0; n < num_entity_blocks; ++n)
      {
        const uint64_t entity_dim = ReadUInt32(file);
        ReadUInt32(file); // entity_tag
        const uint64_t parametric = ReadUInt32(file);
        const uint64_t num_nodes_in_block = ReadSize(file);

        std::vector<uint64_t> node_tags(num_nodes_in_block);
        for (uint64_t i = 0; i < num_nodes_in_block; ++i)
          node_tags[i] = ReadSize(file);

        for (uint64_t i = 0; i < num_nodes_in_block; ++i)
        {
          const double x = ReadReal(file, data_size);
          const double y = ReadReal(file, data_size);
          const double z = ReadReal(file, data_size);
          vertices[node_tags[i] - 1] = {x, y, z};
          if (parametric != 0)
            for (uint64_t p = 0; p < entity_dim; ++p)
              ReadReal(file, data_size);
        }
      }

      ExpectEndSection(file, "$EndNodes");
    }
    else if (element_section_name == line)
    {
      if (not have_entities)
        throw std::logic_error(fname + ": Missing $Entities section for element mapping.");

      const uint64_t num_entity_blocks = ReadSize(file);
      const uint64_t num_elements = ReadSize(file);
      ReadSize(file); // min_element_tag
      ReadSize(file); // max_element_tag

      (void)num_elements;
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

      bool has_3d_elements = false;
      const auto elements_data_pos = file.tellg();
      for (uint64_t n = 0; n < num_entity_blocks; ++n)
      {
        const int entity_dim = static_cast<int>(ReadUInt32(file));
        const int entity_tag = static_cast<int>(ReadUInt32(file));
        const int element_type = static_cast<int>(ReadUInt32(file));
        const uint64_t num_elems_in_block = ReadSize(file);

        (void)entity_dim;
        (void)entity_tag;

        if (element_type == 15)
        {
          const uint64_t bytes_per_elem = sizeof(uint64_t) * 2;
          file.seekg(static_cast<std::streamoff>(bytes_per_elem * num_elems_in_block),
                     std::ios::cur);
          continue;
        }

        if (IsElementType3D(element_type))
          has_3d_elements = true;

        int num_cell_nodes = 0;
        if (element_type == 1)
          num_cell_nodes = 2;
        else if (element_type == 2)
          num_cell_nodes = 3;
        else if (element_type == 3 || element_type == 4)
          num_cell_nodes = 4;
        else if (element_type == 5)
          num_cell_nodes = 8;
        else if (element_type == 6)
          num_cell_nodes = 6;
        else if (element_type == 7)
          num_cell_nodes = 5;

        const uint64_t bytes_per_elem = sizeof(uint64_t) * (1 + num_cell_nodes);
        file.seekg(static_cast<std::streamoff>(bytes_per_elem * num_elems_in_block), std::ios::cur);
      }

      file.clear();
      file.seekg(elements_data_pos);
      mesh_is_2D = !has_3d_elements;

      for (uint64_t n = 0; n < num_entity_blocks; ++n)
      {
        const int entity_dim = static_cast<int>(ReadUInt32(file));
        const int entity_tag = static_cast<int>(ReadUInt32(file));
        const int element_type = static_cast<int>(ReadUInt32(file));
        const uint64_t num_elems_in_block = ReadSize(file);

        if (element_type == 15)
        {
          for (uint64_t e = 0; e < num_elems_in_block; ++e)
          {
            ReadSize(file); // element_tag
            ReadSize(file); // node_tag
          }
          continue;
        }

        if (not IsElementSupported(element_type))
          throw std::logic_error(fname + ": Found unsupported element type.");

        int physical_reg = -1;
        if (entity_dim == 1)
          physical_reg = curve_entities[entity_tag];
        else if (entity_dim == 2)
          physical_reg = surface_entities[entity_tag];
        else if (entity_dim == 3)
          physical_reg = volume_entities[entity_tag];
        if (physical_reg == -1)
        {
          throw std::logic_error(fname + ": Failed to map physical region for block " +
                                 std::to_string(n) + ".");
        }

        int num_cell_nodes = 0;
        if (element_type == 1)
          num_cell_nodes = 2;
        else if (element_type == 2)
          num_cell_nodes = 3;
        else if (element_type == 3 or element_type == 4)
          num_cell_nodes = 4;
        else if (element_type == 5)
          num_cell_nodes = 8;

        for (uint64_t i = 0; i < num_elems_in_block; ++i)
        {
          ReadSize(file); // element_tag

          std::vector<size_t> node_tags(num_cell_nodes);
          for (int j = 0; j < num_cell_nodes; ++j)
            node_tags[j] = static_cast<size_t>(ReadSize(file));

          auto& raw_boundary_cells = mesh->GetRawBoundaryCells();
          auto& raw_cells = mesh->GetRawCells();

          std::shared_ptr<UnpartitionedMesh::LightWeightCell> raw_cell;
          if (mesh_is_2D)
          {
            if (IsElementType1D(element_type))
            {
              raw_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::SLAB,
                                                                              CellType::SLAB);
              raw_boundary_cells.push_back(raw_cell);
            }
            else if (IsElementType2D(element_type))
            {
              raw_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(
                CellType::POLYGON, CellTypeFromMSHTypeID(element_type));
              raw_cells.push_back(raw_cell);
            }
          }
          else
          {
            if (IsElementType2D(element_type))
            {
              raw_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(
                CellType::POLYGON, CellTypeFromMSHTypeID(element_type));
              raw_boundary_cells.push_back(raw_cell);
            }
            else if (IsElementType3D(element_type))
            {
              raw_cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(
                CellType::POLYHEDRON, CellTypeFromMSHTypeID(element_type));
              raw_cells.push_back(raw_cell);
            }
          }

          if (raw_cell == nullptr)
            continue;

          auto& cell = *raw_cell;
          cell.block_id = physical_reg;
          std::vector<uint64_t> nodes(node_tags.size());
          for (size_t k = 0; k < node_tags.size(); ++k)
            nodes[k] = node_tags[k] - 1;
          cell.vertex_ids = nodes;

          if (element_type == 1)
          {
            UnpartitionedMesh::LightWeightFace face0;
            UnpartitionedMesh::LightWeightFace face1;

            face0.vertex_ids = {cell.vertex_ids.at(0)};
            face1.vertex_ids = {cell.vertex_ids.at(1)};

            cell.faces.push_back(face0);
            cell.faces.push_back(face1);
          }
          else if (element_type == 2 or element_type == 3)
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
          else if (element_type == 4)
          {
            auto& v = cell.vertex_ids;
            std::vector<UnpartitionedMesh::LightWeightFace> lw_faces(4);
            lw_faces[0].vertex_ids = {v[0], v[2], v[1]}; // Base face
            lw_faces[1].vertex_ids = {v[0], v[3], v[2]};
            lw_faces[2].vertex_ids = {v[3], v[1], v[2]};
            lw_faces[3].vertex_ids = {v[3], v[0], v[1]};

            for (auto& lw_face : lw_faces)
              cell.faces.push_back(lw_face);
          }
          else if (element_type == 5)
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
          else if (element_type == 6 or element_type == 7)
          {
            throw std::logic_error(fname + ": Polyhedral cells are unsupported for binary reader.");
          }
        }
      }

      ExpectEndSection(file, "$EndElements");
    }
  }

  // create boundary names and IDs
  const unsigned int dimension = (mesh_is_2D) ? 2 : 3;
  mesh->SetDimension(dimension);
  for (const auto& [id, data] : physical_names)
  {
    if (std::cmp_equal(std::get<0>(data), dimension - 1))
    {
      const auto& name = std::get<1>(data);
      mesh->AddBoundary(id, name);
    }
  }

  mesh->SetType(UNSTRUCTURED);
  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // remap boundary cells onto cell faces
  std::map<std::set<uint64_t>, unsigned int> bnd_cell_to_bnd_id_map;
  for (auto& bnd_cell : mesh->GetRawBoundaryCells())
  {
    std::set<uint64_t> key;
    for (auto& vid : bnd_cell->vertex_ids)
      key.insert(vid);
    bnd_cell_to_bnd_id_map[key] = bnd_cell->block_id;
  }
  for (auto& cell_ptr : mesh->GetRawCells())
    for (auto& face : cell_ptr->faces)
      if (not face.has_neighbor)
      {
        std::set<uint64_t> key;
        for (auto& vid : face.vertex_ids)
          key.insert(vid);

        auto it = bnd_cell_to_bnd_id_map.find(key);
        if (it != bnd_cell_to_bnd_id_map.end())
          face.neighbor = it->second;
      }

  log.Log() << "Done processing " << options.file_name << ".\n"
            << "Number of nodes read: " << mesh->GetVertices().size() << "\n"
            << "Number of cells read: " << mesh->GetRawCells().size();

  return mesh;
}

} // namespace opensn
