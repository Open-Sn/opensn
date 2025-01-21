// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <fstream>

namespace opensn
{

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromGmshV41(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::FromGmshV41";
  const std::string format_section_name = "$MeshFormat";
  const std::string node_section_name = "$Nodes";
  const std::string element_section_name = "$Elements";
  const std::string entities_section_name = "$Entities";
  const std::string partitioned_entities_section_name = "$PartitionedEntities";

  // Open file
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
  if (!(iss >> gmsh_version) || gmsh_version != 4.1)
    throw std::logic_error(fname + ": Only Gmsh version 4.1 format is supported.");

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

  // Read $Entities section
  file.seekg(0);
  while (std::getline(file, file_line))
    if (entities_section_name == file_line)
      break;

  std::getline(file, file_line);
  iss = std::istringstream(file_line);
  size_t num_points, num_curves, num_surfaces, num_volumes;
  if (!(iss >> num_points >> num_curves >> num_surfaces >> num_volumes))
    throw std::logic_error(fname + ": Failed to read number of entitities.");

  // Skip point entities
  for (int i = 0; i < num_points; ++i)
    std::getline(file, file_line);

  // Read curve, surface, and volume entities
  auto read_entity = [&](std::map<int, int>& entity_map, size_t num_entities)
  {
    for (size_t i = 0; i < num_entities; ++i)
    {
      int entity_tag, physical_tag;
      size_t num_physical_tags;
      double minx, miny, minz, maxx, maxy, maxz;
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
  size_t num_entity_blocks, num_nodes, min_node_tag, max_node_tag;
  if (not(iss >> num_entity_blocks >> num_nodes >> min_node_tag >> max_node_tag))
    throw std::logic_error(fname + ": Failed to read the number of node entity blocks.");

  auto& vertices = mesh->GetVertices();
  vertices.clear();
  vertices.resize(num_nodes);

  for (int n = 0; n < num_entity_blocks; ++n)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);
    int entity_dim, entity_tag, parametric;
    size_t num_nodes_in_block;
    if (not(iss >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block))
    {
      throw std::logic_error(fname + ": Failed to read number of nodes in block " +
                             std::to_string(n) + ".");
    }

    std::vector<size_t> node_tags(num_nodes_in_block);
    for (int i = 0; i < num_nodes_in_block; ++i)
    {
      std::getline(file, file_line);
      iss = std::istringstream(file_line);
      iss >> node_tags[i];
    }

    for (int i = 0; i < num_nodes_in_block; ++i)
    {
      std::getline(file, file_line);
      iss = std::istringstream(file_line);
      double x, y, z;
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
  size_t num_elements, min_element_tag, max_element_tag;
  if (not(iss >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag))
    throw std::logic_error(fname + ": Failed to read number of element entity blocks.");

  for (int n = 0; n < num_entity_blocks; ++n)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);
    int entity_dim, entity_tag, element_type;
    size_t num_elems_in_block;
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
      log.Log() << "Mesh identified as 3D.";
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

  std::vector<std::tuple<size_t, int, int, std::vector<size_t>>> element_data;
  for (int n = 0; n < num_entity_blocks; ++n)
  {
    std::getline(file, file_line);
    iss = std::istringstream(file_line);
    int entity_dim, entity_tag, element_type;
    size_t num_elems_in_block;
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
      size_t element_tag;
      iss >> element_tag;

      std::vector<size_t> node_tags(num_cell_nodes);
      for (int j = 0; j < num_cell_nodes; ++j)
        iss >> node_tags[j];

      element_data.emplace_back(
        std::make_tuple(element_tag, element_type, physical_reg, node_tags));
    }
  }

  if (element_data.size() != num_elements)
    throw std::logic_error(fname + ": Missing element data.");

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

  // Remap block-ids
  std::set<int> block_ids_set_as_read;
  std::map<int, int> block_mapping;

  for (auto& cell : mesh->GetRawCells())
    block_ids_set_as_read.insert(cell->block_id);

  std::set<int> boundary_ids_set_as_read;
  std::map<int, int> boundary_mapping;

  for (auto& cell : mesh->GetRawBoundaryCells())
    boundary_ids_set_as_read.insert(cell->block_id);

  {
    int m = 0;
    for (const auto& blk_id : block_ids_set_as_read)
      block_mapping.insert(std::make_pair(blk_id, m++));

    int b = 0;
    for (const auto& bndry_id : boundary_ids_set_as_read)
      boundary_mapping.insert(std::make_pair(bndry_id, b++));
  }

  for (auto& cell : mesh->GetRawCells())
    cell->block_id = block_mapping[cell->block_id];

  for (auto& cell : mesh->GetRawBoundaryCells())
    cell->block_id = boundary_mapping[cell->block_id];

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
