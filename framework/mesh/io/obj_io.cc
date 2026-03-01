// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include "framework/utils/utils.h"
#include <fstream>
#include <algorithm>
#include <optional>
#include <stdexcept>

namespace opensn
{

namespace
{

/// extract first part from the `vertex description`
std::string
ExtractFirstPart(std::string_view input)
{
  const auto pos = input.find('/');
  return std::string(input.substr(0, pos));
};

int
ConvertToInt(const std::string& str, const std::string_view file_name, int line_no)
{
  try
  {
    return std::stoi(str);
  }
  catch (const std::invalid_argument& ia)
  {
    std::ostringstream oss;
    oss << "Failed to convert string to a number in " << file_name << ", line " << line_no
        << std::endl;
    throw std::runtime_error(oss.str());
  }
}

double
ConvertToDouble(const std::string& str, const std::string_view file_name, int line_no)
{
  try
  {
    return std::stod(str);
  }
  catch (const std::invalid_argument& ia)
  {
    std::ostringstream oss;
    oss << "Failed to convert string to a floating-point number in " << file_name << ", line "
        << line_no << std::endl;
    throw std::runtime_error(oss.str());
  }
}

} // namespace

std::shared_ptr<UnpartitionedMesh>
MeshIO::FromOBJ(const UnpartitionedMesh::Options& options)
{
  const std::string fname = "MeshIO::ReadFromWavefrontOBJ";

  // Opening the file
  std::ifstream file;
  file.open(options.file_name);
  if (not file.is_open())
  {
    std::ostringstream oss;
    oss << fname << ": Failed to open file: " << options.file_name;
    throw std::runtime_error(oss.str());
  }

  log.Log() << "Making Unpartitioned mesh from wavefront file " << options.file_name;

  std::shared_ptr<UnpartitionedMesh> mesh = std::make_shared<UnpartitionedMesh>();

  struct BlockData
  {
    std::string name;
    std::vector<std::shared_ptr<UnpartitionedMesh::LightWeightCell>> cells;
    std::vector<std::pair<uint64_t, uint64_t>> edges;
  };

  std::vector<BlockData> block_data;
  std::vector<Vector3> file_vertices;

  // Reading every line
  int line_no = 0;
  std::string file_line;
  std::optional<unsigned int> material_id;
  while (std::getline(file, file_line))
  {
    file_line = StringTrim(file_line);
    line_no++;
    auto parts = StringSplit(file_line, " ");
    if (parts.empty())
      continue;

    const auto& first_word = parts[0];
    if (first_word == "o")
    {
      if (parts.size() < 2)
        throw std::runtime_error("Expected block name, but got malformed line");
      auto block_name = parts[1];
      block_data.push_back({block_name, {}});
    }
    else if (first_word == "usemtl")
    {
      log.Log0Verbose1() << "New material at cell count: " << block_data.back().cells.size();

      if (not material_id.has_value())
        material_id = 0;
      else
        ++material_id.value();
    }
    else if (first_word == "v")
    {
      // Keyword "v" for Vertex
      if (parts.size() < 4)
        throw std::runtime_error("Expected line with vertex, but got malformed line");
      Vector3 new_vertex;
      for (int k = 1; k <= 3; ++k)
      {
        auto num_value = ConvertToDouble(parts[k], options.file_name, line_no);
        if (k == 1)
          new_vertex.x = num_value;
        else if (k == 2)
          new_vertex.y = num_value;
        else if (k == 3)
          new_vertex.z = num_value;
      }
      file_vertices.push_back(new_vertex);
    } // if (first_word == "v")
    else if (first_word == "f")
    {
      // Keyword "f" for face
      auto number_of_verts = parts.size() - 1;

      CellType sub_type = CellType::POLYGON;
      if (number_of_verts == 3)
        sub_type = CellType::TRIANGLE;
      else if (number_of_verts == 4)
        sub_type = CellType::QUADRILATERAL;

      auto cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYGON, sub_type);
      cell->block_id = material_id.value_or(std::numeric_limits<unsigned int>::max());

      // Populate vertex-ids
      for (size_t k = 1; k <= number_of_verts; ++k)
      {
        // Extract the vertex ID
        auto vert_word = ExtractFirstPart(parts[k]);
        auto num_value = ConvertToInt(vert_word, options.file_name, line_no);
        cell->vertex_ids.push_back(num_value - 1);
      }

      // Build faces
      const size_t num_verts = cell->vertex_ids.size();
      for (size_t v = 0; v < num_verts; ++v)
      {
        UnpartitionedMesh::LightWeightFace face;

        face.vertex_ids.resize(2);
        face.vertex_ids[0] = cell->vertex_ids[v];
        face.vertex_ids[1] = (v < (num_verts - 1)) ? cell->vertex_ids[v + 1] : cell->vertex_ids[0];

        cell->faces.push_back(std::move(face));
      }

      if (block_data.empty())
        throw std::logic_error(fname + ": Could not add cell to block-data. "
                                       "This normally indicates that the file does not have the "
                                       "\"o Object Name\" entry.");

      block_data.back().cells.push_back(cell);
    } // if (first_word == "f")
    else if (first_word == "l")
    {
      // Keyword "l" for edge
      if (parts.size() < 3)
        throw std::runtime_error("Expected line with vertex, but got malformed line");
      std::pair<uint64_t, uint64_t> edge;
      for (int k = 1; k <= 2; ++k)
      {
        auto vertex_id = ConvertToInt(parts[k], options.file_name, line_no);
        if (k == 1)
          edge.first = vertex_id - 1;
        if (k == 2)
          edge.second = vertex_id - 1;
      } // for k

      if (block_data.empty())
        throw std::logic_error(fname + ": Could not add edge to block-data. "
                                       "This normally indicates that the file does not have the "
                                       "\"o Object Name\" entry.");

      block_data.back().edges.push_back(edge);
    } // if (first_word == "l")
  }
  file.close();
  if (material_id.has_value())
    log.Log0Verbose0() << "Max material id: " << material_id.value();

  // Error checks
  for (const auto& block : block_data)
    for (const auto& cell : block.cells)
      for (const auto vid : cell->vertex_ids)
      {
        OpenSnLogicalErrorIf(
          vid >= file_vertices.size(),
          "Cell vertex id " + std::to_string(vid) +
            " not in list of vertices read (size=" + std::to_string(file_vertices.size()) + ").");
      }

  // Filter blocks
  std::vector<size_t> bndry_block_ids;
  size_t num_cell_blocks = 0;
  size_t main_block_id = 0;
  for (size_t block_id = 0; block_id < block_data.size(); ++block_id)
  {
    if (not block_data[block_id].edges.empty())
      bndry_block_ids.push_back(block_id);
    if (not block_data[block_id].cells.empty())
    {
      ++num_cell_blocks;
      main_block_id = block_id;
    }
  } // for block_id

  if (num_cell_blocks != 1)
    throw std::logic_error(fname +
                           ": More than one cell-block has been read "
                           "from the file. Only a single face-containing object is supported. "
                           "If you exported this mesh from blender, be sure to export "
                           "\"selection only\"");

  // Process blocks
  std::vector<Vector3> cell_vertices;
  {
    // Initial map is straight
    std::vector<uint64_t> vertex_map;
    vertex_map.reserve(file_vertices.size());
    for (size_t m = 0; m < file_vertices.size(); ++m)
      vertex_map.push_back(m);

    // Build set of cell vertices
    std::set<uint64_t> cell_vertex_id_set;
    for (const auto& cell_ptr : block_data.at(main_block_id).cells)
      for (auto vid : cell_ptr->vertex_ids)
        cell_vertex_id_set.insert(vid);

    // Make cell_vertices and edit map
    {
      uint64_t new_id = 0;
      for (auto vid : cell_vertex_id_set)
      {
        cell_vertices.push_back(file_vertices[vid]);
        vertex_map[vid] = new_id;
        ++new_id;
      }
    }

    // Build set of bndry vertices
    std::set<uint64_t> bndry_vertex_id_set;
    for (size_t block_id : bndry_block_ids)
      for (const auto& edge : block_data[block_id].edges)
      {
        bndry_vertex_id_set.insert(edge.first);
        bndry_vertex_id_set.insert(edge.second);
      }

    // Find a match for each boundary vertex and
    // place it in the map
    for (auto bvid : bndry_vertex_id_set)
    {
      const auto& bndry_vertex = file_vertices[bvid];

      bool match_found = false;
      for (size_t cvid = 0; cvid < cell_vertices.size(); ++cvid)
      {
        const auto& cell_vertex = cell_vertices[cvid];

        if ((bndry_vertex - cell_vertex).NormSquare() < 1.0e-12)
        {
          vertex_map[bvid] = cvid;
          match_found = true;
          break;
        }
      } // for cvid
      if (not match_found)
        throw std::logic_error(fname +
                               ": Failed to map a boundary vertex to"
                               "any cell vertex. Check that the edges are conformal with the "
                               "object containing the faces.");
    } // for bvid

    // Change cell and face vertex ids to cell vertex ids
    // using vertex map
    for (auto& cell_ptr : block_data.at(main_block_id).cells)
    {
      for (uint64_t& vid : cell_ptr->vertex_ids)
        vid = vertex_map[vid];

      for (auto& face : cell_ptr->faces)
        for (uint64_t& vid : face.vertex_ids)
          vid = vertex_map[vid];
    }

    // Change edge vertex ids to cell vertex ids using
    // the vertex map
    for (size_t block_id : bndry_block_ids)
      for (auto& edge : block_data[block_id].edges)
      {
        edge.first = vertex_map[edge.first];
        edge.second = vertex_map[edge.second];
      }
  }
  mesh->GetVertices() = cell_vertices;
  mesh->GetRawCells() = block_data[main_block_id].cells;

  // Always do this
  mesh->SetDimension(2);
  mesh->SetType(UNSTRUCTURED);

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // Set boundary ids
  if (not bndry_block_ids.empty())
  {
    std::vector<UnpartitionedMesh::LightWeightFace*> bndry_faces;
    for (auto& cell_ptr : mesh->GetRawCells())
      for (auto& face : cell_ptr->faces)
        if (not face.has_neighbor)
          bndry_faces.push_back(&face);

    size_t bndry_id = 0;
    for (size_t bid : bndry_block_ids)
    {
      const auto& bndry_edges = block_data[bid].edges;

      size_t num_faces_boundarified = 0;
      for (const auto& edge : bndry_edges)
      {
        std::set<uint64_t> edge_vert_id_set({edge.first, edge.second});

        for (auto& face_ptr : bndry_faces)
        {
          const auto& vert_ids = face_ptr->vertex_ids;
          std::set<uint64_t> face_vert_id_set(vert_ids.begin(), vert_ids.end());

          if (face_vert_id_set == edge_vert_id_set)
          {
            face_ptr->neighbor = bndry_id;
            ++num_faces_boundarified;
            break;
          }
        } // for face
      } // for edge

      log.Log() << "UnpartitionedMesh: assigned " << num_faces_boundarified
                << " faces to boundary id " << bndry_id << " with name " << block_data[bid].name;

      mesh->AddBoundary(bndry_id, block_data[bid].name);

      ++bndry_id;
    } // for boundary block
  }

  return mesh;
}

} // namespace opensn
