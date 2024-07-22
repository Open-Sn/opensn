// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/io/mesh_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include <fstream>
#include <algorithm>

namespace opensn
{

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
  std::string file_line;
  std::string delimiter = " ";
  int material_id = -1;
  while (std::getline(file, file_line))
  {
    // Get the first word
    size_t beg_of_word = file_line.find_first_not_of(delimiter);
    size_t end_of_word = file_line.find(delimiter, beg_of_word - beg_of_word);
    std::string first_word = file_line.substr(beg_of_word, end_of_word);
    std::string sub_word;

    if (first_word == "o")
    {
      beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
      end_of_word = file_line.find(delimiter, beg_of_word);
      sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

      std::string block_name = sub_word;
      block_data.push_back({block_name, {}});
    }

    if (first_word == "usemtl")
    {
      log.Log0Verbose1() << "New material at cell count: " << block_data.back().cells.size();
      ++material_id;
    }

    // Keyword "v" for Vertex
    if (first_word == "v")
    {
      Vector3 newVertex;
      for (int k = 1; k <= 3; k++)
      {
        // Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

        // Convert word to number
        try
        {
          double numValue = std::stod(sub_word);

          if (k == 1)
            newVertex.x = numValue;
          else if (k == 2)
            newVertex.y = numValue;
          else if (k == 3)
            newVertex.z = numValue;
        }

        // Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          log.Log0Warning() << "Failed to convert vertex in line " << file_line << std::endl;
        }

        // Stop word extraction on line-end
        if (end_of_word == std::string::npos)
          break;
      }
      file_vertices.push_back(newVertex);
    } // if (first_word == "v")

    // Keyword "f" for face
    if (first_word == "f")
    {
      size_t number_of_verts = std::count(file_line.begin(), file_line.end(), '/') / 2;

      CellType sub_type = CellType::POLYGON;
      if (number_of_verts == 3)
        sub_type = CellType::TRIANGLE;
      else if (number_of_verts == 4)
        sub_type = CellType::QUADRILATERAL;

      auto cell = std::make_shared<UnpartitionedMesh::LightWeightCell>(CellType::POLYGON, sub_type);
      cell->material_id = material_id;

      // Populate vertex-ids
      for (size_t k = 1; k <= number_of_verts; k++)
      {
        // Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

        // Extract locations of hyphens
        size_t first_dash = sub_word.find('/');
        size_t last_dash = sub_word.find_last_of('/');

        // Extract the words ass. vertex and normal
        std::string vert_word = sub_word.substr(0, first_dash - 0);
        std::string norm_word = sub_word.substr(last_dash + 1, sub_word.length() - last_dash - 1);

        // Convert word to number (Vertex)
        try
        {
          int numValue = std::stoi(vert_word);
          cell->vertex_ids.push_back(numValue - 1);
        }
        catch (const std::invalid_argument& ia)
        {
          log.Log0Warning() << "Failed converting word to number in line " << file_line
                            << std::endl;
        }

        // Stop word extraction on line-end
        if (end_of_word == std::string::npos)
        {
          break;
        }
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

    // Keyword "l" for edge
    if (first_word == "l")
    {
      std::pair<uint64_t, uint64_t> edge;
      for (int k = 1; k <= 2; ++k)
      {
        // Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

        // Convert word to number
        try
        {
          int vertex_id = std::stoi(sub_word);
          if (k == 1)
            edge.first = vertex_id - 1;
          if (k == 2)
            edge.second = vertex_id - 1;
        }

        // Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          log.Log0Warning() << "Failed to convert text to integer in line " << file_line
                            << std::endl;
        }
      } // for k

      if (block_data.empty())
        throw std::logic_error(fname + ": Could not add edge to block-data. "
                                       "This normally indicates that the file does not have the "
                                       "\"o Object Name\" entry.");

      block_data.back().edges.push_back(edge);
    } // if (first_word == "l")
  }
  file.close();
  log.Log0Verbose0() << "Max material id: " << material_id;

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
    std::vector<size_t> vertex_map;
    vertex_map.reserve(file_vertices.size());
    for (size_t m = 0; m < file_vertices.size(); ++m)
      vertex_map.push_back(m);

    // Build set of cell vertices
    std::set<size_t> cell_vertex_id_set;
    for (const auto& cell_ptr : block_data.at(main_block_id).cells)
      for (size_t vid : cell_ptr->vertex_ids)
        cell_vertex_id_set.insert(vid);

    // Make cell_vertices and edit map
    {
      size_t new_id = 0;
      for (size_t vid : cell_vertex_id_set)
      {
        cell_vertices.push_back(file_vertices[vid]);
        vertex_map[vid] = new_id;
        ++new_id;
      }
    }

    // Build set of bndry vertices
    std::set<size_t> bndry_vertex_id_set;
    for (size_t block_id : bndry_block_ids)
      for (const auto& edge : block_data[block_id].edges)
      {
        bndry_vertex_id_set.insert(edge.first);
        bndry_vertex_id_set.insert(edge.second);
      }

    // Find a match for each boundary vertex and
    // place it in the map
    for (size_t bvid : bndry_vertex_id_set)
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
        edge.first = static_cast<int>(vertex_map[edge.first]);
        edge.second = static_cast<int>(vertex_map[edge.second]);
      }
  }
  mesh->Vertices() = cell_vertices;
  mesh->RawCells() = block_data[main_block_id].cells;

  // Always do this
  mesh->SetDimension(2);
  mesh->SetType(UNSTRUCTURED);

  mesh->ComputeCentroids();
  mesh->CheckQuality();
  mesh->BuildMeshConnectivity();

  // Set boundary ids
  if (bndry_block_ids.empty())
  {
    mesh->AddBoundary(0, "Default Boundary");
  }
  else
  {
    std::vector<UnpartitionedMesh::LightWeightFace*> bndry_faces;
    for (auto& cell_ptr : mesh->RawCells())
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
        std::set<size_t> edge_vert_id_set({edge.first, edge.second});

        for (auto& face_ptr : bndry_faces)
        {
          const auto& vert_ids = face_ptr->vertex_ids;
          std::set<size_t> face_vert_id_set(vert_ids.begin(), vert_ids.end());

          if (face_vert_id_set == edge_vert_id_set)
          {
            face_ptr->neighbor = bndry_id;
            ++num_faces_boundarified;
            break;
          }
        } // for face
      }   // for edge

      log.Log() << "UnpartitionedMesh: assigned " << num_faces_boundarified
                << " faces to boundary id " << bndry_id << " with name " << block_data[bid].name;

      mesh->AddBoundary(bndry_id, block_data[bid].name);

      ++bndry_id;
    } // for boundary block
  }

  return mesh;
}

} // namespace opensn
