// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <algorithm>
#include <iostream>
#include <fstream>

namespace opensn
{

InputParameters
SurfaceMesh::GetInputParameters()
{
  InputParameters params = Object::GetInputParameters();
  return params;
}

SurfaceMesh::SurfaceMesh(const InputParameters& params) : Object(params)
{
}

SurfaceMesh::SurfaceMesh() = default;

SurfaceMesh::~SurfaceMesh()
{
  poly_faces_.clear();
}

std::ostream&
operator<<(std::ostream& os, SurfaceMesh& that)
{
  std::vector<Face>::const_iterator curface;
  for (curface = that.GetTriangles().begin(); curface != that.GetTriangles().end(); ++curface)
  {
    long index = std::distance(that.GetTriangles().begin(), curface);
    os << "Face " << index << " v:";
    os << curface->v_index[0] << "->";
    os << curface->v_index[1] << "->";
    os << curface->v_index[2] << "  ";

    os << "e:";
    for (int e = 0; e < 3; ++e)
    {
      os << "[" << curface->e_index[e][0] << ",";
      os << curface->e_index[e][1] << ",";
      os << curface->e_index[e][2] << ",";
      os << curface->e_index[e][3] << "]";
    }
    os << std::endl;
  }

  return os;
}

void
SurfaceMesh::UpdateInternalConnectivity()
{
  std::vector<std::vector<size_t>> vertex_subscriptions;

  // Initialize vertex subscription
  size_t num_verts = vertices_.size();
  vertex_subscriptions.resize(num_verts);

  for (auto& vert_sub : vertex_subscriptions)
    vert_sub.reserve(5);

  // Loop over cells and determine which cells subscribe to a vertex
  // Triangles
  size_t num_tri_faces = faces_.size();
  for (size_t tf = 0; tf < num_tri_faces; ++tf)
  {
    auto& try_face = faces_[tf];
    for (int v = 0; v < 3; ++v)
      vertex_subscriptions[v].push_back(tf);
  }
  // Polygons
  size_t num_poly_faces = poly_faces_.size();
  for (size_t pf = 0; pf < num_poly_faces; ++pf)
  {
    auto poly_face = poly_faces_[pf];
    for (auto v : poly_face->v_indices)
      vertex_subscriptions[v].push_back(pf);
  }

  // Loop over cells and determine connectivity
  // Triangles
  for (auto curFace : faces_)
  {
    for (int e = 0; e < 3; ++e)
    {
      auto& curface_edge = curFace.e_index[e];
      int vi = curface_edge[0];
      int vf = curface_edge[1];

      // Search cells subscribing to vi
      for (auto ofi : vertex_subscriptions[vi])
      {
        auto& other_face = faces_[ofi];

        for (size_t e2 = 0; e2 < 3; ++e2)
        {
          if ((curface_edge[0] == other_face.e_index[e2][1]) and
              (curface_edge[1] == other_face.e_index[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for ofi

      // Search cells subscribing to vf
      for (auto ofi : vertex_subscriptions[vf])
      {
        auto& other_face = faces_[ofi];

        for (size_t e2 = 0; e2 < 3; ++e2)
        {
          if ((curface_edge[0] == other_face.e_index[e2][1]) and
              (curface_edge[1] == other_face.e_index[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for ofi
    }     // for current face edges
  }       // for faces

  // Loop over cells and determine connectivity
  // Polygons
  for (const auto& curFace : poly_faces_)
  {
    for (auto& curface_edge : curFace->edges)
    {
      int vi = curface_edge[0];
      int vf = curface_edge[1];

      // Search cells subscribing to vi
      for (auto ofi : vertex_subscriptions[vi])
      {
        auto other_face = poly_faces_[ofi];

        for (size_t e2 = 0; e2 < other_face->edges.size(); ++e2)
        {
          if ((curface_edge[0] == other_face->edges[e2][1]) and
              (curface_edge[1] == other_face->edges[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for ofi

      // Search cells subscribing to vf
      for (auto ofi : vertex_subscriptions[vf])
      {
        auto other_face = poly_faces_[ofi];

        for (size_t e2 = 0; e2 < other_face->edges.size(); ++e2)
        {
          if ((curface_edge[0] == other_face->edges[e2][1]) and
              (curface_edge[1] == other_face->edges[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for other faces
    }     // for current face edges
  }       // for faces
}

int
SurfaceMesh::ImportFromOBJFile(const std::string& fileName, bool as_poly, const Vector3& transform)
{
  const std::string fname = "SurfaceMesh::ImportFromOBJFile";

  // Opening the file
  std::ifstream file;
  file.open(fileName);
  if (not file.is_open())
  {
    std::ostringstream oss;
    oss << fname << ": Failed to open file: " << fileName;
    throw std::runtime_error(oss.str());
  }
  log.Log() << "Loading surface mesh with transform " << transform.PrintStr();

  // Reading every line and determining size
  std::string file_line;
  std::string delimiter = " ";
  while (std::getline(file, file_line))
  {
    // Get the first word
    size_t beg_of_word = file_line.find_first_not_of(delimiter);
    size_t end_of_word = file_line.find(delimiter, beg_of_word - beg_of_word);
    std::string first_word = file_line.substr(beg_of_word, end_of_word);
    std::string sub_word;

    // Keyword "v" for Vertex
    if (first_word == "v")
    {
      Vector3 newVertex;
      for (int k = 1; k <= 3; ++k)
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
            newVertex.x = numValue + transform.x;
          else if (k == 2)
            newVertex.y = numValue + transform.y;
          else if (k == 3)
            newVertex.z = numValue + transform.z;
        }

        // Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          std::cout << "Exception caught!" << std::endl;
        }

        // Stop word extraction on line end
        if (end_of_word == std::string::npos)
        {
          break;
        }
      }
      this->vertices_.push_back(newVertex);
    }

    // Keyword "vt" for Vertex
    if (first_word.compare("vt") == 0)
    {
      Vector3 newVertex;
      for (int k = 1; k <= 2; ++k)
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
          {
            newVertex.x = numValue;
          }
          else if (k == 2)
          {
            newVertex.y = numValue;
          }
          else if (k == 3)
          {
            newVertex.z = numValue;
          }
        }

        // Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          //*newVertex->value[k-1] = 0.0;
          std::cout << "Exception caught!" << std::endl;
        }

        // Stop word extraction on line end
        if (end_of_word == std::string::npos)
        {
          break;
        }
      }
      this->tex_vertices_.push_back(newVertex);
    }

    // Keyword "vn" for normal
    if (first_word == "vn")
    {
      Vector3 newNormal;
      for (int k = 1; k <= 3; ++k)
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
          {
            newNormal.x = numValue;
          }
          else if (k == 2)
          {
            newNormal.y = numValue;
          }
          else if (k == 3)
          {
            newNormal.z = numValue;
          }
        }

        // Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          //*newNormal->value[k-1] = 0.0;
          std::cout << "Exception caught!" << std::endl;
        }

        // Stop word extraction on line end
        if (end_of_word == std::string::npos)
        {
          break;
        }
      }
      this->normals_.push_back(newNormal);
    }

    // Keyword "f" for face
    if (first_word.compare("f") == 0)
    {
      int number_of_verts = std::count(file_line.begin(), file_line.end(), '/') / 2;
      if ((number_of_verts == 3) and (not as_poly))
      {
        auto newFace = std::make_shared<Face>();

        for (int k = 1; k <= 3; ++k)
        {
          // Extract sub word
          beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
          end_of_word = file_line.find(delimiter, beg_of_word);
          sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

          // Extract locations of hiphens
          size_t first_dash = sub_word.find('/');
          size_t last_dash = sub_word.find_last_of('/');

          // Extract the words ass. w vertex and normal
          std::string vert_word = sub_word.substr(0, first_dash - 0);
          std::string norm_word = sub_word.substr(last_dash + 1, sub_word.length() - last_dash - 1);

          // Convert word to number (Vertex)
          try
          {
            int numValue = std::stoi(vert_word);
            newFace->v_index[k - 1] = numValue - 1;
          }
          catch (const std::invalid_argument& ia)
          {
            std::cout << "Exception caught!" << std::endl;
          }

          // Convert word to number (Normal)
          try
          {
            int numValue = std::stoi(norm_word);
            newFace->n_index[k - 1] = numValue - 1;
          }
          catch (const std::invalid_argument& ia)
          {
            std::cout << "Exception caught!" << std::endl;
          }

          // Convert word to number (Texture Vertex)
          if (last_dash > (first_dash + 1))
          {
            std::string tvert_word = sub_word.substr(first_dash + 1, last_dash - first_dash - 1);
            try
            {
              int numValue = std::stoi(tvert_word);
              newFace->vt_index[k - 1] = numValue - 1;
            }
            catch (const std::invalid_argument& ia)
            {
              std::cout << "Exception caught!" << std::endl;
            }
          }

          // Stop word extraction on line end
          if (end_of_word == std::string::npos)
          {
            break;
          }
        }

        // Set edges
        newFace->e_index[0][0] = newFace->v_index[0];
        newFace->e_index[0][1] = newFace->v_index[1];

        newFace->e_index[1][0] = newFace->v_index[1];
        newFace->e_index[1][1] = newFace->v_index[2];

        newFace->e_index[2][0] = newFace->v_index[2];
        newFace->e_index[2][1] = newFace->v_index[0];

        this->faces_.push_back(*newFace);
      }
      else
      {
        auto newFace = std::make_shared<PolyFace>();

        for (int k = 1; k <= number_of_verts; ++k)
        {
          // Extract sub word
          beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
          end_of_word = file_line.find(delimiter, beg_of_word);
          sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

          // Extract locations of hiphens
          size_t first_dash = sub_word.find('/');
          size_t last_dash = sub_word.find_last_of('/');

          // Extract the words ass. w vertex and normal
          std::string vert_word = sub_word.substr(0, first_dash - 0);
          std::string norm_word = sub_word.substr(last_dash + 1, sub_word.length() - last_dash - 1);

          // Convert word to number (Vertex)
          try
          {
            int numValue = std::stoi(vert_word);
            newFace->v_indices.push_back(numValue - 1);
          }
          catch (const std::invalid_argument& ia)
          {
            std::cout << "Exception caught!" << std::endl;
          }

          // Convert word to number (Normal)
          try
          {
            int numValue = std::stoi(norm_word);
          }
          catch (const std::invalid_argument& ia)
          {
            std::cout << "Exception caught!" << std::endl;
          }

          // Convert word to number (Texture Vertex)
          if (last_dash > (first_dash + 1))
          {
            std::string tvert_word = sub_word.substr(first_dash + 1, last_dash - first_dash - 1);
            try
            {
              int numValue = std::stoi(tvert_word);
            }
            catch (const std::invalid_argument& ia)
            {
              std::cout << "Exception caught!" << std::endl;
            }
          }

          // Stop word extraction on line end
          if (end_of_word == std::string::npos)
          {
            break;
          }
        }

        for (int v = 0; v < (newFace->v_indices.size()); ++v)
        {
          std::vector<int> side_indices(4);

          side_indices[0] = newFace->v_indices[v];
          if ((v + 1) >= newFace->v_indices.size())
            side_indices[1] = newFace->v_indices[0];
          else
            side_indices[1] = newFace->v_indices[v + 1];
          side_indices[2] = -1;
          side_indices[3] = -1;

          newFace->edges.push_back(side_indices);
        }

        this->poly_faces_.push_back(newFace);
      }
    }

    // Keyword "l" for line
    if (first_word == "l")
    {
      Edge newEdge;

      // Extract sub word
      beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
      end_of_word = file_line.find(delimiter, beg_of_word);
      sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

      // Convert to number
      try
      {
        int numValue = std::stoi(sub_word);
        newEdge.v_index[0] = numValue - 1;
      }
      catch (const std::invalid_argument& ia)
      {
        std::cout << "Exception caught!" << std::endl;
      }

      // Extract sub word
      beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
      end_of_word = file_line.find(delimiter, beg_of_word);
      sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

      // Convert to number
      try
      {
        int numValue = std::stoi(sub_word);
        newEdge.v_index[1] = numValue - 1;
      }
      catch (const std::invalid_argument& ia)
      {
        std::cout << "Exception caught!" << std::endl;
      }

      newEdge.vertices[0] = this->vertices_.at(newEdge.v_index[0]);
      newEdge.vertices[1] = this->vertices_.at(newEdge.v_index[1]);

      this->lines_.push_back(newEdge);
    }
  }
  file.close();

  // Calculate face properties
  std::vector<Face>::iterator curFace;
  for (curFace = this->faces_.begin(); curFace != this->faces_.end(); ++curFace)
  {
    // Calculate geometrical normal
    Vector3 vA = this->vertices_.at(curFace->v_index[0]);
    Vector3 vB = this->vertices_.at(curFace->v_index[1]);
    Vector3 vC = this->vertices_.at(curFace->v_index[2]);

    Vector3 vAB = vB - vA;
    Vector3 vBC = vC - vB;

    curFace->geometric_normal = vAB.Cross(vBC);
    curFace->geometric_normal = curFace->geometric_normal / curFace->geometric_normal.Norm();

    // Calculate Assigned normal
    Vector3 nA = this->normals_.at(curFace->n_index[0]);
    Vector3 nB = this->normals_.at(curFace->n_index[1]);
    Vector3 nC = this->normals_.at(curFace->n_index[2]);

    Vector3 nAvg = (nA + nB + nC) / 3.0;
    nAvg = nAvg / nAvg.Norm();

    curFace->assigned_normal = nAvg;

    // Compute face center
    curFace->face_centroid = (vA + vB + vC) / 3.0;
  }
  for (auto curPFace = this->poly_faces_.begin(); curPFace != this->poly_faces_.end(); ++curPFace)
  {
    Vector3 centroid;
    int num_verts = (*curPFace)->v_indices.size();
    for (int v = 0; v < num_verts; ++v)
      centroid = centroid + vertices_[(*curPFace)->v_indices[v]];

    centroid = centroid / num_verts;

    (*curPFace)->face_centroid = centroid;

    Vector3 n = (vertices_[(*curPFace)->v_indices[1]] - vertices_[(*curPFace)->v_indices[0]])
                  .Cross(centroid - vertices_[(*curPFace)->v_indices[1]]);
    n = n / n.Norm();

    (*curPFace)->geometric_normal = n;
  }

  UpdateInternalConnectivity();

  // Check each vertex is accounted
  log.Log() << "Surface mesh loaded with " << this->faces_.size() << " triangle faces and "
            << this->poly_faces_.size() << " polygon faces.";

  return 0;
}

int
SurfaceMesh::ImportFromTriangleFiles(const char* fileName, bool as_poly = false)
{
  const std::string fname = "SurfaceMesh::ImportFromTriangleFiles";

  std::string node_filename = std::string(fileName) + std::string(".1.node");
  std::string tria_filename = std::string(fileName) + std::string(".1.ele");

  // Opening the node file
  std::ifstream file;
  file.open(node_filename);
  if (not file.is_open())
  {
    std::ostringstream oss;
    oss << fname << ": Failed to open file: " << node_filename;
    throw std::runtime_error(oss.str());
  }

  int num_verts;
  char line[250];
  file >> num_verts;
  file.getline(line, 250);
  for (int v = 1; v <= num_verts; ++v)
  {
    int vert_index;
    Vector3 vertex;
    file >> vert_index >> vertex.x >> vertex.y;
    file.getline(line, 250);

    vertices_.push_back(vertex);
  }

  file.close();

  // Opening the ele file
  file.open(tria_filename);
  if (not file.is_open())
  {
    std::ostringstream oss;
    oss << fname << ": Failed to open file: " << tria_filename;
    throw std::runtime_error(oss.str());
  }

  int num_tris;

  file >> num_tris;
  file.getline(line, 250);
  for (int v = 1; v <= num_tris; ++v)
  {
    int tri_index;
    auto newFace = std::make_shared<PolyFace>();

    int v0, v1, v2;
    file >> tri_index >> v0 >> v1 >> v2;
    file.getline(line, 250);

    newFace->v_indices.resize(3);
    newFace->v_indices[0] = v0 - 1;
    newFace->v_indices[1] = v1 - 1;
    newFace->v_indices[2] = v2 - 1;

    for (int e = 0; e < 3; ++e)
    {
      std::vector<int> side_indices(4);

      if (e < 2)
      {
        side_indices[0] = newFace->v_indices[e];
        side_indices[1] = newFace->v_indices[e + 1];
        side_indices[2] = -1;
        side_indices[3] = -1;
      }
      else
      {
        side_indices[0] = newFace->v_indices[e];
        side_indices[1] = newFace->v_indices[0];
        side_indices[2] = -1;
        side_indices[3] = -1;
      }
      newFace->edges.push_back(side_indices);
    }

    poly_faces_.push_back(newFace);
  }

  file.close();

  // Calculate face properties
  std::vector<Face>::iterator curFace;
  for (curFace = this->faces_.begin(); curFace != this->faces_.end(); ++curFace)
  {
    // Calculate geometrical normal
    Vector3 vA = this->vertices_.at(curFace->v_index[0]);
    Vector3 vB = this->vertices_.at(curFace->v_index[1]);
    Vector3 vC = this->vertices_.at(curFace->v_index[2]);

    Vector3 vAB = vB - vA;
    Vector3 vBC = vC - vB;

    curFace->geometric_normal = vAB.Cross(vBC);
    curFace->geometric_normal = curFace->geometric_normal / curFace->geometric_normal.Norm();

    // Calculate Assigned normal
    Vector3 nA = this->normals_.at(curFace->n_index[0]);
    Vector3 nB = this->normals_.at(curFace->n_index[1]);
    Vector3 nC = this->normals_.at(curFace->n_index[2]);

    Vector3 nAvg = (nA + nB + nC) / 3.0;
    nAvg = nAvg / nAvg.Norm();

    curFace->assigned_normal = nAvg;

    // Compute face center
    curFace->face_centroid = (vA + vB + vC) / 3.0;
  }
  for (auto curPFace = this->poly_faces_.begin(); curPFace != this->poly_faces_.end(); ++curPFace)
  {
    Vector3 centroid;
    int num_verts = (*curPFace)->v_indices.size();
    for (int v = 0; v < num_verts; ++v)
      centroid = centroid + vertices_[(*curPFace)->v_indices[v]];

    centroid = centroid / num_verts;

    (*curPFace)->face_centroid = centroid;

    Vector3 n = (vertices_[(*curPFace)->v_indices[1]] - vertices_[(*curPFace)->v_indices[0]])
                  .Cross(centroid - vertices_[(*curPFace)->v_indices[1]]);
    n = n / n.Norm();

    (*curPFace)->geometric_normal = n;
  }

  UpdateInternalConnectivity();

  // Check each vertex is accounted
  log.Log() << "Surface mesh loaded with " << this->faces_.size() << " triangle faces and "
            << this->poly_faces_.size() << " polygon faces.";

  return 0;
}

int
SurfaceMesh::ImportFromMshFiles(const char* fileName, bool as_poly = false)
{
  const std::string fname = "SurfaceMesh::ImportFromMshFiles";

  const std::string node_section_name = "$Nodes";
  const std::string elements_section_name = "$Elements";

  std::istringstream iss;
  std::string line;

  std::ifstream file;
  file.open(std::string(fileName));

  if (not file.is_open())
  {
    std::ostringstream oss;
    oss << fname << ": Failed to open file: " << fileName;
    throw std::runtime_error(oss.str());
  }

  // Find section with node information and then read information
  while (std::getline(file, line))
  {
    if (node_section_name == line)
      break;
  }

  std::getline(file, line);
  iss = std::istringstream(line);
  int num_nodes;
  if (not(iss >> num_nodes))
  {
    log.LogAllError() << "Failed while trying to read the number of nodes.\n";
    Exit(EXIT_FAILURE);
  }

  vertices_.resize(num_nodes);

  for (int n = 0; n < num_nodes; ++n)
  {
    std::getline(file, line);
    iss = std::istringstream(line);

    Vector3 vertex;
    int vert_index;
    if (not(iss >> vert_index))
    {
      log.LogAllError() << "Failed to read vertex index.\n";
      Exit(EXIT_FAILURE);
    }

    if (not(iss >> vertex.x >> vertex.y >> vertex.z))
    {
      log.LogAllError() << "Failed while reading the vertex coordinates.\n";
      Exit(EXIT_FAILURE);
    }

    vertices_[vert_index - 1] = vertex;
  }

  // Find the element listing section and first read the boundary data
  file.seekg(0);
  while (std::getline(file, line))
  {
    if (elements_section_name.compare(line) == 0)
      break;
  }

  std::getline(file, line);
  iss = std::istringstream(line);
  int num_elems;
  if (not(iss >> num_elems))
  {
    log.LogAllError() << "Failed to read number of elements.\n";
    Exit(EXIT_FAILURE);
  }

  for (int n = 0; n < num_elems; ++n)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;
    auto newFace = std::make_shared<PolyFace>();
    std::getline(file, line);
    iss = std::istringstream(line);

    if (not(iss >> element_index >> elem_type >> num_tags))
    {
      log.LogAllError() << "Failed while reading element index, element "
                           "type, and number of tags.\n";
      Exit(EXIT_FAILURE);
    }

    if (not(iss >> physical_reg))
    {
      log.LogAllError() << "Failed while reading physical region.\n";
      Exit(EXIT_FAILURE);
    }

    for (int i = 1; i < num_tags; ++i)
      if (not(iss >> tag))
      {
        log.LogAllError() << "Failed when reading tags.\n";
        Exit(EXIT_FAILURE);
      }

    if (elem_type == 2)
    {
      const int num_nodes = 3;

      int nodes[num_nodes];
      for (int i = 0; i < num_nodes; ++i)
        if (not(iss >> nodes[i]))
        {
          log.LogAllError() << "Failed when reading element node index.\n";
          Exit(EXIT_FAILURE);
        }

      newFace->v_indices.resize(num_nodes);
      for (int i = 0; i < num_nodes; ++i)
        newFace->v_indices[i] = nodes[i] - 1;
    }
    else if (elem_type == 3)
    {
      const int num_nodes = 4;

      int nodes[num_nodes];
      for (int& node : nodes)
        if (not(iss >> node))
        {
          log.LogAllError() << "Failed when reading element node index.\n";
          Exit(EXIT_FAILURE);
        }

      newFace->v_indices.resize(num_nodes);
      for (int i = 0; i < num_nodes; ++i)
        newFace->v_indices[i] = nodes[i] - 1;
    }
    else
    {
      continue;
    }

    const size_t total_nodes = newFace->v_indices.size();

    for (size_t e = 0; e < total_nodes; ++e)
    {
      std::vector<int> side_indices(total_nodes);
      side_indices[0] = newFace->v_indices[e];

      if (e < total_nodes - 1)
        side_indices[1] = newFace->v_indices[e + 1];
      else
        side_indices[1] = newFace->v_indices[0];

      side_indices[2] = -1;
      side_indices[3] = -1;

      newFace->edges.push_back(side_indices);
    }

    poly_faces_.push_back(newFace);
    physical_region_map_.push_back(physical_reg);
  }

  file.close();

  // Calculate face properties
  for (const auto& poly_face : poly_faces_)
  {
    Vector3 centroid;
    size_t num_verts = poly_face->v_indices.size();

    for (size_t v = 0; v < num_verts; ++v)
      centroid = centroid + vertices_[poly_face->v_indices[v]];

    centroid = centroid / static_cast<double>(num_verts);

    poly_face->face_centroid = centroid;

    Vector3 n = (vertices_[poly_face->v_indices[1]] - vertices_[poly_face->v_indices[0]])
                  .Cross(centroid - vertices_[poly_face->v_indices[1]]);

    n = n / n.Norm();

    poly_face->geometric_normal = n;
  }

  UpdateInternalConnectivity();

  return 0;
}

void
SurfaceMesh::ComputeLoadBalancing(std::vector<double>& x_cuts, std::vector<double>& y_cuts)
{
  log.Log() << "X-cuts to be logged: " << x_cuts.size();
  //  for (auto& val : x_cuts)
  //    log.Log() << val;
  //
  log.Log() << "Y-cuts to be logged: " << y_cuts.size();
  //  for (auto& val : y_cuts)
  //    log.Log() << val;

  // Sort faces into bins
  size_t I = x_cuts.size();
  size_t J = y_cuts.size();

  std::vector<std::vector<int>> IJ_bins(I + 1, std::vector<int>(J + 1, 0));

  for (auto& poly_face : poly_faces_)
  {
    int ref_i = 0;
    int ref_j = 0;
    for (size_t i = 0; i < I; ++i)
    {
      if (poly_face->face_centroid.x >= x_cuts[i])
        ref_i = i + 1;
    } // for i
    for (size_t j = 0; j < J; ++j)
    {
      if (poly_face->face_centroid.y >= y_cuts[j])
        ref_j = j + 1;
    } // for j

    IJ_bins[ref_i][ref_j] += 1;
  } // for face

  // Determine average and max
  int max_bin_size = 0;
  int tot_bin_size = 0;
  int i_max = 0, j_max = 0;

  for (int i = 0; i < (I + 1); ++i)
  {
    for (int j = 0; j < (J + 1); ++j)
    {
      if (IJ_bins[i][j] > max_bin_size)
      {
        max_bin_size = IJ_bins[i][j];
        i_max = i;
        j_max = j;
      }
      tot_bin_size += IJ_bins[i][j];
    }
  }

  double average = tot_bin_size / ((double)(I + 1) * (J + 1));

  log.Log() << "Average faces per set: " << average;
  log.Log() << "Maximum faces per set: " << max_bin_size << " at (i,j)= ( " << i_max << " , "
            << j_max << " )";

  if (i_max == I)
    log.Log() << "X greater than " << x_cuts[i_max - 1];
  else if (i_max == 0)
    log.Log() << "X less than " << x_cuts[0];
  else
    log.Log() << "X greater than " << x_cuts[i_max - 1] << " and less than " << x_cuts[i_max];

  if (j_max == J)
    log.Log() << "Y greater than " << y_cuts[j_max - 1];
  else if (j_max == 0)
    log.Log() << "Y less than " << y_cuts[0];
  else
    log.Log() << "Y greater than " << y_cuts[j_max - 1] << " and less than " << y_cuts[j_max];

  log.Log() << "Max-to-average ratio: " << max_bin_size / average;
}

} // namespace opensn
