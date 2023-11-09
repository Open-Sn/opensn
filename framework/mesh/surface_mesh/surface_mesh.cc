#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/graphs/directed_graph.h"
#include <algorithm>
#include <iostream>
#include <fstream>

namespace chi_mesh
{

SurfaceMesh::SurfaceMesh()
{
}

SurfaceMesh::~SurfaceMesh()
{
  for (auto poly_face : poly_faces_)
  {
    delete poly_face;
  }

  poly_faces_.clear();
}

std::ostream&
operator<<(std::ostream& os, SurfaceMesh& that)
{
  std::vector<Face>::const_iterator curface;
  for (curface = that.GetTriangles().begin(); curface != that.GetTriangles().end(); curface++)
  {
    long index = std::distance(that.GetTriangles().begin(), curface);
    os << "Face " << index << " v:";
    os << curface->v_index[0] << "->";
    os << curface->v_index[1] << "->";
    os << curface->v_index[2] << "  ";

    os << "e:";
    for (int e = 0; e < 3; e++)
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

bool
SurfaceMesh::CheckNegativeSense(double x, double y, double z)
{
  Vector3 xyz = Vector3(x, y, z);

  // Loop through each face
  std::vector<Face>::iterator cur_face;
  for (cur_face = this->faces_.begin(); cur_face != this->faces_.end(); cur_face++)
  {
    // Get a vertex (first one)
    Vertex p;
    try
    {
      p = this->vertices_.at(cur_face->v_index[0]);
    }
    catch (const std::out_of_range& o)
    {
      std::cout << "Invalid vertex handle" << std::endl;
    }

    // Calculate dot product
    Vector3 p_xyz = xyz - p;
    double dprod = cur_face->assigned_normal.Dot(p_xyz);

    if (dprod < 0.0) { return true; }
  }

  return false;
}

void
SurfaceMesh::ExtractOpenEdgesToObj(const char* fileName)
{
  std::vector<std::pair<int, int>> edges;
  for (auto face : poly_faces_)
  {
    for (auto edge : face->edges)
    {
      if (edge[2] < 0) { edges.push_back(std::pair<int, int>(edge[0], edge[1])); }
    } // for edge
  }   // for face

  std::ofstream outfile;
  outfile.open(fileName);

  if (!outfile.is_open())
  {
    Chi::log.LogAllError() << "In call to SurfaceMesh::ExtractOpenEdgesToObj. Failed"
                           << " to open file: " << std::string(fileName);
    Chi::Exit(EXIT_FAILURE);
  }

  outfile << "# ChiTech open edges file\n";
  outfile << "# Single surface mesh\n";

  for (auto vert_pair : edges)
  {
    Vertex& v0 = vertices_[vert_pair.first];
    Vertex& v1 = vertices_[vert_pair.second];
    outfile << "v " << v0.x << " " << v0.y << " " << v0.z << "\n"
            << "v " << v1.x << " " << v1.y << " " << v1.z << "\n";
  }

  for (size_t e = 0; e < edges.size(); ++e)
  {
    const auto v_count = 2 * e + 1;
    outfile << "l " << v_count << " " << v_count + 1 << "\n";
  }

  outfile.close();
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
  //%%%%%% TRIANGLES %%%%%
  size_t num_tri_faces = faces_.size();
  for (size_t tf = 0; tf < num_tri_faces; tf++)
  {
    auto& try_face = faces_[tf];
    for (int v = 0; v < 3; ++v)
      vertex_subscriptions[v].push_back(tf);
  }
  //%%%%% POLYGONS %%%%%
  size_t num_poly_faces = poly_faces_.size();
  for (size_t pf = 0; pf < num_poly_faces; pf++)
  {
    auto poly_face = poly_faces_[pf];
    for (auto v : poly_face->v_indices)
      vertex_subscriptions[v].push_back(pf);
  }

  // Loop over cells and determine connectivity
  //%%%%%% TRIANGLES %%%%%
  for (auto curFace : faces_)
  {
    for (int e = 0; e < 3; e++)
    {
      int* curface_edge = curFace.e_index[e];
      int vi = curface_edge[0];
      int vf = curface_edge[1];

      // Search cells subscribing to vi
      for (auto ofi : vertex_subscriptions[vi])
      {
        auto& other_face = faces_[ofi];

        for (size_t e2 = 0; e2 < 3; e2++)
        {
          if ((curface_edge[0] == other_face.e_index[e2][1]) &&
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

        for (size_t e2 = 0; e2 < 3; e2++)
        {
          if ((curface_edge[0] == other_face.e_index[e2][1]) &&
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
  //%%%%% POLYGONS %%%%%
  for (auto curFace : poly_faces_)
  {
    for (auto& curface_edge : curFace->edges)
    {
      int vi = curface_edge[0];
      int vf = curface_edge[1];

      // Search cells subscribing to vi
      for (auto ofi : vertex_subscriptions[vi])
      {
        auto other_face = poly_faces_[ofi];

        for (size_t e2 = 0; e2 < other_face->edges.size(); e2++)
        {
          if ((curface_edge[0] == other_face->edges[e2][1]) &&
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

        for (size_t e2 = 0; e2 < other_face->edges.size(); e2++)
        {
          if ((curface_edge[0] == other_face->edges[e2][1]) &&
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

  // Opening the file
  std::ifstream file;
  file.open(fileName);
  if (!file.is_open())
  {
    Chi::log.LogAllError() << "Failed to open file: " << fileName << " in call "
                           << "to ImportFromOBJFile \n";
    Chi::Exit(EXIT_FAILURE);
  }
  Chi::log.Log() << "Loading surface mesh with transform " << transform.PrintStr();

  // Reading every line and determining size
  std::string file_line;
  std::string delimiter = " ";
  int counter = 0;
  while (std::getline(file, file_line))
  {
    counter++;

    // Get the first word
    size_t beg_of_word = file_line.find_first_not_of(delimiter);
    size_t end_of_word = file_line.find(delimiter, beg_of_word - beg_of_word);
    std::string first_word = file_line.substr(beg_of_word, end_of_word);
    std::string sub_word;

    // Keyword "v" for Vertex
    if (first_word == "v")
    {
      Vertex newVertex;
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

          if (k == 1) newVertex.x = numValue + transform.x;
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
        if (end_of_word == std::string::npos) { break; }
      }
      this->vertices_.push_back(newVertex);
    }

    // Keyword "vt" for Vertex
    if (first_word.compare("vt") == 0)
    {
      Vertex newVertex;
      for (int k = 1; k <= 2; k++)
      {
        // Extract sub word
        beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
        end_of_word = file_line.find(delimiter, beg_of_word);
        sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

        // Convert word to number
        try
        {
          double numValue = std::stod(sub_word);
          //*newVertex->value[k-1] = numValue;
          if (k == 1) { newVertex.x = numValue; }
          else if (k == 2) { newVertex.y = numValue; }
          else if (k == 3) { newVertex.z = numValue; }
        }

        // Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          //*newVertex->value[k-1] = 0.0;
          std::cout << "Exception caught!" << std::endl;
        }

        // Stop word extraction on line end
        if (end_of_word == std::string::npos) { break; }
      }
      this->tex_vertices_.push_back(newVertex);
    }

    // Keyword "vn" for normal
    if (first_word.compare("vn") == 0)
    {
      Normal newNormal;
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
          //*newNormal->value[k-1] = numValue;
          if (k == 1) { newNormal.x = numValue; }
          else if (k == 2) { newNormal.y = numValue; }
          else if (k == 3) { newNormal.z = numValue; }
        }

        // Catch conversion error
        catch (const std::invalid_argument& ia)
        {
          //*newNormal->value[k-1] = 0.0;
          std::cout << "Exception caught!" << std::endl;
        }

        // Stop word extraction on line end
        if (end_of_word == std::string::npos) { break; }
      }
      this->normals_.push_back(newNormal);
    }

    // Keyword "f" for face
    if (first_word.compare("f") == 0)
    {
      int number_of_verts = std::count(file_line.begin(), file_line.end(), '/') / 2;
      if ((number_of_verts == 3) && (!as_poly))
      {
        Face* newFace = new Face;

        for (int k = 1; k <= 3; k++)
        {
          // Extract sub word
          beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
          end_of_word = file_line.find(delimiter, beg_of_word);
          sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

          // Extract locations of hiphens
          size_t first_dash = sub_word.find("/");
          size_t last_dash = sub_word.find_last_of("/");

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
          if (end_of_word == std::string::npos) { break; }
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
        PolyFace* newFace = new PolyFace;

        for (int k = 1; k <= number_of_verts; k++)
        {
          // Extract sub word
          beg_of_word = file_line.find_first_not_of(delimiter, end_of_word);
          end_of_word = file_line.find(delimiter, beg_of_word);
          sub_word = file_line.substr(beg_of_word, end_of_word - beg_of_word);

          // Extract locations of hiphens
          size_t first_dash = sub_word.find("/");
          size_t last_dash = sub_word.find_last_of("/");

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
            // newFace->n_indices.push_back(numValue-1);
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
              // newFace->vt_indices.push_back(numValue-1);
            }
            catch (const std::invalid_argument& ia)
            {
              std::cout << "Exception caught!" << std::endl;
            }
          }

          // Stop word extraction on line end
          if (end_of_word == std::string::npos) { break; }
        }

        for (int v = 0; v < (newFace->v_indices.size()); v++)
        {
          int* side_indices = new int[4];

          side_indices[0] = newFace->v_indices[v];
          side_indices[1] = newFace->v_indices[v + 1];
          side_indices[2] = -1;
          side_indices[3] = -1;

          if ((v + 1) >= newFace->v_indices.size()) { side_indices[1] = newFace->v_indices[0]; }

          newFace->edges.push_back(side_indices);
        }

        this->poly_faces_.push_back(newFace);
      }
    }

    // Keyword "l" for line
    if (first_word.compare("l") == 0)
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

      // printf("line %d->%d\n",newEdge.v_index[0],newEdge.v_index[1]);
    }
  }
  file.close();

  // Calculate face properties
  std::vector<Face>::iterator curFace;
  for (curFace = this->faces_.begin(); curFace != this->faces_.end(); curFace++)
  {
    // Calculate geometrical normal
    Vertex vA = this->vertices_.at(curFace->v_index[0]);
    Vertex vB = this->vertices_.at(curFace->v_index[1]);
    Vertex vC = this->vertices_.at(curFace->v_index[2]);

    Vector3 vAB = vB - vA;
    Vector3 vBC = vC - vB;

    curFace->geometric_normal = vAB.Cross(vBC);
    curFace->geometric_normal = curFace->geometric_normal / curFace->geometric_normal.Norm();

    // Calculate Assigned normal
    Vertex nA = this->normals_.at(curFace->n_index[0]);
    Vertex nB = this->normals_.at(curFace->n_index[1]);
    Vertex nC = this->normals_.at(curFace->n_index[2]);

    Vector3 nAvg = (nA + nB + nC) / 3.0;
    nAvg = nAvg / nAvg.Norm();

    curFace->assigned_normal = nAvg;

    // Compute face center
    curFace->face_centroid = (vA + vB + vC) / 3.0;
  }
  std::vector<PolyFace*>::iterator curPFace;
  for (curPFace = this->poly_faces_.begin(); curPFace != this->poly_faces_.end(); curPFace++)
  {
    Vector3 centroid;
    int num_verts = (*curPFace)->v_indices.size();
    for (int v = 0; v < num_verts; v++)
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
  Chi::log.Log() << "Surface mesh loaded with " << this->faces_.size() << " triangle faces and "
                 << this->poly_faces_.size() << " polygon faces.";
  // chi::Exit(EXIT_FAILURE);

  return 0;
}

int
SurfaceMesh::ImportFromTriangleFiles(const char* fileName, bool as_poly = false)
{
  std::string node_filename = std::string(fileName) + std::string(".1.node");
  std::string tria_filename = std::string(fileName) + std::string(".1.ele");

  // Opening the node file
  std::ifstream file;
  file.open(node_filename);
  if (!file.is_open())
  {
    Chi::log.LogAllError() << "Failed to open file: " << node_filename << " in call "
                           << "to ImportFromOBJFile \n";
    Chi::Exit(EXIT_FAILURE);
  }

  int num_verts;
  char line[250];
  file >> num_verts;
  file.getline(line, 250);
  for (int v = 1; v <= num_verts; v++)
  {
    int vert_index;
    Vertex vertex;
    file >> vert_index >> vertex.x >> vertex.y;
    file.getline(line, 250);

    vertices_.push_back(vertex);
  }

  file.close();

  // Opening the ele file
  file.open(tria_filename);
  if (!file.is_open())
  {
    Chi::log.LogAllError() << "Failed to open file: " << tria_filename << " in call "
                           << "to ImportFromOBJFile \n";
    Chi::Exit(EXIT_FAILURE);
  }

  int num_tris;

  file >> num_tris;
  file.getline(line, 250);
  for (int v = 1; v <= num_tris; v++)
  {
    int tri_index;
    PolyFace* newFace = new PolyFace;

    int v0, v1, v2;
    file >> tri_index >> v0 >> v1 >> v2;
    file.getline(line, 250);

    newFace->v_indices.resize(3);
    newFace->v_indices[0] = v0 - 1;
    newFace->v_indices[1] = v1 - 1;
    newFace->v_indices[2] = v2 - 1;

    for (int e = 0; e < 3; e++)
    {
      int* side_indices = new int[4];

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
  for (curFace = this->faces_.begin(); curFace != this->faces_.end(); curFace++)
  {
    // Calculate geometrical normal
    Vertex vA = this->vertices_.at(curFace->v_index[0]);
    Vertex vB = this->vertices_.at(curFace->v_index[1]);
    Vertex vC = this->vertices_.at(curFace->v_index[2]);

    Vector3 vAB = vB - vA;
    Vector3 vBC = vC - vB;

    curFace->geometric_normal = vAB.Cross(vBC);
    curFace->geometric_normal = curFace->geometric_normal / curFace->geometric_normal.Norm();

    // Calculate Assigned normal
    Vertex nA = this->normals_.at(curFace->n_index[0]);
    Vertex nB = this->normals_.at(curFace->n_index[1]);
    Vertex nC = this->normals_.at(curFace->n_index[2]);

    Vector3 nAvg = (nA + nB + nC) / 3.0;
    nAvg = nAvg / nAvg.Norm();

    curFace->assigned_normal = nAvg;

    // Compute face center
    curFace->face_centroid = (vA + vB + vC) / 3.0;
  }
  std::vector<PolyFace*>::iterator curPFace;
  for (curPFace = this->poly_faces_.begin(); curPFace != this->poly_faces_.end(); curPFace++)
  {
    Vector3 centroid;
    int num_verts = (*curPFace)->v_indices.size();
    for (int v = 0; v < num_verts; v++)
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
  Chi::log.Log() << "Surface mesh loaded with " << this->faces_.size() << " triangle faces and "
                 << this->poly_faces_.size() << " polygon faces.";
  // chi::Exit(EXIT_FAILURE);

  return 0;
}

SurfaceMesh*
SurfaceMesh::CreateFromDivisions(std::vector<double>& vertices_1d_x,
                                 std::vector<double>& vertices_1d_y)
{
  // Checks if vertices are empty
  if (vertices_1d_x.empty())
  {
    Chi::log.LogAllError() << "SurfaceMesh::CreateFromDivisions. Empty vertex_x list.";
    Chi::Exit(EXIT_FAILURE);
  }
  if (vertices_1d_y.empty())
  {
    Chi::log.LogAllError() << "SurfaceMesh::CreateFromDivisions. Empty vertex_y list.";
    Chi::Exit(EXIT_FAILURE);
  }

  // Populate 2D vertices
  int Nvx = vertices_1d_x.size();
  int Nvy = vertices_1d_y.size();

  int Ncx = Nvx - 1;
  int Ncy = Nvy - 1;

  std::vector<Vertex> vertices_x;
  std::vector<Vertex> vertices_y;

  vertices_x.reserve(Nvx);
  vertices_y.reserve(Nvy);

  for (double v : vertices_1d_x)
    vertices_x.emplace_back(v, 0.0, 0.0);

  for (double v : vertices_1d_y)
    vertices_y.emplace_back(0.0, v, 0.0);

  // Create surface mesh
  auto surf_mesh = new SurfaceMesh();

  // Populate vertices
  std::vector<std::vector<int>> vert_ij_map(Nvx, std::vector<int>(Nvx, -1));
  for (int i = 0; i < Nvy; ++i)
  {
    for (int j = 0; j < Nvx; ++j)
    {
      surf_mesh->vertices_.push_back(vertices_x[j] + vertices_y[i]);
      vert_ij_map[i][j] = surf_mesh->vertices_.size() - 1;
    } // for j
  }   // for i

  // Populate polyfaces
  for (int i = 0; i < Ncy; ++i)
  {
    for (int j = 0; j < Ncx; ++j)
    {
      auto new_face = new PolyFace();
      new_face->v_indices.push_back(vert_ij_map[i][j]);
      new_face->v_indices.push_back(vert_ij_map[i][j + 1]);
      new_face->v_indices.push_back(vert_ij_map[i + 1][j + 1]);
      new_face->v_indices.push_back(vert_ij_map[i + 1][j]);

      for (int v = 0; v < (new_face->v_indices.size()); v++)
      {
        int* side_indices = new int[4];

        side_indices[0] = new_face->v_indices[v];
        side_indices[1] = new_face->v_indices[v + 1];
        side_indices[2] = -1;
        side_indices[3] = -1;

        if ((v + 1) >= new_face->v_indices.size()) side_indices[1] = new_face->v_indices[0];

        new_face->edges.push_back(side_indices);
      } // for v

      surf_mesh->poly_faces_.push_back(new_face);
    } // for j
  }   // for i

  // Compute normals
  for (auto poly_face : surf_mesh->poly_faces_)
  {
    Vector3 centroid;
    int num_verts = poly_face->v_indices.size();
    for (int v = 0; v < num_verts; v++)
      centroid = centroid + surf_mesh->vertices_[poly_face->v_indices[v]];

    centroid = centroid / num_verts;

    poly_face->face_centroid = centroid;

    Vector3 n = (surf_mesh->vertices_[poly_face->v_indices[1]] -
                 surf_mesh->vertices_[poly_face->v_indices[0]])
                  .Cross(centroid - surf_mesh->vertices_[poly_face->v_indices[1]]);
    n = n / n.Norm();

    poly_face->geometric_normal = n;
  }

  surf_mesh->UpdateInternalConnectivity();

  return surf_mesh;
}

int
SurfaceMesh::ImportFromMshFiles(const char* fileName, bool as_poly = false)
{
  const std::string node_section_name = "$Nodes";
  const std::string elements_section_name = "$Elements";

  std::istringstream iss;
  std::string line;

  std::ifstream file;
  file.open(std::string(fileName));

  if (!file.is_open())
  {
    Chi::log.LogAllError() << "Failed to open file: " << fileName << " in call "
                           << "to ImportFromMshFiles \n";
    Chi::Exit(EXIT_FAILURE);
  }

  // Find section with node information and then read information
  while (std::getline(file, line))
  {
    if (node_section_name.compare(line) == 0) break;
  }

  std::getline(file, line);
  iss = std::istringstream(line);
  int num_nodes;
  if (!(iss >> num_nodes))
  {
    Chi::log.LogAllError() << "Failed while trying to read the number of nodes.\n";
    Chi::Exit(EXIT_FAILURE);
  }

  vertices_.resize(num_nodes);

  for (int n = 0; n < num_nodes; n++)
  {
    std::getline(file, line);
    iss = std::istringstream(line);

    Vertex vertex;
    int vert_index;
    if (!(iss >> vert_index))
    {
      Chi::log.LogAllError() << "Failed to read vertex index.\n";
      Chi::Exit(EXIT_FAILURE);
    }

    if (!(iss >> vertex.x >> vertex.y >> vertex.z))
    {
      Chi::log.LogAllError() << "Failed while reading the vertex coordinates.\n";
      Chi::Exit(EXIT_FAILURE);
    }

    vertices_[vert_index - 1] = vertex;
  }

  // Find the element listing section and first read the boundary data
  file.seekg(0);
  while (std::getline(file, line))
  {
    if (elements_section_name.compare(line) == 0) break;
  }

  std::getline(file, line);
  iss = std::istringstream(line);
  int num_elems;
  if (!(iss >> num_elems))
  {
    Chi::log.LogAllError() << "Failed to read number of elements.\n";
    Chi::Exit(EXIT_FAILURE);
  }

  for (int n = 0; n < num_elems; n++)
  {
    int elem_type, num_tags, physical_reg, tag, element_index;
    PolyFace* newFace = new PolyFace;
    std::getline(file, line);
    iss = std::istringstream(line);

    if (!(iss >> element_index >> elem_type >> num_tags))
    {
      Chi::log.LogAllError() << "Failed while reading element index, element "
                                "type, and number of tags.\n";
      Chi::Exit(EXIT_FAILURE);
    }

    if (!(iss >> physical_reg))
    {
      Chi::log.LogAllError() << "Failed while reading physical region.\n";
      Chi::Exit(EXIT_FAILURE);
    }

    for (int i = 1; i < num_tags; i++)
      if (!(iss >> tag))
      {
        Chi::log.LogAllError() << "Failed when reading tags.\n";
        Chi::Exit(EXIT_FAILURE);
      }

    if (elem_type == 2)
    {
      const int num_nodes = 3;

      int nodes[num_nodes];
      for (int i = 0; i < num_nodes; i++)
        if (!(iss >> nodes[i]))
        {
          Chi::log.LogAllError() << "Failed when reading element node index.\n";
          Chi::Exit(EXIT_FAILURE);
        }

      newFace->v_indices.resize(num_nodes);
      for (int i = 0; i < num_nodes; i++)
        newFace->v_indices[i] = nodes[i] - 1;
    }
    else if (elem_type == 3)
    {
      const int num_nodes = 4;

      int nodes[num_nodes];
      for (int& node : nodes)
        if (!(iss >> node))
        {
          Chi::log.LogAllError() << "Failed when reading element node index.\n";
          Chi::Exit(EXIT_FAILURE);
        }

      newFace->v_indices.resize(num_nodes);
      for (int i = 0; i < num_nodes; i++)
        newFace->v_indices[i] = nodes[i] - 1;
    }
    else { continue; }

    const size_t total_nodes = newFace->v_indices.size();

    for (size_t e = 0; e < total_nodes; e++)
    {
      int* side_indices = new int[total_nodes];
      side_indices[0] = newFace->v_indices[e];

      if (e < total_nodes - 1) side_indices[1] = newFace->v_indices[e + 1];
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

    for (size_t v = 0; v < num_verts; v++)
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
SurfaceMesh::ExportToOBJFile(const char* fileName)
{

  //  if (this->faces.empty())
  //  {
  //    std::cout << "Cannot export empty SurfaceMesh\n";
  //    return;
  //  }
  FILE* outputFile = fopen(fileName, "w");
  if (outputFile == nullptr)
  {
    printf("Error creating file %s!\n", fileName);
    return;
  }

  fprintf(outputFile, "# Exported mesh file from tringulation script\n");
  fprintf(outputFile, "o %s\n", "ChitechTriMesh");

  std::vector<Vertex>::iterator cur_v;
  for (cur_v = this->vertices_.begin(); cur_v != this->vertices_.end(); cur_v++)
  {
    fprintf(outputFile, "v %9.6f %9.6f %9.6f\n", cur_v->x, cur_v->y, cur_v->z);
  }

  for (unsigned ell = 0; ell < this->lines_.size(); ell++)
  {
    fprintf(outputFile, "l %d %d \n", lines_[ell].v_index[0] + 1, lines_[ell].v_index[1] + 1);
  }

  if (!faces_.empty())
  {
    Face first_face = this->faces_.front();
    fprintf(outputFile,
            "vn %.4f %.4f %.4f\n",
            first_face.geometric_normal.x,
            first_face.geometric_normal.y,
            first_face.geometric_normal.z);
    fprintf(outputFile, "s off\n");

    std::vector<Face>::iterator cur_face;
    for (cur_face = this->faces_.begin(); cur_face != this->faces_.end(); cur_face++)
    {
      fprintf(outputFile,
              "f %d//1 %d//1 %d//1\n",
              cur_face->v_index[0] + 1,
              cur_face->v_index[1] + 1,
              cur_face->v_index[2] + 1);
    }
  }
  if (!poly_faces_.empty())
  {
    PolyFace* first_face = this->poly_faces_.front();
    fprintf(outputFile,
            "vn %.4f %.4f %.4f\n",
            first_face->geometric_normal.x,
            first_face->geometric_normal.y,
            first_face->geometric_normal.z);
    fprintf(outputFile, "s off\n");

    for (auto& poly_face : poly_faces_)
    {
      fprintf(outputFile, "f ");
      for (int v_indice : poly_face->v_indices)
      {
        fprintf(outputFile, "%d//1 ", v_indice + 1);
      }
      fprintf(outputFile, "\n");
    }
  }

  fclose(outputFile);
  printf("Exported mesh to %s\n", fileName);
}

void
SurfaceMesh::ExportToPolyFile(const char* fileName)
{
  FILE* outputFile = fopen(fileName, "w");
  if (outputFile == nullptr)
  {
    printf("Error creating file %s!\n", fileName);
    return;
  }

  fprintf(outputFile, "%lu 2 0 0\n", vertices_.size());
  for (int v = 0; v < vertices_.size(); v++)
  {
    fprintf(outputFile, "%d %.15f %.15f 0\n", v + 1, vertices_[v].x, vertices_[v].y);
  }

  fprintf(outputFile, "%lu 0\n", lines_.size());
  for (int e = 0; e < lines_.size(); e++)
  {
    fprintf(outputFile, "%d %d %d\n", e + 1, lines_[e].v_index[0] + 1, lines_[e].v_index[1] + 1);
  }

  fprintf(outputFile, "0");

  fclose(outputFile);
  printf("Exported mesh to %s\n", fileName);
}

void
SurfaceMesh::CheckCyclicDependencies(int num_angles)
{
  double tolerance = 1.0e-8;
  double dvarphi = 2.0 * M_PI / num_angles;
  Vector3 khat(0.0, 0.0, 1.0);

  // Loop over angles
  for (int a = 0; a < num_angles; a++)
  {
    double varphi = 0.5 * dvarphi + a * dvarphi;

    Vector3 omega;
    omega.x = cos(varphi);
    omega.y = sin(varphi);
    omega.z = 0.0;

    // Add all polyfaces to graph
    chi::DirectedGraph G;
    size_t num_loc_cells = poly_faces_.size();
    for (size_t c = 0; c < num_loc_cells; c++)
      G.AddVertex();

    // Now construct dependencies
    for (size_t c = 0; c < num_loc_cells; c++)
    {
      PolyFace* face = poly_faces_[c];

      size_t num_edges = face->edges.size();
      for (size_t e = 0; e < num_edges; e++)
      {
        int v0i = face->edges[e][0];
        int v1i = face->edges[e][1];

        Vector3 v01 = vertices_[v1i] - vertices_[v0i];
        Vector3 n = v01.Cross(khat);
        n = n / n.Norm();

        double mu = omega.Dot(n);
        int neighbor = face->edges[e][2];
        if ((mu > (0.0 + tolerance)) and (neighbor >= 0))
        {
          //          boost::add_edge(c,neighbor,G);
          G.AddEdge(c, neighbor);
        } // if outgoing
      }   // for edge
    }     // for cell

    // Generic topological
    //                                                   sorting
    //    typedef boost::graph_traits<CHI_D_GRAPH>::vertex_descriptor gVertex;
    //
    //    boost::property_map<CHI_D_GRAPH, boost::vertex_index_t>::type
    //      index_map = get(boost::vertex_index, G);
    //
    //    std::vector<gVertex> sorted_list;
    //    try{
    //      boost::topological_sort(G,std::back_inserter(sorted_list));
    //    }
    //    catch (const boost::bad_graph& exc)
    //    {
    //      chi::log.LogAllError()
    //        << "Function CheckCyclicDependencies. Detected cyclic depency.";
    //     chi::Exit(EXIT_FAILURE);
    //    }

    auto topological_order = G.GenerateTopologicalSort();
    if (topological_order.empty())
    {
      Chi::log.LogAllError() << "Function CheckCyclicDependencies. Detected cyclic depency.";
      Chi::Exit(EXIT_FAILURE);
    }

    // Cleanup
    G.Clear();
  } // for angles

  GetMeshStats();

  Chi::log.Log() << "Cyclic dependency check complete. No cycles or "
                 << "bad mesh elements were detected";
}

void
SurfaceMesh::GetMeshStats()
{
  std::vector<double> areas;
  std::vector<double> histo_bins;
  std::vector<int> histo;

  int num_negative_areas = 0;

  // Compute areas for each polyface
  size_t num_loc_cells = poly_faces_.size();
  areas.resize(num_loc_cells);
  double max_area = 0.0;
  for (size_t c = 0; c < num_loc_cells; c++)
  {
    PolyFace* face = poly_faces_[c];

    size_t num_edges = face->edges.size();
    double area = 0.0;
    for (size_t e = 0; e < num_edges; e++)
    {
      int v0i = face->edges[e][0];
      int v1i = face->edges[e][1];

      Vector3 v01 = vertices_[v1i] - vertices_[v0i];
      Vector3 v02 = face->face_centroid - vertices_[v0i];

      // This is essentially the combine of the triangle for each side

      area += 0.5 * (v01.x * v02.y - v01.y * v02.x);
    } // for edge

    areas[c] = area;
    if (area > max_area) max_area = area;

    if (area <= 0.0) num_negative_areas += 1;
  } // for cell

  // Sort the areas
  std::sort(areas.begin(), areas.end(), std::greater<double>());

  // Compute histogram bins
  histo_bins.resize(10);
  histo.resize(10, 0);
  histo_bins[0] = max_area * 1.05;
  for (int i = 1; i < 10; i++)
    histo_bins[i] = histo_bins[i - 1] / 2.0;

  // Polulate histogram
  for (auto area : areas)
  {
    int home_bin = 9;
    for (int i = 0; i < 10; i++)
    {

      if (area <= histo_bins[i]) home_bin = i;

    } // check bins

    histo[home_bin] += 1;
  } // for areas

  std::stringstream output;
  for (int i = 0; i < 10; i++)
  {
    char buff[100];
    snprintf(buff, 100, "%11.3e", histo_bins[i]);

    output << "Areas < " << buff << " = " << histo[i] << "\n";
  }
  output << "Number of negative or zero faces = " << num_negative_areas;

  Chi::log.Log() << output.str();
}

void
SurfaceMesh::ComputeLoadBalancing(std::vector<double>& x_cuts, std::vector<double>& y_cuts)
{
  Chi::log.Log() << "X-cuts to be logged: " << x_cuts.size();
  //  for (auto& val : x_cuts)
  //    chi::log.Log() << val;
  //
  Chi::log.Log() << "Y-cuts to be logged: " << y_cuts.size();
  //  for (auto& val : y_cuts)
  //    chi::log.Log() << val;

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
      if (poly_face->face_centroid.x >= x_cuts[i]) ref_i = i + 1;
    } // for i
    for (size_t j = 0; j < J; ++j)
    {
      if (poly_face->face_centroid.y >= y_cuts[j]) ref_j = j + 1;
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

  Chi::log.Log() << "Average faces per set: " << average;
  Chi::log.Log() << "Maximum faces per set: " << max_bin_size << " at (i,j)= ( " << i_max << " , "
                 << j_max << " )";

  if (i_max == I) Chi::log.Log() << "X greater than " << x_cuts[i_max - 1];
  else if (i_max == 0)
    Chi::log.Log() << "X less than " << x_cuts[0];
  else
    Chi::log.Log() << "X greater than " << x_cuts[i_max - 1] << " and less than " << x_cuts[i_max];

  if (j_max == J) Chi::log.Log() << "Y greater than " << y_cuts[j_max - 1];
  else if (j_max == 0)
    Chi::log.Log() << "Y less than " << y_cuts[0];
  else
    Chi::log.Log() << "Y greater than " << y_cuts[j_max - 1] << " and less than " << y_cuts[j_max];

  Chi::log.Log() << "Max-to-average ratio: " << max_bin_size / average;
}

void
SurfaceMesh::SplitByPatch(std::vector<SurfaceMesh*>& patches)
{
  typedef std::vector<Face> FaceList;
  typedef std::vector<FaceList*> FaceListCollection;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Copy all faces from surface
  FaceList unsorted_faces;
  FaceList::iterator cur_face;
  for (cur_face = this->faces_.begin(); cur_face != this->faces_.end(); cur_face++)
  {
    unsorted_faces.push_back((*cur_face));
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize co-planar
  // collection
  FaceListCollection co_planar_lists;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Seed the first collection
  FaceList* first_list = new FaceList;
  co_planar_lists.push_back(first_list);
  first_list->push_back(unsorted_faces.back());
  unsorted_faces.pop_back();

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reverse iterate over unsorted
  FaceList::reverse_iterator us_face;
  for (us_face = unsorted_faces.rbegin(); us_face != unsorted_faces.rend(); us_face++)
  {
    bool matchFound = false;

    // Check if it can be matched
    FaceListCollection::iterator existing_list;
    for (existing_list = co_planar_lists.begin(); existing_list != co_planar_lists.end();
         existing_list++)
    {
      for (cur_face = (*existing_list)->begin(); cur_face != (*existing_list)->end(); cur_face++)
      {
        Normal n1 = cur_face->geometric_normal;
        Normal n2 = us_face->geometric_normal;

        if (fabs(n1.Dot(n2)) > (1.0 - 1.0e-4))
        {
          (*existing_list)->push_back(unsorted_faces.back());
          unsorted_faces.pop_back();
          matchFound = true;
          break;
        }
      }
      if (matchFound) { break; }
    }

    // If no match found, create new list
    if (!matchFound)
    {
      printf("New list created\n");
      FaceList* new_list = new FaceList;
      new_list->push_back(unsorted_faces.back());
      unsorted_faces.pop_back();
      co_planar_lists.push_back(new_list);
    }
  }

  printf("Number of co-planar collections=%lu\n", co_planar_lists.size());

  FaceListCollection patch_list_collection;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loop over co_planar lists
  FaceListCollection::iterator existing_list;
  for (existing_list = co_planar_lists.begin(); existing_list != co_planar_lists.end();
       existing_list++)
  {
    // Inner patch collection
    FaceListCollection inner_patch_list_collection;

    // Add all the faces of this collection to an unused list
    FaceList unused_faces;
    for (cur_face = (*existing_list)->begin(); cur_face != (*existing_list)->end(); cur_face++)
    {
      unused_faces.push_back((*cur_face));
    }

    // Seed the first patch list
    FaceList* first_patch_list = new FaceList;
    inner_patch_list_collection.push_back(first_patch_list);
    first_patch_list->push_back(unused_faces.back());
    unused_faces.pop_back();

    // Loop over unused faces
    while (unused_faces.size() > 0)
    {
      // printf("Unused faces=%d\n",unused_faces.size());
      bool updateMade = false;
      FaceList::iterator us_face_f; // Forward iterator
      for (us_face_f = unused_faces.begin(); us_face_f != unused_faces.end(); us_face_f++)
      {
        // Try to to find a home
        FaceListCollection::iterator inner_list;
        for (inner_list = inner_patch_list_collection.begin();
             inner_list != inner_patch_list_collection.end();
             inner_list++)
        {
          for (cur_face = (*inner_list)->begin(); cur_face != (*inner_list)->end(); cur_face++)
          {
            // Check if any vertices match
            for (int e = 0; e < 3; e++)
            {
              int vi = (*us_face_f).v_index[e];
              for (int e2 = 0; e2 < 3; e2++)
              {
                int vf = (*cur_face).v_index[e2];

                if (vf == vi)
                {
                  (*inner_list)->push_back(*us_face_f);
                  unused_faces.erase(us_face_f);
                  updateMade = true;
                  break;
                }
              }
              if (updateMade) { break; }
            }
            if (updateMade) { break; }
          }
          if (updateMade) { break; }
        }

        if (updateMade) { break; }
      }

      if (!updateMade)
      {
        FaceList* new_patch_list = new FaceList;
        inner_patch_list_collection.push_back(new_patch_list);
        new_patch_list->push_back(unused_faces.back());
        unused_faces.pop_back();
      }
    }

    // Add inner patch lists to outer
    FaceListCollection::iterator inner_list;
    for (inner_list = inner_patch_list_collection.begin();
         inner_list != inner_patch_list_collection.end();
         inner_list++)
    {
      patch_list_collection.push_back(*inner_list);
    }
  }

  printf("Number of patches = %lu\n", patch_list_collection.size());

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create surfaces for each patch
  FaceListCollection::iterator outer_patch;
  for (outer_patch = patch_list_collection.begin(); outer_patch != patch_list_collection.end();
       outer_patch++)
  {

    SurfaceMesh* new_surface = new SurfaceMesh;

    std::vector<int*> vertex_mapping;

    for (cur_face = (*outer_patch)->begin(); cur_face != (*outer_patch)->end(); cur_face++)
    {
      // Copy the face
      Face newFace = (*cur_face);

      // Copy and map vertices
      for (int e = 0; e < 3; e++)
      {
        int vi = newFace.v_index[e];

        // Check if vertex already there
        bool already_there = false;
        int* already_there_mapping;
        std::vector<int*>::iterator cur_map;
        for (cur_map = vertex_mapping.begin(); cur_map != vertex_mapping.end(); cur_map++)
        {
          if ((*cur_map)[0] == vi)
          {
            already_there = true;
            already_there_mapping = (*cur_map);
            break;
          }
        }

        if (already_there)
        {
          // Update vertex index
          newFace.v_index[e] = already_there_mapping[1];
        }
        else
        {
          // Copy vertex
          Vertex v = this->vertices_.at(vi);
          int* newMapping = new int[2];
          newMapping[0] = vi;
          newMapping[1] = new_surface->vertices_.size();

          new_surface->vertices_.push_back(v);
          vertex_mapping.push_back(newMapping);

          newFace.v_index[e] = newMapping[1];
        }

      } // for e
      newFace.e_index[0][0] = newFace.v_index[0];
      newFace.e_index[0][1] = newFace.v_index[1];

      newFace.e_index[1][0] = newFace.v_index[1];
      newFace.e_index[1][1] = newFace.v_index[2];

      newFace.e_index[2][0] = newFace.v_index[2];
      newFace.e_index[2][1] = newFace.v_index[0];

      for (int e = 0; e < 3; e++)
      {
        newFace.e_index[e][2] = -1;
        newFace.e_index[e][3] = -1;
      }

      new_surface->faces_.push_back(newFace);
    }
    new_surface->UpdateInternalConnectivity();
    // std::cout << (*new_surface);
    patches.push_back(new_surface);

    //    std::string fileName = "ZZMesh";
    //    fileName = fileName + std::to_string((int)patches.size());
    //    fileName = fileName + ".obj";
    //
    //    std::cout <<fileName <<"\n";
    //
    //    new_surface->ExportToOBJFile(fileName.c_str());
  }
}

} // namespace chi_mesh
