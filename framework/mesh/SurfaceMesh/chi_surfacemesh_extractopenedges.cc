#include "chi_surfacemesh.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <set>
#include <fstream>

//###################################################################
/**Extract open edges to wavefront obj format.*/
void
chi_mesh::SurfaceMesh::ExtractOpenEdgesToObj(const char* fileName)
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
    Chi::log.LogAllError() << "In call to chi_mesh::SurfaceMesh::ExtractOpenEdgesToObj. Failed"
                           << " to open file: " << std::string(fileName);
    Chi::Exit(EXIT_FAILURE);
  }

  outfile << "# ChiTech open edges file\n";
  outfile << "# Single surface mesh\n";

  for (auto vert_pair : edges)
  {
    chi_mesh::Vertex& v0 = vertices_[vert_pair.first];
    chi_mesh::Vertex& v1 = vertices_[vert_pair.second];
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