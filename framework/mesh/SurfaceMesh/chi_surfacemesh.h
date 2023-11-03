#pragma once

#include <stdio.h>
#include <vector>

#include "../chi_mesh.h"

//###################################################################
/** Generic surface mesh class.
This class facilitates many functions within the mesh environment including
logically determining volumes.*/
class chi_mesh::SurfaceMesh
{
protected:
  std::vector<chi_mesh::Vertex> vertices_;
  std::vector<chi_mesh::Vertex> tex_vertices_; ///< Texture vertices
  std::vector<chi_mesh::Normal> normals_;
  std::vector<chi_mesh::Face> faces_;
  std::vector<chi_mesh::Edge> lines_;
  std::vector<chi_mesh::PolyFace*> poly_faces_; ///< Polygonal faces

  std::vector<int> physical_region_map_;

public:
  const std::vector<chi_mesh::Vertex>& GetVertices() const { return vertices_; }

  const std::vector<chi_mesh::Face>& GetTriangles() const { return faces_; }

  const std::vector<chi_mesh::PolyFace*>& GetPolygons() const { return poly_faces_; }

  // constrdestr.cc
  SurfaceMesh();
  ~SurfaceMesh();
  friend std::ostream& operator<<(std::ostream& os, SurfaceMesh& dt);
  // loadexport.cc
  int ImportFromOBJFile(const std::string& fileName,
                        bool as_poly = false,
                        const chi_mesh::Vector3& transform = Vector3(0, 0, 0));
  int ImportFromTriangleFiles(const char* fileName, bool as_poly);
  int ImportFromMshFiles(const char* fileName, bool as_poly);
  void ExportToOBJFile(const char* fileName);
  void ExportToPolyFile(const char* fileName);
  static SurfaceMesh* CreateFromDivisions(std::vector<double>& vertices_1d_x,
                                          std::vector<double>& vertices_1d_y);

  // internalconn.cc
  void UpdateInternalConnectivity();

  // checksense.cc
  bool CheckNegativeSense(double x, double y, double z);

  // splitbypatch.cc
  void SplitByPatch(std::vector<chi_mesh::SurfaceMesh*>& patches);

  // extractopenedges.cc
  void ExtractOpenEdgesToObj(const char* fileName);

  // meshstats.cc
  void CheckCyclicDependencies(int num_angles);
  void GetMeshStats();
  void ComputeLoadBalancing(std::vector<double>& x_cuts, std::vector<double>& y_cuts);
};


