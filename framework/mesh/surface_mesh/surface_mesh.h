// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdio>
#include <vector>
#include "framework/mesh/mesh.h"
#include "framework/object.h"

namespace opensn
{

/**
 * Generic surface mesh class.
 * This class facilitates many functions within the mesh environment including logically determining
 * volumes.
 */
class SurfaceMesh : public Object
{
public:
  static InputParameters GetInputParameters();

protected:
  explicit SurfaceMesh(const InputParameters& params);

protected:
  std::vector<Vector3> vertices_;
  /// Texture vertices
  std::vector<Vector3> tex_vertices_;
  std::vector<Vector3> normals_;
  std::vector<Face> faces_;
  std::vector<Edge> lines_;
  /// Polygonal faces
  std::vector<std::shared_ptr<PolyFace>> poly_faces_;

  std::vector<int> physical_region_map_;

public:
  const std::vector<Vector3>& GetVertices() const { return vertices_; }

  const std::vector<Face>& GetTriangles() const { return faces_; }

  SurfaceMesh();
  ~SurfaceMesh() override;

  friend std::ostream& operator<<(std::ostream& os, SurfaceMesh& dt);

  /// Loads a surface mesh from a wavefront .obj file.
  int ImportFromOBJFile(const std::string& fileName,
                        bool as_poly = false,
                        const Vector3& transform = Vector3(0, 0, 0));

  /// Loads a surface mesh from triangle's file format.
  int ImportFromTriangleFiles(const char* fileName, bool as_poly);

  /// Loads a surface mesh from gmsh's file format.
  int ImportFromMshFiles(const char* fileName, bool as_poly);

  /**
   * Runs over the faces of the surfacemesh and determines neighbors. The algorithm first
   * establishes which cells subscribe to each vertex and then loops over faces and edges. For each
   * edge, only the subscribing faces are searched for neighbors. This routine has time complexity
   * O(N).
   */
  void UpdateInternalConnectivity();

  /**
   * Computes load balancing parameters from a set of predictive cuts.
   * Does not actually perform these cuts.
   */
  void ComputeLoadBalancing(std::vector<double>& x_cuts, std::vector<double>& y_cuts);
};

} // namespace opensn
