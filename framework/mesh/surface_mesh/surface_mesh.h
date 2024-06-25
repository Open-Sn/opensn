// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <stdio.h>
#include <vector>
#include "framework/mesh/mesh.h"
#include "framework/object.h"

namespace opensn
{

/**
 * Generic surface mesh class.
 * This class facilitates many functions within the mesh environment including logically
 * determining volumes.
 */
class SurfaceMesh : public Object
{
public:
  static InputParameters GetInputParameters();

protected:
  explicit SurfaceMesh(const InputParameters& params);

protected:
  std::vector<Vertex> vertices_;
  std::vector<Vertex> tex_vertices_; ///< Texture vertices
  std::vector<Normal> normals_;
  std::vector<Face> faces_;
  std::vector<Edge> lines_;
  std::vector<std::shared_ptr<PolyFace>> poly_faces_; ///< Polygonal faces

  std::vector<int> physical_region_map_;

public:
  const std::vector<Vertex>& GetVertices() const { return vertices_; }

  const std::vector<Face>& GetTriangles() const { return faces_; }

  const std::vector<std::shared_ptr<PolyFace>>& GetPolygons() const { return poly_faces_; }

  SurfaceMesh();
  ~SurfaceMesh();
  friend std::ostream& operator<<(std::ostream& os, SurfaceMesh& dt);

  /**
   * Loads a surface mesh from a wavefront .obj file.
   */
  int ImportFromOBJFile(const std::string& fileName,
                        bool as_poly = false,
                        const Vector3& transform = Vector3(0, 0, 0));

  /**
   * Loads a surface mesh from triangle's file format.
   */
  int ImportFromTriangleFiles(const char* fileName, bool as_poly);

  /**
   * Loads a surface mesh from gmsh's file format.
   */
  int ImportFromMshFiles(const char* fileName, bool as_poly);

  /**
   * Exports the triangular faces of a surface mesh to
   * wavefront .obj files.
   */
  void ExportToOBJFile(const char* fileName);

  /**
   * Exports a PSLG to triangle1.6's .poly format.
   */
  void ExportToPolyFile(const char* fileName);

  /**
   * Creates a 2D orthogonal mesh from a set of vertices in x and y.
   * The vertices along a dimension merely represents the divisions. They
   * are not the complete vertices defining a cell. For example:
   * \code
   * std::vector<Vertex> vertices_x = {0.0,1.0,2.0};
   * std::vector<Vertex> vertices_y = {0.0,1.0,2.0};
   * SurfaceMesh::CreateFromDivisions(vertices_x,vertices_y);
   * \endcode
   *
   * This code will create a 2x2 mesh with \f$ \vec{x} \in [0,2]^2 \f$.
   */
  static std::shared_ptr<SurfaceMesh> CreateFromDivisions(std::vector<double>& vertices_1d_x,
                                                          std::vector<double>& vertices_1d_y);

  /**
   * Runs over the faces of the surfacemesh and determines
   * neighbors. The algorithm first establishes which cells subscribe to each
   * vertex and then loops over faces and edges. For each edge, only the
   * subscribing faces are searched for neighbors. This routine has
   * time complexity O(N).
   */
  void UpdateInternalConnectivity();

  bool CheckNegativeSense(double x, double y, double z);

  /**
   * Splits the surface by patch.
   */
  void SplitByPatch(std::vector<std::shared_ptr<SurfaceMesh>>& patches);

  /**
   * Extract open edges to wavefront obj format.
   */
  void ExtractOpenEdgesToObj(const char* fileName);

  /**
   * Checks for cyclic dependencies in this mesh.
   * Transport type sweeps have a step where the inter-cell dependency is acyclically sorted.
   * This step is repeated here.
   */
  void CheckCyclicDependencies(int num_angles);

  /**
   * Gets simple mesh statistics.
   */
  void GetMeshStats();

  /**
   * Computes load balancing parameters from a set of predictive cuts.
   * Does not actually perform these cuts.
   */
  void ComputeLoadBalancing(std::vector<double>& x_cuts, std::vector<double>& y_cuts);
};

} // namespace opensn
