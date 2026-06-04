// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/input_parameters.h"
#include <cstdio>
#include <vector>

namespace opensn
{

/**
 * Generic surface mesh class.
 * This class facilitates many functions within the mesh environment including logically determining
 * volumes.
 */
class SurfaceMesh
{
public:
  /// Structure containing edge properties
  struct Edge
  {
    /// Indices of the vertices
    std::array<int, 2> v_index{-1, -1};
    /// Indices of faces adjoining it
    std::array<int, 4> f_index{-1, -1, -1, -1};
    /// Vector vertices
    std::array<Vector3, 2> vertices{};
  };

  /// Data structure for a triangular face.
  struct Face
  {
    /// vertex indices
    std::array<int, 3> v_index{-1, -1, -1};
    /// normal indices
    std::array<int, 3> n_index{-1, -1, -1};
    /// vertex texture indices
    std::array<int, 3> vt_index{-1, -1, -1};
    /// edge indices
    std::array<std::array<int, 4>, 3> e_index{
      {{-1, -1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, -1, -1}}};

    Vector3 geometric_normal;
    Vector3 assigned_normal;
    Vector3 face_centroid;

    bool invalidated{false};
  };

  /**
   * Data structure for a polygon face.
   *
   * edges
   * An array of 4 integers.
   * [0] = Vertex index of edge start.
   * [1] = Vertex index of edge end.
   * [2] = Index of the face adjoining this edge (not the current face).
   *       -1 if not connected to anything,-1*boundary_index if connected
   *       to a boundary.
   * [3] = Edge number of adjoining face. -1 if not connected
   *       to anything. 0 if a boundary.
   *
   * face_indices
   *  [0] = Index of the adjoining cell. -1 if not connected to anything.
   *        -1*boundary_index if connected to a boundary.
   *  [1] = Face number of adjoining cell. -1 if not connected
   *       to anything. 0 if a boundary.
   *  [2] = Partition ID of adjacent cell.
   */
  struct PolyFace
  {
    std::vector<int> v_indices;
    std::vector<std::vector<int>> edges;
    std::array<int, 3> face_indices{{-1, -1, -1}};

    Vector3 geometric_normal;
    Vector3 face_centroid;

    bool invalidated{false};
  };

  explicit SurfaceMesh(const InputParameters& params);

  const std::vector<Vector3>& GetVertices() const { return vertices_; }

  const std::vector<Face>& GetTriangles() const { return faces_; }

  SurfaceMesh();
  ~SurfaceMesh();

  friend std::ostream& operator<<(std::ostream& os, SurfaceMesh& obj);

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
  static InputParameters GetInputParameters();
  static std::shared_ptr<SurfaceMesh> Create(const ParameterBlock& params);
};

} // namespace opensn
