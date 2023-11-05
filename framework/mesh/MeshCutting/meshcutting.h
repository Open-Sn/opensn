#pragma once

#include "opensn/framework/mesh/chi_mesh.h"
#include "opensn/framework/mesh/Cell/cell.h"

namespace chi_mesh
{
namespace mesh_cutting
{

typedef std::pair<uint64_t, uint64_t> Edge;

/**
 * A general structure storing the two vertex-ids
 * to be cut or being cut. Additionally a cut-point
 * id is stored.
 */
struct ECI
{
  Edge vertex_ids;
  uint64_t cut_point_id = 0;

  ECI() = default;

  explicit ECI(Edge in_edge, uint64_t in_cutpoint_id)
    : vertex_ids(std::move(in_edge)), cut_point_id(in_cutpoint_id)
  {
  }

  static bool Comparator(const ECI& edge_cut_info, const Edge& ref_edge)
  {
    return ref_edge == edge_cut_info.vertex_ids;
  }
};

/**
 * Makes a unique edge from a regular edge. A unique edge is an edge with its 1st vertex-id the
 * smallest of the two vertex-ids.
 */
Edge MakeUniqueEdge(const Edge& edge);

/**
 * Make an edge for a polygon given its edge index.
 */
Edge MakeEdgeFromPolygonEdgeIndex(const std::vector<uint64_t>& vertex_ids, size_t edge_index);

/**
 * Computes the centroid of an edge.
 */
chi_mesh::Vector3 GetEdgeCentroid(const Edge& edge, const chi_mesh::MeshContinuum& grid);

void PopulatePolygonFromVertices(const MeshContinuum& mesh,
                                 const std::vector<uint64_t>& vertex_ids,
                                 chi_mesh::Cell& cell);

/**
 * Performs a quality check of a given polygon. The simple quality requirement
 * is that, if we form a triangle with the cell-centroid and the vertices of
 * each edge (in ccw orientation), no inverted triangles are present.
 */
bool CheckPolygonQuality(const MeshContinuum& mesh,
                         const chi_mesh::Cell& cell,
                         bool check_convexity = false);

/**
 * Cuts a polygon.
 */
void CutPolygon(const std::vector<ECI>& cut_edges,
                const std::set<uint64_t>& cut_vertices,
                const Vector3& plane_point,
                const Vector3& plane_normal,
                MeshContinuum& mesh,
                chi_mesh::Cell& cell);

/**
 * Checks the quality of a polyhedron against the following.
 *
 * Take each face and form triangles using the face edges and the face
 * centroid. Now take each triangle and form a tetrahedron with the
 * cell-centroid as the 4th vertex. The quality requirement for this tetrahedron
 * is that its primary face must not be inverted with respect to the
 * cell-centroid.
 *
 * Another optional requirement is that the cell must be convex. This is checked
 * by using a face centroid CFC and any neighboring face's centroid NFC and
 * normal, n_N. If the dot-product of (NFC-CFC) with n_N is negative, then the
 * cell cannot be classified as convex.
 */
bool CheckPolyhedronQuality(const MeshContinuum& mesh,
                            const chi_mesh::Cell& cell,
                            bool check_convexity = false);

/**
 * Returns the face-indices that have adjacent edges to face with index face_index.
 */
std::set<size_t> FindNeighborFaceIndices(const std::vector<std::vector<uint64_t>>& proxy_faces,
                                         size_t face_index);

/**
 * Finds the non-manifold edges of a collection of faces.
 */
std::vector<Edge> FindNonManifoldEdges(const std::vector<std::vector<uint64_t>>& proxy_faces);

/**
 * Attempts to stitch edges end-to-end.
 */
std::vector<Edge> StitchEdgesEndToEnd(const std::vector<Edge>& edges);

/**
 * Defines a polyhedron based only on its faces.
 */
void PopulatePolyhedronFromFaces(const MeshContinuum& mesh,
                                 const std::vector<std::vector<uint64_t>>& raw_faces,
                                 chi_mesh::Cell& cell);

void Cut3DCell(const std::vector<ECI>& global_cut_edges,
               const std::set<uint64_t>& number,
               const Vector3& plane_point,
               const Vector3& plane_normal,
               double float_compare,
               MeshContinuum& mesh,
               chi_mesh::Cell& cell,
               bool verbose = false);

/**
 * Cuts a mesh with a plane.
 */
void CutMeshWithPlane(MeshContinuum& mesh,
                      const Vector3& plane_point,
                      const Vector3& plane_normal,
                      double merge_tolerance = 1.0e-3,
                      double float_compare = 1.0e-10);

} // namespace mesh_cutting
} // namespace chi_mesh
