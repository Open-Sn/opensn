#pragma once

#include "mesh/chi_mesh.h"
#include "mesh/Cell/cell.h"

namespace chi_mesh
{
namespace mesh_cutting
{
typedef std::pair<uint64_t, uint64_t> Edge;

/**A general structure storing the two vertex-ids
 * to be cut or being cut. Additionally a cut-point
 * id is stored.*/
struct ECI // EdgeCutInfo
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

// 2D_utils
Edge MakeUniqueEdge(const Edge& edge);
Edge MakeEdgeFromPolygonEdgeIndex(const std::vector<uint64_t>& vertex_ids, size_t edge_index);
chi_mesh::Vector3 GetEdgeCentroid(const Edge& edge, const chi_mesh::MeshContinuum& grid);

void PopulatePolygonFromVertices(const MeshContinuum& mesh,
                                 const std::vector<uint64_t>& vertex_ids,
                                 chi_mesh::Cell& cell);

bool CheckPolygonQuality(const MeshContinuum& mesh,
                         const chi_mesh::Cell& cell,
                         bool check_convexity = false);

// 2D_cutcell
void CutPolygon(const std::vector<ECI>& cut_edges,
                const std::set<uint64_t>& cut_vertices,
                const Vector3& plane_point,
                const Vector3& plane_normal,
                MeshContinuum& mesh,
                chi_mesh::Cell& cell);

// 3D_utils
bool CheckPolyhedronQuality(const MeshContinuum& mesh,
                            const chi_mesh::Cell& cell,
                            bool check_convexity = false);

std::set<size_t> FindNeighborFaceIndices(const std::vector<std::vector<uint64_t>>& proxy_faces,
                                         size_t face_index);

std::vector<Edge> FindNonManifoldEdges(const std::vector<std::vector<uint64_t>>& proxy_faces);

std::vector<Edge> StitchEdgesEndToEnd(const std::vector<Edge>& edges);

void PopulatePolyhedronFromFaces(const MeshContinuum& mesh,
                                 const std::vector<std::vector<uint64_t>>& raw_faces,
                                 chi_mesh::Cell& cell);

// 3D_cutcell
void Cut3DCell(const std::vector<ECI>& global_cut_edges,
               const std::set<uint64_t>& number,
               const Vector3& plane_point,
               const Vector3& plane_normal,
               double float_compare,
               MeshContinuum& mesh,
               chi_mesh::Cell& cell,
               bool verbose = false);

// plane
void CutMeshWithPlane(MeshContinuum& mesh,
                      const Vector3& plane_point,
                      const Vector3& plane_normal,
                      double merge_tolerance = 1.0e-3,
                      double float_compare = 1.0e-10);
} // namespace mesh_cutting
} // namespace chi_mesh

