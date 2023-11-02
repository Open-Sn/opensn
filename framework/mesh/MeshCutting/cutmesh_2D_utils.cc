#include "meshcutting.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Makes a unique edge from a regular edge. A unique edge is an edge
 * with its 1st vertex-id the smallest of the two vertex-ids.*/
chi_mesh::mesh_cutting::Edge
chi_mesh::mesh_cutting::MakeUniqueEdge(const Edge& edge)
{
  return std::make_pair(std::min(edge.first, edge.second), std::max(edge.first, edge.second));
}

//###################################################################
/**Make an edge for a polygon given its edge index.*/
std::pair<uint64_t, uint64_t>
chi_mesh::mesh_cutting::MakeEdgeFromPolygonEdgeIndex(const std::vector<uint64_t>& vertex_ids,
                                                     size_t edge_index)
{
  const int e = static_cast<int>(edge_index);
  const int num_verts = static_cast<int>(vertex_ids.size());

  int next_v = (e < (num_verts - 1)) ? e + 1 : 0;
  uint64_t v0_id = vertex_ids[e];
  uint64_t v1_id = vertex_ids[next_v];

  return std::make_pair(v0_id, v1_id);
}

//###################################################################
/**Computes the centroid of an edge.*/
chi_mesh::Vector3
chi_mesh::mesh_cutting::GetEdgeCentroid(const Edge& edge, const chi_mesh::MeshContinuum& grid)
{
  auto& v0 = grid.vertices[edge.first];
  auto& v1 = grid.vertices[edge.second];

  return 0.5 * (v0 + v1);
}

//###################################################################
/***/
void
chi_mesh::mesh_cutting::PopulatePolygonFromVertices(const MeshContinuum& mesh,
                                                    const std::vector<uint64_t>& vertex_ids,
                                                    chi_mesh::Cell& cell)
{
  cell.faces_.clear();
  cell.faces_.reserve(vertex_ids.size());
  cell.vertex_ids_ = vertex_ids;

  cell.centroid_ = chi_mesh::Vector3(0.0, 0.0, 0.0);
  for (uint64_t vid : cell.vertex_ids_)
    cell.centroid_ += mesh.vertices[vid];
  cell.centroid_ /= double(cell.vertex_ids_.size());

  size_t num_verts = vertex_ids.size();
  for (size_t v = 0; v < num_verts; ++v)
  {
    size_t v1_ref = (v < (num_verts - 1)) ? v + 1 : 0;

    uint64_t v0id = cell.vertex_ids_[v];
    uint64_t v1id = cell.vertex_ids_[v1_ref];

    const auto& v0 = mesh.vertices[v0id];
    const auto& v1 = mesh.vertices[v1id];

    chi_mesh::Vector3 v01 = v1 - v0;

    chi_mesh::CellFace face;
    face.vertex_ids_ = {v0id, v1id};
    face.normal_ = v01.Cross(chi_mesh::Normal(0.0, 0.0, 1.0)).Normalized();
    face.centroid_ = (v0 + v1) / 2.0;

    cell.faces_.push_back(face);
  }
}

//###################################################################
/**Performs a quality check of a given polygon. The simple quality requirement
 * is that, if we form a triangle with the cell-centroid and the vertices of
 * each edge (in ccw orientation), no inverted triangles are present.*/
bool
chi_mesh::mesh_cutting::CheckPolygonQuality(const MeshContinuum& mesh,
                                            const chi_mesh::Cell& cell,
                                            const bool check_convexity /*=false*/)
{
  const chi_mesh::Vector3 khat(0.0, 0.0, 1.0);

  auto& v0 = cell.centroid_;
  size_t num_edges = cell.vertex_ids_.size();

  for (size_t e = 0; e < num_edges; ++e)
  {
    auto edge = MakeEdgeFromPolygonEdgeIndex(cell.vertex_ids_, e);

    const auto& v1 = mesh.vertices[edge.first];
    const auto& v2 = mesh.vertices[edge.second];

    auto v01 = v1 - v0;
    auto v02 = v2 - v0;

    if (v01.Cross(v02).Dot(khat) < 0.0) return false;
  } // for edge

  //============================================= Optional convexity check
  if (check_convexity) {}

  return true;
}
