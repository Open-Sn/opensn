#include "meshcutting.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include <algorithm>
#include <queue>

//###################################################################
/**Checks the quality of a polyhedron against the following.
 *
 * Take each face and form triangles using the face edges and the face
 * centroid. Now take each triangle and form a tetrahedron with the
 * cell-centroid as the 4th vertex. The quality requirement for this tetrahedron
 * is that its primary face must not be inverted with respect to the
 * cell-centroid.
 *
 * Another optional requirement is that the cell must be convex. This is checked
 * by using a face centroid CFC and any neighboring face's centroid NFC and normal, n_N. If
 * the dot-product of (NFC-CFC) with n_N is negative, then the
 * cell cannot be classified as convex.*/
bool chi_mesh::mesh_cutting::
  CheckPolyhedronQuality(const MeshContinuum &mesh,
                         const chi_mesh::Cell& cell,
                         const bool check_convexity/*=false*/)
{
  const auto& C = cell.centroid_;

  for (auto& face : cell.faces_)
  {
    const auto& v0 = face.centroid_;
    const size_t num_edges = face.vertex_ids_.size();

    for (size_t e=0; e<num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids_, e);

      const auto& v1 = mesh.vertices[edge.first];
      const auto& v2 = mesh.vertices[edge.second];

      auto v01 = v1-v0;
      auto v02 = v2-v0;

      auto n = v01.Cross(v02);
      auto C_tri = (v0+v1+v2)/3.0;

      auto CC = C_tri - C;

      if (CC.Dot(n)<0.0) return false;
    }
  }//for face

  //============================================= Optional convexity check
  if (check_convexity)
  {
    std::vector<std::vector<uint64_t>> proxy_faces;
    for (const auto& face : cell.faces_)
      proxy_faces.push_back(face.vertex_ids_);

    size_t f=0;
    for (const auto& face : cell.faces_)
    {
      std::set<size_t> neighbor_face_indices =
        FindNeighborFaceIndices(proxy_faces, f);

      for (size_t ofi : neighbor_face_indices)
      {
        auto& other_face = cell.faces_[ofi];
        auto CFC_OFC = other_face.centroid_ - face.centroid_;

        if (CFC_OFC.Dot(other_face.normal_) < 0.0)
          return false;
      }
      ++f;
    }//for f
  }

  return true;
}


//###################################################################
/**Returns the face-indices that have adjacent edges to face with
 * index face_index.*/
std::set<size_t> chi_mesh::mesh_cutting::
  FindNeighborFaceIndices(const std::vector<std::vector<uint64_t>> &proxy_faces,
                          const size_t face_index)
{
  const auto& face = proxy_faces.at(face_index);
  const size_t num_faces = proxy_faces.size();

  std::set<size_t> neighbor_face_indices;
  const size_t face_num_edges = face.size();
  for (size_t e=0; e<face_num_edges; ++e)
  {
    auto face_edge = MakeEdgeFromPolygonEdgeIndex(face,e);
    auto face_uniq_edge = MakeUniqueEdge(face_edge);

    for (size_t of=0; of<num_faces; ++of)//other face
    {
      if (of == face_index) continue;
      const auto& other_face = proxy_faces[of];

      const size_t other_face_num_edges = other_face.size();
      for (size_t ofe=0; ofe<other_face_num_edges; ++ofe)//other face edge
      {
        auto other_face_edge = MakeEdgeFromPolygonEdgeIndex(other_face, ofe);
        auto other_face_uniq_edge = MakeUniqueEdge(other_face_edge);

        if (face_uniq_edge == other_face_uniq_edge)
        {
          neighbor_face_indices.insert(of);
          break; //from ofe-loop
        }
      }//for ofe
    }//for of
  }//for face edge

  return neighbor_face_indices;
}


//###################################################################
/**Finds the non-manifold edges of a collection of faces.*/
std::vector<chi_mesh::mesh_cutting::Edge> chi_mesh::mesh_cutting::
  FindNonManifoldEdges(const std::vector<std::vector<uint64_t>>& proxy_faces)
{
  std::vector<Edge> non_manifold_edges;
  //The reserve amount below is chosen as follows:
  //If we cut a tet the number of non-manifold edges will be max 3.
  //A Hex will be 4. Arbitrary polyhedra as produced by STAR-CCM+ are
  //on average hexagons so lets fluff this by an additional 4 edges
  //and call 10 edges a good estimate.
  non_manifold_edges.reserve(10);

  //=================================== Build edges for each proxy face
  // This saves us some readability later on
  const size_t num_proxy_faces = proxy_faces.size();
  std::vector<std::vector<Edge>> faces_edges;
  faces_edges.reserve(num_proxy_faces);
  for (const auto& proxy_face : proxy_faces)
  {
    std::vector<Edge> face_edges(proxy_face.size());
    for (size_t e=0; e<proxy_face.size(); ++e)
      face_edges[e] = MakeEdgeFromPolygonEdgeIndex(proxy_face, e);

    faces_edges.push_back(std::move(face_edges));
  }//for proxy face

  //=================================== Search each proxy face edge for
  //                                    being attached to an edge of any
  //                                    other proxy face
  for (size_t cf=0; cf<num_proxy_faces; ++cf)
  {
    for (const auto& cur_edge : faces_edges[cf])
    {
      bool is_manifold = false;
      for (size_t af=0; af<num_proxy_faces; ++af)
      {
        if (af == cf) continue;
        for (const auto& alt_edge : faces_edges[af])
        {
          if (MakeUniqueEdge(cur_edge) == MakeUniqueEdge(alt_edge))
          {
            is_manifold = true;
            goto before_next_edge;
          }
        }//for alternate edge
      }//for alternate face

      before_next_edge:
        if (not is_manifold) {non_manifold_edges.push_back(cur_edge);}
    }//for current edge
  }//for current face

  return non_manifold_edges;
}


//###################################################################
/**Attempts to stitch edges end-to-end.*/
std::vector<chi_mesh::mesh_cutting::Edge> chi_mesh::mesh_cutting::
  StitchEdgesEndToEnd(const std::vector<Edge>& edges)
{
  const std::string fname = __FUNCTION__;
  std::vector<Edge> stitched_nm_edges;
  stitched_nm_edges.reserve(edges.size());

  //=================================== Declaring two queues used to
  //                                    stitch
  std::queue<Edge> unused_candidate_edges;
  std::queue<Edge> unused_noncandidate_edges;

  // All the edges are initially candidates
  for (auto& edge : edges)
    unused_candidate_edges.push(edge);

  // Begin the stitch using the first candidate
  // edge. This edge is no longer part of the game
  stitched_nm_edges.push_back(unused_candidate_edges.front());
  unused_candidate_edges.pop();

  //=================================== Attempt to add all candidate edges
  //                                    to the stitch
  const size_t num_candidates = unused_candidate_edges.size();
  for (size_t k=0; k<num_candidates; ++k)
  {
    const auto& edge_previously_added = stitched_nm_edges.back();

    // The following while loop is gauranteed to terminate
    // because we pop the first element each time.
    while (not unused_candidate_edges.empty())
    {
      const auto& edge = unused_candidate_edges.front();

      if (edge.first == edge_previously_added.second)
        stitched_nm_edges.push_back(edge);
      else
        unused_noncandidate_edges.push(edge);

      unused_candidate_edges.pop(); //removes first element
    }

    std::swap(unused_noncandidate_edges, unused_candidate_edges);
    if (unused_candidate_edges.empty()) break;
  }//for k

  if (stitched_nm_edges.size() != edges.size())
    throw std::logic_error(fname + ": stitching failed.");

  return stitched_nm_edges;
}

//###################################################################
/**Defines a polyhedron based only on its faces.*/
void chi_mesh::mesh_cutting::
  PopulatePolyhedronFromFaces(const MeshContinuum &mesh,
                              const std::vector<std::vector<uint64_t>> &raw_faces,
                              chi_mesh::Cell &cell)
{
  std::vector<uint64_t> cell_vertex_ids;

  /**Lamda to add a vertex to a list.*/
  auto AddIDToCellVertexIDs = [&cell_vertex_ids](uint64_t new_id)
  {
    auto find_result = std::find(cell_vertex_ids.begin(),
                                 cell_vertex_ids.end(),
                                 new_id);

    if (find_result == cell_vertex_ids.end())
      cell_vertex_ids.push_back(new_id);
  };

  //======================================== Build faces
  cell.faces_.clear();
  cell.faces_.reserve(raw_faces.size());
  for (auto& raw_face : raw_faces)
  {
    chi_mesh::Vector3 face_centroid;
    for (uint64_t vid : raw_face)
    {
      face_centroid += mesh.vertices[vid];
      AddIDToCellVertexIDs(vid);
    }
    face_centroid /= static_cast<double>(raw_face.size());

    chi_mesh::CellFace new_face;
    new_face.centroid_ = face_centroid;
    new_face.vertex_ids_ = raw_face;

    cell.faces_.push_back(std::move(new_face));
  }//for raw face

  //======================================== Compute cell centroid
  chi_mesh::Vector3 cell_centroid;
  for (uint64_t vid : cell_vertex_ids)
    cell_centroid += mesh.vertices[vid];
  cell_centroid /= static_cast<double>(cell_vertex_ids.size());

  cell.centroid_ = cell_centroid;
  cell.vertex_ids_ = cell_vertex_ids;
}