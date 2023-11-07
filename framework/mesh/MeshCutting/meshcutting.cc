#include "framework/mesh/MeshCutting/meshcutting.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/mesh/Raytrace/raytracer.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include <algorithm>
#include <functional>
#include <queue>

namespace chi_mesh::mesh_cutting
{

Edge
MakeUniqueEdge(const Edge& edge)
{
  return std::make_pair(std::min(edge.first, edge.second), std::max(edge.first, edge.second));
}

std::pair<uint64_t, uint64_t>
MakeEdgeFromPolygonEdgeIndex(const std::vector<uint64_t>& vertex_ids, size_t edge_index)
{
  const int e = static_cast<int>(edge_index);
  const int num_verts = static_cast<int>(vertex_ids.size());

  int next_v = (e < (num_verts - 1)) ? e + 1 : 0;
  uint64_t v0_id = vertex_ids[e];
  uint64_t v1_id = vertex_ids[next_v];

  return std::make_pair(v0_id, v1_id);
}

chi_mesh::Vector3
GetEdgeCentroid(const Edge& edge, const chi_mesh::MeshContinuum& grid)
{
  auto& v0 = grid.vertices[edge.first];
  auto& v1 = grid.vertices[edge.second];

  return 0.5 * (v0 + v1);
}

void
PopulatePolygonFromVertices(const MeshContinuum& mesh,
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

bool
CheckPolygonQuality(const MeshContinuum& mesh,
                    const chi_mesh::Cell& cell,
                    const bool check_convexity)
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

  // Optional convexity check
  if (check_convexity) {}

  return true;
}

void
CutPolygon(const std::vector<ECI>& cut_edges,
           const std::set<uint64_t>& cut_vertices,
           const Vector3& plane_point,
           const Vector3& plane_normal,
           MeshContinuum& mesh,
           chi_mesh::Cell& cell)
{
  const auto& p = plane_point;
  const auto& n = plane_normal;

  /**Utility lambda to check if a vertex is in "cut_vertices" list.*/
  auto VertexIsCut = [&cut_vertices](uint64_t vid)
  {
    auto result = cut_vertices.find(vid);

    if (result != cut_vertices.end()) return true;

    return false;
  };

  /**Utility function to check if an edge is in the "cut_edges" list.*/
  auto EdgeIsCut = [&cut_edges](const Edge& edge)
  {
    Edge edge_set(std::min(edge.first, edge.second), std::max(edge.first, edge.second));

    constexpr auto Arg1 = std::placeholders::_1;
    constexpr auto Comparator = &ECI::Comparator;

    auto result =
      std::find_if(cut_edges.begin(), cut_edges.end(), std::bind(Comparator, Arg1, edge_set));

    if (result != cut_edges.end()) return std::make_pair(true, *result);

    return std::make_pair(false, *result);
  };

  // Create and set vertex and edge cut flags for the current cell. Also populate edge_cut_info
  size_t num_verts = cell.vertex_ids_.size();
  size_t num_edges = num_verts;

  std::vector<bool> vertex_cut_flags(num_verts, false);
  std::vector<bool> edge_cut_flags(num_edges, false);
  std::vector<ECI> edge_cut_info(num_edges);

  for (size_t e = 0; e < num_edges; ++e)
  {
    vertex_cut_flags[e] = VertexIsCut(cell.vertex_ids_[e]);

    auto edge = MakeEdgeFromPolygonEdgeIndex(cell.vertex_ids_, e);
    auto cut_nature = EdgeIsCut(edge);
    edge_cut_flags[e] = cut_nature.first;
    if (cut_nature.first) edge_cut_info[e] = cut_nature.second;

    edge_cut_info[e].vertex_ids = edge;
  } // populate flags

  // Lamda for edge loop
  enum class CurVertex
  {
    AT_FIRST,
    AT_CUT_POINT,
    AT_SECOND,
    NONE
  };

  struct CurCutInfo
  {
    int which_edge = 0;
    CurVertex which_vertex = CurVertex::AT_FIRST;

    CurCutInfo() = default;

    CurCutInfo(int in_which_edge, CurVertex in_which_vertex)
      : which_edge(in_which_edge), which_vertex(in_which_vertex)
    {
    }
  };

  /**This lamda function starts from a current cut-edge, which is either
   * an edge where the first vertex is cut or an edge that is cut
   * somewhere along its length, and then follows the edges in a ccw fashion
   * until it finds another cut. This last cut is just as either an edge
   * cut along its length or cut at the second vertex. This then completes
   * an edge loop that can be used to define another polygon.*/
  auto GetVerticesTillNextCut =
    [&cell, &edge_cut_flags, &edge_cut_info, &VertexIsCut](CurCutInfo start_cut_info)
  {
    size_t num_verts = cell.vertex_ids_.size();
    std::vector<uint64_t> vertex_ids;
    vertex_ids.reserve(num_verts); // Ought to be more than enough

    int e = start_cut_info.which_edge;
    auto end_type = CurVertex::NONE;

    switch (start_cut_info.which_vertex)
    {
      case CurVertex::AT_FIRST:
      {
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.first);
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);

        if (VertexIsCut(edge_cut_info[e].vertex_ids.second))
        {
          end_type = CurVertex::AT_SECOND;
          goto skip_to_return_portion;
        }

        break;
      }
      case CurVertex::AT_CUT_POINT:
      {
        vertex_ids.push_back(edge_cut_info[e].cut_point_id);
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);

        if (VertexIsCut(edge_cut_info[e].vertex_ids.second))
        {
          end_type = CurVertex::AT_SECOND;
          goto skip_to_return_portion;
        }

        break;
      }
      case CurVertex::AT_SECOND:
      case CurVertex::NONE:
      default:
        break;
    } // switch type of starting cut

    // Look at downstream ccw edges and check for
    // edges cut or end-point cuts
    for (int eref = 0; eref < num_verts; ++eref)
    {
      e = (e < (num_verts - 1)) ? e + 1 : 0;

      if (e == start_cut_info.which_edge) break;

      if (edge_cut_flags[e])
      {
        vertex_ids.push_back(edge_cut_info[e].cut_point_id);
        end_type = CurVertex::AT_CUT_POINT;
        break;
      }
      else if (VertexIsCut(edge_cut_info[e].vertex_ids.second))
      {
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);
        end_type = CurVertex::AT_SECOND;
        break;
      }
      else
        vertex_ids.push_back(edge_cut_info[e].vertex_ids.second);
    } // for eref

  skip_to_return_portion:
    CurCutInfo end_cut_info(e, end_type);

    return std::make_pair(vertex_ids, end_cut_info);
  };

  typedef std::pair<std::vector<uint64_t>, CurCutInfo> LoopInfo;

  // Process all edges and create edge loops with associated cells
  std::vector<std::vector<uint64_t>> loops_to_add_to_mesh;
  std::vector<CurCutInfo> cut_history;
  for (size_t e = 0; e < num_edges; ++e)
  {
    LoopInfo loop_info;

    if (vertex_cut_flags[e]) loop_info = GetVerticesTillNextCut(CurCutInfo(e, CurVertex::AT_FIRST));
    else if (edge_cut_flags[e])
      loop_info = GetVerticesTillNextCut(CurCutInfo(e, CurVertex::AT_CUT_POINT));
    else
      continue;

    std::vector<uint64_t> verts_to_next_cut = loop_info.first;
    int end_edge = loop_info.second.which_edge;
    CurVertex end_type = loop_info.second.which_vertex;

    // Notes:
    //  - If the end_edge == e then this means trouble as a polygon cannot
    //    be cut like that.
    //  - If end_edge < e then the end-edge is definitely the edge right before
    //    the first cut edge. We should still process this edge-loop, but stop
    //    searching.

    if (end_edge < e)         // if looped past num_edges
      e = int(num_edges) - 1; // stop search
    else if (end_type == CurVertex::AT_SECOND)
      e = end_edge; // resume search after end_e
    else
      e = end_edge - 1; // resume search at end_e. e will get ++ to end_edge

    loops_to_add_to_mesh.push_back(verts_to_next_cut);
  } // for e

  // Add derivative cells to mesh

  // Take the back cell and paste onto
  // the current cell reference
  if (not loops_to_add_to_mesh.empty())
  {
    auto& back_loop = loops_to_add_to_mesh.back();
    PopulatePolygonFromVertices(mesh, back_loop, cell);

    loops_to_add_to_mesh.pop_back();
  } // if there are cells to add

  // Now push-up new cells from the remainder
  for (auto& loop : loops_to_add_to_mesh)
  {
    auto new_cell = std::make_unique<chi_mesh::Cell>(CellType::POLYGON, CellType::TRIANGLE);
    PopulatePolygonFromVertices(mesh, loop, *new_cell);
    mesh.cells.push_back(std::move(new_cell));
  }
}

bool
CheckPolyhedronQuality(const MeshContinuum& mesh,
                       const chi_mesh::Cell& cell,
                       const bool check_convexity)
{
  const auto& C = cell.centroid_;

  for (auto& face : cell.faces_)
  {
    const auto& v0 = face.centroid_;
    const size_t num_edges = face.vertex_ids_.size();

    for (size_t e = 0; e < num_edges; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids_, e);

      const auto& v1 = mesh.vertices[edge.first];
      const auto& v2 = mesh.vertices[edge.second];

      auto v01 = v1 - v0;
      auto v02 = v2 - v0;

      auto n = v01.Cross(v02);
      auto C_tri = (v0 + v1 + v2) / 3.0;

      auto CC = C_tri - C;

      if (CC.Dot(n) < 0.0) return false;
    }
  } // for face

  // Optional convexity check
  if (check_convexity)
  {
    std::vector<std::vector<uint64_t>> proxy_faces;
    for (const auto& face : cell.faces_)
      proxy_faces.push_back(face.vertex_ids_);

    size_t f = 0;
    for (const auto& face : cell.faces_)
    {
      std::set<size_t> neighbor_face_indices = FindNeighborFaceIndices(proxy_faces, f);

      for (size_t ofi : neighbor_face_indices)
      {
        auto& other_face = cell.faces_[ofi];
        auto CFC_OFC = other_face.centroid_ - face.centroid_;

        if (CFC_OFC.Dot(other_face.normal_) < 0.0) return false;
      }
      ++f;
    } // for f
  }

  return true;
}

std::set<size_t>
FindNeighborFaceIndices(const std::vector<std::vector<uint64_t>>& proxy_faces,
                        const size_t face_index)
{
  const auto& face = proxy_faces.at(face_index);
  const size_t num_faces = proxy_faces.size();

  std::set<size_t> neighbor_face_indices;
  const size_t face_num_edges = face.size();
  for (size_t e = 0; e < face_num_edges; ++e)
  {
    auto face_edge = MakeEdgeFromPolygonEdgeIndex(face, e);
    auto face_uniq_edge = MakeUniqueEdge(face_edge);

    for (size_t of = 0; of < num_faces; ++of) // other face
    {
      if (of == face_index) continue;
      const auto& other_face = proxy_faces[of];

      const size_t other_face_num_edges = other_face.size();
      for (size_t ofe = 0; ofe < other_face_num_edges; ++ofe) // other face edge
      {
        auto other_face_edge = MakeEdgeFromPolygonEdgeIndex(other_face, ofe);
        auto other_face_uniq_edge = MakeUniqueEdge(other_face_edge);

        if (face_uniq_edge == other_face_uniq_edge)
        {
          neighbor_face_indices.insert(of);
          break; // from ofe-loop
        }
      } // for ofe
    }   // for of
  }     // for face edge

  return neighbor_face_indices;
}

std::vector<Edge>
FindNonManifoldEdges(const std::vector<std::vector<uint64_t>>& proxy_faces)
{
  std::vector<Edge> non_manifold_edges;
  // The reserve amount below is chosen as follows:
  // If we cut a tet the number of non-manifold edges will be max 3.
  // A Hex will be 4. Arbitrary polyhedra as produced by STAR-CCM+ are
  // on average hexagons so lets fluff this by an additional 4 edges
  // and call 10 edges a good estimate.
  non_manifold_edges.reserve(10);

  // Build edges for each proxy face
  // This saves us some readability later on
  const size_t num_proxy_faces = proxy_faces.size();
  std::vector<std::vector<Edge>> faces_edges;
  faces_edges.reserve(num_proxy_faces);
  for (const auto& proxy_face : proxy_faces)
  {
    std::vector<Edge> face_edges(proxy_face.size());
    for (size_t e = 0; e < proxy_face.size(); ++e)
      face_edges[e] = MakeEdgeFromPolygonEdgeIndex(proxy_face, e);

    faces_edges.push_back(std::move(face_edges));
  } // for proxy face

  // Search each proxy face edge for being attached to an edge of any other proxy face
  for (size_t cf = 0; cf < num_proxy_faces; ++cf)
  {
    for (const auto& cur_edge : faces_edges[cf])
    {
      bool is_manifold = false;
      for (size_t af = 0; af < num_proxy_faces; ++af)
      {
        if (af == cf) continue;
        for (const auto& alt_edge : faces_edges[af])
        {
          if (MakeUniqueEdge(cur_edge) == MakeUniqueEdge(alt_edge))
          {
            is_manifold = true;
            goto before_next_edge;
          }
        } // for alternate edge
      }   // for alternate face

    before_next_edge:
      if (not is_manifold) { non_manifold_edges.push_back(cur_edge); }
    } // for current edge
  }   // for current face

  return non_manifold_edges;
}

std::vector<Edge>
StitchEdgesEndToEnd(const std::vector<Edge>& edges)
{
  const std::string fname = __FUNCTION__;
  std::vector<Edge> stitched_nm_edges;
  stitched_nm_edges.reserve(edges.size());

  // Declaring two queues used to stitch
  std::queue<Edge> unused_candidate_edges;
  std::queue<Edge> unused_noncandidate_edges;

  // All the edges are initially candidates
  for (auto& edge : edges)
    unused_candidate_edges.push(edge);

  // Begin the stitch using the first candidate
  // edge. This edge is no longer part of the game
  stitched_nm_edges.push_back(unused_candidate_edges.front());
  unused_candidate_edges.pop();

  // Attempt to add all candidate edges to the stitch
  const size_t num_candidates = unused_candidate_edges.size();
  for (size_t k = 0; k < num_candidates; ++k)
  {
    const auto& edge_previously_added = stitched_nm_edges.back();

    // The following while loop is gauranteed to terminate
    // because we pop the first element each time.
    while (not unused_candidate_edges.empty())
    {
      const auto& edge = unused_candidate_edges.front();

      if (edge.first == edge_previously_added.second) stitched_nm_edges.push_back(edge);
      else
        unused_noncandidate_edges.push(edge);

      unused_candidate_edges.pop(); // removes first element
    }

    std::swap(unused_noncandidate_edges, unused_candidate_edges);
    if (unused_candidate_edges.empty()) break;
  } // for k

  if (stitched_nm_edges.size() != edges.size())
    throw std::logic_error(fname + ": stitching failed.");

  return stitched_nm_edges;
}

void
PopulatePolyhedronFromFaces(const MeshContinuum& mesh,
                            const std::vector<std::vector<uint64_t>>& raw_faces,
                            chi_mesh::Cell& cell)
{
  std::vector<uint64_t> cell_vertex_ids;

  /**Lamda to add a vertex to a list.*/
  auto AddIDToCellVertexIDs = [&cell_vertex_ids](uint64_t new_id)
  {
    auto find_result = std::find(cell_vertex_ids.begin(), cell_vertex_ids.end(), new_id);

    if (find_result == cell_vertex_ids.end()) cell_vertex_ids.push_back(new_id);
  };

  // Build faces
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
  } // for raw face

  // Compute cell centroid
  chi_mesh::Vector3 cell_centroid;
  for (uint64_t vid : cell_vertex_ids)
    cell_centroid += mesh.vertices[vid];
  cell_centroid /= static_cast<double>(cell_vertex_ids.size());

  cell.centroid_ = cell_centroid;
  cell.vertex_ids_ = cell_vertex_ids;
}

void
Cut3DCell(const std::vector<ECI>& global_cut_edges,
          const std::set<uint64_t>& global_cut_vertices,
          const Vector3& plane_point,
          const Vector3& plane_normal,
          double float_compare,
          MeshContinuum& mesh,
          chi_mesh::Cell& cell,
          bool verbose)
{
  const auto& p = plane_point;
  const auto& n = plane_normal;
  const size_t cell_num_faces = cell.faces_.size();

  /**Utility lambda to check if a vertex is in generic set.*/
  auto NumberInSet = [](uint64_t number, const std::set<uint64_t>& the_set)
  {
    if (the_set.find(number) != the_set.end()) return true;

    return false;
  };

  /**Utility function to check if an edge is in the "cut_edges" list.*/
  auto EdgeIsCut = [](const Edge& edge, const std::vector<ECI>& cut_edges)
  {
    Edge edge_set(std::min(edge.first, edge.second), std::max(edge.first, edge.second));

    constexpr auto Arg1 = std::placeholders::_1;
    constexpr auto Comparator = &ECI::Comparator;

    auto result =
      std::find_if(cut_edges.begin(), cut_edges.end(), std::bind(Comparator, Arg1, edge_set));

    if (result != cut_edges.end()) return std::make_pair(true, *result);

    return std::make_pair(false, *result);
  };

  /**Utility lambda to check if a vertex is in generic list.*/
  auto NumberInList = [](uint64_t number, const std::vector<uint64_t>& list)
  {
    auto result = std::find(list.begin(), list.end(), number);

    if (result != list.end()) return true;

    return false;
  };
  auto NumberInList_t = [](size_t number, const std::vector<size_t>& list)
  {
    auto result = std::find(list.begin(), list.end(), number);

    if (result != list.end()) return true;

    return false;
  };

  /**Utility lambda to flip and edge.*/
  auto FlipEdge = [](const Edge& in_edge) { return Edge(in_edge.second, in_edge.first); };

  /**Utility lambda to convert a vector of edges to a vertex list.*/
  auto PopulateFragment =
    [NumberInList](std::vector<uint64_t>& frag, const std::vector<Edge>& frag_edges)
  {
    for (auto edge : frag_edges)
    {
      if (not NumberInList(edge.first, frag)) frag.push_back(edge.first);
      if (not NumberInList(edge.second, frag)) frag.push_back(edge.second);
    }
  };

  if (verbose) Chi::log.Log() << "Checkpoint 1";

  // Determine cut-edges relevant to this cell
  std::vector<ECI> local_cut_edges;
  const std::set<uint64_t> local_vertex_id_set(cell.vertex_ids_.begin(), cell.vertex_ids_.end());

  for (auto& eci : global_cut_edges)
    if (NumberInSet(eci.vertex_ids.first, local_vertex_id_set) and
        NumberInSet(eci.vertex_ids.second, local_vertex_id_set))
      local_cut_edges.push_back(eci);

  if (verbose)
  {
    std::stringstream outstr;
    outstr << "Local cut edges:\n";
    for (const auto& eci : local_cut_edges)
      outstr << eci.vertex_ids.first << "->" << eci.vertex_ids.second << " "
             << "cut at vertex " << eci.cut_point_id << " "
             << mesh.vertices[eci.cut_point_id].PrintStr() << "\n";
    Chi::log.Log() << outstr.str();
  }

  // Determine cut- and uncut-faces
  // A face is cut if its vertices have a
  // differing sense wrt to the plane.
  std::vector<size_t> cut_faces_indices;
  std::vector<size_t> uncut_faces_indices_A;
  std::vector<size_t> uncut_faces_indices_B;
  cut_faces_indices.reserve(cell_num_faces);
  {
    size_t face_index = 0;
    for (const auto& face : cell.faces_)
    {
      size_t num_neg_senses = 0;
      size_t num_pos_senses = 0;
      for (auto vid : face.vertex_ids_)
      {
        const auto& x = mesh.vertices[vid];
        double new_sense = n.Dot(x - p);

        if (new_sense < (0.0 - float_compare)) ++num_neg_senses;
        if (new_sense > (0.0 + float_compare)) ++num_pos_senses;

        if (num_neg_senses > 0 && num_pos_senses > 0)
        {
          cut_faces_indices.push_back(face_index);
          break;
        }
      } // for vid

      if (num_neg_senses > 0 and num_pos_senses == 0) uncut_faces_indices_A.push_back(face_index);
      if (num_neg_senses == 0 and num_pos_senses > 0) uncut_faces_indices_B.push_back(face_index);
      ++face_index;
    } // for cell face
  }

  if (verbose)
  {
    Chi::log.Log() << "Checkpoint 2";
    std::stringstream outstr;
    outstr << "Cut faces:\n";
    for (const auto& f : cut_faces_indices)
      outstr << f << " ";
    outstr << "\n";

    outstr << "Uncut faces A:\n";
    for (const auto& f : uncut_faces_indices_A)
      outstr << f << " ";
    outstr << "\n";
    outstr << "Uncut faces B:\n";
    for (const auto& f : uncut_faces_indices_B)
      outstr << f << " ";
    outstr << "\n";

    Chi::log.Log() << outstr.str();
  }

  // Split cut-faces into fragments
  // We create a map here for each cut-face.
  // We map the cut-face-index to a fragment sense-wise.
  std::map<size_t, std::vector<uint64_t>> cut_faces_fragments_A;
  std::map<size_t, std::vector<uint64_t>> cut_faces_fragments_B;

  for (const auto& cut_face_id : cut_faces_indices)
  {
    const auto& cut_face = cell.faces_[cut_face_id];

    if (verbose) Chi::log.Log() << "cut face " << cut_face_id;

    // First extract the edges that form the fragments
    std::vector<Edge> fragment_edges_A;
    std::vector<Edge> fragment_edges_B;
    const size_t face_num_verts = cut_face.vertex_ids_.size();
    for (size_t e = 0; e < face_num_verts; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(cut_face.vertex_ids_, e);

      auto edge_cut_state = EdgeIsCut(edge, local_cut_edges);
      bool edge_is_cut = edge_cut_state.first;
      ECI& edge_cut_info = edge_cut_state.second;

      // Extract edges according to if they are cut
      std::vector<Edge> extracted_edges;
      extracted_edges.reserve(2);
      if (edge_is_cut)
      {
        extracted_edges.push_back(std::make_pair(edge.first, edge_cut_info.cut_point_id));
        extracted_edges.push_back(std::make_pair(edge_cut_info.cut_point_id, edge.second));
      }
      else
        extracted_edges.push_back(edge);

      if (verbose)
      {
        std::stringstream outstr;
        for (const auto& eedge : extracted_edges)
          outstr << eedge.first << "->" << eedge.second << " ";
        Chi::log.Log() << "edge " << e << " " << edge.first << "->" << edge.second << " cut? "
                       << ((edge_is_cut) ? "yes" : "no") << " " << outstr.str();
      }

      // Enlist the edges to the corrected fragment
      for (const auto& extracted_edge : extracted_edges)
      {
        if (n.Dot(GetEdgeCentroid(extracted_edge, mesh) - p) < 0.0 - float_compare)
          fragment_edges_A.push_back(extracted_edge);
        else
          fragment_edges_B.push_back(extracted_edge);
      } // for extracted edge
    }   // for edge e of face

    // Now convert the fragment edges to vertex-id list
    std::vector<uint64_t> fragment_A;
    std::vector<uint64_t> fragment_B;

    fragment_A.reserve(fragment_edges_A.size() + 1);
    fragment_B.reserve(fragment_edges_B.size() + 1);

    PopulateFragment(fragment_A, fragment_edges_A);
    PopulateFragment(fragment_B, fragment_edges_B);

    // Finally map the fragments
    cut_faces_fragments_A[cut_face_id] = fragment_A;
    cut_faces_fragments_B[cut_face_id] = fragment_B;
  } // for cut-face

  if (verbose)
  {
    Chi::log.Log() << "Checkpoint 3";
    std::stringstream outstr;
    outstr << "Cut faces fragments A:\n";
    for (const auto& f_fragment : cut_faces_fragments_A)
    {
      outstr << "  face " << f_fragment.first << ": ";
      for (uint64_t vid : f_fragment.second)
        outstr << vid << " ";
      outstr << "\n";
    }

    outstr << "Cut faces fragments B:\n";
    for (const auto& f_fragment : cut_faces_fragments_B)
    {
      outstr << "  face " << f_fragment.first << ": ";
      for (uint64_t vid : f_fragment.second)
        outstr << vid << " ";
      outstr << "\n";
    }

    Chi::log.Log() << outstr.str();
  }

  // Make proxy-faces from data
  /**This lambda can be applied to uncut faces and cut-face fragments
   * to give a collection of proxy faces.*/
  auto MakeProxyFaces = [NumberInList_t, cut_faces_indices, FlipEdge, PopulateFragment](
                          const chi_mesh::Cell& parent_cell,
                          const std::vector<size_t>& uncut_faces_indices,
                          const std::map<size_t, std::vector<uint64_t>>& cut_face_fragments)
  {
    // Build proxy faces based on uncut-faces and cut faces
    std::vector<std::vector<uint64_t>> proxy_faces;
    proxy_faces.reserve(parent_cell.faces_.size());

    size_t f = 0;
    for (const auto& face : parent_cell.faces_)
    {
      if (NumberInList_t(f, uncut_faces_indices)) proxy_faces.push_back(face.vertex_ids_);
      else if (NumberInList_t(f, cut_faces_indices))
        proxy_faces.push_back(cut_face_fragments.at(f));
      ++f;
    } // for face f

    // Now build the interface proxy-face
    // Find non-manifold edges
    auto non_manifold_edges = FindNonManifoldEdges(proxy_faces);

    // Flip the non-manifold edges
    for (auto& edge : non_manifold_edges)
      edge = FlipEdge(edge);

    // Stitch the edges end-to-end
    auto stitched_nm_edges = StitchEdgesEndToEnd(non_manifold_edges);

    std::vector<uint64_t> interface_proxy_face;
    PopulateFragment(interface_proxy_face, stitched_nm_edges);

    proxy_faces.push_back(std::move(interface_proxy_face));

    return proxy_faces;
  };

  auto proxies_A = MakeProxyFaces(cell, uncut_faces_indices_A, cut_faces_fragments_A);
  auto proxies_B = MakeProxyFaces(cell, uncut_faces_indices_B, cut_faces_fragments_B);

  chi_mesh::Cell cell_A(CellType::POLYHEDRON, CellType::POLYHEDRON);
  auto cell_A_ptr = &cell_A;
  auto cell_B_ptr = std::make_unique<chi_mesh::Cell>(CellType::POLYHEDRON, CellType::POLYHEDRON);

  PopulatePolyhedronFromFaces(mesh, proxies_A, *cell_A_ptr);
  PopulatePolyhedronFromFaces(mesh, proxies_B, *cell_B_ptr);

  cell_A_ptr->local_id_ = cell.local_id_;
  cell_A_ptr->global_id_ = cell.global_id_;

  cell_B_ptr->global_id_ = mesh.local_cells.size();

  cell_A_ptr->material_id_ = cell.material_id_;
  cell_B_ptr->material_id_ = cell.material_id_;

  if (verbose)
  {
    std::set<uint64_t> verts;
    for (uint64_t vid : cell_A_ptr->vertex_ids_)
      verts.insert(vid);
    for (uint64_t vid : cell_B_ptr->vertex_ids_)
      verts.insert(vid);
    for (uint64_t vid : cell.vertex_ids_)
      verts.insert(vid);
    for (uint64_t vid : verts)
      Chi::log.Log() << vid << " " << mesh.vertices[vid].PrintStr();
    Chi::log.Log() << "Checkpoint 4:\n"
                   << cell.ToString() << cell_A_ptr->ToString() << cell_B_ptr->ToString();
  }

  cell = cell_A;
  mesh.cells.push_back(std::move(cell_B_ptr));
}

void
CutMeshWithPlane(MeshContinuum& mesh,
                 const Vector3& plane_point,
                 const Vector3& plane_normal,
                 double merge_tolerance,
                 double float_compare)
{
  const std::string function_name = __FUNCTION__;

  const auto& p = plane_point;
  const auto& n = plane_normal.Normalized();

  Chi::log.Log() << "Cutting mesh with plane. "
                 << "Ref. Point: " << p.PrintS() << " Normal: " << n.PrintS();

  // Snap vertices to plane within merge tolerance to avoid creating small cells or cutting parallel
  // faces Order N, num_vertices
  size_t num_verts_snapped = 0;
  for (auto& id_vertex : mesh.vertices)
  {
    auto& vertex = id_vertex.second;
    double d_from_plane = n.Dot(vertex - p);

    if (std::fabs(d_from_plane) < merge_tolerance)
    {
      vertex -= n * d_from_plane;
      ++num_verts_snapped;
    }
  }

  Chi::log.Log() << "Number of vertices snapped to plane: " << num_verts_snapped;

  // Perform quality checks
  size_t num_bad_quality_cells = 0;
  for (const auto& cell : mesh.local_cells)
  {
    if (cell.Type() == CellType::POLYGON)
    {
      if (not CheckPolygonQuality(mesh, cell)) ++num_bad_quality_cells;
    }
    else if (cell.Type() == CellType::POLYHEDRON)
    {
      if (not CheckPolyhedronQuality(mesh, cell, true)) ++num_bad_quality_cells;
    }
    else
      throw std::logic_error(function_name + ": Called for a mesh containing"
                                             " an unsupported cell-type.");
  }
  if (num_bad_quality_cells > 0)
    throw std::logic_error(function_name + ": Called for a mesh containing " +
                           std::to_string(num_bad_quality_cells) + " bad quality cells.");

  // Determine cells to cut
  // Order N, num_cells
  // A cell is a candidate for cutting if its
  // vertices lay on both sides of the plane.
  // So this algorithm just checks the sense wrt
  // the plane. Works for both 2D and 3D
  std::vector<chi_mesh::Cell*> cells_to_cut;

  for (auto& cell : mesh.local_cells)
  {
    size_t num_neg_senses = 0;
    size_t num_pos_senses = 0;
    for (auto vid : cell.vertex_ids_)
    {
      const auto& x = mesh.vertices[vid];
      double new_sense = n.Dot(x - p);

      if (new_sense < (0.0 - float_compare)) ++num_neg_senses;
      if (new_sense > (0.0 + float_compare)) ++num_pos_senses;

      if (num_neg_senses > 0 && num_pos_senses > 0)
      {
        cells_to_cut.emplace_back(&cell);
        break;
      }
    } // for vid
  }   // for cell
  Chi::log.Log() << "Number of cells to cut: " << cells_to_cut.size();

  // Two-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYGON)
  {
    // Determine cut vertices
    std::set<uint64_t> cut_vertices;
    {
      for (auto& cell_ptr : cells_to_cut)
        for (uint64_t vid : cell_ptr->vertex_ids_)
        {
          const auto vertex = mesh.vertices[vid];
          double dv = std::fabs((vertex - p).Dot(n));
          if (dv < float_compare) cut_vertices.insert(vid);
        } // for vid
    }     // populate cut_vertices

    Chi::log.Log() << "Number of cut vertices: " << cut_vertices.size();

    // Build unique edges
    size_t num_edges_cut = 0;
    std::set<Edge> edges_set;
    for (auto& cell_ptr : cells_to_cut)
    {
      const auto& cell = *cell_ptr;
      const size_t num_edges = cell.vertex_ids_.size();

      for (size_t e = 0; e < num_edges; ++e)
      {
        auto edge = MakeEdgeFromPolygonEdgeIndex(cell.vertex_ids_, e);
        edges_set.insert(
          std::make_pair(std::min(edge.first, edge.second), std::max(edge.first, edge.second)));
      }
    } // for cell - built edges_set

    uint64_t new_vertex_address = mesh.GetGlobalVertexCount();

    // Determine cut edges
    std::vector<ECI> cut_edges;
    {
      for (auto& edge : edges_set)
      {
        const auto& v0 = mesh.vertices[edge.first];
        const auto& v1 = mesh.vertices[edge.second];

        chi_mesh::Vector3 cut_point;

        if (CheckPlaneLineIntersect(n, p, v0, v1, cut_point))
        {
          double dv0 = std::fabs((v0 - p).Dot(n));
          double dv1 = std::fabs((v1 - p).Dot(n));
          if (dv0 > float_compare and dv1 > float_compare)
          {
            mesh.vertices.Insert(new_vertex_address, cut_point);
            cut_edges.emplace_back(edge, new_vertex_address++);
            ++num_edges_cut;
          }
        }
      } // for edge - determine cut
    }   // populate edges cut

    Chi::log.Log() << "Number of cut edges: " << num_edges_cut;

    // Process cells that are cut
    for (auto& cell_ptr : cells_to_cut)
    {
      auto& cell = *cell_ptr;

      CutPolygon(cut_edges, cut_vertices, p, n, mesh, cell);
    } // for cell_ptr
  }   // two-D

  // Three-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYHEDRON)
  {
    // Determine cut vertices
    std::set<uint64_t> cut_vertices;
    {
      for (auto& cell_ptr : cells_to_cut)
        for (uint64_t vid : cell_ptr->vertex_ids_)
        {
          const auto& vertex = mesh.vertices[vid];
          double dv = std::fabs((vertex - p).Dot(n));
          if (dv < float_compare) cut_vertices.insert(vid);
        } // for vid
    }     // populate cut_vertices

    // Build unique edges
    size_t num_edges_cut = 0;
    std::set<Edge> edges_set;
    for (auto& cell_ptr : cells_to_cut)
    {
      const auto& cell = *cell_ptr;

      for (auto& face : cell.faces_)
      {
        const size_t num_edges = face.vertex_ids_.size();

        for (size_t e = 0; e < num_edges; ++e)
        {
          auto edge = MakeEdgeFromPolygonEdgeIndex(face.vertex_ids_, e);
          edges_set.insert(
            std::make_pair(std::min(edge.first, edge.second), std::max(edge.first, edge.second)));
        }
      } // for face
    }   // for cell - built edges_set

    uint64_t new_vertex_address = mesh.GetGlobalVertexCount();

    // Determine cut edges
    std::vector<ECI> cut_edges;
    {
      for (auto& edge : edges_set)
      {
        const auto& v0 = mesh.vertices[edge.first];
        const auto& v1 = mesh.vertices[edge.second];

        chi_mesh::Vector3 cut_point;

        if (CheckPlaneLineIntersect(n, p, v0, v1, cut_point))
        {
          double dv0 = std::fabs((v0 - p).Dot(n));
          double dv1 = std::fabs((v1 - p).Dot(n));
          if (dv0 > float_compare and dv1 > float_compare)
          {
            mesh.vertices.Insert(new_vertex_address, cut_point);
            cut_edges.emplace_back(edge, new_vertex_address++);
            ++num_edges_cut;
          }
        }
      } // for edge - determine cut
    }   // populate edges cut

    Chi::log.Log() << "Number of cut edges: " << num_edges_cut;

    // Process cells that are cut
    for (auto& cell_ptr : cells_to_cut)
    {
      Cut3DCell(cut_edges, cut_vertices, plane_point, plane_normal, float_compare, mesh, *cell_ptr);
    } // for cell_ptr
  }

  Chi::log.Log() << "Done cutting mesh with plane. Num cells = " << mesh.local_cells.size();
}

} // namespace chi_mesh::mesh_cutting
