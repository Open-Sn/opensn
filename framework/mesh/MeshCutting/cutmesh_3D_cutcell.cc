#include "meshcutting.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <queue>
#include <algorithm>
#include <functional>

//###################################################################
/***/
void chi_mesh::mesh_cutting::
  Cut3DCell(const std::vector<ECI> &global_cut_edges,
            const std::set<uint64_t> &global_cut_vertices,
            const Vector3 &plane_point,
            const Vector3 &plane_normal,
            double float_compare,
            MeshContinuum &mesh,
            chi_mesh::Cell &cell,
            bool verbose/*=false*/)
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
    Edge edge_set(std::min(edge.first,edge.second),
                  std::max(edge.first,edge.second));

    constexpr auto Arg1       = std::placeholders::_1;
    constexpr auto Comparator = &ECI::Comparator;

    auto result = std::find_if(cut_edges.begin(),
                               cut_edges.end(),
                               std::bind(Comparator,Arg1,edge_set));

    if (result != cut_edges.end())
      return std::make_pair(true,*result);

    return std::make_pair(false,*result);
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
  auto FlipEdge = [](const Edge& in_edge)
  {
    return Edge(in_edge.second, in_edge.first);
  };

  /**Utility lambda to convert a vector of edges to a vertex list.*/
  auto PopulateFragment = [NumberInList](std::vector<uint64_t>& frag,
                                           const std::vector<Edge>& frag_edges)
  {
    for (auto edge : frag_edges)
    {
      if (not NumberInList(edge.first , frag)) frag.push_back(edge.first );
      if (not NumberInList(edge.second, frag)) frag.push_back(edge.second);
    }
  };

  if (verbose) Chi::log.Log() << "Checkpoint 1";

  //============================================= Determine cut-edges relevant
  //                                              to this cell
  std::vector<ECI> local_cut_edges;
  const std::set<uint64_t> local_vertex_id_set(cell.vertex_ids_.begin(),
                                               cell.vertex_ids_.end());

  for (auto& eci : global_cut_edges)
    if (NumberInSet(eci.vertex_ids.first , local_vertex_id_set) and
        NumberInSet(eci.vertex_ids.second, local_vertex_id_set))
      local_cut_edges.push_back(eci);

  if (verbose)
  {
    std::stringstream outstr;
    outstr << "Local cut edges:\n";
    for (const auto& eci : local_cut_edges)
      outstr << eci.vertex_ids.first << "->"
             << eci.vertex_ids.second << " "
             << "cut at vertex " << eci.cut_point_id << " "
             << mesh.vertices[eci.cut_point_id].PrintStr()
             << "\n";
    Chi::log.Log() << outstr.str();
  }

  //============================================= Determine cut- and uncut-faces
  // A face is cut if its vertices have a
  // differing sense wrt to the plane.
  std::vector<size_t> cut_faces_indices;
  std::vector<size_t> uncut_faces_indices_A;
  std::vector<size_t> uncut_faces_indices_B;
  cut_faces_indices.reserve(cell_num_faces);
  {
    size_t face_index=0;
    for (const auto& face : cell.faces_)
    {
      size_t num_neg_senses = 0;
      size_t num_pos_senses = 0;
      for (auto vid : face.vertex_ids_)
      {
        const auto& x = mesh.vertices[vid];
        double new_sense = n.Dot(x-p);

        if (new_sense < (0.0-float_compare)) ++num_neg_senses;
        if (new_sense > (0.0+float_compare)) ++num_pos_senses;

        if (num_neg_senses>0 && num_pos_senses>0)
        { cut_faces_indices.push_back(face_index); break; }
      }//for vid

      if (num_neg_senses >  0 and num_pos_senses == 0)
        uncut_faces_indices_A.push_back(face_index);
      if (num_neg_senses == 0 and num_pos_senses >  0)
        uncut_faces_indices_B.push_back(face_index);
      ++face_index;
    }//for cell face
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

  //============================================= Split cut-faces into fragments
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
    for (size_t e=0; e<face_num_verts; ++e)
    {
      auto edge = MakeEdgeFromPolygonEdgeIndex(cut_face.vertex_ids_, e);

      auto edge_cut_state = EdgeIsCut(edge, local_cut_edges);
      bool edge_is_cut   = edge_cut_state.first;
      ECI& edge_cut_info = edge_cut_state.second;

      // Extract edges according to if they are cut
      std::vector<Edge> extracted_edges;
      extracted_edges.reserve(2);
      if (edge_is_cut)
      {
        extracted_edges.push_back(
          std::make_pair(edge.first, edge_cut_info.cut_point_id));
        extracted_edges.push_back(
          std::make_pair(edge_cut_info.cut_point_id, edge.second));
      }
      else
        extracted_edges.push_back(edge);

      if (verbose)
      {
        std::stringstream outstr;
        for (const auto& eedge : extracted_edges)
          outstr << eedge.first << "->" << eedge.second << " ";
        Chi::log.Log() << "edge " << e << " "
                      << edge.first << "->" << edge.second
                      << " cut? " << ((edge_is_cut)? "yes" : "no")
                      << " " << outstr.str();
      }

      // Enlist the edges to the corrected fragment
      for (const auto& extracted_edge : extracted_edges)
      {
        if (n.Dot(GetEdgeCentroid(extracted_edge, mesh) - p) < 0.0 - float_compare)
          fragment_edges_A.push_back(extracted_edge);
        else
          fragment_edges_B.push_back(extracted_edge);
      }//for extracted edge
    }// for edge e of face

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
  }//for cut-face

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

  //============================================= Make proxy-faces from data
  /**This lambda can be applied to uncut faces and cut-face fragments
   * to give a collection of proxy faces.*/
  auto MakeProxyFaces = [NumberInList_t, cut_faces_indices, FlipEdge,
                       PopulateFragment]
    (const chi_mesh::Cell& parent_cell,
     const std::vector<size_t>& uncut_faces_indices,
     const std::map<size_t, std::vector<uint64_t>>& cut_face_fragments)
  {
    //================================= Build proxy faces based on uncut-faces
    //                                  and cut faces
    std::vector<std::vector<uint64_t>> proxy_faces;
    proxy_faces.reserve(parent_cell.faces_.size());

    size_t f=0;
    for (const auto& face : parent_cell.faces_)
    {
      if (NumberInList_t(f, uncut_faces_indices))
        proxy_faces.push_back(face.vertex_ids_);
      else if (NumberInList_t(f, cut_faces_indices))
        proxy_faces.push_back(cut_face_fragments.at(f));
      ++f;
    }//for face f

    //================================= Now build the interface proxy-face
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

  auto proxies_A =
    MakeProxyFaces(cell, uncut_faces_indices_A, cut_faces_fragments_A);
  auto proxies_B =
    MakeProxyFaces(cell, uncut_faces_indices_B, cut_faces_fragments_B);

  chi_mesh::Cell cell_A(CellType::POLYHEDRON, CellType::POLYHEDRON);
  auto cell_A_ptr = &cell_A;
  auto cell_B_ptr = std::make_unique<chi_mesh::Cell>(CellType::POLYHEDRON,
                                                     CellType::POLYHEDRON);

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
                  << cell.ToString()
                  << cell_A_ptr->ToString()
                  << cell_B_ptr->ToString();
  }

  cell = cell_A;
  mesh.cells.push_back(std::move(cell_B_ptr));


}