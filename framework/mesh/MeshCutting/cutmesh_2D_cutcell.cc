#include "meshcutting.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include <algorithm>
#include <functional>

//###################################################################
/**Cuts a polygon.*/
void chi_mesh::mesh_cutting::
  CutPolygon(const std::vector<ECI> &cut_edges,
             const std::set<uint64_t> &cut_vertices,
             const Vector3 &plane_point,
             const Vector3 &plane_normal,
             MeshContinuum &mesh,
             chi_mesh::Cell& cell)
{
  const auto& p = plane_point;
  const auto& n = plane_normal;

  /**Utility lambda to check if a vertex is in "cut_vertices" list.*/
  auto VertexIsCut = [&cut_vertices](uint64_t vid)
  {
    auto result = cut_vertices.find(vid);

    if (result != cut_vertices.end())
      return true;

    return false;
  };

  /**Utility function to check if an edge is in the "cut_edges" list.*/
  auto EdgeIsCut = [&cut_edges](const Edge& edge)
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

  //============================================= Create and set vertex and edge
  //                                              cut flags for the current cell.
  //                                              Also populate edge_cut_info
  size_t num_verts = cell.vertex_ids_.size();
  size_t num_edges = num_verts;

  std::vector<bool> vertex_cut_flags(num_verts,false);
  std::vector<bool> edge_cut_flags(num_edges,false);
  std::vector<ECI>  edge_cut_info(num_edges);

  for (size_t e=0; e<num_edges; ++e)
  {
    vertex_cut_flags[e] = VertexIsCut(cell.vertex_ids_[e]);

    auto edge = MakeEdgeFromPolygonEdgeIndex(cell.vertex_ids_, e);
    auto cut_nature = EdgeIsCut(edge);
    edge_cut_flags[e] = cut_nature.first;
    if (cut_nature.first) edge_cut_info[e] = cut_nature.second;

    edge_cut_info[e].vertex_ids = edge;
  }//populate flags

  //============================================= Lamda for edge loop
  enum class CurVertex
  {
    AT_FIRST,
    AT_CUT_POINT,
    AT_SECOND,
    NONE
  };

  struct CurCutInfo
  {
    int which_edge=0;
    CurVertex which_vertex=CurVertex::AT_FIRST;

    CurCutInfo() = default;

    CurCutInfo(int in_which_edge, CurVertex in_which_vertex) :
      which_edge(in_which_edge),
      which_vertex(in_which_vertex) {}
  };

  /**This lamda function starts from a current cut-edge, which is either
   * an edge where the first vertex is cut or an edge that is cut
   * somewhere along its length, and then follows the edges in a ccw fashion
   * until it finds another cut. This last cut is just as either an edge
   * cut along its length or cut at the second vertex. This then completes
   * an edge loop that can be used to define another polygon.*/
  auto GetVerticesTillNextCut =
    [&cell,&edge_cut_flags,&edge_cut_info,&VertexIsCut](
      CurCutInfo start_cut_info)
    {
      size_t num_verts = cell.vertex_ids_.size();
      std::vector<uint64_t> vertex_ids;
      vertex_ids.reserve(num_verts); //Ought to be more than enough

      int  e        = start_cut_info.which_edge;
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
      }//switch type of starting cut

      //Look at downstream ccw edges and check for
      //edges cut or end-point cuts
      for (int eref=0; eref<num_verts; ++eref)
      {
        e = (e<(num_verts-1))? e+1 : 0;

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
      }//for eref

      skip_to_return_portion:
      CurCutInfo end_cut_info(e, end_type);

      return std::make_pair(vertex_ids,end_cut_info);
    };

  typedef std::pair<std::vector<uint64_t>,CurCutInfo> LoopInfo;

  //============================================= Process all edges and create
  //                                              edge loops with associated
  //                                              cells
  std::vector<std::vector<uint64_t>> loops_to_add_to_mesh;
  std::vector<CurCutInfo> cut_history;
  for (size_t e=0; e<num_edges; ++e)
  {
    LoopInfo loop_info;

    if (vertex_cut_flags[e])
      loop_info = GetVerticesTillNextCut(CurCutInfo(e,CurVertex::AT_FIRST));
    else if (edge_cut_flags[e])
      loop_info = GetVerticesTillNextCut(CurCutInfo(e,CurVertex::AT_CUT_POINT));
    else continue;

    std::vector<uint64_t> verts_to_next_cut = loop_info.first;
    int                   end_edge          = loop_info.second.which_edge;
    CurVertex             end_type          = loop_info.second.which_vertex;

    //Notes:
    // - If the end_edge == e then this means trouble as a polygon cannot
    //   be cut like that.
    // - If end_edge < e then the end-edge is definitely the edge right before
    //   the first cut edge. We should still process this edge-loop, but stop
    //   searching.

    if (end_edge < e) //if looped past num_edges
      e = int(num_edges)-1; //stop search
    else if (end_type == CurVertex::AT_SECOND)
      e = end_edge;            //resume search after end_e
    else
      e = end_edge - 1;        //resume search at end_e. e will get ++ to end_edge

    loops_to_add_to_mesh.push_back(verts_to_next_cut);
  }//for e

  //================================================== Add derivative cells to
  //                                                   mesh

  // Take the back cell and paste onto
  // the current cell reference
  if (not loops_to_add_to_mesh.empty())
  {
    auto& back_loop = loops_to_add_to_mesh.back();
    PopulatePolygonFromVertices(mesh, back_loop, cell);

    loops_to_add_to_mesh.pop_back();
  }//if there are cells to add

  // Now push-up new cells from the remainder
  for (auto& loop : loops_to_add_to_mesh)
  {
    auto new_cell = std::make_unique<chi_mesh::Cell>(CellType::POLYGON,
                                                     CellType::TRIANGLE);
    PopulatePolygonFromVertices(mesh, loop, *new_cell);
    mesh.cells.push_back(std::move(new_cell));
  }
}

