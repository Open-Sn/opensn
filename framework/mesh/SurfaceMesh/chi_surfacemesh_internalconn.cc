#include "chi_surfacemesh.h"

#include <chi_log.h>

#include "utils/chi_timer.h"

#include <algorithm>

//#########################################################
/** Runs over the faces of the surfacemesh and determines
 * neighbors. The algorithm first establishes which cells subscribe to each
 * vertex and then loops over faces and edges. For each edge, only the
 * subscribing faces are searched for neighbors. This routine has
 * time complexity O(N).*/
void
chi_mesh::SurfaceMesh::UpdateInternalConnectivity()
{
  std::vector<std::vector<size_t>> vertex_subscriptions;

  //======================================== Initialize vertex subscription
  size_t num_verts = vertices_.size();
  vertex_subscriptions.resize(num_verts);

  for (auto& vert_sub : vertex_subscriptions)
    vert_sub.reserve(5);

  //======================================== Loop over cells and determine
  //                                         which cells subscribe to a vertex
  //%%%%%% TRIANGLES %%%%%
  size_t num_tri_faces = faces_.size();
  for (size_t tf = 0; tf < num_tri_faces; tf++)
  {
    auto& try_face = faces_[tf];
    for (int v = 0; v < 3; ++v)
      vertex_subscriptions[v].push_back(tf);
  }
  //%%%%% POLYGONS %%%%%
  size_t num_poly_faces = poly_faces_.size();
  for (size_t pf = 0; pf < num_poly_faces; pf++)
  {
    auto poly_face = poly_faces_[pf];
    for (auto v : poly_face->v_indices)
      vertex_subscriptions[v].push_back(pf);
  }

  //======================================== Loop over cells and determine
  //                                         connectivity
  //%%%%%% TRIANGLES %%%%%
  for (auto curFace : faces_)
  {
    for (int e = 0; e < 3; e++)
    {
      int* curface_edge = curFace.e_index[e];
      int vi = curface_edge[0];
      int vf = curface_edge[1];

      //=============================== Search cells subscribing to vi
      for (auto ofi : vertex_subscriptions[vi])
      {
        auto& other_face = faces_[ofi];

        for (size_t e2 = 0; e2 < 3; e2++)
        {
          if ((curface_edge[0] == other_face.e_index[e2][1]) &&
              (curface_edge[1] == other_face.e_index[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for ofi

      //=============================== Search cells subscribing to vf
      for (auto ofi : vertex_subscriptions[vf])
      {
        auto& other_face = faces_[ofi];

        for (size_t e2 = 0; e2 < 3; e2++)
        {
          if ((curface_edge[0] == other_face.e_index[e2][1]) &&
              (curface_edge[1] == other_face.e_index[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for ofi
    }     // for current face edges
  }       // for faces

  //======================================== Loop over cells and determine
  //                                         connectivity
  //%%%%% POLYGONS %%%%%
  for (auto curFace : poly_faces_)
  {
    for (auto& curface_edge : curFace->edges)
    {
      int vi = curface_edge[0];
      int vf = curface_edge[1];

      //=============================== Search cells subscribing to vi
      for (auto ofi : vertex_subscriptions[vi])
      {
        auto other_face = poly_faces_[ofi];

        for (size_t e2 = 0; e2 < other_face->edges.size(); e2++)
        {
          if ((curface_edge[0] == other_face->edges[e2][1]) &&
              (curface_edge[1] == other_face->edges[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for ofi

      //=============================== Search cells subscribing to vf
      for (auto ofi : vertex_subscriptions[vf])
      {
        auto other_face = poly_faces_[ofi];

        for (size_t e2 = 0; e2 < other_face->edges.size(); e2++)
        {
          if ((curface_edge[0] == other_face->edges[e2][1]) &&
              (curface_edge[1] == other_face->edges[e2][0]))
          {
            curface_edge[2] = ofi; // cell index
            curface_edge[3] = e2;  // edge index
          }
        } // for e2
      }   // for other faces
    }     // for current face edges
  }       // for faces
}