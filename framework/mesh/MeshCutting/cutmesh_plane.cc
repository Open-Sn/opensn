#include "meshcutting.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/Raytrace/raytracing.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <algorithm>

//###################################################################
/**Cuts a mesh with a plane.*/
void
chi_mesh::mesh_cutting::CutMeshWithPlane(MeshContinuum& mesh,
                                         const Vector3& plane_point,
                                         const Vector3& plane_normal,
                                         double merge_tolerance /*=1.0e-3*/,
                                         double float_compare /*=1.0e-10*/)
{
  const std::string function_name = __FUNCTION__;

  const auto& p = plane_point;
  const auto& n = plane_normal.Normalized();

  Chi::log.Log() << "Cutting mesh with plane. "
                 << "Ref. Point: " << p.PrintS() << " Normal: " << n.PrintS();

  //============================================= Snap vertices to plane within
  //                                              merge tolerance to avoid
  //                                              creating small cells or
  //                                              cutting parallel faces
  // Order N, num_vertices
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

  //============================================= Perform quality checks
  size_t num_bad_quality_cells = 0;
  for (const auto& cell : mesh.local_cells)
  {
    if (cell.Type() == CellType::POLYGON)
    {
      if (not CheckPolygonQuality(mesh, cell)) ++num_bad_quality_cells;
    }
    else if (cell.Type() == CellType::POLYHEDRON)
    {
      if (not CheckPolyhedronQuality(mesh, cell, /*check_convexity*/ true)) ++num_bad_quality_cells;
    }
    else
      throw std::logic_error(function_name + ": Called for a mesh containing"
                                             " an unsupported cell-type.");
  }
  if (num_bad_quality_cells > 0)
    throw std::logic_error(function_name + ": Called for a mesh containing " +
                           std::to_string(num_bad_quality_cells) + " bad quality cells.");

  //============================================= Determine cells to cut
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

  //============================================= Two-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYGON)
  {
    //====================================== Determine cut vertices
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

    //====================================== Build unique edges
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

    //====================================== Determine cut edges
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

    //====================================== Process cells that are cut
    for (auto& cell_ptr : cells_to_cut)
    {
      auto& cell = *cell_ptr;

      CutPolygon(cut_edges, cut_vertices, p, n, mesh, cell);
    } // for cell_ptr
  }   // two-D

  //============================================= Three-D algorithm
  if (mesh.local_cells[0].Type() == CellType::POLYHEDRON)
  {
    //====================================== Determine cut vertices
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

    //====================================== Build unique edges
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

    //====================================== Determine cut edges
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

    //====================================== Process cells that are cut
    for (auto& cell_ptr : cells_to_cut)
    {
      Cut3DCell(cut_edges, cut_vertices, plane_point, plane_normal, float_compare, mesh, *cell_ptr);
    } // for cell_ptr
  }

  Chi::log.Log() << "Done cutting mesh with plane. Num cells = " << mesh.local_cells.size();
}