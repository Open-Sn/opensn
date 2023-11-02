#include "raytracing.h"

#include "mesh/Cell/cell.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Performs raytracing within a 3D Polyhedron.*/
void chi_mesh::RayTracer::TracePolyhedron(const Cell &cell,
                                          Vector3 &pos_i,
                                          Vector3 &omega_i,
                                          bool& intersection_found,
                                          bool& backward_tolerance_hit,
                                          RayTracerOutputInformation& oi)
{
  const auto& grid = Grid();
  chi_mesh::Vector3 ip = pos_i; //Intersection point

  std::vector<RayTracerOutputInformation> triangle_intersections;

  size_t num_faces = cell.faces_.size();
  triangle_intersections.reserve(num_faces*4); //Guessing 4 tris per face
  for (int f=0; f<num_faces; f++)
  {
    const auto& face = cell.faces_[f];

    size_t num_sides = face.vertex_ids_.size();
    for (size_t s=0; s<num_sides; s++)
    {
      size_t v0_index = face.vertex_ids_[s];
      size_t v1_index = (s<(num_sides-1)) ?     //if not last vertex
                        face.vertex_ids_[s + 1] :
                        face.vertex_ids_[0];    //else
      const auto& v0 = grid.vertices[v0_index];
      const auto& v1 = grid.vertices[v1_index];
      const auto& v2 = cell.faces_[f].centroid_;

      auto v01 = v1 - v0;
      auto v02 = v2 - v0;
      auto n_est = v01.Cross(v02);

      if (n_est.Dot(omega_i) < 0.0) continue;

      RayTracerOutputInformation triangle_oi;

      bool intersects =
        chi_mesh::CheckLineIntersectTriangle2(v0,v1,v2,pos_i,omega_i,ip);

      if (intersects)
      {
        triangle_oi.distance_to_surface = (ip - pos_i).Norm();
        triangle_oi.pos_f = ip;

        triangle_oi.destination_face_index = f;
        triangle_oi.destination_face_neighbor = cell.faces_[f].neighbor_id_;

        intersection_found = true;
        triangle_intersections.emplace_back(std::move(triangle_oi));
        if (not perform_concavity_checks_) break;
      }//if intersects
    }//for side

    if (intersection_found and (not perform_concavity_checks_)) break;
  }//for faces

  //======================================== Determine closest intersection
  if (not perform_concavity_checks_ and not triangle_intersections.empty())
    oi = triangle_intersections.back();
  else if (perform_concavity_checks_ and not triangle_intersections.empty())
  {
    auto closest_intersection = &triangle_intersections.back();
    for (auto& intersection : triangle_intersections)
      if (intersection.distance_to_surface <
          closest_intersection->distance_to_surface)
        closest_intersection = &intersection;

    oi = *closest_intersection;
  }
  else
  {
    RayTracerOutputInformation blank_oi;
    blank_oi.distance_to_surface = 1.0e15;
    oi = blank_oi;
  }
}