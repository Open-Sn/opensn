#include "raytracing.h"

#include "mesh/Cell/cell.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Performs raytracing within a 1D-slab.*/
void
chi_mesh::RayTracer::TraceSlab(const Cell& cell,
                               Vector3& pos_i,
                               Vector3& omega_i,
                               bool& intersection_found,
                               bool& backward_tolerance_hit,
                               RayTracerOutputInformation& oi)
{
  const auto& grid = Grid();
  chi_mesh::Vector3 intersection_point;
  std::pair<double, double> weights;

  const double fabs_mu = std::fabs(omega_i.Dot(cell.faces_[0].normal_));

  double d_extend = (fabs_mu < 1.0e-15) ? 1.0e15 : extension_distance_ / fabs_mu;

  chi_mesh::Vector3 pos_f_line = pos_i + omega_i * d_extend;

  int num_faces = 2;
  for (int f = 0; f < num_faces; f++)
  {
    uint64_t fpi = cell.vertex_ids_[f]; // face point index
    chi_mesh::Vertex face_point = grid.vertices[fpi];

    bool intersects = chi_mesh::CheckPlaneLineIntersect(
      cell.faces_[f].normal_, face_point, pos_i, pos_f_line, intersection_point, &weights);

    double D = weights.first * d_extend;

    if ((D > backward_tolerance_) and intersects)
    {
      oi.distance_to_surface = D;
      oi.pos_f = intersection_point;

      oi.destination_face_index = f;
      oi.destination_face_neighbor = cell.faces_[f].neighbor_id_;
      intersection_found = true;
      break;
    }
    if (intersects) backward_tolerance_hit = true;
  } // for faces
}