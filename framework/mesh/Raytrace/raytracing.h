#pragma once

#include "../chi_mesh.h"

namespace chi_mesh
{
/**Data structure to hold output info from the raytracer.*/
struct RayTracerOutputInformation
{
  double distance_to_surface = 0.0;
  Vector3 pos_f;
  unsigned int destination_face_index = 0;
  uint64_t destination_face_neighbor = 0;
  bool particle_lost = false;
  std::string lost_particle_info;
};

//###################################################################
/**Raytracer object.*/
class RayTracer
{
private:
  const chi_mesh::MeshContinuum& reference_grid_;
  std::vector<double> cell_sizes_;
  double epsilon_nudge_ = 1.0e-8;
  double backward_tolerance_ = 1.0e-10;
  double extension_distance_ = 1.0e5;
  bool perform_concavity_checks_ = true;

public:
  explicit RayTracer(const chi_mesh::MeshContinuum& grid,
                     std::vector<double> in_cell_sizes,
                     bool in_perform_concavity_checks = true)
    : reference_grid_(grid),
      cell_sizes_(std::move(in_cell_sizes)),
      perform_concavity_checks_(in_perform_concavity_checks)
  {
  }

private:
  const chi_mesh::MeshContinuum& Grid() const;

  void SetTolerancesFromCellSize(double cell_size)
  {
    epsilon_nudge_ = cell_size * 1.0e-2;
    backward_tolerance_ = cell_size * 1.0e-10;
    extension_distance_ = 3.0 * cell_size;
  }

public:
  /**Traces a ray with an initial position either within the cell or
   * on the cell surface, and with a direction vector pointing inward
   * toward the cell. If the ray-trace fails the particle will be marked
   * as lost.*/
  RayTracerOutputInformation
  TraceRay(const Cell& cell, Vector3& pos_i, Vector3& omega_i, int function_depth = 0);

  /**Traces a ray with an initial position, presumed to be outside the cell,
   * to an incident face.*/
  RayTracerOutputInformation
  TraceIncidentRay(const Cell& cell, const Vector3& pos_i, const Vector3& omega_i);

private:
  void TraceSlab(const Cell& cell,
                 Vector3& pos_i,
                 Vector3& omega_i,
                 bool& intersection_found,
                 bool& backward_tolerance_hit,
                 RayTracerOutputInformation& oi);
  void TracePolygon(const Cell& cell,
                    Vector3& pos_i,
                    Vector3& omega_i,
                    bool& intersection_found,
                    bool& backward_tolerance_hit,
                    RayTracerOutputInformation& oi);
  void TracePolyhedron(const Cell& cell,
                       Vector3& pos_i,
                       Vector3& omega_i,
                       bool& intersection_found,
                       bool& backward_tolerance_hit,
                       RayTracerOutputInformation& oi);
};

bool CheckPlaneLineIntersect(const chi_mesh::Normal& plane_normal,
                             const chi_mesh::Vector3& plane_point,
                             const chi_mesh::Vector3& line_point_0,
                             const chi_mesh::Vector3& line_point_1,
                             chi_mesh::Vector3& intersection_point,
                             std::pair<double, double>* weights = nullptr);

bool CheckLineIntersectStrip(const chi_mesh::Vector3& strip_point0,
                             const chi_mesh::Vector3& strip_point1,
                             const chi_mesh::Vector3& strip_normal,
                             const chi_mesh::Vector3& line_point0,
                             const chi_mesh::Vector3& line_point1,
                             chi_mesh::Vector3& intersection_point,
                             double* distance_to_intersection = nullptr);

bool CheckLineIntersectTriangle2(const chi_mesh::Vector3& tri_point0,
                                 const chi_mesh::Vector3& tri_point1,
                                 const chi_mesh::Vector3& tri_point2,
                                 const chi_mesh::Vector3& ray_posi,
                                 const chi_mesh::Vector3& ray_dir,
                                 chi_mesh::Vector3& intersection_point,
                                 double* distance_to_intersection = nullptr);

bool CheckPointInTriangle(const chi_mesh::Vector3& v0,
                          const chi_mesh::Vector3& v1,
                          const chi_mesh::Vector3& v2,
                          const chi_mesh::Normal& n,
                          const chi_mesh::Vector3& point);

bool CheckPlaneTetIntersect(const chi_mesh::Normal& plane_normal,
                            const chi_mesh::Vector3& plane_point,
                            const std::vector<chi_mesh::Vector3>& tet_points);

void PopulateRaySegmentLengths(const chi_mesh::MeshContinuum& grid,
                               const Cell& cell,
                               const chi_mesh::Vector3& line_point0,
                               const chi_mesh::Vector3& line_point1,
                               const chi_mesh::Vector3& omega,
                               std::vector<double>& segment_lengths);

} // namespace chi_mesh
