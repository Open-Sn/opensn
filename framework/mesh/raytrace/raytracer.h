// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/mesh.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_vector.h"

namespace opensn
{

class Cell;
class MeshContinuum;

/// Data structure to hold output info from the raytracer.
struct RayTracerOutputInformation
{
  double distance_to_surface = 0.0;
  Vector3 pos_f;
  unsigned int destination_face_index = 0;
  uint64_t destination_face_neighbor = 0;
  bool particle_lost = false;
  std::string lost_particle_info;
};

/// Raytracer object.
class RayTracer
{
private:
  const MeshContinuum& reference_grid_;
  std::vector<double> cell_sizes_;
  double epsilon_nudge_ = 1.0e-8;
  double backward_tolerance_ = 1.0e-10;
  double extension_distance_ = 1.0e5;
  bool perform_concavity_checks_ = true;

public:
  explicit RayTracer(const MeshContinuum& grid,
                     std::vector<double> cell_sizes,
                     bool perform_concavity_checks = true)
    : reference_grid_(grid),
      cell_sizes_(std::move(cell_sizes)),
      perform_concavity_checks_(perform_concavity_checks)
  {
  }

private:
  const MeshContinuum& Grid() const;

  void SetTolerancesFromCellSize(double cell_size)
  {
    epsilon_nudge_ = cell_size * 1.0e-2;
    backward_tolerance_ = cell_size * 1.0e-10;
    extension_distance_ = 3.0 * cell_size;
  }

public:
  /**
   * Traces a ray with an initial position either within the cell or on the cell surface, and with a
   * direction vector pointing inward toward the cell. If the ray-trace fails the particle will be
   * marked as lost.
   */
  RayTracerOutputInformation
  TraceRay(const Cell& cell, Vector3& pos_i, Vector3& omega_i, int function_depth = 0);

  /// Traces a ray with an initial position, presumed to be outside the cell, to an incident face.
  RayTracerOutputInformation
  TraceIncidentRay(const Cell& cell, const Vector3& pos_i, const Vector3& omega_i);

private:
  /// Performs raytracing within a 1D-slab.
  void TraceSlab(const Cell& cell,
                 Vector3& pos_i,
                 Vector3& omega_i,
                 bool& intersection_found,
                 bool& backward_tolerance_hit,
                 RayTracerOutputInformation& oi);

  /// Performs raytracing within a 2D Polygon.
  void TracePolygon(const Cell& cell,
                    Vector3& pos_i,
                    Vector3& omega_i,
                    bool& intersection_found,
                    bool& backward_tolerance_hit,
                    RayTracerOutputInformation& oi);

  /// Performs raytracing within a 3D Polyhedron.
  void TracePolyhedron(const Cell& cell,
                       Vector3& pos_i,
                       Vector3& omega_i,
                       bool& intersection_found,
                       bool& backward_tolerance_hit,
                       RayTracerOutputInformation& oi);
};

/**
 * Computes the intersection of a line with a plane.
 *
 * The first step of this algorithm is to compute v0 and v1. These
 * are vectors from the plane's reference point to each of the line-points,
 * respectively.
 * We then take the dot-products of these
 * vectors with the plane normal.
 * We then say that the vectors have a positive sense if the dot-product
 * is positive and a negative sense if the dot-product is negative.
 * If the senses are not equal then the line intersects the plane.
 *
 * Since the face normal is a normalized vector the dot-product of v0 or v1
 * will give the projection of the relevant vector along the normal to the
 * plane. We can use this projection to compute a weight associated with
 * each vector. This also then allows us to compute the intersection point.
 *
 * \param plane_normal The normal associated with the plane
 * \param plane_point The reference point for the plane
 * \param line_point_0 The line's initial point
 * \param line_point_1 The line's destination point
 * \param intersection_point The point to be populated with the intersection
 *                           point
 * \param weights The weights associated with this intersection
 *
 * \return Returns true if the line intersects the plane and false otherwise.
 */
bool CheckPlaneLineIntersect(const Vector3& plane_normal,
                             const Vector3& plane_point,
                             const Vector3& line_point_0,
                             const Vector3& line_point_1,
                             Vector3& intersection_point,
                             std::pair<double, double>* weights = nullptr);

/**
 * Given a strip defined by two points (v0,v1) and a normal, n, (meaning infinite in the direction
 * defined by (v1-v0).cross(n), this function determines if a line, defined from p0 to p1,
 * intersects it. If it does then `true` is returned and `intersection_point` contains the point of
 * intersection. If it does not then `false` is returned and `intersection_point` remains unchanged.
 */
bool CheckLineIntersectStrip(const Vector3& strip_point0,
                             const Vector3& strip_point1,
                             const Vector3& strip_normal,
                             const Vector3& line_point0,
                             const Vector3& line_point1,
                             Vector3& intersection_point,
                             double* distance_to_intersection = nullptr);

/// Given a triangle defined by three points, computes whether a line intersects this triangle.
bool CheckLineIntersectTriangle2(const Vector3& tri_point0,
                                 const Vector3& tri_point1,
                                 const Vector3& tri_point2,
                                 const Vector3& ray_posi,
                                 const Vector3& ray_dir,
                                 Vector3& intersection_point,
                                 double* distance_to_intersection = nullptr);

/// Check whether a point lies in a triangle.
bool CheckPointInTriangle(
  const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& n, const Vector3& point);

/**
 * This functions checks the intersection of a plane with a tetrahedron.
 * The equation of a plane is
 *      nx(x-x0) + ny(y-y0) + nz(z-z0) = 0
 * Where the plane normal is (nx,ny,nz) and the plane point is (x0,y0,z0).
 * If we form a dot product between the normal and a vector (x-x0,y-y0,z-z0) then sign of the result
 * gives the sense to the surface. Therefore, if we encounter differing senses then the plane is
 * indeed intersecting.
 */
bool CheckPlaneTetIntersect(const Vector3& plane_normal,
                            const Vector3& plane_point,
                            const std::vector<Vector3>& tet_points);

/// Populates segment lengths along a ray. Sorted along the direction.
void PopulateRaySegmentLengths(const MeshContinuum& grid,
                               const Cell& cell,
                               const Vector3& line_point0,
                               const Vector3& line_point1,
                               const Vector3& omega,
                               std::vector<double>& segment_lengths);

} // namespace opensn
