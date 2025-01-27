// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/raytrace/raytracer.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/cell/cell.h"
#include "framework/logging/log.h"
#include <algorithm>
#include <set>

namespace opensn
{

const MeshContinuum&
RayTracer::Grid() const
{
  return reference_grid_;
}

RayTracerOutputInformation
RayTracer::TraceRay(const Cell& cell, Vector3& pos_i, Vector3& omega_i, int function_depth)
{
  if (not cell_sizes_.empty())
    SetTolerancesFromCellSize(cell_sizes_[cell.local_id]);

  RayTracerOutputInformation oi;

  bool intersection_found = false;
  bool backward_tolerance_hit = false;

  if (cell.GetType() == CellType::SLAB)
    TraceSlab(cell, pos_i, omega_i, intersection_found, backward_tolerance_hit, oi);
  else if (cell.GetType() == CellType::POLYGON)
    TracePolygon(cell, pos_i, omega_i, intersection_found, backward_tolerance_hit, oi);
  else if (cell.GetType() == CellType::POLYHEDRON)
    TracePolyhedron(cell, pos_i, omega_i, intersection_found, backward_tolerance_hit, oi);
  else
    throw std::logic_error("Unsupported cell type encountered in call to "
                           "RayTrace.");

  if (not intersection_found)
  {
    if (function_depth < 5)
    {
      // Nudge particle towards centroid
      Vector3 v_p_i_cc = (cell.centroid - pos_i);
      Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge_;

      oi = TraceRay(cell, pos_i_nudged, omega_i, function_depth + 1);

      return oi;
    }

    if (function_depth < 7)
    {
      // Nudge particle away from line between location and cell center
      Vector3 v_p_i_cc = (cell.centroid - pos_i).Cross(omega_i);
      Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge_;

      oi = TraceRay(cell, pos_i_nudged, omega_i, function_depth + 1);

      return oi;
    }

    std::stringstream outstr;

    outstr << "Intersection not found at function level " << function_depth << "."
           << ((backward_tolerance_hit) ? " Backward tolerance hit. " : "")
           << "For particle xyz=" << pos_i.PrintStr() << " uvw=" << omega_i.PrintStr() << " "
           << (pos_i + extension_distance_ * omega_i).PrintStr() << " " << extension_distance_
           << " in cell " << cell.global_id << " with vertices: \n";

    const auto& grid = Grid();

    for (auto vi : cell.vertex_ids)
      outstr << grid.vertices[vi].PrintStr() << "\n";

    for (auto& face : cell.faces)
    {
      outstr << "Face with centroid: " << face.centroid.PrintStr() << " ";
      outstr << "n=" << face.normal.PrintStr() << "\n";
      for (auto vi : face.vertex_ids)
        outstr << grid.vertices[vi].PrintStr() << "\n";
    }

    outstr << "o Cell\n";
    for (auto& vid : cell.vertex_ids)
    {
      auto& v = grid.vertices[vid];
      outstr << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    for (auto& face : cell.faces)
    {
      auto& v = face.centroid;
      outstr << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    for (size_t f = 0; f < cell.faces.size(); ++f)
    {
      auto& face = cell.faces[f];
      outstr << "f ";
      for (auto vid : face.vertex_ids)
      {
        size_t ref_cell_id = 0;
        for (size_t cid = 0; cid < cell.vertex_ids.size(); ++cid)
          if (cell.vertex_ids[cid] == vid)
            ref_cell_id = cid + 1;

        outstr << ref_cell_id << "// ";
      }
      outstr << "\n";
    }

    oi.particle_lost = true;
    oi.lost_particle_info = outstr.str();
  }

  return oi;
}

RayTracerOutputInformation
RayTracer::TraceIncidentRay(const Cell& cell, const Vector3& pos_i, const Vector3& omega_i)
{
  const auto cell_type = cell.GetType();
  const double cell_char_length = cell_sizes_[cell.local_id];
  const auto& grid = reference_grid_;

  bool intersects_cell = false;
  Vector3 I;

  size_t f = 0;
  for (const auto& face : cell.faces)
  {
    if (face.normal.Dot(omega_i) > 0.0)
    {
      ++f;
      continue /*the loop*/;
    }

    const auto& p0 = grid.vertices[face.vertex_ids[0]];
    const auto& n = face.normal;

    const auto ppos_i = p0 - pos_i;
    const double d = ppos_i.Dot(omega_i);

    const auto pos_ext = pos_i + (d + cell_char_length) * omega_i;

    if (cell_type == CellType::SLAB)
    {
      intersects_cell = CheckPlaneLineIntersect(n, p0, pos_i, pos_ext, I);
    } // SLAB
    else if (cell_type == CellType::POLYGON)
    {
      const auto& p1 = grid.vertices[face.vertex_ids[1]];
      intersects_cell = CheckLineIntersectStrip(p0, p1, n, pos_i, pos_ext, I);
    } // POLYGON
    else if (cell_type == CellType::POLYHEDRON)
    {
      const auto& vids = face.vertex_ids;
      const size_t num_sides = face.vertex_ids.size();
      for (size_t s = 0; s < num_sides; ++s)
      {
        uint64_t v0i = vids[s];
        uint64_t v1i = (s < (num_sides - 1)) ? vids[s + 1] : vids[0];

        const auto& v0 = grid.vertices[v0i];
        const auto& v1 = grid.vertices[v1i];
        const auto& v2 = face.centroid;

        const auto v01 = v1 - v0;
        const auto v02 = v2 - v0;
        const auto n_est = v01.Cross(v02);

        if (n_est.Dot(omega_i) > 0.0)
          continue;

        intersects_cell = CheckLineIntersectTriangle2(v0, v1, v2, pos_i, omega_i, I);
        if (intersects_cell)
          break;
      } // for side
    }   // POLYHEDRON

    if (intersects_cell)
      break;
    ++f;
  } // for face

  RayTracerOutputInformation oi;
  if (intersects_cell)
  {
    oi.distance_to_surface = (I - pos_i).Norm();
    oi.pos_f = I;
    oi.destination_face_index = f;
    oi.destination_face_neighbor = cell.global_id;
    oi.particle_lost = false;
  }
  else
    oi.particle_lost = true;

  return oi;
}

void
RayTracer::TraceSlab(const Cell& cell,
                     Vector3& pos_i,
                     Vector3& omega_i,
                     bool& intersection_found,
                     bool& backward_tolerance_hit,
                     RayTracerOutputInformation& oi)
{
  const auto& grid = Grid();
  Vector3 intersection_point;
  std::pair<double, double> weights;

  const double fabs_mu = std::fabs(omega_i.Dot(cell.faces[0].normal));

  double d_extend = (fabs_mu < 1.0e-15) ? 1.0e15 : extension_distance_ / fabs_mu;

  Vector3 pos_f_line = pos_i + omega_i * d_extend;

  int num_faces = 2;
  for (int f = 0; f < num_faces; ++f)
  {
    uint64_t fpi = cell.vertex_ids[f]; // face point index
    Vector3 face_point = grid.vertices[fpi];

    bool intersects = CheckPlaneLineIntersect(
      cell.faces[f].normal, face_point, pos_i, pos_f_line, intersection_point, &weights);

    double D = weights.first * d_extend;

    if ((D > backward_tolerance_) and intersects)
    {
      oi.distance_to_surface = D;
      oi.pos_f = intersection_point;

      oi.destination_face_index = f;
      oi.destination_face_neighbor = cell.faces[f].neighbor_id;
      intersection_found = true;
      break;
    }
    if (intersects)
      backward_tolerance_hit = true;
  } // for faces
}

void
RayTracer::TracePolygon(const Cell& cell,
                        Vector3& pos_i,
                        Vector3& omega_i,
                        bool& intersection_found,
                        bool& backward_tolerance_hit,
                        RayTracerOutputInformation& oi)
{
  const auto& grid = Grid();
  Vector3 ip; // intersection point

  const double fabs_mu = std::fabs(omega_i.Dot(cell.faces[0].normal));

  double d_extend = (fabs_mu < 1.0e-15) ? 1.0e15 : extension_distance_ / fabs_mu;

  Vector3 pos_f_line = pos_i + omega_i * d_extend;

  std::vector<RayTracerOutputInformation> face_intersections;

  size_t num_faces = cell.faces.size();
  face_intersections.reserve(num_faces);
  for (int f = 0; f < num_faces; ++f)
  {
    if (cell.faces[f].normal.Dot(omega_i) < 0.0)
      continue;

    RayTracerOutputInformation face_oi;

    uint64_t fpi = cell.faces[f].vertex_ids[0]; // face point index 0
    uint64_t fpf = cell.faces[f].vertex_ids[1]; // face point index 1
    const Vector3& face_point_i = grid.vertices[fpi];
    const Vector3& face_point_f = grid.vertices[fpf];

    bool intersects = CheckLineIntersectStrip(
      face_point_i, face_point_f, cell.faces[f].normal, pos_i, pos_f_line, ip);

    double D = (ip - pos_i).Norm();

    if ((D > backward_tolerance_) and intersects)
    {
      face_oi.distance_to_surface = D;
      face_oi.pos_f = ip;

      face_oi.destination_face_index = f;
      face_oi.destination_face_neighbor = cell.faces[f].neighbor_id;
      intersection_found = true;
      face_intersections.emplace_back(std::move(face_oi));
      if (not perform_concavity_checks_)
        break;
    } // if intersects
    if ((D < backward_tolerance_) and intersects)
      backward_tolerance_hit = true;
  } // for faces

  // Determine closest intersection
  if (not perform_concavity_checks_ and not face_intersections.empty())
    oi = face_intersections.back();
  else if (perform_concavity_checks_ and not face_intersections.empty())
  {
    auto closest_intersection = &face_intersections.back();
    for (auto& intersection : face_intersections)
      if (intersection.distance_to_surface < closest_intersection->distance_to_surface)
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

void
RayTracer::TracePolyhedron(const Cell& cell,
                           Vector3& pos_i,
                           Vector3& omega_i,
                           bool& intersection_found,
                           bool& backward_tolerance_hit,
                           RayTracerOutputInformation& oi)
{
  const auto& grid = Grid();
  Vector3 ip = pos_i; // Intersection point

  std::vector<RayTracerOutputInformation> triangle_intersections;

  size_t num_faces = cell.faces.size();
  triangle_intersections.reserve(num_faces * 4); // Guessing 4 tris per face
  for (int f = 0; f < num_faces; ++f)
  {
    const auto& face = cell.faces[f];

    size_t num_sides = face.vertex_ids.size();
    for (size_t s = 0; s < num_sides; ++s)
    {
      size_t v0_index = face.vertex_ids[s];
      size_t v1_index = (s < (num_sides - 1)) ? // if not last vertex
                          face.vertex_ids[s + 1]
                                              : face.vertex_ids[0]; // else
      const auto& v0 = grid.vertices[v0_index];
      const auto& v1 = grid.vertices[v1_index];
      const auto& v2 = cell.faces[f].centroid;

      auto v01 = v1 - v0;
      auto v02 = v2 - v0;
      auto n_est = v01.Cross(v02);

      if (n_est.Dot(omega_i) < 0.0)
        continue;

      RayTracerOutputInformation triangle_oi;

      bool intersects = CheckLineIntersectTriangle2(v0, v1, v2, pos_i, omega_i, ip);

      if (intersects)
      {
        triangle_oi.distance_to_surface = (ip - pos_i).Norm();
        triangle_oi.pos_f = ip;

        triangle_oi.destination_face_index = f;
        triangle_oi.destination_face_neighbor = cell.faces[f].neighbor_id;

        intersection_found = true;
        triangle_intersections.emplace_back(std::move(triangle_oi));
        if (not perform_concavity_checks_)
          break;
      } // if intersects
    }   // for side

    if (intersection_found and (not perform_concavity_checks_))
      break;
  } // for faces

  // Determine closest intersection
  if (not perform_concavity_checks_ and not triangle_intersections.empty())
    oi = triangle_intersections.back();
  else if (perform_concavity_checks_ and not triangle_intersections.empty())
  {
    auto closest_intersection = &triangle_intersections.back();
    for (auto& intersection : triangle_intersections)
      if (intersection.distance_to_surface < closest_intersection->distance_to_surface)
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

bool
CheckPlaneLineIntersect(const Vector3& plane_normal,
                        const Vector3& plane_point,
                        const Vector3& line_point_0,
                        const Vector3& line_point_1,
                        Vector3& intersection_point,
                        std::pair<double, double>* weights)
{
  Vector3 v0 = line_point_0 - plane_point;
  Vector3 v1 = line_point_1 - plane_point;

  double dotp_0 = plane_normal.Dot(v0);
  double dotp_1 = plane_normal.Dot(v1);

  bool sense_0 = (dotp_0 >= 0.0);
  bool sense_1 = (dotp_1 >= 0.0);

  if (sense_0 != sense_1)
  {
    double dotp_total = std::fabs(dotp_0) + std::fabs(dotp_1);
    double w0 = (std::fabs(dotp_0) / dotp_total);
    double w1 = 1.0 - w0;
    intersection_point = line_point_0 * w1 + line_point_1 * w0;

    if (weights != nullptr)
      *weights = {w0, w1};
    return true;
  }

  return false;
}

bool
CheckLineIntersectStrip(const Vector3& strip_point0,
                        const Vector3& strip_point1,
                        const Vector3& strip_normal,
                        const Vector3& line_point0,
                        const Vector3& line_point1,
                        Vector3& intersection_point,
                        double* distance_to_intersection)
{
  Vector3 plane_intersection_point;
  std::pair<double, double> weights;

  bool intersects_plane = CheckPlaneLineIntersect(
    strip_normal, strip_point0, line_point0, line_point1, plane_intersection_point, &weights);

  if (not intersects_plane)
    return false;

  Vector3 edge_vec = strip_point1 - strip_point0;
  Vector3 ints_vec1 = plane_intersection_point - strip_point0;
  Vector3 ints_vec2 = plane_intersection_point - strip_point1;

  bool sense1 = edge_vec.Dot(ints_vec1) >= 0.0;
  bool sense2 = edge_vec.Dot(ints_vec2) >= 0.0;

  if (distance_to_intersection != nullptr)
    *distance_to_intersection = (plane_intersection_point - line_point0).Norm();

  if (sense1 != sense2)
  {
    intersection_point = plane_intersection_point;

    return true;
  }

  return false;
}

bool
CheckLineIntersectTriangle2(const Vector3& tri_point0,
                            const Vector3& tri_point1,
                            const Vector3& tri_point2,
                            const Vector3& ray_posi,
                            const Vector3& ray_dir,
                            Vector3& intersection_point,
                            double* distance_to_intersection)
{
  double epsilon = 1.0e-12;
  Vector3 edge1 = tri_point1 - tri_point0;
  Vector3 edge2 = tri_point2 - tri_point0;

  // Compute characteristic vector for incident angle
  // This vector becomes perpendicular to the plane
  // when the ray is parallel to triangle
  Vector3 h = ray_dir.Cross(edge2);

  // If h is indeed perpendicular to the plane,
  // the dot product of the other leg with this h will be close
  // to zero.
  double a = edge1.Dot(h);
  if (std::fabs(a) < epsilon)
    return false;

  Vector3 s = ray_posi - tri_point0;

  double f = 1.0 / a;

  // If, s projected onto h, is greater than,
  // v01 projected onto h, or negative, there is now way
  // the ray can intersect
  double u = f * (s.Dot(h));
  if (u < 0.0 or u > 1.0)
    return false;

  Vector3 q = s.Cross(edge1);

  // If, q projected onto omega, is greater than,
  // v01 projected onto h, or negative, there is now way
  // the ray can intersect
  double v = f * ray_dir.Dot(q);
  if (v < 0.0 or (u + v) > 1.0)
    return false;

  double t = f * edge2.Dot(q);

  if (distance_to_intersection != nullptr)
    *distance_to_intersection = t;

  if (t > epsilon and t < (1.0 / epsilon))
  {
    intersection_point = ray_posi + ray_dir * t;
    return true;
  }
  else
  {
    return false;
  }
}

bool
CheckPointInTriangle(
  const Vector3& v0, const Vector3& v1, const Vector3& v2, const Vector3& n, const Vector3& point)
{
  auto v01 = v1 - v0;
  auto v12 = v2 - v1;
  auto v20 = v0 - v2;

  auto v0p = point - v0;
  auto v1p = point - v1;
  auto v2p = point - v2;

  auto vc0 = v01.Cross(v0p);
  auto vc1 = v12.Cross(v1p);
  auto vc2 = v20.Cross(v2p);

  bool dp0 = (vc0.Dot(n) >= 0.0);
  bool dp1 = (vc1.Dot(n) >= 0.0);
  bool dp2 = (vc2.Dot(n) >= 0.0);

  if (dp0 and dp1 and dp2)
    return true;
  else
    return false;
}

bool
CheckPlaneTetIntersect(const Vector3& plane_normal,
                       const Vector3& plane_point,
                       const std::vector<Vector3>& tet_points)
{
  bool current_sense = false;

  size_t num_points = tet_points.size();
  for (size_t i = 0; i < num_points; ++i)
  {
    Vector3 v = tet_points[i] - plane_point;
    double dotp = plane_normal.Dot(v);

    bool new_sense = (dotp >= 0.0);

    if (i == 0)
      current_sense = new_sense;
    else if (new_sense != current_sense)
      return true;
  }
  return false;
}

void
PopulateRaySegmentLengths(const MeshContinuum& grid,
                          const Cell& cell,
                          const Vector3& line_point0,
                          const Vector3& line_point1,
                          const Vector3& omega,
                          std::vector<double>& segment_lengths)
{
  const Vector3 khat(0, 0, 1);
  std::set<double> distance_set;

  double track_length;
  if (segment_lengths.empty())
  {
    track_length = (line_point1 - line_point0).Norm();
    segment_lengths.push_back(track_length);
  }

  track_length = segment_lengths.front();
  distance_set.insert(track_length);

  // Determine intersection points
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
  // Since there are no segments within a slab we will only have
  // a single segment length. It is already pushed

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
  // A polygon can be decomposed into "sides" by means of its
  // edges. Each side comprises a triangle formed by: the two
  // vertices of the associated edge v0 and v1, and the cell
  // centroid vc.
  // Since the triangles all share an edge we only determine
  // segment lengths from the strip defined by v0 to vc.
  if (cell.GetType() == CellType::POLYGON)
  {
    for (auto& face : cell.faces) // edges
    {
      const auto& v0 = grid.vertices[face.vertex_ids[0]];
      const auto& vc = cell.centroid;

      auto n0 = (vc - v0).Cross(khat).Normalized();

      Vector3 intersection_point;
      double d = 0.0;
      bool intersects =
        CheckLineIntersectStrip(v0, vc, n0, line_point0, line_point1, intersection_point, &d);

      if (intersects)
      {
        //        double d = (intersection_point - line_point0).Norm();
        distance_set.insert(d);
      }

    } // for face
  }
  else if (cell.GetType() == CellType::POLYHEDRON)
  {
    auto& vcc = cell.centroid;

    for (auto& face : cell.faces)
    {
      auto& vfc = face.centroid;

      // Face center to vertex segments
      for (auto vi : face.vertex_ids)
      {
        auto& vert = grid.vertices[vi];

        Vector3 intersection_point;

        double d = 0.0;
        bool intersects =
          CheckLineIntersectTriangle2(vert, vfc, vcc, line_point0, omega, intersection_point, &d);

        if (intersects)
        {
          if (d < track_length)
            distance_set.insert(d);
        }
      } // for edge

      // Face edge to cell center segments
      for (int v = 0; v < face.vertex_ids.size(); ++v)
      {
        uint64_t vid_0 = face.vertex_ids[v];
        uint64_t vid_1 =
          (v < (face.vertex_ids.size() - 1)) ? face.vertex_ids[v + 1] : face.vertex_ids[0];

        auto& v0 = grid.vertices[vid_0];
        auto& v1 = grid.vertices[vid_1];
        auto& v2 = vcc;

        Vector3 intersection_point;

        double d = 0.0;
        bool intersects =
          CheckLineIntersectTriangle2(v0, v1, v2, line_point0, omega, intersection_point, &d);

        if (intersects)
        {
          if (d < track_length)
            distance_set.insert(d);
        }
      } // for edge
    }   // for face
  }

  // Populate segment lengths
  // if there are N segments intersected then there will always be
  // N+1 distances.
  segment_lengths.clear();
  double last_distance = 0.0;
  for (double dl : distance_set)
  {
    double new_seg_length = dl - last_distance;
    last_distance = dl;
    segment_lengths.push_back(new_seg_length);
  }
}

} // namespace opensn
