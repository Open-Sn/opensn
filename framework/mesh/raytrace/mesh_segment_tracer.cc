// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/raytrace/mesh_segment_tracer.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/logging/log.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>

namespace opensn
{

namespace
{

/// Directions within this cosine of a plane are treated as parallel to it.
constexpr double parallel_cosine = 1.0e-14;
/// Tolerance of geometric tests relative to the size of the cell or of the segment.
constexpr double relative_tolerance = 1.0e-10;
/// Allowance for the round-off of a plane evaluation relative to the magnitude of the
/// coordinates, about 450 times the machine epsilon.
constexpr double coordinate_roundoff = 1.0e-13;
/// Shortest altitude of a side triangle, relative to the size of its cell, below which the
/// direction of its plane normal, whose error is the round-off of the vertices over the
/// altitude, is too uncertain for the tolerance at the far side of the cell.
constexpr double thin_altitude = 1.0e-3;

/// Largest magnitude of a coordinate of the points.
double
CoordinateScale(std::initializer_list<Vector3> points)
{
  double scale = 0.0;
  for (const auto& p : points)
    scale = std::max({scale, std::abs(p.x), std::abs(p.y), std::abs(p.z)});
  return scale;
}

struct Plane
{
  Vector3 normal;
  double offset;
};

/// Plane through p, q, r with its normal on the side away from inside, or nullopt if the three
/// points are collinear.
std::optional<Plane>
PlaneAwayFrom(const Vector3& p, const Vector3& q, const Vector3& r, const Vector3& inside)
{
  Vector3 normal = (q - p).Cross(r - p);
  const double norm = normal.Norm();
  if (norm == 0.0)
    return std::nullopt;
  normal = normal / norm;
  if (normal.Dot(inside - p) > 0.0)
    normal = -1.0 * normal;
  return Plane{normal, normal.Dot(p)};
}

/// Shortest altitude of the triangle p, q, r: twice its area over its longest edge.
double
ShortestAltitude(const Vector3& p, const Vector3& q, const Vector3& r)
{
  const double longest = std::max({(q - p).Norm(), (r - q).Norm(), (p - r).Norm()});
  return longest > 0.0 ? (q - p).Cross(r - p).Norm() / longest : 0.0;
}

} // namespace

MeshSegmentTracer::MeshSegmentTracer(const Mesh& grid, const bool split_all_cells)
{
  const auto& cells = grid.GetLocalCells();
  const size_t num_cells = cells.size();
  OpenSnLogicalErrorIf(num_cells >= no_neighbor, "MeshSegmentTracer: too many local cells.");
  cell_convex_.assign(num_cells, false);
  cell_size_.assign(num_cells, 0.0);
  cell_piece_offset_.assign(1, 0);
  piece_plane_offset_.assign(1, 0);

  std::unordered_map<uint64_t, std::uint32_t> compact_vertex;
  cell_vertex_offset_.assign(1, 0);

  for (const auto& cell_ptr : cells)
  {
    const Cell& cell = *cell_ptr;
    const auto c = static_cast<std::uint32_t>(cell.local_id);
    OpenSnLogicalErrorIf(c != cell_vertex_offset_.size() - 1,
                         "MeshSegmentTracer: local cells are not ordered by local id.");

    // Compact vertex indices and the cell size.
    Vector3 lo = grid.GlobalVertex(cell.vertex_ids.front());
    Vector3 hi = lo;
    for (const auto vid : cell.vertex_ids)
    {
      const auto [it, inserted] =
        compact_vertex.try_emplace(vid, static_cast<std::uint32_t>(compact_vertex.size()));
      cell_vertices_.push_back(it->second);
      const auto& v = grid.GlobalVertex(vid);
      lo = Vector3(std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z));
      hi = Vector3(std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z));
    }
    cell_vertex_offset_.push_back(static_cast<std::uint32_t>(cell_vertices_.size()));
    const double size = (hi - lo).Norm();
    cell_size_[c] = size;
    const double tolerance = relative_tolerance * std::max(size, 1.0e-300) +
                             coordinate_roundoff * CoordinateScale({lo, hi});

    // The planes of the side triangles of every face, oriented outward, with their neighbors.
    // A planar face contributes its plane once.
    std::vector<Plane> planes;
    std::vector<std::uint32_t> neighbors;
    for (const auto& face : cell.faces)
    {
      std::uint32_t neighbor = no_neighbor;
      if (face.has_neighbor and grid.IsCellLocal(face.neighbor_id))
        neighbor = static_cast<std::uint32_t>(grid.GetGlobalCell(face.neighbor_id).local_id);
      const size_t num_sides = face.vertex_ids.size();
      const size_t first_plane = planes.size();
      // The sides in the order of the shortest altitude of their triangles, so that a planar
      // face keeps the plane of its best conditioned side triangle.
      std::vector<std::pair<double, size_t>> sides;
      for (size_t s = 0; s < num_sides; ++s)
      {
        const auto& a = grid.GlobalVertex(face.vertex_ids[s]);
        const auto& b = grid.GlobalVertex(face.vertex_ids[(s + 1) % num_sides]);
        sides.emplace_back(-ShortestAltitude(a, b, face.centroid), s);
      }
      std::sort(sides.begin(), sides.end());
      for (const auto& [negative_altitude, s] : sides)
      {
        const auto& a = grid.GlobalVertex(face.vertex_ids[s]);
        const auto& b = grid.GlobalVertex(face.vertex_ids[(s + 1) % num_sides]);
        const auto plane = PlaneAwayFrom(a, b, face.centroid, cell.centroid);
        if (not plane)
          continue;
        const bool repeated =
          std::any_of(planes.begin() + static_cast<std::ptrdiff_t>(first_plane),
                      planes.end(),
                      [&](const Plane& other)
                      {
                        return (other.normal - plane->normal).Norm() <= relative_tolerance and
                               std::abs(other.offset - plane->offset) <= tolerance;
                      });
        if (repeated)
          continue;
        if (-negative_altitude < thin_altitude * size)
          ++num_thin_side_triangles_;
        planes.push_back(*plane);
        neighbors.push_back(neighbor);
      }
    }

    // The cell is convex when every vertex and face centroid lies on the inner side of every
    // plane.
    bool convex = true;
    for (const auto& plane : planes)
    {
      for (const auto vid : cell.vertex_ids)
        convex = convex and plane.normal.Dot(grid.GlobalVertex(vid)) - plane.offset <= tolerance;
      for (const auto& face : cell.faces)
        convex = convex and plane.normal.Dot(face.centroid) - plane.offset <= tolerance;
      if (not convex)
        break;
    }
    convex = convex and not split_all_cells;
    cell_convex_[c] = convex;

    if (convex)
    {
      for (size_t k = 0; k < planes.size(); ++k)
      {
        plane_normal_.push_back(planes[k].normal);
        plane_offset_.push_back(planes[k].offset);
        plane_neighbor_.push_back(neighbors[k]);
      }
      piece_plane_offset_.push_back(static_cast<std::uint32_t>(plane_normal_.size()));
    }
    else
    {
      // Side tetrahedra: a side triangle of a face and the cell centroid. Only the side
      // triangle lies on the cell boundary. They fill the cell only if the cell centroid lies on
      // the inner side of every side triangle, which the face vertex order (counterclockwise
      // seen from outside) defines. Otherwise a side tetrahedron lies outside the cell, as does
      // the corresponding sub-element of the piecewise-linear discretizations.
      for (size_t f_index = 0; f_index < cell.faces.size(); ++f_index)
      {
        const auto& face = cell.faces[f_index];
        std::uint32_t neighbor = no_neighbor;
        if (face.has_neighbor and grid.IsCellLocal(face.neighbor_id))
          neighbor = static_cast<std::uint32_t>(grid.GetGlobalCell(face.neighbor_id).local_id);
        const size_t num_sides = face.vertex_ids.size();
        for (size_t s = 0; s < num_sides; ++s)
        {
          const auto& a = grid.GlobalVertex(face.vertex_ids[s]);
          const auto& b = grid.GlobalVertex(face.vertex_ids[(s + 1) % num_sides]);
          const auto& f = face.centroid;
          const auto& m = cell.centroid;
          const Vector3 outward = (a - f).Cross(b - f);
          // The height of the cell centroid above the plane of the side triangle, times the
          // norm of outward, positive on the outer side.
          const double centroid_height = outward.Dot(m - a);
          OpenSnInvalidArgumentIf(centroid_height > tolerance * outward.Norm(),
                                  "MeshSegmentTracer: cell " + std::to_string(cell.global_id) +
                                    " (centroid " + m.PrintStr() + "), face " +
                                    std::to_string(f_index) + ", side " + std::to_string(s) +
                                    ": the cell centroid lies on the outer side of the plane of "
                                    "the side triangle, so the cell is not star-shaped about its "
                                    "centroid.");
          // A side tetrahedron whose apex lies in the plane of its side triangle, or whose side
          // triangle is degenerate, has no volume. The other side tetrahedra fill the cell
          // without it. Its planes would coincide, and a ray in that plane would lie in it
          // without an entry or an exit.
          if (centroid_height >= -tolerance * outward.Norm())
            continue;
          const std::array<std::array<Vector3, 4>, 4> triangles = {
            {{a, b, f, m}, {a, b, m, f}, {a, f, m, b}, {b, f, m, a}}};
          for (size_t k = 0; k < 4; ++k)
          {
            const auto& tri = triangles[k];
            const auto plane = PlaneAwayFrom(tri[0], tri[1], tri[2], tri[3]);
            OpenSnLogicalErrorIf(not plane,
                                 "MeshSegmentTracer: a face of a side tetrahedron with volume "
                                 "is degenerate.");
            if (ShortestAltitude(tri[0], tri[1], tri[2]) < thin_altitude * size)
              ++num_thin_side_triangles_;
            plane_normal_.push_back(plane->normal);
            plane_offset_.push_back(plane->offset);
            plane_neighbor_.push_back(k == 0 ? neighbor : no_neighbor);
          }
          piece_plane_offset_.push_back(static_cast<std::uint32_t>(plane_normal_.size()));
        }
      }
    }
    cell_piece_offset_.push_back(static_cast<std::uint32_t>(piece_plane_offset_.size() - 1));
  }

  // Cells of each vertex.
  vertex_cell_offset_.assign(compact_vertex.size() + 1, 0);
  for (const auto v : cell_vertices_)
    ++vertex_cell_offset_[v + 1];
  for (size_t v = 1; v < vertex_cell_offset_.size(); ++v)
    vertex_cell_offset_[v] += vertex_cell_offset_[v - 1];
  vertex_cells_.resize(cell_vertices_.size());
  std::vector<std::uint32_t> fill(vertex_cell_offset_.begin(), vertex_cell_offset_.end() - 1);
  for (std::uint32_t c = 0; c < num_cells; ++c)
    for (auto k = cell_vertex_offset_[c]; k < cell_vertex_offset_[c + 1]; ++k)
      vertex_cells_[fill[cell_vertices_[k]]++] = c;

  if (num_thin_side_triangles_ > 0)
    log.LogAllWarning() << "MeshSegmentTracer: " << num_thin_side_triangles_
                        << " side triangles of faces are thinner than " << thin_altitude
                        << " of their cell, so the planes computed from them are less accurate "
                           "than the tracing tolerance. A segment that crosses such a plane may "
                           "have a small part assigned to the neighboring cell.";
}

void
MeshSegmentTracer::PieceInterval(const std::uint32_t piece,
                                 const Ray& ray,
                                 const double t,
                                 double& entry,
                                 double& exit,
                                 std::vector<std::uint32_t>* exit_neighbors) const
{
  // Whether the point of the ray at t lies in a half-space is decided from its signed distance
  // to the plane, which rounds off like the coordinates. The value of t at which the ray crosses
  // the plane is that distance divided by the cosine, so for a plane the ray grazes it is far
  // less accurate, and it is only used for planes that the point lies clear of.
  const double tol = ray.tolerance;
  const Vector3 point = ray.start + t * ray.omega;
  entry = -std::numeric_limits<double>::infinity();
  exit = std::numeric_limits<double>::infinity();
  const auto first = piece_plane_offset_[piece];
  const auto last = piece_plane_offset_[piece + 1];
  for (auto k = first; k < last; ++k)
  {
    const double cosine = plane_normal_[k].Dot(ray.omega);
    // The distance from the point to the plane, positive on the inner side.
    const double room = plane_offset_[k] - plane_normal_[k].Dot(point);
    if (room < -tol)
    {
      // The point lies outside the half-space: the ray enters it later, or never.
      if (cosine < -parallel_cosine)
        entry = std::max(entry, t + room / cosine);
      else
      {
        entry = std::numeric_limits<double>::infinity();
        exit = -std::numeric_limits<double>::infinity();
        return;
      }
    }
    else if (room <= tol)
    {
      // The point lies on the plane: the ray leaves or enters the half-space at t, or runs
      // along the plane.
      if (cosine > parallel_cosine)
        exit = std::min(exit, t);
      else if (cosine < -parallel_cosine)
        entry = std::max(entry, t);
    }
    else
    {
      // The point lies inside the half-space.
      if (cosine > parallel_cosine)
        exit = std::min(exit, t + room / cosine);
      else if (cosine < -parallel_cosine)
        entry = std::max(entry, t + room / cosine);
    }
  }
  if (exit_neighbors and exit > entry and std::isfinite(exit))
  {
    // The neighbors across the planes on which the exit point lies.
    const Vector3 exit_point = ray.start + exit * ray.omega;
    for (auto k = first; k < last; ++k)
      if (plane_neighbor_[k] != no_neighbor and
          plane_normal_[k].Dot(ray.omega) > parallel_cosine and
          plane_offset_[k] - plane_normal_[k].Dot(exit_point) <= tol)
        exit_neighbors->push_back(plane_neighbor_[k]);
  }
}

std::optional<double>
MeshSegmentTracer::ExitAfter(const std::uint32_t cell,
                             const Ray& ray,
                             const double t,
                             std::vector<std::uint32_t>* exit_neighbors) const
{
  const double tol = ray.tolerance;
  const auto first = cell_piece_offset_[cell];
  const auto last = cell_piece_offset_[cell + 1];
  if (exit_neighbors)
    exit_neighbors->clear();

  if (cell_convex_[cell])
  {
    double entry = 0.0;
    double exit = 0.0;
    PieceInterval(first, ray, t, entry, exit, exit_neighbors);
    if (entry <= t + tol and exit > t + tol)
      return exit;
    if (exit_neighbors)
      exit_neighbors->clear();
    return std::nullopt;
  }

  // A non-convex cell: join the intervals of the side tetrahedra that continue the ray from t.
  // Each piece is decided at the point that the joined interval has reached, so that a piece
  // entered through a plane that the ray grazes is found from the distance to that plane.
  const auto num_pieces = last - first;
  double reach = t;
  for (std::uint32_t pass = 0;; ++pass)
  {
    OpenSnLogicalErrorIf(pass > num_pieces,
                         "MeshSegmentTracer: the side tetrahedra of a cell did not stop extending "
                         "the interval of a ray.");
    bool extended = false;
    for (std::uint32_t p = 0; p < num_pieces; ++p)
    {
      double entry = 0.0;
      double exit = 0.0;
      PieceInterval(first + p, ray, reach, entry, exit, nullptr);
      if (entry <= reach + tol and exit > reach + tol)
      {
        reach = exit;
        extended = true;
      }
    }
    if (not extended)
      break;
  }
  if (reach <= t + tol)
    return std::nullopt;
  if (exit_neighbors)
    for (std::uint32_t p = 0; p < num_pieces; ++p)
    {
      // The pieces that the reach lies in or on, which the ray leaves there.
      double entry = 0.0;
      double exit = 0.0;
      PieceInterval(first + p, ray, reach, entry, exit, nullptr);
      if (entry <= reach + tol and exit >= reach - tol)
        PieceInterval(first + p, ray, reach, entry, exit, exit_neighbors);
    }
  return reach;
}

std::optional<std::uint32_t>
MeshSegmentTracer::FindNext(const std::uint32_t cell,
                            const Ray& ray,
                            const double t,
                            const std::vector<std::uint32_t>& exit_neighbors,
                            std::vector<std::uint32_t>& visited,
                            double& exit) const
{
  std::optional<std::uint32_t> best;
  double best_exit = t + ray.tolerance;
  visited.clear();
  const auto Consider = [&](const std::uint32_t candidate)
  {
    if (std::find(visited.begin(), visited.end(), candidate) != visited.end())
      return;
    visited.push_back(candidate);
    const auto candidate_exit = ExitAfter(candidate, ray, t, nullptr);
    if (candidate_exit and *candidate_exit > best_exit)
    {
      best = candidate;
      best_exit = *candidate_exit;
    }
  };

  // The neighbors across the faces through which the ray leaves, which is all that is needed
  // unless the ray passes through an edge or vertex.
  for (const auto neighbor : exit_neighbors)
    Consider(neighbor);
  if (not best)
  {
    // A cell that the ray enters at an edge or vertex of this cell shares a vertex with it.
    for (auto k = cell_vertex_offset_[cell]; k < cell_vertex_offset_[cell + 1]; ++k)
    {
      const auto v = cell_vertices_[k];
      for (auto j = vertex_cell_offset_[v]; j < vertex_cell_offset_[v + 1]; ++j)
        Consider(vertex_cells_[j]);
    }
  }
  if (not best)
  {
    // Only a non-conforming mesh needs a second layer. Consider appends to visited, so the
    // first layer is its leading part.
    const size_t first_layer_size = visited.size();
    for (size_t i = 0; i < first_layer_size; ++i)
    {
      const auto layer_cell = visited[i];
      for (auto k = cell_vertex_offset_[layer_cell]; k < cell_vertex_offset_[layer_cell + 1]; ++k)
      {
        const auto v = cell_vertices_[k];
        for (auto j = vertex_cell_offset_[v]; j < vertex_cell_offset_[v + 1]; ++j)
          Consider(vertex_cells_[j]);
      }
    }
  }
  if (best)
    exit = best_exit;
  return best;
}

bool
MeshSegmentTracer::Trace(const std::uint32_t start_cell,
                         const Vector3& start,
                         const Vector3& end,
                         std::vector<Segment>& segments) const
{
  const double length = (end - start).Norm();
  if (length == 0.0)
    return true;
  // Values of t within the tolerance are equal: a cell may be entered or left that much early
  // or late. The tolerance covers the round-off of the plane evaluations at the coordinates of
  // the segment, and a small fraction of the geometry so that slivers and repeated planes are
  // tolerated.
  const double tolerance = relative_tolerance * std::max(length, cell_size_[start_cell]) +
                           coordinate_roundoff * CoordinateScale({start, end});
  const Ray ray{start, (end - start) / length, tolerance};

  // The start must lie in or on the start cell.
  {
    bool contains = false;
    for (auto p = cell_piece_offset_[start_cell];
         p < cell_piece_offset_[start_cell + 1] and not contains;
         ++p)
    {
      bool inside = true;
      for (auto k = piece_plane_offset_[p]; k < piece_plane_offset_[p + 1]; ++k)
        inside = inside and plane_normal_[k].Dot(start) - plane_offset_[k] <= ray.tolerance;
      contains = inside;
    }
    OpenSnLogicalErrorIf(
      not contains, "MeshSegmentTracer: the start of a segment does not lie in its start cell.");
  }

  std::vector<std::uint32_t> exit_neighbors;
  std::vector<std::uint32_t> visited;
  double t = 0.0;
  double exit = 0.0;
  std::uint32_t cell = start_cell;
  if (const auto start_exit = ExitAfter(cell, ray, t, &exit_neighbors))
    exit = *start_exit;
  else
  {
    // The segment only touches the start cell. One that is shorter than the tolerance, and
    // that no cell continues, lies in the start cell up to the tolerance.
    const auto next = FindNext(cell, ray, t, exit_neighbors, visited, exit);
    if (not next)
    {
      if (length <= ray.tolerance)
      {
        segments.push_back({cell, length});
        return true;
      }
      return false;
    }
    cell = *next;
    ExitAfter(cell, ray, t, &exit_neighbors);
  }

  const size_t max_steps = 4 * (piece_plane_offset_.size() + 4);
  for (size_t step = 0; step < max_steps; ++step)
  {
    const double stop = std::min(exit, length);
    segments.push_back({cell, stop - t});
    t = stop;
    if (t >= length - ray.tolerance)
    {
      segments.back().length += length - t;
      return true;
    }
    const auto next = FindNext(cell, ray, t, exit_neighbors, visited, exit);
    if (not next)
      return false;
    cell = *next;
    ExitAfter(cell, ray, t, &exit_neighbors);
  }
  OpenSnLogicalError("MeshSegmentTracer: a segment did not terminate after " +
                     std::to_string(max_steps) + " cells.");
}

std::optional<std::uint32_t>
MeshSegmentTracer::LocateAhead(const std::uint32_t hint_cell,
                               const Vector3& point,
                               const Vector3& omega) const
{
  const double tolerance =
    relative_tolerance * cell_size_[hint_cell] + coordinate_roundoff * CoordinateScale({point});
  const Ray ray{point, omega, tolerance};
  if (ExitAfter(hint_cell, ray, 0.0, nullptr))
    return hint_cell;
  double exit = 0.0;
  std::vector<std::uint32_t> visited;
  return FindNext(hint_cell, ray, 0.0, {}, visited, exit);
}

} // namespace opensn
