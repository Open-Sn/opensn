// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/raytrace/segment_tracer.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>

namespace opensn
{

namespace
{

/// Directions within this cosine of a plane are treated as parallel to it.
constexpr double parallel_cosine = 1.0e-14;
/// Tolerance of geometric tests relative to the size of the cell or of the coordinates.
constexpr double relative_tolerance = 1.0e-10;

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

} // namespace

SegmentTracer::SegmentTracer(const Mesh& grid, const bool split_all_cells)
{
  const auto& cells = grid.GetLocalCells();
  const size_t num_cells = cells.size();
  OpenSnLogicalErrorIf(num_cells >= no_neighbor, "SegmentTracer: too many local cells.");
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
                         "SegmentTracer: local cells are not ordered by local id.");

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
    const double tolerance = relative_tolerance * std::max(size, 1.0e-300);

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
      for (size_t s = 0; s < num_sides; ++s)
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
      // triangle lies on the cell boundary.
      for (const auto& face : cell.faces)
      {
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
          const std::array<std::array<Vector3, 4>, 4> triangles = {
            {{a, b, f, m}, {a, b, m, f}, {a, f, m, b}, {b, f, m, a}}};
          std::array<std::optional<Plane>, 4> tet_planes;
          bool degenerate = false;
          for (size_t k = 0; k < 4; ++k)
          {
            const auto& tri = triangles[k];
            tet_planes[k] = PlaneAwayFrom(tri[0], tri[1], tri[2], tri[3]);
            degenerate = degenerate or not tet_planes[k];
          }
          if (degenerate)
            continue;
          for (size_t k = 0; k < 4; ++k)
          {
            plane_normal_.push_back(tet_planes[k]->normal);
            plane_offset_.push_back(tet_planes[k]->offset);
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
}

void
SegmentTracer::PieceInterval(const std::uint32_t piece,
                             const Ray& ray,
                             double& entry,
                             double& exit,
                             std::vector<std::uint32_t>* exit_neighbors,
                             const double exit_tolerance) const
{
  entry = -std::numeric_limits<double>::infinity();
  exit = std::numeric_limits<double>::infinity();
  const auto first = piece_plane_offset_[piece];
  const auto last = piece_plane_offset_[piece + 1];
  for (auto k = first; k < last; ++k)
  {
    const double cosine = plane_normal_[k].Dot(ray.omega);
    const double room = plane_offset_[k] - plane_normal_[k].Dot(ray.start);
    if (std::abs(cosine) <= parallel_cosine)
    {
      // The ray is parallel to the plane: inside the piece's half-space or never.
      if (room < -ray.tolerance)
      {
        entry = std::numeric_limits<double>::infinity();
        exit = -std::numeric_limits<double>::infinity();
        return;
      }
      continue;
    }
    const double t = room / cosine;
    if (cosine > 0.0)
      exit = std::min(exit, t);
    else
      entry = std::max(entry, t);
  }
  if (exit_neighbors and exit > entry)
    for (auto k = first; k < last; ++k)
    {
      const double cosine = plane_normal_[k].Dot(ray.omega);
      if (cosine > parallel_cosine and plane_neighbor_[k] != no_neighbor)
      {
        const double t = (plane_offset_[k] - plane_normal_[k].Dot(ray.start)) / cosine;
        if (t <= exit + exit_tolerance)
          exit_neighbors->push_back(plane_neighbor_[k]);
      }
    }
}

std::optional<double>
SegmentTracer::ExitAfter(const std::uint32_t cell,
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
    PieceInterval(first, ray, entry, exit, exit_neighbors, tol);
    if (entry <= t + tol and exit > t + tol)
      return exit;
    if (exit_neighbors)
      exit_neighbors->clear();
    return std::nullopt;
  }

  // A non-convex cell: join the intervals of its side tetrahedra that continue from t.
  const auto num_pieces = last - first;
  std::vector<double> entries(num_pieces);
  std::vector<double> exits(num_pieces);
  for (std::uint32_t p = 0; p < num_pieces; ++p)
    PieceInterval(first + p, ray, entries[p], exits[p], nullptr, tol);
  double reach = t;
  bool extended = true;
  while (extended)
  {
    extended = false;
    for (std::uint32_t p = 0; p < num_pieces; ++p)
      if (entries[p] <= reach + tol and exits[p] > reach + tol)
      {
        reach = exits[p];
        extended = true;
      }
  }
  if (reach <= t + tol)
    return std::nullopt;
  if (exit_neighbors)
    for (std::uint32_t p = 0; p < num_pieces; ++p)
      if (exits[p] > entries[p] and std::abs(exits[p] - reach) <= tol)
      {
        double entry = 0.0;
        double exit = 0.0;
        PieceInterval(first + p, ray, entry, exit, exit_neighbors, tol);
      }
  return reach;
}

std::optional<std::uint32_t>
SegmentTracer::FindNext(const std::uint32_t cell,
                        const Ray& ray,
                        const double t,
                        const std::vector<std::uint32_t>& exit_neighbors,
                        double& exit) const
{
  std::optional<std::uint32_t> best;
  double best_exit = t + ray.tolerance;
  std::vector<std::uint32_t> visited;
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
SegmentTracer::Trace(const std::uint32_t start_cell,
                     const Vector3& start,
                     const Vector3& end,
                     std::vector<Segment>& segments) const
{
  const double length = (end - start).Norm();
  if (length == 0.0)
    return true;
  const double scale = std::max({std::abs(start.x),
                                 std::abs(start.y),
                                 std::abs(start.z),
                                 std::abs(end.x),
                                 std::abs(end.y),
                                 std::abs(end.z),
                                 length,
                                 cell_size_[start_cell]});
  const Ray ray{start, (end - start) / length, relative_tolerance * scale};

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
    OpenSnLogicalErrorIf(not contains,
                         "SegmentTracer: the start of a segment does not lie in its start cell.");
  }

  std::vector<std::uint32_t> exit_neighbors;
  double t = 0.0;
  double exit = 0.0;
  std::uint32_t cell = start_cell;
  if (const auto start_exit = ExitAfter(cell, ray, t, &exit_neighbors))
    exit = *start_exit;
  else
  {
    // The segment only touches the start cell.
    const auto next = FindNext(cell, ray, t, exit_neighbors, exit);
    if (not next)
      return false;
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
    const auto next = FindNext(cell, ray, t, exit_neighbors, exit);
    if (not next)
      return false;
    cell = *next;
    ExitAfter(cell, ray, t, &exit_neighbors);
  }
  OpenSnLogicalError("SegmentTracer: a segment did not terminate after " +
                     std::to_string(max_steps) + " cells.");
}

std::optional<std::uint32_t>
SegmentTracer::LocateAhead(const std::uint32_t hint_cell,
                           const Vector3& point,
                           const Vector3& omega) const
{
  const double scale =
    std::max({std::abs(point.x), std::abs(point.y), std::abs(point.z), cell_size_[hint_cell]});
  const Ray ray{point, omega, relative_tolerance * scale};
  if (ExitAfter(hint_cell, ray, 0.0, nullptr))
    return hint_cell;
  double exit = 0.0;
  return FindNext(hint_cell, ray, 0.0, {}, exit);
}

} // namespace opensn
