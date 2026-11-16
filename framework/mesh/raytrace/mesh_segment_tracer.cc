// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/raytrace/mesh_segment_tracer.h"
#include "framework/math/exact_predicates.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <limits>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>

namespace opensn
{

namespace
{

/// Cosine between the line and the plane of its exit triangle below which the crossing of the
/// plane is bounded by the chord of the line over the triangle rather than used as computed.
constexpr double grazing_cosine = 1.0e-6;
/// Fraction of the length of a segment within which it may leave the mesh before an end that
/// lies on a boundary face within the round-off of its coordinates, and the relative round-off
/// allowed in locating that end on the chord over the exit triangle.
constexpr double end_fraction = 1.0e-6;

/**
 * Sign of the orientation predicates of a line against the edges of a triangle when the line
 * threads the triangle in the direction of its outward normal, found once from a reference
 * configuration rather than derived by hand.
 */
int
OutwardThreadingSign()
{
  // Triangle in the plane z = 0, counterclockwise seen from +z, so its outward normal is +z.
  const Vector3 a(0.0, 0.0, 0.0);
  const Vector3 b(1.0, 0.0, 0.0);
  const Vector3 c(0.0, 1.0, 0.0);
  // A line through the triangle along +z.
  const Vector3 p(0.25, 0.25, -1.0);
  const Vector3 q(0.25, 0.25, 1.0);
  const int s1 = Orient3D(p, q, a, b);
  const int s2 = Orient3D(p, q, b, c);
  const int s3 = Orient3D(p, q, c, a);
  OpenSnLogicalErrorIf(s1 == 0 or s1 != s2 or s2 != s3,
                       "MeshSegmentTracer: inconsistent orientation predicates.");
  return s1;
}

const int outward_sign = OutwardThreadingSign();

struct TripleHash
{
  size_t operator()(const std::array<std::uint32_t, 3>& t) const
  {
    size_t h = t[0];
    h = h * 1000003U ^ t[1];
    h = h * 1000003U ^ t[2];
    return h;
  }
};

std::array<std::uint32_t, 3>
SortedTriple(std::uint32_t a, std::uint32_t b, std::uint32_t c)
{
  std::array<std::uint32_t, 3> t = {a, b, c};
  std::sort(t.begin(), t.end());
  return t;
}

} // namespace

MeshSegmentTracer::MeshSegmentTracer(const Mesh& grid, const bool split_all_cells)
{
  const auto& cells = grid.GetLocalCells();
  const auto num_cells = static_cast<std::uint32_t>(cells.size());
  cell_convex_.assign(num_cells, false);
  cell_size_.assign(num_cells, 0.0);
  cell_interior_.assign(num_cells, Vector3());
  cell_piece_offset_.assign(1, 0);
  piece_triangle_offset_.assign(1, 0);
  piece_edge_offset_.assign(1, 0);
  piece_face_offset_.assign(1, 0);

  // Compact mesh vertices.
  std::unordered_map<std::uint64_t, std::uint32_t> compact_vertex;
  for (const auto& cell : cells)
    for (const auto vid : cell->vertex_ids)
    {
      const auto [it, inserted] =
        compact_vertex.try_emplace(vid, static_cast<std::uint32_t>(points_.size()));
      if (inserted)
        points_.push_back(grid.GlobalVertex(vid));
    }
  const auto num_vertices = static_cast<std::uint32_t>(points_.size());

  // Face centroids, one point per face, computed from the vertices in a canonical order so that
  // both cells sharing the face get the same point.
  std::map<std::vector<std::uint64_t>, std::uint32_t> face_centroid;
  const auto FaceCentroid = [&](const std::vector<std::uint64_t>& vertex_ids)
  {
    std::vector<std::uint64_t> key(vertex_ids);
    std::sort(key.begin(), key.end());
    const auto it = face_centroid.find(key);
    if (it != face_centroid.end())
      return it->second;
    Vector3 sum;
    for (const auto vid : key)
      sum = sum + grid.GlobalVertex(vid);
    const Vector3 centroid = sum / static_cast<double>(key.size());
    points_.push_back(centroid);
    const auto id = static_cast<std::uint32_t>(points_.size() - 1);
    face_centroid.emplace(std::move(key), id);
    return id;
  };

  cell_vertex_offset_.assign(1, 0);
  std::vector<std::vector<std::uint32_t>> vertex_cells(num_vertices);
  std::unordered_map<std::array<std::uint32_t, 3>, std::vector<std::uint32_t>, TripleHash>
    triangle_pieces;

  // Boundary triangles of a cell, counterclockwise seen from outside, face by face, and the
  // index of the first triangle of each face.
  std::vector<std::array<std::uint32_t, 3>> cell_triangles;
  std::vector<std::uint32_t> cell_face_starts;
  for (const auto& cell : cells)
  {
    const auto c = static_cast<std::uint32_t>(cell->local_id);
    Vector3 lo = grid.GlobalVertex(cell->vertex_ids.front());
    Vector3 hi = lo;
    for (const auto vid : cell->vertex_ids)
    {
      const auto compact = compact_vertex.at(vid);
      cell_vertices_.push_back(compact);
      vertex_cells[compact].push_back(c);
      const auto& v = grid.GlobalVertex(vid);
      lo = Vector3(std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z));
      hi = Vector3(std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z));
    }
    cell_vertex_offset_.push_back(static_cast<std::uint32_t>(cell_vertices_.size()));
    cell_size_[c] = (hi - lo).Norm();

    cell_triangles.clear();
    cell_face_starts.clear();
    std::vector<std::uint32_t> cell_points(cell_vertices_.begin() + cell_vertex_offset_[c],
                                           cell_vertices_.end());
    for (const auto& face : cell->faces)
    {
      const size_t n = face.vertex_ids.size();
      OpenSnInvalidArgumentIf(n < 3, "MeshSegmentTracer: a face has fewer than three vertices.");
      cell_face_starts.push_back(static_cast<std::uint32_t>(cell_triangles.size()));
      if (n == 3)
        cell_triangles.push_back({compact_vertex.at(face.vertex_ids[0]),
                                  compact_vertex.at(face.vertex_ids[1]),
                                  compact_vertex.at(face.vertex_ids[2])});
      else
      {
        const auto f = FaceCentroid(face.vertex_ids);
        cell_points.push_back(f);
        for (size_t s = 0; s < n; ++s)
          cell_triangles.push_back({compact_vertex.at(face.vertex_ids[s]),
                                    compact_vertex.at(face.vertex_ids[(s + 1) % n]),
                                    f});
      }
    }

    // The cell is convex when every point of it lies in or on the half-space of every triangle.
    bool convex = not split_all_cells;
    for (const auto& tri : cell_triangles)
    {
      for (const auto pt : cell_points)
        if (Orient3D(points_[tri[0]], points_[tri[1]], points_[tri[2]], points_[pt]) < 0)
        {
          convex = false;
          break;
        }
      if (not convex)
        break;
    }

    // An interior point about which the cell is star-shaped: the average of the cell points, or
    // a point between it and a face centroid when the average lies on a triangle plane.
    Vector3 interior;
    for (const auto pt : cell_points)
      interior = interior + points_[pt];
    interior = interior / static_cast<double>(cell_points.size());
    const auto StarShapedAbout = [&](const Vector3& m)
    {
      for (const auto& tri : cell_triangles)
        if (Orient3D(points_[tri[0]], points_[tri[1]], points_[tri[2]], m) <= 0)
          return false;
      return true;
    };
    if (not StarShapedAbout(interior))
    {
      // The points about which the cell is star-shaped form its kernel, where the smallest signed
      // distance to the triangle planes is positive. That smallest distance is a concave function
      // of the point, so its maximum (the center of the largest ball in the kernel) is found by
      // climbing from the average point along the inward normal of the nearest plane and the
      // axes, with a step that halves when no direction improves. The result is verified exactly.
      std::vector<Vector3> inward_normals;
      std::vector<Vector3> anchors;
      for (const auto& tri : cell_triangles)
      {
        const Vector3& a = points_[tri[0]];
        const Vector3 n = (points_[tri[1]] - a).Cross(points_[tri[2]] - a);
        if (n.Norm() > 0.0)
        {
          inward_normals.push_back(-1.0 / n.Norm() * n);
          anchors.push_back(a);
        }
      }
      // The three nearest planes at a point: at a ridge or a corner of the function, the sums of
      // their inward normals are the directions along which it still increases.
      std::array<size_t, 3> nearest = {0, 0, 0};
      const auto SmallestDistance = [&](const Vector3& m)
      {
        std::array<double, 3> smallest = {std::numeric_limits<double>::infinity(),
                                          std::numeric_limits<double>::infinity(),
                                          std::numeric_limits<double>::infinity()};
        for (size_t k = 0; k < inward_normals.size(); ++k)
        {
          const double d = inward_normals[k].Dot(m - anchors[k]);
          for (size_t r = 0; r < 3; ++r)
            if (d < smallest[r])
            {
              for (size_t t = 2; t > r; --t)
              {
                smallest[t] = smallest[t - 1];
                nearest[t] = nearest[t - 1];
              }
              smallest[r] = d;
              nearest[r] = k;
              break;
            }
        }
        return smallest[0];
      };
      Vector3 center = interior;
      double value = SmallestDistance(center);
      double step = 0.25 * cell_size_[c];
      for (int iteration = 0; iteration < 1000 and step > 1.0e-12 * cell_size_[c]; ++iteration)
      {
        const Vector3 n0 = inward_normals[nearest[0]];
        const Vector3 n1 = inward_normals[nearest[1]];
        const Vector3 n2 = inward_normals[nearest[2]];
        bool improved = false;
        for (const Vector3& direction : {n0,
                                         n0 + n1,
                                         n0 + n1 + n2,
                                         n0 + n2,
                                         Vector3(1.0, 0.0, 0.0),
                                         Vector3(-1.0, 0.0, 0.0),
                                         Vector3(0.0, 1.0, 0.0),
                                         Vector3(0.0, -1.0, 0.0),
                                         Vector3(0.0, 0.0, 1.0),
                                         Vector3(0.0, 0.0, -1.0)})
        {
          // Opposite normals (a flat cell) sum to nothing.
          const double norm = direction.Norm();
          if (norm < 1.0e-6)
            continue;
          const Vector3 candidate = center + step / norm * direction;
          const double candidate_value = SmallestDistance(candidate);
          if (candidate_value > value)
          {
            center = candidate;
            value = candidate_value;
            improved = true;
          }
        }
        if (not improved)
          step *= 0.5;
      }
      OpenSnInvalidArgumentIf(not StarShapedAbout(center),
                              "MeshSegmentTracer: cell " + std::to_string(cell->global_id) +
                                " (centroid " + interior.PrintStr() +
                                ") is not star-shaped about any point, so it cannot be traced "
                                "through side tetrahedra.");
      interior = center;
    }
    cell_convex_[c] = convex;
    cell_interior_[c] = interior;

    const auto AddPiece = [&](const std::vector<std::array<std::uint32_t, 3>>& triangles,
                              const std::vector<std::uint32_t>& face_starts)
    {
      const auto piece = static_cast<std::uint32_t>(piece_cell_.size());
      const auto first_edge = static_cast<std::uint32_t>(edges_.size());
      const auto first_triangle = static_cast<std::uint32_t>(triangles_.size());
      const auto EdgeIndex = [&](const std::uint32_t u, const std::uint32_t v)
      {
        const std::pair<std::uint32_t, std::uint32_t> key = {std::min(u, v), std::max(u, v)};
        for (auto k = first_edge; k < edges_.size(); ++k)
          if (edges_[k] == key)
            return static_cast<std::uint32_t>(k - first_edge);
        edges_.push_back(key);
        return static_cast<std::uint32_t>(edges_.size() - 1 - first_edge);
      };
      for (const auto& tri : triangles)
      {
        Triangle t{tri[0], tri[1], tri[2], boundary, {}, {}};
        const std::array<std::pair<std::uint32_t, std::uint32_t>, 3> sides = {
          {{tri[0], tri[1]}, {tri[1], tri[2]}, {tri[2], tri[0]}}};
        for (size_t k = 0; k < 3; ++k)
        {
          t.edge[k] = EdgeIndex(sides[k].first, sides[k].second);
          t.forward[k] = sides[k].first < sides[k].second;
        }
        triangles_.push_back(t);
        triangle_pieces[SortedTriple(tri[0], tri[1], tri[2])].push_back(
          static_cast<std::uint32_t>(triangles_.size() - 1));
      }
      for (const auto start : face_starts)
      {
        const auto& t = triangles_[first_triangle + start];
        const Vector3& a = points_[t.a];
        const Vector3 normal = (points_[t.b] - a).Cross(points_[t.c] - a);
        face_first_triangle_.push_back(first_triangle + start);
        face_planes_.push_back({normal, normal.Dot(a)});
      }
      piece_triangle_offset_.push_back(static_cast<std::uint32_t>(triangles_.size()));
      piece_edge_offset_.push_back(static_cast<std::uint32_t>(edges_.size()));
      piece_face_offset_.push_back(static_cast<std::uint32_t>(face_first_triangle_.size()));
      piece_cell_.push_back(c);
      return piece;
    };

    if (convex)
      AddPiece(cell_triangles, cell_face_starts);
    else
    {
      points_.push_back(interior);
      const auto m = static_cast<std::uint32_t>(points_.size() - 1);
      for (const auto& tri : cell_triangles)
      {
        // The tetrahedron of the triangle and the interior point, every face oriented outward.
        std::vector<std::array<std::uint32_t, 3>> tet = {
          tri, {tri[0], m, tri[1]}, {tri[1], m, tri[2]}, {tri[2], m, tri[0]}};
        const std::array<std::uint32_t, 4> opposite = {m, tri[2], tri[0], tri[1]};
        for (size_t k = 0; k < 4; ++k)
          if (Orient3D(
                points_[tet[k][0]], points_[tet[k][1]], points_[tet[k][2]], points_[opposite[k]]) <
              0)
            std::swap(tet[k][1], tet[k][2]);
        AddPiece(tet, {0, 1, 2, 3});
      }
    }
    cell_piece_offset_.push_back(static_cast<std::uint32_t>(piece_cell_.size()));
  }

  // Pair the triangles: the two pieces on either side of a triangle hold the same three points.
  for (const auto& [key, tris] : triangle_pieces)
  {
    OpenSnLogicalErrorIf(tris.size() > 2,
                         "MeshSegmentTracer: a triangle is shared by more than two pieces.");
    if (tris.size() == 2)
    {
      const auto PieceOf = [&](const std::uint32_t t)
      {
        return static_cast<std::uint32_t>(
          std::upper_bound(piece_triangle_offset_.begin(), piece_triangle_offset_.end(), t) -
          piece_triangle_offset_.begin() - 1);
      };
      // The two pieces lie on opposite sides of the triangle, so they orient it oppositely.
      const auto& t0 = triangles_[tris[0]];
      const auto& t1 = triangles_[tris[1]];
      const bool opposite = (t1.a == t0.a and t1.b == t0.c and t1.c == t0.b) or
                            (t1.a == t0.b and t1.b == t0.a and t1.c == t0.c) or
                            (t1.a == t0.c and t1.b == t0.b and t1.c == t0.a);
      OpenSnInvalidArgumentIf(not opposite,
                              "MeshSegmentTracer: cells " +
                                std::to_string(cells[piece_cell_[PieceOf(tris[0])]]->global_id) +
                                " and " +
                                std::to_string(cells[piece_cell_[PieceOf(tris[1])]]->global_id) +
                                " lie on the same side of a face they share: the mesh is "
                                "inverted or self-overlapping.");
      triangles_[tris[0]].across = PieceOf(tris[1]);
      triangles_[tris[1]].across = PieceOf(tris[0]);
    }
  }

  for (const auto& point : points_)
    scale_ = std::max({scale_, std::abs(point.x), std::abs(point.y), std::abs(point.z)});

  vertex_cell_offset_.assign(1, 0);
  for (const auto& list : vertex_cells)
  {
    vertex_cells_.insert(vertex_cells_.end(), list.begin(), list.end());
    vertex_cell_offset_.push_back(static_cast<std::uint32_t>(vertex_cells_.size()));
  }
}

const MeshSegmentTracer::Triangle*
MeshSegmentTracer::Threaded(const std::uint32_t piece,
                            const Vector3& p,
                            const Vector3& q,
                            const bool outward) const
{
  // The faces of the piece in the order of the floating-point crossing parameters of their planes
  // along the line: for the exit, the faces whose plane the line crosses outward, first crossing
  // first; for the entry, the faces whose plane it crosses inward, last crossing first. The
  // triangles of the faces are verified exactly in that order, and the first that the line
  // threads in the wanted sense is the answer, which is unique for a convex piece.
  thread_local std::vector<std::pair<double, std::uint32_t>> candidates;
  thread_local std::vector<signed char> signs;
  const Vector3 d = q - p;
  const int wanted = outward ? outward_sign : -outward_sign;
  candidates.clear();
  const auto first_face = piece_face_offset_[piece];
  const auto last_face = piece_face_offset_[piece + 1];
  for (auto f = first_face; f < last_face; ++f)
  {
    const auto& plane = face_planes_[f];
    const double denominator = plane.normal.Dot(d);
    double s = std::numeric_limits<double>::infinity();
    if (outward ? denominator > 0.0 : denominator < 0.0)
    {
      s = (plane.offset - plane.normal.Dot(p)) / denominator;
      if (not std::isfinite(s))
        s = std::numeric_limits<double>::infinity();
      else if (not outward)
        s = -s;
    }
    candidates.emplace_back(s, f);
  }
  std::sort(candidates.begin(), candidates.end());
  const auto first_edge = piece_edge_offset_[piece];
  constexpr signed char unknown = 2;
  signs.assign(piece_edge_offset_[piece + 1] - first_edge, unknown);
  const auto Sign = [&](const Triangle& tri, const size_t e) -> int
  {
    auto& sign = signs[tri.edge[e]];
    if (sign == unknown)
    {
      const auto& edge = edges_[first_edge + tri.edge[e]];
      const int value = Orient3D(p, q, points_[edge.first], points_[edge.second]);
      if (value == 0)
        throw Degenerate{};
      sign = static_cast<signed char>(value);
    }
    return tri.forward[e] ? sign : -sign;
  };
  for (const auto& [s, f] : candidates)
  {
    const auto first = face_first_triangle_[f];
    const auto last =
      f + 1 < last_face ? face_first_triangle_[f + 1] : piece_triangle_offset_[piece + 1];
    for (auto t = first; t < last; ++t)
    {
      const auto& tri = triangles_[t];
      if (Sign(tri, 0) == wanted and Sign(tri, 1) == wanted and Sign(tri, 2) == wanted)
        return &tri;
    }
  }
  return nullptr;
}

const MeshSegmentTracer::Triangle&
MeshSegmentTracer::ExitTriangle(const std::uint32_t piece, const Vector3& p, const Vector3& q) const
{
  const Triangle* exit = Threaded(piece, p, q, true);
  OpenSnLogicalErrorIf(exit == nullptr,
                       "MeshSegmentTracer: a line leaves a convex piece through no triangle.");
  return *exit;
}

bool
MeshSegmentTracer::WithinRoundoff(const Triangle& triangle, const Vector3& point) const
{
  const auto& a = points_[triangle.a];
  const auto& b = points_[triangle.b];
  const auto& c = points_[triangle.c];
  const Vector3 normal = (b - a).Cross(c - a);
  const double norm = normal.Norm();
  if (norm == 0.0)
    return false;
  return std::abs(normal.Dot(point - a)) / norm <= std::max(Roundoff(a, point), Roundoff(b, c));
}

std::optional<std::uint32_t>
MeshSegmentTracer::FindStartPiece(const std::uint32_t cell,
                                  const Vector3& p,
                                  const Vector3& q,
                                  const bool search_neighbors) const
{
  // The segment travels in a piece that contains p when p lies strictly inside the plane of the
  // triangle through which the line leaves it. When p lies on that plane, the segment leaves the
  // piece at once and continues in the piece across the triangle, which also contains p.
  const auto TryPiece = [&](std::uint32_t piece) -> std::optional<std::uint32_t>
  {
    // The line crosses the convex piece between the triangles it threads inward and outward; p
    // lies in or on the piece exactly when it lies on the inner side of both of their planes.
    const Triangle* entry = Threaded(piece, p, q, false);
    if (entry == nullptr or
        Orient3D(points_[entry->a], points_[entry->b], points_[entry->c], p) < 0)
      return std::nullopt;
    for (size_t hop = 0; hop <= piece_cell_.size(); ++hop)
    {
      const Triangle* exit = Threaded(piece, p, q, true);
      if (exit == nullptr)
        return std::nullopt;
      const int side = Orient3D(points_[exit->a], points_[exit->b], points_[exit->c], p);
      if (side < 0)
        return std::nullopt;
      if (side > 0 or exit->across == boundary)
        return piece;
      // p lies on the exit plane: the segment leaves this piece at once and continues in the piece
      // across the triangle, which has p on its entry triangle.
      piece = exit->across;
    }
    OpenSnLogicalError(
      "MeshSegmentTracer: the start of a segment does not leave a cycle of pieces.");
  };
  const auto TryCell = [&](const std::uint32_t c) -> std::optional<std::uint32_t>
  {
    for (auto piece = cell_piece_offset_[c]; piece < cell_piece_offset_[c + 1]; ++piece)
      if (const auto found = TryPiece(piece))
        return found;
    return std::nullopt;
  };
  if (const auto piece = TryCell(cell))
    return piece;
  if (not search_neighbors)
    return std::nullopt;
  // The cells within two layers of the cell.
  std::vector<std::uint32_t> visited = {cell};
  std::vector<std::uint32_t> layer = {cell};
  for (int depth = 0; depth < 2; ++depth)
  {
    std::vector<std::uint32_t> next;
    for (const auto c : layer)
      for (auto k = cell_vertex_offset_[c]; k < cell_vertex_offset_[c + 1]; ++k)
      {
        const auto v = cell_vertices_[k];
        for (auto j = vertex_cell_offset_[v]; j < vertex_cell_offset_[v + 1]; ++j)
        {
          const auto candidate = vertex_cells_[j];
          if (std::find(visited.begin(), visited.end(), candidate) != visited.end())
            continue;
          visited.push_back(candidate);
          next.push_back(candidate);
          if (const auto piece = TryCell(candidate))
            return piece;
        }
      }
    layer = std::move(next);
  }
  return std::nullopt;
}

MeshSegmentTracer::Outcome
MeshSegmentTracer::Walk(const std::uint32_t start_cell,
                        const Vector3& p,
                        const Vector3& q,
                        std::vector<Segment>& segments,
                        const bool search_neighbors) const
{
  const auto start = FindStartPiece(start_cell, p, q, search_neighbors);
  if (not start)
    return Outcome::no_start;
  const double length = (q - p).Norm();
  const auto Add = [&](const std::uint32_t cell, const double piece_length)
  {
    if (not segments.empty() and segments.back().cell_local_id == cell)
      segments.back().length += piece_length;
    else
      segments.push_back({cell, piece_length});
  };
  std::uint32_t piece = *start;
  double t = 0.0;
  const size_t max_steps = piece_cell_.size() + 1;
  for (size_t step = 0; step < max_steps; ++step)
  {
    const auto cell = piece_cell_[piece];
    const auto& exit = ExitTriangle(piece, p, q);
    // The crossing of the plane of the exit triangle. For a triangle that the line grazes it is
    // ill-conditioned, but it lies on the chord of the line through the prism over the triangle,
    // which the walls through the edges along the triangle normal bound well; otherwise the
    // crossing is accurate and lies on the chord anyway.
    const Vector3& ea = points_[exit.a];
    const Vector3& eb = points_[exit.b];
    const Vector3& ec = points_[exit.c];
    const Vector3 normal = (eb - ea).Cross(ec - ea);
    const Vector3 d = q - p;
    const double denominator = normal.Dot(d);
    const double crossing =
      denominator == 0.0 ? t : std::clamp(normal.Dot(ea - p) / denominator * length, t, length);
    double lo = crossing;
    double hi = crossing;
    const bool grazing = std::abs(denominator) < grazing_cosine * normal.Norm() * length;
    // An end within the round-off of its coordinates of the plane of a boundary triangle, as an
    // end that the caller computed to lie on a boundary face can be.
    const bool end_near_boundary = exit.across == boundary and WithinRoundoff(exit, q);
    double chord_hi = length;
    if (grazing or end_near_boundary)
    {
      // The chord of the line through the prism over the triangle, bounded by the walls through
      // the edges along the triangle normal. For a grazed triangle the crossing is ill-conditioned
      // but lies on the chord; otherwise the crossing is accurate and lies on the chord anyway.
      double chord_lo = t;
      const Vector3 direction = d / length;
      const std::array<std::pair<const Vector3*, const Vector3*>, 3> edges = {
        {{&ea, &eb}, {&eb, &ec}, {&ec, &ea}}};
      const std::array<const Vector3*, 3> opposite = {&ec, &ea, &eb};
      for (size_t e = 0; e < 3; ++e)
      {
        Vector3 wall = (*edges[e].second - *edges[e].first).Cross(normal);
        if (wall.Dot(*opposite[e] - *edges[e].first) > 0.0)
          wall = -1.0 * wall;
        const double slope = wall.Dot(direction);
        const double offset = wall.Dot(p - *edges[e].first);
        // Inside the wall where offset + s * slope <= 0. A wall that the line runs along within
        // round-off bounds nothing; its quotient would be noise.
        const double parallel = 1.0e-9 * wall.Norm();
        if (slope > parallel)
          chord_hi = std::min(chord_hi, -offset / slope);
        else if (slope < -parallel)
          chord_lo = std::max(chord_lo, -offset / slope);
      }
      if (not std::isfinite(chord_lo) or not std::isfinite(chord_hi) or chord_hi < chord_lo)
        chord_hi = chord_lo = std::clamp(std::isfinite(chord_lo) ? chord_lo : t, t, length);
      if (grazing)
      {
        lo = chord_lo;
        hi = chord_hi;
      }
    }
    // The line entered the piece before the end of the segment, so the piece, convex, holds the
    // end exactly when the end lies on or inside the plane of the exit triangle. This decision is
    // exact: the floating-point crossing can fall short of an end that lies inside the piece by
    // its round-off, which exceeds the length of a short segment.
    if (Orient3D(ea, eb, ec, q) >= 0)
    {
      Add(cell, length - t);
      return Outcome::reached;
    }
    const double stop = std::clamp(std::isfinite(crossing) ? crossing : lo, lo, hi);
    if (end_near_boundary and (chord_hi >= length * (1.0 - end_fraction) or
                               length - stop <= std::max(end_fraction * length, Roundoff(p, q))))
    {
      // The end lies outside the plane of a boundary triangle by round-off, and either over that
      // triangle, since the chord of the line over it reaches the end, or within round-off of the
      // crossing point along the line, which lies in the triangle: the segment ends where it
      // leaves. Nearness to the plane alone is not enough: a line that crosses the plane at a
      // shallow angle stays within round-off of it far beyond the triangle. The chord alone is not
      // enough either: the crossing point can lie next to the edge of the triangle with the end
      // over the neighboring triangle of the same face.
      Add(cell, length - t);
      return Outcome::reached;
    }
    Add(cell, stop - t);
    t = stop;
    if (exit.across == boundary)
      return end_near_boundary ? Outcome::left_near_end : Outcome::left;
    piece = exit.across;
  }
  OpenSnLogicalError("MeshSegmentTracer: a segment did not terminate after " +
                     std::to_string(max_steps) + " pieces.");
}

double
MeshSegmentTracer::Roundoff(const Vector3& p, const Vector3& q) const
{
  double scale = scale_;
  for (const double c : {p.x, p.y, p.z, q.x, q.y, q.z})
    scale = std::max(scale, std::abs(c));
  return 64.0 * std::ldexp(1.0, std::ilogb(std::max(scale, 1.0e-300)) - 52);
}

Vector3
MeshSegmentTracer::Displacement(const Vector3& p,
                                const Vector3& q,
                                const int direction,
                                const bool far_end) const
{
  const double sx = (direction & 1) != 0 ? -1.0 : 1.0;
  const double sy = (direction & 2) != 0 ? -1.0 : 1.0;
  const double sz = (direction & 4) != 0 ? -1.0 : 1.0;
  const Vector3 unit = far_end ? Vector3(sx * 3.0, sy * 1.0, sz * 2.0).Normalized()
                               : Vector3(sx * 1.0, sy * 2.0, sz * 3.0).Normalized();
  // Twice (three times for the far end) the round-off allowed to a point computed by a caller,
  // so that a start or end that lies outside the mesh by that round-off is moved inside for some
  // direction.
  return (far_end ? 3.0 : 2.0) * Roundoff(p, q) * unit;
}

std::optional<bool>
MeshSegmentTracer::Resolve(const std::uint32_t start_cell,
                           const Vector3& start,
                           const Vector3& end,
                           std::vector<Segment>& segments) const
{
  if ((end - start).Norm() == 0.0)
    return true;
  std::vector<Segment> traced;
  std::optional<std::vector<Segment>> exact_left;
  try
  {
    // The start cell, or the cell across one of its faces, is expected to hold the start.
    Outcome outcome = Walk(start_cell, start, end, traced, false);
    if (outcome == Outcome::no_start)
    {
      // A start computed by the caller to lie on a face of the start cell can fall outside it by
      // round-off: move it toward the interior of the start cell by the round-off allowance.
      const Vector3 inward = (cell_interior_[start_cell] - start).Normalized();
      const Vector3 p = start + 2.0 * Roundoff(start, end) * inward;
      traced.clear();
      outcome = Walk(start_cell, p, end, traced, false);
      const double factor = (end - start).Norm() / (end - p).Norm();
      for (auto& segment : traced)
        segment.length *= factor;
    }
    if (outcome == Outcome::no_start)
    {
      // A start cell given by a caller that located it approximately: the cells near it.
      traced.clear();
      outcome = Walk(start_cell, start, end, traced, true);
    }
    // The segment is not degenerate: its outcome is final, unless it left the mesh through a
    // boundary triangle whose plane holds the end within round-off. A segment that runs along a
    // boundary face within round-off of it can fall outside the triangles of that face, whose
    // centroid is rounded, over a long stretch; it reaches its end if some displacement of its
    // ends by that round-off keeps it inside, as decided below, and otherwise left as traced.
    if (outcome != Outcome::no_start and outcome != Outcome::left_near_end)
    {
      segments.insert(segments.end(), traced.begin(), traced.end());
      return outcome == Outcome::reached;
    }
    if (outcome == Outcome::left_near_end)
      exact_left = traced;
    traced.clear();
  }
  catch (const Degenerate&)
  {
    // The line meets an edge or a vertex, or lies in a face plane: the traversals of the displaced
    // segment below decide.
    traced.clear();
  }
  // The line of the segment meets an edge or a vertex, or lies in a face plane. The ends are
  // displaced independently, by a few units in the last place, in every pair of octant
  // directions, which moves the line off the degeneracy. A displaced end may fall outside the
  // mesh where the segment ends on the boundary, so the segment reaches its end if any
  // displacement does, and otherwise leaves the mesh where the displacement that stays inside
  // longest does: a segment on the boundary counts as inside.
  std::optional<std::vector<Segment>> left;
  double left_length = 0.0;
  for (int i = 0; i < num_directions; ++i)
    for (int j = 0; j < num_directions; ++j)
    {
      const Vector3 p = start + Displacement(start, end, i, false);
      const Vector3 q = end + Displacement(start, end, j, true);
      traced.clear();
      Outcome outcome = Outcome::no_start;
      try
      {
        outcome = Walk(start_cell, p, q, traced, true);
      }
      catch (const Degenerate&)
      {
        continue;
      }
      // The displaced segment differs in length by a few units in the last place.
      const double factor = (end - start).Norm() / (q - p).Norm();
      double traced_length = 0.0;
      for (auto& segment : traced)
      {
        segment.length *= factor;
        traced_length += segment.length;
      }
      if (outcome == Outcome::reached)
      {
        segments.insert(segments.end(), traced.begin(), traced.end());
        return true;
      }
      if ((outcome == Outcome::left or outcome == Outcome::left_near_end) and
          (not left or traced_length > left_length))
      {
        left = traced;
        left_length = traced_length;
      }
    }
  if (exact_left)
  {
    // No displacement keeps the segment inside: it left where the exact traversal says.
    segments.insert(segments.end(), exact_left->begin(), exact_left->end());
    return false;
  }
  if (left)
  {
    segments.insert(segments.end(), left->begin(), left->end());
    return false;
  }
  return std::nullopt;
}

bool
MeshSegmentTracer::Trace(const std::uint32_t start_cell,
                         const Vector3& start,
                         const Vector3& end,
                         std::vector<Segment>& segments) const
{
  const auto reached = Resolve(start_cell, start, end, segments);
  OpenSnLogicalErrorIf(not reached,
                       "MeshSegmentTracer: the start " + start.PrintStr() + " of the segment to " +
                         end.PrintStr() + " does not lie in, on, or near its start cell " +
                         std::to_string(start_cell) + ".");
  return *reached;
}

} // namespace opensn
