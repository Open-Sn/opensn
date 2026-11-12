// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace opensn
{

class Mesh;

/**
 * Exact traversal of straight segments through the local cells of a mesh.
 *
 * Every cell is the union of its side tetrahedra, each formed by a side triangle of a face
 * (two consecutive face vertices and the face centroid) and the cell centroid. A cell for which
 * all of its vertices and face centroids lie on the inner side of the plane of every side
 * triangle is convex, and is the intersection of the half-spaces of its side triangles. Any
 * other cell is traced through its side tetrahedra, which are convex.
 *
 * A segment from start to end is parameterized by its distance t from start along its fixed
 * direction. For each convex piece the tracer computes the exact interval of t in which the line
 * lies in the piece, from the planes of the piece alone. It never moves the ray or uses a point
 * computed in another cell. A cell continues the traversal at t when its interval begins at t,
 * within a tolerance scaled to the coordinates, and extends beyond it. The next cell is looked
 * for first among the neighbors across the faces through which the ray leaves, then among all
 * cells that share a vertex with the current cell, and then among the cells that share a vertex
 * with those. A ray that leaves a cell exactly through an edge or vertex, or that only touches a
 * cell, is therefore continued in the cell it enters. Each step advances t, and the lengths of
 * a completed traversal sum to the segment length exactly.
 *
 * The tracer works on the local cells of a mesh and treats a face whose neighbor is not local as
 * a boundary face.
 */
class SegmentTracer
{
public:
  /// A part of a segment in one local cell.
  struct Segment
  {
    std::uint32_t cell_local_id;
    double length;
  };

  /**
   * Builds the planes of every local cell and the vertex-to-cell adjacency. O(cells).
   *
   * With split_all_cells, every cell is traced through its side tetrahedra, as a non-convex cell
   * is. The result is the same up to round-off; the option exists to test that path on meshes of
   * convex cells.
   */
  explicit SegmentTracer(const Mesh& grid, bool split_all_cells = false);

  /**
   * Traces the segment from start to end, beginning in the local cell start_cell, which must
   * contain start or have it on its boundary. Appends the cells that the segment crosses, in
   * order, with the length of the segment in each, to segments.
   *
   * Returns true when the segment reaches end. Returns false when the segment leaves the local
   * cells first, through a boundary face, a face whose neighbor is not local, or a gap in a
   * non-conforming mesh; segments then holds the part traced before the segment left. Throws
   * if start does not lie in or on start_cell.
   */
  bool Trace(std::uint32_t start_cell,
             const Vector3& start,
             const Vector3& end,
             std::vector<Segment>& segments) const;

  /**
   * Returns the local cell that contains the points just after point along the unit direction
   * omega: hint_cell if it does, else one of the cells near it. Returns std::nullopt if no local
   * cell does, for example when the ray leaves the mesh at point.
   */
  std::optional<std::uint32_t>
  LocateAhead(std::uint32_t hint_cell, const Vector3& point, const Vector3& omega) const;

  /// Whether the cell is traced as one convex piece rather than through its side tetrahedra.
  bool IsConvexCell(std::uint32_t cell_local_id) const { return cell_convex_[cell_local_id]; }

private:
  /// The line start + t omega, with the absolute tolerance used to compare values of t.
  struct Ray
  {
    Vector3 start;
    Vector3 omega;
    double tolerance = 0.0;
  };

  /**
   * Returns the value of t at which the ray leaves the cell, if the ray lies in the cell just
   * after t, and appends the local ids of the neighbors across the planes through which it
   * leaves to exit_neighbors.
   */
  std::optional<double> ExitAfter(std::uint32_t cell,
                                  const Ray& ray,
                                  double t,
                                  std::vector<std::uint32_t>* exit_neighbors) const;

  /// Interval of t in which the ray lies in a convex piece; empty when exit <= entry.
  void PieceInterval(std::uint32_t piece,
                     const Ray& ray,
                     double& entry,
                     double& exit,
                     std::vector<std::uint32_t>* exit_neighbors,
                     double exit_tolerance) const;

  /// Returns the cell that continues the ray at t, searching outward from cell.
  std::optional<std::uint32_t> FindNext(std::uint32_t cell,
                                        const Ray& ray,
                                        double t,
                                        const std::vector<std::uint32_t>& exit_neighbors,
                                        double& exit) const;

  static constexpr std::uint32_t no_neighbor = std::numeric_limits<std::uint32_t>::max();

  std::vector<bool> cell_convex_;
  std::vector<double> cell_size_;
  /// Pieces of each cell: [cell_piece_offset_[c], cell_piece_offset_[c + 1]).
  std::vector<std::uint32_t> cell_piece_offset_;
  /// Planes of each piece: [piece_plane_offset_[p], piece_plane_offset_[p + 1]).
  std::vector<std::uint32_t> piece_plane_offset_;
  /// Outward unit normal n and offset d of each plane: the piece lies where n . x <= d.
  std::vector<Vector3> plane_normal_;
  std::vector<double> plane_offset_;
  /// Local id of the cell across each plane, or no_neighbor for a boundary face, a face whose
  /// neighbor is not local, or a plane inside a cell.
  std::vector<std::uint32_t> plane_neighbor_;
  /// Local cells that share each vertex, and vertices of each cell, as compact vertex indices.
  std::vector<std::uint32_t> vertex_cell_offset_;
  std::vector<std::uint32_t> vertex_cells_;
  std::vector<std::uint32_t> cell_vertex_offset_;
  std::vector<std::uint32_t> cell_vertices_;
};

} // namespace opensn
