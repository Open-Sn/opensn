// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/vector3.h"
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace opensn
{

class Mesh;

/**
 * Traversal of straight segments through the local cells of a mesh, decided by exact geometric
 * predicates.
 *
 * Every cell is traced as a set of convex pieces whose boundaries are triangles: a triangular
 * face is one triangle, and a face with more vertices is the fan of triangles from its sides to
 * its centroid. The centroid of a face is computed from the vertices in a canonical order, so
 * that the two cells sharing a face hold bit-identical triangles. A cell that is convex, with
 * every vertex of it on the inner side of every triangle, is one piece. Any other cell is split
 * into its side tetrahedra, each formed by a boundary triangle and an interior point of the
 * cell, which requires the cell to be star-shaped about that point; the tracer rejects a cell
 * that is not.
 *
 * A segment leaves a convex piece through the one boundary triangle that its line threads
 * outward, decided from the signs of the orientation predicate of the segment against the three
 * edges of the triangle (Platis and Theoharis 2003, Lagae and Dutre 2008). The predicates are
 * evaluated exactly (Shewchuk 1997), and the two pieces on either side of a triangle test the
 * same edges from the same points, so the traversal never finds a gap or an overlap between
 * pieces: it cannot lose a segment in the interior of a conforming mesh. Whether the segment
 * ends in a piece is decided exactly as well. Only the distances along the segment at which it
 * crosses the triangles are computed in floating point, and they are kept monotone, so the
 * lengths of a traversal are nonnegative and sum to the segment length.
 *
 * A segment whose line passes exactly through an edge or a vertex, or lies in the plane of a
 * face, makes a predicate zero. The tracer then traces the segment translated by a few units in
 * the last place of the coordinates, in several directions, which keeps the lengths within that
 * round-off, and reports the segment as reaching its end if the segment translated in any
 * direction does.
 *
 * The tracer works on the local cells of a mesh and treats a face whose neighbor is not local as
 * a boundary face.
 */
class MeshSegmentTracer
{
public:
  /// A part of a segment in one local cell.
  struct Segment
  {
    std::uint32_t cell_local_id;
    double length;
  };

  /**
   * Builds the pieces of every local cell and the vertex-to-cell adjacency. O(cells). A cell that
   * is not convex is traced through the tetrahedra of its triangles and an interior point about
   * which it is star-shaped, searched for from the average of its points. Throws
   * std::invalid_argument for a cell that is not star-shaped about any point.
   *
   * With split_all_cells, every cell is traced through its side tetrahedra, as a non-convex cell
   * is. The result is the same up to round-off; the option exists to test that path on meshes of
   * convex cells.
   */
  explicit MeshSegmentTracer(const Mesh& grid, bool split_all_cells = false);

  /**
   * Traces the segment from start to end, beginning in the local cell start_cell, which
   * contains start, has it on its boundary, or is next to the cell that does. A start that lies
   * outside start_cell by the round-off of its coordinates, as a point computed to lie on one of
   * its faces can, is moved into it. Appends the cells that the segment crosses, in order, with
   * the length of the segment in each, to segments.
   *
   * Returns true when the segment reaches end, which includes an end that lies outside a
   * boundary face by the round-off of its coordinates. Returns false when the segment leaves the
   * local cells first, through a boundary face, a face whose neighbor is not local, or a gap in
   * a non-conforming mesh; segments then holds the part traced before the segment left. Throws
   * if start does not lie in, on, or near start_cell.
   */
  bool Trace(std::uint32_t start_cell,
             const Vector3& start,
             const Vector3& end,
             std::vector<Segment>& segments) const;

  /// Whether the cell is traced as one convex piece rather than through its side tetrahedra.
  bool IsConvexCell(std::uint32_t cell_local_id) const { return cell_convex_[cell_local_id]; }

private:
  /// A boundary triangle of a piece, counterclockwise seen from outside the piece.
  struct Triangle
  {
    std::uint32_t a;
    std::uint32_t b;
    std::uint32_t c;
    /// The piece across the triangle, or boundary.
    std::uint32_t across;
    /// The edges (a, b), (b, c), (c, a) as indices into the edges of the piece, and whether each
    /// runs from the smaller to the larger point index, as the piece edge is stored.
    std::array<std::uint32_t, 3> edge;
    std::array<bool, 3> forward;
  };

  /// Thrown inside a traversal when a predicate is zero.
  struct Degenerate
  {
  };

  /// The outcome of one traversal of a (possibly translated) segment.
  enum class Outcome
  {
    reached,
    left,
    /// Left through a boundary triangle whose plane holds the end within round-off, so that
    /// displaced ends may keep the line inside.
    left_near_end,
    no_start
  };

  /**
   * The boundary triangle that the line from p to q threads outward (the exit of the line from
   * the piece) or inward (its entry), by the exact signs of the orientation predicate of the line
   * against its three edges, or nullptr when the line misses the piece. The triangles are tried in
   * the order of the floating-point crossings of their planes, so that the answer is usually the
   * first tried; the order decides nothing. Throws Degenerate when the line meets the line of an
   * edge that is tested.
   */
  const Triangle*
  Threaded(std::uint32_t piece, const Vector3& p, const Vector3& q, bool outward) const;

  /// The exit of the line from p to q from the piece, which must exist.
  const Triangle& ExitTriangle(std::uint32_t piece, const Vector3& p, const Vector3& q) const;

  /**
   * The round-off allowed to the ends of a segment that a caller computed, such as the crossing
   * of a ray with a reflecting plane, along each axis: 64 units in the last place of the largest
   * coordinate of the mesh and of the two ends along that axis. The allowances differ between axes
   * when the mesh is far from the origin along one of them and thin along another. The extent of
   * the mesh enters because the caller computes such points from coordinates of the mesh: a mesh
   * whose extent along an axis exceeds the size of its cells by more than about 1e13 cannot be
   * traced to the cell size along that axis.
   */
  Vector3 Roundoff(const Vector3& p, const Vector3& q) const;

  /// Whether the point lies in or on the piece (exact).
  bool Contains(std::uint32_t piece, const Vector3& point) const;

  /**
   * The local cells that contain the point, in or on: among the cell near_cell and the cells
   * within two layers of it, or, without near_cell, among all local cells. O(cells) in the second
   * form.
   */
  std::vector<std::uint32_t> CellsContaining(const Vector3& point,
                                             std::optional<std::uint32_t> near_cell) const;

  /// The round-off of a point with the allowances along a direction (a unit vector).
  static double Along(const Vector3& allowance, const Vector3& direction);

  /// Whether the point lies within the given allowance of the plane of the triangle.
  bool
  WithinRoundoff(const Triangle& triangle, const Vector3& point, const Vector3& allowance) const;

  /**
   * The piece in which the segment from p to q starts: a piece that contains p and that the
   * segment travels in, found from the pieces of the cell by following the triangles through
   * which the line leaves the pieces that have p on their boundary. With search_neighbors, the
   * cells within two layers of the cell are tried after it.
   */
  std::optional<std::uint32_t> FindStartPiece(std::uint32_t cell,
                                              const Vector3& p,
                                              const Vector3& q,
                                              bool search_neighbors) const;

  /// One traversal, which may throw Degenerate.
  Outcome Walk(std::uint32_t start_cell,
               const Vector3& p,
               const Vector3& q,
               std::vector<Segment>& segments,
               bool search_neighbors) const;

  /**
   * The traversal of the segment, with the displacements of a degenerate segment. Returns
   * whether the segment reaches its end, or std::nullopt when its start does not lie in, on, or
   * near start_cell.
   */
  std::optional<bool> Resolve(std::uint32_t start_cell,
                              const Vector3& start,
                              const Vector3& end,
                              std::vector<Segment>& segments) const;

  /**
   * Displacement of an end of a degenerate segment: a few units in the last place of the
   * coordinates along one of the eight directions with components along every axis, one per
   * octant, or toward the given interior point of the start cell, which points into the mesh at a
   * boundary point however thin the cell is there. The far end uses other directions and a larger
   * step, so that the displaced line is never parallel to the line of the segment and leaves
   * every degeneracy of it.
   */
  Vector3 Displacement(const Vector3& point,
                       const Vector3& other,
                       int direction,
                       bool far_end,
                       const Vector3& interior) const;

  static constexpr std::uint32_t boundary = std::numeric_limits<std::uint32_t>::max();
  /// Directions of displacement tried for each end of a degenerate segment: the eight octants,
  /// the direction toward the interior of the start cell, and the direction along the segment.
  static constexpr int num_directions = 10;

  /**
   * Displacements of an end of a segment into a cell that contains it: toward the interior point
   * of the cell and toward the midpoints between that point and the vertices of the cell, which
   * lie inside a convex piece and break the symmetry that keeps a displaced line on a face plane.
   */
  std::vector<Vector3> InwardDisplacements(const Vector3& point,
                                           const Vector3& other,
                                           bool far_end,
                                           std::uint32_t cell) const;

  /// Mesh vertices, face centroids, and cell centroids used by the pieces.
  std::vector<Vector3> points_;
  /// Largest magnitude of a coordinate of the mesh along each axis.
  Vector3 scale_;
  /// The plane of the first triangle of a face of a piece, for ordering the exit candidates.
  struct Plane
  {
    Vector3 normal;
    double offset = 0.0;
  };

  /// Triangles of each piece: [piece_triangle_offset_[p], piece_triangle_offset_[p + 1]).
  std::vector<Triangle> triangles_;
  std::vector<std::uint32_t> piece_triangle_offset_;
  /// The faces of each piece, [piece_face_offset_[p], piece_face_offset_[p + 1]), each the
  /// triangles from face_first_triangle_[f] to the next face's first triangle or the end of the
  /// piece, with the plane of its first triangle.
  std::vector<std::uint32_t> face_first_triangle_;
  std::vector<Plane> face_planes_;
  std::vector<std::uint32_t> piece_face_offset_;
  /// Distinct edges of each piece, from the smaller to the larger point index:
  /// [piece_edge_offset_[p], piece_edge_offset_[p + 1]). A line is tested against each edge once.
  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges_;
  std::vector<std::uint32_t> piece_edge_offset_;
  /// Cell of each piece, and pieces of each cell: [cell_piece_offset_[c], cell_piece_offset_[c +
  /// 1]).
  std::vector<std::uint32_t> piece_cell_;
  std::vector<std::uint32_t> cell_piece_offset_;
  std::vector<bool> cell_convex_;
  /// Diagonal of the bounding box of each cell.
  std::vector<double> cell_size_;
  /// The interior point of each cell about which it is star-shaped.
  std::vector<Vector3> cell_interior_;
  /// Local cells that share each vertex, and vertices of each cell, as compact vertex indices.
  std::vector<std::uint32_t> cell_vertices_;
  std::vector<std::uint32_t> cell_vertex_offset_;
  std::vector<std::uint32_t> vertex_cells_;
  std::vector<std::uint32_t> vertex_cell_offset_;
};

} // namespace opensn
