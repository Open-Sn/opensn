// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/raytrace/segment_tracer.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/runtime.h"
#include "test/unit/common/mesh_builders.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace opensn;

namespace
{

std::vector<std::uint32_t>
ContainingCells(const Mesh& grid, const Vector3& point)
{
  std::vector<std::uint32_t> cells;
  for (const auto& cell : grid.GetLocalCells())
    if (grid.CheckPointInsideCell(*cell, point))
      cells.push_back(cell->local_id);
  return cells;
}

/// Points of the generalized golden-ratio sequence in [0, 1)^6, spread over the cube.
std::array<double, 6>
SequencePoint(int k)
{
  // phi_6, the root of x^7 = x + 1.
  const double g = 1.11257658;
  std::array<double, 6> point{};
  double a = 1.0;
  for (auto& x : point)
  {
    a /= g;
    x = std::fmod(0.5 + k * a, 1.0);
  }
  return point;
}

/**
 * Checks a traced segment against an independent oracle. The segment is cut at its crossings
 * with the plane of every side triangle of every cell. Between two crossings it lies in one
 * cell, or on a face shared by two, so the cell that the tracer assigns to the midpoint of every
 * piece must contain that midpoint according to Mesh::CheckPointInsideCell. The traced lengths
 * must sum to the segment length.
 */
void
CheckSegment(const Mesh& grid, const SegmentTracer& tracer, const Vector3& a, const Vector3& b)
{
  const auto starts = ContainingCells(grid, a);
  ASSERT_FALSE(starts.empty()) << "start " << a.PrintStr() << " is not in the mesh";
  const double length = (b - a).Norm();
  const Vector3 omega = (b - a) / length;
  std::vector<double> cuts = {0.0, length};
  for (const auto& cell : grid.GetLocalCells())
    for (const auto& face : cell->faces)
    {
      const size_t num_sides = face.vertex_ids.size();
      for (size_t s = 0; s < num_sides; ++s)
      {
        const auto& v0 = grid.GlobalVertex(face.vertex_ids[s]);
        const auto& v1 = grid.GlobalVertex(face.vertex_ids[(s + 1) % num_sides]);
        const Vector3 normal = (v1 - v0).Cross(face.centroid - v0);
        const double cosine = normal.Dot(omega);
        if (std::abs(cosine) <= 1.0e-14 * normal.Norm())
          continue;
        const double t = normal.Dot(v0 - a) / cosine;
        if (t > 0.0 and t < length)
          cuts.push_back(t);
      }
    }
  std::sort(cuts.begin(), cuts.end());

  // Start in up to four of the cells that contain the start, so that rays that only touch the
  // start cell are tested.
  for (size_t k = 0; k < std::min<size_t>(starts.size(), 4); ++k)
  {
    std::vector<SegmentTracer::Segment> segments;
    ASSERT_TRUE(tracer.Trace(starts[k], a, b, segments))
      << "segment " << a.PrintStr() << " -> " << b.PrintStr() << " left the mesh";
    double total = 0.0;
    std::vector<double> ends;
    for (const auto& segment : segments)
    {
      EXPECT_GE(segment.length, 0.0);
      total += segment.length;
      ends.push_back(total);
    }
    EXPECT_NEAR(total, length, 1.0e-12 * length);
    for (size_t c = 0; c + 1 < cuts.size(); ++c)
    {
      if (cuts[c + 1] - cuts[c] <= 1.0e-7 * length)
        continue;
      const double t = 0.5 * (cuts[c] + cuts[c + 1]);
      const auto piece = std::upper_bound(ends.begin(), ends.end(), t) - ends.begin();
      ASSERT_LT(static_cast<size_t>(piece), segments.size());
      const auto& cell = grid.GetLocalCell(segments[piece].cell_local_id);
      EXPECT_TRUE(grid.CheckPointInsideCell(cell, a + t * omega))
        << "segment " << a.PrintStr() << " -> " << b.PrintStr() << " assigns t=" << t << " to cell "
        << cell.global_id << ", which does not contain it";
    }
  }
}

/// Segments between generic points, between mesh vertices (through vertices, along edges, and in
/// face planes), and from vertices to generic points, in the box [lo, hi]^3.
void
CheckMesh(const Mesh& grid, const SegmentTracer& tracer, double lo, double hi)
{
  std::vector<Vector3> vertices;
  for (const auto& cell : grid.GetLocalCells())
    for (const auto vid : cell->vertex_ids)
      vertices.push_back(grid.GlobalVertex(vid));
  const double span = hi - lo;
  for (int k = 1; k <= 60; ++k)
  {
    const auto p = SequencePoint(k);
    const Vector3 a(lo + span * p[0], lo + span * p[1], lo + span * p[2]);
    const Vector3 b(lo + span * p[3], lo + span * p[4], lo + span * p[5]);
    CheckSegment(grid, tracer, a, b);
    const auto& va = vertices[static_cast<size_t>(p[0] * static_cast<double>(vertices.size()))];
    const auto& vb = vertices[static_cast<size_t>(p[3] * static_cast<double>(vertices.size()))];
    if ((vb - va).Norm() > 1.0e-6 * span)
    {
      CheckSegment(grid, tracer, va, vb);
      CheckSegment(grid, tracer, va, b);
    }
    // A segment that lies in a coordinate plane through a vertex: on faces of an orthogonal
    // mesh, and through many edges and vertices.
    const Vector3 c(va.x, b.y, b.z);
    if ((c - va).Norm() > 1.0e-6 * span)
      CheckSegment(grid, tracer, va, c);
  }
}

std::filesystem::path
Asset(const std::string& name)
{
  return std::filesystem::path(OPENSN_TEST_ROOT) / "assets" / "mesh" / name;
}

/**
 * Writes the polygon (counterclockwise in the xy-plane) extruded from z = 0 to z = 1 as a single
 * VTK polyhedron to a temporary VTU file, reads it, and removes the file.
 */
std::shared_ptr<Mesh>
BuildExtrudedPolygonCell(const std::vector<std::array<double, 2>>& polygon)
{
  const size_t n = polygon.size();
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_segment_tracer_cell_" + std::to_string(now) + ".vtu");
  {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n<Piece NumberOfPoints=\"" << 2 * n << "\" NumberOfCells=\"1\">\n"
        << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const double z : {0.0, 1.0})
      for (const auto& [x, y] : polygon)
        out << x << " " << y << " " << z << "\n";
    out << "</DataArray>\n</Points>\n<Cells>\n"
        << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t i = 0; i < 2 * n; ++i)
      out << i << " ";
    out << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n"
        << 2 * n << "\n</DataArray>\n"
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n42\n</DataArray>\n";
    // Faces with their vertices counterclockwise seen from outside: the bottom, the top, and the
    // sides.
    std::vector<size_t> faces = {n + 2, n};
    for (size_t i = 0; i < n; ++i)
      faces.push_back(n - 1 - i);
    faces.push_back(n);
    for (size_t i = 0; i < n; ++i)
      faces.push_back(n + i);
    for (size_t i = 0; i < n; ++i)
      for (const size_t v : {size_t{4}, i, (i + 1) % n, n + (i + 1) % n, n + i})
        faces.push_back(v);
    out << "<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n";
    for (const auto v : faces)
      out << v << " ";
    out << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">\n"
        << faces.size() << "\n</DataArray>\n</Cells>\n<CellData>\n"
        << "<DataArray type=\"Int32\" Name=\"BlockID\" format=\"ascii\">\n0\n</DataArray>\n"
        << "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  }
  std::shared_ptr<Mesh> grid;
  try
  {
    grid = BuildMeshFromFile(path);
  }
  catch (...)
  {
    std::filesystem::remove(path);
    throw;
  }
  std::filesystem::remove(path);
  return grid;
}

double
TracedLength(const SegmentTracer& tracer, const Vector3& a, const Vector3& b, bool& reached_end)
{
  std::vector<SegmentTracer::Segment> segments;
  reached_end = tracer.Trace(0, a, b, segments);
  double total = 0.0;
  for (const auto& segment : segments)
  {
    EXPECT_EQ(segment.cell_local_id, 0);
    total += segment.length;
  }
  return total;
}

} // namespace

TEST(SegmentTracer, OrthogonalHexahedra)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  for (const bool split : {false, true})
    CheckMesh(*grid, SegmentTracer(*grid, split), 0.0, 1.0);
}

TEST(SegmentTracer, Tetrahedra)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildMeshFromFile(Asset("cube3.2.msh"));
  for (const bool split : {false, true})
    CheckMesh(*grid, SegmentTracer(*grid, split), 0.0, 0.032);
}

// segment_tracer_prisms.vtu: the unit square meshed by Gmsh 4.15.2 with 0.22 triangles and
// extruded to z = 1 in 4 layers (264 triangular prisms), written as VTK wedges.
TEST(SegmentTracer, Prisms)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildMeshFromFile(Asset("segment_tracer_prisms.vtu"));
  for (const bool split : {false, true})
    CheckMesh(*grid, SegmentTracer(*grid, split), 0.0, 1.0);
}

TEST(SegmentTracer, LeavesMesh)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  const SegmentTracer tracer(*grid);
  // From inside to beyond the face x = 1 through the edge x = 1, y = 0.5.
  const Vector3 a(0.3, 0.2, 0.6);
  const Vector3 b(1.4, 0.2 + 1.1 * 0.3 / 0.7, 0.6);
  std::vector<SegmentTracer::Segment> segments;
  EXPECT_FALSE(tracer.Trace(ContainingCells(*grid, a).front(), a, b, segments));
  double total = 0.0;
  for (const auto& segment : segments)
    total += segment.length;
  EXPECT_NEAR(total, 0.7 / 1.1 * (b - a).Norm(), 1.0e-12);
}

TEST(SegmentTracer, LocateAheadAtVertex)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  const SegmentTracer tracer(*grid);
  // From the vertex (0.5, 0.5, 0.5), shared by eight cells, into each octant.
  const Vector3 vertex(0.5, 0.5, 0.5);
  const auto hint = ContainingCells(*grid, vertex).front();
  for (const double sx : {-1.0, 1.0})
    for (const double sy : {-1.0, 1.0})
      for (const double sz : {-1.0, 1.0})
      {
        const Vector3 omega = Vector3(sx * 0.3, sy * 0.5, sz * 0.8).Normalized();
        const auto cell = tracer.LocateAhead(hint, vertex, omega);
        if (cell)
          EXPECT_TRUE(grid->CheckPointInsideCell(grid->GetLocalCell(*cell), vertex + 0.01 * omega));
        else
          ADD_FAILURE() << "no cell ahead of the vertex along " << omega.PrintStr();
      }
}

// A square with a V-shaped notch in its top edge, extruded in z. The cell is not convex, but its
// centroid (1, 9.6 / 7, 0.5) lies on the inner side of every side triangle, so its side
// tetrahedra fill it. A segment at y = 1.9 leaves the cell where it meets the notch, at
// x = 0.8 + 0.5 (2 - 1.9) = 0.85.
TEST(SegmentTracer, NonConvexCell)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildExtrudedPolygonCell(
    {{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {1.2, 2.0}, {1.0, 1.6}, {0.8, 2.0}, {0.0, 2.0}});
  const SegmentTracer tracer(*grid);
  EXPECT_FALSE(tracer.IsConvexCell(0));
  bool reached_end = false;
  EXPECT_NEAR(
    TracedLength(tracer, Vector3(0.2, 1.5, 0.5), Vector3(1.8, 1.5, 0.5), reached_end), 1.6, 1e-12);
  EXPECT_TRUE(reached_end);
  EXPECT_NEAR(
    TracedLength(tracer, Vector3(0.2, 1.9, 0.5), Vector3(1.8, 1.9, 0.5), reached_end), 0.65, 1e-12);
  EXPECT_FALSE(reached_end);
}

// A U-shaped cell, extruded in z. Its centroid (1.5, 1.25, 0.5) lies in the notch, outside the
// cell, so some of its side tetrahedra lie outside it. The tracer rejects the cell.
TEST(SegmentTracer, RejectsCellNotStarShapedAboutCentroid)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildExtrudedPolygonCell({{0.0, 0.0},
                                              {3.0, 0.0},
                                              {3.0, 2.0},
                                              {2.0, 2.0},
                                              {2.0, 1.0},
                                              {1.0, 1.0},
                                              {1.0, 2.0},
                                              {0.0, 2.0}});
  EXPECT_THROW(SegmentTracer tracer(*grid), std::invalid_argument);
}
