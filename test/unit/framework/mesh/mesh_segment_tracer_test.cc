// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/raytrace/mesh_segment_tracer.h"
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
#include <limits>
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
CheckSegment(const Mesh& grid, const MeshSegmentTracer& tracer, const Vector3& a, const Vector3& b)
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
    std::vector<MeshSegmentTracer::Segment> segments;
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
CheckMesh(const Mesh& grid, const MeshSegmentTracer& tracer, double lo, double hi)
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
 * Writes the cells, each a list of faces given by point indices counterclockwise seen from
 * outside, as VTK polyhedra to a temporary VTU file, reads it, and removes the file.
 */
std::shared_ptr<Mesh>
BuildPolyhedraMesh(const std::vector<Vector3>& points,
                   const std::vector<std::vector<std::vector<size_t>>>& cells)
{
  const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
  const auto path = std::filesystem::temp_directory_path() /
                    ("opensn_mesh_segment_tracer_cell_" + std::to_string(now) + ".vtu");
  {
    std::ofstream out(path);
    out << std::setprecision(17);
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n<Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\""
        << cells.size() << "\">\n"
        << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& p : points)
      out << p.x << " " << p.y << " " << p.z << "\n";
    out << "</DataArray>\n</Points>\n<Cells>\n"
        << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
    std::vector<size_t> offsets;
    std::vector<size_t> faces;
    std::vector<size_t> face_offsets;
    size_t num_cell_points = 0;
    for (const auto& cell : cells)
    {
      std::vector<size_t> cell_points;
      for (const auto& face : cell)
        for (const auto v : face)
          if (std::find(cell_points.begin(), cell_points.end(), v) == cell_points.end())
            cell_points.push_back(v);
      std::sort(cell_points.begin(), cell_points.end());
      for (const auto v : cell_points)
        out << v << " ";
      num_cell_points += cell_points.size();
      offsets.push_back(num_cell_points);
      faces.push_back(cell.size());
      for (const auto& face : cell)
      {
        faces.push_back(face.size());
        faces.insert(faces.end(), face.begin(), face.end());
      }
      face_offsets.push_back(faces.size());
    }
    out << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
    for (const auto o : offsets)
      out << o << " ";
    out << "\n</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (size_t c = 0; c < cells.size(); ++c)
      out << "42 ";
    out << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"faces\" format=\"ascii\">\n";
    for (const auto v : faces)
      out << v << " ";
    out << "\n</DataArray>\n<DataArray type=\"Int64\" Name=\"faceoffsets\" format=\"ascii\">\n";
    for (const auto o : face_offsets)
      out << o << " ";
    out << "\n</DataArray>\n</Cells>\n<CellData>\n"
        << "<DataArray type=\"Int32\" Name=\"BlockID\" format=\"ascii\">\n";
    for (size_t c = 0; c < cells.size(); ++c)
      out << "0 ";
    out << "\n</DataArray>\n</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
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

/// The polygon (counterclockwise in the xy-plane) extruded from z = 0 to z = 1 as a single cell.
std::shared_ptr<Mesh>
BuildExtrudedPolygonCell(const std::vector<std::array<double, 2>>& polygon)
{
  const size_t n = polygon.size();
  std::vector<Vector3> points;
  for (const double z : {0.0, 1.0})
    for (const auto& [x, y] : polygon)
      points.emplace_back(x, y, z);
  // The bottom, the top, and the sides.
  std::vector<std::vector<size_t>> faces(2);
  for (size_t i = 0; i < n; ++i)
  {
    faces[0].push_back(n - 1 - i);
    faces[1].push_back(n + i);
  }
  for (size_t i = 0; i < n; ++i)
    faces.push_back({i, (i + 1) % n, n + (i + 1) % n, n + i});
  return BuildPolyhedraMesh(points, {faces});
}

/**
 * The unit cube with its top face replaced by a pyramidal dent whose apex is the center of the
 * cube. The average of the nine vertices is the apex, so the cell centroid lies in the plane of
 * each of the four dent triangles. The cell is the part of the cube below the cone
 * z - 1/2 = max(|x - 1/2|, |y - 1/2|) and is star-shaped about the apex.
 */
std::shared_ptr<Mesh>
BuildDentedCube()
{
  const std::vector<Vector3> points = {{0.0, 0.0, 0.0},
                                       {1.0, 0.0, 0.0},
                                       {1.0, 1.0, 0.0},
                                       {0.0, 1.0, 0.0},
                                       {0.0, 0.0, 1.0},
                                       {1.0, 0.0, 1.0},
                                       {1.0, 1.0, 1.0},
                                       {0.0, 1.0, 1.0},
                                       {0.5, 0.5, 0.5}};
  const std::vector<std::vector<size_t>> faces = {{0, 3, 2, 1},
                                                  {0, 1, 5, 4},
                                                  {1, 2, 6, 5},
                                                  {2, 3, 7, 6},
                                                  {3, 0, 4, 7},
                                                  {4, 5, 8},
                                                  {5, 6, 8},
                                                  {6, 7, 8},
                                                  {7, 4, 8}};
  return BuildPolyhedraMesh(points, {faces});
}

double
TracedLength(const MeshSegmentTracer& tracer, const Vector3& a, const Vector3& b, bool& reached_end)
{
  std::vector<MeshSegmentTracer::Segment> segments;
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

TEST(MeshSegmentTracer, OrthogonalHexahedra)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  for (const bool split : {false, true})
    CheckMesh(*grid, MeshSegmentTracer(*grid, split), 0.0, 1.0);
}

TEST(MeshSegmentTracer, Tetrahedra)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildMeshFromFile(Asset("cube3.2.msh"));
  for (const bool split : {false, true})
    CheckMesh(*grid, MeshSegmentTracer(*grid, split), 0.0, 0.032);
}

// mesh_segment_tracer_prisms.vtu: the unit square meshed by Gmsh 4.15.2 with 0.22 triangles and
// extruded to z = 1 in 4 layers (264 triangular prisms), written as VTK wedges.
TEST(MeshSegmentTracer, Prisms)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildMeshFromFile(Asset("mesh_segment_tracer_prisms.vtu"));
  for (const bool split : {false, true})
    CheckMesh(*grid, MeshSegmentTracer(*grid, split), 0.0, 1.0);
}

// A segment that crosses an interior face within 1e-7 radians of its plane. The two cells
// compute the crossing at values of t that differ by their round-off divided by the cosine, which
// exceeds the tolerance, so the crossing must be decided from the distances to the plane. The
// segment is symmetric about the face centroid, so each cell receives half of it.
TEST(MeshSegmentTracer, GrazingFaces)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildMeshFromFile(Asset("cube3.2.msh"));
  const MeshSegmentTracer tracer(*grid);
  constexpr double grazing_slope = 1.0e-7;
  size_t count = 0;
  for (const auto& cell : grid->GetLocalCells())
    for (const auto& face : cell->faces)
    {
      if (not face.has_neighbor or face.neighbor_id < cell->global_id)
        continue;
      const auto& f = face.centroid;
      const size_t num_sides = face.vertex_ids.size();
      double edge_distance = std::numeric_limits<double>::infinity();
      for (size_t s = 0; s < num_sides; ++s)
      {
        const auto& v0 = grid->GlobalVertex(face.vertex_ids[s]);
        const auto& v1 = grid->GlobalVertex(face.vertex_ids[(s + 1) % num_sides]);
        edge_distance = std::min(edge_distance, (v1 - v0).Cross(f - v0).Norm() / (v1 - v0).Norm());
      }
      const double half_length = 0.3 * edge_distance;
      const Vector3 u = (grid->GlobalVertex(face.vertex_ids[0]) - f).Normalized();
      const Vector3 step = half_length * (u + grazing_slope * face.normal);
      const Vector3 a = f - step;
      const Vector3 b = f + step;
      // The oracle costs a pass over the mesh per segment, so it checks a sample of the faces.
      if (count % 32 == 0)
        CheckSegment(*grid, tracer, a, b);

      std::vector<MeshSegmentTracer::Segment> segments;
      ASSERT_TRUE(tracer.Trace(cell->local_id, a, b, segments))
        << "segment " << a.PrintStr() << " -> " << b.PrintStr() << " left the mesh";
      const auto neighbor = grid->GetGlobalCell(face.neighbor_id).local_id;
      double in_cell = 0.0;
      double in_neighbor = 0.0;
      for (const auto& segment : segments)
      {
        if (segment.cell_local_id == cell->local_id)
          in_cell += segment.length;
        else if (segment.cell_local_id == neighbor)
          in_neighbor += segment.length;
        else
          ADD_FAILURE() << "segment " << a.PrintStr() << " -> " << b.PrintStr()
                        << " assigns length to cell "
                        << grid->GetLocalCell(segment.cell_local_id).global_id;
      }
      EXPECT_NEAR(in_cell, step.Norm(), 1.0e-6 * step.Norm());
      EXPECT_NEAR(in_neighbor, step.Norm(), 1.0e-6 * step.Norm());
      ++count;
    }
  EXPECT_GT(count, 100);
}

// Every point of a cell traced to a point source, as the uncollided problem does: from the cell
// centroid, from points toward every vertex, and from points near every face, to a generic
// source, to a vertex, to a point on an edge, and to a point on a face. No ray may leave the
// convex mesh, the lengths must sum to the distance, and the ray must end in a cell that
// contains the source.
TEST(MeshSegmentTracer, AllRaysToSource)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildMeshFromFile(Asset("cube3.2.msh"));
  const Vector3 generic_source(0.0102586, 0.0114131, 0.0146416);
  const auto source_cell = ContainingCells(*grid, generic_source);
  ASSERT_EQ(source_cell.size(), 1U);
  const auto& cell_of_source = grid->GetLocalCell(source_cell.front());
  // An interior face of the source cell, with a vertex and an edge of it.
  const CellFace* interior_face = nullptr;
  for (const auto& face : cell_of_source.faces)
    if (face.has_neighbor)
      interior_face = &face;
  ASSERT_NE(interior_face, nullptr);
  const auto& v0 = grid->GlobalVertex(interior_face->vertex_ids[0]);
  const auto& v1 = grid->GlobalVertex(interior_face->vertex_ids[1]);
  const std::vector<Vector3> sources = {
    generic_source, v0, 0.5 * (v0 + v1), interior_face->centroid};

  for (const bool split : {false, true})
  {
    const MeshSegmentTracer tracer(*grid, split);
    size_t traced = 0;
    size_t lost = 0;
    size_t wrong_end = 0;
    double worst_length_error = 0.0;
    std::vector<MeshSegmentTracer::Segment> segments;
    for (const auto& cell : grid->GetLocalCells())
    {
      std::vector<Vector3> starts = {cell->centroid};
      for (const auto vid : cell->vertex_ids)
        starts.push_back(0.5 * (cell->centroid + grid->GlobalVertex(vid)));
      for (const auto& face : cell->faces)
        starts.push_back(cell->centroid + 0.999 * (face.centroid - cell->centroid));
      for (const auto& start : starts)
        for (const auto& source : sources)
        {
          const double distance = (source - start).Norm();
          if (distance == 0.0)
            continue;
          segments.clear();
          ++traced;
          if (not tracer.Trace(cell->local_id, start, source, segments))
          {
            ++lost;
            continue;
          }
          double total = 0.0;
          for (const auto& segment : segments)
          {
            EXPECT_GE(segment.length, 0.0);
            total += segment.length;
          }
          worst_length_error = std::max(worst_length_error, std::abs(total - distance) / distance);
          if (not grid->CheckPointInsideCell(grid->GetLocalCell(segments.back().cell_local_id),
                                             source))
            ++wrong_end;
        }
    }
    EXPECT_GT(traced, 100000U);
    EXPECT_EQ(lost, 0U) << "split=" << split;
    EXPECT_EQ(wrong_end, 0U) << "split=" << split;
    EXPECT_LT(worst_length_error, 1.0e-12) << "split=" << split;
  }
}

// A hexahedron whose top face is a strip 1e-4 wide and 1 long. The side triangles of that face
// are thinner than 1e-3 of the cell, so the tracer counts them. Ordinary meshes have none.
TEST(MeshSegmentTracer, ThinSideTriangles)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const double w = 1.0e-4;
  const std::vector<Vector3> points = {{0, 0, 0},
                                       {1, 0, 0},
                                       {1, 1, 0},
                                       {0, 1, 0},
                                       {0, 0.5 - 0.5 * w, 1},
                                       {1, 0.5 - 0.5 * w, 1},
                                       {1, 0.5 + 0.5 * w, 1},
                                       {0, 0.5 + 0.5 * w, 1}};
  const std::vector<std::vector<size_t>> hexahedron = {
    {0, 3, 2, 1}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}, {4, 5, 6, 7}};
  const auto grid = BuildPolyhedraMesh(points, {hexahedron});
  const MeshSegmentTracer tracer(*grid);
  // The top face: its four side triangles have altitudes of about w / 2 or less, and one of them
  // is kept as the plane of the face.
  EXPECT_EQ(tracer.NumThinSideTriangles(), 1U);
  bool reached_end = false;
  const Vector3 a(0.2, 0.3, 0.1);
  const Vector3 b(0.7, 0.52, 0.9);
  EXPECT_NEAR(TracedLength(tracer, a, b, reached_end), (b - a).Norm(), 1.0e-12);
  EXPECT_TRUE(reached_end);

  EXPECT_EQ(MeshSegmentTracer(*BuildBoxMesh(1.0, 2, 0.0)).NumThinSideTriangles(), 0U);
  EXPECT_EQ(MeshSegmentTracer(*BuildMeshFromFile(Asset("cube3.2.msh"))).NumThinSideTriangles(), 0U);
}

TEST(MeshSegmentTracer, LeavesMesh)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  const MeshSegmentTracer tracer(*grid);
  // From inside to beyond the face x = 1 through the edge x = 1, y = 0.5.
  const Vector3 a(0.3, 0.2, 0.6);
  const Vector3 b(1.4, 0.2 + 1.1 * 0.3 / 0.7, 0.6);
  std::vector<MeshSegmentTracer::Segment> segments;
  EXPECT_FALSE(tracer.Trace(ContainingCells(*grid, a).front(), a, b, segments));
  double total = 0.0;
  for (const auto& segment : segments)
    total += segment.length;
  EXPECT_NEAR(total, 0.7 / 1.1 * (b - a).Norm(), 1.0e-12);
}

TEST(MeshSegmentTracer, LocateAheadAtVertex)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  const MeshSegmentTracer tracer(*grid);
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
TEST(MeshSegmentTracer, NonConvexCell)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildExtrudedPolygonCell(
    {{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {1.2, 2.0}, {1.0, 1.6}, {0.8, 2.0}, {0.0, 2.0}});
  const MeshSegmentTracer tracer(*grid);
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
TEST(MeshSegmentTracer, RejectsCellNotStarShapedAboutCentroid)
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
  EXPECT_THROW(MeshSegmentTracer tracer(*grid), std::invalid_argument);
}

// The dented cube's centroid lies in the planes of the four dent triangles, so the side
// tetrahedra of those faces are flat and have no volume. The tracer leaves them out; the side
// tetrahedra of the other faces fill the cell. Were they kept, their four planes would coincide,
// and a segment in the plane of a dent triangle would lie in them without an exit. The plane of
// the dent triangle on the edge y = 0, z = 1 is z = 1 - y.
TEST(MeshSegmentTracer, FlatSideTetrahedra)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildDentedCube();
  ASSERT_EQ(grid->GetLocalCells().size(), 1);
  const auto& centroid = grid->GetLocalCell(0).centroid;
  EXPECT_NEAR((centroid - Vector3(0.5, 0.5, 0.5)).Norm(), 0.0, 1.0e-15);
  const MeshSegmentTracer tracer(*grid);
  EXPECT_FALSE(tracer.IsConvexCell(0));
  bool reached_end = false;
  // In the plane z = 1 - y, strictly inside the cell until it leaves through the face x = 1 at
  // half its length.
  const Vector3 a(0.7, 0.4, 0.6);
  const Vector3 b(1.3, 0.1, 0.9);
  EXPECT_NEAR(TracedLength(tracer, a, b, reached_end), 0.5 * (b - a).Norm(), 1.0e-12);
  EXPECT_FALSE(reached_end);
  // In the same plane, lying on the dent triangle for part of its length, then leaving through
  // the face x = 1 at three quarters of its length.
  const Vector3 c(0.1, 0.4, 0.6);
  const Vector3 d(1.3, 0.1, 0.9);
  EXPECT_NEAR(TracedLength(tracer, c, d, reached_end), 0.75 * (d - c).Norm(), 1.0e-12);
  EXPECT_FALSE(reached_end);
  // Through the apex, which is the common vertex of every side tetrahedron.
  const Vector3 e(0.2, 0.2, 0.35);
  const Vector3 f(0.8, 0.8, 0.65);
  EXPECT_NEAR(TracedLength(tracer, e, f, reached_end), (f - e).Norm(), 1.0e-12);
  EXPECT_TRUE(reached_end);
}

// A mesh far from the origin: coordinates of 1e6 with cells of size 1. The tolerance stays far
// below the cell size, so a short segment that ends 2.5e-5 before the face x = 1 is traced in
// its cell, and the general segments are traced as near the origin.
TEST(MeshSegmentTracer, TranslatedMesh)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const double z0 = 1.0e6;
  const auto grid = BuildOrthogonalMesh({{0.0, 1.0}, {0.0, 1.0}, {z0, z0 + 1.0}});
  const MeshSegmentTracer tracer(*grid);
  const Vector3 a(0.99995, 0.5, z0 + 0.5);
  const Vector3 b(0.999975, 0.5, z0 + 0.5);
  std::vector<MeshSegmentTracer::Segment> segments;
  EXPECT_TRUE(tracer.Trace(0, a, b, segments));
  ASSERT_EQ(segments.size(), 1);
  EXPECT_EQ(segments[0].cell_local_id, 0);
  EXPECT_NEAR(segments[0].length, 2.5e-5, 1.0e-12);
  EXPECT_EQ(tracer.LocateAhead(0, a, Vector3(1.0, 0.0, 0.0)), 0);

  const auto box = BuildBoxMesh(1.0, 2, z0);
  for (const bool split : {false, true})
    CheckMesh(*box, MeshSegmentTracer(*box, split), z0, z0 + 1.0);
}
