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
#include <map>
#include <memory>
#include <optional>
#include <sstream>
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
  for (std::uint32_t c = 0; c < grid.GetLocalCellCount(); ++c)
    if (grid.CheckPointInsideCell(c, point))
      cells.push_back(c);
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
  for (std::uint64_t cell_local_id = 0; cell_local_id < grid.GetLocalCellCount(); ++cell_local_id)
  {
    const auto& cell = grid.GetLocalCell(cell_local_id);
    auto cell_faces = grid.GetCellFaces(cell_local_id);
    for (std::uint32_t fi = 0; fi < cell_faces.size(); ++fi)
    {
      const auto& face = cell_faces[fi];
      auto face_vertex_ids = grid.GetCellFaceConnectivity(cell_local_id, fi);
      const size_t num_sides = face_vertex_ids.size();
      for (size_t s = 0; s < num_sides; ++s)
      {
        const auto& v0 = grid.GlobalVertex(face_vertex_ids[s]);
        const auto& v1 = grid.GlobalVertex(face_vertex_ids[(s + 1) % num_sides]);
        const Vector3 normal = (v1 - v0).Cross(face.centroid - v0);
        const double cosine = normal.Dot(omega);
        if (std::abs(cosine) <= 1.0e-14 * normal.Norm())
          continue;
        const double t = normal.Dot(v0 - a) / cosine;
        if (t > 0.0 and t < length)
          cuts.push_back(t);
      }
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
      EXPECT_TRUE(grid.CheckPointInsideCell(segments[piece].cell_local_id, a + t * omega))
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
  for (std::uint32_t cell_local_id = 0; cell_local_id < grid.GetLocalCellCount(); ++cell_local_id)
  {
    auto cell_vertex_ids = grid.GetCellConnectivity(cell_local_id);
    for (const auto vid : cell_vertex_ids)
      vertices.push_back(grid.GlobalVertex(vid));
  }
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
  for (std::uint32_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount(); ++cell_local_id)
  {
    const auto& cell = grid->GetLocalCell(cell_local_id);
    auto cell_faces = grid->GetCellFaces(cell_local_id);
    for (std::uint32_t fi = 0; fi < cell_faces.size(); ++fi)
    {
      const auto& face = cell_faces[fi];
      if (not face.has_neighbor or face.neighbor_id < cell.global_id)
        continue;
      const auto& f = face.centroid;
      auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, fi);
      const size_t num_sides = face_vertex_ids.size();
      double edge_distance = std::numeric_limits<double>::infinity();
      for (size_t s = 0; s < num_sides; ++s)
      {
        const auto& v0 = grid->GlobalVertex(face_vertex_ids[s]);
        const auto& v1 = grid->GlobalVertex(face_vertex_ids[(s + 1) % num_sides]);
        edge_distance = std::min(edge_distance, (v1 - v0).Cross(f - v0).Norm() / (v1 - v0).Norm());
      }
      const double half_length = 0.3 * edge_distance;
      const Vector3 u = (grid->GlobalVertex(face_vertex_ids[0]) - f).Normalized();
      const Vector3 step = half_length * (u + grazing_slope * face.normal);
      const Vector3 a = f - step;
      const Vector3 b = f + step;
      // The oracle costs a pass over the mesh per segment, so it checks a sample of the faces.
      if (count % 32 == 0)
        CheckSegment(*grid, tracer, a, b);

      std::vector<MeshSegmentTracer::Segment> segments;
      ASSERT_TRUE(tracer.Trace(cell_local_id, a, b, segments))
        << "segment " << a.PrintStr() << " -> " << b.PrintStr() << " left the mesh";
      const auto neighbor = grid->MapCellGlobalID2LocalID(face.neighbor_id);
      double in_cell = 0.0;
      double in_neighbor = 0.0;
      for (const auto& segment : segments)
      {
        if (segment.cell_local_id == cell_local_id)
          in_cell += segment.length;
        else if (segment.cell_local_id == neighbor)
          in_neighbor += segment.length;
        else
          ADD_FAILURE() << "segment " << a.PrintStr() << " -> " << b.PrintStr()
                        << " assigns length to cell "
                        << grid->GetLocalCell(segment.cell_local_id).global_id;
      }
      // The segment passes exactly through the face centroid, a point of every triangle of the
      // face, so the tracer traces it displaced by the round-off allowance; at the grazing slope
      // the crossing moves by that displacement over the slope, about 1e-5 of the half length.
      EXPECT_NEAR(in_cell, step.Norm(), 1.0e-4 * step.Norm());
      EXPECT_NEAR(in_neighbor, step.Norm(), 1.0e-4 * step.Norm());
      ++count;
    }
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
  std::span<const std::uint64_t> interior_face_vertex_ids;
  auto cell_of_source_faces = grid->GetCellFaces(source_cell.front());
  for (std::uint32_t face_idx = 0; face_idx < cell_of_source_faces.size(); ++face_idx)
  {
    const auto& face = cell_of_source_faces[face_idx];
    if (face.has_neighbor)
    {
      interior_face = &face;
      interior_face_vertex_ids = grid->GetCellFaceConnectivity(source_cell.front(), face_idx);
    }
  }
  ASSERT_NE(interior_face, nullptr);
  const auto& v0 = grid->GlobalVertex(interior_face_vertex_ids[0]);
  const auto& v1 = grid->GlobalVertex(interior_face_vertex_ids[1]);
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
    for (std::uint32_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount();
         ++cell_local_id)
    {
      const auto& cell = grid->GetLocalCell(cell_local_id);
      auto cell_vertex_ids = grid->GetCellConnectivity(cell_local_id);
      std::vector<Vector3> starts = {cell.centroid};
      for (const auto vid : cell_vertex_ids)
        starts.push_back(0.5 * (cell.centroid + grid->GlobalVertex(vid)));
      auto cell_faces = grid->GetCellFaces(cell_local_id);
      for (const auto& face : cell_faces)
        starts.push_back(cell.centroid + 0.999 * (face.centroid - cell.centroid));
      for (const auto& start : starts)
        for (const auto& source : sources)
        {
          const double distance = (source - start).Norm();
          if (distance == 0.0)
            continue;
          segments.clear();
          ++traced;
          if (not tracer.Trace(cell_local_id, start, source, segments))
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
          if (not grid->CheckPointInsideCell(segments.back().cell_local_id, source))
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
  // A hexahedron whose top face is a strip 1e-4 wide: its side triangles are thin, but exact
  // predicates do not depend on their conditioning.
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
  bool reached_end = false;
  const Vector3 a(0.2, 0.3, 0.1);
  const Vector3 b(0.7, 0.52, 0.9);
  EXPECT_NEAR(TracedLength(tracer, a, b, reached_end), (b - a).Norm(), 1.0e-12);
  EXPECT_TRUE(reached_end);
  // Through the strip to the top: the segment ends on the top face.
  const Vector3 c(0.5, 0.5, 1.0);
  EXPECT_NEAR(TracedLength(tracer, a, c, reached_end), (c - a).Norm(), 1.0e-12);
  EXPECT_TRUE(reached_end);
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

TEST(MeshSegmentTracer, EndsOnBoundaryNearEdge)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  // Segments that end on the boundary face z = 0 just past the edge x = 0.25 of that face. The
  // ray crosses the face x = 0.25 at a height below the tracing tolerance, so the cell that it
  // enters there sees it on the boundary plane, although the part of the segment that remains is
  // longer than the tolerance. At the shallow slope, the crossing of the boundary plane is too
  // uncertain to use, so the end is found from the cell that holds the rest of the segment.
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  const MeshSegmentTracer tracer(*grid);
  for (const double slope : {0.1, 1.0e-4})
  {
    const double past_edge = 3.0e-11 / slope;
    const Vector3 end(0.25 + past_edge, 0.1, 0.0);
    const Vector3 direction = Vector3(1.0, 0.0, -slope).Normalized();
    const Vector3 start = end - 0.2 * direction;
    std::vector<MeshSegmentTracer::Segment> segments;
    EXPECT_TRUE(tracer.Trace(ContainingCells(*grid, start).front(), start, end, segments))
      << "slope " << slope;
    double total = 0.0;
    for (const auto& segment : segments)
      total += segment.length;
    EXPECT_NEAR(total, (end - start).Norm(), 1.0e-12) << "slope " << slope;
  }
}

TEST(MeshSegmentTracer, SliverNearInteriorEdge)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  // Three prisms around the edge x = y = 0, with faces in the directions 0, 150 and 240 degrees
  // from the edge. A segment that crosses the first face at a distance d from the edge, in the
  // direction 165 degrees, runs 1.93 d through the second prism before it enters the third. For
  // d between 1.1 and 2 tolerances, the second prism sees the segment on its exit face while the
  // third sees it clearly outside its face in the direction 240 degrees, so neither continues
  // the segment at the crossing unless the crossing of the exit face is computed. The distances
  // scanned cover that window for the tolerance of this mesh.
  const auto Rim = [](const double degrees)
  {
    return std::array<double, 2>{std::cos(degrees * M_PI / 180.0),
                                 std::sin(degrees * M_PI / 180.0)};
  };
  std::vector<Vector3> points;
  for (const double z : {0.0, 1.0})
  {
    points.emplace_back(0.0, 0.0, z);
    for (const double degrees : {0.0, 150.0, 240.0})
      points.emplace_back(Rim(degrees)[0], Rim(degrees)[1], z);
  }
  // Prism s spans the rim points s + 1 and s + 2 (cyclic) with the edge, counterclockwise.
  std::vector<std::vector<std::vector<size_t>>> cells;
  for (size_t s = 0; s < 3; ++s)
  {
    const size_t a = 1 + s;
    const size_t b = 1 + (s + 1) % 3;
    cells.push_back(
      {{0, b, a}, {4, 4 + a, 4 + b}, {0, a, 4 + a, 4}, {a, b, 4 + b, 4 + a}, {b, 0, 4, 4 + b}});
  }
  const auto grid = BuildPolyhedraMesh(points, cells);
  const MeshSegmentTracer tracer(*grid);
  const Vector3 direction(std::cos(165.0 * M_PI / 180.0), std::sin(165.0 * M_PI / 180.0), 0.0);
  for (int k = 0; k <= 40; ++k)
  {
    const double d = 1.0e-11 * std::pow(10.0, 2.0 * k / 40.0);
    const Vector3 crossing(d, 0.0, 0.5);
    const Vector3 start = crossing - 0.3 * direction;
    const Vector3 end = crossing + 0.3 * direction;
    std::vector<MeshSegmentTracer::Segment> segments;
    EXPECT_TRUE(tracer.Trace(ContainingCells(*grid, start).front(), start, end, segments))
      << "d = " << d;
    double total = 0.0;
    for (const auto& segment : segments)
      total += segment.length;
    EXPECT_NEAR(total, (end - start).Norm(), 1.0e-12) << "d = " << d;
    ASSERT_FALSE(segments.empty());
    EXPECT_EQ(segments.back().cell_local_id, 1) << "d = " << d;
  }
}

namespace
{

/// Segments traced near a feature, with the first few that the tracer lost.
struct NearTraces
{
  size_t traced = 0;
  size_t lost = 0;
  std::vector<std::string> examples;
};

/**
 * Traces the segment from a to b, which must lie in the mesh, starting in a cell that contains a
 * among the cells within two layers of the vertices given, or any cell. Counts the segment as
 * lost when the trace does not reach b or the lengths do not sum to the length of the segment.
 */
void
TraceNear(const Mesh& grid,
          const MeshSegmentTracer& tracer,
          const std::vector<std::vector<std::uint32_t>>& vertex_cells,
          const std::vector<std::uint64_t>& near_vertices,
          const Vector3& a,
          const Vector3& b,
          const std::string& kind,
          NearTraces& traces)
{
  std::optional<std::uint32_t> start;
  std::vector<std::uint32_t> candidates;
  for (const auto vid : near_vertices)
    candidates.insert(candidates.end(), vertex_cells[vid].begin(), vertex_cells[vid].end());
  const size_t first_layer = candidates.size();
  for (size_t i = 0; i < first_layer; ++i)
  {
    auto cell_vertex_ids = grid.GetCellConnectivity(candidates[i]);
    for (const auto vid : cell_vertex_ids)
      candidates.insert(candidates.end(), vertex_cells[vid].begin(), vertex_cells[vid].end());
  }
  for (const auto c : candidates)
    if (grid.CheckPointInsideCell(c, a))
    {
      start = c;
      break;
    }
  if (not start)
  {
    const auto cells = ContainingCells(grid, a);
    if (cells.empty())
      return;
    start = cells.front();
  }
  ++traces.traced;
  std::vector<MeshSegmentTracer::Segment> segments;
  const bool reached = tracer.Trace(*start, a, b, segments);
  double total = 0.0;
  for (const auto& segment : segments)
    total += segment.length;
  const double length = (b - a).Norm();
  if (not reached or std::abs(total - length) > 1.0e-12 * length)
  {
    ++traces.lost;
    if (traces.examples.size() < 5)
    {
      std::ostringstream text;
      text << std::setprecision(17) << kind << ": " << a.x << " " << a.y << " " << a.z << " -> "
           << b.x << " " << b.y << " " << b.z << (reached ? " traced " : " lost after ") << total
           << " of " << length;
      traces.examples.push_back(text.str());
    }
  }
}

} // namespace

TEST(MeshSegmentTracer, NearEdgesAndVerticesAtTolerance)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  // Segments that pass edges and vertices of a tetrahedral mesh at distances from a hundredth to
  // a hundred times the tracing tolerance, in many directions, and segments that end on the
  // boundary just past a boundary edge at slopes down to 1e-5. Near an edge, the distances of a
  // point to the planes through the edge straddle the tolerance, so the rules that decide the
  // cell the segment continues in must agree with each other at every distance.
  const auto grid = BuildMeshFromFile(Asset("cube3.2.msh"));
  Vector3 lo = grid->GlobalVertex(grid->GetCellConnectivity(0)[0]);
  Vector3 hi = lo;
  std::vector<std::vector<std::uint32_t>> vertex_cells(grid->GetGlobalVertexCount());
  std::map<std::pair<std::uint64_t, std::uint64_t>, bool> edges;
  for (std::uint32_t cell_local_id = 0; cell_local_id < grid->GetLocalCellCount(); ++cell_local_id)
  {
    auto cell_vertex_ids = grid->GetCellConnectivity(cell_local_id);
    const auto& cell = grid->GetLocalCell(cell_local_id);
    for (const auto vid : cell_vertex_ids)
    {
      vertex_cells[vid].push_back(cell_local_id);
      const auto& v = grid->GlobalVertex(vid);
      lo = Vector3(std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z));
      hi = Vector3(std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z));
    }
    auto cell_faces = grid->GetCellFaces(cell_local_id);
    for (std::uint32_t fi = 0; fi < cell_faces.size(); ++fi)
    {
      const auto face = cell_faces[fi];
      auto face_vertex_ids = grid->GetCellFaceConnectivity(cell_local_id, fi);
      const size_t n = face_vertex_ids.size();
      for (size_t s = 0; s < n; ++s)
      {
        const auto p = face_vertex_ids[s];
        const auto q = face_vertex_ids[(s + 1) % n];
        auto& boundary = edges[{std::min(p, q), std::max(p, q)}];
        boundary = boundary or not face.has_neighbor;
      }
    }
  }
  const double box = (hi - lo).Norm();
  const double half_length = 0.05 * box;
  // The tracing tolerance of these segments, up to the size of the start cell.
  const double tolerance = 1.0e-10 * 2.0 * half_length + 1.0e-13 * std::max(hi.Norm(), lo.Norm());
  const auto Inside = [&](const Vector3& p)
  {
    const double margin = 1.0e-3 * box;
    return p.x > lo.x + margin and p.x < hi.x - margin and p.y > lo.y + margin and
           p.y < hi.y - margin and p.z > lo.z + margin and p.z < hi.z - margin;
  };
  const auto Perpendicular = [](const Vector3& e)
  {
    const Vector3 trial = std::abs(e.x) < 0.9 ? Vector3(1.0, 0.0, 0.0) : Vector3(0.0, 1.0, 0.0);
    return e.Cross(trial).Normalized();
  };
  std::vector<double> offsets;
  for (int k = 0; k <= 16; ++k)
    offsets.push_back(tolerance * std::pow(10.0, -2.0 + 4.0 * k / 16.0));

  for (const bool split : {false, true})
  {
    const MeshSegmentTracer tracer(*grid, split);
    NearTraces traces;
    size_t edge_count = 0;
    for (const auto& [edge, boundary] : edges)
    {
      if (edge_count++ % 4 != 0)
        continue;
      const auto& v0 = grid->GlobalVertex(edge.first);
      const auto& v1 = grid->GlobalVertex(edge.second);
      const Vector3 m = 0.5 * (v0 + v1);
      const Vector3 e = (v1 - v0).Normalized();
      const Vector3 u = Perpendicular(e);
      const Vector3 v = e.Cross(u);
      for (const double offset : offsets)
        for (int k = 0; k < 6; ++k)
        {
          const auto q = SequencePoint(static_cast<int>(edge_count) * 6 + k);
          const double alpha = 2.0 * M_PI * q[0];
          const double beta = M_PI * (q[1] - 0.5);
          const Vector3 side = std::cos(alpha) * u + std::sin(alpha) * v;
          const Vector3 across = e.Cross(side);
          const Vector3 direction = std::cos(beta) * across + std::sin(beta) * e;
          const Vector3 p = m + offset * side;
          const Vector3 a = p - half_length * direction;
          const Vector3 b = p + half_length * direction;
          if (Inside(a) and Inside(b))
            TraceNear(*grid, tracer, vertex_cells, {edge.first, edge.second}, a, b, "edge", traces);
        }
      if (not boundary)
        continue;
      // Segments that end on the boundary face through the edge, just past the edge, at slopes
      // from steep to grazing.
      for (const int axis : {0, 1, 2})
        for (const double plane : {lo[axis], hi[axis]})
        {
          if (v0[axis] != plane or v1[axis] != plane)
            continue;
          const double into = plane == lo[axis] ? 1.0 : -1.0;
          const Vector3 inward(
            axis == 0 ? into : 0.0, axis == 1 ? into : 0.0, axis == 2 ? into : 0.0);
          const Vector3 along = e.Cross(inward);
          for (const double sign : {-1.0, 1.0})
            for (const double offset : offsets)
              for (const double slope : {0.3, 1.0e-2, 1.0e-3, 1.0e-5})
                for (int k = 0; k < 3; ++k)
                {
                  const auto q = SequencePoint(static_cast<int>(edge_count) * 3 + k);
                  const double alpha = 0.8 * M_PI * (q[2] - 0.5);
                  const Vector3 end = m + sign * offset * along;
                  const Vector3 in_plane = std::cos(alpha) * along + std::sin(alpha) * e;
                  const Vector3 direction = (-sign * in_plane + slope * inward).Normalized();
                  const Vector3 start = end + half_length * direction;
                  // The end must lie on the boundary face, inside the box.
                  bool on_face = true;
                  for (const int other : {0, 1, 2})
                    if (other != axis)
                      on_face = on_face and end[other] > lo[other] and end[other] < hi[other];
                  if (on_face and Inside(start))
                    TraceNear(*grid,
                              tracer,
                              vertex_cells,
                              {edge.first, edge.second},
                              start,
                              end,
                              "boundary end",
                              traces);
                }
        }
    }
    size_t vertex_count = 0;
    for (std::uint64_t vid = 0; vid < grid->GetGlobalVertexCount(); ++vid)
    {
      if (vertex_cells[vid].empty() or vertex_count++ % 4 != 0)
        continue;
      const auto& vertex = grid->GlobalVertex(vid);
      for (const double offset : offsets)
        for (int k = 0; k < 6; ++k)
        {
          const auto q = SequencePoint(static_cast<int>(vertex_count) * 6 + k);
          const Vector3 side = Vector3(q[0] - 0.5, q[1] - 0.5, q[2] - 0.5).Normalized();
          const Vector3 direction = Perpendicular(side).Cross(side);
          const Vector3 turned =
            (std::cos(M_PI * q[3]) * direction + std::sin(M_PI * q[3]) * Perpendicular(side))
              .Normalized();
          const Vector3 p = vertex + offset * side;
          const Vector3 a = p - half_length * turned;
          const Vector3 b = p + half_length * turned;
          if (Inside(a) and Inside(b))
            TraceNear(*grid, tracer, vertex_cells, {vid}, a, b, "vertex", traces);
        }
    }
    EXPECT_GT(traces.traced, 100000U);
    EXPECT_EQ(traces.lost, 0U) << "split=" << split;
    for (const auto& example : traces.examples)
      ADD_FAILURE() << "split=" << split << " " << example;
  }
}

TEST(MeshSegmentTracer, GrazingBoundaryEndFarFromOrigin)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  // An orthogonal mesh of 20 x 23 x 15 cells with every third layer 1e-4 times as thick, scaled
  // to 1e5 and placed at 1e8. The segment ends on the boundary face y = max, approaching it at a
  // cosine of 1.5e-5, and crosses an interior face 2e-4 before its end. At these coordinates the
  // round-off of a plane evaluation, divided by the cosine, exceeds the tolerance, so the crossing
  // of the boundary plane computed along the ray lands before the end of the segment, where no
  // cell continues the ray: the segment must still reach its end.
  const auto nodes = [](int n)
  {
    std::vector<double> x = {0.0};
    for (int i = 0; i < n; ++i)
      x.push_back(x.back() + (i % 3 == 1 ? 1.0e-4 : 1.0));
    for (auto& v : x)
      v = 1.0e8 + 1.0e5 * v / x.back();
    return x;
  };
  const auto grid = BuildOrthogonalMesh({nodes(20), nodes(23), nodes(15)});
  const MeshSegmentTracer tracer(*grid);
  const Vector3 start(100055724.88939425, 100099999.69010857, 100049999.74423848);
  const Vector3 end(100046153.66863474, 100100000.0, 100068826.83724038);
  const auto cells = ContainingCells(*grid, start);
  ASSERT_FALSE(cells.empty());
  std::vector<MeshSegmentTracer::Segment> segments;
  EXPECT_TRUE(tracer.Trace(cells.front(), start, end, segments));
  double total = 0.0;
  for (const auto& segment : segments)
    total += segment.length;
  EXPECT_NEAR(total, (end - start).Norm(), 1.0e-9 * (end - start).Norm());
}

TEST(MeshSegmentTracer, AlongSharedFaceFarFromOrigin)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  // Two tetrahedra of size 0.03 that share a face, placed 3000 away from the origin, so that the
  // round-off of the face normals computed from different side triangles (about 1e-11) is far
  // above any fixed cosine threshold. Segments in the plane of the shared face, traced from
  // either cell, must not see the ray leaving both cells.
  const Vector3 o(1000.0, -2000.0, 3000.0);
  const std::vector<Vector3> points = {o + Vector3(0.0, 0.0, 0.0),
                                       o + Vector3(0.031, 0.002, 0.0),
                                       o + Vector3(0.012, 0.029, 0.0),
                                       o + Vector3(0.011, 0.009, 0.027),
                                       o + Vector3(0.013, 0.012, -0.025)};
  // Faces counterclockwise seen from outside. Cell 0 has apex 3 above the shared face (0, 1, 2),
  // cell 1 has apex 4 below it.
  const std::vector<std::vector<std::vector<size_t>>> cells = {
    {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}, {{0, 1, 2}, {0, 4, 1}, {1, 4, 2}, {2, 4, 0}}};
  const auto grid = BuildPolyhedraMesh(points, cells);
  const MeshSegmentTracer tracer(*grid);
  for (int k = 0; k < 24; ++k)
  {
    const auto q = SequencePoint(k);
    // Two points in the shared face, with barycentric coordinates away from its edges.
    const auto in_face = [&](double u, double v)
    {
      const double w0 = 0.1 + 0.8 * u * (1.0 - v);
      const double w1 = 0.1 + 0.8 * (1.0 - u) * (1.0 - v);
      const double w2 = 1.0 - w0 - w1;
      return w0 * points[0] + w1 * points[1] + w2 * points[2];
    };
    const Vector3 a = in_face(q[0], 0.5 * q[1]);
    const Vector3 b = in_face(q[2], 0.5 * q[3]);
    if ((b - a).Norm() < 1.0e-3)
      continue;
    for (const std::uint32_t hint : {0U, 1U})
    {
      std::vector<MeshSegmentTracer::Segment> segments;
      EXPECT_TRUE(tracer.Trace(hint, a, b, segments)) << "k=" << k << " hint " << hint;
      double total = 0.0;
      for (const auto& segment : segments)
        total += segment.length;
      EXPECT_NEAR(total, (b - a).Norm(), 1.0e-9 * (b - a).Norm()) << "k=" << k;
    }
  }
}

TEST(MeshSegmentTracer, FromVertexIntoEveryOctant)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const auto grid = BuildBoxMesh(1.0, 4, 0.0);
  const MeshSegmentTracer tracer(*grid);
  // From the vertex (0.5, 0.5, 0.5), shared by eight cells, into each octant, with any of the
  // eight cells as the start cell: the segment begins in the cell of its octant.
  const Vector3 vertex(0.5, 0.5, 0.5);
  for (const auto hint : ContainingCells(*grid, vertex))
    for (const double sx : {-1.0, 1.0})
      for (const double sy : {-1.0, 1.0})
        for (const double sz : {-1.0, 1.0})
        {
          const Vector3 omega = Vector3(sx * 0.3, sy * 0.5, sz * 0.8).Normalized();
          std::vector<MeshSegmentTracer::Segment> segments;
          ASSERT_TRUE(tracer.Trace(hint, vertex, vertex + 0.1 * omega, segments))
            << "hint " << hint << " along " << omega.PrintStr();
          // A start at a vertex is traced with displaced ends, which may enter the other cells at
          // the vertex over lengths of the order of the displacement.
          ASSERT_FALSE(segments.empty()) << "hint " << hint << " along " << omega.PrintStr();
          const auto longest =
            std::max_element(segments.begin(),
                             segments.end(),
                             [](const auto& s, const auto& t) { return s.length < t.length; });
          EXPECT_TRUE(grid->CheckPointInsideCell(longest->cell_local_id, vertex + 0.01 * omega))
            << "hint " << hint << " along " << omega.PrintStr();
          EXPECT_NEAR(longest->length, 0.1, 1.0e-12);
          for (const auto& segment : segments)
            if (&segment != &*longest)
              EXPECT_LE(segment.length, 1.0e-12)
                << "hint " << hint << " along " << omega.PrintStr();
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
  // The dented cube: the average of its vertices is the apex of the dent, which lies in the
  // plane of every dent triangle, so the tracer must choose another interior point for the side
  // tetrahedra. The coordinates are exactly representable, so that the segments below lie
  // exactly in the plane z = 1 - y of one dent triangle: a segment on the boundary counts as
  // inside the mesh.
  const auto grid = BuildDentedCube();
  const auto& cell = grid->GetLocalCell(0);
  EXPECT_NEAR((cell.centroid - Vector3(0.5, 0.5, 0.5)).Norm(), 0.0, 1.0e-15);
  const MeshSegmentTracer tracer(*grid);
  EXPECT_FALSE(tracer.IsConvexCell(0));
  bool reached_end = false;
  // Along the dent plane from inside the cube out through the face x = 1 at half the length.
  const Vector3 a(0.75, 0.375, 0.625);
  const Vector3 b(1.25, 0.125, 0.875);
  EXPECT_NEAR(TracedLength(tracer, a, b, reached_end), 0.5 * (b - a).Norm(), 1.0e-12);
  EXPECT_FALSE(reached_end);
  // The same plane, crossing the crease between two dent triangles on the way, out through
  // x = 1 at 7/9 of the length.
  const Vector3 c(0.125, 0.375, 0.625);
  const Vector3 d(1.25, 0.125, 0.875);
  EXPECT_NEAR(TracedLength(tracer, c, d, reached_end), 7.0 / 9.0 * (d - c).Norm(), 1.0e-12);
  EXPECT_FALSE(reached_end);
  // Below the dent from side to side.
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

  const auto box = BuildBoxMesh(1.0, 2, z0);
  for (const bool split : {false, true})
    CheckMesh(*box, MeshSegmentTracer(*box, split), z0, z0 + 1.0);
}

// Segments of length 1e-12 that end on a tilted shared face 3000 away from the origin, where the
// end rounds to within a few 1e-13 of the face plane and the floating-point crossing of the plane
// is uncertain by as much: whether a segment ends in a piece must be decided exactly, in both
// directions.
TEST(MeshSegmentTracer, TinySegmentsAcrossTiltedFace)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const Vector3 o(1000.0, -2000.0, 3000.0);
  // The shared face (0, 1, 2) is tilted, so that a point computed on it rounds to either side.
  const std::vector<Vector3> points = {o + Vector3(0.0, 0.0, 0.0),
                                       o + Vector3(0.031, 0.002, 0.004),
                                       o + Vector3(0.012, 0.029, -0.003),
                                       o + Vector3(0.011, 0.009, 0.027),
                                       o + Vector3(0.013, 0.012, -0.025)};
  const std::vector<std::vector<std::vector<size_t>>> cells = {
    {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}, {{0, 1, 2}, {0, 4, 1}, {1, 4, 2}, {2, 4, 0}}};
  const auto grid = BuildPolyhedraMesh(points, cells);
  const MeshSegmentTracer tracer(*grid);
  const Vector3 normal = (points[1] - points[0]).Cross(points[2] - points[0]).Normalized();
  for (int k = 0; k < 200; ++k)
  {
    const auto q = SequencePoint(k);
    const double w0 = 0.1 + 0.8 * q[0] * (1.0 - q[1]);
    const double w1 = 0.1 + 0.8 * (1.0 - q[0]) * (1.0 - q[1]);
    const Vector3 x = w0 * points[0] + w1 * points[1] + (1.0 - w0 - w1) * points[2];
    // A direction with a component across the face between 0.05 and 1.
    Vector3 d = Vector3(q[2] - 0.5, q[3] - 0.5, q[4] - 0.5);
    if (d.Norm() < 1.0e-3)
      d = normal;
    d = d.Normalized();
    if (std::abs(d.Dot(normal)) < 0.05)
      d = (d + normal).Normalized();
    const Vector3 a = x - 1.0e-12 * d;
    for (const auto& [start, end] : {std::pair{a, x}, std::pair{x, a}})
      for (const std::uint32_t hint : {0U, 1U})
      {
        std::vector<MeshSegmentTracer::Segment> segments;
        EXPECT_TRUE(tracer.Trace(hint, start, end, segments)) << "k=" << k << " hint " << hint;
        double total = 0.0;
        for (const auto& segment : segments)
          total += segment.length;
        // The ends round to the coordinate grid, about 5e-13 here, so the length is not 1e-12.
        EXPECT_NEAR(total, (end - start).Norm(), 1.0e-9 * (end - start).Norm())
          << "k=" << k << " hint " << hint;
      }
  }
}

// A segment of length 1e-3 whose end lies 3 units in the last place outside a boundary face,
// approached at a shallow angle (cosine 0.02), 1e5 away from the origin. The end is within
// round-off of the face, so the segment reaches it; its distance along the line past the crossing
// is 150 units in the last place, which is not the measure of that round-off.
TEST(MeshSegmentTracer, ObliqueEndOutsideBoundaryByRoundoff)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const double o = 1.0e5;
  const auto grid = BuildBoxMesh(1.0e5, 2, o);
  const double top = o + 1.0e5;
  const double ulp = std::ldexp(1.0, std::ilogb(top) - 52);
  for (const bool split : {false, true})
  {
    const MeshSegmentTracer tracer(*grid, split);
    for (int k = 0; k < 50; ++k)
    {
      const auto q = SequencePoint(k);
      // The end just outside the top face, the start 1e-3 away inside at a shallow angle.
      const Vector3 end(o + 2.0e4 + 6.0e4 * q[0], o + 2.0e4 + 6.0e4 * q[1], top + 3.0 * ulp);
      const double angle = 2.0 * M_PI * q[2];
      const Vector3 along(std::cos(angle), std::sin(angle), 0.0);
      const Vector3 start =
        end - 1.0e-3 * (std::sqrt(1.0 - 0.02 * 0.02) * along + 0.02 * Vector3(0.0, 0.0, 1.0));
      const auto hint = ContainingCells(*grid, start).front();
      std::vector<MeshSegmentTracer::Segment> segments;
      EXPECT_TRUE(tracer.Trace(hint, start, end, segments)) << "k=" << k << " split " << split;
      double total = 0.0;
      for (const auto& segment : segments)
        total += segment.length;
      EXPECT_NEAR(total, (end - start).Norm(), 1.0e-9 * (end - start).Norm())
        << "k=" << k << " split " << split;
    }
  }
}

// The orthogonal mesh of a fuzzing failure: layers of relative thickness 1 and 1e-4 over 1e5,
// placed at 1e8. A segment in an internal face plane (traced with displaced ends) that runs 7.5e-7
// below the top boundary, within the round-off allowance of that plane, must put every piece of
// its length in the cell that contains it: a displaced line that leaves through the top at a
// shallow angle is not an end on the boundary.
TEST(MeshSegmentTracer, ShallowLineBelowTopBoundaryFarFromOrigin)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  std::vector<std::vector<double>> sets;
  for (const int n : {23, 23, 3})
  {
    std::vector<double> nodes = {0.0};
    double x = 0.0;
    for (int i = 0; i < n; ++i)
    {
      x += (i % 3 == 1) ? 1.0e-4 : 1.0;
      nodes.push_back(x);
    }
    const double len = nodes.back();
    for (auto& v : nodes)
      v = 1.0e5 * v / len + 1.0e8;
    sets.push_back(nodes);
  }
  const auto grid = BuildOrthogonalMesh(sets);
  const MeshSegmentTracer tracer(*grid);
  const Vector3 a(100073333.73753235, 100039999.86667378, 100099999.99999925);
  const Vector3 b(100073332.44027151, 100039999.86667378, 100099999.99999925);
  const auto hint = ContainingCells(*grid, a).front();
  std::vector<MeshSegmentTracer::Segment> segments;
  ASSERT_TRUE(tracer.Trace(hint, a, b, segments));
  // Three cells along x, the middle one a thin layer.
  EXPECT_EQ(segments.size(), 3);
  double t = 0.0;
  const double length = (b - a).Norm();
  for (const auto& segment : segments)
  {
    const Vector3 midpoint = a + (t + 0.5 * segment.length) / length * (b - a);
    EXPECT_TRUE(grid->CheckPointInsideCell(segment.cell_local_id, midpoint))
      << "cell " << segment.cell_local_id;
    t += segment.length;
  }
  EXPECT_NEAR(t, length, 1.0e-9 * length);
}

// An L-shaped polyhedron (two merged boxes) from a jittered grid, with non-planar faces. It is not
// star-shaped about the average of its points, nor about the points between that average and its
// face centroids, but about a small kernel of points elsewhere (radius about 0.01 of a cell of
// size 0.35), which the tracer must find to trace the cell through side tetrahedra.
TEST(MeshSegmentTracer, JitteredLCellStarShapedOffCentroid)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const std::vector<Vector3> points = {
    Vector3(0.12651251618053377, 0.5633351652742743, 0.7494488069072369),
    Vector3(0.1250037856818183, 0.5627314233111153, 0.8867048961242754),
    Vector3(0.14697618953061467, 0.7390600820631072, 0.7239082929689022),
    Vector3(0.17547996125567442, 0.7228486043478811, 0.8356187176368864),
    Vector3(0.29571629714748804, 0.5721436989020123, 0.7277635453510685),
    Vector3(0.3020344146589612, 0.5407426350683723, 0.8815673391169403),
    Vector3(0.2730620141126997, 0.6855197990453811, 0.7427591413241965),
    Vector3(0.27535813620574706, 0.7072208689177655, 0.8549516075604895),
    Vector3(0.28075116172628106, 0.8498276674777979, 0.68929688443365),
    Vector3(0.25814434970112043, 0.8528660477413939, 0.8295267619961877),
    Vector3(0.4204496635896583, 0.606325876127527, 0.7471416998626698),
    Vector3(0.4217131677578642, 0.6016808544199302, 0.8630290050224434),
    Vector3(0.4276577790740371, 0.6857787674637315, 0.7434771999925175),
    Vector3(0.4295667088502454, 0.7496755381894493, 0.8292072727371338),
    Vector3(0.44069437967411224, 0.8702730603089949, 0.7469207830866319),
    Vector3(0.46426634445300086, 0.838343275280225, 0.8817461954348187),
  };
  const std::vector<std::vector<std::vector<size_t>>> cells = {{{0, 2, 6, 4},
                                                                {1, 5, 7, 3},
                                                                {0, 4, 5, 1},
                                                                {6, 2, 3, 7},
                                                                {2, 0, 1, 3},
                                                                {4, 6, 12, 10},
                                                                {5, 11, 13, 7},
                                                                {4, 10, 11, 5},
                                                                {10, 12, 13, 11},
                                                                {6, 8, 14, 12},
                                                                {7, 13, 15, 9},
                                                                {12, 14, 15, 13},
                                                                {14, 8, 9, 15},
                                                                {8, 6, 7, 9}}};
  const auto grid = BuildPolyhedraMesh(points, cells);
  std::unique_ptr<MeshSegmentTracer> tracer;
  ASSERT_NO_THROW(tracer = std::make_unique<MeshSegmentTracer>(*grid));
  // Segments between points of the kernel region and points near the vertices of each box.
  const Vector3 kernel(0.297, 0.666, 0.765);
  for (int k = 0; k < 40; ++k)
  {
    const auto q = SequencePoint(k);
    const Vector3 target = 0.9 * points[static_cast<size_t>(q[0] * 16.0) % 16] + 0.1 * kernel;
    const Vector3 a = kernel + 0.3 * (q[1] - 0.5) * (target - kernel);
    std::vector<MeshSegmentTracer::Segment> segments;
    EXPECT_TRUE(tracer->Trace(0, a, target, segments)) << "k=" << k;
    double total = 0.0;
    for (const auto& segment : segments)
      total += segment.length;
    EXPECT_NEAR(total, (target - a).Norm(), 1.0e-9 * (target - a).Norm()) << "k=" << k;
  }
}

// A pentagonal prism 1e9 away from the origin whose pentagonal boundary face is planar only to
// round-off: one rim vertex lies two units in the last place beyond the others, as the vertices
// of a cell clipped to a box do. The face centroid, the rounded average, lies with the lower
// vertices, so four of the five fan triangles lie two units inside the box plane. A segment that
// runs along the face one unit inside the box plane then lies outside those triangles, but it
// lies on the boundary within round-off and must reach its end in both directions.
TEST(MeshSegmentTracer, AlongBoundaryFaceWithRoundedCentroidFarFromOrigin)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const double o = 1.0e9;
  const double ulp = std::ldexp(1.0, std::ilogb(o) - 52);
  for (int trial = 0; trial < 4; ++trial)
  {
    const double X = o + 1000.0 + 7.0 * trial * ulp;
    std::vector<Vector3> points;
    for (const double x : {X - 2000.0, X})
      for (int k = 0; k < 5; ++k)
      {
        const double angle = 2.0 * M_PI * k / 5.0;
        const double rim = (x == X and k != 0) ? X - 2.0 * ulp : x;
        points.emplace_back(rim, o + 1000.0 * std::cos(angle), o + 1000.0 * std::sin(angle));
      }
    // Faces counterclockwise seen from outside: the two pentagons and five quads.
    std::vector<std::vector<std::vector<size_t>>> cells = {{{0, 4, 3, 2, 1}, {5, 6, 7, 8, 9}}};
    for (size_t k = 0; k < 5; ++k)
      cells[0].push_back({k, (k + 1) % 5, 5 + (k + 1) % 5, 5 + k});
    const auto grid = BuildPolyhedraMesh(points, cells);
    const MeshSegmentTracer tracer(*grid);
    const double x = X - ulp;
    // Segments across the pentagon, one unit in the last place inside the box plane.
    for (int k = 0; k < 10; ++k)
    {
      const auto q = SequencePoint(k);
      const double angle = 2.0 * M_PI * q[0];
      const Vector3 a(x, o + 600.0 * std::cos(angle), o + 600.0 * std::sin(angle));
      const Vector3 b(x, o - 600.0 * std::cos(angle), o - 600.0 * std::sin(angle));
      for (const auto& [start, end] : {std::pair{a, b}, std::pair{b, a}})
      {
        std::vector<MeshSegmentTracer::Segment> segments;
        EXPECT_TRUE(tracer.Trace(0, start, end, segments)) << "trial " << trial << " k=" << k;
        double total = 0.0;
        for (const auto& segment : segments)
          total += segment.length;
        EXPECT_NEAR(total, (end - start).Norm(), 1.0e-9 * (end - start).Norm())
          << "trial " << trial << " k=" << k;
      }
    }
  }
}

// A box 1e-4 wide in x, 1 in y and z, placed at z = 1e12. The round-off allowance along z is
// about 1e-4, the whole width of the box, while along x it is 1e-20: a degenerate segment (through
// the centers of the faces normal to y) displaced by the z allowance along x would leave the box
// for every direction. Allowances are per axis.
TEST(MeshSegmentTracer, AnisotropicBoxFarFromOriginAlongOneAxis)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const double z0 = 1.0e12;
  const auto grid =
    BuildOrthogonalMesh({{0.0, 0.5e-4, 1.0e-4}, {0.0, 0.5, 1.0}, {z0, z0 + 0.5, z0 + 1.0}});
  for (const bool split : {false, true})
  {
    const MeshSegmentTracer tracer(*grid, split);
    const Vector3 a(0.5e-4, 0.25, z0 + 0.5);
    const Vector3 b(0.5e-4, 0.75, z0 + 0.5);
    const auto hint = ContainingCells(*grid, a).front();
    std::vector<MeshSegmentTracer::Segment> segments;
    EXPECT_TRUE(tracer.Trace(hint, a, b, segments)) << "split " << split;
    double total = 0.0;
    for (const auto& segment : segments)
      total += segment.length;
    EXPECT_NEAR(total, 0.5, 1.0e-9) << "split " << split;
  }
}

// Tetrahedra whose solid angle at the vertex (0, 0, 0) is thin, with half-widths h down to 1e-6,
// around the diagonal edge from that vertex to (1, 1, 1). A segment along that edge meets edges
// exactly, and every octant displacement of its start (a few units in the last place along each
// axis, all of extent 1) lies outside the cone once h is small enough; the displacement toward
// the interior of the cell does not.
TEST(MeshSegmentTracer, AlongEdgeOfThinTetrahedron)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  for (const double h : {1.0e-2, 1.0e-4, 1.0e-6})
  {
    const Vector3 u(1.0, 1.0, 1.0);
    const std::vector<Vector3> points = {
      Vector3(0.0, 0.0, 0.0), u, u + h * Vector3(1.0, -1.0, 0.0), u + h * Vector3(1.0, 1.0, -2.0)};
    // Faces counterclockwise seen from outside.
    const std::vector<std::vector<std::vector<size_t>>> cells = {
      {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}}};
    const auto grid = BuildPolyhedraMesh(points, cells);
    const MeshSegmentTracer tracer(*grid);
    const Vector3 a(0.0, 0.0, 0.0);
    const Vector3 b = 0.5 * u;
    for (const auto& [start, end] : {std::pair{a, b}, std::pair{b, a}})
    {
      std::vector<MeshSegmentTracer::Segment> segments;
      EXPECT_TRUE(tracer.Trace(0, start, end, segments)) << "h " << h;
      double total = 0.0;
      for (const auto& segment : segments)
        total += segment.length;
      EXPECT_NEAR(total, b.Norm(), 1.0e-12) << "h " << h;
    }
  }
}

// Three convex cells: two tetrahedra whose cones at the apexes (0, 0, 0) and (1, 0, 0) are thin
// (height 1e-10 over a base of 0.5), joined by the convex hull of their bases and four far points.
// The segment between the apexes runs along the axis of both cones, through vertices and in face
// planes. Displacements of the far end toward the interior of the start cell would miss the cone
// at the far end, and octant displacements leave it; the displacement along the segment stays in
// the closed mesh.
TEST(MeshSegmentTracer, BetweenApexesOfThinCones)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const double h = 1.0e-10;
  const std::vector<Vector3> points = {Vector3(0.0, 0.0, 0.0),
                                       Vector3(0.25, -0.25, 0.0),
                                       Vector3(0.25, 0.25, 0.0),
                                       Vector3(0.25, 0.0, h),
                                       Vector3(0.75, -0.25, 0.0),
                                       Vector3(0.75, 0.25, 0.0),
                                       Vector3(0.75, 0.0, -h),
                                       Vector3(0.5, -1.0, -1.0),
                                       Vector3(0.5, -1.0, 1.0),
                                       Vector3(0.5, 1.0, -1.0),
                                       Vector3(0.5, 1.0, 1.0),
                                       Vector3(1.0, 0.0, 0.0)};
  // Faces counterclockwise seen from outside each cell; the caps (1, 2, 3) and (4, 5, 6) are
  // shared between a tetrahedron and the hull.
  const std::vector<std::vector<std::vector<size_t>>> cells = {
    {{3, 1, 2}, {1, 3, 0}, {2, 1, 0}, {3, 2, 0}},
    {{8, 1, 7},
     {7, 4, 8},
     {9, 2, 10},
     {10, 5, 9},
     {2, 1, 3},
     {3, 1, 8},
     {10, 3, 8},
     {10, 2, 3},
     {9, 6, 7},
     {6, 4, 7},
     {9, 5, 6},
     {5, 4, 6},
     {7, 1, 2, 9},
     {8, 4, 5, 10}},
    {{6, 4, 5}, {4, 6, 11}, {5, 4, 11}, {6, 5, 11}}};
  const auto grid = BuildPolyhedraMesh(points, cells);
  const MeshSegmentTracer tracer(*grid);
  // The apexes exactly, and the apexes rounded outward by one unit in the last place, as a caller
  // computing them could produce: the cells that hold such an end within round-off must still
  // supply the inward directions.
  const double outside = std::nextafter(1.0, 2.0) - 1.0;
  for (const auto& [a, b] :
       {std::pair{Vector3(0.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0)},
        std::pair{Vector3(-outside, 0.0, 0.0), Vector3(1.0 + outside, 0.0, 0.0)}})
    for (const auto& [start, end, hint] : {std::tuple{a, b, 0U}, std::tuple{b, a, 2U}})
    {
      std::vector<MeshSegmentTracer::Segment> segments;
      EXPECT_TRUE(tracer.Trace(hint, start, end, segments)) << "from " << start.PrintStr();
      double total = 0.0;
      for (const auto& segment : segments)
        total += segment.length;
      EXPECT_NEAR(total, (b - a).Norm(), 1.0e-12) << "from " << start.PrintStr();
    }
}

// A unit cube whose top is dented to a pyramid with its apex at z = 0.1: five planar square
// faces, whose fans share a plane four times, and four triangles. The average of its points lies
// in the cavity, outside the cell, so the kernel point (below the apex) must be found by the
// search, whose ascent at a ridge between a square and a dent plane needs the two distinct planes.
TEST(MeshSegmentTracer, DentedCubeKernelSearch)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const std::vector<Vector3> points = {Vector3(0.0, 0.0, 0.0),
                                       Vector3(1.0, 0.0, 0.0),
                                       Vector3(1.0, 1.0, 0.0),
                                       Vector3(0.0, 1.0, 0.0),
                                       Vector3(0.0, 0.0, 1.0),
                                       Vector3(1.0, 0.0, 1.0),
                                       Vector3(1.0, 1.0, 1.0),
                                       Vector3(0.0, 1.0, 1.0),
                                       Vector3(0.5, 0.5, 0.1)};
  // Faces counterclockwise seen from outside.
  const std::vector<std::vector<std::vector<size_t>>> cells = {{{0, 3, 2, 1},
                                                                {0, 1, 5, 4},
                                                                {1, 2, 6, 5},
                                                                {2, 3, 7, 6},
                                                                {3, 0, 4, 7},
                                                                {4, 5, 8},
                                                                {5, 6, 8},
                                                                {6, 7, 8},
                                                                {7, 4, 8}}};
  const auto grid = BuildPolyhedraMesh(points, cells);
  std::unique_ptr<MeshSegmentTracer> tracer;
  ASSERT_NO_THROW(tracer = std::make_unique<MeshSegmentTracer>(*grid));
  for (const auto& [a, b] : {std::pair{Vector3(0.1, 0.1, 0.05), Vector3(0.9, 0.9, 0.05)},
                             std::pair{Vector3(0.05, 0.5, 0.9), Vector3(0.05, 0.5, 0.02)},
                             std::pair{Vector3(0.1, 0.9, 0.8), Vector3(0.1, 0.1, 0.8)}})
  {
    std::vector<MeshSegmentTracer::Segment> segments;
    EXPECT_TRUE(tracer->Trace(0, a, b, segments)) << a.PrintStr();
    double total = 0.0;
    for (const auto& segment : segments)
      total += segment.length;
    EXPECT_NEAR(total, (b - a).Norm(), 1.0e-12 * (b - a).Norm()) << a.PrintStr();
  }
}

// A triangular prism whose bottom face is replaced by three triangles meeting at the origin, far
// up inside the prism: the cell is star-shaped about w / 32, w = (1, 1, 8), but not about the
// average of its points, and at the corner of the search function reached from that average every
// ascent direction that the climb tries decreases one of the three active distances. The complete
// fallback must find the kernel.
TEST(MeshSegmentTracer, DentedPrismKernelFallback)
{
  if (opensn::mpi_comm.size() != 1)
    return;
  const Vector3 w(1.0, 1.0, 8.0);
  const std::vector<Vector3> rim = {Vector3(-1.0, -3.0, -20.0) / 16.0,
                                    Vector3(-3.0, -1.0, -20.0) / 16.0,
                                    Vector3(-4.0, -4.0, -24.0) / 16.0};
  std::vector<Vector3> points = rim;
  for (const auto& v : rim)
    points.push_back(v + 0.25 * w);
  points.emplace_back(0.0, 0.0, 0.0);
  // Faces counterclockwise seen from outside: top, three sides, three triangles of the dent.
  const std::vector<std::vector<std::vector<size_t>>> cells = {
    {{3, 4, 5}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {1, 0, 6}, {2, 1, 6}, {0, 2, 6}}};
  const auto grid = BuildPolyhedraMesh(points, cells);
  std::unique_ptr<MeshSegmentTracer> tracer;
  ASSERT_NO_THROW(tracer = std::make_unique<MeshSegmentTracer>(*grid));
  // Segments along and near the axis of the prism, above the apex of the dent.
  const Vector3 kernel = w / 32.0;
  for (const auto& [a, b] : {std::pair{Vector3(0.0, 0.0, 0.02), Vector3(0.0625, 0.0625, 0.5)},
                             std::pair{kernel, Vector3(0.03, 0.02, 0.3)},
                             std::pair{Vector3(0.07, 0.07, 0.55), Vector3(0.005, 0.01, 0.1)}})
  {
    std::vector<MeshSegmentTracer::Segment> segments;
    EXPECT_TRUE(tracer->Trace(0, a, b, segments)) << a.PrintStr();
    double total = 0.0;
    for (const auto& segment : segments)
      total += segment.length;
    EXPECT_NEAR(total, (b - a).Norm(), 1.0e-12 * (b - a).Norm()) << a.PrintStr();
  }
}
