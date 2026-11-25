// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/error.h"
#include <algorithm>
#include <array>
#include <map>
#include <cmath>
#include <sstream>
#include <string>

namespace opensn
{

namespace
{

/**
 * A cell whose area or volume is at most this fraction of the magnitude of the terms summed to
 * compute it is zero to within rounding, and degenerate. It is about 4500 times the unit roundoff.
 */
constexpr double degenerate_measure_tolerance = 1.0e-12;

/// Returns the center of a face, the average of its vertices, relative to `origin`.
Vector3
FaceCenter(const UnpartitionedMesh::LightWeightFace& face,
           const std::vector<Vector3>& vertices,
           const Vector3& origin = Vector3())
{
  Vector3 center;
  for (const auto vid : face.vertex_ids)
    center += vertices[vid] - origin;
  return center / static_cast<double>(face.vertex_ids.size());
}

/**
 * Returns the area vector of a face of a polygon or polyhedron, oriented by the order of its
 * vertices. For a polygon, whose faces are edges in the x-y plane, it is the edge length times the
 * normal that points to the right of the edge. For a polyhedron, it is the sum of the area vectors
 * of the triangles that join each edge of the face to its center. The faces of a valid cell are
 * ordered so that this vector points out of the cell. The vertices are taken relative to `origin`,
 * which changes the vector only by rounding.
 */
Vector3
FaceAreaVector(const CellType cell_type,
               const UnpartitionedMesh::LightWeightFace& face,
               const std::vector<Vector3>& vertices,
               const Vector3& origin = Vector3())
{
  if (cell_type == CellType::POLYGON)
  {
    const auto& v0 = vertices[face.vertex_ids[0]];
    const auto& v1 = vertices[face.vertex_ids[1]];
    return (v1 - v0).Cross(Vector3(0.0, 0.0, 1.0));
  }

  const auto center = FaceCenter(face, vertices, origin);
  const size_t num_vertices = face.vertex_ids.size();
  Vector3 area;
  for (size_t i = 0; i < num_vertices; ++i)
  {
    const auto v0 = vertices[face.vertex_ids[i]] - origin;
    const auto v1 = vertices[face.vertex_ids[(i + 1) % num_vertices]] - origin;
    area += 0.5 * (v0 - center).Cross(v1 - center);
  }
  return area;
}

/// Signed area of a polygon or signed volume of a polyhedron, with the scale of its rounding error.
struct CellMeasure
{
  /**
   * Positive when the faces are ordered outward, negative when the cell is inverted, and zero when
   * it is degenerate.
   */
  double value = 0.0;
  /**
   * Sum of the magnitudes of the terms summed for the value. The rounding error of the value is
   * a small multiple of the unit roundoff times this scale. It is much smaller than the cube of
   * the longest edge for a needle-shaped cell, whose faces are thin.
   */
  double scale = 0.0;
};

/**
 * Returns the signed area of a polygon or the signed volume of a polyhedron, from the divergence
 * theorem over its faces. Concave cells are measured exactly.
 *
 * The faces are measured relative to a vertex of the cell. Far from the origin, absolute
 * coordinates round with errors much larger than the volume of a thin cell, while the difference
 * of two nearby coordinates is exact. The divergence theorem holds for any reference point, so the
 * vertex is also the point that the face arms start from. A cell without vertices has measure
 * zero.
 */
CellMeasure
SignedCellMeasure(const UnpartitionedMesh::LightWeightCell& cell,
                  const std::vector<Vector3>& vertices)
{
  const double dimension = cell.type == CellType::POLYGON ? 2.0 : 3.0;
  CellMeasure measure;
  if (cell.vertex_ids.empty())
    return measure;
  const Vector3& origin = vertices[cell.vertex_ids.front()];
  for (const auto& face : cell.faces)
  {
    const Vector3 arm = FaceCenter(face, vertices, origin);
    const Vector3 area = FaceAreaVector(cell.type, face, vertices, origin);
    measure.value += arm.Dot(area);
    measure.scale += arm.Norm() * area.Norm();
  }
  measure.value /= dimension;
  measure.scale /= dimension;
  return measure;
}

/**
 * Two faces closer than this fraction of the mesh's bounding-box diagonal, in their plane offset
 * and in their overlap, coincide. Vertex coordinates written with 10 to 12 digits round to well
 * below it.
 */
constexpr double coincidence_tolerance = 1.0e-8;

/// Faces of the mesh that overlap a face of another cell without sharing its vertices.
struct HangingFaces
{
  size_t count = 0;
  uint64_t first_cell = 0;
  Vector3 first_center;
};

/// Returns whether point p lies strictly inside the polygon `polygon` of the plane, by winding.
bool
StrictlyInsidePolygon(const std::vector<std::array<double, 2>>& polygon,
                      const std::array<double, 2>& p,
                      const double tolerance)
{
  // The point must be inside, and farther than the tolerance from every edge.
  bool inside = false;
  const size_t n = polygon.size();
  for (size_t i = 0, j = n - 1; i < n; j = i++)
  {
    const auto& a = polygon[i];
    const auto& b = polygon[j];
    const double ex = b[0] - a[0];
    const double ey = b[1] - a[1];
    const double length = std::hypot(ex, ey);
    if (length > 0.0)
    {
      const double s =
        std::clamp(((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / (length * length), 0.0, 1.0);
      if (std::hypot(p[0] - a[0] - s * ex, p[1] - a[1] - s * ey) <= tolerance)
        return false;
    }
    if (((a[1] > p[1]) != (b[1] > p[1])) and
        (p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]))
      inside = not inside;
  }
  return inside;
}

/**
 * Finds the faces without a neighbor that overlap, with positive length or area, a face without a
 * neighbor of another cell lying in the same line or plane and facing the other way. In a
 * conforming mesh two such faces are the same face and are connected. When the vertices of one
 * cell lie inside an edge or face of its neighbor (a hanging node or T-junction), the two faces
 * are not equal, stay unconnected, and would become boundaries inside the domain.
 */
HangingFaces
FindHangingFaces(const std::vector<std::shared_ptr<UnpartitionedMesh::LightWeightCell>>& cells,
                 const std::vector<Vector3>& vertices)
{
  HangingFaces result;
  if (vertices.empty())
    return result;
  Vector3 low = vertices.front();
  Vector3 high = vertices.front();
  for (const auto& v : vertices)
  {
    low = Vector3(std::min(low.x, v.x), std::min(low.y, v.y), std::min(low.z, v.z));
    high = Vector3(std::max(high.x, v.x), std::max(high.y, v.y), std::max(high.z, v.z));
  }
  const double tolerance = coincidence_tolerance * (high - low).Norm();
  if (not(tolerance > 0.0))
    return result;

  // Unconnected faces, grouped by line or plane: the normal with its first significant component
  // made positive, and the offset of the plane along it.
  struct Face
  {
    uint64_t cell = 0;
    const UnpartitionedMesh::LightWeightFace* face = nullptr;
    bool is_polygon = false;
    int orientation = 1;
    Vector3 center;
    Vector3 normal;
  };
  std::map<std::array<int64_t, 4>, std::vector<Face>> planes;
  for (uint64_t c = 0; c < cells.size(); ++c)
  {
    const auto& cell = *cells[c];
    if (cell.type != CellType::POLYGON and cell.type != CellType::POLYHEDRON)
      continue;
    for (const auto& face : cell.faces)
    {
      if (face.has_neighbor)
        continue;
      const Vector3 area = FaceAreaVector(cell.type, face, vertices);
      const double magnitude = area.Norm();
      if (not(magnitude > 0.0))
        continue;
      Vector3 normal = area / magnitude;
      int orientation = 1;
      for (const double component : {normal.x, normal.y, normal.z})
        if (std::fabs(component) > 1.0e-6)
        {
          orientation = component > 0.0 ? 1 : -1;
          break;
        }
      normal = normal * static_cast<double>(orientation);
      const Vector3 center = FaceCenter(face, vertices);
      // The normal is rounded to bins of 1e-6, so that an axis normal and its rounded copy share a
      // bin. The offset is binned with floor, and the search also looks in the next bins.
      const auto Round = [](const double value) { return std::llround(value / 1.0e-6); };
      planes[{Round(normal.x),
              Round(normal.y),
              Round(normal.z),
              static_cast<int64_t>(std::floor(normal.Dot(center) / (4.0 * tolerance)))}]
        .push_back({c, &face, cell.type == CellType::POLYGON, orientation, center, normal});
    }
  }

  // Two faces of one plane, facing opposite ways, that overlap in more than the tolerance.
  const auto Overlap = [&](const Face& a, const Face& b)
  {
    if (std::fabs(a.normal.Dot(a.center) - a.normal.Dot(b.center)) > tolerance)
      return false;
    // A basis of the plane. For a polygon's edge, which lies in the x-y plane, the edge direction.
    Vector3 t1 = a.normal.Cross(Vector3(0.0, 0.0, 1.0));
    if (t1.Norm() < 0.5)
      t1 = a.normal.Cross(Vector3(1.0, 0.0, 0.0));
    t1 = t1 / t1.Norm();
    const Vector3 t2 = a.normal.Cross(t1);
    const auto Project = [&](const Face& f)
    {
      std::vector<std::array<double, 2>> points;
      for (const auto vid : f.face->vertex_ids)
      {
        const Vector3 d = vertices[vid] - a.center;
        points.push_back({d.Dot(t1), d.Dot(t2)});
      }
      return points;
    };
    const auto pa = Project(a);
    const auto pb = Project(b);
    if (a.is_polygon)
    {
      // Edges: overlapping intervals along the edge direction.
      const double a0 = std::min(pa[0][0], pa[1][0]);
      const double a1 = std::max(pa[0][0], pa[1][0]);
      const double b0 = std::min(pb[0][0], pb[1][0]);
      const double b1 = std::max(pb[0][0], pb[1][0]);
      return std::min(a1, b1) - std::max(a0, b0) > tolerance;
    }
    // Polygons: a vertex or the center of one strictly inside the other.
    const auto Center = [](const std::vector<std::array<double, 2>>& points)
    {
      std::array<double, 2> center = {0.0, 0.0};
      for (const auto& p : points)
        for (int i = 0; i < 2; ++i)
          center[i] += p[i] / static_cast<double>(points.size());
      return center;
    };
    if (StrictlyInsidePolygon(pb, Center(pa), tolerance) or
        StrictlyInsidePolygon(pa, Center(pb), tolerance))
      return true;
    for (const auto& p : pa)
      if (StrictlyInsidePolygon(pb, p, tolerance))
        return true;
    for (const auto& p : pb)
      if (StrictlyInsidePolygon(pa, p, tolerance))
        return true;
    return false;
  };

  // Each face against the opposite faces of its offset bin and of the two next to it, which can
  // lie within the tolerance.
  for (const auto& [key, faces] : planes)
    for (const auto& a : faces)
    {
      bool hanging = false;
      for (const int64_t shift : {-1, 0, 1})
      {
        const auto other = planes.find({key[0], key[1], key[2], key[3] + shift});
        if (other == planes.end())
          continue;
        for (const auto& b : other->second)
          if (b.cell != a.cell and b.orientation != a.orientation and Overlap(a, b))
          {
            hanging = true;
            break;
          }
        if (hanging)
          break;
      }
      if (not hanging)
        continue;
      if (result.count == 0)
      {
        result.first_cell = a.cell;
        result.first_center = a.center;
      }
      ++result.count;
    }
  return result;
}

} // namespace

UnpartitionedMesh::UnpartitionedMesh()
  : dim_(0), coord_sys_(CoordinateSystemType::CARTESIAN), mesh_type_(UNSTRUCTURED), extruded_(false)
{
}

UnpartitionedMesh::~UnpartitionedMesh()
{
  log.Log0Verbose1() << "Unpartitioned Mesh destructor called";
}

void
UnpartitionedMesh::ComputeBoundingBox()
{
  for (size_t p = 0; p < vertices_.size(); ++p)
  {
    const auto& vertex = vertices_[p];
    if (p == 0)
      bound_box_ = {vertex.x, vertex.x, vertex.y, vertex.y, vertex.z, vertex.z};
    else
    {
      bound_box_.xmin = std::min(bound_box_.xmin, vertex.x);
      bound_box_.xmax = std::max(bound_box_.xmax, vertex.x);
      bound_box_.ymin = std::min(bound_box_.ymin, vertex.y);
      bound_box_.ymax = std::max(bound_box_.ymax, vertex.y);
      bound_box_.zmin = std::min(bound_box_.zmin, vertex.z);
      bound_box_.zmax = std::max(bound_box_.zmax, vertex.z);
    }
  }
}

void
UnpartitionedMesh::ComputeCentroids()
{
  log.Log0Verbose1() << "Computing cell-centroids.";
  for (auto& cell : raw_cells_)
  {
    cell->centroid = Vector3(0.0, 0.0, 0.0);
    for (auto vid : cell->vertex_ids)
      cell->centroid += vertices_[vid];

    cell->centroid = cell->centroid / static_cast<double>(cell->vertex_ids.size());
  }
  log.Log0Verbose1() << "Done computing cell-centroids.";
}

void
UnpartitionedMesh::CheckQuality()
{
  log.Log0Verbose1() << "Checking cell-center-to-face orientations";
  const Vector3 khat(0.0, 0.0, 1.0);
  size_t num_negative_volume_elements = 0;
  for (const auto& cell : raw_cells_)
  {
    if (cell->type == CellType::POLYGON)
    {
      // Form triangles
      size_t num_verts = cell->vertex_ids.size();
      for (size_t v = 0; v < num_verts; ++v)
      {
        size_t vp1 = (v < (num_verts - 1)) ? v + 1 : 0;

        const auto& v0 = vertices_[cell->vertex_ids[v]];
        const auto& v1 = vertices_[cell->vertex_ids[vp1]];

        auto E01 = v1 - v0;
        auto n = E01.Cross(khat).Normalized();

        if (n.Dot(v0 - cell->centroid) < 0.0)
          ++num_negative_volume_elements;
      } // for v
    }
    else if (cell->type == CellType::POLYHEDRON)
    {
      for (auto& face : cell->faces)
      {
        if (face.vertex_ids.size() < 2)
          throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                 ": cell-center-to-face check encountered face "
                                 "with less than 2 vertices on a face, making "
                                 "normal computation impossible.");

        // Compute centroid
        Vector3 face_centroid;
        for (uint64_t vid : face.vertex_ids)
          face_centroid += vertices_[vid];
        face_centroid /= static_cast<double>(face.vertex_ids.size());

        // Form tets for each face edge
        size_t num_face_verts = face.vertex_ids.size();
        for (size_t fv = 0; fv < face.vertex_ids.size(); ++fv)
        {
          size_t fvp1 = (fv < (num_face_verts - 1)) ? fv + 1 : 0;

          const auto& fv1 = vertices_[face.vertex_ids[fv]];
          const auto& fv2 = vertices_[face.vertex_ids[fvp1]];

          auto E0 = fv1 - face_centroid;
          auto E1 = fv2 - face_centroid;
          auto n = E0.Cross(E1).Normalized();

          if (n.Dot(fv1 - cell->centroid) < 0.0)
            ++num_negative_volume_elements;
        }
      } // for face
    } // if polyhedron
  } // for cell in raw_cells

  log.Log0Verbose1() << "Checking face sizes";
  size_t cell_id = 0;
  for (const auto& cell : raw_cells_)
  {
    if (cell->type == CellType::POLYGON)
    {
      size_t f = 0;
      for (const auto& face : cell->faces)
      {
        const auto& v0 = vertices_.at(face.vertex_ids[0]);
        const auto& v1 = vertices_.at(face.vertex_ids[1]);
        OpenSnLogicalErrorIf((v1 - v0).Norm() < 1.0e-12,
                             "Cell " + std::to_string(cell_id) +
                               " (centroid=" + cell->centroid.PrintStr() + ") face " +
                               std::to_string(f) + ": Face has length < 1.0e-12.");
        ++f;
      }
    } // if polygon
    if (cell->type == CellType::POLYHEDRON)
    {
      size_t f = 0;
      for (const auto& face : cell->faces)
      {
        size_t num_face_verts = face.vertex_ids.size();
        for (size_t s = 0; s < face.vertex_ids.size(); ++s)
        {
          size_t fvp1 = (s < (num_face_verts - 1)) ? s + 1 : 0;

          const auto& v0 = vertices_.at(face.vertex_ids[s]);
          const auto& v1 = vertices_.at(face.vertex_ids[fvp1]);

          OpenSnLogicalErrorIf((v1 - v0).Norm() < 1.0e-12,
                               "Cell " + std::to_string(cell_id) + " (centroid=" +
                                 cell->centroid.PrintStr() + ") face " + std::to_string(f) +
                                 " side " + std::to_string(s) + ": Side has length < 1.0e-12.");
        }

        ++f;
      }
    } // if polyhedron
    ++cell_id;
  } // for cell in raw_cells

  if (num_negative_volume_elements > 0)
    log.LogAllWarning() << "Cell quality checks detected " << num_negative_volume_elements
                        << " negative volume sub-elements (sub-triangle or sub-tetrahedron)."
                        << " This issue could result in incorrect quantities"
                        << " under some circumstances.";
  log.Log() << "Done checking cell-center-to-face orientations";

  // An inverted cell has its faces ordered inward, so its face normals point into it. A degenerate
  // cell has no volume. Neither can be transported through.
  log.Log0Verbose1() << "Checking cell volumes";
  size_t num_inverted = 0;
  size_t num_degenerate = 0;
  size_t first_bad_cell = 0;
  double first_bad_measure = 0.0;
  for (size_t c = 0; c < raw_cells_.size(); ++c)
  {
    const auto& cell = *raw_cells_[c];
    if (cell.type != CellType::POLYGON and cell.type != CellType::POLYHEDRON)
      continue;
    const auto [measure, scale] = SignedCellMeasure(cell, vertices_);
    const bool degenerate = std::fabs(measure) <= degenerate_measure_tolerance * scale;
    if (not degenerate and measure > 0.0)
      continue;
    if (num_inverted + num_degenerate == 0)
    {
      first_bad_cell = c;
      first_bad_measure = measure;
    }
    if (degenerate)
      ++num_degenerate;
    else
      ++num_inverted;
  }
  if (num_inverted + num_degenerate > 0)
  {
    std::ostringstream message;
    message << "The mesh has " << num_inverted << " inverted cells and " << num_degenerate
            << " cells of zero volume. The first is cell " << first_bad_cell << " at centroid "
            << raw_cells_[first_bad_cell]->centroid.PrintStr() << ", with volume "
            << std::scientific << first_bad_measure
            << ". A cell is inverted when its faces are ordered inward, for example a polygon "
               "whose vertices are ordered clockwise.";
    OpenSnLogicalError(message.str());
  }
}

uint64_t
UnpartitionedMesh::MakeBoundaryID(const std::string& boundary_name)
{
  if (boundary_id_map_.empty())
    return 0;

  for (const auto& [id, name] : boundary_id_map_)
    if (boundary_name == name)
      return id;

  uint64_t max_id = 0;
  for (const auto& [id, name] : boundary_id_map_)
    max_id = std::max(id, max_id);

  return max_id + 1;
}

void
UnpartitionedMesh::AddBoundary(uint64_t id, const std::string& name)
{
  boundary_id_map_[id] = name;
  boundary_name_map_[name] = id;
}

void
UnpartitionedMesh::SetOrthoAttributes(size_t nx, size_t ny, size_t nz)
{
  mesh_type_ = ORTHOGONAL;
  ortho_attrs_.Nx = nx;
  ortho_attrs_.Ny = ny;
  ortho_attrs_.Nz = nz;
}

void
UnpartitionedMesh::BuildMeshConnectivity()
{
  const size_t num_raw_cells = raw_cells_.size();
  const size_t num_raw_vertices = vertices_.size();

  // Reset all cell neighbors
  int num_bndry_faces = 0;
  for (auto& cell : raw_cells_)
    for (auto& face : cell->faces)
      if (not face.has_neighbor)
        ++num_bndry_faces;

  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Number of unconnected faces before connectivity: " << num_bndry_faces;

  log.Log() << program_timer.GetTimeString() << " Establishing cell connectivity.";

  // Establish internal connectivity
  // Populate vertex subscriptions to internal cells
  vertex_cell_subscriptions_.resize(num_raw_vertices);
  {
    uint64_t cur_cell_id = 0;
    for (const auto& cell : raw_cells_)
    {
      for (auto vid : cell->vertex_ids)
        vertex_cell_subscriptions_.at(vid).insert(cur_cell_id);
      ++cur_cell_id;
    }
  }

  log.Log() << program_timer.GetTimeString() << " Vertex cell subscriptions complete.";

  // Process raw cells. The two cells of an interior face must lie on opposite sides of it, so their
  // copies of the face have opposite area vectors. When they have the same direction, the cells lie
  // on the same side and overlap: the mesh is folded.
  //
  // A face must also belong to at most two cells. Cells are processed in order of ID and their
  // candidate neighbors are searched in order of ID, so the cell with the lowest ID among those
  // sharing a face connects it first. Every later cell with the same face then finds that copy
  // already connected before any other unconnected copy, whether or not the cells after it pair up
  // among themselves.
  // Copies of each folded face, which are disconnected again before the error is raised.
  std::vector<std::pair<LightWeightFace*, LightWeightFace*>> folded_faces;
  uint64_t first_folded_cell = 0;
  uint64_t first_folded_neighbor = 0;
  Vector3 first_folded_center;
  size_t num_shared_faces = 0;
  uint64_t first_shared_cell = 0;
  Vector3 first_shared_center;
  {
    uint64_t aux_counter = 0;
    uint64_t cur_cell_id = 0;
    for (auto& cell : raw_cells_)
    {
      for (auto& cur_cell_face : cell->faces)
      {
        if (cur_cell_face.has_neighbor)
          continue;
        const std::set<uint64_t> cfvids(cur_cell_face.vertex_ids.begin(),
                                        cur_cell_face.vertex_ids.end());

        std::set<size_t> cells_to_search;
        for (uint64_t vid : cfvids)
          for (uint64_t cell_id : vertex_cell_subscriptions_.at(vid))
            if (cell_id != cur_cell_id)
              cells_to_search.insert(cell_id);

        for (uint64_t adj_cell_id : cells_to_search)
        {
          auto adj_cell = raw_cells_.at(adj_cell_id);

          for (auto& adj_cell_face : adj_cell->faces)
          {
            // Every vertex of an equal face is a vertex of the current face. Testing that first
            // rejects most faces without building their vertex sets.
            const auto& adj_vids = adj_cell_face.vertex_ids;
            if (not std::all_of(adj_vids.begin(),
                                adj_vids.end(),
                                [&cfvids](const uint64_t vid) { return cfvids.contains(vid); }))
              continue;
            const std::set<uint64_t> afvids(adj_vids.begin(), adj_vids.end());
            if (cfvids != afvids)
              continue;

            // The face already connects two other cells.
            if (adj_cell_face.has_neighbor)
            {
              if (num_shared_faces == 0)
              {
                first_shared_cell = cur_cell_id;
                first_shared_center = FaceCenter(cur_cell_face, vertices_);
              }
              ++num_shared_faces;
              goto face_neighbor_found;
            }

            cur_cell_face.neighbor = adj_cell_id;
            adj_cell_face.neighbor = cur_cell_id;

            cur_cell_face.has_neighbor = true;
            adj_cell_face.has_neighbor = true;

            if ((cell->type == CellType::POLYGON or cell->type == CellType::POLYHEDRON) and
                FaceAreaVector(cell->type, cur_cell_face, vertices_)
                    .Dot(FaceAreaVector(adj_cell->type, adj_cell_face, vertices_)) >= 0.0)
            {
              if (folded_faces.empty())
              {
                first_folded_cell = cur_cell_id;
                first_folded_neighbor = adj_cell_id;
                first_folded_center = FaceCenter(cur_cell_face, vertices_);
              }
              folded_faces.emplace_back(&cur_cell_face, &adj_cell_face);
            }

            goto face_neighbor_found;
          } // for adjacent cell face
        }
      face_neighbor_found:;
      } // for face

      ++cur_cell_id;
      const double fraction_complete =
        static_cast<double>(cur_cell_id) / static_cast<double>(num_raw_cells);
      if (fraction_complete >= static_cast<double>(aux_counter + 1) * 0.1)
      {
        log.Log() << program_timer.GetTimeString() << " Surpassing cell " << cur_cell_id << " of "
                  << num_raw_cells << " (" << (aux_counter + 1) * 10 << "%)";
        ++aux_counter;
      }
    } // for cell
  }

  if (not folded_faces.empty())
  {
    // Connected faces are skipped when the connectivity is built again, so a folded pair left
    // connected would then be accepted. Faces of more than two cells need no reset, because their
    // extra copies were never connected.
    for (const auto& [face, adj_face] : folded_faces)
      for (auto* folded_face : {face, adj_face})
      {
        folded_face->has_neighbor = false;
        folded_face->neighbor = LightWeightFace().neighbor;
      }
    OpenSnLogicalError(
      "The mesh is folded: " + std::to_string(folded_faces.size()) +
      " interior faces have both of their cells on the same side, so the cells overlap. The first "
      "is the face between cells " +
      std::to_string(first_folded_cell) + " and " + std::to_string(first_folded_neighbor) + " at " +
      first_folded_center.PrintStr() + ".");
  }

  if (num_shared_faces > 0)
    OpenSnLogicalError("The mesh is not conforming: " + std::to_string(num_shared_faces) +
                       " cell faces coincide with a face that two other cells already share. The "
                       "first is a face of cell " +
                       std::to_string(first_shared_cell) + " at " + first_shared_center.PrintStr() +
                       ".");

  if (const auto hanging = FindHangingFaces(raw_cells_, vertices_); hanging.count > 0)
    OpenSnLogicalError(
      "The mesh is not conforming: " + std::to_string(hanging.count) +
      " cell faces overlap a face of another cell without sharing its vertices (a hanging node "
      "or T-junction). The two cells are then not connected, and both faces would become "
      "boundaries inside the domain. The first is a face of cell " +
      std::to_string(hanging.first_cell) + " at " + hanging.first_center.PrintStr() + ".");

  log.Log() << program_timer.GetTimeString() << " Establishing cell boundary connectivity.";

  // Establish boundary connectivity
  // Make list of internal cells on the boundary
  std::vector<std::shared_ptr<LightWeightCell>> internal_cells_on_boundary;
  for (auto& cell : raw_cells_)
  {
    bool cell_on_boundary = false;
    for (auto& face : cell->faces)
      if (not face.has_neighbor)
      {
        cell_on_boundary = true;
        break;
      }

    if (cell_on_boundary)
      internal_cells_on_boundary.push_back(cell);
  }

  // Populate vertex subscriptions to boundary cells
  std::vector<std::set<uint64_t>> vertex_bndry_cell_subscriptions(vertices_.size());
  {
    uint64_t cur_cell_id = 0;
    for (auto& cell : raw_boundary_cells_)
    {
      for (auto vid : cell->vertex_ids)
        vertex_bndry_cell_subscriptions.at(vid).insert(cur_cell_id);
      ++cur_cell_id;
    }
  }

  // Process boundary cells
  for (auto& cell : internal_cells_on_boundary)
    for (auto& face : cell->faces)
    {
      if (face.has_neighbor)
        continue;
      std::set<uint64_t> cfvids(face.vertex_ids.begin(), face.vertex_ids.end());

      std::set<size_t> cells_to_search;
      for (uint64_t vid : face.vertex_ids)
        for (uint64_t cell_id : vertex_bndry_cell_subscriptions[vid])
          cells_to_search.insert(cell_id);

      for (uint64_t adj_cell_id : cells_to_search)
      {
        auto& adj_cell = raw_boundary_cells_[adj_cell_id];

        std::set<uint64_t> afvids(adj_cell->vertex_ids.begin(), adj_cell->vertex_ids.end());

        if (cfvids == afvids)
        {
          face.neighbor = adj_cell_id;
          break;
        }
      } // for adj_cell_id
    } // for face

  num_bndry_faces = 0;
  for (const auto& cell : raw_cells_)
    for (auto& face : cell->faces)
      if (not face.has_neighbor)
        ++num_bndry_faces;

  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Number of boundary faces "
                        "after connectivity: "
                     << num_bndry_faces;

  log.Log() << program_timer.GetTimeString() << " Done establishing cell connectivity.";
}

} // namespace opensn
