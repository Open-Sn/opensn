// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/utils/error.h"
#include <algorithm>
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

std::span<const uint64_t>
UnpartitionedMesh::GetCellConnectivity(std::uint32_t cell_global_id) const
{
  if (cell_global_id < connect_ofst_.size() - 1)
  {
    auto first = connect_ofst_[cell_global_id];
    auto last = connect_ofst_[cell_global_id + 1];
    return std::span{connect_ids_.data() + first, connect_ids_.data() + last};
  }
  else
    throw std::out_of_range("Cell local id out of range");
}

std::span<const CellFace>
UnpartitionedMesh::GetCellFaces(std::uint32_t cell_global_id) const
{
  if (cell_global_id < connect_ofst_.size() - 1)
  {
    auto first = face_connect_ofst_[cell_global_id];
    auto last = face_connect_ofst_[cell_global_id + 1];
    return std::span{faces_.data() + first, faces_.data() + last};
  }
  else
    throw std::out_of_range("Cell global id out of range");
}

std::span<CellFace>
UnpartitionedMesh::GetCellFaces(std::uint32_t cell_global_id)
{
  if (cell_global_id < connect_ofst_.size() - 1)
  {
    auto first = face_connect_ofst_[cell_global_id];
    auto last = face_connect_ofst_[cell_global_id + 1];
    return std::span{faces_.data() + first, faces_.data() + last};
  }
  else
    throw std::out_of_range("Cell global id out of range");
}

std::uint64_t
UnpartitionedMesh::GetCellFaceCount(std::uint32_t cell_global_id) const
{
  return cell_face_connectivity_[cell_global_id].size();
}

const CellFace&
UnpartitionedMesh::GetCellFace(std::uint32_t cell_global_id, std::uint32_t face_idx) const
{
  if (cell_global_id < face_connect_ofst_.size() - 1)
  {
    auto ofst = face_connect_ofst_[cell_global_id] + face_idx;
    return faces_[ofst];
  }
  else
    throw std::out_of_range("Cell local id out of range");
}

CellFace&
UnpartitionedMesh::GetCellFace(std::uint32_t cell_global_id, std::uint32_t face_idx)
{
  if (cell_global_id < face_connect_ofst_.size() - 1)
  {
    auto ofst = face_connect_ofst_[cell_global_id] + face_idx;
    return faces_[ofst];
  }
  else
    throw std::out_of_range("Cell local id out of range");
}

void
UnpartitionedMesh::ComputeCentroids()
{
  log.Log0Verbose1() << "Computing cell-centroids.";
  for (std::size_t cell_id = 0; cell_id < cells_.size(); ++cell_id)
  {
    auto cell_vertex_ids = GetCellConnectivity(cell_id);
    auto& cell = cells_[cell_id];
    cell.centroid = Vector3(0.0, 0.0, 0.0);
    for (auto vid : cell_vertex_ids)
      cell.centroid += vertices_[vid];

    cell.centroid = cell.centroid / static_cast<double>(cell_vertex_ids.size());
  }
  log.Log0Verbose1() << "Done computing cell-centroids.";
}

void
UnpartitionedMesh::CheckQuality()
{
  log.Log0Verbose1() << "Checking cell-center-to-face orientations";
  const Vector3 khat(0.0, 0.0, 1.0);
  size_t num_negative_volume_elements = 0;
  for (std::size_t cell_id = 0; cell_id < cells_.size(); ++cell_id)
  {
    auto& cell = cells_[cell_id];
    if (cell.GetType() == CellType::POLYGON)
    {
      auto cell_vertex_ids = GetCellConnectivity(cell_id);
      // Form triangles
      size_t num_verts = cell_vertex_ids.size();
      for (size_t v = 0; v < num_verts; ++v)
      {
        size_t vp1 = (v < (num_verts - 1)) ? v + 1 : 0;

        const auto& v0 = vertices_[cell_vertex_ids[v]];
        const auto& v1 = vertices_[cell_vertex_ids[vp1]];

        auto E01 = v1 - v0;
        auto n = E01.Cross(khat).Normalized();

        if (n.Dot(v0 - cell.centroid) < 0.0)
          ++num_negative_volume_elements;
      } // for v
    }
    else if (cell.GetType() == CellType::POLYHEDRON)
    {
      for (size_t f = 0; f < cell_face_connectivity_[cell_id].size(); ++f)
      {
        const auto& face_vertex_ids = cell_face_connectivity_[cell_id][f];
        if (face_vertex_ids.size() < 2)
          throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                                 ": cell-center-to-face check encountered face "
                                 "with less than 2 vertices on a face, making "
                                 "normal computation impossible.");

        // Compute centroid
        Vector3 face_centroid;
        for (uint64_t vid : face_vertex_ids)
          face_centroid += vertices_[vid];
        face_centroid /= static_cast<double>(face_vertex_ids.size());

        // Form tets for each face edge
        size_t num_face_verts = face_vertex_ids.size();
        for (size_t fv = 0; fv < face_vertex_ids.size(); ++fv)
        {
          size_t fvp1 = (fv < (num_face_verts - 1)) ? fv + 1 : 0;

          const auto& fv1 = vertices_[face_vertex_ids[fv]];
          const auto& fv2 = vertices_[face_vertex_ids[fvp1]];

          auto E0 = fv1 - face_centroid;
          auto E1 = fv2 - face_centroid;
          auto n = E0.Cross(E1).Normalized();

          if (n.Dot(fv1 - cell.centroid) < 0.0)
            ++num_negative_volume_elements;
        }
      } // for face
    } // if polyhedron
  } // for cell in raw_cells

  log.Log0Verbose1() << "Checking face sizes";
  size_t cell_id = 0;
  for (const auto& cell : cells_)
  {
    if (cell.GetType() == CellType::POLYGON)
    {
      for (size_t f = 0; f < cell_face_connectivity_[cell_id].size(); ++f)
      {
        const auto& face_vertex_ids = cell_face_connectivity_[cell_id][f];
        const auto& v0 = vertices_.at(face_vertex_ids[0]);
        const auto& v1 = vertices_.at(face_vertex_ids[1]);
        OpenSnLogicalErrorIf((v1 - v0).Norm() < 1.0e-12,
                             "Cell " + std::to_string(cell_id) +
                               " (centroid=" + cell.centroid.PrintStr() + ") face " +
                               std::to_string(f) + ": Face has length < 1.0e-12.");
      }
    } // if polygon
    if (cell.GetType() == CellType::POLYHEDRON)
    {
      for (size_t f = 0; f < cell_face_connectivity_[cell_id].size(); ++f)
      {
        const auto& face_vertex_ids = cell_face_connectivity_[cell_id][f];
        size_t num_face_verts = face_vertex_ids.size();
        for (size_t s = 0; s < face_vertex_ids.size(); ++s)
        {
          size_t fvp1 = (s < (num_face_verts - 1)) ? s + 1 : 0;

          const auto& v0 = vertices_.at(face_vertex_ids[s]);
          const auto& v1 = vertices_.at(face_vertex_ids[fvp1]);

          OpenSnLogicalErrorIf((v1 - v0).Norm() < 1.0e-12,
                               "Cell " + std::to_string(cell_id) + " (centroid=" +
                                 cell.centroid.PrintStr() + ") face " + std::to_string(f) +
                                 " side " + std::to_string(s) + ": Side has length < 1.0e-12.");
        }
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
UnpartitionedMesh::SetCells(std::vector<Cell>&& cells,
                            const std::vector<std::vector<std::uint64_t>>& cell_connectivity)
{
  cells_ = std::move(cells);
  // Assign global numbering. Otherwise partitioning fails, since it now requires valid
  // Cell::global_id for indexing into cell faces.
  std::uint64_t gid = 0;
  for (auto& cell : cells_)
    cell.global_id = gid++;

  std::size_t total_len = 0;
  for (const auto& cell : cell_connectivity)
    total_len += cell.size();
  connect_ids_.reserve(total_len);

  connect_ofst_.reserve(cell_connectivity.size() + 1);
  for (const auto& cell : cell_connectivity)
  {
    connect_ofst_.push_back(connect_ids_.size());
    for (const auto& id : cell)
      connect_ids_.push_back(id);
  }
  connect_ofst_.push_back(total_len);
}

void
UnpartitionedMesh::SetCellFaces(
  std::vector<CellFace>&& faces,
  const std::vector<std::vector<std::vector<std::uint64_t>>>& cell_face_connectivity)
{
  faces_ = std::move(faces);
  face_connect_ofst_.push_back(0);
  for (const auto& cell : cells_)
  {
    try
    {
      const auto n_faces = cell_face_connectivity.at(cell.global_id).size();
      face_connect_ofst_.push_back(face_connect_ofst_.back() + n_faces);
    }
    catch (const std::out_of_range& e)
    {
      throw std::out_of_range("map::at: key not found for local cell.global_id = " +
                              std::to_string(cell.global_id));
    }
  }
  cell_face_connectivity_ = cell_face_connectivity;
}

void
UnpartitionedMesh::BuildMeshConnectivity()
{
  const size_t num_raw_cells = cells_.size();
  const size_t num_raw_vertices = vertices_.size();

  // Reset all cell neighbors
  int num_bndry_faces = 0;
  for (auto& face : faces_)
    if (not face.has_neighbor)
      ++num_bndry_faces;

  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Number of unconnected faces before connectivity: " << num_bndry_faces;

  log.Log() << program_timer.GetTimeString() << " Establishing cell connectivity.";

  // Establish internal connectivity
  // Populate vertex subscriptions to internal cells
  vertex_cell_subscriptions_.resize(num_raw_vertices);
  {
    for (std::size_t cell_id = 0; cell_id < cells_.size(); ++cell_id)
    {
      auto cell_vertex_ids = GetCellConnectivity(cell_id);
      for (auto vid : cell_vertex_ids)
        vertex_cell_subscriptions_.at(vid).insert(cell_id);
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
    for (auto& cell : cells_)
    {
      auto cell_faces = GetCellFaces(cur_cell_id);
      for (size_t f = 0; f < cell_face_connectivity_[cur_cell_id].size(); ++f)
      {
        auto& cur_cell_face = cell_faces[f];
        if (cur_cell_face.has_neighbor)
          continue;
        const auto& cfvids_vec = cell_face_connectivity_[cur_cell_id][f];
        const std::set<uint64_t> cfvids(cfvids_vec.begin(), cfvids_vec.end());

        std::set<size_t> cells_to_search;
        for (uint64_t vid : cfvids)
          for (uint64_t cell_id : vertex_cell_subscriptions_.at(vid))
            if (cell_id != cur_cell_id)
              cells_to_search.insert(cell_id);

        for (uint64_t adj_cell_id : cells_to_search)
        {
          auto& adj_cell = cells_.at(adj_cell_id);

          for (size_t af = 0; af < cell_face_connectivity_[adj_cell_id].size(); ++af)
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

  log.Log() << program_timer.GetTimeString() << " Establishing cell boundary connectivity.";

  // Establish boundary connectivity
  // Make list of internal cells on the boundary
  std::vector<Cell*> internal_cells_on_boundary;
  for (auto& cell : cells_)
  {
    bool cell_on_boundary = false;
    auto cell_faces = GetCellFaces(cell.global_id);
    for (auto& face : cell_faces)
      if (not face.has_neighbor)
      {
        cell_on_boundary = true;
        break;
      }

    if (cell_on_boundary)
      internal_cells_on_boundary.push_back(&cell);
  }

  num_bndry_faces = 0;
  for (const auto& face : faces_)
    if (not face.has_neighbor)
      ++num_bndry_faces;

  log.Log0Verbose1() << program_timer.GetTimeString()
                     << " Number of boundary faces "
                        "after connectivity: "
                     << num_bndry_faces;

  log.Log() << program_timer.GetTimeString() << " Done establishing cell connectivity.";
}

} // namespace opensn
