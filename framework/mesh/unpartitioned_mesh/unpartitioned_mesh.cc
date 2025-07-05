// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/mesh/cell/cell.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include <algorithm>

namespace opensn
{

UnpartitionedMesh::UnpartitionedMesh() : dim_(0), mesh_type_(UNSTRUCTURED), extruded_(false)
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
                     << " Number of unconnected faces "
                        "before connectivity: "
                     << num_bndry_faces;

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

  // Process raw cells
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
            if (adj_cell_face.has_neighbor)
              continue;
            const std::set<uint64_t> afvids(adj_cell_face.vertex_ids.begin(),
                                            adj_cell_face.vertex_ids.end());

            if (cfvids == afvids)
            {
              cur_cell_face.neighbor = adj_cell_id;
              adj_cell_face.neighbor = cur_cell_id;

              cur_cell_face.has_neighbor = true;
              adj_cell_face.has_neighbor = true;

              goto face_neighbor_found;
            }
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
          face.neighbor = adj_cell->block_id;
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
