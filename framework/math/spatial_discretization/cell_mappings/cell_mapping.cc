// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include <utility>

namespace opensn
{

CellMapping::CellMapping(const MeshContinuum& grid,
                         const Cell& cell,
                         size_t num_nodes,
                         std::vector<Vector3> node_locations,
                         std::vector<std::vector<int>> face_node_mappings,
                         const VandAFunction& volume_area_function)
  : ref_grid_(grid),
    cell_(cell),
    num_nodes_(num_nodes),
    node_locations_(std::move(node_locations)),
    face_node_mappings_(std::move(face_node_mappings))
{
  volume_area_function(ref_grid_, cell, volume_, areas_);
}

const Cell&
CellMapping::ReferenceCell() const
{
  return cell_;
}

const MeshContinuum&
CellMapping::ReferenceGrid() const
{
  return ref_grid_;
}

size_t
CellMapping::NumNodes() const
{
  return num_nodes_;
}

size_t
CellMapping::NumFaceNodes(size_t face_index) const
{
  return face_node_mappings_.at(face_index).size();
}

const std::vector<std::vector<int>>&
CellMapping::GetFaceNodeMappings() const
{
  return face_node_mappings_;
}

double
CellMapping::CellVolume() const
{
  return volume_;
}

double
CellMapping::FaceArea(size_t face_index) const
{
  return areas_[face_index];
}

int
CellMapping::MapFaceNode(size_t face_index, size_t face_node_index) const
{
  try
  {
    return face_node_mappings_.at(face_index).at(face_node_index);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range("chi_math::CellMapping::MapFaceNode: "
                            "Either face_index or face_node_index is out of range");
  }
}

void
CellMapping::ComputeCellVolumeAndAreas(const MeshContinuum& grid,
                                       const Cell& cell,
                                       double& volume,
                                       std::vector<double>& areas)
{
  switch (cell.Type())
  {
    case CellType::SLAB:
    {
      const auto& v0 = grid.vertices[cell.vertex_ids_[0]];
      const auto& v1 = grid.vertices[cell.vertex_ids_[1]];

      volume = (v1 - v0).Norm();
      areas = {1.0, 1.0};
      break;
    }
    case CellType::POLYGON:
    {
      volume = 0.0;
      const auto& v2 = cell.centroid_;

      size_t num_faces = cell.faces_.size();
      areas.reserve(num_faces);

      for (size_t f = 0; f < num_faces; ++f)
      {
        const uint64_t v0i = cell.faces_[f].vertex_ids_[0];
        const uint64_t v1i = cell.faces_[f].vertex_ids_[1];

        const auto& v0 = grid.vertices[v0i];
        const auto& v1 = grid.vertices[v1i];

        areas.push_back((v1 - v0).Norm());

        const Vector3 sidev01 = v1 - v0;
        const Vector3 sidev02 = v2 - v0;

        double sidedetJ = ((sidev01.x) * (sidev02.y) - (sidev02.x) * (sidev01.y));

        volume += sidedetJ / 2.0;
      } // for face

      break;
    }
    case CellType::POLYHEDRON:
    {
      volume = 0.0;
      const auto& vcc = cell.centroid_;

      size_t num_faces = cell.faces_.size();
      areas.assign(num_faces, 0.0);
      for (size_t f = 0; f < num_faces; f++)
      {
        const auto& face = cell.faces_[f];
        const size_t num_edges = face.vertex_ids_.size();
        for (size_t e = 0; e < num_edges; ++e)
        {
          size_t ep1 = (e < (num_edges - 1)) ? e + 1 : 0;
          uint64_t v0i = face.vertex_ids_[e];
          uint64_t v1i = face.vertex_ids_[ep1];

          const auto& v0 = grid.vertices[v0i];
          const auto& v1 = cell.faces_[f].centroid_;
          const auto& v2 = grid.vertices[v1i];
          const auto& v3 = vcc;

          const auto sidev01 = v1 - v0;
          const auto sidev02 = v2 - v0;
          const auto sidev03 = v3 - v0;

          Matrix3x3 J;

          J.SetColJVec(0, sidev01);
          J.SetColJVec(1, sidev02);
          J.SetColJVec(2, sidev03);

          areas[f] += (sidev01.Cross(sidev02)).Norm() / 2.0;
          volume += J.Det() / 6.0;
        } // for edge
      }   // for face
      break;
    }
    default:
      throw std::logic_error("chi_math::CellMapping::ComputeCellVolume: "
                             "Unsupported cell type.");
  }
}

const std::vector<Vector3>&
CellMapping::GetNodeLocations() const
{
  return node_locations_;
}

} // namespace opensn
