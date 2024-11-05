// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/petsc_utils/petsc_utils.h"
#include "framework/logging/log.h"

namespace opensn
{

SpatialDiscretization::SpatialDiscretization(const MeshContinuum& grid,
                                             CoordinateSystemType cs_type,
                                             SpatialDiscretizationType sdm_type)
  : UNITARY_UNKNOWN_MANAGER({std::make_pair(UnknownType::SCALAR, 0)}),
    ref_grid_(grid),
    coord_sys_type_(cs_type),
    type_(sdm_type)
{
}

const CellMapping&
SpatialDiscretization::GetCellMapping(const Cell& cell) const
{
  constexpr std::string_view fname = "spatial_discretization::"
                                     "GetCellMapping";
  try
  {
    if (Grid().IsCellLocal(cell.global_id))
      return *cell_mappings_.at(cell.local_id);
    else
      return *nb_cell_mappings_.at(cell.global_id);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(std::string(fname) + ": Failed to obtain cell mapping.");
  }
}

SpatialDiscretizationType
SpatialDiscretization::Type() const
{
  return type_;
}

const MeshContinuum&
SpatialDiscretization::Grid() const
{
  return ref_grid_;
}

CoordinateSystemType
SpatialDiscretization::GetCoordinateSystemType() const
{
  return coord_sys_type_;
}

size_t
SpatialDiscretization::GetNumLocalDOFs(const UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return local_base_block_size_ * N;
}

size_t
SpatialDiscretization::GetNumGlobalDOFs(const UnknownManager& unknown_manager) const
{
  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  return globl_base_block_size_ * N;
}

size_t
SpatialDiscretization::GetNumLocalAndGhostDOFs(const UnknownManager& unknown_manager) const
{
  return GetNumLocalDOFs(unknown_manager) + GetNumGhostDOFs(unknown_manager);
}

size_t
SpatialDiscretization::GetCellNumNodes(const Cell& cell) const
{
  return GetCellMapping(cell).NumNodes();
}

const std::vector<Vector3>&
SpatialDiscretization::GetCellNodeLocations(const Cell& cell) const
{
  return GetCellMapping(cell).NodeLocations();
}

std::pair<std::set<uint32_t>, std::set<uint32_t>>
SpatialDiscretization::MakeCellInternalAndBndryNodeIDs(const Cell& cell) const
{
  const auto& cell_mapping = GetCellMapping(cell);
  const size_t num_faces = cell.faces.size();
  const size_t num_nodes = cell_mapping.NumNodes();

  // Determine which nodes are on the boundary
  std::set<uint32_t> boundary_nodes;
  for (size_t f = 0; f < num_faces; ++f)
  {
    if (not cell.faces[f].has_neighbor)
    {
      const size_t num_face_nodes = cell_mapping.NumFaceNodes(f);
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
        boundary_nodes.insert(cell_mapping.MapFaceNode(f, fi));
    }
  } // for f

  // Determine non-boundary nodes
  std::set<uint32_t> internal_nodes;
  for (size_t i = 0; i < num_nodes; ++i)
    if (boundary_nodes.find(i) == boundary_nodes.end())
      internal_nodes.insert(i);

  return {internal_nodes, boundary_nodes};
}

std::vector<std::vector<std::vector<int>>>
SpatialDiscretization::MakeInternalFaceNodeMappings(const double tolerance) const
{
  const auto& grid = this->ref_grid_;

  std::vector<std::vector<std::vector<int>>> cell_adj_mapping;
  for (const auto& cell : grid.local_cells)
  {
    const auto& cell_mapping = this->GetCellMapping(cell);
    const auto& node_locations = cell_mapping.NodeLocations();
    const size_t num_faces = cell.faces.size();

    std::vector<std::vector<int>> per_face_adj_mapping;

    for (size_t f = 0; f < num_faces; ++f)
    {
      const auto& face = cell.faces[f];
      const auto num_face_nodes = cell_mapping.NumFaceNodes(f);
      std::vector<int> face_adj_mapping(num_face_nodes, -1);
      if (face.has_neighbor)
      {
        const auto& adj_cell = grid.cells[face.neighbor_id];
        const auto& adj_cell_mapping = this->GetCellMapping(adj_cell);
        const auto& adj_node_locations = adj_cell_mapping.NodeLocations();
        const size_t adj_num_nodes = adj_cell_mapping.NumNodes();

        for (size_t fi = 0; fi < num_face_nodes; ++fi)
        {
          const int i = cell_mapping.MapFaceNode(f, fi);
          const auto& ivec3 = node_locations[i];

          for (size_t ai = 0; ai < adj_num_nodes; ++ai)
          {
            const auto& aivec3 = adj_node_locations[ai];
            if ((ivec3 - aivec3).NormSquare() < tolerance)
            {
              face_adj_mapping[fi] = static_cast<int>(ai);
              break;
            }
          } // for ai
          if (face_adj_mapping[fi] < 0)
            throw std::logic_error("Face node mapping failed");
        } // for fi
      }   // if internal face

      per_face_adj_mapping.push_back(std::move(face_adj_mapping));
    } // for face

    cell_adj_mapping.push_back(std::move(per_face_adj_mapping));
  } // for cell

  return cell_adj_mapping;
}

void
SpatialDiscretization::CopyVectorWithUnknownScope(const std::vector<double>& from_vector,
                                                  std::vector<double>& to_vector,
                                                  const UnknownManager& from_vec_uk_structure,
                                                  const unsigned int from_vec_uk_id,
                                                  const UnknownManager& to_vec_uk_structure,
                                                  const unsigned int to_vec_uk_id) const
{
  const std::string fname = "spatial_discretization::"
                            "CopyVectorWithUnknownScope";
  const auto& ukmanF = from_vec_uk_structure;
  const auto& ukmanT = to_vec_uk_structure;
  const auto& ukidF = from_vec_uk_id;
  const auto& ukidT = to_vec_uk_id;
  try
  {
    const auto& ukA = from_vec_uk_structure.unknowns.at(from_vec_uk_id);
    const auto& ukB = to_vec_uk_structure.unknowns.at(to_vec_uk_id);

    if (ukA.num_components != ukB.num_components)
      throw std::logic_error(fname + " Unknowns do not have the "
                                     "same number of components");

    const size_t num_comps = ukA.num_components;

    for (const auto& cell : ref_grid_.local_cells)
    {
      const auto& cell_mapping = this->GetCellMapping(cell);
      const size_t num_nodes = cell_mapping.NumNodes();

      for (size_t i = 0; i < num_nodes; ++i)
      {
        for (size_t c = 0; c < num_comps; ++c)
        {
          const int64_t fmap = MapDOFLocal(cell, i, ukmanF, ukidF, c);
          const int64_t imap = MapDOFLocal(cell, i, ukmanT, ukidT, c);

          to_vector[imap] = from_vector[fmap];
        } // for component c
      }   // for node i
    }     // for cell
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range(fname + ": either from_vec_uk_id or to_vec_uk_id is "
                                    "out of range for its respective "
                                    "unknown manager.");
  }
}

void
SpatialDiscretization::LocalizePETScVector(Vec petsc_vector,
                                           std::vector<double>& local_vector,
                                           const UnknownManager& unknown_manager) const
{
  size_t num_local_dofs = GetNumLocalDOFs(unknown_manager);

  CopyVecToSTLvector(petsc_vector, local_vector, num_local_dofs);
}

void
SpatialDiscretization::LocalizePETScVectorWithGhosts(Vec petsc_vector,
                                                     std::vector<double>& local_vector,
                                                     const UnknownManager& unknown_manager) const
{
  size_t num_local_dofs = GetNumLocalAndGhostDOFs(unknown_manager);

  CopyVecToSTLvectorWithGhosts(petsc_vector, local_vector, num_local_dofs);
}

double
SpatialDiscretization::CartesianSpatialWeightFunction(const Vector3& point)
{
  return 1.0;
}

double
SpatialDiscretization::CylindricalRZSpatialWeightFunction(const Vector3& point)
{
  return 2.0 * M_PI * point[0];
}

double
SpatialDiscretization::Spherical1DSpatialWeightFunction(const Vector3& point)
{
  const double r = point[2];
  return 4.0 * M_PI * r * r;
}

SpatialDiscretization::SpatialWeightFunction
SpatialDiscretization::GetSpatialWeightingFunction() const
{
  switch (coord_sys_type_)
  {
    case CoordinateSystemType::CARTESIAN:
      return CartesianSpatialWeightFunction;
    case CoordinateSystemType::CYLINDRICAL:
      return CylindricalRZSpatialWeightFunction;
    case CoordinateSystemType::SPHERICAL:
      return Spherical1DSpatialWeightFunction;
    case CoordinateSystemType::UNDEFINED:
    default:
      OpenSnLogicalError("Coordinate system undefined.");
  }
}

} // namespace opensn
