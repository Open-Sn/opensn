// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_mapping/mesh_mapping.h"
#include "framework/logging/log.h"
#include "framework/mesh/cell/cell.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/math/spatial_discretization/finite_element/piecewise_linear/piecewise_linear_continuous.h"

#include <sstream>

namespace opensn
{

const std::size_t MeshMapping::invalid_face_index = std::numeric_limits<std::size_t>::max();

MeshMapping::CoarseMapping::CoarseMapping(const Cell& coarse_cell)
  : fine_faces{coarse_cell.faces.size()}
{
}

MeshMapping::FineMapping::FineMapping(const Cell& fine_cell)
  : coarse_cell(nullptr), coarse_faces(fine_cell.faces.size(), MeshMapping::invalid_face_index)
{
}

void
MeshMapping::Build(const std::shared_ptr<MeshContinuum> fine_grid,
                   const std::shared_ptr<MeshContinuum> coarse_grid)

{
  if (opensn::mpi_comm.size() > 1)
    OpenSnLogicalError("MeshMapping is not currently supported in parallel.");
  if (fine_grid->GetDimension() != coarse_grid->GetDimension())
    OpenSnLogicalError("Grid dimensions are not equal for mapping. Fine dimension = " +
                       std::to_string(fine_grid->GetDimension()) +
                       ", coarse dimension = " + std::to_string(coarse_grid->GetDimension()) + ".");

  coarse_to_fine_.clear();
  fine_to_coarse_.clear();

  // Instantiate the maps; constructors take the cell to size the face maps.
  for (const auto& coarse_cell : coarse_grid->local_cells)
    coarse_to_fine_.emplace(&coarse_cell, coarse_cell);
  for (const auto& fine_cell : fine_grid->local_cells)
    fine_to_coarse_.emplace(&fine_cell, fine_cell);

  // Volumetric mapping; find the coarse cell that contains a fine cell centroid
  for (auto& [fine_cell_ptr, fine_mapping] : fine_to_coarse_)
  {
    const auto& fine_cell = *fine_cell_ptr;
    for (const auto& coarse_cell : coarse_grid->local_cells)
      if (coarse_grid->CheckPointInsideCell(coarse_cell, fine_cell.centroid))
      {
        fine_mapping.coarse_cell = &coarse_cell;
        break;
      }

    if (!fine_mapping.coarse_cell)
    {
      std::ostringstream oss;
      oss << "Failed to find a corresponding coarse cell for fine cell " << fine_cell.global_id
          << " with centroid " << fine_cell.centroid.PrintStr() << ".";
      OpenSnLogicalError(oss.str());
    }

    coarse_to_fine_.at(fine_mapping.coarse_cell).fine_cells.push_back(fine_cell_ptr);
  }

  // Ensure that coarse cell volume is equal to the sum of the fine cell volumes contained within it
  auto fine_sdm_ptr = PieceWiseLinearContinuous::New(fine_grid);
  auto& fine_sdm = *fine_sdm_ptr;
  auto coarse_sdm_ptr = PieceWiseLinearContinuous::New(coarse_grid);
  auto& coarse_sdm = *coarse_sdm_ptr;
  for (const auto& [coarse_cell_ptr, coarse_mapping] : coarse_to_fine_)
  {
    const auto& coarse_cell = *coarse_cell_ptr;
    const auto& coarse_cell_mapping = coarse_sdm.GetCellMapping(coarse_cell);
    const auto coarse_cell_volume = coarse_cell_mapping.GetCellVolume();
    double total_fine_volume = 0;
    for (const auto fine_cell_ptr : coarse_mapping.fine_cells)
    {
      const auto& fine_cell = *fine_cell_ptr;
      const auto& fine_cell_mapping = fine_sdm.GetCellMapping(fine_cell);
      total_fine_volume += fine_cell_mapping.GetCellVolume();
    }
    if (std::abs(total_fine_volume - coarse_cell_volume) > 1.e-6)
    {
      std::ostringstream oss;
      oss << "Coarse cell " << coarse_cell.global_id << " with centroid "
          << coarse_cell.centroid.PrintStr() << " volumetric mapping failed.";
      OpenSnLogicalError(oss.str());
    }
  }

  // Surface mapping; find the coarse cell face that contains a fine cell face centroid
  for (auto& [fine_cell_ptr, fine_mapping] : fine_to_coarse_)
  {
    const auto& fine_cell = *fine_cell_ptr;
    const auto& coarse_cell = *fine_mapping.coarse_cell;
    auto& coarse_mapping = coarse_to_fine_.at(&coarse_cell);
    for (size_t fine_face_i = 0; fine_face_i < fine_cell.faces.size(); ++fine_face_i)
    {
      const auto& fine_face = fine_cell.faces[fine_face_i];
      for (size_t coarse_face_i = 0; coarse_face_i < coarse_cell.faces.size(); ++coarse_face_i)
      {
        if (coarse_grid->CheckPointInsideCellFace(coarse_cell, coarse_face_i, fine_face.centroid))
        {
          coarse_mapping.fine_faces[coarse_face_i].emplace_back(fine_cell_ptr, fine_face_i);
          fine_mapping.coarse_faces[fine_face_i] = coarse_face_i;
          break;
        }
      }
    }
  }

  // Ensure that coarse cell area is equal to the sum of the fine cell areas contained within it
  for (const auto& [coarse_cell_ptr, coarse_mapping] : coarse_to_fine_)
  {
    const auto& coarse_cell = *coarse_cell_ptr;
    for (size_t coarse_face_i = 0; coarse_face_i < coarse_cell.faces.size(); ++coarse_face_i)
    {
      double total_fine_face_area = 0;
      const auto& fine_faces = coarse_mapping.fine_faces[coarse_face_i];
      for (const auto& [fine_cell_ptr, fine_face_i] : fine_faces)
      {
        const auto& fine_face = fine_cell_ptr->faces[fine_face_i];
        total_fine_face_area += fine_face.ComputeFaceArea(fine_grid.get());
      }
      const auto& coarse_face = coarse_cell.faces[coarse_face_i];
      const auto coarse_face_area = coarse_face.ComputeFaceArea(coarse_grid.get());
      if (std::abs(total_fine_face_area - coarse_face_area) > 1.e-6)
      {
        std::ostringstream oss;
        oss << "Coarse cell " << coarse_cell.global_id << " face " << coarse_face_i
            << " with centroid " << coarse_face.centroid.PrintStr() << " surface mapping failed.";
        OpenSnLogicalError(oss.str());
      }
    }
  }
}

const MeshMapping::CoarseMapping&
MeshMapping::GetCoarseMapping(const Cell& coarse_cell) const
{
  const auto it = coarse_to_fine_.find(&coarse_cell);
  if (it == coarse_to_fine_.end())
    OpenSnLogicalError("MeshMapping::GetCoarseMapping(): Coarse cell not found in mapping.");
  return it->second;
}

const MeshMapping::FineMapping&
MeshMapping::GetFineMapping(const Cell& fine_cell) const
{
  const auto it = fine_to_coarse_.find(&fine_cell);
  if (it == fine_to_coarse_.end())
    OpenSnLogicalError("MeshMapping::GetFineMapping(): Fine cell not found in mapping.");
  return it->second;
}

} // namespace opensn
