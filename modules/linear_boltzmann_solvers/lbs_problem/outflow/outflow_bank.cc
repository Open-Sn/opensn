// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/outflow_bank.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

namespace
{

inline std::uint64_t
Merge(const std::uint32_t& high, const std::uint32_t& low)
{
  return (static_cast<std::uint64_t>(high) << 32) | low;
}

} // namespace

OutflowBank::OutflowBank(const MeshContinuum& grid, unsigned int num_groups)
  : views_(grid.local_cells.size())
{
  for (const auto& cell : grid.local_cells)
    views_[cell.local_id] = CellOutflowView(cell.faces.size(), num_groups);

  std::vector<std::int64_t> cell_offsets(grid.local_cells.size());
  cell_offsets.assign(cell_offsets.size(), static_cast<std::int64_t>(-1));

  std::size_t num_boundary_faces = 0;
  for (const auto& cell : grid.local_cells)
  {
    std::size_t face_offset_count = 0;
    for (std::size_t f = 0; f < cell.faces.size(); ++f)
    {
      const auto& face = cell.faces[f];
      if (face.has_neighbor)
        continue;

      auto& cell_view = views_[cell.local_id];
      if (not cell_view.HasFaceOffsets())
      {
        cell_view.InitializeFaceOffsets();
        cell_offsets[cell.local_id] = static_cast<std::int64_t>(num_boundary_faces);
      }
      cell_view.SetFaceOffset(f, static_cast<std::int64_t>(face_offset_count * num_groups));
      ++face_offset_count;
      cellface_map_[Merge(cell.local_id, f)] = num_boundary_faces * num_groups;
      num_boundary_faces++;
    }
  }

  const std::size_t size = num_boundary_faces * num_groups;
  outflows_.reserve(size);
  outflows_.assign(size, 0.0);
  for (std::size_t cell_local_idx = 0; cell_local_idx < grid.local_cells.size(); ++cell_local_idx)
  {
    const auto& cell_offset = cell_offsets[cell_local_idx];
    if (cell_offset < 0)
      continue;
    views_[cell_local_idx].Assign(outflows_.data() + cell_offset * num_groups);
  }
}

std::uint64_t
OutflowBank::GetOffset(std::uint32_t cell_local_idx, std::uint32_t face_idx) const
{
  return cellface_map_.at(Merge(cell_local_idx, face_idx));
}

} // namespace opensn
