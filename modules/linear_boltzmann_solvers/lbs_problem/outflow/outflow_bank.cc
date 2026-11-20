// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/lbs_problem/outflow/outflow_bank.h"
#include "framework/mesh/mesh/mesh.h"

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

OutflowBank::OutflowBank(const Mesh& grid, unsigned int num_groups, bool include_internal_faces)
  : views_(grid.GetLocalCellCount())
{
  for (std::uint32_t cell_local_id = 0; cell_local_id < grid.GetLocalCellCount(); ++cell_local_id)
  {
    const auto& cell = grid.GetLocalCell(cell_local_id);
    views_[cell_local_id] = CellOutflowView(grid.GetCellFaceCount(cell_local_id), num_groups);
  }

  std::vector<std::int64_t> cell_offsets(grid.GetLocalCellCount(), -1);

  std::size_t num_faces = 0;
  for (std::uint32_t cell_local_id = 0; cell_local_id < grid.GetLocalCellCount(); ++cell_local_id)
  {
    const auto& cell = grid.GetLocalCell(cell_local_id);
    auto& cell_view = views_[cell_local_id];
    cell_view.InitializeFaceOffsets();
    std::size_t cell_num_stored_faces = 0;

    auto cell_faces = grid.GetCellFaces(cell_local_id);
    for (std::size_t f = 0; f < cell_faces.size(); ++f)
    {
      const auto& face = cell_faces[f];
      if (face.has_neighbor and not include_internal_faces)
        continue;

      if (cell_num_stored_faces == 0)
        cell_offsets[cell_local_id] = static_cast<std::int64_t>(num_faces);

      cell_view.SetFaceOffset(f, static_cast<std::int64_t>(cell_num_stored_faces) * num_groups);
      cellface_map_[Merge(cell_local_id, f)] = num_faces * num_groups;
      ++num_faces;
      ++cell_num_stored_faces;
    }
  }

  const std::size_t size = num_faces * num_groups;
  outflows_.reserve(size);
  outflows_.assign(size, 0.0);
  for (std::size_t cell_local_idx = 0; cell_local_idx < grid.GetLocalCellCount(); ++cell_local_idx)
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

bool
OutflowBank::HasOffset(std::uint32_t cell_local_idx, std::uint32_t face_idx) const
{
  return cellface_map_.count(Merge(cell_local_idx, face_idx)) != 0;
}

} // namespace opensn
