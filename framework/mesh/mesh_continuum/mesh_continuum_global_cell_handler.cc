// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

void
GlobalCellHandler::PushBack(std::unique_ptr<Cell> new_cell)
{
  if (new_cell->partition_id == static_cast<uint64_t>(opensn::mpi_comm.rank()))
  {
    new_cell->local_id = local_cells_ref_.size();
    local_cells_ref_.push_back(std::move(new_cell));
    global_to_local_map_[local_cells_ref_.back()->global_id] = local_cells_ref_.size() - 1;
  }
  else
  {
    ghost_cells_ref_.push_back(std::move(new_cell));
    global_to_ghost_map_[ghost_cells_ref_.back()->global_id] = ghost_cells_ref_.size() - 1;
  }
}

Cell&
GlobalCellHandler::operator[](uint64_t cell_global_index)
{
  auto local_it = global_to_local_map_.find(cell_global_index);
  if (local_it != global_to_local_map_.end())
    return *local_cells_ref_[local_it->second];

  auto ghost_it = global_to_ghost_map_.find(cell_global_index);
  if (ghost_it != global_to_ghost_map_.end())
    return *ghost_cells_ref_[ghost_it->second];

  throw std::out_of_range("Cell with global ID " + std::to_string(cell_global_index) +
                          " not found.");
}

const Cell&
GlobalCellHandler::operator[](uint64_t cell_global_index) const
{
  auto local_it = global_to_local_map_.find(cell_global_index);
  if (local_it != global_to_local_map_.end())
    return *local_cells_ref_[local_it->second];

  auto ghost_it = global_to_ghost_map_.find(cell_global_index);
  if (ghost_it != global_to_ghost_map_.end())
    return *ghost_cells_ref_[ghost_it->second];

  throw std::out_of_range("Cell with global ID " + std::to_string(cell_global_index) +
                          " not found.");
}

std::vector<uint64_t>
GlobalCellHandler::GetGhostGlobalIDs() const
{
  std::vector<uint64_t> ids;
  ids.reserve(GhostCellCount());

  for (auto& cell : ghost_cells_ref_)
    ids.push_back(cell->global_id);

  return ids;
}

uint64_t
GlobalCellHandler::GetGhostLocalID(uint64_t cell_global_index) const
{
  auto foreign_location = global_to_ghost_map_.find(cell_global_index);

  if (foreign_location != global_to_ghost_map_.end())
    return foreign_location->second;

  throw std::out_of_range("Cell with global ID " + std::to_string(cell_global_index) +
                          " not found.");
}

} // namespace opensn
