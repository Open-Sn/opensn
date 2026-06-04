// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

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

} // namespace opensn
