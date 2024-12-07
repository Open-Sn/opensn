// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/mesh/mesh_continuum/mesh_continuum.h"

namespace opensn
{

Cell&
LocalCellHandler::operator[](uint64_t cell_local_index)
{
  if (local_cells.empty())
  {
    throw std::out_of_range("LocalCells attempted to access local cell " +
                            std::to_string(cell_local_index) +
                            ". Local cells are empty. Check the partitioning.");
  }

  if (cell_local_index >= local_cells.size())
  {
    throw std::out_of_range("Local cell index out of range: " + std::to_string(cell_local_index) +
                            " (max: " + std::to_string(local_cells.size() - 1) + ").");
  }

  return *local_cells[cell_local_index];
}

const Cell&
LocalCellHandler::operator[](uint64_t cell_local_index) const
{
  if (local_cells.empty())
  {
    throw std::out_of_range("LocalCells attempted to access local cell " +
                            std::to_string(cell_local_index) +
                            ". Local cells are empty. Check the partitioning.");
  }

  if (cell_local_index >= local_cells.size())
  {
    throw std::out_of_range("Local cell index out of range: " + std::to_string(cell_local_index) +
                            " (max: " + std::to_string(local_cells.size() - 1) + ").");
  }

  return *local_cells[cell_local_index];
}

} // namespace opensn
