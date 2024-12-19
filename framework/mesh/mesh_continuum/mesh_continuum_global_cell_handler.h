// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/cell/cell.h"
#include <map>

namespace opensn
{

/// Handles all global index queries.
class GlobalCellHandler
{
  friend class MeshContinuum;

private:
  std::vector<std::shared_ptr<Cell>>& local_cells_ref_;
  std::vector<std::shared_ptr<Cell>>& ghost_cells_ref_;

  /// Global to local ID map
  std::map<uint64_t, uint64_t>& global_to_local_map_;
  /// Global to ghost ID map
  std::map<uint64_t, uint64_t>& global_to_ghost_map_;

private:
  explicit GlobalCellHandler(std::vector<std::shared_ptr<Cell>>& native_cells,
                             std::vector<std::shared_ptr<Cell>>& foreign_cells,
                             std::map<uint64_t, uint64_t>& global_to_local_map,
                             std::map<uint64_t, uint64_t>& global_to_ghost_map)
    : local_cells_ref_(native_cells),
      ghost_cells_ref_(foreign_cells),
      global_to_local_map_(global_to_local_map),
      global_to_ghost_map_(global_to_ghost_map)
  {
  }

public:
  /**
   * Adds a new cell to the appropriate category (local or ghost).
   * @param new_cell The cell to add.
   */
  void PushBack(std::shared_ptr<Cell> new_cell);

  /// Returns a reference to a cell given its global cell index.
  Cell& operator[](uint64_t cell_global_index);

  /// Returns a const reference to a cell given its global cell index.
  const Cell& operator[](uint64_t cell_global_index) const;

  /// Returns the the total number of ghost cells
  size_t GhostCellCount() const { return global_to_ghost_map_.size(); }

  /**
   * Returns the cell global ids of all ghost cells. These are cells that neighbors to this
   * partition's cells but are on a different partition.
   */
  std::vector<uint64_t> GetGhostGlobalIDs() const;

  /**
   * Returns the local storage address of a ghost cell. If the ghost is not truly a ghost then -1 is
   * returned, but is wasteful and therefore the user of this function should implement code to
   * prevent it.
   */
  uint64_t GetGhostLocalID(uint64_t cell_global_index) const;
};

} // namespace opensn
