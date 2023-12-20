#pragma once

#include "framework/mesh/cell/cell.h"

#include <map>

namespace opensn
{

/**Handles all global index queries.*/
class GlobalCellHandler
{
  friend class MeshContinuum;

private:
  std::vector<std::unique_ptr<Cell>>& local_cells_ref_;
  std::vector<std::unique_ptr<Cell>>& ghost_cells_ref_;

  std::map<uint64_t, uint64_t>& global_cell_id_to_native_id_map;
  std::map<uint64_t, uint64_t>& global_cell_id_to_foreign_id_map;

private:
  explicit GlobalCellHandler(std::vector<std::unique_ptr<Cell>>& in_native_cells,
                             std::vector<std::unique_ptr<Cell>>& in_foreign_cells,
                             std::map<uint64_t, uint64_t>& in_global_cell_id_to_native_id_map,
                             std::map<uint64_t, uint64_t>& in_global_cell_id_to_foreign_id_map)
    : local_cells_ref_(in_native_cells),
      ghost_cells_ref_(in_foreign_cells),
      global_cell_id_to_native_id_map(in_global_cell_id_to_native_id_map),
      global_cell_id_to_foreign_id_map(in_global_cell_id_to_foreign_id_map)
  {
  }

public:
  /**Adds a new cell to grid registry.*/
  void push_back(std::unique_ptr<Cell> new_cell);
  /**Returns a reference to a cell given its global cell index.*/
  Cell& operator[](uint64_t cell_global_index);
  /**Returns a const reference to a cell given its global cell index.*/
  const Cell& operator[](uint64_t cell_global_index) const;

  size_t GetNumGhosts() const { return global_cell_id_to_foreign_id_map.size(); }

  /**Returns the cell global ids of all ghost cells. These are cells that
   * neighbors to this partition's cells but are on a different
   * partition.*/
  std::vector<uint64_t> GetGhostGlobalIDs() const;

  /**Returns the local storage address of a ghost cell. If the
   * ghost is not truly a ghost then -1 is returned, but is wasteful and
   * therefore the user of this function should implement code
   * to prevent it.*/
  uint64_t GetGhostLocalID(uint64_t cell_global_index) const;
};

} // namespace opensn
