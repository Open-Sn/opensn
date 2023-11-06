#pragma once

#include "framework/mesh/Cell/cell.h"

#include <map>

namespace chi_mesh
{
// ##################################################
/**Handles all global index queries.*/
class GlobalCellHandler
{
  friend class MeshContinuum;

private:
  std::vector<std::unique_ptr<chi_mesh::Cell>>& local_cells_ref_;
  std::vector<std::unique_ptr<chi_mesh::Cell>>& ghost_cells_ref_;

  std::map<uint64_t, uint64_t>& global_cell_id_to_native_id_map;
  std::map<uint64_t, uint64_t>& global_cell_id_to_foreign_id_map;

private:
  explicit GlobalCellHandler(std::vector<std::unique_ptr<chi_mesh::Cell>>& in_native_cells,
                             std::vector<std::unique_ptr<chi_mesh::Cell>>& in_foreign_cells,
                             std::map<uint64_t, uint64_t>& in_global_cell_id_to_native_id_map,
                             std::map<uint64_t, uint64_t>& in_global_cell_id_to_foreign_id_map)
    : local_cells_ref_(in_native_cells),
      ghost_cells_ref_(in_foreign_cells),
      global_cell_id_to_native_id_map(in_global_cell_id_to_native_id_map),
      global_cell_id_to_foreign_id_map(in_global_cell_id_to_foreign_id_map)
  {
  }

public:
  void push_back(std::unique_ptr<chi_mesh::Cell> new_cell);
  chi_mesh::Cell& operator[](uint64_t cell_global_index);
  const chi_mesh::Cell& operator[](uint64_t cell_global_index) const;

  size_t GetNumGhosts() const { return global_cell_id_to_foreign_id_map.size(); }

  std::vector<uint64_t> GetGhostGlobalIDs() const;

  uint64_t GetGhostLocalID(uint64_t cell_global_index) const;
};

} // namespace chi_mesh
