#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/chi_runtime.h"
#include "framework/mpi/chi_mpi.h"
#include "framework/logging/chi_log.h"

/**Adds a new cell to grid registry.*/
void
chi_mesh::GlobalCellHandler::push_back(std::unique_ptr<chi_mesh::Cell> new_cell)
{
  if (new_cell->partition_id_ == static_cast<uint64_t>(Chi::mpi.location_id))
  {
    new_cell->local_id_ = local_cells_ref_.size();

    local_cells_ref_.push_back(std::move(new_cell));

    const auto& cell = local_cells_ref_.back();

    global_cell_id_to_native_id_map.insert(
      std::make_pair(cell->global_id_, local_cells_ref_.size() - 1));
  }
  else
  {
    ghost_cells_ref_.push_back(std::move(new_cell));

    const auto& cell = ghost_cells_ref_.back();

    global_cell_id_to_foreign_id_map.insert(
      std::make_pair(cell->global_id_, ghost_cells_ref_.size() - 1));
  }
}

/**Returns a reference to a cell given its global cell index.*/
chi_mesh::Cell&
chi_mesh::GlobalCellHandler::operator[](uint64_t cell_global_index)
{
  auto native_location = global_cell_id_to_native_id_map.find(cell_global_index);

  if (native_location != global_cell_id_to_native_id_map.end())
    return *local_cells_ref_[native_location->second];
  else
  {
    auto foreign_location = global_cell_id_to_foreign_id_map.find(cell_global_index);
    if (foreign_location != global_cell_id_to_foreign_id_map.end())
      return *ghost_cells_ref_[foreign_location->second];
  }

  std::stringstream ostr;
  ostr << "chi_mesh::MeshContinuum::cells. Mapping error."
       << "\n"
       << cell_global_index;

  throw std::invalid_argument(ostr.str());
}

/**Returns a const reference to a cell given its global cell index.*/
const chi_mesh::Cell&
chi_mesh::GlobalCellHandler::operator[](uint64_t cell_global_index) const
{
  auto native_location = global_cell_id_to_native_id_map.find(cell_global_index);

  if (native_location != global_cell_id_to_native_id_map.end())
    return *local_cells_ref_[native_location->second];
  else
  {
    auto foreign_location = global_cell_id_to_foreign_id_map.find(cell_global_index);
    if (foreign_location != global_cell_id_to_foreign_id_map.end())
      return *ghost_cells_ref_[foreign_location->second];
  }

  std::stringstream ostr;
  ostr << "chi_mesh::MeshContinuum::cells. Mapping error."
       << "\n"
       << cell_global_index;

  throw std::invalid_argument(ostr.str());
}

/**Returns the cell global ids of all ghost cells. These are cells that
 * neighbors to this partition's cells but are on a different
 * partition.*/
std::vector<uint64_t>
chi_mesh::GlobalCellHandler::GetGhostGlobalIDs() const
{
  std::vector<uint64_t> ids;
  ids.reserve(GetNumGhosts());

  for (auto& cell : ghost_cells_ref_)
    ids.push_back(cell->global_id_);

  return ids;
}

/**Returns the local storage address of a ghost cell. If the
 * ghost is not truly a ghost then -1 is returned, but is wasteful and
 * therefore the user of this function should implement code
 * to prevent it.*/
uint64_t
chi_mesh::GlobalCellHandler::GetGhostLocalID(uint64_t cell_global_index) const
{
  auto foreign_location = global_cell_id_to_foreign_id_map.find(cell_global_index);

  if (foreign_location != global_cell_id_to_foreign_id_map.end()) return foreign_location->second;

  std::stringstream ostr;
  ostr << "Grid GetGhostLocalID failed to find cell " << cell_global_index;

  throw std::invalid_argument(ostr.str());
}
