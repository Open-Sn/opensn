#include "FiniteVolume.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "mesh/Cell/cell.h"

#include "chi_runtime.h"
#include "chi_log.h"

#define sc_int64 static_cast<int64_t>

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Maps a finite volume degree of freedom using an unknown manager.*/
int64_t FiniteVolume::MapDOF(
  const chi_mesh::Cell& cell,
  const unsigned int,
  const chi_math::UnknownManager& unknown_manager,
  const unsigned int unknown_id,
  const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type_;

  const size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  const size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  const size_t num_local_cells = ref_grid_.local_cells.size();

  if (component >= num_unknowns) return -1;

  int64_t address = -1;
  if (cell.partition_id_ == Chi::mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = sc_int64(local_block_address_) * num_unknowns +
                num_local_cells * block_id + cell.local_id_;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = sc_int64(local_block_address_) * num_unknowns +
                cell.local_id_ * num_unknowns + block_id;
  }
  else
  {
    const uint64_t ghost_local_id =
      neighbor_cell_local_ids_.at(cell.global_id_);

    if (storage == chi_math::UnknownStorageType::BLOCK)
      address =
        sc_int64(locJ_block_address_[cell.partition_id_]) * num_unknowns +
        locJ_block_size_[cell.partition_id_] * block_id + ghost_local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address =
        sc_int64(locJ_block_address_[cell.partition_id_]) * num_unknowns +
        ghost_local_id * num_unknowns + block_id;
  }

  return address;
}

// ###################################################################
/**Maps a finite volume degree of freedom to a local address using
 * an unknown manager.*/
int64_t
FiniteVolume::MapDOFLocal(
  const chi_mesh::Cell& cell,
  const unsigned int,
  const chi_math::UnknownManager& unknown_manager,
  const unsigned int unknown_id,
  const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type_;

  const size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  const size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  const size_t num_local_cells = ref_grid_.local_cells.size();

  if (component >= num_unknowns) return -1;

  int64_t address = -1;
  if (cell.partition_id_ == Chi::mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = sc_int64(num_local_cells) * block_id + cell.local_id_;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = sc_int64(cell.local_id_) * num_unknowns + block_id;
  }
  else
  {
    const size_t num_local_dofs = GetNumLocalDOFs(unknown_manager);
    const size_t num_ghost_nodes = GetNumGhostDOFs(UNITARY_UNKNOWN_MANAGER);
    const uint64_t ghost_local_id =
      ref_grid_.cells.GetGhostLocalID(cell.global_id_);

    if (storage == chi_math::UnknownStorageType::BLOCK)
      address = sc_int64(num_local_dofs) +
                sc_int64(num_ghost_nodes) * block_id + ghost_local_id;
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = sc_int64(num_local_dofs) +
                num_unknowns * sc_int64(ghost_local_id) + block_id;
  }

  return address;
}

} // namespace chi_math::spatial_discretization
