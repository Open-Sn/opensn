#include "LagrangeContinuous.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#define sc_int64 static_cast<int64_t>

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Maps a vertex id according to a developed node ordering.*/
int64_t LagrangeContinuous::MapDOF(
  const chi_mesh::Cell& cell,
  const unsigned int node,
  const chi_math::UnknownManager& unknown_manager,
  const unsigned int unknown_id,
  const unsigned int component /*=0*/) const
{
  const uint64_t vertex_id = cell.vertex_ids_[node];

  ChiLogicalErrorIf(node_mapping_.count(vertex_id) == 0,
                    std::string("Bad trouble mapping vertex ") +
                      std::to_string(vertex_id));
  const int64_t global_id = node_mapping_.at(vertex_id);

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  auto storage = unknown_manager.dof_storage_type_;

  int64_t address = -1;
  if (storage == chi_math::UnknownStorageType::BLOCK)
  {
    for (int locJ = 0; locJ < Chi::mpi.process_count; ++locJ)
    {
      const int64_t local_id = global_id - sc_int64(locJ_block_address_[locJ]);

      if (local_id < 0 or local_id >= locJ_block_size_[locJ]) continue;

      address = sc_int64(locJ_block_address_[locJ] * num_unknowns) +
                sc_int64(locJ_block_size_[locJ] * block_id) + local_id;
      break;
    }
  }
  else if (storage == chi_math::UnknownStorageType::NODAL)
    address = global_id * sc_int64(num_unknowns) + sc_int64(block_id);

  return address;
}

// ###################################################################
/**Maps a vertex id according to a developed node ordering.*/
int64_t LagrangeContinuous::MapDOFLocal(
  const chi_mesh::Cell& cell,
  const unsigned int node,
  const chi_math::UnknownManager& unknown_manager,
  const unsigned int unknown_id,
  const unsigned int component /*=0*/) const
{
  const uint64_t vertex_id = cell.vertex_ids_[node];

  ChiLogicalErrorIf(node_mapping_.count(vertex_id) == 0, "Bad trouble");
  const int64_t node_global_id = node_mapping_.at(vertex_id);

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id = unknown_manager.MapUnknown(unknown_id, component);
  auto storage = unknown_manager.dof_storage_type_;

  const int64_t local_id = node_global_id - sc_int64(local_block_address_);
  const bool is_local = not(local_id < 0 or local_id >= local_base_block_size_);

  int64_t address = -1;
  if (is_local)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
    {
      address = sc_int64(local_base_block_size_ * block_id) + local_id;
    }
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address = local_id * sc_int64(num_unknowns) + sc_int64(block_id);
  } // if is_local
  else
  {
    const size_t num_local_dofs = GetNumLocalDOFs(unknown_manager);
    int64_t ghost_local_node_id = -1;
    int64_t counter = 0;
    for (const auto& vid_gnid : ghost_node_mapping_)
    {
      if (node_global_id == vid_gnid.second)
      {
        ghost_local_node_id = counter;
        break;
      }
      ++counter;
    }
    if (storage == chi_math::UnknownStorageType::BLOCK)
    {
      address =
        sc_int64(ghost_node_mapping_.size() * block_id) + ghost_local_node_id;
    }
    else if (storage == chi_math::UnknownStorageType::NODAL)
      address =
        ghost_local_node_id * sc_int64(num_unknowns) + sc_int64(block_id);

    address += sc_int64(num_local_dofs);
  }

  return address;
}

} // namespace chi_math::spatial_discretization
