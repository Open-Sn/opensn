#include "PieceWiseLinearDiscontinuous.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#define sc_int64 static_cast<int64_t>

namespace chi_math::spatial_discretization
{
// ###################################################################
/**Provides a mapping of cell's DOF from a DFEM perspective.*/
int64_t PieceWiseLinearDiscontinuous::MapDOF(
  const chi_mesh::Cell& cell,
  const unsigned int node,
  const chi_math::UnknownManager& unknown_manager,
  const unsigned int unknown_id,
  const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type_;

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id = unknown_manager.MapUnknown(unknown_id, component);

  if (cell.partition_id_ == Chi::mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
    {
      int64_t address = sc_int64(local_block_address_ * num_unknowns) +
                        cell_local_block_address_[cell.local_id_] +
                        local_base_block_size_ * block_id + node;
      return address;
    }
    else if (storage == chi_math::UnknownStorageType::NODAL)
    {
      int64_t address =
        sc_int64(local_block_address_ * num_unknowns) +
        cell_local_block_address_[cell.local_id_] * num_unknowns +
        node * num_unknowns + block_id;
      return address;
    }
  }
  else
  {
    int index = 0;
    bool found = false;
    for (auto neighbor_info : neighbor_cell_block_address_)
    {
      if (neighbor_info.first == cell.global_id_)
      {
        found = true;
        break;
      }
      ++index;
    }

    if (!found)
    {
      Chi::log.LogAllError()
        << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
        << "with global index " << cell.global_id_ << " and partition-ID "
        << cell.partition_id_;
      Chi::Exit(EXIT_FAILURE);
    }

    if (storage == chi_math::UnknownStorageType::BLOCK)
    {
      int64_t address = sc_int64(neighbor_cell_block_address_[index].second) +
                        locJ_block_size_[cell.partition_id_] * block_id + node;
      return address;
    }
    else if (storage == chi_math::UnknownStorageType::NODAL)
    {
      int64_t address =
        sc_int64(neighbor_cell_block_address_[index].second * num_unknowns) +
        node * num_unknowns + block_id;
      return address;
    }
  }

  return -1;
}

// ###################################################################
/**Provides a mapping of cell's DOF from a DFEM perspective.*/
int64_t PieceWiseLinearDiscontinuous::MapDOFLocal(
  const chi_mesh::Cell& cell,
  const unsigned int node,
  const chi_math::UnknownManager& unknown_manager,
  const unsigned int unknown_id,
  const unsigned int component) const
{
  auto storage = unknown_manager.dof_storage_type_;

  size_t num_unknowns = unknown_manager.GetTotalUnknownStructureSize();
  size_t block_id = unknown_manager.MapUnknown(unknown_id, component);

  if (cell.partition_id_ == Chi::mpi.location_id)
  {
    if (storage == chi_math::UnknownStorageType::BLOCK)
    {
      int64_t address = sc_int64(cell_local_block_address_[cell.local_id_]) +
                        local_base_block_size_ * block_id + node;
      return address;
    }
    else if (storage == chi_math::UnknownStorageType::NODAL)
    {
      int64_t address =
        sc_int64(cell_local_block_address_[cell.local_id_] * num_unknowns) +
        node * num_unknowns + block_id;
      return address;
    }
  }
  else
  {
    int index = 0;
    bool found = false;
    for (auto neighbor_info : neighbor_cell_block_address_)
    {
      if (neighbor_info.first == cell.global_id_)
      {
        found = true;
        break;
      }
      ++index;
    }

    if (!found)
    {
      Chi::log.LogAllError()
        << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
        << "with global index " << cell.global_id_ << " and partition-ID "
        << cell.partition_id_;
      Chi::Exit(EXIT_FAILURE);
    }

    if (storage == chi_math::UnknownStorageType::BLOCK)
    {
      int64_t address = sc_int64(neighbor_cell_block_address_[index].second) +
                        locJ_block_size_[cell.partition_id_] * block_id + node;
      return address;
    }
    else if (storage == chi_math::UnknownStorageType::NODAL)
    {
      int64_t address =
        sc_int64(neighbor_cell_block_address_[index].second * num_unknowns) +
        node * num_unknowns + block_id;
      return address;
    }
  }

  return -1;
}

} // namespace chi_math::spatial_discretization
