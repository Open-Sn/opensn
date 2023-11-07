#include "framework/math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"
#include "framework/math/UnknownManager/unknown_manager.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"
#include "framework/utils/chi_timer.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/mpi/chi_mpi.h"
#include "framework/mpi/chi_mpi_utils.h"

#define sc_int64 static_cast<int64_t>

namespace chi_math::spatial_discretization
{

PieceWiseLinearDiscontinuous::PieceWiseLinearDiscontinuous(const chi_mesh::MeshContinuum& grid,
                                                           chi_math::QuadratureOrder q_order,
                                                           chi_math::CoordinateSystemType cs_type)
  : PieceWiseLinearBase(grid, q_order, SDMType::PIECEWISE_LINEAR_DISCONTINUOUS, cs_type)
{
  CreateCellMappings();

  OrderNodes();
}

std::shared_ptr<PieceWiseLinearDiscontinuous>
PieceWiseLinearDiscontinuous::New(const chi_mesh::MeshContinuum& grid,
                                  QuadratureOrder q_order /*=QuadratureOrder::SECOND*/,
                                  CoordinateSystemType cs_type /*=CoordinateSystemType::CARTESIAN*/)

{
  const auto PWLD = SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;
  // First try to find an existing spatial discretization that matches the
  // one requested.
  for (auto& sdm : Chi::sdm_stack)
    if (sdm->Type() == PWLD and std::addressof(sdm->Grid()) == std::addressof(grid) and
        sdm->GetCoordinateSystemType() == cs_type)
    {
      auto fe_ptr = std::dynamic_pointer_cast<FiniteElementBase>(sdm);

      ChiLogicalErrorIf(not fe_ptr, "Casting failure to FE");

      if (fe_ptr->GetQuadratureOrder() != q_order) break;

      auto sdm_ptr = std::dynamic_pointer_cast<PieceWiseLinearDiscontinuous>(fe_ptr);

      ChiLogicalErrorIf(not sdm_ptr, "Casting failure");

      return sdm_ptr;
    }

  auto new_sdm = std::shared_ptr<PieceWiseLinearDiscontinuous>(
    new PieceWiseLinearDiscontinuous(grid, q_order, cs_type));

  Chi::sdm_stack.push_back(new_sdm);

  return new_sdm;
}

void
PieceWiseLinearDiscontinuous::OrderNodes()
{
  const std::string fname = __FUNCTION__;
  chi::Timer t_stage[6];

  t_stage[0].Reset();
  // Check cell views avail
  size_t num_loc_cells = ref_grid_.local_cells.size();

  // Get local DOF count and set cell_local_block_address
  cell_local_block_address_.resize(num_loc_cells, 0);

  uint64_t local_node_count = 0;
  for (const auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);
    cell_local_block_address_[cell.local_id_] = static_cast<int64_t>(local_node_count);
    local_node_count += cell_mapping.NumNodes();
  }

  // Allgather node_counts
  locJ_block_size_.assign(Chi::mpi.process_count, 0);
  MPI_Allgather(&local_node_count, // sendbuf
                1,
                MPI_UNSIGNED_LONG_LONG,  // sendcount, sendtype
                locJ_block_size_.data(), // recvbuf
                1,
                MPI_UNSIGNED_LONG_LONG, // recvcount, recvtype
                Chi::mpi.comm);         // comm

  // Assign
  // local_block_address
  uint64_t running_block_address = 0;
  for (int locI = 0; locI < Chi::mpi.process_count; ++locI)
  {
    if (locI == Chi::mpi.location_id)
      local_block_address_ = static_cast<int64_t>(running_block_address);

    running_block_address += locJ_block_size_[locI];
  }
  const uint64_t global_node_count = running_block_address;

  local_base_block_size_ = local_node_count;
  globl_base_block_size_ = global_node_count;

  // Collect ghost cell ids needing block addresses
  std::map<int, std::vector<uint64_t>> ghost_cell_ids_consolidated;

  for (uint64_t global_id : ref_grid_.cells.GetGhostGlobalIDs())
  {
    const auto& cell = ref_grid_.cells[global_id];
    const int locI = static_cast<int>(cell.partition_id_);

    std::vector<uint64_t>& locI_cell_id_list = ghost_cell_ids_consolidated[locI];

    locI_cell_id_list.push_back(cell.global_id_);
  }

  // AllToAll to get query cell-ids
  const std::map<int, std::vector<uint64_t>> query_ghost_cell_ids_consolidated =
    chi_mpi_utils::MapAllToAll(ghost_cell_ids_consolidated, MPI_UNSIGNED_LONG_LONG);

  // Map all query cell-ids
  std::map<int, std::vector<uint64_t>> mapped_ghost_cell_ids_consolidated;
  for (const auto& [pid, cell_id_list] : query_ghost_cell_ids_consolidated)
  {
    std::vector<uint64_t>& map_list = mapped_ghost_cell_ids_consolidated[pid];

    for (uint64_t cell_global_id : cell_id_list)
    {
      const auto& cell = ref_grid_.cells[cell_global_id];

      const uint64_t cell_block_address =
        local_block_address_ + cell_local_block_address_[cell.local_id_];
      map_list.push_back(cell_block_address);
    }
  }

  // Communicate back the mapping
  const std::map<int, std::vector<uint64_t>> global_id_mapping =
    chi_mpi_utils::MapAllToAll(mapped_ghost_cell_ids_consolidated, MPI_UNSIGNED_LONG_LONG);

  // Process global id mapping
  for (const auto& [pid, mapping_list] : global_id_mapping)
  {
    const auto& global_id_list = ghost_cell_ids_consolidated.at(pid);

    if (mapping_list.size() != global_id_list.size())
      throw std::logic_error(fname + ": Ghost cell mapping error.");

    const size_t list_size = mapping_list.size();
    for (size_t k = 0; k < list_size; ++k)
      neighbor_cell_block_address_.emplace_back(global_id_list[k],
                                                static_cast<int64_t>(mapping_list[k]));
  }

  // Print info
  Chi::log.LogAllVerbose2() << "Local dof count, start, total " << local_node_count << " "
                            << local_block_address_ << " " << global_node_count;
}

void
PieceWiseLinearDiscontinuous::BuildSparsityPattern(
  std::vector<int64_t>& nodal_nnz_in_diag,
  std::vector<int64_t>& nodal_nnz_off_diag,
  const chi_math::UnknownManager& unknown_manager) const
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL CONNECTIVITY
  size_t local_dof_count = local_base_block_size_;

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_dof_count, 0);
  nodal_nnz_off_diag.resize(local_dof_count, 0);

  int lc = 0;
  for (const auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);
    size_t num_nodes = cell_mapping.NumNodes();

    // Self connection
    for (int i = 0; i < num_nodes; ++i)
    {
      int64_t ir = cell_local_block_address_[lc] + i;
      nodal_nnz_in_diag[ir] += static_cast<int64_t>(num_nodes);
    }

    // Local adjacent cell connections
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_ and face.IsNeighborLocal(ref_grid_))
      {
        const auto& adj_cell = ref_grid_.cells[face.neighbor_id_];
        const auto& adj_cell_mapping = GetCellMapping(adj_cell);

        for (int i = 0; i < num_nodes; ++i)
        {
          int64_t ir = cell_local_block_address_[lc] + i;
          nodal_nnz_in_diag[ir] += static_cast<int64_t>(adj_cell_mapping.NumNodes());
        }
      }
    } // for face
    ++lc;
  } // for local cell

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEIGHBORING CONNECTIVITY
  lc = 0;
  for (auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);

    // Local adjacent cell connections
    for (auto& face : cell.faces_)
    {
      if (face.has_neighbor_ and (not face.IsNeighborLocal(ref_grid_)))
      {
        const auto& adj_cell = ref_grid_.cells[face.neighbor_id_];
        const auto& adj_cell_mapping = GetCellMapping(adj_cell);

        for (int i = 0; i < cell_mapping.NumNodes(); ++i)
        {
          int64_t ir = cell_local_block_address_[lc] + i;
          nodal_nnz_off_diag[ir] += static_cast<int64_t>(adj_cell_mapping.NumNodes());
        }
      }
    }
    ++lc;
  } // for local cell

  // Spacing according to unknown manager
  auto backup_nnz_in_diag = nodal_nnz_in_diag;
  auto backup_nnz_off_diag = nodal_nnz_off_diag;

  unsigned int N = unknown_manager.GetTotalUnknownStructureSize();

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_base_block_size_ * N, 0);
  nodal_nnz_off_diag.resize(local_base_block_size_ * N, 0);

  if (unknown_manager.dof_storage_type_ == chi_math::UnknownStorageType::NODAL)
  {
    int ir = -1;
    for (int i = 0; i < local_base_block_size_; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        ++ir;
        nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
        nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
      } // for j
    }   // for i
  }
  else if (unknown_manager.dof_storage_type_ == chi_math::UnknownStorageType::BLOCK)
  {
    int ir = -1;
    for (int j = 0; j < N; ++j)
    {
      for (int i = 0; i < local_base_block_size_; ++i)
      {
        ++ir;
        nodal_nnz_in_diag[ir] = backup_nnz_in_diag[i];
        nodal_nnz_off_diag[ir] = backup_nnz_off_diag[i];
      } // for i
    }   // for j
  }

  Chi::mpi.Barrier();
}

int64_t
PieceWiseLinearDiscontinuous::MapDOF(const chi_mesh::Cell& cell,
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
      int64_t address = sc_int64(local_block_address_ * num_unknowns) +
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
      Chi::log.LogAllError() << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
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
      int64_t address = sc_int64(neighbor_cell_block_address_[index].second * num_unknowns) +
                        node * num_unknowns + block_id;
      return address;
    }
  }

  return -1;
}

int64_t
PieceWiseLinearDiscontinuous::MapDOFLocal(const chi_mesh::Cell& cell,
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
      int64_t address = sc_int64(cell_local_block_address_[cell.local_id_] * num_unknowns) +
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
      Chi::log.LogAllError() << "SpatialDiscretization_PWL::MapDFEMDOF. Mapping failed for cell "
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
      int64_t address = sc_int64(neighbor_cell_block_address_[index].second * num_unknowns) +
                        node * num_unknowns + block_id;
      return address;
    }
  }

  return -1;
}

size_t
PieceWiseLinearDiscontinuous::GetNumGhostDOFs(const UnknownManager& unknown_manager) const
{
  return 0;
}

std::vector<int64_t>
PieceWiseLinearDiscontinuous::GetGhostDOFIndices(const UnknownManager& unknown_manager) const
{
  return {};
}

} // namespace chi_math::spatial_discretization
