#include "PieceWiseLinearDiscontinuous.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
#include "chi_mpi.h"

#include "utils/chi_timer.h"

#include "chi_mpi_utils.h"

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Reorders the nodes for parallel computation in a Continuous
 * Finite Element calculation.*/
void PieceWiseLinearDiscontinuous::OrderNodes()
{
  const std::string fname = __FUNCTION__;
  chi::Timer t_stage[6];

  t_stage[0].Reset();
  //================================================== Check cell views avail
  size_t num_loc_cells = ref_grid_.local_cells.size();

  //================================================== Get local DOF count
  //                                                   and set
  //                                                   cell_local_block_address
  cell_local_block_address_.resize(num_loc_cells, 0);

  uint64_t local_node_count = 0;
  for (const auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);
    cell_local_block_address_[cell.local_id_] =
      static_cast<int64_t>(local_node_count);
    local_node_count += cell_mapping.NumNodes();
  }

  //================================================== Allgather node_counts
  locJ_block_size_.assign(Chi::mpi.process_count, 0);
  MPI_Allgather(&local_node_count, // sendbuf
                1,
                MPI_UNSIGNED_LONG_LONG,  // sendcount, sendtype
                locJ_block_size_.data(), // recvbuf
                1,
                MPI_UNSIGNED_LONG_LONG, // recvcount, recvtype
                Chi::mpi.comm);         // comm

  //================================================== Assign
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

  //================================================== Collect ghost cell ids
  //                                                   needing block addresses
  std::map<int, std::vector<uint64_t>> ghost_cell_ids_consolidated;

  for (uint64_t global_id : ref_grid_.cells.GetGhostGlobalIDs())
  {
    const auto& cell = ref_grid_.cells[global_id];
    const int locI = static_cast<int>(cell.partition_id_);

    std::vector<uint64_t>& locI_cell_id_list =
      ghost_cell_ids_consolidated[locI];

    locI_cell_id_list.push_back(cell.global_id_);
  }

  //================================================== AllToAll to get query
  //                                                   cell-ids
  const std::map<int, std::vector<uint64_t>> query_ghost_cell_ids_consolidated =
    chi_mpi_utils::MapAllToAll(ghost_cell_ids_consolidated,
                               MPI_UNSIGNED_LONG_LONG);

  //================================================== Map all query cell-ids
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

  //================================================== Communicate back the
  // mapping
  const std::map<int, std::vector<uint64_t>> global_id_mapping =
    chi_mpi_utils::MapAllToAll(mapped_ghost_cell_ids_consolidated,
                               MPI_UNSIGNED_LONG_LONG);

  //================================================== Process global id mapping
  for (const auto& [pid, mapping_list] : global_id_mapping)
  {
    const auto& global_id_list = ghost_cell_ids_consolidated.at(pid);

    if (mapping_list.size() != global_id_list.size())
      throw std::logic_error(fname + ": Ghost cell mapping error.");

    const size_t list_size = mapping_list.size();
    for (size_t k = 0; k < list_size; ++k)
      neighbor_cell_block_address_.emplace_back(
        global_id_list[k], static_cast<int64_t>(mapping_list[k]));
  }

  //================================================== Print info
  Chi::log.LogAllVerbose2()
    << "Local dof count, start, total " << local_node_count << " "
    << local_block_address_ << " " << global_node_count;
}

} // namespace chi_math::spatial_discretization
