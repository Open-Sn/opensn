#include "FiniteVolume.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
#include "chi_mpi.h"

#include "mpi/chi_mpi_utils_map_all2all.h"

#define MappingError                                                                               \
  "chi_math::SpatialDiscretization_FV::OrderNodes: "                                               \
  "Error mapping neighbor cells"

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Develops node ordering per location.*/
void
chi_math::spatial_discretization::FiniteVolume::OrderNodes()
{
  //============================================= Communicate node counts
  const uint64_t local_num_nodes = ref_grid_.local_cells.size();
  locJ_block_size_.assign(Chi::mpi.process_count, 0);
  MPI_Allgather(&local_num_nodes, // sendbuf
                1,
                MPI_UINT64_T,            // sendcount, sendtype
                locJ_block_size_.data(), // recvbuf
                1,
                MPI_UINT64_T,   // recvcount, recvtype
                Chi::mpi.comm); // comm

  //============================================= Build block addresses
  locJ_block_address_.assign(Chi::mpi.process_count, 0);
  uint64_t global_num_nodes = 0;
  for (int j = 0; j < Chi::mpi.process_count; ++j)
  {
    locJ_block_address_[j] = global_num_nodes;
    global_num_nodes += locJ_block_size_[j];
  }

  local_block_address_ = locJ_block_address_[Chi::mpi.location_id];

  local_base_block_size_ = local_num_nodes;
  globl_base_block_size_ = global_num_nodes;

  //============================================= Sort neigbor ids
  const auto neighbor_gids = ref_grid_.cells.GetGhostGlobalIDs();
  std::map<uint64_t, std::vector<uint64_t>> sorted_nb_gids;
  for (uint64_t gid : neighbor_gids)
  {
    const auto& cell = ref_grid_.cells[gid];
    sorted_nb_gids[cell.partition_id_].push_back(gid);
  }

  //============================================= Communicate neighbor ids
  //                                              requiring mapping
  const auto query_nb_gids = chi_mpi_utils::MapAllToAll(sorted_nb_gids, // map
                                                        MPI_UINT64_T,   // datatype
                                                        Chi::mpi.comm); // comm

  //============================================= Map the ids
  std::map<uint64_t, std::vector<uint64_t>> mapped_query_nb_gids;
  for (const auto& pid_list_pair : query_nb_gids)
  {
    const uint64_t pid = pid_list_pair.first;
    const auto& gids = pid_list_pair.second;

    auto& local_ids = mapped_query_nb_gids[pid];
    local_ids.reserve(gids.size());
    for (uint64_t gid : gids)
    {
      if (not ref_grid_.IsCellLocal(gid)) throw std::logic_error(MappingError);

      const auto& local_cell = ref_grid_.cells[gid];
      local_ids.push_back(local_cell.local_id_);
    } // for gid
  }   // for pid_list_pair

  //============================================= Communicate back the mapped
  //                                              ids
  const auto mapped_nb_gids = chi_mpi_utils::MapAllToAll(mapped_query_nb_gids, // map
                                                         MPI_UINT64_T,         // datatype
                                                         Chi::mpi.comm);       // comm

  //============================================= Create the neighbor cell
  //                                              mapping
  neighbor_cell_local_ids_.clear();
  for (const auto& pid_list_pair : sorted_nb_gids)
  {
    try
    {
      const auto& pid = pid_list_pair.first;
      const auto& gid_list = pid_list_pair.second;
      const auto& lid_list = mapped_nb_gids.at(pid);

      if (gid_list.size() != lid_list.size())
        throw std::logic_error(MappingError + std::string(" Size-mismatch."));

      for (size_t i = 0; i < gid_list.size(); ++i)
        neighbor_cell_local_ids_.insert(std::make_pair(gid_list[i], lid_list[i]));
    }
    catch (const std::out_of_range& oor)
    {
      throw std::logic_error(MappingError + std::string(" OOR."));
    }
  } // for pid_list_pair

  local_base_block_size_ = ref_grid_.local_cells.size();
  globl_base_block_size_ = ref_grid_.GetGlobalNumberOfCells();
}

} // namespace chi_math::spatial_discretization
