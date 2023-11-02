#include "LagrangeContinuous.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include <algorithm>

#define sc_int64 static_cast<int64_t>

namespace chi_math::spatial_discretization
{

// ###################################################################
/**Builds the sparsity pattern for a Continuous Finite Element Method.*/
void LagrangeContinuous::BuildSparsityPattern(
  std::vector<int64_t>& nodal_nnz_in_diag,
  std::vector<int64_t>& nodal_nnz_off_diag,
  const chi_math::UnknownManager& unknown_manager) const
{
  //**************************************** DEFINE UTILITIES

  //=================================== Dof-handler
  /**Utility mappings*/
  struct DOFHandler
  {
    const int64_t local_block_start = 0;
    const int64_t local_block_end = 0;
    const std::vector<uint64_t>& locI_block_addr;

    DOFHandler(int64_t block_start,
               int64_t block_end,
               const std::vector<uint64_t>& in_locJ_block_address)
      : local_block_start(block_start),
        local_block_end(block_end),
        locI_block_addr(in_locJ_block_address)
    {
    }

    bool IsMapLocal(int64_t ir) const
    {
      return (ir >= local_block_start and ir < local_block_end);
    }

    int64_t MapIRLocal(int64_t ir) const { return ir - local_block_start; }

    int64_t GetLocFromIR(int64_t ir)
    {
      int64_t locI =
        std::upper_bound(locI_block_addr.begin(), locI_block_addr.end(), ir) -
        locI_block_addr.begin() - 1;
      return locI;
    }

  } dof_handler(sc_int64(local_block_address_),
                sc_int64(local_block_address_ + local_base_block_size_),
                locJ_block_address_);

  // Writes a message on ir error
  auto IR_MAP_ERROR = []()
  {
    Chi::log.LogAllError() << "PWL-MapCFEMDOF: ir Mapping error node ";
    Chi::Exit(EXIT_FAILURE);
  };

  // Writes a message on jr error
  auto JR_MAP_ERROR = []()
  {
    Chi::log.LogAllError() << "PWL-MapCFEMDOF: jr Mapping error node ";
    Chi::Exit(EXIT_FAILURE);
  };

  // Checks whether an integer is already in a vector
  auto IS_VALUE_IN_VECTOR = [](const std::vector<int64_t>& vec, int64_t val)
  {
    bool already_there = false;
    for (auto check_val : vec)
      if (check_val == val)
      {
        already_there = true;
        break;
      }
    return already_there;
  };

  //**************************************** END OF UTILITIES

  //======================================== Build local sparsity pattern
  Chi::log.Log0Verbose1() << "Building local sparsity pattern.";
  std::vector<std::vector<int64_t>> nodal_connections(local_base_block_size_);

  nodal_nnz_in_diag.clear();
  nodal_nnz_off_diag.clear();

  nodal_nnz_in_diag.resize(local_base_block_size_, 0);
  nodal_nnz_off_diag.resize(local_base_block_size_, 0);

  for (auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);
    for (unsigned int i = 0; i < cell_mapping.NumNodes(); ++i)
    {
      const int64_t ir = MapDOF(cell, i);
      if (ir < 0) IR_MAP_ERROR();

      if (dof_handler.IsMapLocal(ir))
      {
        const int64_t il = dof_handler.MapIRLocal(ir);
        std::vector<int64_t>& node_links = nodal_connections[il];

        for (unsigned int j = 0; j < cell_mapping.NumNodes(); ++j)
        {
          const int64_t jr = MapDOF(cell, j);
          if (jr < 0) JR_MAP_ERROR();

          if (IS_VALUE_IN_VECTOR(node_links, jr)) continue;

          node_links.push_back(jr);
          if (dof_handler.IsMapLocal(jr)) nodal_nnz_in_diag[il] += 1;
          else
            nodal_nnz_off_diag[il] += 1;
        } // for j
      }   // if i local
    }     // for i
  }       // for cell

  //======================================== Build non-local sparsity pattern
  Chi::log.Log0Verbose1() << "Building non-local sparsity pattern.";

  // In this process we build a list
  // of ir-nodes that are not local. Each ir-node needs to
  // be furnished with the jr-nodes it links to.

  typedef std::pair<int64_t, std::vector<int64_t>> ROWJLINKS;
  std::vector<ROWJLINKS> ir_links;

  for (auto& cell : ref_grid_.local_cells)
  {
    const auto& cell_mapping = GetCellMapping(cell);

    for (unsigned int i = 0; i < cell_mapping.NumNodes(); ++i)
    {
      const int64_t ir = MapDOF(cell, i);
      if (ir < 0) IR_MAP_ERROR();

      if (not dof_handler.IsMapLocal(ir))
      {
        ROWJLINKS new_ir_link;
        ROWJLINKS* cur_ir_link = &new_ir_link;

        //============================= Check if ir already there
        bool ir_already_there = false;
        for (auto& ir_link : ir_links)
          if (ir == ir_link.first)
          {
            ir_already_there = true;
            cur_ir_link = &ir_link;
            break;
          }

        //============================= Now add links
        auto& node_links = cur_ir_link->second;
        for (unsigned int j = 0; j < cell_mapping.NumNodes(); ++j)
        {
          const int64_t jr = MapDOF(cell, j);
          if (jr < 0) JR_MAP_ERROR();

          if (IS_VALUE_IN_VECTOR(node_links, jr)) continue;
          else
            node_links.push_back(jr);
        } // for j

        //============================= If its not add it
        if (not ir_already_there)
        {
          cur_ir_link->first = ir;
          ir_links.push_back(*cur_ir_link);
        }

      } // if i not local
    }   // for i
  }     // for cell

  //======================================== Build communication structure
  Chi::log.Log0Verbose1() << "Building communication structure.";

  //=================================== Step 1
  // We now serialize the non-local data
  std::vector<std::vector<int64_t>> locI_serialized(Chi::mpi.process_count);

  for (const auto& ir_linkage : ir_links)
  {
    const int64_t locI = dof_handler.GetLocFromIR(ir_linkage.first);

    locI_serialized[locI].push_back(
      sc_int64(ir_linkage.second.size())); // row cols amount
    locI_serialized[locI].push_back(sc_int64(ir_linkage.first)); // row num
    for (int64_t jr : ir_linkage.second)
      locI_serialized[locI].push_back(jr); // col num
  }

  //=================================== Step 2
  // Establish the size of the serialized data
  // to send to each location and communicate
  // to get receive count.
  std::vector<int> sendcount(Chi::mpi.process_count, 0);
  std::vector<int> recvcount(Chi::mpi.process_count, 0);
  int locI = 0;
  for (const auto& locI_data : locI_serialized)
  {
    sendcount[locI] = static_cast<int>(locI_data.size());

    if (Chi::mpi.location_id == 0)
      Chi::log.LogAllVerbose1()
        << "To send to " << locI << " = " << sendcount[locI];

    ++locI;
  }

  MPI_Alltoall(
    sendcount.data(), 1, MPI_INT, recvcount.data(), 1, MPI_INT, Chi::mpi.comm);

  //=================================== Step 3
  // We now establish send displacements and
  // receive displacements.
  std::vector<int> send_displs(Chi::mpi.process_count, 0);
  std::vector<int> recv_displs(Chi::mpi.process_count, 0);

  int send_displ_c = 0;
  int recv_displ_c = 0;

  int c = 0;
  for (int send_count : sendcount)
  {
    send_displs[c++] = send_displ_c;
    send_displ_c += send_count;
  }
  c = 0;
  for (int recv_count : recvcount)
  {
    recv_displs[c++] = recv_displ_c;
    recv_displ_c += recv_count;
  }

  //======================================== Communicate data
  Chi::log.Log0Verbose1() << "Communicating non-local rows.";

  // We now initialize the buffers and
  // communicate the data
  std::vector<int64_t> sendbuf;
  std::vector<int64_t> recvbuf;

  sendbuf.reserve(send_displ_c);
  recvbuf.resize(recv_displ_c, 0);

  for (const auto& serial_block : locI_serialized)
    for (int64_t data_val : serial_block)
      sendbuf.push_back(data_val);

  MPI_Alltoallv(sendbuf.data(),
                sendcount.data(),
                send_displs.data(),
                MPI_INT64_T,
                recvbuf.data(),
                recvcount.data(),
                recv_displs.data(),
                MPI_INT64_T,
                Chi::mpi.comm);

  //======================================== Deserialze data
  Chi::log.Log0Verbose1() << "Deserialize data.";

  std::vector<ROWJLINKS> foreign_ir_links;

  for (size_t k = 0; k < recvbuf.size();)
  {
    const int64_t num_values = recvbuf[k++];
    const int64_t ir = recvbuf[k++];

    ROWJLINKS new_links;
    new_links.first = ir;
    new_links.second.reserve(num_values);
    for (int i = 0; i < num_values; ++i)
      new_links.second.push_back(recvbuf[k++]);

    foreign_ir_links.push_back(std::move(new_links));
  }

  //======================================== Adding to sparsity pattern
  for (const auto& ir_linkage : foreign_ir_links)
  {
    const int64_t ir = ir_linkage.first;

    if (not dof_handler.IsMapLocal(ir)) IR_MAP_ERROR();

    int64_t il = dof_handler.MapIRLocal(ir);
    std::vector<int64_t>& node_links = nodal_connections[il];

    for (int64_t jr : ir_linkage.second)
    {
      if (IS_VALUE_IN_VECTOR(node_links, jr)) continue;

      node_links.push_back(jr);
      if (dof_handler.IsMapLocal(jr)) nodal_nnz_in_diag[il] += 1;
      else
        nodal_nnz_off_diag[il] += 1;
    }
  }

  Chi::mpi.Barrier();

  //======================================== Spacing according to unknown
  //                                         manager
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
  else if (unknown_manager.dof_storage_type_ ==
           chi_math::UnknownStorageType::BLOCK)
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
}

} // namespace chi_math::spatial_discretization