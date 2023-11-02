#include "chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "chi_mpi.h"

//###################################################################
/**Builds and returns a vector of unique boundary id's present in
 * the mesh.*/
std::vector<uint64_t>
chi_mesh::MeshContinuum::GetDomainUniqueBoundaryIDs() const
{
  Chi::mpi.Barrier();
  Chi::log.Log() << "Identifying unique boundary-ids.";

  //====================================== Develop local bndry-id set
  std::set<uint64_t> local_bndry_ids_set;
  for (auto& cell : local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_) local_bndry_ids_set.insert(face.neighbor_id_);

  //====================================== Vectorify it and get local count
  std::vector<uint64_t> local_bndry_ids(local_bndry_ids_set.begin(), local_bndry_ids_set.end());
  int local_num_bndry_ids = (int)local_bndry_ids.size();

  //====================================== Everyone now tells everyone
  //                                       how many bndry-ids they have
  std::vector<int> locI_bndry_count(Chi::mpi.process_count, 0);

  MPI_Allgather(&local_num_bndry_ids,    // sendbuf
                1,                       // sendcount
                MPI_INT,                 // sendtype
                locI_bndry_count.data(), // recvbuf
                1,                       // recvcount
                MPI_INT,                 // recvtype
                Chi::mpi.comm);          // communicator

  //====================================== Build a displacement list, in prep
  //                                       for gathering all bndry-ids
  std::vector<int> locI_bndry_ids_displs(Chi::mpi.process_count, 0);
  size_t total_num_global_bndry_ids = locI_bndry_count[0];
  for (int locI = 1; locI < Chi::mpi.process_count; ++locI)
  {
    locI_bndry_ids_displs[locI] = locI_bndry_ids_displs[locI - 1] + locI_bndry_count[locI - 1];
    total_num_global_bndry_ids += locI_bndry_count[locI];
  }

  //====================================== Everyone now sends everyone
  //                                       they're boundary-ids
  std::vector<uint64_t> globl_bndry_ids(total_num_global_bndry_ids);

  MPI_Allgatherv(local_bndry_ids.data(),       // sendbuf
                 local_num_bndry_ids,          // sendcount
                 MPI_UNSIGNED_LONG_LONG,       // sendtype
                 globl_bndry_ids.data(),       // recvbuf
                 locI_bndry_count.data(),      // recvcounts
                 locI_bndry_ids_displs.data(), // displs
                 MPI_UNSIGNED_LONG_LONG,       // recvtype
                 Chi::mpi.comm);               // communicator

  std::set<uint64_t> globl_bndry_ids_set(globl_bndry_ids.begin(), globl_bndry_ids.end());

  std::vector<uint64_t> unique_bdnry_ids(globl_bndry_ids_set.begin(), globl_bndry_ids_set.end());
  return unique_bdnry_ids;
}