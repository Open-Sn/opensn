#include "framework/mpi/chi_mpi_utils.h"

namespace chi_mpi_utils
{

int
GetLocationID(MPI_Comm mpi_comm)
{
  int location_id;
  MPI_Comm_rank(mpi_comm, &location_id);

  return location_id;
}

int
GetProcessCount(MPI_Comm mpi_comm)
{
  int process_count;
  MPI_Comm_size(mpi_comm, &process_count);

  return process_count;
}

std::vector<uint64_t>
BuildLocationExtents(uint64_t local_size, MPI_Comm comm)
{
  const int process_count = GetProcessCount(comm);
  // Get the local vector sizes per process
  std::vector<uint64_t> local_sizes(process_count, 0);
  MPI_Allgather(&local_size, 1, MPI_UINT64_T, local_sizes.data(), 1, MPI_UINT64_T, comm);

  // With the vector sizes per processor, now the offsets for each
  // processor can be defined using a cumulative sum per processor.
  // This allows for the determination of whether a global index is
  // locally owned or not.
  std::vector<uint64_t> extents(process_count + 1, 0);
  for (size_t locJ = 1; locJ < process_count; ++locJ)
    extents[locJ] = extents[locJ - 1] + local_sizes[locJ - 1];
  extents[process_count] = extents[process_count - 1] + local_sizes.back();

  return extents;
}
} // namespace chi_mpi_utils
