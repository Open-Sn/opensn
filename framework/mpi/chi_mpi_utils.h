#pragma once

#include "opensn/framework/mpi/chi_mpi_utils_map_all2all.h"

namespace chi_mpi_utils
{
/**Returns the current rank on the specified communicator.*/
int GetLocationID(MPI_Comm mpi_comm);
/**Returns the total number of ranks on the specified communicator.*/
int GetProcessCount(MPI_Comm mpi_comm);

/**Given each location's local size (of items), builds a vector (dimension
 * comm-size plus 1) of where each location's global indices start and end.
 * Example: location i starts at extents[i] and ends at extents[i+1]*/
std::vector<uint64_t> BuildLocationExtents(uint64_t local_size, MPI_Comm comm);
} // namespace chi_mpi_utils
