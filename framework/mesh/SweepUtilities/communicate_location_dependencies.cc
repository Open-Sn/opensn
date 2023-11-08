#include "framework/mesh/SweepUtilities/sweep_namespace.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "framework/mpi/chi_mpi.h"

void
chi_mesh::sweep_management::CommunicateLocationDependencies(
  const std::vector<int>& location_dependencies, std::vector<std::vector<int>>& global_dependencies)
{
  int P = Chi::mpi.process_count;

  // Communicate location dep
  // counts
  std::vector<int> depcount_per_loc(P, 0);
  int current_loc_dep_count = location_dependencies.size();
  MPI_Allgather(
    &current_loc_dep_count, 1, MPI_INT, depcount_per_loc.data(), 1, MPI_INT, Chi::mpi.comm);

  // Broadcast dependencies
  std::vector<int> raw_depvec_displs(P, 0);
  int recv_buf_size = depcount_per_loc[0];
  for (int locI = 1; locI < P; ++locI)
  {
    raw_depvec_displs[locI] = raw_depvec_displs[locI - 1] + depcount_per_loc[locI - 1];
    recv_buf_size += depcount_per_loc[locI];
  }

  std::vector<int> raw_dependencies(recv_buf_size, 0);

  MPI_Allgatherv(location_dependencies.data(),
                 int(location_dependencies.size()),
                 MPI_INT,
                 raw_dependencies.data(),
                 depcount_per_loc.data(),
                 raw_depvec_displs.data(),
                 MPI_INT,
                 Chi::mpi.comm);

  for (int locI = 0; locI < P; ++locI)
  {
    global_dependencies[locI].resize(depcount_per_loc[locI], 0);
    for (int c = 0; c < depcount_per_loc[locI]; ++c)
    {
      int addr = raw_depvec_displs[locI] + c;
      global_dependencies[locI][c] = raw_dependencies[addr];
    }
  }
}
