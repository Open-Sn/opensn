// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

void
CommunicateLocationDependencies(const std::vector<int>& location_dependencies,
                                std::vector<std::vector<int>>& global_dependencies)
{
  CALI_CXX_MARK_FUNCTION;

  int P = opensn::mpi_comm.size();

  // Communicate location dep counts
  std::vector<int> depcount_per_loc;
  int current_loc_dep_count = location_dependencies.size();
  mpi_comm.all_gather(current_loc_dep_count, depcount_per_loc);

  // Broadcast dependencies
  std::vector<int> raw_depvec_displs(P, 0);
  for (int locI = 1; locI < P; ++locI)
    raw_depvec_displs[locI] = raw_depvec_displs[locI - 1] + depcount_per_loc[locI - 1];

  std::vector<int> raw_dependencies;
  mpi_comm.all_gather(location_dependencies, raw_dependencies, depcount_per_loc, raw_depvec_displs);

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

} // namespace opensn
