#include "lbsadj_solver.h"

#include "mesh/LogicalVolume/LogicalVolume.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs
{

void
DiscreteOrdinatesAdjointSolver::InitQOIs()
{
  //============================================= Initialize QOIs
  for (auto& qoi_pair : response_functions_)
  {
    const auto& qoi_designation = qoi_pair.first;
    auto& qoi_cell_subscription = qoi_pair.second;

    for (const auto& cell : grid_ptr_->local_cells)
      if (qoi_designation.logical_volume->Inside(cell.centroid_))
        qoi_cell_subscription.push_back(cell.local_id_);

    size_t num_local_subs = qoi_cell_subscription.size();
    size_t num_globl_subs = 0;

    MPI_Allreduce(&num_local_subs, // sendbuf
                  &num_globl_subs, // recvbuf
                  1,
                  MPI_UNSIGNED_LONG_LONG, // count + datatype
                  MPI_SUM,                // operation
                  Chi::mpi.comm);         // communicator

    Chi::log.Log() << "LBAdjointSolver: Number of cells subscribed to " << qoi_designation.name
                   << " = " << num_globl_subs;
  }
}

} // namespace lbs