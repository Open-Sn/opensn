#include "lbs_discrete_ordinates_solver.h"

#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include "chi_runtime.h"
#include "console/chi_console.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include <iomanip>

//###################################################################
/**Clears all the sweep orderings for a groupset in preperation for
 * another.*/
void
lbs::DiscreteOrdinatesSolver::ResetSweepOrderings(LBSGroupset& groupset)
{
  Chi::log.Log0Verbose1() << "Resetting SPDS and FLUDS";

  groupset.angle_agg_->angle_set_groups.clear();

  Chi::mpi.Barrier();

  Chi::log.Log() << "SPDS and FLUDS reset complete.            Process memory = "
                 << std::setprecision(3) << chi::Console::GetMemoryUsageInMB() << " MB";

  double local_app_memory = Chi::log.ProcessEvent(chi::ChiLog::StdTags::MAX_MEMORY_USAGE,
                                                  chi::ChiLog::EventOperation::MAX_VALUE);
  double total_app_memory = 0.0;
  MPI_Allreduce(&local_app_memory, &total_app_memory, 1, MPI_DOUBLE, MPI_SUM, Chi::mpi.comm);
  double max_proc_memory = 0.0;
  MPI_Allreduce(&local_app_memory, &max_proc_memory, 1, MPI_DOUBLE, MPI_MAX, Chi::mpi.comm);

  Chi::log.Log() << "\n"
                 << std::setprecision(3)
                 << "           Total application memory (max): " << total_app_memory / 1024.0
                 << " GB\n"
                 << "           Maximum process memory        : " << max_proc_memory / 1024.0
                 << " GB\n\n";
}