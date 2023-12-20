#include "framework/mesh/sweep_utilities/sweep_boundary/sweep_boundary.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mpi/mpi.h"

namespace opensn
{

double*
SweepBoundary::HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                        unsigned int face_num,
                                        unsigned int fi,
                                        unsigned int angle_num,
                                        int group_num,
                                        size_t gs_ss_begin)
{
  Chi::log.LogAllError()
    << "HeterogeneousPsiIncoming call made to boundary that has no such information.";
  Chi::Exit(EXIT_FAILURE);
  return nullptr;
}

double*
SweepBoundary::HeterogeneousPsiOutgoing(uint64_t cell_local_id,
                                        unsigned int face_num,
                                        unsigned int fi,
                                        unsigned int angle_num,
                                        size_t gs_ss_begin)
{
  Chi::log.LogAllError()
    << "HeterogeneousPsiOutgoing call made to boundary that has no such information.";
  Chi::Exit(EXIT_FAILURE);
  return nullptr;
}

} // namespace opensn
