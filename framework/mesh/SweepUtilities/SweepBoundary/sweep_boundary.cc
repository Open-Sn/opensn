#include "sweep_boundaries.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Returns a pointer to a heterogeneous flux storage location.*/
double*
chi_mesh::sweep_management::SweepBoundary::HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                                                    unsigned int face_num,
                                                                    unsigned int fi,
                                                                    unsigned int angle_num,
                                                                    int group_num,
                                                                    size_t gs_ss_begin)
{
  Chi::log.LogAllError() << "HeterogeneousPsiIncoming call made to boundary "
                            "that has no such information.";
  Chi::Exit(EXIT_FAILURE);
  return nullptr;
}

//###################################################################
/**Returns a pointer to a heterogeneous flux storage location.*/
double*
chi_mesh::sweep_management::SweepBoundary::HeterogeneousPsiOutgoing(uint64_t cell_local_id,
                                                                    unsigned int face_num,
                                                                    unsigned int fi,
                                                                    unsigned int angle_num,
                                                                    size_t gs_ss_begin)
{
  Chi::log.LogAllError() << "HeterogeneousPsiOutgoing call made to boundary "
                            "that has no such information.";
  Chi::Exit(EXIT_FAILURE);
  return nullptr;
}